!**********************************************************************
module rns_ls
!**********************************************************************
use rns_base_ls
!!$if($CYL_SKEW_LS)
!!use cyl_skew_base_ls, only : tr_t, clindx_to_loc_id, brindx_to_loc_id, ntree
!!use cyl_skew_ls, only : cyl_skew_fill_tree_array_ls, ngen, ngen_reslv
!!!use cyl_skew_ls, only : cyl_skew_fill_cl_ref_plane_array_ls
!!!use cyl_skew_ls, only : cyl_skew_get_branch_id_ls
!!$endif

implicit none

save
private

public :: rns_init_ls !, rns_u_write_ls
public :: rns_CD_ls
public :: rns_forcing_ls ! Apply forcing

character (*), parameter :: mod_name = 'rns_ls'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************
use messages
use param, only : USE_MPI, coord
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : tree, tr_t, ntree
use cyl_skew_ls, only : cyl_skew_fill_tree_array_ls
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'

integer :: nt, np, ng, nc

integer :: ncluster_tot, ncount

! Load brindx
call brindx_init()
! Load filtered indicator function (chi)
call chi_init()

$if($CYL_SKEW_LS)
if(coord == 0) call mesg ( sub_name, 'filling tree array' )
call cyl_skew_fill_tree_array_ls()
if(coord == 0) call mesg ( sub_name, 'tree array filled' )
$endif

!if(coord == 0) write(*,*) 'ncluster_tot : ', ncluster_tot

!if(use_main_tree_ref) then
!  !ncluster_ref = tr_t(1) % ncluster
!  ntree_ref = 1
!else
!  !ncluster_ref = ncluster_tot
!  ntree_ref = ntree
!endif
if ( rns_ntree > ntree ) call error ( sub_name, 'rns_ntree > ntree ')
!  this is used to map the brindx to correct rns tree
allocate( rns_tree_iarray( ntree ) )

if(rns_tree_layout == 1) then 
  do nt = 1, ntree
  
    if( nt < rns_ntree ) then
        rns_tree_iarray( nt ) = nt
    else
        rns_tree_iarray(nt) = rns_ntree
    endif
    
  enddo
endif


!if(coord == 0) write(*,*) 'ncluster_ref, ntree_ref : ', ncluster_ref, ntree_ref

ncluster_reslv = 0
do nt = 1, rns_ntree
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      if(ng <= tr_t(nt) % ngen_reslv) ncluster_reslv = ncluster_reslv + 1   
    enddo    
  enddo
enddo

!  Get the total number of clusters 
ncluster_tot = 0
do nt = 1, ntree
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      ncluster_tot = ncluster_tot + 1   
    enddo    
  enddo
enddo

allocate( rns_reslv_cl_iarray( ncluster_tot ) )
ncount=0
do nt = 1, rns_ntree
  do ng = 1, tr_t(nt) % ngen_reslv
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      ncount = ncount + 1
      rns_reslv_cl_iarray( tr_t(nt) % gen_t (ng) % cl_t (nc) %indx ) = ncount 
    enddo    
  enddo
enddo

!ncluster_unreslv = ncluster_tot - ncluster_reslv

if(coord == 0) write(*,*) 'ncluster_reslv : ', ncluster_reslv

if ( .not. use_beta_sub_regions ) then
  nbeta = rns_ntree
else
  call error(sub_name, 'use_beta_sub_regions must be false')  
endif

!if(use_main_tree_ref) then
!  ncluster_unreslv_ref = ncluster_ref - ncluster_reslv_ref
!else
!  ncluster_unreslv_ref = ncluster_unreslv
!endif

!if(coord == 0) write(*,*) 'ncluster_unreslv : ', ncluster_unreslv

!if(coord == 0) then

!  write(*,*) 'ncluster_tot : ', ncluster_tot
!  write(*,*) 'ncluster_reslv : ', ncluster_reslv
!  write(*,*) 'ncluster_unreslv : ', ncluster_unreslv
!  
!endif

!!  Set the number of resolved clusters
!if(use_main_tree_only) then

!  ncluster_reslv = 0
!  
!  tr_t_p => tr_t(1)
!  
!  do ng = 1, ngen_reslv
!  
!    do nc = 1, tr_t_p % gen_t(ng) % ncluster
!    
!      ncluster_reslv = ncluster_reslv + 1
!      
!    enddo
!    
!  enddo
!  
!  nullify(tr_t_p)
!  
!else

!  call error(sub_name, 'Not set up yet for multiple trees')
!  
!endif

if(clforce_calc) then
  !  Create cluster reference value plane
  call rns_fill_ref_plane_array_ls()
  !  Create cluster index array
  call rns_fill_indx_array_ls()
  !  Set parent used for assigning CD's on unresolved branches
  !if(coord == 0) write(*,*) ' calling rns_set_cl_parent_ls'
  !call rns_set_cl_parent_ls()
  !if(coord == 0) write(*,*) 'called rns_set_cl_parent_ls'
  
  allocate( clforce_t ( ncluster_reslv ) ) 
  clforce_t = force(parent = 0, CD = 0._rprec, fD = 0._rprec, kappa=0._rprec)
  
  allocate( beta_force_t( nbeta ) )
  beta_force_t = force(parent = 0, CD = 0._rprec, fD = 0._rprec, kappa=0._rprec)

endif

write(*,*) 'exiting rns_init_ls'
 
return
end subroutine rns_init_ls

!**********************************************************************
subroutine rns_fill_ref_plane_array_ls()
!**********************************************************************
!  This subroutine fills all reference plane arrays
!
use types, only : rprec, vec2d
use param, only : dy, dz, USE_MPI, coord
use messages
$if($CYL_SKEW_LS) 
use cyl_skew_base_ls, only : tr_t
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_fill_cl_ref_plane_array_ls'
real(rprec), parameter :: alpha=1._rprec
real(rprec), parameter :: alpha_rbeta = 1.25_rprec
real(rprec), parameter :: alpha_beta = alpha_rbeta

integer :: nt, ng, nc, nb, rns_clindx

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)

integer, pointer :: clindx_p, nbranch_p
real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p
real(rprec), pointer, dimension(:) :: origin_p

type(vec2d) :: rvec_t

nullify(d_p, l_p, skew_angle_p, clindx_p)

allocate( cl_ref_plane_t( ncluster_reslv ) )
allocate( rbeta_ref_plane_t ( rns_ntree ) )

!if(ntree_ref < 1) call error( sub_name, 'ntree_ref not specified correctly')

do nt=1, rns_ntree

  do ng=1, tr_t(nt) % ngen_reslv
  
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
    
      nbranch_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
      
      clindx_p => rns_reslv_cl_iarray( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx )
      
      rns_clindx = rns_clindx + 1
      
      !if(clindx_p > ncluster_ref) then
      !  call mesg(sub_name, 'clindx_p : ', clindx_p)
      !  call mesg(sub_name, 'ncluster_ref : ', ncluster_ref)
      !  call error(sub_name, 'clindx_p > ncluster_ref')
      !endif
      
      h_m = 0._rprec
      area_proj = 0._rprec
      
      do nb = 1, nbranch_p

        d_p          => tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t(nb)%d
        l_p          => tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t(nb)%l
        skew_angle_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t(nb)%skew_angle
        
        h         = l_p * dcos(skew_angle_p)
        h_m       = h_m + h
        area_proj = area_proj + d_p * h
        
        nullify(d_p, l_p, skew_angle_p)
      
      enddo
      
      !  Mean height of branch cluster  and height of reference area
      h_m = h_m / nbranch_p
      !  width of reference area
      w   = area_proj / h_m     

      cl_ref_plane_t(clindx_p) % area = area_proj
      !  These are defined to be x - planes (no not the NASA experimental planes)
      cl_ref_plane_t(clindx_p) % nzeta = ceiling( w / dy + 1)
      cl_ref_plane_t(clindx_p) % neta  = ceiling( h_m / dz + 1)
      
      origin_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%origin
      
      !  Offset in the upstream x-direction
      zeta_c = origin_p + (/ -alpha * w, 0._rprec, 0._rprec /)
      
      cl_ref_plane_t(clindx_p) % p1    = zeta_c 
      cl_ref_plane_t(clindx_p) % p1(2) = cl_ref_plane_t(clindx_p) % p1(2) + w / 2._rprec
      
      cl_ref_plane_t(clindx_p) % p2    = cl_ref_plane_t(clindx_p) % p1
      cl_ref_plane_t(clindx_p) % p2(2) = cl_ref_plane_t(clindx_p) % p2(2) - w
      
      cl_ref_plane_t(clindx_p) % p3    = cl_ref_plane_t(clindx_p) % p2
      cl_ref_plane_t(clindx_p) % p3(3) = cl_ref_plane_t(clindx_p) % p3(3) + h_m
      
      nullify(nbranch_p, clindx_p, origin_p)
      
    enddo
    
  enddo
 
enddo

!  Need to set rbeta_ref_plane_t
do nt = 1, rns_ntree
  
   hbot_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen_reslv ) % bplane
   htop_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen ) % tplane
   
   origin_p => tr_t(nt) % origin
   
   h = htop_p - hbot_p
   
   !  Let w = 2 * (2D distance) from the tree origin to the center of a top most cluster
   rvec_t % xy = tr_t(nt) % gen_t ( tr_t(nt) % ngen ) % cl_t(1) % origin(1:2)
   
   rvec_t % xy = rvec_t % xy - origin_p(1:2)
   
   call vector_magnitude_2d( rvec_t % xy, rvec_t % mag )
   
   w = 2._rprec * rvec_t % mag
   
   rbeta_ref_plane_t(nt) % area = h * w 
   
   rbeta_ref_plane_t( nt ) % nzeta = ceiling( w / dy + 1)
   rbeta_ref_plane_t( nt ) % neta  = ceiling( h / dz + 1)
      
      !  Offset in the upstream x-direction and displace in z
    zeta_c = origin_p + (/ -alpha_rbeta * rvec_t % mag, 0._rprec, hbot_p /)  
    
    rbeta_ref_plane_t(nt) % p1    = zeta_c 
    rbeta_ref_plane_t(nt) % p1(2) = rbeta_ref_plane_t(nt) % p1(2) + w / 2._rprec
      
    rbeta_ref_plane_t(nt) % p2    = rbeta_ref_plane_t(nt) % p1
    rbeta_ref_plane_t(nt) % p2(2) = rbeta_ref_plane_t(nt) % p2(2) - w
      
    rbeta_ref_plane_t(nt) % p3    = rbeta_ref_plane_t(nt) % p2
    rbeta_ref_plane_t(nt) % p3(3) = rbeta_ref_plane_t(nt) % p3(3) + h
    
    nullify(hbot_p, htop_p, origin_p)
    
enddo 

!  Now need to define beta_ref_plane_t; needs to be consistent with corresponding volume of beta_indx_array
if(.not. use_beta_sub_regions) then

  !  Define entire unresolved region of tree as a beta region
  allocate( beta_ref_plane_t( nbeta ) )
  
  do nt = 1, nbeta

   hbot_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen_reslv ) % tplane
   htop_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen ) % tplane
   
   origin_p => tr_t(nt) % origin
   
   h = htop_p - hbot_p
   
   !  Let w = 2 * (2D distance) from the tree origin to the center of a top most cluster
   rvec_t % xy = tr_t(nt) % gen_t ( tr_t(nt) % ngen ) % cl_t(1) % origin(1:2)
   
   rvec_t % xy = rvec_t % xy - origin_p(1:2)
   
   call vector_magnitude_2d( rvec_t % xy, rvec_t % mag )
   
   w = 2._rprec * rvec_t % mag
   
   beta_ref_plane_t(nt) % area = h * w 
   
   beta_ref_plane_t( nt ) % nzeta = ceiling( w / dy + 1)
   beta_ref_plane_t( nt ) % neta  = ceiling( h / dz + 1)
      
      !  Offset in the upstream x-direction and displace in z
    zeta_c = origin_p + (/ -alpha_beta * rvec_t % mag, 0._rprec, hbot_p /)  
    
    beta_ref_plane_t(nt) % p1    = zeta_c 
    beta_ref_plane_t(nt) % p1(2) = beta_ref_plane_t(nt) % p1(2) + w / 2._rprec
      
    beta_ref_plane_t(nt) % p2    = beta_ref_plane_t(nt) % p1
    beta_ref_plane_t(nt) % p2(2) = beta_ref_plane_t(nt) % p2(2) - w
      
    beta_ref_plane_t(nt) % p3    = beta_ref_plane_t(nt) % p2
    beta_ref_plane_t(nt) % p3(3) = beta_ref_plane_t(nt) % p3(3) + h
    
    nullify(hbot_p, htop_p, origin_p)
    
  enddo
  
endif
    
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) 'Reference Plane Values For All Trees : '
  do nt=1, rns_ntree

    do ng = 1, tr_t(nt)%ngen_reslv
      do nc = 1, tr_t(nt)%gen_t(ng)%ncluster

        write(*,*) '-------------------------'
        write(*,*) 'nt, ng, nc : ', nt, ng, nc  
        write(*,*) 'nzeta, neta : ', cl_ref_plane_t(rns_reslv_cl_iarray( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx )) % nzeta, &
           cl_ref_plane_t(rns_reslv_cl_iarray ( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % neta
         write(*,*) 'p1 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p1
         write(*,*) 'p2 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p2
         write(*,*) 'p3 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p3
         write(*,*) 'area : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % area
        
         write(*,*) '-------------------------'
      enddo
    enddo
    enddo

endif
          
      
return

end subroutine rns_fill_ref_plane_array_ls

!**********************************************************************
subroutine rns_fill_indx_array_ls()
!**********************************************************************
!  This subroutine sets the indx_array for both resolved and unresolved 
!  regions
!
use types, only : rprec
use param, only : nx,ny,nz, coord
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : ngen, ngen_reslv, brindx_to_loc_id, tr_t
$endif
use level_set_base, only : phi
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_fill_indx_array_ls'

integer :: i,j,k, nc, np, nb
integer, pointer :: clindx_p, brindx_p
integer, pointer :: nt_p, ng_p, nc_p

integer, pointer, dimension(:) :: br_loc_id_p
type(indx_array), pointer, dimension(:) :: cl_pre_indx_array_t, beta_pre_indx_array_t

if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

! ---- Nullify all pointers ----
nullify(clindx_p, br_loc_id_p, brindx_p)
nullify(nt_p, ng_p, nc_p)
nullify(cl_pre_indx_array_t, beta_pre_indx_array_t)

! ---- Allocate all arrays ----
allocate(cl_pre_indx_array_t( ncluster_reslv ) )

if (.not. use_beta_sub_regions ) then
  allocate(beta_pre_indx_array_t ( nbeta ) )
endif

do nc=1, ncluster_reslv
  allocate(cl_pre_indx_array_t(nc) % iarray(3,nx*ny*(nz-1)))
enddo

do nb=1, nbeta
  allocate(beta_pre_indx_array_t(nb) % iarray(3,nx*ny*(nz-1)))
enddo

!  Intialize the number of points assigned to the cluster
do nc=1, ncluster_reslv
  cl_pre_indx_array_t(nc) % npoint = 0
enddo

do nb=1, nbeta
  beta_pre_indx_array_t(nb) % npoint = 0
enddo

if(.not. chi_initialized) call error(sub_name, 'chi not initialized')

do k=1, nz - 1

  do j=1, ny

    do i = 1, nx
    
      ! map brindx to clindx
      brindx_p => brindx(i,j,k)
      
      !write(*,*) 'brindx_p : ', brindx_p
      
      if ( brindx_p > 0 ) then
      
        br_loc_id_p => brindx_to_loc_id(:, brindx_p)
        
        nt_p => rns_tree_iarray( br_loc_id_p(1) ) ! Map unique tree it to wrapped tree in rns domain
        ng_p => br_loc_id_p(2)
        nc_p => br_loc_id_p(3)
      
        if( nt_p < 1) then ! nt < 1
          call error(sub_name, 'brindx(i,j,k) : ', brindx_p)
          call mesg(sub_name, 'coord : ', coord)
          call error(sub_name, 'br_loc_id_p(1) : ', br_loc_id_p(1))
        endif

        !  Setting cluster id it belongs to
        clindx_p => rns_reslv_cl_iarray( tr_t( nt_p ) % gen_t( ng_p ) % cl_t ( ng_p ) % indx )
      
        if( ng_p <= tr_t( nt_p ) % ngen_reslv ) then
          !  Use only inside points
          if ( phi(i,j,k) <= 0._rprec ) then 
          
            !write(*,*) 'setting resolved point'
        
            cl_pre_indx_array_t(clindx_p) % npoint = cl_pre_indx_array_t(clindx_p) % npoint + 1
        
            cl_pre_indx_array_t(clindx_p) % iarray(1, cl_pre_indx_array_t(clindx_p) % npoint) = i
            cl_pre_indx_array_t(clindx_p) % iarray(2, cl_pre_indx_array_t(clindx_p) % npoint) = j
            cl_pre_indx_array_t(clindx_p) % iarray(3, cl_pre_indx_array_t(clindx_p) % npoint) = k
          
          endif
        
        elseif ( ng_p <= tr_t( nt_p ) % ngen ) then

          !write(*,*) 'chi(i,j,k) : ', chi(i,j,k)
          
          if( chi(i,j,k) > chi_cutoff) then
          
            !write(*,*) 'setting unresolved point'
            if( .not. use_beta_sub_regions) then
            
              beta_pre_indx_array_t( nt_p ) % npoint = beta_pre_indx_array_t( nt_p ) % npoint + 1
          
              beta_pre_indx_array_t( nt_p ) % iarray(1, beta_pre_indx_array_t( nt_p ) % npoint) = i
              beta_pre_indx_array_t( nt_p ) % iarray(2, beta_pre_indx_array_t( nt_p ) % npoint) = j
              beta_pre_indx_array_t( nt_p ) % iarray(3, beta_pre_indx_array_t( nt_p ) % npoint) = k 
              
            endif
          
          endif
      
        else
          
          call error(sub_name, 'Total number of generations in tree exceeded.')
        
        endif
      
        nullify(br_loc_id_p, clindx_p, brindx_p)
        nullify(nt_p, ng_p, nc_p)
        
      endif
      
      nullify(brindx_p)
      
    enddo
    
  enddo
  
enddo

!  Allocate true indx_array
allocate(cl_indx_array_t( ncluster_reslv ) )


do nc=1, ncluster_reslv
  allocate(cl_indx_array_t(nc) % iarray(3, cl_pre_indx_array_t(nc) % npoint))
enddo


do nc=1, ncluster_reslv
  
  cl_indx_array_t(nc) % npoint = cl_pre_indx_array_t(nc) % npoint
  
  do np = 1, cl_indx_array_t(nc) % npoint
  
    cl_indx_array_t(nc) % iarray(:,np) = cl_pre_indx_array_t(nc) % iarray(:,np)
    
  enddo
  
enddo

allocate(beta_indx_array_t ( nbeta ) )

do nb = 1, nbeta
  allocate( beta_indx_array_t(nb) % iarray(3, beta_pre_indx_array_t(nb) % npoint))
enddo

do nb = 1, nbeta
  beta_indx_array_t(nb) % npoint = beta_pre_indx_array_t(nb) % npoint
  
  do np = 1, beta_indx_array_t(nb) % npoint
    beta_indx_array_t(nb) % iarray(:,np) = beta_pre_indx_array_t(nb) % iarray(:,np)
  enddo
  
enddo

!if(.not. USE_MPI .or. (USE_MPI .and. coord == 2)) then
!  write(*,*) 'Number of points per cluster : '
!  do nc = 1, ncluster_tot
!    write(*,*) '-------------------------'
!    write(*,*) 'nc, cl_indx_array_t(nc) % npoint : ', nc, cl_indx_array_t(nc) % npoint
!    write(*,*) '-------------------------'
!  enddo
!endif

!  No longer needed
deallocate(cl_pre_indx_array_t)
deallocate(beta_pre_indx_array_t)

!!  Sort each cl_indx_array_t into column major order on the iarray output
!do nc=1, ncluster_tot

!  call isortcm(cl_indx_array_t(nc) % iarray, 3, cl_indx_array_t(nc) % npoint)
!  
!enddo

!if(coord == 0) call mesg(sub_name, 'Exiting ' // sub_name)

return
end subroutine rns_fill_indx_array_ls

!!**********************************************************************
!subroutine rns_set_cl_parent_ls()
!!**********************************************************************
!!  This subroutine sets the parent (a resolved cluster) to each unresolved
!!  cluster which will be used for the CD calculations
!!
!!use types, only : rprec
!!$if($CYL_SKEW_LS)
!!use cyl_skew_base_ls, only : tree, cluster, tr_t, ntree
!!$endif
!!use messages
!implicit none

!character (*), parameter :: sub_name = mod_name // '.rns_set_cl_parent_ls'

!integer :: nt, ng, nc

!type(tree), pointer :: tr_t_p
!type(cluster), pointer :: cl_t_p

!nullify(tr_t_p, cl_t_p)

!!  Allocate and initialize forces
!!if(.not. use_main_tree_only) call error(sub_name,'use_main_tree_only must be true')
!allocate( clforce_t ( ncluster_tot ) ) 
!clforce_t = force(parent = 0, CD = 0._rprec, fD = 0._rprec, kappa=0._rprec)

!do nt = 1, ntree

!  tr_t_p => tr_t(nt) 

!  do ng = 1, tr_t_p % ngen
!  
!    do nc = 1, tr_t_p % gen_t(ng) % ncluster
!    
!      cl_t_p  => tr_t_p % gen_t(ng) % cl_t(nc)
!    
!      if(ng <= tr_t_p % ngen_reslv) then
!       !  is itself
!        clforce_t(cl_t_p % indx) % parent = cl_t_p % indx
!      
!      elseif (ng <= tr_t_p % ngen) then
!    
!        clforce_t(cl_t_p % indx) % parent = clforce_t(cl_t_p % parent) % parent
!      
!      else
!     
!        call error(sub_name, 'Generation number not computed correctly.')
!      
!      endif
!      
!      nullify(cl_t_p)
!      
!    enddo
!  
!  enddo
!  
!  nullify(tr_t_p)
!  
!enddo

!return
!end subroutine rns_set_cl_parent_ls

!**********************************************************************
subroutine rns_CD_ls()
!**********************************************************************
!  This subroutine handles all CD calculation within the RNS module; 
!  all CD and force calculations associated with
!
!  tree -> generation -> cluster -> branch
!
!  are handled here
!
use param, only : jt, USE_MPI, coord
use messages
!!$if($CYL_SKEW_LS)
!!use cyl_skew_base_ls, only : ngen, ngen_reslv
!!$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_CD_ls'

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

!if(clforce_calc) then

  call rns_cl_force_ls()     !  Get CD, force etc for resolved regions
  call rns_beta_force_ls()  !  Get force of 
    
  !if(ngen > ngen_reslv) call rns_cl_unreslv_CD_ls()
    
  if(modulo (jt, clforce_nskip) == 0) then
    
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then

      call rns_write_cl_CD_ls()
      call rns_write_cl_fD_ls()
      call rns_write_cl_vel_ls()
      
      call rns_write_beta_CD_ls()
      call rns_write_beta_fD_ls()
      call rns_write_beta_vel_ls()
      call rns_write_beta_kappa_ls()
      
    endif
    
  endif
  
!endif

!if(coord == 0) call mesg(sub_name, 'Exiting ' // sub_name)

return
end subroutine rns_CD_ls


!**********************************************************************
subroutine rns_cl_force_ls()
!**********************************************************************
!  This subroutine computes the CD of the branch cluster (cl) associated
!  with each region dictated by the brindx value. The cl is mapped from 
!  brindex
!
use types, only : rprec
!!use param, only : nx, ny, nz, dx, dy, dz
!!use param, only : USE_MPI, coord
!!$if($MPI)
!!use param, only : MPI_RPREC, MPI_SUM, rank_of_coord, comm, ierr
!!$endif
!!use sim_param, only : u
!!use functions, only : plane_avg_3D
!!use level_set_base, only : phi
!!use immersedbc, only : fx
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : ntree, tr_t
$endif
use messages
use param, only : nx, ny, nz, dx, dy, dz, coord
$if($MPI)
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif
use sim_param, only : u
use functions, only : plane_avg_3D
use immersedbc, only : fx
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_cl_force_ls'

integer, pointer :: clindx_p
integer, pointer :: npoint_p
integer, pointer, dimension(:,:) :: iarray_p
integer, pointer :: i, j, k
real(rprec), pointer :: fD_p

integer :: ncluster_tot
integer :: nt, ng, nc, np

$if ($MPI)
real(rprec) :: cl_fD
$endif

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

!  Comment starts here 
nullify(clindx_p)
nullify(npoint_p, iarray_p)
nullify(i,j,k)
nullify(fD_p)

!!$if ($MPI)
!!allocate (cl_fD ( ncluster_reslv_ref ) )
!!cl_fD = 0._rprec
!!$endif

do nt = 1, rns_ntree

  do ng = 1, tr_t(nt) % ngen_reslv
  
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster
    
      clindx_p => rns_reslv_cl_iarray(tr_t(nt) % gen_t(ng) % cl_t(nc) % indx)
   
      cl_ref_plane_t(clindx_p) % u = plane_avg_3D( u(1:nx,:,1:nz), cl_ref_plane_t(clindx_p) % p1, cl_ref_plane_t(clindx_p) % p2, &
        cl_ref_plane_t(clindx_p) % p3, cl_ref_plane_t(clindx_p) % nzeta, cl_ref_plane_t(clindx_p) % neta )
     
      npoint_p => cl_indx_array_t(clindx_p) % npoint
      iarray_p => cl_indx_array_t(clindx_p) % iarray
  
      $if($MPI)
      cl_fD = 0._rprec
      $endif
  !clforce_t(clindx_p) % fD = 0._rprec
  
      fD_p => clforce_t( clindx_p ) % fD
      fD_p = 0._rprec
  
      do np=1, npoint_p
  
        i => iarray_p(1,np)
        j => iarray_p(2,np)
        k => iarray_p(3,np)
  
        $if($MPI)
        cl_fD = cl_fD + fx(i,j,k) * dx * dy * dz
        $else
        fD_p = fD_p + fx(i,j,k) * dx * dy * dz
        $endif
    
        nullify(i,j,k)
    
      enddo
  
      $if($MPI)
      call mpi_allreduce (cl_fD, fD_p, 1, MPI_RPREC, MPI_SUM, comm, ierr)
      $endif
  
      clforce_t(clindx_p) % CD = -fD_p / (0.5_rprec * cl_ref_plane_t(clindx_p)%area * (cl_ref_plane_t(clindx_p)%u)**2)
  !if(clforce_t(clindx_p) % CD <= 0._rprec) call mesg(sub_name,'CD < 0')
  
      nullify(clindx_p)
      nullify(npoint_p, iarray_p)
      nullify(fD_p)
 
      
    enddo
    
  enddo
 
enddo

!  Comment ends here

!do nc = 1, ncluster_reslv_ref
!  reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
!  clindx_p => tr_t(reslv_cl_loc_id_p(1)) % gen_t(reslv_cl_loc_id_p(2)) % cl_t(reslv_cl_loc_id_p(3)) % indx
!  clforce_t(clindx_p) % CD = 1._rprec
!  nullify(reslv_cl_loc_id_p, clindx_p)
!      
!enddo

!if(use_main_tree_ref) then
!!  Need to put CD on other resolved clusters (on other trees)
!  do nt = 2, ntree
!    do ng = 1, tr_t(nt) % ngen_reslv
!      do nc = 1, tr_t(nt) % gen_t (ng) % ncluster

!        clindx_p       => tr_t(1) % gen_t(ng) % cl_t(nc) % indx
!        clindx_other_p => tr_t(nt) % gen_t(ng) % cl_t(nc) % indx
!        
!        clforce_t(clindx_other_p) % CD = clforce_t(clindx_p) % CD
!        
!        nullify(clindx_p, clindx_other_p)
!        
!      enddo
!    enddo
!  enddo
!  
!endif

!if(coord == 0) call mesg(sub_name, 'Exiting ' // sub_name)

return
end subroutine rns_cl_force_ls

!**********************************************************************
subroutine rns_beta_force_ls()
!**********************************************************************
!  This subroutine computes the CD of the branch cluster (cl) associated
!  with each region dictated by the brindx value. 
!
use types, only : rprec
use param, only : dx, dy, dz, nx, ny, nz, jt, coord
use messages
use sim_param, only : u
use immersedbc, only : fx
use functions, only : plane_avg_3D
$if($MPI)
use mpi
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : generation, tr_t
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_beta_force_ls'

integer :: nb, np, nt, nc

integer, pointer :: i,j,k
integer, pointer :: npoint_p
integer, pointer :: clindx_p

real(rprec), pointer, dimension(:) :: p1_p, p2_p, p3_p   
integer, pointer :: nzeta_p, neta_p 
real(rprec), pointer :: area_p, u_p
real(rprec), pointer :: kappa_p, CD_p

real(rprec) :: sigma
real(rprec), allocatable, dimension(:) :: fD_dir

real(rprec) :: CD_num, CD_denom, fD_tot, CD, Lint
$if($MPI)
real(rprec) :: Lint_global
$endif

type(generation), pointer :: gen_t_p

$if($MPI)
real(rprec) :: fD
$endif

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

nullify(i,j,k)
nullify(npoint_p)
nullify(clindx_p)

nullify(p1_p, p2_p, p3_p)
nullify(nzeta_p, neta_p)
nullify(area_p, u_p)
nullify(kappa_p, CD_p)

nullify(gen_t_p)

allocate(fD_dir(nbeta))

!  Compute total drag force all unresolved (beta) regions
!  Need more work to have beta as sub regions
do nb = 1, nbeta 
 
    !  Loop over number of points used in beta region
  npoint_p => beta_indx_array_t( nb ) % npoint
  
  $if($MPI)
  fD = 0._rprec
  $endif
  
  beta_force_t(nb) % fD = 0._rprec
  
  do np = 1, npoint_p
  
    i => beta_indx_array_t( nb ) % iarray(1,np)
    j => beta_indx_array_t( nb ) % iarray(2,np)
    k => beta_indx_array_t( nb ) % iarray(3,np)
    
    $if($MPI)
    fD = fD + fx(i,j,k) * dx * dy * dz
    $else    
    beta_force_t(nb) % fD = beta_force_t(nb) % fD + fx(i,j,k) * dx * dy * dz
    $endif
 
    nullify(i,j,k)
    
  enddo
  
  $if($MPI)
  call mpi_allreduce (fD, beta_force_t(nb) % fD, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  $endif
  
  nullify(npoint_p)

  if(beta_force_t(nb) % fD < 0._rprec) then

    fD_dir(nb) = -1._rprec

  else

    fD_dir(nb) = 1._rprec

  endif
  
enddo

if(use_single_beta_CD) then
  
  CD_num = 0._rprec
  CD_denom = 0._rprec

  do nt = 1, rns_ntree
  
    fD_tot = 0._rprec
  
    p1_p    => rbeta_ref_plane_t (nt) % p1
    p2_p    => rbeta_ref_plane_t (nt) % p2
    p3_p    => rbeta_ref_plane_t (nt) % p3
    nzeta_p => rbeta_ref_plane_t (nt) % nzeta
    neta_p  => rbeta_ref_plane_t (nt) % neta
    area_p  => rbeta_ref_plane_t (nt) % area
    u_p     => rbeta_ref_plane_t (nt) % u
    
    u_p = plane_avg_3D(u(1:nx,1:ny,1:nz), p1_p, p2_p, p3_p, nzeta_p, neta_p)
  
    nullify(p1_p, p2_p, p3_p, nzeta_p, neta_p)  
    
    gen_t_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen_reslv )
    
    do nc = 1, gen_t_p % ncluster
      
      clindx_p =>  rns_reslv_cl_iarray( gen_t_p % cl_t(nc) % indx )
      
      fD_tot = fD_tot + clforce_t(clindx_p) % fD
      
      nullify(clindx_p)
      
    enddo
    
    nullify(gen_t_p)
    
    !  Should sum over all beta regions of tree
    fD_tot = fD_tot + beta_force_t(nt) % fD ! Computed in rns_beta_force_ls
    
    CD_num = CD_num + fD_tot * u_p * dabs( u_p ) * area_p
    CD_denom = CD_denom + u_p * u_p * u_p * u_p * area_p**2
    
    nullify(area_p, u_p)
    
  enddo

  !  Compute CD
  CD = -2._rprec * CD_num / CD_denom
  
  if( jt < CD_ramp_nstep ) CD = dble(jt)/dble(CD_ramp_nstep) * CD
  
  !  This CD goes with the regions beta on all trees ! and is the new CD
  do nb = 1, nbeta
  
    beta_force_t( nb ) % CD = CD
    
  enddo
  
  !  Compute kappa
  !  Compute Lint over each region beta
  do nb = 1, nbeta 
  
    p1_p    => beta_ref_plane_t (nb) % p1
    p2_p    => beta_ref_plane_t (nb) % p2
    p3_p    => beta_ref_plane_t (nb) % p3
    nzeta_p => beta_ref_plane_t (nb) % nzeta
    neta_p  => beta_ref_plane_t (nb) % neta
    area_p  => beta_ref_plane_t (nb) % area
    u_p     => beta_ref_plane_t (nb) % u
    
    u_p = plane_avg_3D(u(1:nx,1:ny,1:nz), p1_p, p2_p, p3_p, nzeta_p, neta_p)
  
    nullify(p1_p, p2_p, p3_p, nzeta_p, neta_p)  
 
    !  Loop over number of points used in beta region
    npoint_p => beta_indx_array_t( nb ) % npoint
    
    Lint = 0._rprec
    
    $if($MPI)
    Lint_global = 0._rprec
    $endif
  
    do np = 1, npoint_p
  
      i => beta_indx_array_t( nb ) % iarray(1,np)
      j => beta_indx_array_t( nb ) % iarray(2,np)
      k => beta_indx_array_t( nb ) % iarray(3,np)
    
      Lint = Lint + dabs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
      nullify(i,j,k)
      
    enddo
    
    nullify( npoint_p )
    
    $if($MPI)
    call mpi_allreduce (Lint, Lint_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
    Lint = Lint_global
    $endif
    
    kappa_p => beta_force_t(nb) % kappa
    CD_p    => beta_force_t(nb) % CD
    
    kappa_p = CD_p * dabs ( u_p ) * area_p * u_p / ( 2._rprec * Lint * dx * dy * dz )
    
    if(coord == 0 .and. (modulo (jt, clforce_nskip) == 0)) write(*,'(1a,i,3f12.6)') 'nb, kappa, CD, Lint : ', nb, kappa_p, CD_p, Lint
    
    nullify(kappa_p, CD_p)
    nullify(u_p, area_p)
        
  enddo
   
  
else

  call error(sub_name, 'use_single_beta_CD must be true')
  
endif
  
deallocate(fD_dir)

end subroutine rns_beta_force_ls

!**********************************************************************
subroutine rns_forcing_ls()
!**********************************************************************
!  This subroutine computes the forces on the unresolved branches
!
use types, only : rprec
use sim_param, only : u
use immersedbc, only : fx
$if($MPI)
use mpi
use param, only : MPI_RPREC, up, down, comm, status, ierr, ld, ny, nz, nproc
$endif
use param, only : dx, dy, dz, coord, jt

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_forcing_ls'

integer :: nb, np

integer, pointer :: i, j, k

integer, pointer :: npoint_p

real(rprec), pointer :: kappa_p

nullify(i,j,k)
nullify(npoint_p)
nullify(kappa_p)


do nb = 1, nbeta 
 
    !  Loop over number of points used in beta region
    npoint_p => beta_indx_array_t( nb ) % npoint
    
    kappa_p  => beta_force_t( nb ) % kappa
  
    do np = 1, npoint_p
  
      i => beta_indx_array_t( nb ) % iarray(1,np)
      j => beta_indx_array_t( nb ) % iarray(2,np)
      k => beta_indx_array_t( nb ) % iarray(3,np)
    
      fx(i,j,k) = - kappa_p * dabs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
      nullify(i,j,k)
      
    enddo
    
    nullify( npoint_p, kappa_p )
    
enddo

$if($MPI)
!  Sync data at overlap nodes
if(coord < nproc - 1) then
  call mpi_recv (fx(:,:,nz), ld*ny, MPI_RPREC, up, 1, comm, status, ierr)
endif
if(coord > 0) then
  call mpi_send (fx(:,:,1), ld*ny, MPI_RPREC, down, 1, comm, ierr)
endif
$endif

return

end subroutine rns_forcing_ls


!**********************************************************************
subroutine rns_write_cl_CD_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only :  reslv_clindx_to_loc_id
$endif

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_CD_ls'
character(*), parameter :: fname = path // 'output/rns_cl_CD_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count
integer, pointer, dimension(:) :: reslv_cl_loc_id_p => null()

!  Write cluster force (CD) for all trees + time step

nvar = ncluster_reslv + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, reslv_cl_loc_id_p(1))
    call strcat(var_list, ',')
    call strcat(var_list, reslv_cl_loc_id_p(2))
    call strcat(var_list, ',')
    call strcat(var_list, reslv_cl_loc_id_p(3))
    call strcat(var_list, '</sub>"')
  enddo
  nullify(reslv_cl_loc_id_p)
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, clforce_t(1:nvar-1)%CD /))

return
end subroutine rns_write_cl_CD_ls



!**********************************************************************
subroutine rns_write_cl_vel_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only :  reslv_clindx_to_loc_id
$endif

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_vel_ls'
character(*), parameter :: fname = path // 'output/rns_cl_vel_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count
integer, pointer, dimension(:) :: reslv_cl_loc_id_p => null()

!  Write cluster force (CD) for all trees + time step

nvar = ncluster_reslv + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, reslv_cl_loc_id_p(1))
    call strcat(var_list, ',')
    call strcat(var_list, reslv_cl_loc_id_p(2))
    call strcat(var_list, ',')
    call strcat(var_list, reslv_cl_loc_id_p(3))
    call strcat(var_list, '</sub>"')
  enddo
  nullify(reslv_cl_loc_id_p)
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, cl_ref_plane_t(1:nvar-1)%u /))

return
end subroutine rns_write_cl_vel_ls

!**********************************************************************
subroutine rns_write_beta_vel_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_beta_vel_ls'
character(*), parameter :: fname = path // 'output/rns_beta_vel_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count

!  Write cluster force (CD) for all trees + time step
nvar = nbeta + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, nc)
    call strcat(var_list, ',')
    call strcat(var_list, 1) ! Need ability to specify beta region
    call strcat(var_list, '</sub>"')
  enddo

  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, beta_ref_plane_t(1:nvar-1) % u /))

return
end subroutine rns_write_beta_vel_ls

!**********************************************************************
subroutine rns_write_cl_fD_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only :  reslv_clindx_to_loc_id
$endif

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_fD_ls'
character(*), parameter :: fname = path // 'output/rns_cl_fD_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count
integer, pointer, dimension(:) :: reslv_cl_loc_id_p => null()

!  Write cluster force (CD) for all trees + time step

nvar = ncluster_reslv + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
    !  Create variable list name:
    call strcat(var_list, ',"fD<sub>')
    call strcat(var_list, reslv_cl_loc_id_p(1))
    call strcat(var_list, ',')
    call strcat(var_list, reslv_cl_loc_id_p(2))
    call strcat(var_list, ',')
    call strcat(var_list, reslv_cl_loc_id_p(3))
    call strcat(var_list, '</sub>"')
  enddo
  nullify(reslv_cl_loc_id_p)
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, clforce_t(1:nvar-1)%fD /))

return
end subroutine rns_write_cl_fD_ls

!**********************************************************************
subroutine rns_write_beta_fD_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_beta_fD_ls'
character(*), parameter :: fname = path // 'output/rns_beta_fD_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count

!  Write cluster force (CD) for all trees + time step
nvar = nbeta + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    !  Create variable list name:
    call strcat(var_list, ',"fD<sub>')
    call strcat(var_list, nc)
    call strcat(var_list, ',')
    call strcat(var_list, 1) ! Need ability to specify beta region
    call strcat(var_list, '</sub>"')
  enddo

  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, beta_force_t(1:nvar-1) % fD /))

return
end subroutine rns_write_beta_fD_ls

!**********************************************************************
subroutine rns_write_beta_CD_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_beta_CD_ls'
character(*), parameter :: fname = path // 'output/rns_beta_CD_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count

!  Write cluster force (CD) for all trees + time step
nvar = nbeta + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, nc)
    call strcat(var_list, ',')
    call strcat(var_list, 1) ! Need ability to specify beta region
    call strcat(var_list, '</sub>"')
  enddo

  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, beta_force_t(1:nvar-1) % CD /))

return
end subroutine rns_write_beta_CD_ls

!**********************************************************************
subroutine rns_write_beta_kappa_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_beta_kappa_ls'
character(*), parameter :: fname = path // 'output/rns_beta_kappa_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count

!  Write cluster force (CD) for all trees + time step
nvar = nbeta + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    !  Create variable list name:
    call strcat(var_list, ',"<greek>k</greek><sub>')
    call strcat(var_list, nc)
    call strcat(var_list, ',')
    call strcat(var_list, 1) ! Need ability to specify beta region
    call strcat(var_list, '</sub>"')
  enddo

  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, beta_force_t(1:nvar-1) % kappa /))

return
end subroutine rns_write_beta_kappa_ls

!**********************************************************************
subroutine brindx_init ()
!**********************************************************************
use param, only : iBOGUS, coord, ld, ny, nz
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.brindx_init'
character (*), parameter :: fbrindx_in = 'brindx.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fbrindx_in_MPI
$endif

integer :: ip, i,j,k, brindx_max

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fbrindx_in_MPI, '(a,a,i0)') fbrindx_in, MPI_suffix, coord
    
  inquire (file=fbrindx_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fbrindx_in_MPI)

  open (1, file=fbrindx_in_MPI, action='read', position='rewind',  &
         form='unformatted')
  read (1) brindx
  close (1)

  brindx(:, :, nz) = iBOGUS

$else

  inquire (file=fbrindx_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fbrindx_in)

  open (1, file=fbrindx_in, action='read', position='rewind',  &
         form='unformatted')
  read (1) brindx
  close (1)

$endif

!brindx_max = -100
!do k=0,nz
!  do j=1,ny
!    do i=1,ld
!    
!      if(brindx(i,j,k) > brindx_max) brindx_max = brindx(i,j,k)

!    enddo
!  enddo
!enddo

brindx_max = maxval(brindx)
write(*,*) 'coord, brindx_max : ', coord, brindx_max



brindx_initialized = .true.

end subroutine brindx_init

!**********************************************************************
subroutine chi_init ()
!**********************************************************************
use param, only : iBOGUS, coord
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.chi_init'
character (*), parameter :: fchi_in = 'chi.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fchi_in_MPI
$endif

integer :: ip
real(rprec) :: chi_max

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fchi_in_MPI, '(a,a,i0)') fchi_in, MPI_suffix, coord
    
  inquire (file=fchi_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fchi_in_MPI)

  open (1, file=fchi_in_MPI, action='read', position='rewind',  &
         form='unformatted')
  read (1) chi
  close (1)

  chi(:, :, nz) = iBOGUS

$else

  inquire (file=fchi_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fchi_in)

  open (1, file=fchi_in, action='read', position='rewind',  &
         form='unformatted')
  read (1) chi
  close (1)

$endif

chi_max = maxval(chi)
write(*,*) 'coord, chi_max : ', coord, chi_max

chi_initialized = .true.

end subroutine chi_init

end module rns_ls


