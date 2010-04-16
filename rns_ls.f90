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
public :: rns_forcing_ls

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

type(tree), pointer :: tr_t_p

nullify(tr_t_p)

! Load brindx
call brindx_init()
! Load filtered indicator function (chi)
call chi_init()

$if($CYL_SKEW_LS)
if(coord == 0) call mesg ( sub_name, 'filling tree array' )
call cyl_skew_fill_tree_array_ls()
if(coord == 0) call mesg ( sub_name, 'tree array filled' )
$endif

!  Get the total number of clusters 
ncluster_tot = 0
do nt = 1, ntree
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      ncluster_tot = ncluster_tot + 1   
    enddo    
  enddo
enddo

if(use_main_tree_only) then
  ncluster_ref = tr_t(1) % ncluster
  ntree_ref = 1
else
  ncluster_ref = ncluster_tot
  ntree_ref = ntree
endif

ncluster_reslv = 0
do nt = 1, ntree
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      if(ng <= tr_t(nt) % ngen_reslv) ncluster_reslv = ncluster_reslv + 1   
    enddo    
  enddo
enddo

ncluster_unreslv = ncluster_tot - ncluster_reslv

ncluster_reslv_ref = 0
do nt = 1, ntree_ref
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      if(ng <= tr_t(nt) % ngen_reslv) ncluster_reslv_ref = ncluster_reslv_ref + 1   
    enddo    
  enddo
enddo

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
  call rns_fill_cl_ref_plane_array_ls()
  !  Create cluster index array
  call rns_fill_cl_indx_array_ls()
  !  Set parent used for assigning CD's on unresolved branches
  if(coord == 0) write(*,*) ' calling rns_set_cl_parent_ls'
  call rns_set_cl_parent_ls()
  if(coord == 0) write(*,*) 'called rns_set_cl_parent_ls'

endif

write(*,*) 'exiting rns_init_ls'
 
return
end subroutine rns_init_ls

!**********************************************************************
subroutine rns_fill_cl_ref_plane_array_ls()
!**********************************************************************
use types, only : rprec
use param, only : dy, dz, USE_MPI, coord
use messages
$if($CYL_SKEW_LS) 
use cyl_skew_base_ls, only : tr_t
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_fill_cl_ref_plane_array_ls'
real(rprec), parameter :: alpha=1._rprec

integer :: nt, ng, nc, nb

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)

integer, pointer :: clindx_p, nbranch_p
real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer, dimension(:) :: origin_p

nullify(d_p, l_p, skew_angle_p, clindx_p)

allocate(cl_ref_plane_t( ncluster_ref ))

call mesg(sub_name, 'ncluster_ref - tr_t(1) % ncluster : ', ncluster_ref - tr_t(1) % ncluster)

if(ntree_ref < 1) call error( sub_name, 'ntree_ref not specified correctly')

do nt=1, ntree_ref

  do ng=1, tr_t(nt)%ngen
  
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
    
      nbranch_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
      
      clindx_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%indx
      
      if(clindx_p > ncluster_ref) then
        call mesg(sub_name, 'clindx_p : ', clindx_p)
        call mesg(sub_name, 'ncluster_ref : ', ncluster_ref)
        call error(sub_name, 'clindx_p > ncluster_ref')
      endif
      
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

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) 'Reference Plane Values for Tree 1 : '
  nt=1
    do ng = 1, tr_t(nt)%ngen
      do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
        write(*,*) '-------------------------'
        write(*,*) 'nt, ng, nc : ', nt, ng, nc
        write(*,*) 'nzeta, neta : ', cl_ref_plane_t(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx) % nzeta, &
          cl_ref_plane_t(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx) % neta
        write(*,*) 'p1 : ', cl_ref_plane_t(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx) % p1
        write(*,*) 'p2 : ', cl_ref_plane_t(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx) % p2
        write(*,*) 'p3 : ', cl_ref_plane_t(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx) % p3
        write(*,*) 'area : ', cl_ref_plane_t(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx) % area
        write(*,*) '-------------------------'
      enddo
    enddo

endif
          
      
return

end subroutine rns_fill_cl_ref_plane_array_ls

!**********************************************************************
subroutine rns_fill_cl_indx_array_ls()
!**********************************************************************
!  This subroutine sets the indx_array for both resolved and unresolved 
!  branches
!
use types, only : rprec
use param, only : nx,ny,nz, coord
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : ngen, ngen_reslv, brindx_to_loc_id, tr_t
$endif
use level_set_base, only : phi
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_fill_cl_indx_array_ls'

integer :: i,j,k, nc, np
integer, pointer :: clindx_p, ng_p, brindx_p
integer, pointer, dimension(:) :: br_loc_id_p
type(indx_array), pointer, dimension(:) :: cl_pre_indx_array

nullify(ng_p, clindx_p, br_loc_id_p, cl_pre_indx_array)

allocate(cl_pre_indx_array( ncluster_tot ) )

do nc=1, ncluster_tot !  same for all clusters
  allocate(cl_pre_indx_array(nc) % iarray(3,nx*ny*(nz-1)))
enddo

!  Intialize the number of points assigned to the cluster
cl_pre_indx_array % npoint = 0

if(.not. chi_initialized) call error(sub_name, 'chi not initialized')

do k=1, nz - 1

  do j=1, ny

    do i = 1, nx
    
      ! map brindx to clindx
      brindx_p => brindx(i,j,k)
      
      !write(*,*) 'brindx_p : ', brindx_p
      
      if ( brindx_p > 0 ) then
      
        br_loc_id_p => brindx_to_loc_id(:, brindx_p)
      
        if(br_loc_id_p(1) < 1) then
          call error(sub_name, 'brindx(i,j,k) : ', brindx_p)
          call mesg(sub_name, 'coord : ', coord)
          call error(sub_name, 'br_loc_id_p(1) : ', br_loc_id_p(1))
        endif
      
        clindx_p => tr_t(br_loc_id_p(1)) % gen_t(br_loc_id_p(2)) % cl_t (br_loc_id_p(3)) % indx
 
        ng_p => br_loc_id_p(2)
        
        !write(*,*) 'coord, br_loc_id_p(1), br_loc_id_p(2), br_loc_id_p(3) : ', coord, br_loc_id_p(1), br_loc_id_p(2), br_loc_id_p(3)
      
        if( ng_p <= ngen_reslv ) then
      
          if ( phi(i,j,k) <= 0._rprec ) then 
          
            write(*,*) 'setting resolved point'
        
            cl_pre_indx_array(clindx_p) % npoint = cl_pre_indx_array(clindx_p) % npoint + 1
        
            cl_pre_indx_array(clindx_p) % iarray(1, cl_pre_indx_array(clindx_p) % npoint) = i
            cl_pre_indx_array(clindx_p) % iarray(2, cl_pre_indx_array(clindx_p) % npoint) = j
            cl_pre_indx_array(clindx_p) % iarray(3, cl_pre_indx_array(clindx_p) % npoint) = k
          
          endif
        
        elseif ( ng_p <= ngen ) then
        
          write(*,*) 'setting unresolved point'
                  
          !write(*,*) 'chi(i,j,k) : ', chi(i,j,k)
          
          if( chi(i,j,k) > chi_cutoff) then
        
            cl_pre_indx_array(clindx_p) % npoint = cl_pre_indx_array(clindx_p) % npoint + 1
          
            cl_pre_indx_array(clindx_p) % iarray(1, cl_pre_indx_array(clindx_p) % npoint) = i
            cl_pre_indx_array(clindx_p) % iarray(2, cl_pre_indx_array(clindx_p) % npoint) = j
            cl_pre_indx_array(clindx_p) % iarray(3, cl_pre_indx_array(clindx_p) % npoint) = k          
          
          endif
      
        else
      
          call error(sub_name, 'Generation number not computed correctly.')
        
        endif
      
        nullify(br_loc_id_p, clindx_p, ng_p)
        
      endif
      
      nullify(brindx_p)
      
    enddo
    
  enddo
  
enddo

!  Allocate true indx_array
allocate(cl_indx_array(ncluster_tot))

do nc=1, ncluster_tot
  allocate(cl_indx_array(nc) % iarray(3, cl_pre_indx_array(nc) % npoint))
enddo

do nc=1, ncluster_tot
  
  cl_indx_array(nc) % npoint = cl_pre_indx_array(nc) % npoint
  
  do np = 1, cl_indx_array(nc) % npoint
  
    cl_indx_array(nc) % iarray(:,np) = cl_pre_indx_array(nc) % iarray(:,np)
    
  enddo
  
enddo

!do nc = 1, ncluster_tot

!  write(*,*) 'cl_indx_array(nc) % npoint : ', cl_indx_array(nc) % npoint
!enddo

!  No longer needed
deallocate(cl_pre_indx_array)

!!  Sort each cl_indx_array into column major order on the iarray output
!do nc=1, ncluster_tot

!  call isortcm(cl_indx_array(nc) % iarray, 3, cl_indx_array(nc) % npoint)
!  
!enddo

return
end subroutine rns_fill_cl_indx_array_ls

!**********************************************************************
subroutine rns_set_cl_parent_ls()
!**********************************************************************
!  This subroutine sets the parent (a resolved cluster) to each unresolved
!  cluster which will be used for the CD calculations
!
use types, only : rprec
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : tree, cluster, tr_t, ntree
$endif
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_set_cl_parent_ls'

integer :: nt, ng, nc

type(tree), pointer :: tr_t_p
type(cluster), pointer :: cl_t_p

nullify(tr_t_p, cl_t_p)

!if(.not. use_main_tree_only) call error(sub_name,'use_main_tree_only must be true')
allocate( clforce_t ( ncluster_tot ) ) 

clforce_t = force(parent = 0, CD = 0._rprec, fD = 0._rprec)

do nt = 1, ntree

  tr_t_p => tr_t(nt) 

  do ng = 1, tr_t_p % ngen
  
    do nc = 1, tr_t_p % gen_t(ng) % ncluster
    
      cl_t_p   => tr_t_p % gen_t(ng) % cl_t(nc)
    
      if(ng <= tr_t_p % ngen_reslv) then
       !  is itself
        clforce_t(cl_t_p % indx) % parent = cl_t_p % indx
      
      elseif (ng <= tr_t_p % ngen) then
    
        clforce_t(cl_t_p % indx) % parent = clforce_t(cl_t_p % parent) % parent
      
      else
     
        call error(sub_name, 'Generation number not computed correctly.')
      
      endif
      
      nullify(cl_t_p)
      
    enddo
  
  enddo
  
  nullify(tr_t_p)
  
enddo

return
end subroutine rns_set_cl_parent_ls

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
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : ngen, ngen_reslv
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'

if(clforce_calc) then

  call rns_cl_reslv_CD_ls()
    
  if(ngen > ngen_reslv) call rns_cl_unreslv_CD_ls()
    
  if(modulo (jt, clforce_nskip) == 0) then
    
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then
      call rns_write_cl_CD_ls()
      if(clforce_vel_write) call rns_write_cl_vel_ls()
    endif
    
  endif
  
endif

return
end subroutine rns_CD_ls


!**********************************************************************
subroutine rns_cl_reslv_CD_ls()
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
use cyl_skew_base_ls, only : reslv_clindx_to_loc_id, ntree, tr_t
$endif
use messages
use param, only : nx, ny, nz, dx, dy, dz
$if($MPI)
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif
use sim_param, only : u
use functions, only : plane_avg_3D
use immersedbc, only : fx
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_cl_reslv_CD_ls'

integer, pointer, dimension(:) :: reslv_cl_loc_id_p 
integer, pointer :: clindx_p, clindx_other_p
integer, pointer :: npoint_p
integer, pointer, dimension(:,:) :: iarray_p
integer, pointer :: i, j, k

integer :: ncluster_tot
integer :: nt, ng, nc, np
$if ($MPI)
real(rprec), pointer, dimension(:) :: cl_fD
$endif

nullify(reslv_cl_loc_id_p, clindx_p)
nullify(clindx_other_p)
nullify(npoint_p, iarray_p)
nullify(i,j,k)

$if ($MPI)
allocate (cl_fD ( ncluster_reslv_ref ) )
cl_fD = 0._rprec
$endif

!  Get reference velocity for only the resolved reference planes reference planes
do nc = 1, ncluster_reslv_ref

  reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
  clindx_p => tr_t(reslv_cl_loc_id_p(1)) % gen_t(reslv_cl_loc_id_p(2)) % cl_t(reslv_cl_loc_id_p(3)) % indx
  
  cl_ref_plane_t(clindx_p) % u = plane_avg_3D( u(1:nx,:,1:nz), cl_ref_plane_t(clindx_p) % p1, cl_ref_plane_t(nc) % p2, &
    cl_ref_plane_t(clindx_p) % p3, cl_ref_plane_t(clindx_p) % nzeta, cl_ref_plane_t(clindx_p) % neta )
    
  nullify(reslv_cl_loc_id_p, clindx_p)
enddo

clforce_t % fD = 0._rprec

do nc = 1, ncluster_reslv_ref

  reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
  clindx_p => tr_t(reslv_cl_loc_id_p(1)) % gen_t(reslv_cl_loc_id_p(2)) % cl_t(reslv_cl_loc_id_p(3)) % indx
  
  npoint_p => cl_indx_array(clindx_p) % npoint
  iarray_p => cl_indx_array(clindx_p) % iarray
  
  do np=1, npoint_p
  
    i => iarray_p(1,np)
    j => iarray_p(2,np)
    k => iarray_p(3,np)
  
    $if($MPI)
    cl_fD(clindx_p) = cl_fD(clindx_p) - fx(i,j,k) * dx * dy * dz
    $else
    clforce_t(clindx_p)%fD = clforce_t(clindx_p)%fD - fx(i,j,k) * dx * dy * dz
    $endif
    
    nullify(i,j,k)
    
  enddo
  
  nullify(reslv_cl_loc_id_p, clindx_p)
  nullify(npoint_p, iarray_p)
  
enddo
 
$if($MPI)
!  Need to sum forces over all processors
do nc=1, ncluster_reslv_ref

  reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
  clindx_p => tr_t(reslv_cl_loc_id_p(1)) % gen_t(reslv_cl_loc_id_p(2)) % cl_t(reslv_cl_loc_id_p(3)) % indx
  
  call mpi_allreduce (cl_fD(clindx_p), clforce_t(clindx_p)%fD, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  
  nullify(reslv_cl_loc_id_p, clindx_p)
enddo

deallocate(cl_fD)
$endif 

do nc = 1, ncluster_reslv_ref

  reslv_cl_loc_id_p => reslv_clindx_to_loc_id(:,nc)
  clindx_p => tr_t(reslv_cl_loc_id_p(1)) % gen_t(reslv_cl_loc_id_p(2)) % cl_t(reslv_cl_loc_id_p(3)) % indx

  clforce_t(clindx_p) % CD = clforce_t(clindx_p)%fD / (0.5_rprec * cl_ref_plane_t(clindx_p)%area * (cl_ref_plane_t(clindx_p)%u)**2)
  
  nullify(reslv_cl_loc_id_p, clindx_p)
      
enddo

if(use_main_tree_only) then
!  Need to put CD on other resolved clusters (on other trees)
  do nt = 2, ntree
    do ng = 1, tr_t(nt) % ngen_reslv
      do nc = 1, tr_t(nt) % gen_t (ng) % ncluster

        clindx_p       => tr_t(1) % gen_t(ng) % cl_t(nc) % indx
        clindx_other_p => tr_t(nt) % gen_t(ng) % cl_t(nc) % indx
        
        clforce_t(clindx_other_p) % CD = clforce_t(clindx_p) % CD
        
        nullify(clindx_p, clindx_other_p)
        
      enddo
    enddo
  enddo
  
endif
        

return
end subroutine rns_cl_reslv_CD_ls

!**********************************************************************
subroutine rns_cl_unreslv_CD_ls()
!**********************************************************************
!  This subroutine computes the CD of the branch cluster (cl) associated
!  with each region dictated by the brindx value. 
!
use types, only : rprec
use param, only : coord
use messages
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : unreslv_clindx_to_loc_id, tr_t
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_cl_unreslv_CD_ls'

integer, pointer :: clindx_p
integer, pointer, dimension(:) :: unreslv_cl_loc_id_p
integer, pointer :: parent_p

integer :: nc

nullify(unreslv_cl_loc_id_p, clindx_p, parent_p)

!if(coord == 0) write(*,*) 'called rns_cl_unreslv_CD_ls '

!if(coord == 0) write(*,*) 'ncluster_unreslv : ', ncluster_unreslv

!  All unresolved branches get a CD 
do nc = 1, ncluster_unreslv

  unreslv_cl_loc_id_p => unreslv_clindx_to_loc_id(:, nc)

  clindx_p            => tr_t(unreslv_cl_loc_id_p(1)) % gen_t(unreslv_cl_loc_id_p(2)) % cl_t(unreslv_cl_loc_id_p(3)) % indx
  
 ! if(coord == 0) call mesg(sub_name, 'clindx_p ', clindx_p)
  
  parent_p => clforce_t(clindx_p) % parent
  
 ! if(coord == 0) call mesg(sub_name, 'parent_p ', parent_p)

  clforce_t(clindx_p) % CD = clforce_t( parent_p ) % CD
  
  nullify(unreslv_cl_loc_id_p, clindx_p)
  nullify(parent_p)
  
enddo

return
end subroutine rns_cl_unreslv_CD_ls

!**********************************************************************
subroutine rns_forcing_ls()
!**********************************************************************
!  This subroutine computes the forces on the unresolved branches
!
use types, only : rprec
use sim_param, only : u
use immersedbc, only : fx
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : ngen, ngen_reslv, tr_t, unreslv_clindx_to_loc_id
$endif
use param, only : dx, dy, dz, coord

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_forcing_ls'

integer :: nc, np

integer, pointer :: i, j, k

integer, pointer :: clindx_p, npoint_p
integer, pointer, dimension(:) :: unreslv_cl_loc_id_p

nullify(i,j,k)
nullify(unreslv_cl_loc_id_p, clindx_p, npoint_p)



!  Compute force due to unresolved clusters
do nc = 1, ncluster_unreslv

  !  Get the global cluster index
  unreslv_cl_loc_id_p => unreslv_clindx_to_loc_id(:, nc)
  clindx_p => tr_t(unreslv_cl_loc_id_p(1)) % gen_t(unreslv_cl_loc_id_p(2)) % cl_t(unreslv_cl_loc_id_p(3)) % indx
  
  !  Loop over number of points used in cluster calc
  npoint_p => cl_indx_array( clindx_p ) % npoint
  
  !write(*,*) 'coord, npoint_p : ', coord, npoint_p
  
  do np = 1, npoint_p
  
    !write(*,*) 'coord, nc, np : ', coord, nc, np
  
    i => cl_indx_array( clindx_p ) % iarray(1,np)
    j => cl_indx_array( clindx_p ) % iarray(2,np)
    k => cl_indx_array( clindx_p ) % iarray(3,np)
    
    fx(i,j,k) = - 0.5_rprec * clforce_t( clindx_p ) % CD * abs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) / (dx*dy*dz)
 
    nullify(i,j,k)
    
  enddo
  
  nullify(unreslv_cl_loc_id_p, clindx_p, npoint_p)
  
enddo

return

end subroutine rns_forcing_ls


!**********************************************************************
subroutine rns_write_cl_CD_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only :  clindx_to_loc_id
$endif

implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_CD_ls'
character(*), parameter :: fname = path // 'output/rns_cl_CD_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count
integer, pointer, dimension(:) :: cl_loc_id_p => null()

real(rprec), pointer, dimension(:) :: cl_CD_p

!  Write cluster force (CD) for all trees + time step
nvar = ncluster_ref + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    cl_loc_id_p => clindx_to_loc_id(:,nc)
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, cl_loc_id_p(1))
    call strcat(var_list, ',')
    call strcat(var_list, cl_loc_id_p(2))
    call strcat(var_list, ',')
    call strcat(var_list, cl_loc_id_p(3))
    call strcat(var_list, '</sub>"')
  enddo
  nullify(cl_loc_id_p)
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
use cyl_skew_base_ls, only :  clindx_to_loc_id
$endif
implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_vel_ls'
character(*), parameter :: fname = path // 'output/rns_cl_vel_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count
integer, pointer, dimension(:) :: cl_loc_id_p => null()

real(rprec), pointer, dimension(:) :: cl_CD_p

!  Write cluster velocity for all trees + time step
nvar = ncluster_ref + 1

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do nc = 1, nvar-1
  
    cl_loc_id_p => clindx_to_loc_id(:,nc)
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, cl_loc_id_p(1))
    call strcat(var_list, ',')
    call strcat(var_list, cl_loc_id_p(2))
    call strcat(var_list, ',')
    call strcat(var_list, cl_loc_id_p(3))
    call strcat(var_list, '</sub>"')
  enddo
  nullify(cl_loc_id_p)
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', nvar, (/ jt_total*dt, cl_ref_plane_t(1:nvar-1)%u /))

return
end subroutine rns_write_cl_vel_ls

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


