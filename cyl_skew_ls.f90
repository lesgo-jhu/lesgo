!**********************************************************************
module cylinder_skew_ls
!**********************************************************************
use types, only : rprec
use param
use cylinder_skew_base_ls
implicit none

save
private

public :: cylinder_skew_init_ls, cylinder_skew_CD_ls
public :: cylinder_skew_fill_tree_array_ls

$if($RNS_LS)
public :: cylinder_skew_fill_ref_plane_array_ls
$endif

contains

!**********************************************************************
subroutine cylinder_skew_init_ls()
!**********************************************************************

implicit none

character(64) :: fname, temp
integer :: i,j,k,ng

!  Open file which to write global data
fname = path // 'cylinder_skew_gen_ls.out'
$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

!  Read in cylinder_skew_gen.dat file
open (unit = 2,file = fname, status='old',form='formatted', &
  action='read',position='rewind')

allocate(igen(ngen))
allocate(kbottom_inside(ngen))
allocate(kbottom(ngen))
allocate(dz_bottom(ngen))
allocate(ktop_inside(ngen))
allocate(ktop(ngen))
allocate(dz_top(ngen))
allocate(itype(nx+2,ny,nz))

do ng=1,ngen
  read(2,*) igen(ng), kbottom_inside(ng), kbottom(ng), &
    dz_bottom(ng), ktop_inside(ng), ktop(ng), &
    dz_top(ng)
enddo
close(2)

!  Check 1st generation only for need ground association
if(igen(1) /= -1) then
  !  Open file which to write global data
  fname = path // 'cylinder_skew_point_ls.out'
  $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
  $endif

  !  Read in cylinder_skew_gen.dat file
  open (unit = 2,file = fname, status='old',form='formatted', &
    action='read',position='rewind')
  do k=1,nz
    do j = 1,ny
      do i = 1,nx+2
        read(2,*) itype(i,j,k)
      enddo
    enddo
  enddo
   
  close(2)
endif

return
end subroutine cylinder_skew_init_ls

!**********************************************************************
subroutine cylinder_skew_CD_ls ()
!**********************************************************************
use immersedbc, only : fx
use sim_param, only : u
use io, only : jt_total
use messages
use grid_defs, only : zw
implicit none

character (*), parameter :: sub_name = mod_name // '.cylinder_skew_CD_ls'
character (*), parameter :: fCD_out = 'output/cylinder_skew_CD_ls.dat'
character(64) :: fname, temp

integer, parameter :: lun = 991  !--keep open between calls
integer, parameter :: n_calc_CD = 10  !--# t-steps between updates

real (rprec), parameter :: Ap = 1._rprec !--projected area

logical, save, dimension(10) :: file_init=.false. !  May want to change this to allocatable to
                                                  !  match the generation number
logical :: opn, exst

real (rprec) :: CD
real (rprec) :: UAinf, Ainf   !--velocity scale used in calculation of CD
real (rprec) :: fD     !--drag, lift force
real (rprec) :: Uinf_global, UAinf_global, Ainf_global

integer :: i,j,k,n,ng
integer :: kstart, kend

real(rprec) :: dz_start, dz_end
real(rprec) :: dz_p
real(rprec) :: gen_thck

!---------------------------------------------------------------------

if (modulo (jt, n_calc_CD) /= 0) return  !--do nothing

!  Perform an area average for infinite velocity
if(z_bottom_surf <= zw(1)) then
  !  The the global Uinf from the inlet plane; average for proc
  UAinf = sum (u(1, :, 1:nz-1)) * dy * dz
  Ainf  = ny * (nz - 1) * dy * dz
!  Check if z_bottom_surf is associated with proc domain
elseif(zw(1) <= z_bottom_surf .and. z_bottom_surf < zw(nz-1)) then
  UAinf = 0.
  Ainf = 0.
  do k=2,nz-1
    if(zw(k-1) <= z_bottom_surf .and. z_bottom_surf < zw(k)) then
      dz_p = zw(k) - z_bottom_surf
      UAinf = UAinf + sum(u(1,:,k)) * dy * dz_p
      Ainf = Ainf + ny * dy * dz_p     
    elseif(z_bottom_surf <= zw(k-1)) then
      UAinf = UAinf + sum(u(1,:,k)) * dy * dz
      Ainf  = Ainf + ny * dy * dz
    endif
  enddo
endif

$if ($MPI)

  !  Sum Uinf from all procs and redestribute
  call mpi_allreduce(UAinf, UAinf_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce(Ainf, Ainf_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  !  Average over all procs; assuming distribution is even
  Uinf_global = UAinf_global / Ainf_global
  
$else

  Uinf_global = UAinf/Ainf

$endif

do ng=1,ngen
  if(igen(ng) /= -1) then
    
    fD=0.

    !  Check if bottom is in proc domain
    if(kbottom_inside(ng) == 1) then
      kstart = kbottom(ng)
      dz_start = dz_bottom(ng)
    else
      kstart = 1
      dz_start = dz
    endif

    if(ktop_inside(ng) == 1) then
      kend = ktop(ng)
      dz_end = dz_top(ng)
    else
      kend = nz-1 ! -1 to avoid interprocessor overlap
      dz_end = dz
    endif

     !--(-) since want force ON cylinder
     !--dx*dy*dz is since force is per cell (unit volume)
     !--may want to restrict this sum to points with phi < 0.
    do k=kstart,kend
      if(k==kstart) then
        dz_p = dz_start
      elseif(k==kend) then
        dz_p = dz_end
      else
        dz_p = dz
      endif
      if(ng==1) then !  Want to check with ground association
        do j=1,ny
          do i=1,nx
            if(itype(i,j,k) /= 0) fD = fD - fx(i,j,k) * dx * dy * dz_p
          enddo
        enddo
      else
        fD = fD - sum(fx(1:nx, :, k)) * dx * dy * dz_p
      endif
    enddo

    CD = fD / (0.5_rprec * Ap * Uinf_global**2)

    inquire (lun, exist=exst, opened=opn)

    if (.not. exst) call error (sub_name, 'unit', lun, ' nonexistant')
    if (opn) call error (sub_name, 'unit', lun, ' is already open')

    !  Open file which to write global data
    fname = trim(fCD_out)
    write(temp,'(".g",i0)') ng
    fname = trim (fname) // temp

    $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
    $endif

    if (.not. file_init(ng)) then  !--set up file for output
    !  Compute thickness of generation associated with proc
      open (lun, file=fname, position='rewind')

      gen_thck=0._rprec
      do k=kstart,kend
        if(k==kstart) then
          dz_p = dz_start
        elseif(k==kend) then
          dz_p = dz_end
        else
          dz_p = dz
        endif
        gen_thck = gen_thck + dz_p
      enddo
    !--write a header
      write (lun, '(a,es12.5)') '# Ap = ', Ap
      write (lun, '(a,es12.5)') '# Gen thickness = ', gen_thck
      write (lun, '(a)') '# t, CD, fD, Uinf' 

      file_init(ng) = .true.
    else
   
      open (lun, file=fname, position='append')

    end if


    !--output to file
    write (lun, '(4(es12.5,1x))') (jt_total * dt), CD, fD, Uinf_global

    close (lun)  !--only do this to force a flush

  end if

enddo

return
end subroutine cylinder_skew_CD_ls

$if ($DEVEL)
!**********************************************************************
subroutine cylinder_skew_RNS_CD_ls ()
!**********************************************************************
use immersedbc, only : fx
use sim_param, only : u
use io, only : jt_total
use messages
use grid_defs, only : zw
implicit none

character (*), parameter :: sub_name = mod_name // '.cylinder_skew_CD_ls'
character (*), parameter :: fCD_out = 'output/cylinder_skew_CD_ls.dat'
character(64) :: fname, temp

integer, parameter :: lun = 991  !--keep open between calls
integer, parameter :: n_calc_CD = 10  !--# t-steps between updates

real (rprec), parameter :: Ap = 1._rprec !--projected area

logical, save, dimension(10) :: file_init=.false. !  May want to change this to allocatable to
                                                  !  match the generation number
logical :: opn, exst

real (rprec) :: CD
real (rprec) :: UAinf, Ainf   !--velocity scale used in calculation of CD
real (rprec) :: fD     !--drag, lift force
real (rprec) :: Uinf_global, UAinf_global, Ainf_global

integer :: i,j,k,n,ng
integer :: kstart, kend

real(rprec) :: dz_start, dz_end
real(rprec) :: dz_p
real(rprec) :: gen_thck

!---------------------------------------------------------------------

if (modulo (jt, n_calc_CD) /= 0) return  !--do nothing

!  Perform an area average for infinite velocity
if(z_bottom_surf <= zw(1)) then
  !  The the global Uinf from the inlet plane; average for proc
  UAinf = sum (u(1, :, 1:nz-1)) * dy * dz
  Ainf  = ny * (nz - 1) * dy * dz
!  Check if z_bottom_surf is associated with proc domain
elseif(zw(1) <= z_bottom_surf .and. z_bottom_surf < zw(nz-1)) then
  UAinf = 0.
  Ainf = 0.
  do k=2,nz-1
    if(zw(k-1) <= z_bottom_surf .and. z_bottom_surf < zw(k)) then
      dz_p = zw(k) - z_bottom_surf
      UAinf = UAinf + sum(u(1,:,k)) * dy * dz_p
      Ainf = Ainf + ny * dy * dz_p     
    elseif(z_bottom_surf <= zw(k-1)) then
      UAinf = UAinf + sum(u(1,:,k)) * dy * dz
      Ainf  = Ainf + ny * dy * dz
    endif
  enddo
endif

$if ($MPI)

  !  Sum Uinf from all procs and redestribute
  call mpi_allreduce(UAinf, UAinf_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce(Ainf, Ainf_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  !  Average over all procs; assuming distribution is even
  Uinf_global = UAinf_global / Ainf_global
  
$else

  Uinf_global = UAinf/Ainf

$endif

do ng=1,ngen
  if(igen(ng) /= -1) then
    
    fD=0.

    !  Check if bottom is in proc domain
    if(kbottom_inside(ng) == 1) then
      kstart = kbottom(ng)
      dz_start = dz_bottom(ng)
    else
      kstart = 1
      dz_start = dz
    endif

    if(ktop_inside(ng) == 1) then
      kend = ktop(ng)
      dz_end = dz_top(ng)
    else
      kend = nz-1 ! -1 to avoid interprocessor overlap
      dz_end = dz
    endif

     !--(-) since want force ON cylinder
     !--dx*dy*dz is since force is per cell (unit volume)
     !--may want to restrict this sum to points with phi < 0.
    do k=kstart,kend
      if(k==kstart) then
        dz_p = dz_start
      elseif(k==kend) then
        dz_p = dz_end
      else
        dz_p = dz
      endif
      if(ng==1) then !  Want to check with ground association
        do j=1,ny
          do i=1,nx
            if(itype(i,j,k) /= 0) fD = fD - fx(i,j,k) * dx * dy * dz_p
          enddo
        enddo
      else
        fD = fD - sum(fx(1:nx, :, k)) * dx * dy * dz_p
      endif
    enddo

    CD = fD / (0.5_rprec * Ap * Uinf_global**2)

    inquire (lun, exist=exst, opened=opn)

    if (.not. exst) call error (sub_name, 'unit', lun, ' nonexistant')
    if (opn) call error (sub_name, 'unit', lun, ' is already open')

    !  Open file which to write global data
    fname = trim(fCD_out)
    write(temp,'(".g",i0)') ng
    fname = trim (fname) // temp

    $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
    $endif

    if (.not. file_init(ng)) then  !--set up file for output
    !  Compute thickness of generation associated with proc
      open (lun, file=fname, position='rewind')

      gen_thck=0._rprec
      do k=kstart,kend
        if(k==kstart) then
          dz_p = dz_start
        elseif(k==kend) then
          dz_p = dz_end
        else
          dz_p = dz
        endif
        gen_thck = gen_thck + dz_p
      enddo
    !--write a header
      write (lun, '(a,es12.5)') '# Ap = ', Ap
      write (lun, '(a,es12.5)') '# Gen thickness = ', gen_thck
      write (lun, '(a)') '# t, CD, fD, Uinf' 

      file_init(ng) = .true.
    else
   
      open (lun, file=fname, position='append')

    end if


    !--output to file
    write (lun, '(4(es12.5,1x))') (jt_total * dt), CD, fD, Uinf_global

    close (lun)  !--only do this to force a flush

  end if

enddo

return
end subroutine cylinder_skew_RNS_CD_ls
$endif

!**********************************************************************
subroutine cylinder_skew_fill_tree_array_ls()
!**********************************************************************
!
!  This subroutine sets all values for the tree struct - tr_t; it also
!  defines the key arrays clindex_to_loc_id and brindx_to_loc_id 
!  which are used to get the local id of a global index. This subroutine
!  should be called any time the tree struct and its settings as defined
!  in cylinder_skew_base_ls are needed (only once though)
!
implicit none

integer :: nt, ng, nc, nb, nc_g1
integer :: clindx, brindx

integer, pointer, dimension(:) :: clindx_loc_id

real(rprec) :: angle
real(rprec) :: gen_scale_fact

write(*,*) 'ntree : ', ntree

allocate(tr_t(ntree))

!  Set the number of generations in the tree
tr_t%ngen = ngen
tr_t%ngen_reslv = ngen_reslv

!  Allocate the number of clusters in the generation
do nt=1,ntree
  
  allocate(tr_t(nt)%gen_t( tr_t(nt)%ngen ))

  do ng=1, tr_t(nt)%ngen
    
    !  Set the number of clusters for the generation
    tr_t(nt)%gen_t(ng)%ncluster = nbranch**(ng - 1)
    
    !  Allocate space for the clusters    
    allocate( tr_t(nt)%gen_t(ng)%cl_t( tr_t(nt)%gen_t(ng)%ncluster ))
    
    do nc=1, tr_t(nt)%gen_t(ng)%ncluster
    
        !  Set the number of branches for the cluster - here they are all the same
        tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch = nbranch
        !  Allocate space for the branches
        allocate( tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t( tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch ))
     
    enddo
    
  enddo
 
enddo

!  ----- Start setting values for the tree struc -----

!  Must allocate all arrays before setting data; may be that
!  memory addresses move as needed during the allocation of 
!  sub types
do nt=1, ntree
  call set_tree_origin(nt,tr_t(nt)%origin)
enddo

do nt = 1, ntree
    
  do ng = 1, tr_t(nt)%ngen
    
    gen_scale_fact = scale_fact**(ng-1)
          
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster

      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % offset = gen_scale_fact*offset
           
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % d = gen_scale_fact*d
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % l = gen_scale_fact*l
            
      ! Ellipse minor axis
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % b = &
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % d / 2._rprec 
            
      ! Ellipse major axis
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % a = &
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % b/dcos(skew_angle) 
            
      do nb = 1, tr_t(nt) % gen_t(ng) % cl_t(nc) % nbranch

        angle =  zrot_angle + &
          2.*pi*(nb-1)/(tr_t(nt) % gen_t(ng) % cl_t(nc) % nbranch) + &
          (ng - 1)*pi ! Rotate 180 degrees for each generation
                    
        tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % angle = angle
                   
                
        tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % skew_axis = &
          (/ dcos(angle +pi/2.), dsin(angle + pi/2.), 0._rprec/)
                        
        tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % skew_angle = skew_angle         
                  
                 
                 
      enddo
        
    enddo
        
  enddo
    
enddo

!  Get a global count of clusters and branches (all may not be the same)
clindx = 0
brindx = 0
do nt = 1, ntree
  do ng = 1, tr_t(nt)%ngen
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
      clindx = clindx + 1
      do nb = 1, tr_t(1)%gen_t(ng)%cl_t(nc)%nbranch
        brindx = brindx + 1
      enddo
    enddo
  enddo
enddo

!  Set the total number of clusters and branches of the tree
tr_t = tree ( ncluster = clindx, nbranch = brindx )

allocate(clindx_to_loc_id(3,clindx))
allocate(brindx_to_loc_id(4,brindx))

!  Initialize global indexes for clusters and branches
clindx = 0
brindx = 0
do nt = 1, ntree

    do ng=1, tr_t(nt)%ngen
 
        !  Set cluster id for ng+1 generation
        nc_g1 = 0
        
        do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
        
            !  Update global cluster index
            clindx = clindx + 1  
            tr_t(nt)%gen_t(ng)%cl_t(nc)%indx = clindx
            clindx_to_loc_id(:,clindx) = (/ nt, ng, nc /)
            nullify(clindx_loc_id)
            clindx_loc_id => clindx_to_loc_id(:,clindx)
            
            !  Set cluster origin to tree origin
            if( ng == 1 ) tr_t(nt) % gen_t(ng) % cl_t(nc) % origin = tr_t(nt) % origin
            
            do nb = 1, tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
                !  Update global branch index
                brindx = brindx + 1
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % indx = brindx
                brindx_to_loc_id(:,brindx) = (/ nt, ng, nc, nb /)
                
                !  Update cluster id for ng+1 generation
                nc_g1 = nc_g1 + 1
                    
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % origin
                
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot(1) = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % origin(1) + &          
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % offset * &
                    dcos(tr_t(nt) % gen_t(ng) %cl_t(nc) % br_t(nb) % angle)
                
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot(2) = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % origin(2) + &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % offset * &
                    dsin(tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % angle)
                   
                call rotation_axis_vector_3d( tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % skew_axis, &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % skew_angle, &
                    (/ 0._rprec, 0._rprec, tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % l /), &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % top )
                
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % top = & 
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % top + &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot
                
            !  Now set the cluster origin of the ng+1 cluster (with nc = nb)
                if ( ng < tr_t(nt)%ngen ) then
                    tr_t(nt) % gen_t(ng+1) % cl_t(nc_g1) % origin = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % top
                endif

 
            enddo
                   
        enddo
        
        !  Set the top and bottom plane of the generation - this assumes that all 
        !  branches are the same height!      
        tr_t(nt) % gen_t(ng) % bplane = tr_t(nt) % gen_t(ng) % cl_t(1) % br_t(1) % bot(3)
        tr_t(nt) % gen_t(ng) % tplane = tr_t(nt) % gen_t(ng) % cl_t(1) % br_t(1) % top(3)
    
    enddo
    
enddo

return
end subroutine cylinder_skew_fill_tree_array_ls

$if ($RNS_LS)
!**********************************************************************
subroutine cylinder_skew_fill_ref_plane_array_ls()
!**********************************************************************

use rns_base_ls, only : ref_plane_t
use param, only : dy, dz

implicit none

real(rprec), parameter :: alpha=1._rprec

integer :: nt, ng, nc, nb

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)

integer, pointer :: clindx_p, nbranch_p
real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer, dimension(:) :: origin_p

nullify(d_p, l_p, skew_angle_p, clindx_p)

!  allocate ref_array_t
allocate(ref_plane_t( size(clindx_to_loc_id,2) ))

do nt=1, ntree

  do ng=1, tr_t(nt)%ngen
  
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
    
      nbranch_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
      
      clindx_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%indx
           
      do nb = 1, nbranch

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

      ref_plane_t(clindx_p) % area = area_proj
      !  These are defined to be x - planes (no not the NASA experimental planes)
      ref_plane_t(clindx_p) % nzeta = ceiling( w / dy + 1)
      ref_plane_t(clindx_p) % neta  = ceiling( h_m / dz + 1)
      
      origin_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%origin
      
      !  Offset in the upstream x-direction
      zeta_c = origin_p + (/ -alpha * w, 0._rprec, 0._rprec /)
      
      ref_plane_t(clindx_p) % p1    = zeta_c 
      ref_plane_t(clindx_p) % p1(2) = ref_plane_t(clindx_p) % p1(2) + w / 2._rprec
      
      ref_plane_t(clindx_p) % p2    = ref_plane_t(clindx_p) % p1
      ref_plane_t(clindx_p) % p2(2) = ref_plane_t(clindx_p) % p2(2) - w
      
      ref_plane_t(clindx_p) % p3    = ref_plane_t(clindx_p) % p2
      ref_plane_t(clindx_p) % p3(3) = ref_plane_t(clindx_p) % p3(3) + h_m
      
      nullify(clindx_p, origin_p)
      
    enddo
    
  enddo
 
enddo
      
return

end subroutine cylinder_skew_fill_ref_plane_array_ls
$endif

!**********************************************************************
subroutine cylinder_skew_branch_id_ls( indx_based , indx,  branch_id)
!**********************************************************************
implicit none

logical, intent(in) :: indx_based
integer, intent(in) :: indx
integer, dimension

return
end subroutine cylinder_skew_branch_id_ls

end module cylinder_skew_ls
