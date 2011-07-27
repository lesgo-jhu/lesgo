$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif


!**********************************************************************
module cyl_skew_pre_base_ls
!**********************************************************************
use types, only : rprec, vec3d
use param, only : pi, BOGUS
use cyl_skew_base_ls
use cyl_skew_ls, only : fill_tree_array_ls

implicit none

save 

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

!  Defined local processor definitions
$if($MPI)
integer :: nx_proc
integer :: nproc_csp, global_rank_csp
integer :: stride
$endif

!  cs{0,1} all correspond to vectors with the origin at the
!  corresponding coordinate system
type(cs0), target, allocatable, dimension(:,:,:) :: gcs_t
type(cs1) :: lcs_t, slcs_t, sgcs_t, ecs_t

!  vectors do not have starting point a origin of corresponding
!  coordinate system
type(vec3d) :: vgcs_t

logical :: DIST_CALC=.true.

real(rprec), parameter :: eps = 1.e-12
real(rprec), parameter, dimension(3) :: zrot_axis = (/0.,0.,1./)

logical :: in_cir, in_cyl
logical :: in_cyl_top, in_cyl_bottom
logical :: above_cyl, below_cyl
logical :: in_bottom_surf, btw_planes

integer, dimension(3) :: cyl_loc

end module cyl_skew_pre_base_ls

!**************************************************************
program cyl_skew_pre_ls
!***************************************************************
use cyl_skew_pre_base_ls, only : DIST_CALC, ntree
use cyl_skew_base_ls, only : ngen, ngen_reslv, filter_chi
use messages
implicit none

character (*), parameter :: prog_name = 'cyl_skew_pre_ls'
integer :: nt

if(ngen_reslv > ngen) call error(prog_name, ' ngen_reslv > ngen ')    

call initialize()
!  Loop over all trees
do nt = 1, ntree
  
  $if ($DEBUG)
  if(coord == 0) write(*,*) 'Tree id : ', nt
  $endif
    
  if(DIST_CALC) call main_loop(nt)
  
enddo 

!  Uses global tree info from ebgcs_t and etgcs_t to

if( filter_chi ) call compute_chi()

call finalize()

write(*,*) 'Program completed successfully.'
stop

end program cyl_skew_pre_ls

!**********************************************************************
subroutine initialize()
!**********************************************************************
use param, only : nx, nz_tot, BOGUS
use cyl_skew_pre_base_ls, only : gcs_t, ntree
$if($MPI)
use mpi
use param, only : ierr
use cyl_skew_pre_base_ls, only : global_rank_csp, nproc_csp, stride, nx_proc
$endif
use cyl_skew_base_ls, only : use_bottom_surf, z_bottom_surf, ngen, tr_t
use cyl_skew_ls, only : fill_tree_array_ls

implicit none

integer :: ng, k

$if($MPI)
integer :: nx_proc_sum
integer :: nx_remain, nx_extra
$endif

$if ($MPI)
!call initialize_mpi ()
call initialize_mpi_csp ()

! Load balancing: any left over get consumed in rank order, one-by-one
nx_remain = modulo( nx , nproc_csp )

nx_extra = 0
! Spread the wealth around
if( global_rank_csp < nx_remain )  nx_extra = 1

! Set nx for the given processor
nx_proc = nx / nproc_csp + nx_extra

! Set stride to serve as offset for x-indexing; includes accumulated
! extras from load balacing
stride = (nx / nproc_csp)*global_rank_csp + min( global_rank_csp, nx_remain )

! Write info to screen
write(*,*) 'ID, nx_proc, stride : ', global_rank_csp, nx_proc, stride

! Sanity check for loadbalancing
call mpi_allreduce(nx_proc, nx_proc_sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
if(global_rank_csp == 0 .and. nx_proc_sum /= nx) then
  write(*,*) 'Error in x decomposition - nx_proc_sum, nx : ', nx_proc_sum, nx
  stop
endif
  

$endif

call allocate_arrays()
if( ntree > 0 ) call fill_tree_array_ls()
call generate_grid()

!  Initialize the distance function
gcs_t(:,:,:)%phi = -BOGUS
gcs_t(:,:,:)%brindx=0
gcs_t(:,:,:) % clindx=0

!  Initialize the iset flag
gcs_t(:,:,:)%iset=0

!  Initialize the point to surface association
gcs_t(:,:,:)%itype=-1 !  0 - bottom, 1 - elsewhere

if(use_bottom_surf) then
  gcs_t(:,:,:)%itype=0
!  Loop over all global coordinates
  do k=$lbz,nz_tot
    gcs_t(:,:,k)%phi = gcs_t(:,:,k)%xyz(3) - z_bottom_surf   
    if(gcs_t(1,1,k)%phi <= 0.) then 
        !if(.not. brindx_set) then
        !  gcs_t(:,:,k)%brindx = -1
        !  gcs_t(:,:,k)%clindx = -1
        !endif
    endif
   
  enddo
endif
  
!  Top and bottom z-plane in gcs (same for all cylinders in generation)
do ng=1,ngen

  $if($MPI)
  if(global_rank_csp == 0) then
    write(*,*) 'generation # : ', ng
    write(*,*) 'bplane and tplane = ', tr_t(1)%gen_t(ng)%bplane, tr_t(1)%gen_t(ng)%tplane
  endif
  $else
  write(*,*) 'generation # : ', ng
  write(*,*) 'bplane and tplane = ', tr_t(1)%gen_t(ng)%bplane, tr_t(1)%gen_t(ng)%tplane
  $endif

enddo


return 

contains

$if($MPI)
!**********************************************************************
subroutine initialize_mpi_csp()
!**********************************************************************
use mpi
use types, only : rprec
use param, only : ierr, MPI_RPREC
implicit none

!--check for consistent preprocessor & param.f90 definitions of 
!  MPI and $MPI
!if (.not. USE_MPI) then
!  write (*, *) 'inconsistent use of USE_MPI and $MPI'
!  stop
!end if

call mpi_init (ierr)
call mpi_comm_size (MPI_COMM_WORLD, nproc_csp, ierr)
call mpi_comm_rank (MPI_COMM_WORLD, global_rank_csp, ierr)

!--set the MPI_RPREC variable
if (rprec == kind (1.e0)) then
  MPI_RPREC = MPI_REAL
else if (rprec == kind (1.d0)) then
  MPI_RPREC = MPI_DOUBLE_PRECISION
else
  write (*, *) 'error defining MPI_RPREC'
  stop
end if

return
end subroutine initialize_mpi_csp
$endif

!**********************************************************************
subroutine allocate_arrays()
!**********************************************************************
use param, only : ny, nz_tot
$if($MPI)
use cyl_skew_pre_base_ls, only : nx_proc
$else
use param, only : nx_proc => nx
$endif
implicit none

!  Allocate x,y,z for all coordinate systems
allocate(gcs_t(nx_proc,ny,$lbz:nz_tot))

return
end subroutine allocate_arrays

!**********************************************************************
subroutine generate_grid()
!**********************************************************************
! This subroutine generates the xyz values on all the points in the domain
! (global coordinate system) in gcs_t using the grid generation routine
! grid_build()
!
use types, only : rprec
use param, only : ny,nz_tot,dx,dy,dz
$if($MPI)
use cyl_skew_pre_base_ls, only : nx_proc, stride
$else
use param, only : nx_proc => nx
$endif

implicit none

integer :: i,j,k

do k=$lbz,nz_tot
  do j=1,ny
    do i=1,nx_proc
      $if($MPI)
      gcs_t(i,j,k)%xyz(1) = (i - 1 +  stride )*dx
      $else
      gcs_t(i,j,k)%xyz(1) = (i - 1)*dx
      $endif
      gcs_t(i,j,k)%xyz(2) = (j - 1)*dy
      gcs_t(i,j,k)%xyz(3) = (k - 0.5_rprec) * dz
    enddo
  enddo
enddo
     
return
end subroutine generate_grid

end subroutine initialize

!**********************************************************************
subroutine main_loop(nt)
!**********************************************************************
use types, only : rprec
$if($MPI)
use mpi
use param, only : ierr
use cyl_skew_pre_base_ls, only : nx_proc, global_rank_csp
$else
use param, only : nx_proc => nx
$endif
use param, only : ny, nz_tot
use cyl_skew_base_ls, only : tr_t
use cyl_skew_pre_base_ls, only : gcs_t

implicit none

integer, intent(IN) :: nt

integer :: ng, nc, nb,i,j,k
!  Loop over all global coordinates

$if ($MPI)
integer :: dumb_indx
$endif

do k=$lbz,nz_tot

  do j=1,ny

    do i=1,nx_proc

      $if ($MPI)
      !  To keep mpi stuff flowing during bad load balancing runs
      call mpi_allreduce(global_rank_csp, dumb_indx, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      $endif
        
      do ng = 1, tr_t(nt) % ngen_reslv
        
        do nc = 1, tr_t(nt)%gen_t(ng)%ncluster

          do nb=1, tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch

            if(gcs_t(i,j,k)%phi > 0._rprec) then
              call pt_loc(nt,ng,nc,nb,i,j,k)
              call point_dist(nt,ng,nc,nb,i,j,k)
              call set_sign(i,j,k)
            endif
                
          enddo
            
        enddo
      enddo
      
    enddo
  enddo
enddo



return
end subroutine main_loop


!**********************************************************************
subroutine pt_loc(nt,ng,nc,nb,i,j,k)
!**********************************************************************
use types, only : rprec
use cyl_skew_base_ls, only : tr_t, branch
use cyl_skew_base_ls, only : z_bottom_surf, use_bottom_surf
use cyl_skew_pre_base_ls, only : btw_planes, in_cir, &
  in_cyl, in_bottom_surf, in_cyl_top, in_cyl_bottom, above_cyl, below_cyl
use cyl_skew_pre_base_ls, only : gcs_t, vgcs_t, lcs_t, zrot_axis, ecs_t
implicit none

integer, intent(IN) :: nt,ng,nc,nb,i,j,k

real(rprec) :: circk, eck
real(rprec), pointer :: a => null(), b=> null()
real(rprec), pointer :: bplane => null(), tplane=> null(), skw_angle => null(), angle=>null()
real(rprec), pointer, dimension(:) :: bot=>null(), top=>null(), skw_axis=> null()
type(branch), pointer :: br_t_p => null()

!  Intialize flags
btw_planes=.false.
in_cir=.false.
in_cyl=.false.
in_bottom_surf = .false.
in_cyl_top=.false.
in_cyl_bottom=.false.

!  Set branch pointer to correct branch
br_t_p => tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb)

a         => br_t_p % a
b         => br_t_p % b
bplane    => br_t_p % bot(3)
tplane    => br_t_p % top(3)
bot       => br_t_p % bot
top       => br_t_p % top
skw_angle => br_t_p % skew_angle
skw_axis  => br_t_p % skew_axis
angle     => br_t_p % angle

!  Also check if point is below bottom surface
if(use_bottom_surf .and. ng == 1) then
  if(gcs_t(i,j,k)%xyz(3) <= z_bottom_surf) in_bottom_surf = .true.
endif

!  First check if points are between the top and bottom planes in the z - gcs
if(gcs_t(i,j,k)%xyz(3) >= bplane .and. gcs_t(i,j,k)%xyz(3) <= tplane) then
  btw_planes=.true.
elseif(gcs_t(i,j,k)%xyz(3) > tplane) then
!  Check if point is below bottom ellipse
  above_cyl = .true.
elseif(gcs_t(i,j,k)%xyz(3) < bplane) then
!  Check if point is below bottom ellipse
  below_cyl = .true.
else
  write(*,*) 'Error in pt_loc: cannot be anywhere else'
  stop
endif
      
!  Compute vector to point in the gcs from the lcs 
vgcs_t%xyz = gcs_t(i,j,k)%xyz - bot
!  Rotate gcs vector into local coordinate system
call rotation_axis_vector_3d(skw_axis, -skw_angle, vgcs_t%xyz, lcs_t%xyz)

!  Check if the point lies in the cylinder circle
 circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
if(circk <= (tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % d / 2._rprec)**2) in_cir = .true.

!  Check if point is in cylinder
if(btw_planes .and. in_cir) in_cyl = .true.

!  Check if point lies in top ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - top
call rotation_axis_vector_3d(zrot_axis, -angle, vgcs_t%xyz, ecs_t%xyz)
  
eck = ecs_t%xyz(1)**2/(a**2) + ecs_t%xyz(2)**2/(b**2)
    
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) > (tplane + bplane)/2.) in_cyl_top=.true. !  Could be below or above

!  Check if point lies in bottom ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - bot
call rotation_axis_vector_3d(zrot_axis, -angle, vgcs_t%xyz, ecs_t%xyz)
eck = ecs_t%xyz(1)**2/(a**2) + ecs_t%xyz(2)**2/(b**2)
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) < (tplane + bplane)/2.) in_cyl_bottom=.true. !  Could be below or above

!  Nullify pointers
nullify(a,b,bplane,tplane,bot,top,skw_angle,skw_axis)

return
end subroutine pt_loc

!**********************************************************************
subroutine point_dist(nt,ng,nc,nb,i,j,k)
!**********************************************************************
use types, only : rprec
use cyl_skew_base_ls, only : branch, tr_t
use cyl_skew_pre_base_ls, only : lcs_t, slcs_t, vgcs_t, sgcs_t, gcs_t
use cyl_skew_pre_base_ls, only : in_cyl_top, zrot_axis, ecs_t, eps
use cyl_skew_pre_base_ls, only : in_cyl_bottom
implicit none

integer, intent(IN) :: nt,ng,nc,nb,i,j,k
real(rprec) :: atan4, theta, dist
integer, pointer :: brindx_p => null(), clindx_p => null()
real(rprec), pointer :: a => null(), b=> null()
real(rprec), pointer :: bplane => null(), tplane=> null(), skw_angle => null(), angle=>null()
real(rprec), pointer, dimension(:) :: bot=>null(), top=>null(), skw_axis=> null()

type(branch), pointer :: br_t_p => null()

!  Set branch pointer to correct branch
br_t_p => tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb)

a         => br_t_p % a
b         => br_t_p % b
bplane    => br_t_p % bot(3)
tplane    => br_t_p % top(3)
bot       => br_t_p % bot
top       => br_t_p % top
skw_angle => br_t_p % skew_angle
skw_axis  => br_t_p % skew_axis
angle     => br_t_p % angle
brindx_p  => br_t_p % indx

clindx_p  => tr_t(nt) % gen_t(ng) % cl_t(nc) % indx

!  Compute theta value on lcs using geometry.atan4
theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

slcs_t%xyz(1) = b*dcos(theta)
slcs_t%xyz(2) = b*dsin(theta)
slcs_t%xyz(3) = lcs_t%xyz(3)

!  Rotate the surface vector in the lcs back into the gcs
call rotation_axis_vector_3d(skw_axis, skw_angle, slcs_t%xyz, vgcs_t%xyz)

sgcs_t%xyz = vgcs_t%xyz + bot !  Vector now corresponds with origin of gcs

!  Check if point on cylinder surface is between cutting planes
if(sgcs_t%xyz(3) >= bplane .and. sgcs_t%xyz(3) <= tplane) then

  call vector_magnitude_3d(lcs_t%xyz - slcs_t%xyz,dist)

  if(dist <= dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    !if(.not. brindx_set) then
    !  gcs_t(i,j,k)%brindx = brindx_p
    !  gcs_t(i,j,k)%clindx = clindx_p
    !endif
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif
else
  if(sgcs_t%xyz(3) >= tplane .and. .not. in_cyl_top) then

    vgcs_t%xyz = gcs_t(i,j,k)%xyz - top

  !  Get vector in ellipse coordinate system
    call rotation_axis_vector_3d(zrot_axis, -angle, vgcs_t%xyz, ecs_t%xyz)

    call ellipse_point_dist_2D_3(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

    call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

    if(dist <= dabs(gcs_t(i,j,k)%phi)) then
      gcs_t(i,j,k)%phi = dist
      !if(.not. brindx_set) then
      !  gcs_t(i,j,k)%brindx = brindx_p
      !  gcs_t(i,j,k)%clindx = clindx_p
      !endif
      gcs_t(i,j,k)%itype = 1
      call set_iset(i,j,k)
    endif

  elseif(sgcs_t%xyz(3) <= bplane .and. .not. in_cyl_bottom) then
    vgcs_t%xyz = gcs_t(i,j,k)%xyz - bot

  !  Get vector in ellipse coordinate system
    call rotation_axis_vector_3d(zrot_axis, -angle, vgcs_t%xyz, ecs_t%xyz)

    call ellipse_point_dist_2D_3(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

    call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

    if(dist <= dabs(gcs_t(i,j,k)%phi)) then
      gcs_t(i,j,k)%phi = dist
      !if(.not. brindx_set) then
      !  gcs_t(i,j,k)%brindx = brindx_p
      !  gcs_t(i,j,k)%clindx = clindx_p
      !endif
      gcs_t(i,j,k)%itype = 1
      call set_iset(i,j,k)
    endif

  endif 

endif

!  Check also if the point lies on the ellipses
if(in_cyl_top) then
  dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
  if(dist <= dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    !if(.not. brindx_set) then
    !  gcs_t(i,j,k)%brindx = brindx_p
    !  gcs_t(i,j,k)%clindx = clindx_p
    !endif
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif
endif

if(in_cyl_bottom) then
  dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane)
  if(dist <= dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    !if(.not. brindx_set) then
    !  gcs_t(i,j,k)%brindx = brindx_p
    !  gcs_t(i,j,k)%clindx = clindx_p
    !endif
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif
endif

!  Nullify pointers
nullify(br_t_p, a,b,bplane,tplane,bot,top,skw_angle,skw_axis,brindx_p, clindx_p)

return
end subroutine point_dist

!**********************************************************************
subroutine set_iset(i,j,k)
!**********************************************************************
use cyl_skew_pre_base_ls, only : gcs_t

implicit none

logical, parameter :: VERBOSE=.false.
integer, intent(IN) :: i,j,k

if(gcs_t(i,j,k)%iset == 1) then
  if(VERBOSE) write(*,*) 'iset already 1 - resetting phi at i,j,k : ', i,j,k
else
  gcs_t(i,j,k)%iset = 1
endif

return
end subroutine set_iset

!**********************************************************************
subroutine set_sign(i,j,k)
!**********************************************************************
use cyl_skew_pre_base_ls, only : in_cyl, in_bottom_surf, gcs_t
implicit none

integer, intent(IN) :: i,j,k
!if(gcs_t(i,j,k)%phi > 0) then
  if(in_cyl .or. in_bottom_surf) then
    gcs_t(i,j,k)%phi = -dabs(gcs_t(i,j,k)%phi)
 !   gcs_t(i,j,k)%brindx = 1
  else    
 !   gcs_t(i,j,k)%brindx = 0
  endif
!endif
return
end subroutine set_sign

!######################################################################

!**********************************************************************
subroutine compute_chi()
!**********************************************************************
!  This subroutine filters the indicator function chi
use types, only : rprec
$if($MPI)
use mpi
use param, only : ierr
use cyl_skew_pre_base_ls, only : nx_proc, global_rank_csp
$else
use param, only : nx_proc => nx
$endif
use param, only : ny, nz_tot, dz
use messages
use cyl_skew_pre_base_ls, only : gcs_t
use cyl_skew_base_ls, only : tr_t, ngen, filt_width, brindx_to_loc_id


implicit none
character (*), parameter :: sub_name = 'compute_chi'
real(rprec), dimension(:), allocatable :: z_w ! Used for checking vertical locations

integer :: i,j,k
integer :: n, nf
integer :: gen_cell_bot, gen_cell_top, gen_id
real(rprec) :: chi
real(rprec) :: z_star
real(rprec) :: zcell_bot, zcell_top
real(rprec), pointer :: bplane_p => null(), tplane_p => null()
integer, pointer, dimension(:) :: brindx_loc_id_p

$if ($MPI)
integer :: dumb_indx
$endif


nullify(brindx_loc_id_p)

allocate(z_w($lbz:nz_tot))

!  Create w-grid (physical grid)
do k=$lbz,nz_tot
  z_w(k) = gcs_t(1,1,k)%xyz(3) - dz/2.
enddo

gcs_t(:,:,:) % chi = 0._rprec

!  Do not set top most chi value; for MPI jobs
!  this is the overlap node and must be sync'd
do k=$lbz,nz_tot
  do j=1,ny

    $if ($MPI)
    !  To keep mpi stuff flowing during bad load balancing runs
    call mpi_allreduce(global_rank_csp, dumb_indx, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    $endif

    do i=1,nx_proc
          
      zcell_bot = gcs_t(i,j,k) % xyz(3) - dz/2.
      zcell_top = gcs_t(i,j,k) % xyz(3) + dz/2.
      
      call find_assoc_gen( zcell_bot, gen_cell_bot )
      call find_assoc_gen( zcell_top, gen_cell_top )

      !write(*,*) 'gen_cell_bot : ', gen_cell_bot
      !write(*,*) 'gen_cell_top : ', gen_cell_top
      
      if( gen_cell_bot == -1 .and. gen_cell_top == -1) then
      
        gcs_t(i,j,k)%chi = 0.
        
      else
      
        if( gen_cell_bot == gen_cell_top ) then
          
          call filter_chi(gcs_t(i,j,k) % xyz, gen_cell_bot, filt_width, gcs_t(i,j,k)%chi, gcs_t(i,j,k) % brindx)      

        elseif( gen_cell_bot == -1 .and. gen_cell_top .ne. -1 ) then
   
          !  Filter at ending generation
          z_star = tr_t(1)%gen_t(gen_cell_top)%bplane
          call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_cell_top, filt_width, gcs_t(i,j,k)%chi, gcs_t(i,j,k) % brindx)
          gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi * (zcell_top - z_star) / dz
          
          nf = gen_cell_top - 1 
          
          !  Filter over intermediate generations
          do n=1, nf
          
            gen_id = n
            tplane_p => tr_t(1)%gen_t(gen_id)%tplane
            bplane_p => tr_t(1)%gen_t(gen_id)%bplane
            
            z_star = 0.5_rprec * ( tplane_p + bplane_p )

            call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_id, filt_width, chi, gcs_t(i,j,k) % brindx)
            gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi + chi * (tplane_p - bplane_p) / dz
            
            nullify(tplane_p, bplane_p)
            
          enddo          
          
        elseif( gen_cell_bot .ne. -1 .and. gen_cell_top == -1 ) then 

          !  Filter at beginning generation
          z_star = tr_t(1)%gen_t(gen_cell_bot)%tplane
          call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_cell_bot, filt_width, gcs_t(i,j,k)%chi, gcs_t(i,j,k) % brindx)
          gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi * (z_star - zcell_bot) / dz
          
          nf = ngen - gen_cell_bot
          
          !  Filter over intermediate generations
          do n=1, nf
          
            gen_id = gen_cell_bot + n
            tplane_p => tr_t(1)%gen_t(gen_id)%tplane
            bplane_p => tr_t(1)%gen_t(gen_id)%bplane
            
            z_star = 0.5_rprec * ( tplane_p + bplane_p )

            call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_id, filt_width, chi, gcs_t(i,j,k) % brindx)
            gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi + chi * (tplane_p - bplane_p) / dz
            
            nullify(tplane_p, bplane_p)
            
          enddo  
          
          !nullify(tplane_p)
         
        else
        
          !  Filter at beginning generation
          z_star = tr_t(1)%gen_t(gen_cell_bot)%tplane
          call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_cell_bot, filt_width, gcs_t(i,j,k)%chi, gcs_t(i,j,k) % brindx)
          gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi * (z_star - zcell_bot) / dz
          
          !  Filter at ending generation
          z_star = tr_t(1)%gen_t(gen_cell_top)%bplane
          call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_cell_top, filt_width, chi, gcs_t(i,j,k) % brindx)
          gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi + chi * (zcell_top - z_star) / dz
          
          nf = ( gen_cell_top - gen_cell_bot + 1) - 2
          
          !  Filter over intermediate generations
          do n=1, nf
          
            gen_id = gen_cell_bot + n
            tplane_p => tr_t(1)%gen_t(gen_id)%tplane
            bplane_p => tr_t(1)%gen_t(gen_id)%bplane
            
            !  Filter at mid-height of generation
            z_star = 0.5_rprec * ( tplane_p + bplane_p )

            call filter_chi((/ gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), z_star /), gen_id, filt_width, chi, gcs_t(i,j,k) % brindx)
            gcs_t(i,j,k)%chi = gcs_t(i,j,k)%chi + chi * (tplane_p - bplane_p) / dz
            
            nullify(tplane_p, bplane_p)
            
          enddo
          
        endif
        
      endif    

    enddo
  enddo
enddo

deallocate(z_w)

!  Ensure all pointers are nullified
nullify(bplane_p,tplane_p)

!  Set clindx based on brindx
do k=$lbz,nz_tot
  do j=1,ny
    do i=1,nx_proc
    
      if(gcs_t(i,j,k) % brindx > 0) then
        brindx_loc_id_p => brindx_to_loc_id(:, gcs_t(i,j,k) % brindx )

        gcs_t(i,j,k) % clindx = tr_t(brindx_loc_id_p(1)) % gen_t(brindx_loc_id_p(2)) % cl_t(brindx_loc_id_p(3)) % indx
        nullify(brindx_loc_id_p)
      endif
    enddo
    
  enddo
  
enddo

return

end subroutine compute_chi

!**********************************************************************
subroutine find_assoc_gen(z, gen_id)
!**********************************************************************
!  This subroutine finds the generation associated with a given point
!  on the uv grid (i.e. cell center). This routine biases the bottom
!  generation when an interior interface falls within a cell

use types, only : rprec
use cyl_skew_base_ls, only : ngen, tr_t
implicit none

real(rprec), intent(in) :: z ! on uv grid

integer, intent(out) :: gen_id

integer :: ng
real(rprec), pointer :: bplane_p => null(), tplane_p => null() 

gen_id = -1

isearch_gen : do ng=1,ngen
    !  Assume all trees have same bplane and tplane values for all generations
  bplane_p => tr_t(1)%gen_t(ng)%bplane
  tplane_p => tr_t(1)%gen_t(ng)%tplane
    
  if( bplane_p <= z .and. z <= tplane_p ) then
  
    gen_id = ng
    exit isearch_gen
    
  endif
  
  nullify( bplane_p, tplane_p )

enddo isearch_gen

nullify(bplane_p, tplane_p)

return
end subroutine find_assoc_gen

!**********************************************************************
subroutine filter_chi(xyz, id_gen, delta, chi, brindx)
!**********************************************************************
!  This subroutine performs filtering in the horizontal planes
!
!  delta - filter width
!  chi   - filtered indicator function
!
use cyl_skew_pre_base_ls
implicit none

real(rprec), intent(in), dimension(3) :: xyz
integer,     intent(in)               :: id_gen
real(rprec), intent(in)               :: delta

real(rprec), intent(out)              :: chi
integer, intent(out)                  :: brindx

integer :: nt, nc

real(rprec) :: delta2, brdist

integer, pointer :: ncluster_p
type(cluster), pointer, dimension(:) :: cl_t_p

chi=0.

delta2 = delta*delta

brdist = huge(1.)
brindx = -1

do nt=1, ntree
  
  !  Loop over all branch clusters
  cl_t_p => tr_t(nt) % gen_t(id_gen) % cl_t
  ncluster_p => tr_t(nt)%gen_t(id_gen)%ncluster
  
  do nc = 1, ncluster_p
  
    !  Filter over cluster
    call filter_cl_chi(xyz, cl_t_p(nc), delta, chi, brdist, brindx) 

  enddo
  
  nullify(cl_t_p, ncluster_p)
  
enddo

!  Normalize after completing integration
!chi = chi/(2._rprec*pi*delta2)
chi = 6._rprec/(pi*delta2) * chi

return

end subroutine filter_chi

!**********************************************************************
subroutine filter_cl_chi(xyz, cl_t, delta, chi, brdist, brindx)
!**********************************************************************
!  Also assigns branch index consistent with filtering of chi
!
use types, only : rprec, vec3d
use cyl_skew_base_ls, only : cluster, branch, point_2d
use cyl_skew_pre_base_ls, only : zrot_axis
implicit none

real(rprec),  intent(in), dimension(3) :: xyz
type(cluster), target, intent(in)      :: cl_t
real(rprec), intent(in)                :: delta
real(rprec), intent(inout)             :: chi
real(rprec), intent(inout)             :: brdist
integer, intent(inout)                 :: brindx

integer     :: nb
real(rprec) :: chi_int, ds
real(rprec) :: brdist_check

integer, pointer                   :: nbranch_p
real(rprec), pointer               :: skew_angle_p, angle_p
real(rprec), pointer, dimension(:) :: angle2_p
real(rprec), pointer, dimension(:) :: a_p, b_p
real(rprec), pointer, dimension(:) :: bot_p, top_p, cl_origin_p
integer, pointer                   :: indx_p


type(point_2d), allocatable, dimension(:) :: lpnt_t, cpnt_t
real(rprec), dimension(3) :: xyz_rot !  point used for rotatations 

type(vec3d) :: svec_t

type(branch), pointer, dimension(:) :: br_t_p

! Nullify all pointers
nullify(nbranch_p, cl_origin_p, br_t_p)
nullify(top_p, bot_p, skew_angle_p, angle_p)
nullify(indx_p)

nbranch_p   => cl_t % nbranch
cl_origin_p => cl_t % origin
br_t_p      => cl_t % br_t
      
allocate(lpnt_t(nbranch_p), cpnt_t(nbranch_p))

!  Set all branch settings
do nb = 1, nbranch_p

  top_p        => br_t_p(nb) % top
  bot_p        => br_t_p(nb) % bot
  skew_angle_p => br_t_p(nb) % skew_angle
  angle_p      => br_t_p(nb) % angle
  indx_p       => br_t_p(nb) % indx

  svec_t % xyz = top_p - bot_p

  call vector_magnitude_3d(svec_t % xyz, svec_t % mag)
  
  ! find the projection length along the branch axis      
  ds = ( xyz(3) - bot_p(3) ) / (svec_t % mag * dcos( skew_angle_p )) 
  
  ! center point of area to integrate over in z-plane
  cpnt_t(nb) % xy = bot_p(1:2) + ds * svec_t%xyz(1:2) 
    
  !    write(*,'(1a,3f12.6)') 'xyz ', xyz
  !    write(*,'(1a,3f12.6)') 'xyz_c : ', xyz_c
   
  !  Compute local vector to branch coordinate system (must be in 3D); project to z=0
  xyz_rot = (/ xyz(1), xyz(2), 0._rprec /) - (/ cpnt_t(nb) % xy(1), cpnt_t(nb) % xy(2), 0._rprec /)
  
  !  Compute magnitude of distance from point to center of ellipse
  call vector_magnitude_3d(xyz_rot, brdist_check)
  !write(*,*) 'brdist_check : ', brdist_check
  if( brdist_check < brdist ) then
    brdist = brdist_check
    brindx = indx_p ! This is the new closest branch
  endif
 
  !   write(*,'(1a,f12.6)') 'zrot_t(id_gen)%angle(n)*180/pi : ', zrot_t(id_gen)%angle(n)*180./pi
  !  Perform rotation of local vector about z-axis into ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -angle_p, xyz_rot, xyz_rot)
  
  ! 2d local point relative to ellipse coordinate system
  lpnt_t(nb) % xy = xyz_rot(1:2) 
  
  nullify(top_p, bot_p, skew_angle_p, angle_p)
  nullify(indx_p)
  
enddo

a_p      => br_t_p % a
b_p      => br_t_p % b
angle2_p => br_t_p % angle

!  Perform weighted integration over branch cluster
call weighted_cl_chi_int(a_p, b_p, angle2_p, cpnt_t, lpnt_t, nbranch_p, delta, chi_int)

nullify(a_p, b_p, angle2_p)
deallocate(cpnt_t, lpnt_t)

chi = chi + chi_int

return

end subroutine filter_cl_chi

!**********************************************************************
subroutine weighted_cl_chi_int(a, b, angle, cpnt_t, lpnt_t, nbranch, delta, chi)
!**********************************************************************
use types, only : rprec
use cyl_skew_base_ls, only : point_2d, point_3d
use cyl_skew_pre_base_ls, only : zrot_axis
!  Does not normalize
implicit none

integer, intent(in) :: nbranch
real(rprec), dimension(nbranch), intent(in) :: a, b, angle
type(point_2d), dimension(nbranch), intent(in) :: cpnt_t ! center of ellipse (in global coordinate system)
type(point_2d), dimension(nbranch), intent(in) :: lpnt_t ! local point relative to ellipse
real(rprec), intent(in) :: delta
real(rprec), intent(out) :: chi

integer, parameter :: Nx=20, Ny=20

logical :: inside_self, inside_other

integer :: i, j, n, nm

real(rprec) :: delta2, dist, dist2
real(rprec), allocatable, dimension(:) :: dx, dy

type(point_2d) :: cell_t
type(point_3d) :: cell_gcs_t
type(point_3d) :: test_point_t ! used for test if integrand points lie in multiple elipses


! we assume a, b, etc. are the same length

allocate(dx(nbranch), dy(nbranch)) ! Each branch recieves its own
!allocate(cell_center_t(nbranch)) ! center of integration cells

do n = 1, nbranch
  dx(n) = 2._rprec*a(n) / Nx
  dy(n) = 2._rprec*b(n) / Ny 
enddo

delta2 = delta*delta

chi=0.

!  Now for each branch check if inside itself and other
do n = 1, nbranch

  do j = 1, ny

  !  y-value of i,j cell center for each branch (relative to ellipse coordinate system)
    cell_t % xy(2) = -b(n) + (j - 0.5)*dy(n)
  
    do i = 1, nx
  
      cell_t % xy(1) = -a(n) + (i - 0.5)*dx(n) ! (relative to ellipse coordinate system)

      inside_self = .false.
      inside_other = .false.

      call ellipse_contains_point_2d(a(n), b(n), cell_t % xy, (/ 0._rprec, 0._rprec /), inside_self)
      !if((cell_t%xy(1)/a(n))**2 + (cell_t%xy(2)/b(n))**2 <= 1._rprec) inside_self = .true.
      
      if(inside_self) then
        
        if(n > 1) then
        
          nm = n ! should skip inside_other_chk if n = 1
          
          ! Need to rotate cell_center_t vector into global coordinate system
          call rotation_axis_vector_3d(zrot_axis, angle(n), &
            (/ cell_t % xy(1), cell_t %xy(2), 0._rprec /), cell_gcs_t % xyz)
             
          inside_other_chk : do while ( nm > 1 )
        
            nm = nm - 1
                        
            !  xc, yc is relative to each ellipse center (cpnt_t)
            !  compute center relative to previous ellipse center
            test_point_t % xyz(1:2) = cell_gcs_t % xyz(1:2) + cpnt_t(n) % xy - cpnt_t(nm) % xy 
            test_point_t % xyz(3) = 0._rprec 
          
          !  Rotate test point in to ellipse nm's coordinate system
            call rotation_axis_vector_3d(zrot_axis, -angle(nm), test_point_t % xyz, test_point_t % xyz)

            call ellipse_contains_point_2d(a(nm), b(nm), test_point_t % xyz(1:2), (/ 0._rprec, 0._rprec /), inside_other)
            !if((test_point_t % xyz(1)/a(nm))**2 + (test_point_t % xyz(2)/b(nm))**2 <= 1._rprec) inside_other = .true.       
          
            if(inside_other) then

              exit inside_other_chk
            endif
          
          enddo inside_other_chk
          
        endif
       
        if (.not. inside_other) then
        
        !  distance from cell center and specified point
          call vector_magnitude_2d((/ lpnt_t(n)%xy(1) - cell_t % xy(1), &
            lpnt_t(n)%xy(2) - cell_t % xy(2) /), dist)
                   
          dist2 = dist*dist

          chi = chi + dx(n)*dy(n)*exp(-6._rprec*dist2/delta2)
          
        endif
        
      endif
      
    enddo
    
  enddo
  
enddo

return
end subroutine weighted_cl_chi_int


!############################################################################################################

!############################################################################################################
!**********************************************************************
subroutine finalize()
!**********************************************************************
$if ($MPI)
use mpi_defs
use param, only : ierr
$endif
use cyl_skew_pre_base_ls
use param, only : nproc, coord

implicit none

if(DIST_CALC) call write_output()

$if ($MPI)
!  Finalize mpi communication
call MPI_FINALIZE(ierr)
$endif

return
contains

!**********************************************************************
subroutine write_output()
!**********************************************************************
$if($MPI)
use mpi
$endif
use param, only : ld, nx, ny, nz, nz_tot,dx,dy
$if($MPI)
use param, only :  MPI_RPREC,ierr, nproc, status
use cyl_skew_pre_base_ls, only : nx_proc
$endif
!use cyl_skew_base_ls, only : filter_chi, brindx_to_loc_id, tr_t
use strmod, only : numtostr
implicit none
include 'tecio.h'

character (64) :: fname, fname_phi, fname_brindx, fname_clindx, fname_chi, temp
integer :: i,j,k,n

integer :: istart, iend
integer :: kstart, kend

integer, pointer, dimension(:,:,:) :: brindx, clindx
real(rprec), pointer, dimension(:,:,:) :: phi, chi

$if ($MPI)
integer, allocatable, dimension(:,:,:) :: brindx_proc, clindx_proc
real(rprec), allocatable, dimension(:,:,:) :: rbrindx_proc, rclindx_proc
real(rprec), allocatable, dimension(:,:,:) :: phi_proc, chi_proc
integer :: sendcnt, recvcnt
$endif

real(rprec), allocatable, dimension(:) :: x, y, z
integer, pointer, dimension(:) :: br_loc_id_p

$if($MPI)
sendcnt = nx_proc * ny * (nz_tot - $lbz + 1)
!gcs_t(:,:,$lbz)%phi = BOGUS
$endif

$if($MPI)
!  Gather data one by one in rank order
if( global_rank_csp > 0 ) then

  !write(*,*) 'global_rank_csp, sendcnt : ',  global_rank_csp, sendcnt
    !  Open file which to write global data
    write (fname,*) 'x.dat'
    fname = trim(adjustl(fname)) 

    write (temp, '(".c",i0)') global_rank_csp
    fname = trim (fname) // temp  
    open (1, file=fname, form='formatted',status='unknown',position='rewind')
    do n=1,nx_proc
      write(1,*) gcs_t(n,1,1)%xyz(1)
    enddo
    close (1)

  call mpi_ssend( gcs_t(:,:,:)%phi, sendcnt, MPI_RPREC, 0, 1, &
    MPI_COMM_WORLD, ierr)
  call mpi_ssend( gcs_t(:,:,:)%chi, sendcnt, MPI_RPREC, 0, 2, &
    MPI_COMM_WORLD, ierr)
  call mpi_ssend( gcs_t(:,:,:)%brindx, sendcnt, MPI_INTEGER, 0, 3, &
    MPI_COMM_WORLD, ierr)
  call mpi_ssend( gcs_t(:,:,:)%clindx, sendcnt, MPI_INTEGER, 0, 4, &
    MPI_COMM_WORLD, ierr)    
  call mpi_send( gcs_t(:,1,1)%xyz(1), nx_proc, MPI_RPREC, 0, 5, &
    MPI_COMM_WORLD, ierr)


else

  allocate(phi(ld,ny,$lbz:nz_tot))
  allocate(chi(ld,ny,$lbz:nz_tot))
  allocate(brindx(ld,ny,$lbz:nz_tot))
  allocate(clindx(ld,ny,$lbz:nz_tot))
  allocate(x(nx+1),y(ny+1),z($lbz:nz_tot))
 
  !  Initialize
  phi=BOGUS
  chi=BOGUS
  brindx=0
  clindx=0
  x=-1._rprec
  y=-1._rprec
  z=-1._rprec

  !  Get from proc 0
  phi(1:nx_proc,:,:) = gcs_t(:,:,:)%phi
  chi(1:nx_proc,:,:) = gcs_t(:,:,:)%chi
  brindx(1:nx_proc,:,:) = gcs_t(:,:,:)%brindx
  clindx(1:nx_proc,:,:) = gcs_t(:,:,:)%clindx
  x(1:nx_proc) = gcs_t(:,1,1)%xyz(1)
  y(1:ny) = gcs_t(1,:,1)%xyz(2)
  z = gcs_t(1,1,:)%xyz(3)
        
  !  Get from all other procs
  do n=1,nproc_csp-1

    istart = nx_proc*n + 1
    if( n == nproc_csp - 1 ) then
      iend = nx
    else                                         
      iend = nx_proc*(n+1)
    endif

    recvcnt = (iend-istart+1) * ny * (nz_tot - $lbz + 1)
    write(*,*) 'n, istart, iend, recvcnt : ', n, istart, iend, recvcnt
    ierr=1
    call mpi_recv( phi(istart:iend,:,:), recvcnt, MPI_RPREC, n, 1, &
      MPI_COMM_WORLD, status, ierr)
    call mpi_recv( chi(istart:iend,:,:), recvcnt, MPI_RPREC, n, 2, &
      MPI_COMM_WORLD, status, ierr)
    call mpi_recv( brindx(istart:iend,:,:), recvcnt, MPI_INTEGER, n, 3, &
      MPI_COMM_WORLD, status, ierr)            
    call mpi_recv( clindx(istart:iend,:,:), recvcnt, MPI_INTEGER, n, 4, &
      MPI_COMM_WORLD, status, ierr)
    call mpi_recv( x(istart:iend), iend-istart+1, MPI_RPREC, n, 5, &
      MPI_COMM_WORLD, status, ierr)

    if( ierr /= 0 ) then
       write(*,*) 'mpi_recv error wiht proc : ',ierr, n
    endif

  enddo

  !  Set buffer x and y points for periodicity
  x(nx+1) = x(nx)+dx
  y(ny+1) = y(ny)+dy

endif  
write(*,*) 'Finalized local to global send/receive'

$else


allocate(phi(ld,ny,nz))
allocate(chi(ld,ny,nz))
allocate(brindx(ld,ny,nz))
allocate(clindx(ld,ny,nz))
allocate(x(nx+1),y(ny+1),z(nz))

phi(1:nx,:,:) = gcs_t(1:nx,:,:)%phi
chi(1:nx,:,:) = gcs_t(1:nx,:,:)%chi
brindx(1:nx,:,:) = gcs_t(1:nx,:,:)%brindx
clindx(1:nx,:,:) = gcs_t(1:nx,:,:)%clindx
x(1:nx) = gcs_t(1:nx,1,1)%xyz(1)
y(1:nx) = gcs_t(1,1:ny,1)%xyz(2)
z(:) = gcs_t(1,1,:)%xyz(3)

!  Set buffer x and y points for periodicity
x(nx+1) = x(nx)+dx
y(ny+1) = y(ny)+dy

$endif

deallocate( gcs_t )

$if($MPI)
if( global_rank_csp == 0 ) then

!  Write processor files for lesgo
  allocate(phi_proc(ld,ny,$lbz:nz))
  allocate(rbrindx_proc(ld,ny,$lbz:nz))
  allocate(rclindx_proc(ld,ny,$lbz:nz))
  allocate(chi_proc(ld,ny,$lbz:nz))

  allocate(brindx_proc(ld,ny,$lbz:nz))
  allocate(clindx_proc(ld,ny,$lbz:nz))
  
  !  Write data files
  do n=0,nproc-1

    kstart = (nz-1)*n
    kend = kstart + nz

    write(*,*) 'n, nz, nz_tot, kstart, kend :', n, nz, nz_tot, kstart, kend    

    phi_proc = BOGUS
    rbrindx_proc=0._rprec
    rclindx_proc=0._rprec    
    brindx_proc = 0
    clindx_proc = 0
    chi_proc = BOGUS

    if(kend-kstart+1 /= nz-$lbz+1) then
      write(*,*) 'z dimension for proc ',n,' not specified correctly'
    endif
  
    phi_proc(:,:,:) = phi(:,:,kstart:kend)
    rbrindx_proc(:,:,:) = 1.*brindx(:,:,kstart:kend)
    rclindx_proc(:,:,:) = 1.*clindx(:,:,kstart:kend)    
    brindx_proc(:,:,:) = brindx(:,:,kstart:kend)
    clindx_proc(:,:,:) = clindx(:,:,kstart:kend) 
    chi_proc(:,:,:) = chi(:,:,kstart:kend) 

    !  Open file which to write global data
    write (fname,*) 'cyl_skew_ls.dat'
    fname = trim(adjustl(fname)) 

    write (temp, '(".c",i0)') n
    fname = trim (fname) // temp

    call write_tecplot_header_ND(fname, 'rewind', 7, &
      (/ Nx+1, Ny+1, Nz-$lbz+1 /), &
      '"x", "y", "z", "phi", "brindx", "clindx", "chi"', &
      numtostr(n,6), 2)

    call write_real_data_3D( fname, 'append', 'formatted', 4, Nx, Ny, nz-$lbz+1,&
      (/phi_proc(1:nx,:,:), rbrindx_proc(1:nx,:,:), rclindx_proc(1:nx,:,:), chi_proc(1:nx,:,:)/),&
      4, x, y, z(kstart:kend))    
 
    !  Open file which to write global data
    write (fname_phi,*) 'phi.out'
    fname_phi = trim(adjustl(fname_phi)) 
    write (fname_brindx,*) 'brindx.out'
    fname_brindx = trim(adjustl(fname_brindx)) 
    write (fname_clindx,*) 'clindx.out'
    fname_clindx = trim(adjustl(fname_clindx))     
    write (fname_chi,*) 'chi.out'
    fname_chi = trim(adjustl(fname_chi)) 
 

    write (temp, '(".c",i0)') n
    fname_phi = trim (fname_phi) // temp
    write (temp, '(".c",i0)') n
    fname_brindx = trim (fname_brindx) // temp
    write (temp, '(".c",i0)') n
    fname_clindx = trim (fname_clindx) // temp
    write (temp, '(".c",i0)') n
    fname_chi = trim (fname_chi) // temp 

    !  Write binary data for lesgo
    $if ($WRITE_BIG_ENDIAN)
    open (1, file=fname_phi, form='unformatted', convert='big_endian')
    $elseif ($WRITE_LITTLE_ENDIAN)
    open (1, file=fname_phi, form='unformatted', convert='little_endian')
    $else
    open (1, file=fname_phi, form='unformatted')
    $endif
    write(1) phi_proc
    close (1)

    !  Write binary data for lesgo
    $if ($WRITE_BIG_ENDIAN)
    open (1, file=fname_brindx, form='unformatted', convert='big_endian')
    $elseif ($WRITE_LITTLE_ENDIAN)
    open (1, file=fname_brindx, form='unformatted', convert='little_endian')
    $else
    open (1, file=fname_brindx, form='unformatted')
    $endif
    write(1) brindx_proc
    close (1) 

    !  Write binary data for lesgo
    $if ($WRITE_BIG_ENDIAN)
    open (1, file=fname_clindx, form='unformatted', convert='big_endian')
    $elseif ($WRITE_LITTLE_ENDIAN)
    open (1, file=fname_clindx, form='unformatted', convert='little_endian')
    $else
    open (1, file=fname_clindx, form='unformatted')
    $endif
    write(1) clindx_proc
    close (1)  

    !  Write binary data for lesgo
    $if ($WRITE_BIG_ENDIAN)
    open (1, file=fname_chi, form='unformatted', convert='big_endian')
    $elseif ($WRITE_LITTLE_ENDIAN)
    open (1, file=fname_chi, form='unformatted', convert='little_endian')
    $else
    open (1, file=fname_chi, form='unformatted')
    $endif
    write(1) chi_proc
    close (1)    

  enddo

  deallocate(phi,phi_proc)
  deallocate(rbrindx_proc, rclindx_proc)
  deallocate(brindx,brindx_proc)
  deallocate(clindx,clindx_proc)
  deallocate(chi,chi_proc)
 

endif  

$else
  
  !write(*,*) 'No output yet for single processor'
  
  !  Open file which to write global data
  write (fname,*) 'cyl_skew_ls.dat'
  fname = trim(adjustl(fname)) 

  call write_tecplot_header_ND(fname, 'rewind', 7, &
    (/ Nx+1, Ny+1, Nz /), &
    '"x", "y", "z", "phi", "brindx", "clindx", "chi"', &
    numtostr(n,6), 2)

  call write_real_data_3D( fname, 'append', 'formatted', 4, Nx, Ny, Nz,&
    (/phi(1:nx,:,1:Nz), real(brindx(1:nx,:,1:Nz),rprec), real(clindx(1:nx,:,1:Nz),rprec), chi(1:nx,:,1:Nz)/),&
    4, x, y, z(1:Nz))    

  !  Open file which to write global data
  write (fname_phi,*) 'phi.out'
  fname_phi = trim(adjustl(fname_phi)) 
  write (fname_brindx,*) 'brindx.out'
  fname_brindx = trim(adjustl(fname_brindx)) 
  write (fname_clindx,*) 'clindx.out'
  fname_clindx = trim(adjustl(fname_clindx))     
  write (fname_chi,*) 'chi.out'
  fname_chi = trim(adjustl(fname_chi)) 
 
  !  Write binary data for lesgo
  $if ($WRITE_BIG_ENDIAN)
  open (1, file=fname_phi, form='unformatted', convert='big_endian')
  $elseif ($WRITE_LITTLE_ENDIAN)
  open (1, file=fname_phi, form='unformatted', convert='little_endian')
  $else
  open (1, file=fname_phi, form='unformatted')
  $endif
  write(1) phi
  close (1)

  !  Write binary data for lesgo
  $if ($WRITE_BIG_ENDIAN)
  open (1, file=fname_brindx, form='unformatted', convert='big_endian')
  $elseif ($WRITE_LITTLE_ENDIAN)
  open (1, file=fname_brindx, form='unformatted', convert='little_endian')
  $else
  open (1, file=fname_brindx, form='unformatted')
  $endif
  write(1) brindx
  close (1) 

  !  Write binary data for lesgo
    $if ($WRITE_BIG_ENDIAN)
    open (1, file=fname_clindx, form='unformatted', convert='big_endian')
    $elseif ($WRITE_LITTLE_ENDIAN)
    open (1, file=fname_clindx, form='unformatted', convert='little_endian')
    $else
    open (1, file=fname_clindx, form='unformatted')
    $endif
    write(1) clindx
    close (1)  

    !  Write binary data for lesgo
    $if ($WRITE_BIG_ENDIAN)
    open (1, file=fname_chi, form='unformatted', convert='big_endian')
    $elseif ($WRITE_LITTLE_ENDIAN)
    open (1, file=fname_chi, form='unformatted', convert='little_endian')
    $else
    open (1, file=fname_chi, form='unformatted')
    $endif
    write(1) chi
    close (1)    

  deallocate(phi)
  deallocate(brindx)
  deallocate(clindx)
  deallocate(chi)  

$endif

!  Generate generation associations to be used in drag force calculations
!  for each generation
!call gen_assoc() !  Generation data from the last tree must match that of the first

return
end subroutine write_output

!**********************************************************************
subroutine gen_assoc()
!**********************************************************************
!
!  This subroutine is used to find where each generation lives at. The
!  information created from this routine is used in lesgo for computing
!  drag force data for individual generations. In order to use this
!  capability the Makefile flag should be set to USE_CYLINDER_SKEW=yes
!           
use param, only : nz,dz

implicit none
character(64) :: fname, temp
integer :: ng, k

real(rprec) :: bplane, tplane
integer, dimension(:), allocatable :: igen, kbottom, kbottom_inside, ktop, ktop_inside
real(rprec), dimension(:), allocatable :: gcs_w, dz_bottom, dz_top

allocate(gcs_w(nz))
allocate(igen(ngen))
allocate(kbottom(ngen), kbottom_inside(ngen))
allocate(ktop(ngen), ktop_inside(ngen))
allocate(dz_bottom(ngen), dz_top(ngen))


!  Create w-grid (physical grid)
do k=1,nz
  gcs_w(k) = gcs_t(1,1,k)%xyz(3) - dz/2.
enddo

do ng=1,ngen

    !  Assume all trees have same bplane and tplane values for all generations
    bplane = tr_t(1)%gen_t(ng)%bplane
    tplane = tr_t(1)%gen_t(ng)%tplane

  if(bplane < gcs_w(1) .and. tplane < gcs_w(1)) then
    igen(ng) = -1
    kbottom(ng) = -1
    ktop(ng) = -1
    kbottom_inside(ng) = 0
    ktop_inside(ng) = 0
    dz_bottom(ng) = 0.
    dz_top(ng) = 0.
  elseif(bplane > gcs_w(nz) .and. tplane > gcs_w(nz)) then
    igen(ng) = -1
    kbottom(ng) = -1
    ktop(ng) = -1
    kbottom_inside(ng) = 0
    ktop_inside(ng) = 0
    dz_bottom(ng) = 0.
    dz_top(ng) = 0.
  else
    igen(ng) = ng
  !  Perform kbottom, kbottom_inside, ktop, ktop_inside search
    if(bplane < gcs_w(1)) then
      kbottom(ng) = -1
      kbottom_inside(ng) = 0
      dz_bottom(ng) = 0.
    else
      isearch_bottom: do k=2,nz
        if(gcs_w(k) > bplane) then
          kbottom(ng) = k-1
          kbottom_inside(ng) = 1
          dz_bottom(ng) = gcs_w(k) - bplane
          exit isearch_bottom
        endif
      enddo isearch_bottom
    endif
    if(tplane > gcs_w(nz)) then
      ktop(ng) = -1
      ktop_inside(ng) = 0
      dz_top(ng) = 0.
    else
      isearch_top: do k=2,nz
        if(gcs_w(k) >= tplane) then
          ktop(ng) = k-1
          ktop_inside(ng) = 1
          dz_top(ng) = tplane - gcs_w(k-1)
          exit isearch_top
        endif
      enddo isearch_top
    endif
    if(ng == 1) call point_assoc() !  For gen-1 only check point association with the ground
  endif
enddo

!  Open file which to write global data
write (fname,*) 'cyl_skew_gen_ls.out'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')
do ng=1,ngen
  write(2,*) igen(ng), kbottom_inside(ng), kbottom(ng), dz_bottom(ng), ktop_inside(ng), ktop(ng), dz_top(ng)
enddo
close(2)

deallocate(gcs_w)
deallocate(igen)
deallocate(kbottom, kbottom_inside)
deallocate(ktop, ktop_inside)
deallocate(dz_bottom, dz_top)

return

end subroutine gen_assoc

!**********************************************************************
subroutine point_assoc()
!**********************************************************************
use param, only : nx, ny, nz
$if ($MPI)
use param, only : nproc, coord
$endif
implicit none

character(64) :: fname, temp
integer :: i,j,k

!  Open file which to write global data
write (fname,*) 'cyl_skew_point_ls.out'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')
do k=1,nz
  do j = 1,ny
    do i = 1,nx+2
      write(2,*) gcs_t(i,j,k)%itype
    enddo
  enddo
enddo
   
close(2)

return
end subroutine point_assoc

end subroutine finalize
!############################################################################################################
