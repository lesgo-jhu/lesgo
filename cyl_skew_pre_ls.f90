!**********************************************************************
module cyl_skew_pre_base_ls
!**********************************************************************
use types, only : rprec, vec3d
use param, only : pi
use cyl_skew_base_ls
use cyl_skew_ls, only : fill_tree_array_ls
use io, only : write_tecplot_header_xyline, write_tecplot_header_ND
use io, only : write_real_data, write_real_data_1D, write_real_data_2D, write_real_data_3D

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

!  cs{0,1} all correspond to vectors with the origin at the
!  corresponding coordinate system
type(cs0), target, allocatable, dimension(:,:,:) :: gcs_t
type(cs1) :: lcs_t, slcs_t, sgcs_t, ecs_t
!type(cs2), allocatable, dimension(:) :: lgcs_t
!type(cs2), allocatable, dimension(:,:) ::  ebgcs_t, etgcs_t ! Shape (ntree, ngen)
!type(rot), allocatable, dimension(:) :: zrot_t
!  vectors do not have starting point a origin of corresponding
!  coordinate system
type(vec3d) :: vgcs_t

logical, parameter :: mpi_split = .true.
integer, parameter :: mpi_split_nproc = 16

logical :: DIST_CALC=.true.

real(rprec), parameter :: BOGUS = 1234567890._rprec
real(rprec), parameter :: iBOGUS = 1234567890
real(rprec), parameter :: eps = 1.e-12
real(rprec), parameter, dimension(3) :: zrot_axis = (/0.,0.,1./)

!real(rprec), dimension(3,ntree) :: origin

logical :: in_cir, in_cyl
logical :: in_cyl_top, in_cyl_bottom
logical :: above_cyl, below_cyl
logical :: in_bottom_surf, btw_planes

integer, dimension(3) :: cyl_loc

!logical :: brindx_set = .false. !  Used for finding branch closest to each point

!integer, allocatable, dimension(:) :: gen_ntrunk, gen_ncluster
!real(rprec), allocatable, dimension(:) :: crad, clen, rad_offset

end module cyl_skew_pre_base_ls

!**************************************************************
program cyl_skew_pre_ls
!***************************************************************
$if ($MPI)
use param, only : coord
$endif
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
$if ($MPI)
use mpi_defs
$endif

use param, only : nz, coord
use cyl_skew_pre_base_ls, only : gcs_t, BOGUS
use cyl_skew_base_ls, only : use_bottom_surf, z_bottom_surf, ngen, tr_t
use cyl_skew_ls, only : fill_tree_array_ls

implicit none

integer :: ng,i,j,k

$if ($MPI)
call initialize_mpi ()
$endif
call allocate_arrays()
call fill_tree_array_ls()
call generate_grid()

!!  Allocate x,y,z for all coordinate systems
!allocate(gcs_t(nx+2,ny,$lbz:nz))

!  Initialize the distance function
gcs_t(:,:,:)%phi = BOGUS
!!  Set lower level
!gcs_t(:,:,0)%phi = -BOGUS

gcs_t(:,:,:)%brindx=0
gcs_t % clindx = 0


!  Initialize the iset flag
gcs_t(:,:,:)%iset=0

!  Initialize the point to surface association
gcs_t(:,:,:)%itype=-1 !  0 - bottom, 1 - elsewhere

if(use_bottom_surf) then
  gcs_t(:,:,:)%itype=0
!  Loop over all global coordinates
  do k=$lbz,Nz
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

  if(coord == 0) then
    write(*,*) 'generation # : ', ng
    write(*,*) 'bplane and tplane = ', tr_t(1)%gen_t(ng)%bplane, tr_t(1)%gen_t(ng)%tplane
  endif

enddo


return 

contains

!**********************************************************************
subroutine allocate_arrays()
!**********************************************************************
use param, only : ld, nx, ny
implicit none

!  Allocate x,y,z for all coordinate systems
allocate(gcs_t(ld,ny,$lbz:nz))

return
end subroutine allocate_arrays

!**********************************************************************
subroutine generate_grid()
!**********************************************************************
! This subroutine generates the xyz values on all the points in the domain
! (global coordinate system) in gcs_t using the grid generation routine
! grid_build()
!
use param, only : nproc, coord, nx, ny, nz
use grid_defs

implicit none

if(.not. grid_built) call grid_build()

do k=$lbz,nz
  do j=1,ny
    do i=1,nx
      gcs_t(i,j,k)%xyz(1) = x(i)
      gcs_t(i,j,k)%xyz(2) = y(j)
      gcs_t(i,j,k)%xyz(3) = z(k)
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
use param, only : nx, ny, nz
use cyl_skew_base_ls, only : tr_t
use cyl_skew_pre_base_ls, only : gcs_t
implicit none

integer, intent(IN) :: nt

integer :: ng, nc, nb,i,j,k
!  Loop over all global coordinates


do k=$lbz,Nz

  do j=1,ny

    do i=1,nx+2
        
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
use param, only : nx, ny, nz, dz
use messages
use cyl_skew_pre_base_ls, only : gcs_t
use cyl_skew_base_ls, only : tr_t, ngen, filt_width, brindx_to_loc_id
$if($MPI)
use mpi
use param, only : coord, nproc, comm, ierr
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none
character (*), parameter :: sub_name = 'compute_chi'
real(rprec), dimension(:), allocatable :: z_w ! Used for checking vertical locations

integer :: i,j,k,ubz
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

allocate(z_w($lbz:nz))

$if ($MPI)
  !if(coord == nproc - 1) then
  !  ubz = nz
  !else
    ubz = nz - 1
  !endif
$else
  ubz = nz
$endif

!  Create w-grid (physical grid)
do k=$lbz,nz
  z_w(k) = gcs_t(1,1,k)%xyz(3) - dz/2.
enddo

gcs_t(:,:,:) % chi = 0._rprec

!  Do not set top most chi value; for MPI jobs
!  this is the overlap node and must be sync'd
do k=1,ubz
  do j=1,ny
    do i=1,nx
    
      $if ($MPI)
      !  To keep mpi stuff flowing during bad load balancing runs
      call mpi_allreduce(coord, dumb_indx, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      $endif
      
      zcell_bot = gcs_t(i,j,k) % xyz(3) - dz/2.
      zcell_top = gcs_t(i,j,k) % xyz(3) + dz/2.
      
      call find_assoc_gen( zcell_bot, gen_cell_bot )
      call find_assoc_gen( zcell_top, gen_cell_top )

      write(*,*) 'gen_cell_bot : ', gen_cell_bot
      write(*,*) 'gen_cell_top : ', gen_cell_top
      
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

!  Now must sync all overlapping nodes
$if ($MPI)
call mpi_sync_real_array(gcs_t(:,:,:)%chi, MPI_SYNC_DOWNUP)
$endif

!  Ensure all pointers are nullified
nullify(bplane_p,tplane_p)

!  Set clindx based on brindx
do k=$lbz,ubz
  do j=1,ny
    do i=1,nx
    
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
use grid_defs
use param, only : ld, Nx, Ny, Nz
use cyl_skew_base_ls, only : filter_chi
implicit none

character (64) :: fname, temp
integer :: i,j,k

integer, pointer, dimension(:,:,:) :: brindx, clindx
real(rprec), pointer, dimension(:,:,:) :: phi, chi

integer, pointer, dimension(:) :: br_loc_id_p

nullify(phi,brindx)
nullify(br_loc_id_p)
allocate(phi(ld,ny,$lbz:nz))
allocate(brindx(ld,ny,$lbz:nz))
allocate(clindx(ld,ny,$lbz:nz))
allocate(chi(ld,ny,$lbz:nz))

if(nproc > 1 .and. coord == 0) gcs_t(:,:,$lbz)%phi = -BOGUS

!  Open file which to write global data
write (fname,*) 'cyl_skew_ls.dat'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif
!  Create tecplot formatted phi and brindx field file
open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')

write(2,*) 'variables = "x", "y", "z", "phi", "brindx", "clindx", "itype", "chi"';

write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &

1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz-$lbz+1

write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''

do k=$lbz,nz
  do j=1,ny
    do i=1,nx

      !  Update clindx based on brindx
      if(gcs_t(i,j,k)%brindx > 0 ) then
        br_loc_id_p => brindx_to_loc_id(:, gcs_t(i,j,k)%brindx)
     
        gcs_t(i,j,k)%clindx = tr_t(br_loc_id_p(1)) % gen_t(br_loc_id_p(2)) % cl_t(br_loc_id_p(3)) % indx
        nullify( br_loc_id_p )
      else
      
        gcs_t(i,j,k)%clindx = -1
        
      endif
      
      write(2,*) gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), gcs_t(i,j,k)%xyz(3), &
        gcs_t(i,j,k)%phi, gcs_t(i,j,k)%brindx, gcs_t(i,j,k)%clindx, &
        gcs_t(i,j,k)%itype, gcs_t(i,j,k)%chi
        
    enddo
  enddo
enddo
close(2)

do k=$lbz,nz
  do j = 1,ny
    do i = 1,ld
      phi(i,j,k) = gcs_t(i,j,k)%phi
      brindx(i,j,k) = gcs_t(i,j,k)%brindx
      clindx(i,j,k) = gcs_t(i,j,k)%clindx
     
      !if(brindx(i,j,k) > 0 ) then
      !  br_loc_id_p => brindx_to_loc_id(:, brindx(i,j,k))
     
      !  clindx(i,j,k) = tr_t(br_loc_id_p(1)) % gen_t(br_loc_id_p(2)) % cl_t(br_loc_id_p(3)) % indx
      !else
      
      !  clindx(i,j,k) = -1
      !  
      !endif
      
      chi(i,j,k) = gcs_t(i,j,k)%chi
    enddo
  enddo
enddo

$if ($MPI) 
  if(coord == 0) phi(:,:,$lbz) = -BOGUS
$endif

!  Open file which to write global data
write (fname,*) 'phi.out'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif
!  Write binary data for lesgo
$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, form='unformatted', convert='little_endian')
$else
open (1, file=fname, form='unformatted')
$endif
write(1) phi
close (1)

!  Open file which to write global data
write (fname,*) 'brindx.out'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, form='unformatted', convert='little_endian')
$else
open (1, file=fname, form='unformatted')
$endif
write(1) brindx
close (1)

!  Open file which to write global data
write (fname,*) 'clindx.out'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, form='unformatted', convert='little_endian')
$else
open (1, file=fname, form='unformatted')
$endif
write(1) clindx
close (1)

if( filter_chi ) then
write (fname,*) 'chi.out'
fname = trim(adjustl(fname)) 

if(nproc > 1) then
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
endif

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, form='unformatted', convert='little_endian')
$else
open (1, file=fname, form='unformatted')
$endif
write(1) chi
close (1)
endif

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
