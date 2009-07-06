!**********************************************************************
module cylinder_skew_defs
!**********************************************************************

implicit none

save
public

type cs0
     integer :: brindex, iset, itype
     double precision :: phi
     double precision, dimension(3) :: xyz
end type cs0

type cs1
    double precision, dimension(3) :: xyz
end type cs1

!type vector0
!  double precision, dimension(:), pointer :: xyz
!end type vector0

type vector
  double precision, dimension(3) :: xyz
end type vector

!  cs{0,1} all correspond to vectors with the origin at the
!  corresponding coordinate system
type(cs0), target, allocatable, dimension(:,:,:) :: gcs_t
type(cs1) :: lcs_t, slcs_t, sgcs_t, ecs_t
type(cs1), allocatable, dimension(:) :: lgcs_t, ebgcs_t, etgcs_t
!  vectors do not have starting point a origin of corresponding
!  coordinate system
type(vector) :: vgcs_t

integer, parameter :: Nx=64, Ny=64, Nz=57
double precision, parameter :: Lx = 4., dx=Lx/(Nx-1)
double precision, parameter :: Ly = 4., dy=Ly/(Ny-1)
double precision, parameter :: Lz = 3.587301587301587302, dz = Lz/(Nz-1./2.)

double precision, parameter :: pi = dacos(-1.)
double precision, parameter :: BOGUS = 1234567890.
double precision, parameter :: iBOGUS = 1234567890
double precision, parameter :: eps = 1.e-12
double precision, parameter, dimension(3) :: zrot_axis = (/0.,0.,1./)
double precision, parameter :: skew_angle = 30.*pi/180.
double precision, parameter :: thresh = 0.D+00

integer, parameter :: ntrunk = 1
integer, parameter :: ngen = 1
double precision, parameter :: d = 1., l = 1.5
double precision, parameter :: offset = 0.1
double precision, parameter :: scale_fact = 0.5

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
double precision, parameter :: z_bottom_surf = 1.*dz
double precision, dimension(3), parameter :: origin=(/ Lx/2., Ly/2., z_bottom_surf /)

logical :: DEBUG=.true.

logical :: in_cir, in_cyl
logical :: in_cyl_top, in_cyl_bottom
logical :: above_cyl, below_cyl
logical :: in_bottom_surf, btw_planes

integer :: gen_ntrunk
integer, dimension(3) :: cyl_loc

double precision, allocatable, dimension(:) :: zrot_angle
double precision, allocatable, dimension(:,:) :: axis

double precision :: crad, clen, rad_offset
double precision :: circk, dist, theta
double precision :: a,b
double precision :: eck
double precision :: tplane, bplane

end module cylinder_skew_defs

!**************************************************************
program cylinder_skew
!***************************************************************
use cylinder_skew_defs
implicit none

integer :: nt,ng

call generate_grid()
call initialize() 

do ng=1,ngen
  if(ng > 1) call gen_update(ng)

  do nt=1,gen_ntrunk
    call main_loop(nt)
  enddo

enddo

call write_output()

stop

end program cylinder_skew

!**********************************************************************
subroutine generate_grid()
!**********************************************************************
use cylinder_skew_defs

implicit none

integer :: i,j,k

!  Allocate x,y,z for all coordinate systems
allocate(gcs_t(nx+2,ny,0:nz))

!  Create grid in the global coordinate system
do k=0,Nz
  do j=1,ny
    do i=1,nx+2
      gcs_t(i,j,k)%xyz(1)=(i-1)*dx 
      gcs_t(i,j,k)%xyz(2)=(j-1)*dy 
      gcs_t(i,j,k)%xyz(3)=(k-1./2.)*dz 
    enddo
  enddo
enddo

return
end subroutine generate_grid

!**********************************************************************
subroutine initialize()
!**********************************************************************
use cylinder_skew_defs

implicit none

integer :: i,j,k

!  Initialize the distance function
gcs_t(:,:,:)%phi = BOGUS
!  Set lower level
gcs_t(:,:,0)%phi = -BOGUS
gcs_t(:,:,:)%brindex=1

!  Initialize the iset flag
gcs_t(:,:,:)%iset=0

!  Initialize the point to surface association
gcs_t(:,:,:)%itype=-1 !  0 - bottom, 1 - top, 2 - side

if(use_bottom_surf) then
  gcs_t(:,:,:)%itype=-1
!  Loop over all global coordinates
  do k=1,Nz
    gcs_t(:,:,k)%phi = gcs_t(:,:,k)%xyz(3) - z_bottom_surf
    if(gcs_t(1,1,k)%phi <= 0.) gcs_t(:,:,k)%brindex = -1
  enddo
endif

!  Set initial generation
allocate(lgcs_t(ntrunk), ebgcs_t(ntrunk), etgcs_t(ntrunk))
allocate(zrot_angle(ntrunk))
allocate(axis(3,ntrunk))

!  Update parameters
crad = d/2
clen = l
rad_offset = offset

gen_ntrunk = ntrunk

if(DEBUG) then
  write(*,*) 'ntrunk 	 : ', ntrunk
  write(*,*) 'crad 	 : ', crad
  write(*,*) 'clen 	 : ', clen
  write(*,*) 'rad_offset : ', rad_offset
endif

!  Set rotation angle about z-axis with which the skew angle is applied 
zrot_angle=0.
do i=1,ntrunk
  zrot_angle(i) =  (360./ntrunk)*(i-1)*pi/180.
  if(DEBUG) write(*,*) 'zrot_angle(i) : ', zrot_angle(i)*180./pi
  axis(:,i)=(/dcos(zrot_angle(i)+pi/2.),dsin(zrot_angle(i)+pi/2.),0./)
enddo

do i=1,ntrunk
  lgcs_t(i)%xyz = origin
  lgcs_t(i)%xyz(1) = lgcs_t(i)%xyz(1) + rad_offset*dcos(zrot_angle(i))
  lgcs_t(i)%xyz(2) = lgcs_t(i)%xyz(2) + rad_offset*dsin(zrot_angle(i))

  if(DEBUG) write(*,*) 'lgcs_t(i)%xyz : ', lgcs_t(i)%xyz

enddo

!  Set top and bottom of the cylinder
do i=1,ntrunk
  !  Set the center point of the bottom ellipse
  ebgcs_t(i)%xyz=lgcs_t(i)%xyz
  !  Compute the center point of the top ellipse in the gcs
  call rotation_axis_vector_3d (axis(:,i), skew_angle, (/0., 0., clen/),etgcs_t(i)%xyz)
  etgcs_t(i)%xyz = etgcs_t(i)%xyz + ebgcs_t(i)%xyz
  
  if(DEBUG) then
    write(*,*) 'etgcs_t(i)%xyz : ', etgcs_t(i)%xyz
    write(*,*) 'ebgcs_t(i)%xyz : ', ebgcs_t(i)%xyz
  endif

enddo 

!  Top and bottom z-plane in gcs (same for all cylinders in generation)
bplane=ebgcs_t(1)%xyz(3)
tplane=etgcs_t(1)%xyz(3)
write(*,*) 'generation # : ', 1
write(*,*) 'tplane and bplane = ', tplane, bplane

return 
end subroutine initialize


!**********************************************************************
subroutine gen_update(ng)
!**********************************************************************

use cylinder_skew_defs

implicit none

integer, intent(IN) :: ng

integer :: i,j
integer :: istart, iend
integer :: gen_ntrunk_old

type(cs1), allocatable, dimension(:) :: etgcs_old_t


if(DEBUG) write(*,*) 'gen_update called. ng : ', ng

gen_ntrunk_old = ntrunk**(ng - 1)
allocate(etgcs_old_t(gen_ntrunk_old))

do j=1,gen_ntrunk_old
  etgcs_old_t(j)%xyz = etgcs_t(j)%xyz
enddo

if(allocated(lgcs_t)) then
  deallocate(lgcs_t,ebgcs_t,etgcs_t)
  deallocate(zrot_angle,axis)
endif

gen_ntrunk = ntrunk**ng

allocate(lgcs_t(gen_ntrunk), ebgcs_t(gen_ntrunk), etgcs_t(gen_ntrunk))
allocate(zrot_angle(gen_ntrunk))
allocate(axis(3,gen_ntrunk))

!  Update parameters
crad = (scale_fact**ng)*d/2
clen = (scale_fact**ng)*l
rad_offset = (scale_fact**ng)*offset

!  Set rotation angle about z-axis with which the skew angle is applied 
zrot_angle=0.
do i=1,ntrunk
  zrot_angle(i) =  (360./ntrunk)*(i-1)*pi/180.
  if(DEBUG) write(*,*) 'zrot_angle(i) : ', zrot_angle(i)*180./pi
  axis(:,i)=(/dcos(zrot_angle(i)+pi/2.),dsin(zrot_angle(i)+pi/2.),0./)
enddo

istart=0
iend=0
!  Set the lgcs for the new generation
do j=1,gen_ntrunk_old

    istart = (j-1)*ntrunk + 1
    iend   = istart + (ntrunk -1)

    do i=istart,iend
      lgcs_t(i)%xyz = etgcs_old_t(j)%xyz
      lgcs_t(i)%xyz(1) = lgcs_t(i)%xyz(1) + rad_offset*dcos(zrot_angle(i))
      lgcs_t(i)%xyz(2) = lgcs_t(i)%xyz(2) + rad_offset*dsin(zrot_angle(i))
    enddo

enddo
  
!  Set top and bottom of the cylinder
do i=1,gen_ntrunk
  !  Set the center point of the bottom ellipse
  ebgcs_t(i)%xyz=lgcs_t(i)%xyz
  !  Compute the center point of the top ellipse in the gcs
  call rotation_axis_vector_3d (axis(:,i), skew_angle, (/0., 0., clen/),etgcs_t(i)%xyz)
  etgcs_t(i)%xyz = etgcs_t(i)%xyz + ebgcs_t(i)%xyz
enddo

!  Top and bottom z-plane in gcs (same for all cylinders in generation)
bplane=ebgcs_t(1)%xyz(3)
tplane=etgcs_t(1)%xyz(3)
write(*,*) 'generation # : ', ng
write(*,*) 'tplane and bplane = ', tplane, bplane
  
return
end subroutine gen_update


!**********************************************************************
subroutine main_loop(nt)
!**********************************************************************
use cylinder_skew_defs

implicit none

integer, intent(IN) :: nt
integer :: i,j,k
!  Loop over all global coordinates
do k=1,Nz

  do j=1,ny

    do i=1,nx+2

!  Get all prelimenary information about the point to cylinder location
    if(gcs_t(i,j,k)%phi > 0) then

      call pt_loc(nt,i,j,k)

      call assoc_cyl_loc(nt,i,j,k)

      call set_sign(i,j,k)

    endif
    enddo

  enddo

enddo

return
end subroutine main_loop


!**********************************************************************
subroutine pt_loc(nt,i,j,k)
!**********************************************************************

use cylinder_skew_defs

implicit none

integer, intent(IN) :: nt,i,j,k

!  Intialize flags
btw_planes=.false.
in_cir=.false.
in_cyl=.false.
in_bottom_surf = .false.
in_cyl_top=.false.
in_cyl_bottom=.false.

!  Also check if point is below bottom surface
if(use_bottom_surf) then
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
endif
      
!  Compute vector to point in the gcs from the lcs 
vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t(nt)%xyz
!  Rotate gcs vector into local coordinate system
call rotation_axis_vector_3d(axis(:,nt),-skew_angle,vgcs_t%xyz,lcs_t%xyz)
!  Check if the point lies in the cylinder circle
 circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
if(circk <= crad**2) in_cir = .true.
!  Check if point is in cylinder
if(btw_planes .and. in_cir) in_cyl = .true.

!  Check if point lies in top ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t(nt)%xyz
call rotation_axis_vector_3d(zrot_axis, -zrot_angle(nt), vgcs_t%xyz, &
        ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) > bplane) in_cyl_top=.true. !  Could be below or above

!  Check if point lies in bottom ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(nt)%xyz
call rotation_axis_vector_3d(zrot_axis, -zrot_angle(nt), vgcs_t%xyz, &
        ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) < tplane) in_cyl_bottom=.true. !  Could be below or above

return
end subroutine pt_loc

!**********************************************************************
subroutine assoc_cyl_loc(nt,i,j,k)
!**********************************************************************
use cylinder_skew_defs

implicit none

integer, intent(IN) :: nt,i,j,k
double precision :: atan4

!  Compute theta value on lcs using geometry.atan4
theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

slcs_t%xyz(1) = crad*dcos(theta)
slcs_t%xyz(2) = crad*dsin(theta)
slcs_t%xyz(3) = lcs_t%xyz(3)

!  Rotate the surface vector in the lcs back into the gcs
call rotation_axis_vector_3d(axis(3,nt),skew_angle,slcs_t%xyz,vgcs_t%xyz)

sgcs_t%xyz = vgcs_t%xyz + lgcs_t(nt)%xyz !  Vector now corresponds with origin of gcs

!  Check if point on cylinder surface is between cutting planes
if(sgcs_t%xyz(3) >= bplane .and. sgcs_t%xyz(3) <= tplane) then

  call vector_magnitude_3d(lcs_t%xyz - slcs_t%xyz,dist)

  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 2
    call set_iset(i,j,k)
  endif
endif
!else
!    if(below_cyl .and. in_cyl_bottom) then

if(in_cyl_bottom) then
!  Perform bottom ellipse stuff
  vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(nt)%xyz

    !  Get vector in ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -zrot_angle(nt), vgcs_t%xyz, ecs_t%xyz)

  call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

  call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

!   if(gcs_t(i,j,k)%itype == -1 .and. sgcs_t%xyz(3) <= bplane) then
!     gcs_t(i,j,k)%phi = dist
!     gcs_t(i,j,k)%itype = 0
!     call set_iset(i,j,k)      
!   elseif(dist < dabs(gcs_t(i,j,k)%phi)) then
  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 0
    call set_iset(i,j,k)
  endif
endif

  !elseif(sgcs_t%xyz(3) >= tplane .and. .not. in_cyl_top) then
if(sgcs_t%xyz(3) >= tplane .and. .not. in_cyl_top) then
  vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t(nt)%xyz

  !  Get vector in ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -zrot_angle(nt), vgcs_t%xyz, ecs_t%xyz)

  call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

  call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif

endif

!  Check also if the point lies on the ellipses
if(in_cyl_top) then
  dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif
endif

return
end subroutine assoc_cyl_loc

!**********************************************************************
subroutine set_iset(i,j,k)
!**********************************************************************
use cylinder_skew_defs, only : gcs_t

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
use cylinder_skew_defs

implicit none

integer, intent(IN) :: i,j,k
!if(gcs_t(i,j,k)%phi > 0) then
  if(in_cyl .or. in_bottom_surf) then
    gcs_t(i,j,k)%phi = -dabs(gcs_t(i,j,k)%phi)
    gcs_t(i,j,k)%brindex = -1
  else    
    gcs_t(i,j,k)%brindex = 1
  endif
!endif
return
end subroutine set_sign

!**********************************************************************
subroutine write_output()
!**********************************************************************
use cylinder_skew_defs

implicit none

integer :: i,j,k

integer, pointer, dimension(:,:,:) :: brindex
double precision, pointer, dimension(:,:,:) :: phi

!  Create tecplot formatted phi and brindex field file
open (unit = 2,file = 'cylinder_skew.dat', status='unknown',form='formatted', &
  action='write',position='rewind')

write(2,*) 'variables = "x", "y", "z", "phi", "brindex"';

write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &

1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz

write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''

do k=1,nz
  do j=1,ny
    do i=1,nx
      write(2,*) gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), gcs_t(i,j,k)%xyz(3), gcs_t(i,j,k)%phi, gcs_t(i,j,k)%brindex
    enddo
  enddo
enddo
close(2)

nullify(phi,brindex)
allocate(phi(nx+2,ny,1:nz))
allocate(brindex(nx+2,ny,1:nz))
do k=1,nz
  do j = 1,ny
    do i = 1,nx+2
      phi(i,j,k) = gcs_t(i,j,k)%phi
      brindex(i,j,k) = gcs_t(i,j,k)%brindex
    enddo
  enddo
enddo

!  Write binary data for lesgo
open (1, file='phi.out', form='unformatted')
write(1) phi
close (1)

open (1, file='brindex.out', form='unformatted')
write(1) brindex
close (1)

return
end subroutine write_output

