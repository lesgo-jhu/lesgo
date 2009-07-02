!**********************************************************************
module cylinder_skew_defs
!**********************************************************************

implicit none

save
public

logical :: in_cir, in_cyl, in_cyl_top, in_cyl_bottom, above_cyl, below_cyl, in_bottom_surf, btw_planes

type cs0
     integer :: brindex, iset
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
type(cs1) :: lcs_t, lgcs_t, slcs_t, sgcs_t, ecs_t, ebgcs_t, etgcs_t
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
double precision, parameter :: skew_angle=60.*pi/180. !  In radians
double precision, parameter :: thresh = 0.D+00
double precision, dimension(3) :: axis
integer, dimension(3) :: cyl_loc


double precision :: zrot_angle

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
double precision :: z_bottom_surf = 10.*dz

double precision :: tplane, bplane

double precision :: crad, clen

double precision :: circk, dist, theta

double precision :: a,b

double precision :: eck



end module cylinder_skew_defs

!**************************************************************
program cylinder_skew
!***************************************************************
use cylinder_skew_defs
implicit none


integer :: i,j,k

integer, pointer, dimension(:,:,:) :: brindex

double precision :: atan4
double precision, pointer, dimension(:,:,:) :: phi

zrot_angle = 180.*pi/180.

axis=(/dcos(zrot_angle+pi/2.),dsin(zrot_angle+pi/2.),0./)

 crad = 0.5 !  Cylinder radius
 clen=3./dcos(skew_angle) !  Cylinder length
a=crad/cos(skew_angle); b=crad

!  Check if axis has component in z-direction
if(axis(3) .ne. 0.) then
  write(*,*) 'Error: axis cannot have a z component!'
  stop
endif

!  Specify global vector to origin of lcs for the cylinder
lgcs_t%xyz=(/ 2.825, 2., z_bottom_surf /)
!  Set the center point of the bottom ellipse
ebgcs_t%xyz=lgcs_t%xyz
!  Compute the center point of the top ellipse in the gcs
call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen/),etgcs_t%xyz)
etgcs_t%xyz = etgcs_t%xyz + ebgcs_t%xyz

!  Top and bottom z-plane in gcs
bplane=ebgcs_t%xyz(3)
tplane=etgcs_t%xyz(3)

write(*,*) 'tplane and bplane = ', tplane, bplane

call generate_grid()
call initialize() 

!  Loop over all global coordinates
do k=1,Nz

  do j=1,ny

    do i=1,nx+2

    !  Intialize flags
      btw_planes=.false.
      in_cir=.false.
      in_cyl=.false.
      in_bottom_surf = .false.
      in_cyl_top=.false.
      in_cyl_bottom=.false.

!  Get all prelimenary information about the point to cylinder location
      call pt_to_cyl_loc(i,j,k)

      call assoc_cyl_loc(i,j,k)


!  Check also if the point lies on the ellipses
      if(in_cyl_top) then
        dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
        if(dist < dabs(gcs_t(i,j,k)%phi)) then
          gcs_t(i,j,k)%phi = dist
          call set_iset(i,j,k)
        endif
      endif

     if(in_cyl .or. in_bottom_surf) then
       gcs_t(i,j,k)%phi = -dabs(gcs_t(i,j,k)%phi)
       gcs_t(i,j,k)%brindex = -1
     else     
       gcs_t(i,j,k)%brindex = 1
     endif

!      if(.not. in_cyl .and. .not. in_cyl_bottom) then
! !  Perform check for creating bottom surface
!        dist = dabs(gcs_t(i,j,k)%xyz(3) - z_bottom_surf)
!        if(dist < dabs(gcs_t(i,j,k)%phi)) then
!          gcs_t(i,j,k)%phi = dist
!          call set_iset(i,j,k)
! !  Assign appropriate sign
!          if(in_bottom_surf) then
!            gcs_t(i,j,k)%phi = -dabs(gcs_t(i,j,k)%phi)
!            gcs_t(i,j,k)%brindex = -1
!          else
!            gcs_t(i,j,k)%brindex = 1
!          endif
!        endif
!      endif

    enddo

  enddo

enddo

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

stop

end program cylinder_skew

!**********************************************************************
subroutine pt_to_cyl_loc(i,j,k)
!**********************************************************************

use cylinder_skew_defs

implicit none

integer, intent(IN) :: i,j,k

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
vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t%xyz
!  Rotate gcs vector into local coordinate system
call rotation_axis_vector_3d(axis,-skew_angle,vgcs_t%xyz,lcs_t%xyz)
!  Check if the point lies in the cylinder circle
 circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
if(circk <= crad*crad) in_cir = .true.
!  Check if point is in cylinder
if(btw_planes .and. in_cir) in_cyl = .true.

!  Check if point lies in top ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) > bplane) in_cyl_top=.true. !  Could be below or above

!  Check if point lies in bottom ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) < tplane) in_cyl_bottom=.true. !  Could be below or above

return
end subroutine pt_to_cyl_loc

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

if(use_bottom_surf) then
!  Loop over all global coordinates
  do k=1,Nz
    gcs_t(:,:,k)%phi = gcs_t(:,:,k)%xyz(3) - z_bottom_surf
    if(gcs_t(1,1,k)%phi <= 0.) gcs_t(:,:,k)%brindex = -1
  enddo
endif

return 
end subroutine initialize

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
subroutine assoc_cyl_loc(i,j,k)
!**********************************************************************
use cylinder_skew_defs

implicit none

integer, intent(IN) :: i,j,k
double precision :: atan4

!  Compute theta value on lcs using geometry.atan4
theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

slcs_t%xyz(1) = crad*dcos(theta)
slcs_t%xyz(2) = crad*dsin(theta)
slcs_t%xyz(3) = lcs_t%xyz(3)

!  Rotate the surface vector in the lcs back into the gcs
call rotation_axis_vector_3d(axis,skew_angle,slcs_t%xyz,vgcs_t%xyz)

sgcs_t%xyz = vgcs_t%xyz + lgcs_t%xyz !  Vector now corresponds with origin of gcs

!  Check if point on cylinder surface is between cutting planes
if(sgcs_t%xyz(3) >= bplane .and. sgcs_t%xyz(3) <= tplane) then

  call vector_magnitude_3d(lcs_t%xyz - slcs_t%xyz,dist)

  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    call set_iset(i,j,k)
  endif

else
    if(below_cyl .and. in_cyl_bottom) then
!  Perform bottom ellipse stuff
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz

    !  Get vector in ellipse coordinate system
    call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

    call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

    call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

    if(dist < dabs(gcs_t(i,j,k)%phi)) then
      gcs_t(i,j,k)%phi = dist
      call set_iset(i,j,k)
    endif

  elseif(sgcs_t%xyz(3) >= tplane .and. .not. in_cyl_top) then

    vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz

    !  Get vector in ellipse coordinate system
    call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

    call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

    call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

    if(dist < dabs(gcs_t(i,j,k)%phi)) then
      gcs_t(i,j,k)%phi = dist
      call set_iset(i,j,k)
    endif

  endif

endif


return
end subroutine assoc_cyl_loc