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
type(cs1), allocatable, dimension(:,:) :: lgcs_t, ebgcs_t, etgcs_t
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
double precision :: skew_angle
double precision, parameter :: thresh = 0.D+00
double precision, dimension(3) :: axis
integer, dimension(3) :: cyl_loc
integer, parameter :: ntrunk = 3
integer, parameter :: ngen = 5
double precision, parameter :: d = 0.5, l = 1.5

logical :: in_cir, in_cyl, in_cyl_top, in_cyl_bottom, above_cyl, below_cyl, in_bottom_surf, btw_planes

double precision, allocatable, dimension(:,:) :: zrot_angle, rad_offset

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
double precision :: z_bottom_surf = 1.*dz

double precision, allocatable, dimension(:,:) :: tplane, bplane

double precision, allocatable, dimension(:,:) :: crad, clen

double precision :: circk, dist, theta

double precision :: a,b

double precision :: eck



end module cylinder_skew_defs

!**************************************************************
program cylinder_skew
!***************************************************************
use cylinder_skew_defs
implicit none


integer :: nt,ng,i,j,k

integer, pointer, dimension(:,:,:) :: brindex

double precision :: atan4
double precision, pointer, dimension(:,:,:) :: phi

! skew_angle=30.*pi/180. !  In radians
! 
! zrot_angle = 0.*pi/180.
! 
! axis=(/dcos(zrot_angle+pi/2.),dsin(zrot_angle+pi/2.),0./)
! 
!  crad = 0.5 !  Cylinder radius
!  clen=1.5/dcos(skew_angle) !  Cylinder length
! a=crad/cos(skew_angle); b=crad
! 
! !  Check if axis has component in z-direction
! if(axis(3) .ne. 0.) then
!   write(*,*) 'Error: axis cannot have a z component!'
!   stop
! endif
! 
! !  Specify global vector to origin of lcs for the cylinder
! lgcs_t%xyz=(/ 2.+0.1*dcos(zrot_angle), 2.+0.1*dsin(zrot_angle), z_bottom_surf /)
! !  Set the center point of the bottom ellipse
! ebgcs_t%xyz=lgcs_t%xyz
! !  Compute the center point of the top ellipse in the gcs
! call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen/),etgcs_t%xyz)
! etgcs_t%xyz = etgcs_t%xyz + ebgcs_t%xyz
! 
! !  Top and bottom z-plane in gcs
! bplane=ebgcs_t%xyz(3)
! tplane=etgcs_t%xyz(3)
! 
! write(*,*) 'tplane and bplane = ', tplane, bplane

call generate_grid()

call initialize() 


do ng=1,ngen
  do nt=1,ntrunk
    call main_loop(nt,ng)
  enddo
enddo

call write_output()



stop

end program cylinder_skew

!**********************************************************************
subroutine pt_loc(nt,ng,i,j,k)
!**********************************************************************

use cylinder_skew_defs

implicit none

integer, intent(IN) :: nt,ng,i,j,k

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
if(gcs_t(i,j,k)%xyz(3) >= bplane(nt,ng) .and. gcs_t(i,j,k)%xyz(3) <= tplane(nt,ng)) then
  btw_planes=.true.
elseif(gcs_t(i,j,k)%xyz(3) > tplane(nt,ng)) then
!  Check if point is below bottom ellipse
  above_cyl = .true.
elseif(gcs_t(i,j,k)%xyz(3) < bplane(nt,ng)) then
!  Check if point is below bottom ellipse
  below_cyl = .true.
endif
      
!  Compute vector to point in the gcs from the lcs 
vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t(nt,ng)%xyz
!  Rotate gcs vector into local coordinate system
call rotation_axis_vector_3d(axis,-skew_angle,vgcs_t%xyz,lcs_t%xyz)
!  Check if the point lies in the cylinder circle
 circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
if(circk <= crad(nt,ng)**2) in_cir = .true.
!  Check if point is in cylinder
if(btw_planes .and. in_cir) in_cyl = .true.

!  Check if point lies in top ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t(nt,ng)%xyz
call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) > bplane(nt,ng)) in_cyl_top=.true. !  Could be below or above

!  Check if point lies in bottom ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(nt,ng)%xyz
call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) < tplane(nt,ng)) in_cyl_bottom=.true. !  Could be below or above

return
end subroutine pt_loc

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

integer :: nt,ng,i,j,k

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

allocate(crad(ntrunk,ngen),clen(ntrunk,ngen))
allocate(lgcs_t(ntrunk,ngen), ebgcs_t(ntrunk,ngen), etgcs_t(ntrunk,ngen))
allocate(zrot_angle(ntrunk,ngen),rad_offset(ntrunk,ngen))
allocate(tplane(ntrunk,ngen),bplane(ntrunk,ngen))

 crad(:,1) = d/2.
 clen(:,1) = l
 skew_angle=30.*pi/180. !  In radians

 rad_offset = 0.1
 
 zrot_angle(1,1) = 0.*pi/180.
 zrot_angle(2,1) = zrot_angle(1,1) + 120.
 zrot_angle(3,1) = zrot_angle(2,1) + 120.

do ng=1,ngen
  do nt=1,ntrunk

    if(ng > 1) then
      zrot_angle(nt,ng) = zrot_angle(nt,ng-1) + 180.
      axis=(/dcos(zrot_angle(nt,ng)+pi/2.),dsin(zrot_angle(nt,ng)+pi/2.),0./)
    endif

    if(ng > 1) then
      crad(nt,ng) = crad(nt,ng-1)/2.
      clen(nt,ng) = clen(nt,ng-1)/2.
      a=crad(nt,ng)/cos(skew_angle); b=crad(nt,ng)

      rad_offset(nt,ng) = rad_offset(nt,ng-1)/2.

    endif
  
!  Check if axis has component in z-direction
    if(axis(3) .ne. 0.) then
      write(*,*) 'Error: axis cannot have a z component!'
      stop
    endif

    if(ng == 1) then
 !  Specify global vector to origin of lcs for the cylinder
      lgcs_t(nt,ng)%xyz=(/ 2.+rad_offset(nt,ng)*dcos(zrot_angle(nt,ng)), 2.+ &
        rad_offset(nt,ng)*dsin(zrot_angle(nt,ng)), z_bottom_surf /)
    else
      lgcs_t(nt,ng)%xyz = etgcs_t(nt,ng-1)%xyz
      lgcs_t(nt,ng)%xyz(1) = lgcs_t(nt,ng)%xyz(1) + rad_offset(nt,ng)*dcos(zrot_angle(nt,ng))
      lgcs_t(nt,ng)%xyz(2) = lgcs_t(nt,ng)%xyz(2) + rad_offset(nt,ng)*dsin(zrot_angle(nt,ng))
    endif

    !  Set the center point of the bottom ellipse
    ebgcs_t(nt,ng)%xyz=lgcs_t(nt,ng)%xyz
    !  Compute the center point of the top ellipse in the gcs
    call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen(nt,ng)/),etgcs_t(nt,ng)%xyz)
    etgcs_t(nt,ng)%xyz = etgcs_t(nt,ng)%xyz + ebgcs_t(nt,ng)%xyz

    !  Top and bottom z-plane in gcs
    bplane(nt,ng)=ebgcs_t(nt,ng)%xyz(3)
    tplane(nt,ng)=etgcs_t(nt,ng)%xyz(3)

    write(*,*) 'tplane(nt,ng) and bplane(nt,ng) = ', tplane(nt,ng), bplane(nt,ng)

  enddo
enddo
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
subroutine assoc_cyl_loc(nt,ng,i,j,k)
!**********************************************************************
use cylinder_skew_defs

implicit none

integer, intent(IN) :: nt,ng,i,j,k
double precision :: atan4

!  Compute theta value on lcs using geometry.atan4
theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

slcs_t%xyz(1) = crad(nt,ng)*dcos(theta)
slcs_t%xyz(2) = crad(nt,ng)*dsin(theta)
slcs_t%xyz(3) = lcs_t%xyz(3)

!  Rotate the surface vector in the lcs back into the gcs
call rotation_axis_vector_3d(axis,skew_angle,slcs_t%xyz,vgcs_t%xyz)

sgcs_t%xyz = vgcs_t%xyz + lgcs_t(nt,ng)%xyz !  Vector now corresponds with origin of gcs

!  Check if point on cylinder surface is between cutting planes
if(sgcs_t%xyz(3) >= bplane(nt,ng) .and. sgcs_t%xyz(3) <= tplane(nt,ng)) then

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
  vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(nt,ng)%xyz

    !  Get vector in ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

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
if(sgcs_t%xyz(3) >= tplane(nt,ng) .and. .not. in_cyl_top) then
  vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t(nt,ng)%xyz

  !  Get vector in ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

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
  dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane(nt,ng))
  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif
endif

return
end subroutine assoc_cyl_loc

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
subroutine main_loop(nt,ng)
!**********************************************************************
use cylinder_skew_defs

implicit none

integer, intent(IN) :: nt,ng
integer :: i,j,k
!  Loop over all global coordinates
do k=1,Nz

  do j=1,ny

    do i=1,nx+2

!  Get all prelimenary information about the point to cylinder location
    if(gcs_t(i,j,k)%phi > 0) then
 
      call pt_loc(nt,ng,i,j,k)

      call assoc_cyl_loc(nt,ng,i,j,k)

      call set_sign(i,j,k)

    endif
    enddo

  enddo

enddo

return
end subroutine main_loop

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

