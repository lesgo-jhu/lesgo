!***************************************************************
program cylinder_skew
!***************************************************************

implicit none

type cs0
     integer :: brindex
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

integer, parameter :: Nx=64, Ny=64, Nz=114
double precision, parameter :: Lx = 4., dx=Lx/(Nx-1)
double precision, parameter :: Ly = 4., dy=Ly/(Ny-1)
double precision, parameter :: Lz = 3.587301587301587302, dz = Lz/(Nz-1./2.)

logical :: use_bottom_surf = .true. !  True for making a bottom surface
double precision, parameter :: z_bottom_surf = 1.*dz

double precision, parameter :: pi = dacos(-1.)
double precision, parameter :: BOGUS = 1234567890.
double precision, parameter :: iBOGUS = 1234567890
double precision, parameter :: eps = 1.e-12
double precision, parameter, dimension(3) :: zrot_axis = (/0.,0.,1./)
double precision, parameter :: skew_angle=-30.*pi/180. !  In radians
double precision, parameter :: thresh = 0.D+00


double precision :: zrot_angle
double precision :: crad, clen
double precision, dimension(3) :: axis

logical :: incir, insolid, intop, inbottom, below_top
double precision :: tplane, bplane

double precision :: circk, dist, theta
integer :: i,j,k

double precision :: a,b
double precision :: eck, atan4
integer, pointer, dimension(:,:,:) :: brindex
double precision, pointer, dimension(:,:,:) :: phi



zrot_angle = 0.*pi/180.

axis=(/dcos(zrot_angle+pi/2.),dsin(zrot_angle+pi/2.),0./)

 crad = 0.5 !  Cylinder radius
 clen=1./dcos(skew_angle) !  Cylinder length
a=crad/cos(skew_angle); b=crad

!  Check if axis has component in z-direction
if(axis(3) .ne. 0.) then
  write(*,*) 'Error: axis cannot have a z component!'
  stop
endif

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

!  Initialize the distance function
gcs_t(:,:,:)%phi = BOGUS
!  Set lower level
gcs_t(:,:,0)%phi = -BOGUS
gcs_t(:,:,:)%brindex=iBOGUS

!  Loop over all global coordinates
do k=1,Nz

  do j=1,ny

    do i=1,nx+2

    !  Intialize flags
      below_top=.false.
      incir=.false.
      insolid=.false.
      intop=.false.
      inbottom=.false.

!  Also check if point is below bottom surface
      if(use_bottom_surf) then
        if(gcs_t(i,j,k)%xyz(3) <= z_bottom_surf) insolid = .true.
      endif

!  First check if points are between the top and bottom planes in the z - gcs
      if(gcs_t(i,j,k)%xyz(3) < tplane) below_top=.true.
!  Compute vector to point in the gcs from the lcs 
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t%xyz
!  Rotate gcs vector into local coordinate system
      call rotation_axis_vector_3d(axis,-skew_angle,vgcs_t%xyz,lcs_t%xyz)
!  Check if the point lies in the cylinder circle
      circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
      if(circk <= crad*crad) incir = .true.
!  Check if point is in cylinder
      if(below_top .and. incir) insolid = .true.

!  Check if point lies in top ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
      call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
      eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
      if(eck <= 1) intop=.true.

!  Check if point lies in bottom ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
      call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
      eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
      if(eck <= 1) inbottom=.true.

!  Compute theta value on lcs using geometry.atan4
      theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

      !theta = datan2(lcs_t%xyz(2),lcs_t%xyz(1))
      slcs_t%xyz(1) = crad*dcos(theta)
      slcs_t%xyz(2) = crad*dsin(theta)
      slcs_t%xyz(3) = lcs_t%xyz(3)

      !  Rotate the surface vector in the lcs back into the gcs
      call rotation_axis_vector_3d(axis,skew_angle,slcs_t%xyz,vgcs_t%xyz)

      sgcs_t%xyz = vgcs_t%xyz + lgcs_t%xyz !  Vector now corresponds with origin of gcs

!  Check if between cutting planes
      if(sgcs_t%xyz(3) > bplane .and. sgcs_t%xyz(3) < tplane) then
        call vector_magnitude_3d(lcs_t%xyz - slcs_t%xyz,dist)

        if(dist < dabs(gcs_t(i,j,k)%phi)) gcs_t(i,j,k)%phi = dist

      else

!        if(sgcs_t%xyz(3) <= bplane .and. .not. inbottom) then
!        if(sgcs_t%xyz(3) <= bplane) then

!  Perform bottom ellipse stuff
          vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz

          !  Get vector in ellipse coordinate system
          call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

          call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

          if(dist < dabs(gcs_t(i,j,k)%phi)) gcs_t(i,j,k)%phi = dist

!        elseif(sgcs_t%xyz(3) >= tplane .and. .not. intop) then
        if(sgcs_t%xyz(3) >= tplane .and. .not. intop) then

          vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz

          !  Get vector in ellipse coordinate system
          call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

          call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

          if(dist < dabs(gcs_t(i,j,k)%phi)) gcs_t(i,j,k)%phi = dist

        endif



      endif
!  Check also if the point lies on the ellipses

      if(intop) then
        dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
        if(dist < dabs(gcs_t(i,j,k)%phi)) then
          gcs_t(i,j,k)%phi = dist
        endif
      endif

!       if(inbottom) then
!         dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane)
!         if(dabs(dist) < dabs(gcs_t(i,j,k)%phi)) then
!           gcs_t(i,j,k)%phi = dist
!         endif
!       endif

     if(insolid) then
       gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
       gcs_t(i,j,k)%brindex = -1
     else
       gcs_t(i,j,k)%brindex = 1
     endif

     if(.not. inbottom) then
!  Perform check for creating bottom surface
       dist = dabs(gcs_t(i,j,k)%xyz(3) - z_bottom_surf)
       if(dist < dabs(gcs_t(i,j,k)%phi)) then
         gcs_t(i,j,k)%phi = dist
!  Check the sign of 
         if(gcs_t(i,j,k)%xyz(3) > z_bottom_surf) then
           gcs_t(i,j,k)%brindex = 1
         else
           gcs_t(i,j,k)%brindex = -1
           gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
         endif
       endif
     else
       if(gcs_t(i,j,k)%xyz(3) < z_bottom_surf .and. .not. insolid) then
         gcs_t(i,j,k)%brindex = -1
         gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
       endif
     endif

    enddo

  enddo

enddo

!********************* 2nd Cylinder ***********************************
zrot_angle = 180.*pi/180.

axis=(/dcos(zrot_angle+pi/2.),dsin(zrot_angle+pi/2.),0./)

 crad = 0.25 !  Cylinder radius
 clen=0.5/dcos(skew_angle) !  Cylinder length
a=crad/cos(skew_angle); b=crad

!  Check if axis has component in z-direction
if(axis(3) .ne. 0.) then
  write(*,*) 'Error: axis cannot have a z component!'
  stop
endif

!  Specify global vector to origin of lcs for the cylinder
lgcs_t%xyz=etgcs_t%xyz
!  Set the center point of the bottom ellipse
ebgcs_t%xyz=lgcs_t%xyz
!  Compute the center point of the top ellipse in the gcs
call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen/),etgcs_t%xyz)
etgcs_t%xyz = etgcs_t%xyz + ebgcs_t%xyz

!  Top and bottom z-plane in gcs
bplane=ebgcs_t%xyz(3)
tplane=etgcs_t%xyz(3)

write(*,*) 'tplane and bplane = ', tplane, bplane

!  Loop over all global coordinates
do k=1,Nz

  do j=1,ny

    do i=1,nx+2

if(gcs_t(i,j,k)%brindex > 0) then
    !  Intialize flags
      below_top=.false.
      incir=.false.
      insolid=.false.
      intop=.false.
      inbottom=.false.

!  First check if points are between the top and bottom planes in the z - gcs
!      if(gcs_t(i,j,k)%xyz(3) > bplane .and. gcs_t(i,j,k)%xyz(3) < tplane) below_top=.true.
      if(gcs_t(i,j,k)%xyz(3) < tplane) below_top=.true.

      !  Compute vector to point in the gcs from the lcs 
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t%xyz

      !  Rotate gcs vector into local coordinate system
      call rotation_axis_vector_3d(axis,-skew_angle,vgcs_t%xyz,lcs_t%xyz)

!  Check if the point lies in the cylinder circle
      circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
      if(circk <= crad*crad) incir = .true.
!  Check if point is in cylinder
      if(below_top .and. incir) insolid = .true.

!  Check if point lies in top ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
      call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
      eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
      if(eck <= 1) intop=.true.

!  Check if point lies in bottom ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
      call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, &
        ecs_t%xyz)
      eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2
      if(eck <= 1) inbottom=.true.

!  Compute theta value on lcs using geometry.atan4
      theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

      !theta = datan2(lcs_t%xyz(2),lcs_t%xyz(1))
      slcs_t%xyz(1) = crad*dcos(theta)
      slcs_t%xyz(2) = crad*dsin(theta)
      slcs_t%xyz(3) = lcs_t%xyz(3)

      !  Rotate the surface vector in the lcs back into the gcs
      call rotation_axis_vector_3d(axis,skew_angle,slcs_t%xyz,vgcs_t%xyz)

      sgcs_t%xyz = vgcs_t%xyz + lgcs_t%xyz !  Vector now corresponds with origin of gcs

!  Check if between cutting planes
      if(sgcs_t%xyz(3) > bplane .and. sgcs_t%xyz(3) < tplane) then
        call vector_magnitude_3d(lcs_t%xyz - slcs_t%xyz,dist)

        if(dist < dabs(gcs_t(i,j,k)%phi)) gcs_t(i,j,k)%phi = dist

      else

!        if(sgcs_t%xyz(3) <= bplane .and. .not. inbottom) then
!        if(sgcs_t%xyz(3) <= bplane) then

!  Perform bottom ellipse stuff
          vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz

          !  Get vector in ellipse coordinate system
          call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

          call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

          if(dist < dabs(gcs_t(i,j,k)%phi)) gcs_t(i,j,k)%phi = dist

!        elseif(sgcs_t%xyz(3) >= tplane .and. .not. intop) then
        if(sgcs_t%xyz(3) >= tplane .and. .not. intop) then

          vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz

          !  Get vector in ellipse coordinate system
          call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)

          call ellipse_point_dist_2D_2(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

          if(dist < dabs(gcs_t(i,j,k)%phi)) gcs_t(i,j,k)%phi = dist

        endif



      endif
!  Check also if the point lies on the ellipses

      if(intop) then
        dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
        if(dist < dabs(gcs_t(i,j,k)%phi)) then
          gcs_t(i,j,k)%phi = dist
        endif
      endif

!       if(inbottom) then
!         dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane)
!         if(dabs(dist) < dabs(gcs_t(i,j,k)%phi)) then
!           gcs_t(i,j,k)%phi = dist
!         endif
!       endif

     if(insolid) then
       gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
       gcs_t(i,j,k)%brindex = -1
     else
       gcs_t(i,j,k)%brindex = 1
     endif

     if(.not. inbottom) then
!  Perform check for creating bottom surface
       dist = dabs(gcs_t(i,j,k)%xyz(3) - z_bottom_surf)
       if(dist < dabs(gcs_t(i,j,k)%phi)) then
         gcs_t(i,j,k)%phi = dist
!  Check the sign of 
         if(gcs_t(i,j,k)%xyz(3) > z_bottom_surf) then
           gcs_t(i,j,k)%brindex = 1
         else
           gcs_t(i,j,k)%brindex = -1
           gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
         endif
       endif
     else
       if(gcs_t(i,j,k)%xyz(3) < z_bottom_surf .and. .not. insolid) then
         gcs_t(i,j,k)%brindex = -1
         gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
       endif
     endif
endif
    enddo

  enddo

enddo

!********************* 2nd Cylinder ***********************************


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
! do k=1,nz
!   do j=1,ny
!     do i=1,nx+2
!       write (1) gcs_t(i, j, k)%phi
!     enddo
!   enddo
! enddo
!write (1) gcs_t(:, :, 1:nz)%phi
write(1) phi
close (1)

open (1, file='brindex.out', form='unformatted')
! do k=1,nz
!   do j=1,ny
!     do i=1,nx+2
!       write (1) gcs_t(i, j, k)%brindex
!     enddo
!   enddo
! enddo

!write (1) gcs_t(:, :, 1:nz)%brindex
write(1) brindex
! close (1)
stop

end program cylinder_skew


