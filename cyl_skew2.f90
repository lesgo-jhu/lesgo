program cylinder_skew

implicit none

type cs0
     double precision :: phi, brindex
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
type(cs0), allocatable, dimension(:,:,:) :: gcs_t
type(cs1) :: lcs_t, lgcs_t, slcs_t, sgcs_t, ecs_t, ebgcs_t, etgcs_t
!  vectors do not have starting point a origin of corresponding
!  coordinate system
type(vector) :: vgcs_t

integer, parameter :: Nx=128, Ny=128, Nz=128
double precision, parameter :: pi = dacos(-1.)
double precision, parameter :: eps = 1.e-6
double precision :: zrot_angle 
double precision, parameter, dimension(3) :: zrot_axis = (/0.,0.,1./)
double precision, parameter :: skew_angle=30.*pi/180. !  In radians
double precision :: crad !  Cylinder radius
double precision :: clen !  Cylinder length
double precision, dimension(3) :: axis

logical :: inside, incir, incyl, inte, inbe, btplanes
double precision :: tplane, bplane
double precision, parameter :: thresh = 1.e-12
double precision :: circk, dist, theta
integer :: i,j,k,nf

double precision :: eck

double precision, parameter :: xmin=0., xmax=1., dx=(xmax-xmin)/(Nx-1)
double precision, parameter :: ymin=0., ymax=1., dy=(ymax-ymin)/(Ny-1)
double precision, parameter :: zmin=0., zmax=1., dz=(zmax-zmin)/(Nz-1)
double precision :: a,b

!  Parameters 
crad = 0.1
clen = 0.25
a=crad/cos(skew_angle)
b=crad

!  Specify global vector to origin of lcs 
lgcs_t%xyz=(/ .333, .5, 0.25 /)
zrot_angle = 0.*pi/180.
axis=(/dcos(zrot_angle+pi/2.),dsin(zrot_angle+pi/2.),0./)
write(*,*) 'skew_angle = ', skew_angle


!  Check if axis has component in z-direction
if(axis(3) .ne. 0.) then
  write(*,*) 'Error: axis cannot have a z component!'
  stop
endif

!  Allocate x,y,z for all coordinate systems
allocate(gcs_t(nx,ny,nz))

!  Create grid in the global coordinate system
do k=1,Nz
  do j=1,ny
    do i=1,nx
      gcs_t(i,j,k)%xyz(1)=(i-1)*dx
      gcs_t(i,j,k)%xyz(2)=(j-1)*dy
      gcs_t(i,j,k)%xyz(3)=(k-1)*dz
    enddo
  enddo
enddo

!  Set the center point of the bottom ellipse
ebgcs_t%xyz=lgcs_t%xyz
!  Compute the center point of the top ellipse
call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen/),etgcs_t%xyz)
etgcs_t%xyz = etgcs_t%xyz + ebgcs_t%xyz

!  Top and bottom plane in gcs
bplane=ebgcs_t%xyz(3)
tplane=etgcs_t%xyz(3)

write(*,*) 'tplane and bplane = ', tplane, bplane

!  Loop over all global coordinates
do k=1,Nz
   
  do j=1,ny
  
    do i=1,nx

!  Initialize the distance function
	  gcs_t(i,j,k)%phi = 1.
      gcs_t(i,j,k)%brindex=0
	
	!  Intialize flags
	  btplanes=.false.
	  incir=.false.
	  incyl=.false.
	  inte=.false.
	  inbe=.false.
	
!!  First check if points are between the top and bottom planes
      if(gcs_t(i,j,k)%xyz(3) > bplane .and. gcs_t(i,j,k)%xyz(3) < tplane) btplanes=.true.

	  !  Compute vector to point from lcs in the gcs
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t%xyz
	  
	  !  Rotate gcs vector into local coordinate system
	  call rotation_axis_vector_3d(axis,-skew_angle,vgcs_t%xyz,lcs_t%xyz)
	  
!  Check if the point lies in the cylinder circle
	  circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
	  if(circk < crad**2) incir = .true.
!  Check if point is in cylinder	  
	  if(btplanes .and. incir) incyl = .true.	 
	  
!  Check if point lies in top ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
	  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
	  eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2 
      if(eck <= 1) inte=.true.
	  
!  Check if point lies in bottom ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
	  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
	  eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2 
      if(eck <= 1) inbe=.true.	  
	  
	  !  Compute the location on the cylinder surface that corresponds
	  !  to the minimum distance.
	  theta = datan2(lcs_t%xyz(2),lcs_t%xyz(1))
	  slcs_t%xyz(1) = crad*dcos(theta)
	  slcs_t%xyz(2) = crad*dsin(theta)
	  slcs_t%xyz(3) = lcs_t%xyz(3)

	  !  Rotate the surface vector in the lcs back into the gcs
	  call rotation_axis_vector_3d(axis,skew_angle,slcs_t%xyz,vgcs_t%xyz)
	  
	  sgcs_t%xyz = vgcs_t%xyz + lgcs_t%xyz !  Vector now corresponds with origin of gcs
	  
!  Check if between cutting planes
      if(sgcs_t%xyz(3) > bplane .and. sgcs_t%xyz(3) < tplane) then
		dist = magnitude_vector_3d(lcs_t%xyz - slcs_t%xyz)		

		if(dabs(dist) < dabs(gcs_t(i,j,k)%phi)) then
		  gcs_t(i,j,k)%phi = dist
		endif
	  else

	    if(sgcs_t%xyz(3) <= bplane .and. .not. inbe) then
		  
  		  vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
		  
		  !  Get vector in ellipse coordinate system		  
 		  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
 		
          call min_dist_to_ellipse(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          dist = dsqrt(dist**2 + ecs_t%xyz(3)**2)

		  if(dist < dabs(gcs_t(i,j,k)%phi)) then
  		    gcs_t(i,j,k)%phi = dist 
            gcs_t(i,j,k)%brindex = 1.
          endif
		
		elseif(sgcs_t%xyz(3) >= tplane .and. .not. inte) then
		  
		  vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
		  
		  		  !  Get vector in ellipse coordinate system		  
 		  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
 		
          call min_dist_to_ellipse(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          dist = dsqrt(dist**2 + ecs_t%xyz(3)**2)

		  if(dist < dabs(gcs_t(i,j,k)%phi)) then
  		    gcs_t(i,j,k)%phi = dist 
            gcs_t(i,j,k)%brindex = 1.
          endif
			
		endif


	
	  endif
!  Check also if the point lies on the ellipses
          
      if(inte) then
	    dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
	    if(dist < dabs(gcs_t(i,j,k)%phi)) then
  		  gcs_t(i,j,k)%phi = dist
		endif
	  endif
		  
      if(inbe) then
		dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane)
	    if(dabs(dist) < dabs(gcs_t(i,j,k)%phi)) then
  		  gcs_t(i,j,k)%phi = dist
		endif
	  endif	
	  
	 if(incyl) then
	   gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
	   gcs_t(i,j,k)%brindex = 1
	 endif

    enddo

  enddo

enddo

! ----------------------------------------------------------
crad = .05
clen = 0.25
!  Specify global vector to origin of lcs 
lgcs_t%xyz=etgcs_t%xyz
zrot_angle = 180.*pi/180.
axis=(/dcos(zrot_angle+pi/2.),dsin(zrot_angle+pi/2.),0./)

!  Set the center point of the bottom ellipse
ebgcs_t%xyz=lgcs_t%xyz
!  Compute the center point of the top ellipse
call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen/),etgcs_t%xyz)
etgcs_t%xyz = etgcs_t%xyz + ebgcs_t%xyz

!  Top and bottom plane in gcs
bplane=ebgcs_t%xyz(3)
tplane=etgcs_t%xyz(3)

write(*,*) 'tplane and bplane = ', tplane, bplane

!  Loop over all global coordinates
do k=1,Nz
   
  do j=1,ny
  
    do i=1,nx

!!  Initialize the distance function
!	  gcs_t(i,j,k)%phi = 1.
!      gcs_t(i,j,k)%brindex=0
	
	!  Intialize flags
	  btplanes=.false.
	  incir=.false.
	  incyl=.false.
	  inte=.false.
	  inbe=.false.
	
!!  First check if points are between the top and bottom planes
      if(gcs_t(i,j,k)%xyz(3) > bplane .and. gcs_t(i,j,k)%xyz(3) < tplane) btplanes=.true.

	  !  Compute vector to point from lcs in the gcs
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t%xyz
	  
	  !  Rotate gcs vector into local coordinate system
	  call rotation_axis_vector_3d(axis,-skew_angle,vgcs_t%xyz,lcs_t%xyz)
	  
!  Check if the point lies in the cylinder circle
	  circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
	  if(circk < crad**2) incir = .true.
!  Check if point is in cylinder	  
	  if(btplanes .and. incir) incyl = .true.	 
	  
!  Check if point lies in top ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
	  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
	  eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2 
      if(eck <= 1) inte=.true.
	  
!  Check if point lies in bottom ellipse
      vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
	  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
	  eck = ecs_t%xyz(1)**2/a**2 + ecs_t%xyz(2)**2/b**2 
      if(eck <= 1) inbe=.true.	  
	  
	  !  Compute the location on the cylinder surface that corresponds
	  !  to the minimum distance.
	  theta = datan2(lcs_t%xyz(2),lcs_t%xyz(1))
	  slcs_t%xyz(1) = crad*dcos(theta)
	  slcs_t%xyz(2) = crad*dsin(theta)
	  slcs_t%xyz(3) = lcs_t%xyz(3)

	  !  Rotate the surface vector in the lcs back into the gcs
	  call rotation_axis_vector_3d(axis,skew_angle,slcs_t%xyz,vgcs_t%xyz)
	  
	  sgcs_t%xyz = vgcs_t%xyz + lgcs_t%xyz !  Vector now corresponds with origin of gcs
	  
!  Check if between cutting planes
      if(sgcs_t%xyz(3) > bplane .and. sgcs_t%xyz(3) < tplane) then
		dist = magnitude_vector_3d(lcs_t%xyz - slcs_t%xyz)		

		if(dabs(dist) < dabs(gcs_t(i,j,k)%phi)) then
		  gcs_t(i,j,k)%phi = dist
		endif
	  else

	    if(sgcs_t%xyz(3) <= bplane .and. .not. inbe) then
		  
  		  vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t%xyz
		  
		  !  Get vector in ellipse coordinate system		  
 		  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
 		
          call min_dist_to_ellipse(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          dist = dsqrt(dist**2 + ecs_t%xyz(3)**2)

		  if(dist < dabs(gcs_t(i,j,k)%phi)) then
  		    gcs_t(i,j,k)%phi = dist 
            gcs_t(i,j,k)%brindex = 1.
          endif
		
		elseif(sgcs_t%xyz(3) >= tplane .and. .not. inte) then
		  
		  vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t%xyz
		  
		  		  !  Get vector in ellipse coordinate system		  
 		  call rotation_axis_vector_3d(zrot_axis, -zrot_angle, vgcs_t%xyz, ecs_t%xyz)
 		
          call min_dist_to_ellipse(a,b,ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

          dist = dsqrt(dist**2 + ecs_t%xyz(3)**2)

		  if(dist < dabs(gcs_t(i,j,k)%phi)) then
  		    gcs_t(i,j,k)%phi = dist 
            gcs_t(i,j,k)%brindex = 1.
          endif
			
		endif


	
	  endif
!  Check also if the point lies on the ellipses
          
      if(inte) then
	    dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane)
	    if(dist < dabs(gcs_t(i,j,k)%phi)) then
  		  gcs_t(i,j,k)%phi = dist
		endif
	  endif
		  
      if(inbe) then
		dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane)
	    if(dabs(dist) < dabs(gcs_t(i,j,k)%phi)) then
  		  gcs_t(i,j,k)%phi = dist
		endif
	  endif	
	  
	 if(incyl) then
	   gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
	   gcs_t(i,j,k)%brindex = 1
	 endif

    enddo

  enddo

enddo

! -----------------------------------------

nf=1

!  Create tecplot formatted velocity field file  
open (unit = 2,file = 'cylinder_skew2.dat', status='unknown',form='formatted', &
  action='write',position='rewind')

write(2,*) 'variables = "x", "y", "z", "phi", "brindex"'; 

write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &

nf,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz

write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''

do k=1,nz
  do j=1,ny
    do i=1,nx
      write(2,*) gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), gcs_t(i,j,k)%xyz(3), gcs_t(i,j,k)%phi, gcs_t(i,j,k)%brindex
    enddo
  enddo
enddo
close(2)
!  Create plt file the hard way
!    write(fpreplt,*) 'preplot ',ftec,' && rm -v ', ftec
!    call system(fpreplt);
!write(*,*) 'gcs_t%phi = ', gcs_t%phi 

stop

contains

double precision function magnitude_vector_3d(vector)
  double precision, dimension(3), intent(IN) :: vector
  magnitude_vector_3d = dsqrt(vector(1)*vector(1) + vector(2)*vector(2) + vector(3)*vector(3))
  return
end function magnitude_vector_3d

double precision function magnitude_vector_2d(vector)
  double precision, dimension(2), intent(IN) :: vector
  magnitude_vector_2d = dsqrt(vector(1)*vector(1) + vector(2)*vector(2))
  return
end function magnitude_vector_2d



end program cylinder_skew


