program cylinder_skew

implicit none

type cs0
     double precision :: phi
	 double precision, dimension(3) :: xyz
end type cs0

type cs1
	double precision, dimension(3) :: xyz
end type cs1

type vector0
  double precision, dimension(:), pointer :: xyz
end type vector0

type vector1
  double precision, dimension(3) :: xyz
end type vector1

type(cs0), allocatable, target, dimension(:,:,:) :: gcs_t
type(cs1) :: lcs_t
type(cs1) :: tcs_t
type(vector1) :: vtcp_t, vbcp_t

integer, parameter :: Nx=64, Ny=64, Nz=64

double precision, parameter :: skew_angle=30.*3.14/180. !  In radians
double precision, parameter :: crad = 0.1 !  Cylinder radius
double precision, parameter :: clen=0.5 !  Cylinder length
!double precision, parameter, dimension(3) :: axis=(/-cos(90*3.14/180),sin(90.*3.14/180),0/)
double precision, parameter, dimension(3) :: axis=(/1.,0.,0./)
double precision :: tplane,bplane
double precision, parameter :: thresh = 1.e-9
double precision :: echeck,dist,rlcs(3),rrlcs(3)
integer :: lcase,i,j,k,nf

type(vector0) :: vgp_t
type(vector1) :: vlcs_t, vp_t, veck_t
double precision :: eck

double precision, parameter :: xmin=0., xmax=1., dx=(xmax-xmin)/(Nx-1)
double precision, parameter :: ymin=0., ymax=1., dy=(ymax-ymin)/(Ny-1)
double precision, parameter :: zmin=0., zmax=1., dz=(zmax-zmin)/(Nz-1)
double precision, parameter :: a=crad/cos(skew_angle), b=crad

!  Allocate x,y,z for all coordinate systems
allocate(gcs_t(nx,ny,nz))
nullify(vgp_t%xyz)

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

!  Specify global vector to origin of lcs 
vlcs_t%xyz=(/ 0.333, 0.5, 0.25 /)
!  Set the center point of the bottom ellipse
vbcp_t%xyz=vlcs_t%xyz
!  Compute the center point of the top ellipse
call rotation_axis_vector_3d (axis, skew_angle, (/0., 0., clen/),vtcp_t%xyz)
vtcp_t%xyz = vtcp_t%xyz + vbcp_t%xyz

!  Top and bottom plane in gcs
bplane=vlcs_t%xyz(3)
tplane=clen*dcos(skew_angle) + bplane

write(*,*) 'tplane and bplane = ', tplane, bplane

!  Loop over all global coordinates
do k=1,Nz
  do j=1,ny
    do i=1,nx
	  gcs_t(i,j,k)%phi = 1.
	  vgp_t%xyz => gcs_t(i,j,k)%xyz
!  Compute vector to point from lcs in the gcs
      vp_t%xyz = vgp_t%xyz - vlcs_t%xyz
!  Check if below cutting planes
      if(vgp_t%xyz(3) > bplane .and. vgp_t%xyz(3) < tplane) then
        call rotation_axis_vector_3d ( axis, -skew_angle, vp_t%xyz, vp_t%xyz )
		gcs_t(i,j,k)%phi = magnitude_vector_2d(vp_t%xyz(1:2)) - crad
      elseif(vgp_t%xyz(3) > vtcp_t%xyz(3)) then
        veck_t%xyz = vgp_t%xyz - vtcp_t%xyz
        eck = veck_t%xyz(1)**2/a**2 + veck_t%xyz(2)**2/b**2 
        if(eck < 1) then
!  Lies inside of ellipse
          gcs_t(i,j,k)%phi = veck_t%xyz(3)
        endif
      elseif(vgp_t%xyz(3) < vbcp_t%xyz(3)) then
        veck_t%xyz = vgp_t%xyz - vbcp_t%xyz
        eck = veck_t%xyz(1)**2/a**2 + veck_t%xyz(2)**2/b**2 
        if(eck < 1) then
!  Lies inside of ellipse
          gcs_t(i,j,k)%phi = -veck_t%xyz(3)
        endif         
	  endif
    enddo
  enddo
enddo
nf=1
!  Create tecplot formatted velocity field file  
    open (unit = 2,file = 'cylinder_skew.dat', status='unknown',form='formatted', &
      action='write',position='rewind')
    write(2,*) 'variables = "x", "y", "z", "phi"'; 
    write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
      nf,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
    write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE)')//''

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(2,*) gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), gcs_t(i,j,k)%xyz(3), gcs_t(i,j,k)%phi
        enddo
      enddo
    enddo
    close(2)
!  Create plt file the hard way
!    write(fpreplt,*) 'preplot ',ftec,' && rm -v ', ftec
!    call system(fpreplt);


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


