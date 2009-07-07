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

type cs2
    double precision, allocatable, dimension(:,:) :: xyz
end type cs2

type rot
  double precision, allocatable, dimension(:) :: angle
  double precision, allocatable, dimension(:,:) :: axis
end type rot

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
type(cs2), allocatable, dimension(:) :: lgcs_t, ebgcs_t, etgcs_t
type(rot), allocatable, dimension(:) :: zrot_t
!  vectors do not have starting point a origin of corresponding
!  coordinate system
type(vector) :: vgcs_t

integer, parameter :: Nx=64, Ny=64, Nz=64
double precision, parameter :: Lx = 4., dx=Lx/(Nx-1)
double precision, parameter :: Ly = 4., dy=Ly/(Ny-1)
!double precision, parameter :: Lz = 3.587301587301587302, dz = Lz/(Nz-1./2.)
double precision, parameter :: Lz = 4., dz = Lz/(Nz-1./2.)

double precision, parameter :: pi = dacos(-1.)
double precision, parameter :: BOGUS = 1234567890.
double precision, parameter :: iBOGUS = 1234567890
double precision, parameter :: eps = 1.e-12
double precision, parameter, dimension(3) :: zrot_axis = (/0.,0.,1./)
double precision, parameter :: skew_angle = 30.*pi/180.
double precision, parameter :: thresh = 0.D+00

integer, parameter :: ntrunk = 3
integer, parameter :: ngen = 2
double precision, parameter :: d = 0.5, l = 1.
double precision, parameter :: offset = 0.2
double precision, parameter :: scale_fact = 0.5

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
double precision, parameter :: z_bottom_surf = 1.*dz
double precision, dimension(3), parameter :: origin=(/ Lx/2., Ly/2., z_bottom_surf /)

logical :: DEBUG=.true.

logical :: in_cir, in_cyl
logical :: in_cyl_top, in_cyl_bottom
logical :: above_cyl, below_cyl
logical :: in_bottom_surf, btw_planes

integer, dimension(3) :: cyl_loc

integer, allocatable, dimension(:) :: gen_ntrunk
double precision, allocatable, dimension(:) :: crad, clen, rad_offset, a,b, tplane, bplane
double precision, allocatable, dimension(:,:) :: zrot_angle
double precision, allocatable, dimension(:,:,:) :: axis

double precision :: circk, dist, theta
double precision :: eck

end module cylinder_skew_defs

!**************************************************************
program cylinder_skew
!***************************************************************
use cylinder_skew_defs
implicit none

integer :: nt,ng

call generate_grid()

call initialize() 

call main_loop()

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

integer :: ng,nt,i,j,k,istart,iend
double precision :: gen_scale_fact

call allocate_arrays()

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

!  Set cylinder parameters for all generations
do ng=1,ngen
  gen_scale_fact = scale_fact**(ng-1)
  if(DEBUG) write(*,*) 'gen_scale_fact : ', gen_scale_fact
  crad(ng) = gen_scale_fact*d/2.
  clen(ng) = gen_scale_fact*l
  rad_offset = gen_scale_fact*offset
enddo
a = crad/dcos(skew_angle)
b = crad

if(DEBUG) then
  write(*,*) 'ntrunk 	 : ', ntrunk
  write(*,*) 'crad 	 : ', crad
  write(*,*) 'clen 	 : ', clen
  write(*,*) 'rad_offset : ', rad_offset
endif

!  Set rotation angle about z-axis with which the skew angle is applied 
do ng=1,ngen
  zrot_t(ng)%angle(:)=0.
enddo

ng=1 !  Do for the 1st generation (ng = 1)
do nt=1,ntrunk
  zrot_t(ng)%angle(nt) = (360./ntrunk)*(nt-1)*pi/180.
  zrot_t(ng)%axis(:,nt) = (/dcos(zrot_t(ng)%angle(nt)+pi/2.),dsin(zrot_t(ng)%angle(nt)+pi/2.),0./)
  if(DEBUG) then
    write(*,*) 'zrot_t(1)%angle(nt) : ', zrot_t(ng)%angle(nt)*180./pi
    write(*,*) 'zrot_t(1)%axis(:,nt) : ', zrot_t(ng)%axis(:,nt)
  endif
enddo

ng=1 !  Do for the 1st generation (ng = 1)
do nt=1,ntrunk
!  Set the local coordinate system
  lgcs_t(ng)%xyz(:,nt) = origin
  lgcs_t(ng)%xyz(1,nt) = lgcs_t(ng)%xyz(1,nt) + rad_offset(ng)*dcos(zrot_t(ng)%angle(nt))
  lgcs_t(ng)%xyz(2,nt) = lgcs_t(ng)%xyz(2,nt) + rad_offset(ng)*dsin(zrot_t(ng)%angle(nt))

  if(DEBUG) write(*,*) 'lgcs_t(ng)%xyz(:,nt) : ', lgcs_t(ng)%xyz(:,nt)

  !  Set the center point of the bottom ellipse
  ebgcs_t(ng)%xyz(:,nt)=lgcs_t(ng)%xyz(:,nt)
  !  Compute the center point of the top ellipse in the gcs
  call rotation_axis_vector_3d(zrot_t(ng)%axis(:,nt), &
    skew_angle, &
    (/0., 0., clen(ng)/), &
    etgcs_t(ng)%xyz(:,nt))
  etgcs_t(ng)%xyz(:,nt) = etgcs_t(ng)%xyz(:,nt) + ebgcs_t(ng)%xyz(:,nt)
  
  if(DEBUG) then
    write(*,*) 'ebgcs_t(ng)%xyz(:,nt) : ', ebgcs_t(ng)%xyz(:,nt)
    write(*,*) 'etgcs_t(ng)%xyz(:,nt) : ', etgcs_t(ng)%xyz(:,nt)
  endif

enddo

if(ngen > 1) then

  !  Set the lgcs for the new generation
  do ng=2,ngen
  !  Set the rotation angle for each trunk. The first trunk has the same as the first generation
    i=1
    do nt=1,gen_ntrunk(ng)
      zrot_t(ng)%angle(nt) = zrot_t(1)%angle(i)
      if(mod(ng,2)==0) zrot_t(ng)%angle(nt) = zrot_t(ng)%angle(nt) + pi
      zrot_t(ng)%axis(:,nt) = (/dcos(zrot_t(ng)%angle(nt)+pi/2.),dsin(zrot_t(ng)%angle(nt)+pi/2.),0./)
      i = i + 1
      if(i > ntrunk) i = 1
    enddo

  !  Set the local origin in the gcs
  !  Set the lgcs for the new generation
    do j=1,gen_ntrunk(ng - 1)

      istart = (j-1)*ntrunk + 1
      iend   = istart + (ntrunk -1)

      if(DEBUG) then
	write(*,*) 
	write(*,*) 'istart : ', istart
	write(*,*) 'iend   : ', iend
      endif

      do i=istart,iend
	lgcs_t(ng)%xyz(:,i) = etgcs_t(ng-1)%xyz(:,j)
	lgcs_t(ng)%xyz(1,i) = lgcs_t(ng)%xyz(1,i) + rad_offset(ng)*dcos(zrot_t(ng)%angle(i))
	lgcs_t(ng)%xyz(2,i) = lgcs_t(ng)%xyz(2,i) + rad_offset(ng)*dsin(zrot_t(ng)%angle(i))
      enddo

    enddo
  
    !  Set top and bottom of the cylinder
    do nt=1,gen_ntrunk(ng)
      !  Set the center point of the bottom ellipse
      ebgcs_t(ng)%xyz(:,nt)=lgcs_t(ng)%xyz(:,nt)
      !  Compute the center point of the top ellipse in the gcs
      call rotation_axis_vector_3d (zrot_t(ng)%axis(:,nt), &
	skew_angle, &
	(/0., 0., clen(ng)/),&
	etgcs_t(ng)%xyz(:,nt))
      etgcs_t(ng)%xyz(:,nt) = etgcs_t(ng)%xyz(:,nt) + ebgcs_t(ng)%xyz(:,nt)
    enddo

  enddo


endif

!  Top and bottom z-plane in gcs (same for all cylinders in generation)
do ng=1,ngen
  bplane(ng)=ebgcs_t(ng)%xyz(3,1)
  tplane(ng)=etgcs_t(ng)%xyz(3,1)
  write(*,*) 'generation # : ', ng
  write(*,*) 'tplane and bplane = ', bplane(ng), tplane(ng)
enddo

return 

contains

!**********************************************************************
subroutine allocate_arrays()
!**********************************************************************

implicit none

allocate(lgcs_t(ngen), &
  ebgcs_t(ngen), &
  etgcs_t(ngen), &
  zrot_t(ngen))

allocate(a(ngen), &
  b(ngen), &
  crad(ngen),&
  clen(ngen), &
  rad_offset(ngen), &
  tplane(ngen), &
  bplane(ngen), &
  gen_ntrunk(ngen))

do ng=1,ngen
  gen_ntrunk(ng) = ntrunk**ng
enddo

do ng=1,ngen
  allocate(lgcs_t(ng)%xyz(3,gen_ntrunk(ng)))
  allocate(ebgcs_t(ng)%xyz(3,gen_ntrunk(ng)))
  allocate(etgcs_t(ng)%xyz(3,gen_ntrunk(ng)))
  allocate(zrot_t(ng)%angle(gen_ntrunk(ng)))
  allocate(zrot_t(ng)%axis(3,gen_ntrunk(ng)))
enddo

return
end subroutine allocate_arrays

end subroutine initialize


! !**********************************************************************
! subroutine gen_update(ng)
! !**********************************************************************
! 
! use cylinder_skew_defs
! 
! implicit none
! 
! integer, intent(IN) :: ng
! 
! integer :: i,j
! integer :: istart, iend
! integer :: gen_ntrunk_old
! 
! double precision :: alt_angle=pi
! type(cs1), allocatable, dimension(:) :: etgcs_old_t
! 
! 
! if(DEBUG) write(*,*) 'gen_update called. ng : ', ng
! 
! gen_ntrunk_old = ntrunk**(ng - 1)
! allocate(etgcs_old_t(gen_ntrunk_old))
! 
! do j=1,gen_ntrunk_old
!   etgcs_old_t(j)%xyz = etgcs_t(j)%xyz
! enddo
! 
! if(allocated(lgcs_t)) then
!   deallocate(lgcs_t,ebgcs_t,etgcs_t)
!   deallocate(zrot_angle,axis)
! endif
! 
! gen_ntrunk = ntrunk**ng
! 
! allocate(lgcs_t(gen_ntrunk), ebgcs_t(gen_ntrunk), etgcs_t(gen_ntrunk))
! allocate(zrot_angle(gen_ntrunk))
! allocate(axis(3,gen_ntrunk))
! 
! !  Update parameters
!  crad = scale_fact*crad
!  clen = scale_fact*clen
! a=crad/cos(skew_angle); b=crad
! 
! rad_offset = scale_fact*rad_offset
! 
! if(DEBUG) then
!   write(*,*) 'gen_trunk  : ', gen_ntrunk
!   write(*,*) 'crad 	 : ', crad
!   write(*,*) 'clen 	 : ', clen
!   write(*,*) 'rad_offset : ', rad_offset
! endif
! 
! !  Set rotation angle about z-axis with which the skew angle is applied 
! zrot_angle=0.
! do i=1,gen_ntrunk
! 
!   if (mod(ng,2)==0) then
!     zrot_angle(i) =  (360./ntrunk)*(i-1)*pi/180. + alt_angle
!   else
!     zrot_angle(i) =  (360./ntrunk)*(i-1)*pi/180.
!   endif
!   
!   if(DEBUG) write(*,*) 'zrot_angle(i) : ', zrot_angle(i)*180./pi
!   axis(:,i)=(/dcos(zrot_angle(i)+pi/2.),dsin(zrot_angle(i)+pi/2.),0./)
! enddo
! 
! istart=0
! iend=0
! !  Set the lgcs for the new generation
! do j=1,gen_ntrunk_old
! 
!     istart = (j-1)*ntrunk + 1
!     iend   = istart + (ntrunk -1)
! 
!     if(DEBUG) then
!       write(*,*) 
!       write(*,*) 'istart : ', istart
!       write(*,*) 'iend   : ', iend
!     endif
! 
!     do i=istart,iend
!       lgcs_t(i)%xyz = etgcs_old_t(j)%xyz
!       lgcs_t(i)%xyz(1) = lgcs_t(i)%xyz(1) + rad_offset*dcos(zrot_angle(i))
!       lgcs_t(i)%xyz(2) = lgcs_t(i)%xyz(2) + rad_offset*dsin(zrot_angle(i))
!     enddo
! 
! enddo
!   
! !  Set top and bottom of the cylinder
! do i=1,gen_ntrunk
!   !  Set the center point of the bottom ellipse
!   ebgcs_t(i)%xyz=lgcs_t(i)%xyz
!   !  Compute the center point of the top ellipse in the gcs
!   call rotation_axis_vector_3d (axis(:,i), skew_angle, (/0., 0., clen/),etgcs_t(i)%xyz)
!   etgcs_t(i)%xyz = etgcs_t(i)%xyz + ebgcs_t(i)%xyz
! enddo
! 
! !  Top and bottom z-plane in gcs (same for all cylinders in generation)
! bplane=ebgcs_t(1)%xyz(3)
! tplane=etgcs_t(1)%xyz(3)
! write(*,*) 'generation # : ', ng
! write(*,*) 'tplane and bplane = ', tplane, bplane
!   
! return
! end subroutine gen_update


!**********************************************************************
subroutine main_loop()
!**********************************************************************
use cylinder_skew_defs

implicit none

integer :: ng,nt,i,j,k
!  Loop over all global coordinates
do k=1,Nz

  do j=1,ny

    do i=1,nx+2


      do ng=1,ngen
        do nt=1,gen_ntrunk(ng)

	  if(gcs_t(i,j,k)%phi > 0) then

	    call pt_loc(ng,nt,i,j,k)

	    call assoc_cyl_loc(ng,nt,i,j,k)

	    call set_sign(i,j,k)

	  endif

	enddo
      enddo
    enddo

  enddo

enddo

return
end subroutine main_loop


!**********************************************************************
subroutine pt_loc(ng,nt,i,j,k)
!**********************************************************************

use cylinder_skew_defs

implicit none

integer, intent(IN) :: ng,nt,i,j,k

!  Intialize flags
btw_planes=.false.
in_cir=.false.
in_cyl=.false.
in_bottom_surf = .false.
in_cyl_top=.false.
in_cyl_bottom=.false.

!  Also check if point is below bottom surface
if(use_bottom_surf .and. ng == 1) then
  if(gcs_t(i,j,k)%xyz(3) <= z_bottom_surf) in_bottom_surf = .true.
endif

!  First check if points are between the top and bottom planes in the z - gcs
if(gcs_t(i,j,k)%xyz(3) >= bplane(ng) .and. gcs_t(i,j,k)%xyz(3) <= tplane(ng)) then
  btw_planes=.true.
elseif(gcs_t(i,j,k)%xyz(3) > tplane(ng)) then
!  Check if point is below bottom ellipse
  above_cyl = .true.
elseif(gcs_t(i,j,k)%xyz(3) < bplane(ng)) then
!  Check if point is below bottom ellipse
  below_cyl = .true.
else
  write(*,*) 'Error in pt_loc: cannot be anywhere else'
  stop
endif
      
!  Compute vector to point in the gcs from the lcs 
vgcs_t%xyz = gcs_t(i,j,k)%xyz - lgcs_t(ng)%xyz(:,nt)
!  Rotate gcs vector into local coordinate system
call rotation_axis_vector_3d(zrot_t(ng)%axis(:,nt), &
  -skew_angle, &
  vgcs_t%xyz,lcs_t%xyz)

!  Check if the point lies in the cylinder circle
 circk = lcs_t%xyz(1)**2 + lcs_t%xyz(2)**2
if(circk <= crad(ng)**2) in_cir = .true.

!  Check if point is in cylinder
if(btw_planes .and. in_cir) in_cyl = .true.

!  Check if point lies in top ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t(ng)%xyz(:,nt)
call rotation_axis_vector_3d(zrot_axis, &
  -zrot_t(ng)%angle(nt), &
  vgcs_t%xyz, &
  ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a(ng)**2 + ecs_t%xyz(2)**2/b(ng)**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) > bplane(ng)) in_cyl_top=.true. !  Could be below or above

!  Check if point lies in bottom ellipse
vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(ng)%xyz(:,nt)
call rotation_axis_vector_3d(zrot_axis, &
  -zrot_t(ng)%angle(nt), &
  vgcs_t%xyz, &
  ecs_t%xyz)
eck = ecs_t%xyz(1)**2/a(ng)**2 + ecs_t%xyz(2)**2/b(ng)**2
if(eck <= 1 .and. gcs_t(i,j,k)%xyz(3) < tplane(ng)) in_cyl_bottom=.true. !  Could be below or above

return
end subroutine pt_loc

!**********************************************************************
subroutine assoc_cyl_loc(ng,nt,i,j,k)
!**********************************************************************
use cylinder_skew_defs

implicit none

integer, intent(IN) :: ng,nt,i,j,k
double precision :: atan4

!  Compute theta value on lcs using geometry.atan4
theta = atan4(lcs_t%xyz(2),lcs_t%xyz(1))

slcs_t%xyz(1) = crad(ng)*dcos(theta)
slcs_t%xyz(2) = crad(ng)*dsin(theta)
slcs_t%xyz(3) = lcs_t%xyz(3)

!  Rotate the surface vector in the lcs back into the gcs
call rotation_axis_vector_3d(zrot_t(ng)%axis(:,nt), & 
  skew_angle,slcs_t%xyz,vgcs_t%xyz)

sgcs_t%xyz = vgcs_t%xyz + lgcs_t(ng)%xyz(:,nt) !  Vector now corresponds with origin of gcs

!  Check if point on cylinder surface is between cutting planes
if(sgcs_t%xyz(3) >= bplane(ng) .and. sgcs_t%xyz(3) <= tplane(ng)) then

  call vector_magnitude_3d(lcs_t%xyz - slcs_t%xyz,dist)

  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 2
    call set_iset(i,j,k)
  endif
endif
!else
!    if(below_cyl .and. in_cyl_bottom) then
if(use_bottom_surf .and. ng==1 .and. ebgcs_t(ng)%xyz(3,nt) == z_bottom_surf) then
  if(in_cyl_bottom .or. sgcs_t%xyz(3) <= bplane(ng)) then
!  Perform bottom ellipse stuff
  vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(ng)%xyz(:,nt)

    !  Get vector in ellipse coordinate system
    call rotation_axis_vector_3d(zrot_axis, -zrot_t(ng)%angle(nt), vgcs_t%xyz, ecs_t%xyz)

    call ellipse_point_dist_2D_2(a(ng),b(ng),ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

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
elseif(sgcs_t%xyz(3) <= bplane(ng) .and. .not. in_cyl_bottom) then
  vgcs_t%xyz = gcs_t(i,j,k)%xyz - ebgcs_t(ng)%xyz(:,nt)

  !  Get vector in ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -zrot_t(ng)%angle(nt), vgcs_t%xyz, ecs_t%xyz)

  call ellipse_point_dist_2D_2(a(ng),b(ng),ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

  call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif

endif

  !elseif(sgcs_t%xyz(3) >= tplane .and. .not. in_cyl_top) then
if(sgcs_t%xyz(3) >= tplane(ng) .and. .not. in_cyl_top) then

  vgcs_t%xyz = gcs_t(i,j,k)%xyz - etgcs_t(ng)%xyz(:,nt)

  !  Get vector in ellipse coordinate system
  call rotation_axis_vector_3d(zrot_axis, -zrot_t(ng)%angle(nt), vgcs_t%xyz, ecs_t%xyz)

  call ellipse_point_dist_2D_2(a(ng),b(ng),ecs_t%xyz(1),ecs_t%xyz(2),eps, dist)

  call vector_magnitude_2d((/dist, ecs_t%xyz(3) /), dist)

  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif

endif

!  Check also if the point lies on the ellipses
if(in_cyl_top) then
  dist = dabs(gcs_t(i,j,k)%xyz(3) - tplane(ng))
  if(dist < dabs(gcs_t(i,j,k)%phi)) then
    gcs_t(i,j,k)%phi = dist
    gcs_t(i,j,k)%itype = 1
    call set_iset(i,j,k)
  endif
endif

if(ng == 1) then
  if(use_bottom_surf .and. ebgcs_t(ng)%xyz(3,nt) .ne. z_bottom_surf) then
    dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane(ng))
    if(dist < dabs(gcs_t(i,j,k)%phi)) then
      gcs_t(i,j,k)%phi = dist
      gcs_t(i,j,k)%itype = 0
      call set_iset(i,j,k)
    endif
  elseif(.not. use_bottom_surf) then
    dist = dabs(gcs_t(i,j,k)%xyz(3) - bplane(ng))
    if(dist < dabs(gcs_t(i,j,k)%phi)) then
      gcs_t(i,j,k)%phi = dist
      gcs_t(i,j,k)%itype = 0
      call set_iset(i,j,k)
    endif
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

write(2,*) 'variables = "x", "y", "z", "phi", "brindex", "itype"';

write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &

1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz

write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''

do k=1,nz
  do j=1,ny
    do i=1,nx
      write(2,*) gcs_t(i,j,k)%xyz(1), gcs_t(i,j,k)%xyz(2), gcs_t(i,j,k)%xyz(3), gcs_t(i,j,k)%phi, gcs_t(i,j,k)%brindex, gcs_t(i,j,k)%itype
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

