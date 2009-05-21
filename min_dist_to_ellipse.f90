subroutine min_dist_to_ellipse(a,b,xy,dist)
!  This subroutine computes the minimum distance to an ellipse with origin
!  x0,y0 = 0,0. 
use error
use root
implicit none

integer :: i
integer, parameter :: rmax = 4  ! maximum number of roots of ellipse
double precision, intent(IN) :: a,b,xy(2)
double precision, intent(OUT) :: dist
integer :: icount
double precision :: aa,bb,cc,dd,ee,dist_chk,rx,ry,rtmax
double complex, dimension(rmax) :: rt
!double precision, dimension(rmax) :: drx, dry
double precision, parameter :: thresh = 1.e-12

!  Initialize distance value
dist = huge(1.)
!!  Set up coefficients for minimum dist calculations
!aa=b*b*(b*b-a*a)
!bb=2.*a**2*b**2*(b*b-a*a)*xy(1)
!cc=a**4*b**2*xy(1)**2+a**2*b**4*xy(2)**2 - a**2*b**2*(b**2-a**2)
!dd=-2.*a**4*b**2*(b**2-a**2)*xy(1)
!ee=-a**6*b**2*xy(1)**2

!  Check if x < thresh
if(xy(1) <= thresh) then
  dist = dabs(xy(2) - b)
elseif(xy(2) <= thresh) then
  dist = dabs(xy(1) - a)
else
aa=2*b**2 + 2*a**2
bb=b**4 + 4*a**2*b**2 + a**4 - a**2*xy(1)**2 - b**2*xy(2)**2
cc = 2*a**2*b**4 + 2*a**4*b**2 - 2*a**2*b**2*xy(1)**2 - 2*a**2*b**2*xy(2)**2
dd = a**4*b**4 - a**2*b**4*xy(1)**2 - a**4*b**2*xy(2)**2
!call RootPol(bb/aa,cc/aa,dd/aa,ee/aa,rx(1),rx(2),rx(3),rx(4))
call RootPol(aa,bb,cc,dd,rt(1),rt(2),rt(3),rt(4))

rtmax=maxval(cabs(rt))
!write(*,*) 'rt = ', rtmax
!  Compute corresponding y values for the roots
rx = xy(1)/(1+rtmax/a**2)
ry = xy(2)/(1+rtmax/b**2)

!drx = dble(rx)
!dry = dble(ry)
!  Initialize counter for finding exclusively real roots
!icount=0
!do i=1,rmax
!  if(dimag(rx(i)) <= thresh) then
!    icount=icount+1
    !dist_chk = magnitude_vector_2d((/drx(i) - xy(1),dry(i) - xy(2)/))
!!    write(*,*) 'dist_chk = ', dist_chk
!    if(dist_chk < dist) dist = dist_chk
!  endif
!enddo
dist = magnitude_vector_2d((/rx - xy(1),ry - xy(2)/))



!if(icount == 0) then
!  write(*,*) 'Error: No exclusively real roots found!'
!  stop
!endif
endif
return

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
end subroutine min_dist_to_ellipse
