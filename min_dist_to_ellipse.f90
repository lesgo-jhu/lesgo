subroutine min_dist_to_ellipse(a,b,uin,vin,eps,dist)
!  This subroutine computes the minimum distance to an ellipse with origin
!  x0,y0 = 0,0. 
use error
use root
implicit none

integer :: i
logical :: inside
integer, parameter :: rmax = 4  ! maximum number of roots of ellipse
double precision, intent(IN) :: a,b,eps,uin, vin
double precision, intent(OUT) :: dist
double precision :: u,v,x,y,t,tdiff,func,funcp
integer :: icount
double precision :: aa,bb,cc,dd,ee,dist_chk,rx,ry,rtmax,eck
double complex, dimension(rmax) :: rt
!double precision, dimension(rmax) :: drx, dry
double precision, parameter :: thresh = 1.e-12

!  Because ellipse is symmetrical
u=dabs(uin)
v=dabs(vin)

!  Initialize distance value
dist = huge(1.)
inside=.false.

!  Check if inside
eck = (u/a)**2 + (v/b)**2
if(eck <= 1.) inside=.true.

!  Go through various cases
if(u <= thresh .and. v <= thresh) then
  dist = b
elseif(u <= a - b**2/a .and. v <= thresh) then
  dist = b*dsqrt(1 - u**2/(a**2 - b**2))
elseif(a - b**2/a < u .and. u < a .and. v <= thresh) then
  dist = a - u
elseif(u <= thresh .and. v <= b) then
  dist = b - v
elseif( u > 0. .and. v > 0.) then
!  Employ newtons method

  t = b*v - b**2
  tdiff = huge(1.)
  do while (tdiff > eps)
    func = (a*u / (t + a**2))**2 + (b*v / (t + b**2))**2 - 1
	funcp = -2*( (a*u)**2 / ( t + a**2)**3 + (b*v)**2 / (t + b**2)**3)
	tdiff = -func/funcp
	t = t + tdiff
  enddo	
  
  x = a**2*u/(t + a**2)
  y = b**2*v/(t + b**2)
  
  dist = magnitude_vector_2d((/ u - x, v - y /))
  
endif

if(inside) dist = -dist

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
