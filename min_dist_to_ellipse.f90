subroutine min_dist_to_ellipse(a,b,u,v,dist)
!  This subroutine computes the minimum distance to an ellipse with origin
!  x0,y0 = 0,0. 
use error
use root
implicit none

integer :: i
logical :: inside
integer, parameter :: rmax = 4  ! maximum number of roots of ellipse
double precision, intent(IN) :: a,b
double precision, intent(OUT) :: dist
double precision :: u,v
integer :: icount
double precision :: aa,bb,cc,dd,ee,dist_chk,rx,ry,rtmax,eck
double complex, dimension(rmax) :: rt
!double precision, dimension(rmax) :: drx, dry
double precision, parameter :: thresh = 1.e-12

!!  Because ellipse is symmetrical
!xy=dabs(xy)

write(*,*) 'xy = ', xy
!  Initialize distance value
dist = huge(1.)
inside=.false.

!  Check if inside
eck = (xy(1)/a)**2 + (xy(2)/b)**2
if(eck <= 1.) inside=.true.

!  Check if x < thresh
if(dabs(xy(1)) <= thresh) then
  dist = xy(2) - b
elseif(dabs(xy(2)) <= thresh) then
  dist = xy(1) - a
elseif(dabs(xy(1)) <= thresh .and. dabs(xy(2)) <= thresh) then
!  At the origin
  dist = a
else

  aa=2.*(a**2 + b**2)
  bb=b**4 + 4*a**2*b**2 + a**4 - a**2*xy(1)**2 - b**2*xy(2)**2
  cc = 2*(a**2*b**4 + a**4*b**2 - a**2*b**2*xy(1)**2 - a**2*b**2*xy(2)**2)
  dd = a**4*b**4 - a**2*b**4*xy(1)**2 - a**4*b**2*xy(2)**2
  call RootPol(aa,bb,cc,dd,rt(1),rt(2),rt(3),rt(4))

!  Look for real roots
  rtmax=-huge(1.)
  icount=0
  do i=1,rmax
    if(dimag(rt(i)) <= thresh) then
      if(dble(rt(i)) > rtmax) then
  	    rx = xy(1)/(1+rtmax/a**2)
        ry = xy(2)/(1+rtmax/b**2)
	  endif
	  icount = 1
    endif
  enddo

  dist = magnitude_vector_2d((/rx - xy(1),ry - xy(2)/))
  if(inside) dist = -dist

  if(icount == 0) then
    write(*,*) 'Error: No exclusively real roots found!'
    stop
  endif

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
