program test

implicit none

integer, parameter :: Nx=32, Ny=32
double precision 
double precision :: x,y,dist,theta
double precision, parameter :: pi = dacos(-1.)
!double precision, parameter :: skew_angle=30.*pi/180. !  In radians
!double precision, parameter :: crad = 0.1 !  Cylinder radius
double precision, parameter :: a=2., b=1.

!x=-sqrt(3.)/2.
!y=-1./2.
x=2.*crad
y=0.
write(*,*) 'a, b = ', a,b
write(*,*) 'circle check = ', dsqrt(x**2 + y**2) - crad
call min_dist_to_ellipse(a,b,(/x,y/), dist)

theta = datan2(y,x)
write(*,*) 'theta = ', theta*180./pi

write(*,*) 'dist = ', dist

stop
end program test
