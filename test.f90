program test

implicit none

double precision :: x,y,dist,theta
double precision, parameter :: pi = dacos(-1.)
double precision, parameter :: skew_angle=15.*pi/180. !  In radians
double precision, parameter :: crad = 0.1 !  Cylinder radius
double precision, parameter :: a=crad/cos(skew_angle), b=crad

x=-sqrt(3.)/2.
y=-1./2.
write(*,*) 'a, b = ', a,b
write(*,*) 'circle check = ', dsqrt(x**2 + y**2) - crad
call min_dist_to_ellipse(a,b,(/x,y/), dist)

theta = datan2(y,x)
write(*,*) 'theta = ', theta*180./pi

write(*,*) 'dist = ', dist

stop
end program test
