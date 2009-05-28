program test

implicit none

integer :: i,j
integer, parameter :: Nx=256, Ny=256
double precision, parameter :: xmin = -4, xmax=4
double precision, parameter :: ymin = -4, ymax=4
double precision :: x(Nx), y(Ny), dx, dy
double precision :: phi(Nx,Ny),inside(Nx,Ny)
double precision :: eck
double precision, parameter :: pi = dacos(-1.)
double precision, parameter :: eps = 1.e-9
double precision, parameter :: a=2., b=2.

inside=0.
phi=huge(1.)

dx = (xmax - xmin)/(Nx-1)
dy = (ymax - ymin)/(Ny-1)
x(1) = xmin
y(1) = ymin

do i=1,Nx-1
  x(i+1) = x(i) + dx
enddo
do j=1,Ny-1
  y(j+1) = y(j) + dy
enddo

do j=1,Ny
  do i=1,Nx
    eck = (x(i)/a)**2 + (y(j)/b)**2
	if(eck <= 1.) inside(i,j) = 1.
    call min_dist_to_ellipse(a,b,x(i),y(j), eps,phi(i,j))
	!phi(i,j) = dabs(phi(i,j))
	!if(inside(i,j) > 0.) phi(i,j) = 0.
  enddo
enddo
!  Create tecplot formatted velocity field file  
open (unit = 2,file = 'test.dat', status='unknown',form='formatted', &
  action='write',position='rewind')

write(2,*) 'variables = "x", "y", "phi", "inside"'; 

write(2,"(1a,i3,1a,i3)") 'ZONE T="1", DATAPACKING=POINT, i=', Nx,', j=',Ny

write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE)')//''

  do j=1,ny
    do i=1,nx
      write(2,*) x(i), y(j), phi(i,j), inside(i,j)
    enddo
  enddo

close(2)

stop
end program test
