!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by: 
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   The Actuator Turbine Model is free software: you can redistribute it 
!!   and/or modify it under the terms of the GNU General Public License as 
!!   published by the Free Software Foundation, either version 3 of the 
!!   License, or (at your option) any later version.
!!
!!   The Actuator Turbine Model is distributed in the hope that it will be 
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module atm_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module provides basic functionalities to the actuator turbine model
! Real precision variable and dynamic allocation types are stored here

implicit none

! Precision of real numbers
integer, parameter :: rprec = kind (1.d0)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine error (msg)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
character (*), intent (in) :: msg

write (*, '(1x,a)') '  '
write (*, '(1x,a)') '*****ERROR*****'
write (*, '(1x,a)') 'In atm_input_util:'
write (*, '(1x,a)') trim (msg)
write (*, '(1x,a)') '***************'
write (*, '(1x,a)') 'Program aborted'
write (*, '(1x,a)') '  '

stop
end subroutine error

!-------------------------------------------------------------------------------
function interpolate(xp,x,y)
! This function interpolates xp from x and y 
!-------------------------------------------------------------------------------
real(rprec), dimension(:), intent(in) :: x,y
real(rprec), intent(in) ::  xp
real(rprec) :: xa, xb, ya, yb
integer :: i,p
real(rprec) :: interpolate
p=size(x)

if (xp <= x(1)) then 
    interpolate=y(1)
else if (xp>=x(p)) then
    interpolate=y(p)
else
    do i=2,p
        if ( ( xp .ge. x(i-1) ) .and. ( xp .le. x(i) ) ) then
            xa=x(i-1)
            xb=x(i)
            ya=y(i-1)
            yb=y(i)
            interpolate = ya + (yb-ya) * (xp-xa) / (xb-xa) 
        endif
    enddo
endif
return
end function interpolate

!-------------------------------------------------------------------------------
integer function interpolate_i(xp,x,y)
! This function interpolates xp from x and y 
!-------------------------------------------------------------------------------
real(rprec), dimension(:), intent(in) :: x
integer, dimension(:), intent(in) :: y
real(rprec), intent(in) ::  xp
real(rprec) :: xa,xb,ya,yb
integer :: i,p
p=size(x)

if (xp .lt. x(1)) then
    interpolate_i=y(1) 
    else if (xp .gt. x(p)) then
        interpolate_i=y(p)
    else
        do i=2,p
            if ( ( xp .ge. x(i-1) ) .and. ( xp .le. x(i) ) ) then
                xa=x(i-1)
                xb=x(i)
                ya=real(y(i-1),rprec)
                yb=real(y(i),rprec)
                interpolate_i=nint( ya + (yb-ya) * (xp-xa) / (xb-xa) )
!write(*,*) 'Interpolation = ', nint( ya + (yb-ya) * (xp-xa) / (yb-ya) )
                endif
        enddo
endif
!write(*,*) 'Value of y = ', ya, yb
!write(*,*) 'Value of x = ', xa, xb
!write(*,*) 'Value of xp = ', xp
!write(*,*) 'Value of p = ', p
!write(*,*) 'Interpolated Value = ', interpolate_i
return
end function interpolate_i

!-------------------------------------------------------------------------------
function vector_add(a,b)
! This function adds 2 vectors (arrays real(rprec), dimension(3))
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a,b
real(rprec), dimension(3) :: vector_add
vector_add(1)=a(1)+b(1)
vector_add(2)=a(2)+b(2)
vector_add(3)=a(3)+b(3)
return
end function vector_add

!-------------------------------------------------------------------------------
function vector_divide(a,b)
! This function divides one vector (array real(rprec), dimension(3) by a number)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec), intent(in) :: b
real(rprec), dimension(3) :: vector_divide
vector_divide(1)=a(1)/b
vector_divide(2)=a(2)/b
vector_divide(3)=a(3)/b
return
end function vector_divide

!-------------------------------------------------------------------------------
function vector_multiply(a,b)
! This function multiplies one vector (array real(rprec), dimension(3) by
! a real(rprec) number)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec), intent(in) :: b
real(rprec), dimension(3) :: vector_multiply
vector_multiply(1)=a(1)*b
vector_multiply(2)=a(2)*b
vector_multiply(3)=a(3)*b
return
end function vector_multiply

!-------------------------------------------------------------------------------
function vector_mag(a)
! This function calculates the magnitude of a vector
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec) :: vector_mag
vector_mag=abs(sqrt(a(1)**2+a(2)**2+a(3)**2))
return
end function vector_mag

!-------------------------------------------------------------------------------
function rotatePoint(point_in, rotationPoint, axis, angle)
! This function performs rotation of a point with respect to an axis or rotation
! and a certain angle
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: point_in
real(rprec), dimension(3), intent(in) :: rotationPoint
real(rprec), dimension(3), intent(in) :: axis
real(rprec), intent(in) :: angle
real(rprec), dimension(3,3) :: RM ! Rotation Matrix tensor
real(rprec), dimension(3) :: rotatePoint, point

point=point_in

RM(1,1) = axis(1)**2 + (1.0 - axis(1)**2) * cos(angle)
RM(1,2) = axis(1) * axis(2) * (1.0 - cos(angle)) - axis(3) * sin(angle)
RM(1,3) = axis(1) * axis(3) * (1.0 - cos(angle)) + axis(2) * sin(angle)
RM(2,1) = axis(1) * axis(2) * (1.0 - cos(angle)) + axis(3) * sin(angle)
RM(2,2) = axis(2)**2 + (1.0 - axis(2)**2) * cos(angle)
RM(2,3) = axis(2) * axis(3) * (1.0 - cos(angle)) - axis(1) * sin(angle)
RM(3,1) = axis(1) * axis(3) * (1.0 - cos(angle)) - axis(2) * sin(angle)
RM(3,2) = axis(2) * axis(3) * (1.0 - cos(angle)) + axis(1) * sin(angle)
RM(3,3) = axis(3)**2 + (1.0 - axis(3)**2) * cos(angle)

! Rotation matrices make a rotation about the origin, so need to subtract 
! rotation point off the point to be rotated
point=vector_add(point,-rotationPoint)

! Perform rotation (multiplication matrix and vector)
point=matrix_vector(RM,point)

! Return the rotated point to its new location relative to the rotation point
rotatePoint = point + rotationPoint

return 
end function rotatePoint

!-------------------------------------------------------------------------------
function matrix_vector(RM,point)
! This function multiplies a matrix and a vector
!-------------------------------------------------------------------------------
real(rprec), dimension(3,3), intent(in) :: RM ! Matrix
real(rprec), dimension(3), intent(in) :: point ! vector point
real(rprec), dimension(3) :: matrix_vector
! Perform rotation
matrix_vector(1)=RM(1,1)*point(1)+RM(1,2)*point(2)+RM(1,3)*point(3)
matrix_vector(2)=RM(2,1)*point(1)+RM(2,2)*point(2)+RM(2,3)*point(3)
matrix_vector(3)=RM(3,1)*point(1)+RM(3,2)*point(2)+RM(3,3)*point(3)
return
end function matrix_vector

!-------------------------------------------------------------------------------
function cross_product(u,v)
! This function calculates the cross product of 2 vectors
!-------------------------------------------------------------------------------
real(rprec), intent(in) :: u(3),v(3)
real(rprec) :: cross_product(3)
cross_product(1) = u(2)*v(3)-u(3)*v(2)
cross_product(2) = u(3)*v(1)-u(1)*v(3)
cross_product(3) = u(1)*v(2)-u(2)*v(1)
return
end function cross_product

!-------------------------------------------------------------------------------
function distance(a,b)
! This function calculates the distance between a(1,2,3) and b(1,2,3)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a,b
real(rprec) :: distance
distance=abs(sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2))
return
end function distance

!-------------------------------------------------------------------------------
character(len=20) function int2str(k)
!-------------------------------------------------------------------------------
! This function converts an integer to string
! http://stackoverflow.com/questions/1262695/
!converting-integers-to-strings-in-fortran
integer, intent(in) :: k
write (int2str, *) k
int2str = adjustl(int2str)
end function int2str





end module atm_base
