!!
!!  Copyright (C) 2016-2016  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

module util
    
use types, only : rprec
use param, only : pi, CHAR_BUFF_LENGTH
use messages
use string_util
implicit none

public

interface interpolate
    module procedure :: interpolate_scalar
    module procedure :: interpolate_array
end interface interpolate

contains

! linear interpolation between 1D-arrays
subroutine interpolate_scalar(x, y, xi, yi)
    implicit none
    real(rprec), dimension(:), intent(in)  :: x, y
    real(rprec), intent(in)                :: xi
    real(rprec), intent(out)               :: yi
    real(rprec), dimension(:), allocatable :: xi_array, yi_array
    
    allocate( xi_array(1) )
    allocate( yi_array(1) )
    xi_array(1) = xi
    call interpolate(x, y, xi_array, yi_array)
    yi = yi_array(1)
    
end subroutine interpolate_scalar

! linear interpolation between 1D-arrays
subroutine interpolate_array(x, y, xi, yi)
    implicit none
    real(rprec), dimension(:), intent(in)  :: x, y, xi
    real(rprec), dimension(:), intent(out) :: yi
    integer     :: i, j   
    real(rprec) :: dx, t 
    
    if ( size(x) /= size(y) .or. size(xi) /= size(yi)) then
        call error('util.inteprolate','Interpolation pairs must be of equal size.')
    end if
    if ( size(x) /= size(y) ) then
        call error('util.inteprolate','Interpolation pairs must be of equal size.')
    end if
    do i = 2, size(x)
        if ( x(i-1) >= x(i) ) then
            call error('util.interpolate','array x must be monotonically increasing')
        end if
    end do
    do i = 2, size(xi)
        if ( xi(i-1) >= xi(i) ) then
            call error('util.interpolate','array xi must be monotonically increasing')
        end if
    end do

    j = 1
    do i = 1, size(xi)
        if ( xi(i) <= x(1) ) then
            yi(i) = y(1)
        else if ( xi(i) >= x( size(x) ) ) then
            yi(i) = y( size(x) )
        else
            do while ( xi(i) > x(j+1) )
                j = j + 1
            end do
            dx = x(j+1) - x(j)
            t = ( xi(i) - x(j) ) / dx
            yi(i) = (1-t) * y(j) + t * y(j+1)
        end if
    end do
    
end subroutine interpolate_array

end module util