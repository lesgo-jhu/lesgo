!!  Copyright (C) 2009-2013  Johns Hopkins University
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

!! This subroutine performs a linear regression analysis for a set of data 
!! given as (x,y) pairs. The output from the subroutine is the slope and 
!! y-intercept of the least squares best fit straight line through the 
!! data points. Adapted from: www.pgccphy.net/Linreg/linreg_f90.txt

subroutine linreg(x,y,b,m)

use types, only: rprec
implicit none

real(rprec), dimension(:), intent(in) :: x, y
real(rprec) :: b, m



end subroutine linreg
