!!
!!  Copyright 2009,2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

module integrate
use types, only : rp => rprec
implicit none

public

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--a, b do NOT need to line up with grid points, but this means
!  additional accuracy is required to do the cut "cells"
!--x is assumed to equi-spaced
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function integrate_1d (periodic, ain, bin, dx, x, f)
implicit none

real (rp) :: integrate_1d

logical, intent (in) :: periodic  !--whether or not periodic

real (rp), intent (in) :: ain, bin, dx
real (rp), intent (in) :: x(:), f(:)

real (rp), parameter :: thresh = 1.e-2_rp

integer :: sgn
integer :: n
integer :: imin, imin_m1, imin_p1, imin_0
integer :: imax, imax_m1, imax_p1, imax_0
integer :: i, ii

real (rp) :: a, b
real (rp) :: c1, c2, c3
real (rp) :: left_part, right_part, middle
real (rp) :: lfr, rfr
real (rp) :: xmin, xmax

!---------------------------------------------------------------------

n = size (f)

if (n /= size (x)) then
  write (*, *) 'integrate_1d: sizes of x, f are not the same'
  stop
end if

!--make it so we can assume a <= b
if (ain <= bin) then
  a = ain
  b = bin
  sgn = 1
else
  a = bin
  b = ain
  sgn = -1
end if

if (.not. periodic) then
  if ((a < x(1)) .or. (x(n) < b)) then
    write (*, *) 'a = ', a
    write (*, *) 'b = ', b
    write (*, *) 'x(1) = ', x(1)
    write (*, *) 'x(n) = ', x(n)
    write (*, *) 'integrate_1d: a or b is out of range of x'
    stop
  end if
end if

!--find outermost grid pts inside interval (a, b)
imin_0 = ceiling ((a - x(1)) / dx + 1)
imax_0 = floor ((b - x(1)) / dx + 1)

!--we will need to apply bdry treatment at some point
if (periodic) then

  imin = modulo (imin_0 - 1, n) + 1
  imax = modulo (imax_0 - 1, n) + 1
  
else

  if ((imin_0 < 1) .or. (imax_0 > n)) then
    write (*, *) 'integrate_1d: imin or imax out of range'
    stop
  else

    imin = imin_0
    imax = imax_0

  end if

end if

xmin = x(imin)
xmax = x(imax)

!--find fractional part of cut cells
lfr = (xmin - a) / dx
rfr = (b - xmax) / dx

if (lfr > thresh) then  !--perform integration over left cut cell
  if (periodic) then

    imin_p1 = modulo (imin_0, n) + 1
    imin_m1 = modulo (imin_0 - 2, n) + 1

  else

    !--we have a problem when near the boundary!
    imin_p1 = imin_0 + 1
    imin_m1 = imin_0 - 1

    if ((imin_p1 > n) .or. (imin_m1 < 1)) then
      write (*, *) 'integrate_1d: trying to integrate outside of domain'
      write (*, *) 'n = ', n
      write (*, *) 'imin_p1 = ', imin_p1
      write (*, *) 'imin_m1 = ', imin_m1
      write (*, *) 'imin_0 = ', imin_0
      stop
    end if
  
  end if

  !--coefficients in a + b (x - xmin) + c (x - xmin)^2
  c1 = f(imin)
  c2 = (f(imin_p1) - f(imin_m1)) / (2._rp * dx)
  c3 = (f(imin_p1) - 2._rp * f(imin) + f(imin_m1)) / (2._rp * dx**2)

  !--integral from xmin - lfr * dx to xmin
  left_part = c1 * lfr * dx - (c2 / 2._rp ) * (lfr * dx)**2 +  &
              (c3 / 3._rp) * (lfr * dx)**3

else

  left_part = 0._rp

end if

if (rfr > thresh) then  !--perform integration over right cut cell

  if (periodic) then

    imax_p1 = modulo (imax_0, n) + 1
    imax_m1 = modulo (imax_0 - 2, n) + 1

  else

    !--we have a problem when near the boundary!
    imax_p1 = imax_0 + 1
    imax_m1 = imax_0 - 1

    if ((imax_p1 > n) .or. (imax_m1 < 1)) then
      write (*, *) 'integrate_1d: trying to integrate outside of domain'
      write (*, *) 'n = ', n
      write (*, *) 'imax_p1 = ', imax_p1
      write (*, *) 'imax_m1 = ', imax_m1
      write (*, *) 'imax_0 = ', imax_0
      stop
    end if
  
  end if

  c1 = f(imax)
  c2 = (f(imax_p1) - f(imax_m1)) / (2._rp * dx)
  c3 = (f(imax_p1) - 2._rp * f(imax) + f(imax_m1)) / (2._rp * dx **2)

  !--integral from xmax to xmax + rfr * dx
  right_part = c1 * rfr * dx + (c2 / 2._rp) * (rfr * dx)**2 +  &
               (c3 / 3._rp) * (rfr * dx)**3
    
else

  right_part = 0._rp

end if

!--perform integration over middle cells (Trapezoidal rule for now)
middle = 0.5_rp * (f(imin) + f(imax))

do i = imin_0 + 1, imax_0 - 1

  if (periodic) then
    ii = modulo (i - 1, n) + 1
  else
    ii = i
  end if
  
  middle = middle + f(ii)
  
end do

middle = middle * dx

!--add up pieces to get integral
integrate_1d = sgn * (left_part + right_part + middle)

end function integrate_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!program test
!use types, only : rp => rprec
!use integrate
!implicit none
!
!integer, parameter :: n = 8
!
!integer :: j
!
!real (rp) :: integral
!real (rp) :: x(n), f(n)
!real (rp) :: a, b
!real (rp) :: lmin, lmax
!real (rp) :: dx
!
!!---------------------------------------------------------------------
!
!lmin = 0._rp
!lmax = 1._rp
!dx = (lmax - lmin) / (n - 1)
!
!a = 0._rp
!b = 1._rp
!
!x = (/ ( lmin + (j - 1) * dx, j=1, n ) /)
!
!f = x**2
!
!integral = integrate_1d (.false., a, b, dx, x, f)
!
!write (*, *) 'integral = ', integral
!
!end program test
