!!
!!  Copyright (C) 2011-2013  Johns Hopkins University
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

!*********************************************************************
module emul_complex
!*********************************************************************
! 
! The purpose of this module is to provide methods for performing
! complex operations against real arrays emulating complex arrays.  The
! real array is to contain interleaved complex information.
!  
! Written by : 
!     Jason Graham <jgraha8@gmail.com>
!

use types, only : rprec
use messages, only : error
implicit none

save
private 

! public :: mult_real_complex, &
!      mult_real_complex_imag, &
!      mult_real_complex_real 

public :: operator( .MUL. ), &
     operator( .MULI. ), &
     operator( .MULR. )

character (*), parameter :: mod_name = 'emul_complex'

!///////////////////////////////////////
!/// OPERATORS                       ///
!///////////////////////////////////////

! REAL X COMPLEX
interface operator (.MUL.) 
   module procedure &
        mul_real_complex_2D
end interface

! REAL X IMAG(COMPLEX)
interface operator (.MULI.) 
   module procedure &
        mul_real_complex_imag_scalar, &
        mul_real_complex_imag_2D 
end interface

! REAL X REAL(COMPLEX)
interface operator (.MULR.)
   module procedure &
        mul_real_complex_real_2D
end interface

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function mul_real_complex_imag_scalar( a, a_c ) result(b)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This function emulates the multiplication of two complex scalars
!  by emulating the input real vector (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0).
!
!  Input:
!  
!    a (real,size(2,1))  - input real vector
!    a_c (real)          - input imaginary part of complex scalar
!
!  Output:
!
!    b (real, size(2,1)) - output real vector
!
implicit none

character (*), parameter :: sub_name = mod_name // '.mul_real_complex_imag'

real(rprec), dimension(2), intent(in) :: a
real(rprec), intent(in) :: a_c

real(rprec), dimension(2) :: b

!  Cached variables
real(rprec) :: a_c_i
  
!  Cache multi-usage variables
a_c_i = a_c

!  Perform multiplication
b(1) = - a(2) * a_c_i
b(2) =  a(1) * a_c_i

return

end function mul_real_complex_imag_scalar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function mul_real_complex_2D( a, a_c ) result(b)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This function emulates the multiplication of two complex 2D array by
!  emulating the input real array (a) as a complex type.
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (complex,size(nx_c,ny))- input complex array 
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!
!  Note: nx_c must be nx_r/2
!  
implicit none

character (*), parameter :: sub_name = mod_name // '.mul_real_complex_2D'

real(rprec), dimension( :, :), intent(in) :: a
complex(rprec), dimension( :, : ), intent(in) :: a_c

real(rprec), allocatable, dimension(:, :) :: b

!  Cached variables
real(rprec) ::  a_r, a_i, a_c_r, a_c_i

integer :: i,j,ir,ii
integer :: nx_r, nx_c, ny

! Get the size of the incoming arrays
nx_r = size(a,1)
ny   = size(a,2)

nx_c = size(a_c,1)

$if ($SAFETYMODE)
if ( nx_r .NE. 2*nx_c .OR. &
     ny .NE. size(a_c,2) ) call error( sub_name, 'Mismatch in input array sizes')
$endif

! Allocate returned array
allocate( b(nx_r, ny) )

!  Emulate complex multiplication
!  Using outer loop to get contiguous memory access
do j=1, ny
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)
    a_c_r = real(a_c(i,j),kind=rprec)
    $if($DBLPREC)
    a_c_i = dimag(a_c(i,j))
    $else
    a_c_i = aimag(a_c(i,j))
    $endif
    
    !  Perform multiplication
    b(ir,j) = a_r * a_c_r - a_i * a_c_i
    b(ii,j) = a_r * a_c_i + a_i * a_c_r

  enddo
enddo

return

end function mul_real_complex_2D


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function mul_real_complex_imag_2D( a, a_c ) result(b)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This function emulates the multiplication of two complex 2D array by
!  emulating the input real array (a) as a complex type. This subroutine
!  ignores the real part of a_c (e.g. would use this when real(a_c) = 0)
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (real,size(nx_c,ny))   - input imaginary part of complex array 
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!
!  Note: nx_c must be nx_r/2
!
implicit none

character (*), parameter :: sub_name = mod_name // '.mul_real_complex_imag_2D'

real(rprec), dimension( :, : ), intent(in) :: a
real(rprec), dimension( :, : ), intent(in) :: a_c

real(rprec), allocatable, dimension(:, :) :: b

!  Cached variables
real(rprec) ::  a_c_i, cache

integer :: i,j,ii,ir
integer :: nx_r, nx_c, ny

! Get the size of the incoming arrays
nx_r = size(a,1)
ny   = size(a,2)

nx_c = size(a_c,1)

$if ($SAFETYMODE)
if ( nx_r .NE. 2*nx_c .OR. &
     ny .NE. size(a_c,2) ) call error( sub_name, 'Mismatch in array sizes')
$endif

! Allocate the returned array
allocate( b(nx_r, ny ) )

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contiguous memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_c_i = a_c(i,j)

    !  Perform multiplication (cache data to ensure sequential access)
    cache = a(ir,j) * a_c_i
    b(ir,j) = - a(ii,j) * a_c_i
    b(ii,j) =  cache

  enddo
enddo

return

end function mul_real_complex_imag_2D


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function mul_real_complex_real_2D( a, a_c ) result(b)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This function emulates the multiplication of two complex 2D array by
!  emulating the input real array (a) as a complex type. This subroutine
!  ignores the imaginary part of a_c (e.g. would use this when imag(a_c)
!  = 0).
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input/output real array
!    a_c (real,size(nx_c,ny))   - input real part of complex array 
!
!  Output:
!
!    b (real, size(nx_r,ny))    - output real array
!
!  Note: nx_c must be nx_r/2
!
use types, only : rprec
implicit none

character (*), parameter :: sub_name = mod_name // '.mul_real_complex_real_2D'

real(rprec), dimension( :, : ), intent(in) :: a
real(rprec), dimension( :, : ), intent(in) :: a_c

real(rprec), allocatable, dimension(:, :) :: b

integer :: i,j,ii,ir
integer :: nx_r, nx_c, ny

! Get the size of the incoming arrays
nx_r = size(a,1)
ny   = size(a,2)

nx_c = size(a_c,1)

$if ($SAFETYMODE)
if ( nx_r .NE. 2*nx_c .OR. &
     ny .NE. size(a_c,2) ) call error( sub_name, 'Mismatch in array sizes')
! Allocate the returned array
$endif

allocate(b(nx_r,ny))

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contiguous memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1

    !  Perform multiplication
    b(ir:ii,j) = a_c(i,j)*a(ir:ii,j)

  enddo
enddo

return

end function mul_real_complex_real_2D

end module emul_complex
