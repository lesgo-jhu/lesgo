!!
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

module mpi_transpose_mod
use types, only : rp => rprec
use param, only : np => nproc, comm_cart => comm, coord, ierr, MPI_CPREC
use mpi
implicit none

save
private
public :: mpi_transpose

interface mpi_transpose
  module procedure mpi_transpose_c
end interface
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--the sizes of a, b are assumed shape here to allow us to use arrays
!  dimensioned for the Nyquist frequency in x-direction
!--only a(1:mx, 1:my, 1:mz) & b(1:mz*np, 1:my, 1:mx/np) are used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mpi_transpose_c (mx, my, mz, a, b)
implicit none

integer, intent (in) :: mx, my, mz  !--declare as arguments so untransposing
                                    !  does not require additional code
                                    !--may need these to be parameters, so
                                    !  further changes may be required
complex (rp), intent (in) :: a(:, :, :)
complex (rp), intent (out) :: b(:, :, :)
!complex (rp), intent (in) :: a(mx, my, mz)
!complex (rp), intent (out) :: b(mz*np, my, mx/np)

integer :: bs
integer :: ip
integer :: up, down
integer :: status(MPI_STATUS_SIZE)

integer :: i, k, jx, jz

logical, save :: init = .false.

complex (rp) :: tmpout(mx/np, my, mz), tmpin(mx/np, my, mz)

!integer, save :: src(np-1), dest(np-1)
integer, save, allocatable, dimension(:) :: src, dest

logical, save :: arrays_allocated = .false.

!---------------------------------------------------------------------

if( .not. arrays_allocated ) then
   allocate(src(np-1))
   allocate(dest(np-1))
   
   arrays_allocated = .true.

endif

if (.not. init) then
  do ip = 1, np-1
    !--this is a bit awkward: really want periodic topology so can
    !  use cart_shift here but this is not good for the finite
    !  differences, so perhaps create a new cartisian topology from 
    !  existing one, but we no not allow reordering and we do allow 
    !  periodicity
    up = modulo (coord + ip, np)  !--corresponds to dest(ip)
    down = modulo (coord - ip, np)  !--corresponds to src(ip)

    call MPI_cart_rank (comm_cart, (/ up /), dest(ip), ierr)
    call MPI_cart_rank (comm_cart, (/ down /), src(ip), ierr)
  end do

  init = .true.
end if

!--block size
bs = mx*my*mz/np

do ip = 1, np-1

  up = modulo (coord + ip, np)  !--corresponds to dest(ip)
  down = modulo (coord - ip, np)  !--corresponds to src(ip)

  !--copy chunk "up" into buffer (no local transpose)
  do jz = 1, mz
    do jx = 1, mx/np      
      tmpout(jx, :, jz) = a(up*mx/np+jx, :, jz)
    end do
  end do

  !write (*, *) $str($context_doc), coord, ' before sendrecv'
  call MPI_sendrecv (tmpout(1,1,1), bs, MPI_CPREC, dest(ip), ip,  &
                     tmpin(1,1,1), bs, MPI_CPREC, src(ip), ip,    &
                     comm_cart, status, ierr)
  !write (*, *) $str($context_doc), coord, 'after sendrecv'

  !--copy chunk "down" from tmpin to b, in transposed order
  !--may want to keep in natural order and then tranpose whole array 
  !  after all mpi communication is complete
  do i = 1, mx/np
    do k = 1, mz
      jz = down*mz + k
      b(jz, :, i) = tmpin(i, :, k)
    end do
  end do
  
end do

!--local transpose on non-transferred data
!--chunk 'rank' should not have have been sent/received
do i = 1, mx/np
  do k = 1, mz
    jx = coord*mx/np + i
    jz = coord*mz + k
    b(jz, :, i) = a(jx, :, k)
  end do
end do

end subroutine mpi_transpose_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module mpi_transpose_mod
