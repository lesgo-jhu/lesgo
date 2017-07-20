!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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

!*******************************************************************************
module mpi_transpose_mod
!*******************************************************************************
use types, only : rprec
use param, only : np => nproc, comm_cart => comm, coord, ierr, MPI_CPREC
use mpi
implicit none

save
private
public :: mpi_transpose

contains

!*******************************************************************************
subroutine mpi_transpose (mx, my, mz, a, b)
!*******************************************************************************
!
!--the sizes of a, b are assumed shape here to allow us to use arrays
!  dimensioned for the Nyquist frequency in x-direction
!--only a(1:mx, 1:my, 1:mz) & b(1:mz*np, 1:my, 1:mx/np) are used
!
implicit none

! declare as arguments so untransposing does not require additional code
integer, intent (in) :: mx, my, mz
complex(rprec), intent (in) :: a(:, :, :)
complex(rprec), intent (out) :: b(:, :, :)
integer :: bs
integer :: ip
integer :: up, down
integer :: status(MPI_STATUS_SIZE)
integer :: i, k, jx, jz
logical, save :: init = .false.
complex(rprec) :: tmpout(mx/np, my, mz), tmpin(mx/np, my, mz)
integer, save, allocatable, dimension(:) :: src, dest
logical, save :: arrays_allocated = .false.

if( .not. arrays_allocated ) then
    allocate(src(np-1))
    allocate(dest(np-1))
    arrays_allocated = .true.
endif

if (.not. init) then
    do ip = 1, np-1
        ! this is a bit awkward: really want periodic topology so can
        ! use cart_shift here but this is not good for the finite
        ! differences, so perhaps create a new cartisian topology from
        ! existing one, but we no not allow reordering and we do allow
        ! periodicity
        up = modulo (coord + ip, np)    ! corresponds to dest(ip)
        down = modulo (coord - ip, np)  ! corresponds to src(ip)

        call MPI_cart_rank (comm_cart, (/ up /), dest(ip), ierr)
        call MPI_cart_rank (comm_cart, (/ down /), src(ip), ierr)
    end do

    init = .true.
end if

!--block size
bs = mx*my*mz/np

do ip = 1, np-1
    up = modulo (coord + ip, np)    ! corresponds to dest(ip)
    down = modulo (coord - ip, np)  ! corresponds to src(ip)

    ! copy chunk "up" into buffer (no local transpose)
    do jz = 1, mz
    do jx = 1, mx/np
        tmpout(jx, :, jz) = a(up*mx/np+jx, :, jz)
    end do
    end do

    call MPI_sendrecv (tmpout(1,1,1), bs, MPI_CPREC, dest(ip), ip,             &
        tmpin(1,1,1), bs, MPI_CPREC, src(ip), ip, comm_cart, status, ierr)

    ! copy chunk "down" from tmpin to b, in transposed order
    do i = 1, mx/np
    do k = 1, mz
        jz = down*mz + k
        b(jz, :, i) = tmpin(i, :, k)
    end do
    end do
end do

! local transpose on non-transferred data
! chunk 'rank' should not have have been sent/received
do i = 1, mx/np
do k = 1, mz
    jx = coord*mx/np + i
    jz = coord*mz + k
    b(jz, :, i) = a(jx, :, k)
end do
end do

end subroutine mpi_transpose

end module mpi_transpose_mod
