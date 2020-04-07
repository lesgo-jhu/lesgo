!!
!!  Copyright (C) 2009-2020  Johns Hopkins University
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
module fringe
!*******************************************************************************
use types, only : rprec
use mpi_defs
use messages

implicit none

private
public fringe_t

type fringe_t
    ! number of x grid points in fringe region
    integer :: nx
    ! Start of fringe region
    integer :: istart
    ! Start of plateau of fringe region
    integer :: iplateau
    ! End of fringe region
    integer :: iend
    ! End of fringe region as a fraction of L_x
    real(rprec) :: l_end
    ! Length of fringe region as a fraction of L_x
    real(rprec) :: l_len
    ! Wrapped locations
    integer, allocatable, dimension(:) :: iwrap
    ! Weighting functions
    real(rprec), allocatable, dimension(:) :: alpha, beta
contains

end type fringe_t

interface fringe_t
    module procedure :: constructor
end interface fringe_t

contains


!*******************************************************************************
function constructor(l_end, l_len) result(this)
!*******************************************************************************
use param, only : nx, pi
type(fringe_t) :: this
real(rprec), intent(in) :: l_end, l_len
integer :: i

! Copy from input arguments
this%l_end = l_end
this%l_len = l_len

! Location of start, end, and plateau locations, and size of fringe
this%iend = floor(this%l_end * nx + 1.0_rprec)
this%iplateau = floor(( this%l_end - this%l_len / 4 ) * nx + 1.0_rprec)
this%istart = floor((this%l_end - this%l_len) * nx + 1.0_rprec)
this%nx = this%iend - this%istart

! Weighting functions
allocate(this%alpha(this%nx))
allocate(this%beta(this%nx))
do i = 1, this%nx
    if (i+this%istart > this%iplateau) then
        this%beta(i) = 1.0_rprec
    else
        this%beta(i) = 0.5_rprec * (1.0_rprec-cos(pi*real(i,rprec)             &
            / (this%iplateau-this%istart)))
    endif
    this%alpha(i) = 1._rprec - this%beta(i)
end do

! Allocate and assign wrapped index and fringe weights
allocate(this%iwrap(this%nx))
do i = 1, this%nx
    this%iwrap(i) = modulo( this%istart+i-1, nx ) + 1
enddo

end function constructor

end module fringe
