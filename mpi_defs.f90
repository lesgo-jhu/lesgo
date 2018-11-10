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
module mpi_defs
!*******************************************************************************
use mpi
implicit none

save
private

public :: initialize_mpi, mpi_sync_real_array
public :: MPI_SYNC_DOWN, MPI_SYNC_UP, MPI_SYNC_DOWNUP

character (*), parameter :: mod_name = 'mpi_defs'

integer, parameter :: MPI_SYNC_DOWN=1
integer, parameter :: MPI_SYNC_UP=2
integer, parameter :: MPI_SYNC_DOWNUP=3

#ifdef PPCGNS
integer, public :: cgnsParallelComm
#endif

contains

!*******************************************************************************
subroutine initialize_mpi()
!*******************************************************************************
use types, only : rprec
use param
#ifdef PPCPS
use concurrent_precursor
#endif
#ifdef PPCGNS
use cgns
#endif
implicit none

! Set the local communicator
#ifdef PPCPS
    ! Create the local communicator (split from MPI_COMM_WORLD)
    ! This also sets the globally defined intercommunicator (bridge)
    call create_mpi_comms_cps( globalComm )
#else
    globalComm = MPI_COMM_WORLD
#endif

! call mpi_comm_size (globalComm, nproc, ierr)
call mpi_comm_rank (globalComm, global_rank, ierr)

! set the MPI_RPREC variable
if (rprec == kind (1.e0)) then
    MPI_RPREC = MPI_REAL
    MPI_CPREC = MPI_COMPLEX
else if (rprec == kind (1.d0)) then
    MPI_RPREC = MPI_DOUBLE_PRECISION
    MPI_CPREC = MPI_DOUBLE_COMPLEX
else
    write (*, *) 'error defining MPI_RPREC/MPI_CPREC'
    stop
end if

end subroutine initialize_mpi

!*******************************************************************************
subroutine mpi_sync_real_array( var, lbz, isync )
!*******************************************************************************
!
! This subroutine provides a generic method for syncing arrays in
! lesgo. This method applies to arrays indexed in the direction starting
! from both 0 and 1. For arrays starting from index of 1, only the
! SYNC_DOWN procedure may be performed. No assumption is made about the
! dimensions of the other directions (x and y) and the bounds of these
! indices are obtained directly.
!
! The synchronization is provided according to the following rules:
!
! SYNC_DOWN : Send data from k = 1 at coord+1 to k=nz at coord
! SYNC_UP   : Send data from k = nz-1 at coord to k=0 at coord+1
!
! Arguments:
!
! var   : three dimensional array to be sync'd accross processors
! lbz   : the lower bound of var for the z index; its specification resolves
!         descrepencies between arrays indexed starting at 0 from those at 1
! isync : flag for determining the type of synchronization and takes on values,
!         MPI_SYNC_DOWN, MPI_SYNC_UP, or MPI_SYNC_DOWNUP from the MPI_DEFS
!         module.
!
use types, only : rprec
use mpi
use param, only : MPI_RPREC, down, up, comm, status, ierr, nz
use messages

implicit none

character (*), parameter :: sub_name = mod_name // '.mpi_sync_real_array'

real(rprec), dimension(:,:,lbz:), intent(INOUT) :: var
integer, intent(in) :: lbz
integer, intent(in) :: isync

integer :: sx, sy
integer :: ubz
integer :: mpi_datasize

! Get size of var array in x and y directions
sx = size(var,1)
sy = size(var,2)
! Get upper bound of z index; the lower bound is specified
ubz = ubound(var,3)

! We are assuming that the array var has nz as the upper bound - checking this
if( ubz .ne. nz ) call error( sub_name, 'Input array must lbz:nz z dimensions.')

!  Set mpi data size
mpi_datasize = sx*sy

if (isync == MPI_SYNC_DOWN) then
    call sync_down()
else if( isync == MPI_SYNC_UP) then
    if( lbz /= 0 ) call error( sub_name,                                       &
        'Cannot SYNC_UP with variable with non-zero starting index')
    call sync_up()
else if( isync == MPI_SYNC_DOWNUP) then
    if( lbz /= 0 ) call error( sub_name,                                       &
        'Cannot SYNC_DOWNUP with variable with non-zero starting index')
    call sync_down()
    call sync_up()
else
   call error( sub_name, 'isync not specified properly')
end if

if(ierr .ne. 0) call error( sub_name,                                          &
    'Error occured during mpi sync; recieved mpi error code :', ierr)

! Enforce globally synchronous MPI behavior. Most likely safe to comment
! out, but can be enabled to ensure absolute safety.
!call mpi_barrier( comm, ierr )

contains

!*******************************************************************************
subroutine sync_down()
!*******************************************************************************
implicit none

call mpi_sendrecv (var(:,:,1), mpi_datasize, MPI_RPREC, down, 1,               &
    var(:,:,ubz), mpi_datasize, MPI_RPREC, up, 1, comm, status, ierr)

end subroutine sync_down

!*******************************************************************************
subroutine sync_up()
!*******************************************************************************
implicit none

call mpi_sendrecv (var(:,:,ubz-1), mpi_datasize, MPI_RPREC, up, 2,             &
    var(:,:,0), mpi_datasize, MPI_RPREC, down, 2, comm, status, ierr)

end subroutine sync_up

end subroutine mpi_sync_real_array

end module mpi_defs
