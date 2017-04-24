!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
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
module concurrent_precursor
!*******************************************************************************
use types, only : rprec
use fringe_region
use mpi
implicit none

save
private

public :: interComm, color, RED, BLUE
public :: create_mpi_comms, &
          cps_init, &
          synchronize_cps

character (*), parameter :: mod_name = 'concurrent_precursor'

integer, parameter :: RED=0  ! Upstream domain   (producer)
integer, parameter :: BLUE=1 ! Downstream domain (consumer)

integer :: interComm, color

contains

!*******************************************************************************
subroutine create_mpi_comms(localComm)
!*******************************************************************************
!
! This subroutine does two things. It first splits the MPI_COMM_WORLD
! communicator into two communicators (localComm). The two new
! communicators are then bridged to create an intercommunicator
! (interComm).
!
use param, only : ierr, use_cps, inflow_cond
use messages
implicit none

integer, intent(out) :: localComm

integer :: world_np, world_rank
integer :: remoteLeader
integer :: memberKey
integer :: share_int = 0
integer :: dummy

character (*), parameter :: sub_name = mod_name // '.create_mpi_comms'

! Get number of processors in world comm
call mpi_comm_size(MPI_COMM_WORLD, world_np, ierr)
call mpi_comm_rank(MPI_COMM_WORLD, world_rank, ierr)

! Sum together the use_cps values
if (use_cps) share_int = 1
    call MPI_Allreduce(share_int, dummy, 1, MPI_INTEGER, MPI_SUM,              &
    MPI_COMM_WORLD, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    share_int = dummy
if (share_int == 0) then
    ! Do not use CPS
    localComm = MPI_COMM_WORLD
else if (share_int == world_np) then
    ! Use CPS
    write(*,*) "Using Concurrent Precursor method"

    ! Set color
    if (inflow_cond == 5) then
        color = BLUE
    else
        color = RED
    end if

    ! Make sure there are exactly half with inflow conditions
    call MPI_Allreduce(color, share_int, 1, MPI_INTEGER, MPI_SUM,              &
        MPI_COMM_WORLD, ierr)
    if (share_int /= world_np / 2) then
        call error(sub_name ,'There must be one process with inflow conditions')
    end if

    ! set the remote leader for each intercommunicator interComm
    share_int = color
    call MPI_Bcast(share_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (color == share_int) then
        remoteLeader = 0
    else
        remoteLeader = world_np / 2
    end if

    ! Generate member key
    memberKey=modulo(world_rank, world_np / 2)

    ! Split the world communicator into intracommunicators localComm
    call MPI_Comm_split(MPI_COMM_WORLD, color, memberKey, localComm, ierr)

    ! Create intercommunicator interComm
    call mpi_intercomm_create( localComm, 0, MPI_COMM_WORLD, remoteLeader, 1,  &
        interComm, ierr)
else
    call error(sub_name, 'Inconsistent use of use_cps')
end if

end subroutine create_mpi_comms

!*******************************************************************************
subroutine cps_init()
!*******************************************************************************
use param, only : nx, ny, nz
use param, only : coord, rank_of_coord, status, ierr, recycl_region_end
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.cps_init'

integer :: fr_num

if (color == BLUE) then
    ! The fringe region should already be created, but we have to send the size
    ! of the sample block to upstream domain (RED)
    fr_num = fringe%n
    call mpi_send(fr_num, 1, MPI_INTEGER, &
        rank_of_coord(coord), 1, interComm, ierr )
else if (color == RED) then
   ! Receive from downstream domain (BLUE) the length of the sample block
   call mpi_recv(fr_num, 1, MPI_INTEGER, &
        rank_of_coord(coord), 1, interComm, status, ierr)
   recycl = fringe_region_t(recycl_region_end, fr_num)
else
    call error(sub_name,'Erroneous color specification')
endif

end subroutine cps_init

!*******************************************************************************
subroutine synchronize_cps()
!*******************************************************************************
use types, only : rprec
use messages
use param, only : ny, nz
use param, only : coord, rank_of_coord, status, ierr, MPI_RPREC
implicit none

character (*), parameter :: sub_name = mod_name // '.synchronize_cps'

integer :: size

if( color == BLUE ) then
    ! Recieve sampled velocities from upstream (RED)
    size = fringe%n * ny * nz
    call mpi_recv(fringe%u(1,1,1), size, MPI_RPREC, &
        rank_of_coord(coord), 1, interComm, status, ierr)
    call mpi_recv(fringe%v(1,1,1), size, MPI_RPREC, &
        rank_of_coord(coord), 2, interComm, status, ierr)
    call mpi_recv(fringe%w(1,1,1), size, MPI_RPREC, &
        rank_of_coord(coord), 3, interComm, status, ierr)
elseif( color == RED ) then
    ! Send sampled velocities to downstream domain (BLUE)
    call recycl%sample_vel()
    size = recycl%n * ny * nz
    call mpi_send(recycl%u(1,1,1), size, MPI_RPREC, &
        rank_of_coord(coord), 1, interComm, ierr )
    call mpi_send(recycl%v(1,1,1), size, MPI_RPREC, &
        rank_of_coord(coord), 2, interComm, ierr )
    call mpi_send(recycl%w(1,1,1), size, MPI_RPREC, &
        rank_of_coord(coord), 3, interComm, ierr )
else
    call error( sub_name, 'Erroneous color specification')
endif

end subroutine synchronize_cps

end module concurrent_precursor
