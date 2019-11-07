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
use mpi_defs
implicit none

save
private

public :: vel_sample_t
public :: initialize_cps, synchronize_cps, inflow_cond_cps

character(*), parameter :: mod_name = 'concurrent_precursor'

type vel_sample_type
    integer :: nx
    integer :: istart
    integer :: iplateau
    integer :: iend
    integer, allocatable, dimension(:) :: iwrap
    real(rprec), allocatable, dimension(:,:,:) :: u, v, w
#ifdef PPSCALARS
    real(rprec), allocatable, dimension(:,:,:) :: theta
#endif
end type vel_sample_type

type(vel_sample_type), target :: vel_sample_t

! Weights used in fringe region
real(rprec), allocatable, dimension(:) :: alpha, beta

contains

!*******************************************************************************
subroutine initialize_cps()
!*******************************************************************************
use param, only : nx, ny, nz
use param, only : coord, rank_of_coord, status, ierr
use messages
use mpi
use fringe_util, only : fringe_init, fringe_weighting
implicit none

character (*), parameter :: sub_name = mod_name // '.initialize_cps'
integer :: i, index
integer, pointer :: nx_p, istart_p, iplateau_p, iend_p

nullify( nx_p, istart_p, iplateau_p, iend_p )

istart_p   => vel_sample_t % istart
iplateau_p => vel_sample_t % iplateau
iend_p     => vel_sample_t % iend
nx_p       => vel_sample_t % nx

if( color == BLUE ) then
    call fringe_init( istart_p, iplateau_p, iend_p )

    ! Sample size same as buffer region (omitting istart from the block
    ! since velocity is already set there)
    nx_p = iend_p - istart_p

    ! Send size of the sample block to upstream domain (RED)
    call mpi_send( nx_p , 1, MPI_INTEGER, &
    rank_of_coord(coord), 1, interComm, ierr )

    ! Now setup fringe weights
    allocate( alpha( nx_p ), beta( nx_p ) )
    index = 0
    do i = istart_p + 1, iend_p
        index=index+1
        beta(index) = fringe_weighting( i, istart_p, iplateau_p )
    enddo
    alpha = 1.0_rprec - beta
elseif( color == RED ) then
    ! Receive from downstream domain (BLUE) the length of the sample block
    call mpi_recv( nx_p , 1, MPI_INTEGER, &
    rank_of_coord(coord), 1, interComm, status, ierr)

    ! Should end up as nx + 1 (this eventually gets wrapped)
    iend_p = nx + 1
    ! Plateau location not used since no fringe treatment on the RED domain,
    ! but setting so it is at least initialized.
    iplateau_p = iend_p
    ! Set istart based on the size of the sample block
    istart_p = iend_p - nx_p
else
    call error(sub_name,'Erroneous color specification')
endif

! Allocate and assign wrapped index and fringe weights
allocate( vel_sample_t % iwrap( nx_p ) )
index = 0
do i = istart_p + 1, iend_p
    index=index+1
    vel_sample_t % iwrap(index) = modulo( i - 1, nx ) + 1
enddo

! Allocate the sample block
allocate( vel_sample_t % u( nx_p, ny, nz ) )
allocate( vel_sample_t % v( nx_p, ny, nz ) )
allocate( vel_sample_t % w( nx_p, ny, nz ) )
#ifdef PPSCALARS
allocate( vel_sample_t % theta( nx_p, ny, nz ) )
#endif

nullify( nx_p, istart_p, iplateau_p, iend_p )

end subroutine initialize_cps

!*******************************************************************************
subroutine synchronize_cps()
!*******************************************************************************
use types, only : rprec
use messages
use param, only : ny, nz
use param, only : coord, rank_of_coord, status, ierr, MPI_RPREC
use sim_param, only : u,v,w
#ifdef PPSCALARS
use scalars, only : theta
#endif
use coriolis, only : coriolis_forcing, alpha, G
implicit none

character (*), parameter :: sub_name = mod_name // '.synchronize_cps'
integer, pointer :: nx_p
integer, pointer, dimension(:) :: iwrap_p
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
#ifdef PPSCALARS
real(rprec), pointer, dimension(:,:,:) :: theta_p
#endif
integer :: sendsize, recvsize

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif
nullify( nx_p, iwrap_p )

iwrap_p  => vel_sample_t % iwrap
nx_p     => vel_sample_t % nx
u_p      => vel_sample_t % u
v_p      => vel_sample_t % v
w_p      => vel_sample_t % w
#ifdef PPSCALARS
theta_p      => vel_sample_t % theta
#endif

sendsize = nx_p * ny * nz
recvsize = sendsize

if( color == BLUE ) then
    ! Recieve sampled velocities from upstream (RED)
    call mpi_recv( u_p(1,1,1) , recvsize, MPI_RPREC,                           &
        rank_of_coord(coord), 1, interComm, status, ierr)
    call mpi_recv( v_p(1,1,1) , recvsize, MPI_RPREC,                           &
        rank_of_coord(coord), 2, interComm, status, ierr)
    call mpi_recv( w_p(1,1,1) , recvsize, MPI_RPREC,                           &
        rank_of_coord(coord), 3, interComm, status, ierr)
#ifdef PPSCALARS
    call mpi_recv( theta_p(1,1,1) , recvsize, MPI_RPREC,                       &
        rank_of_coord(coord), 4, interComm, status, ierr)
#endif
if (coriolis_forcing>0) then
    call mpi_recv( G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm, status, ierr )
    call mpi_recv( alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm, status, ierr )
end if

elseif( color == RED ) then
    ! Sample velocity and copy to buffers
    u_p(:,:,:) = u(iwrap_p(:),1:ny,1:nz)
    v_p(:,:,:) = v(iwrap_p(:),1:ny,1:nz)
    w_p(:,:,:) = w(iwrap_p(:),1:ny,1:nz)
#ifdef PPSCALARS
    theta_p(:,:,:) = theta(iwrap_p(:),1:ny,1:nz)
#endif

    ! Send sampled velocities to downstream domain (BLUE)
    call mpi_send( u_p(1,1,1), sendsize, MPI_RPREC,                            &
        rank_of_coord(coord), 1, interComm, ierr )
    call mpi_send( v_p(1,1,1), sendsize, MPI_RPREC,                            &
        rank_of_coord(coord), 2, interComm, ierr )
    call mpi_send( w_p(1,1,1), sendsize, MPI_RPREC,                            &
        rank_of_coord(coord), 3, interComm, ierr )
#ifdef PPSCALARS
    call mpi_send( theta_p(1,1,1), sendsize, MPI_RPREC,                        &
        rank_of_coord(coord), 4, interComm, ierr )
#endif
if (coriolis_forcing>0) then
    call mpi_send( G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm, ierr )
    call mpi_send( alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm, ierr )
end if

else
   call error( sub_name, 'Erroneous color specification')
endif

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif
nullify( nx_p, iwrap_p )

end subroutine synchronize_cps

!*******************************************************************************
subroutine inflow_cond_cps ()
!*******************************************************************************
!
!  Enforces prescribed inflow condition from an inlet velocity field
!  generated from a precursor simulation. The inflow condition is
!  enforced by direct modulation on the velocity in the fringe region.
!
use types, only : rprec
use param, only : nx, ny, nz
use sim_param, only : u, v, w
#ifdef PPSCALARS
use scalars, only : theta
#endif
use messages, only : error
implicit none

character (*), parameter :: sub_name = 'inflow_cond_cps'
integer :: j,k
integer :: istart_wrap
integer, pointer :: istart_p
integer, pointer, dimension(:) :: iwrap_p
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
#ifdef PPSCALARS
real(rprec), pointer, dimension(:,:,:) :: theta_p
#endif

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif
nullify( istart_p, iwrap_p )

u_p        => vel_sample_t % u
v_p        => vel_sample_t % v
w_p        => vel_sample_t % w
#ifdef PPSCALARS
theta_p        => vel_sample_t % theta
#endif
istart_p   => vel_sample_t % istart
iwrap_p    => vel_sample_t % iwrap

istart_wrap = modulo( istart_p - 1, nx ) + 1

do k = 1, nz
do j = 1, ny
    u(iwrap_p(:),j,k) = alpha(:) * u(istart_wrap,j,k) + beta(:) * u_p(:,j,k)
    v(iwrap_p(:),j,k) = alpha(:) * v(istart_wrap,j,k) + beta(:) * v_p(:,j,k)
    w(iwrap_p(:),j,k) = alpha(:) * w(istart_wrap,j,k) + beta(:) * w_p(:,j,k)
#ifdef PPSCALARS
    theta(iwrap_p(:),j,k) = alpha(:) * theta(istart_wrap,j,k) + beta(:) * theta_p(:,j,k)
#endif
enddo
enddo

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif
nullify( istart_p, iwrap_p )

end subroutine inflow_cond_cps

end module concurrent_precursor
