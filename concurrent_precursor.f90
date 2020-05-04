!!
!!  Copyright (C) 2011-2020  Johns Hopkins University
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
use fringe
implicit none

save
private

public :: vel_sample_t
public :: initialize_cps, synchronize_cps, inflow_cps

character(*), parameter :: mod_name = 'concurrent_precursor'

type vel_sample_type
    real(rprec), allocatable, dimension(:,:,:) :: u, v, w
#ifdef PPSCALARS
    real(rprec), allocatable, dimension(:,:,:) :: theta
#endif
end type vel_sample_type

type(vel_sample_type), target :: vel_sample_t
type(fringe_t) :: cps_fringe

contains

!*******************************************************************************
subroutine initialize_cps()
!*******************************************************************************
use param, only : nx, ny, nz, dx, L_x, coord, rank_of_coord, status, ierr
use param, only : fringe_region_end, fringe_region_len, sampling_region_end
use messages
use mpi
implicit none

character (*), parameter :: sub_name = mod_name // '.initialize_cps'
integer :: i

if( color == BLUE ) then
    cps_fringe = fringe_t(fringe_region_end, fringe_region_len)
    call mpi_send(cps_fringe%nx , 1, MPI_INTEGER,                              &
        rank_of_coord(coord), 1, interComm, ierr )
elseif( color == RED) then
    call mpi_recv(i , 1, MPI_INTEGER,                                          &
        rank_of_coord(coord), 1, interComm, status, ierr)
    fringe_region_len = (i - 0.5_rprec)*dx/L_x
    cps_fringe = fringe_t(sampling_region_end, fringe_region_len)
else
   call error(sub_name, 'Erroneous color specification')
endif

! Allocate the sample block
allocate(vel_sample_t%u(cps_fringe%nx, ny, nz ))
allocate(vel_sample_t%v(cps_fringe%nx, ny, nz ))
allocate(vel_sample_t%w(cps_fringe%nx, ny, nz ))
#ifdef PPSCALARS
allocate(vel_sample_t%theta(cps_fringe%nx, ny, nz))
#endif

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
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
#ifdef PPSCALARS
real(rprec), pointer, dimension(:,:,:) :: theta_p
#endif
integer :: sendsize, recvsize

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif

u_p => vel_sample_t%u
v_p => vel_sample_t%v
w_p => vel_sample_t%w
#ifdef PPSCALARS
theta_p => vel_sample_t%theta
#endif

sendsize = cps_fringe%nx * ny * nz
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
    call mpi_recv(G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm,         &
        status, ierr)
    call mpi_recv(alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm,     &
        status, ierr)
end if

elseif( color == RED ) then
    ! Sample velocity and copy to buffers
    u_p(:,:,:) = u(cps_fringe%iwrap(:),1:ny,1:nz)
    v_p(:,:,:) = v(cps_fringe%iwrap(:),1:ny,1:nz)
    w_p(:,:,:) = w(cps_fringe%iwrap(:),1:ny,1:nz)
#ifdef PPSCALARS
    theta_p(:,:,:) = theta(cps_fringe%iwrap(:),1:ny,1:nz)
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
    call mpi_send(G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm, ierr)
    call mpi_send(alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm, ierr)
end if

else
   call error(sub_name, 'Erroneous color specification')
endif

nullify(u_p, v_p, w_p)
#ifdef PPSCALARS
nullify(theta_p)
#endif

end subroutine synchronize_cps

!*******************************************************************************
subroutine inflow_cps ()
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
integer :: i, i_w
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
#ifdef PPSCALARS
real(rprec), pointer, dimension(:,:,:) :: theta_p
#endif

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif

u_p => vel_sample_t%u
v_p => vel_sample_t%v
w_p => vel_sample_t%w
#ifdef PPSCALARS
theta_p => vel_sample_t%theta
#endif

do i = 1, cps_fringe%nx
    i_w = cps_fringe%iwrap(i)
    u(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * u(i_w,1:ny,1:nz)                  &
        + cps_fringe%beta(i) * u_p(i,1:ny,1:nz)
    v(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * v(i_w,1:ny,1:nz)                  &
        + cps_fringe%beta(i) * v_p(i,1:ny,1:nz)
    w(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * w(i_w,1:ny,1:nz)                  &
        + cps_fringe%beta(i) * w_p(i,1:ny,1:nz)
#ifdef PPSCALARS
    theta(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * theta(i_w,1:ny,1:nz)          &
        + cps_fringe%beta(i) * theta_p(i,1:ny,1:nz)
#endif
end do

nullify(u_p, v_p, w_p)
#ifdef PPSCALARS
nullify(theta_p)
#endif

end subroutine inflow_cps

end module concurrent_precursor
