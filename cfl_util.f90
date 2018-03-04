!!
!!  Copyright (C) 2010-2013  Johns Hopkins University
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

!**********************************************************************
module cfl_util
!**********************************************************************
!
! This module provides the subroutines/functions for getting CFL related
! quantities 
!
save 
private

public get_max_cfl, get_cfl_dt

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function get_max_cfl() result(cfl)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This function provides the value of the maximum CFL in the entire 
! domain
!
use types, only : rprec
use param, only : dt, dx, dy, dz, nx, ny, nz
use sim_param, only : u,v,w

$if($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif

implicit none
real(rprec) :: cfl

real(rprec) :: cfl_u, cfl_v, cfl_w

$if($MPI)
real(rprec) :: cfl_buf
$endif

cfl_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
cfl_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
cfl_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz

cfl = dt * maxval( (/ cfl_u, cfl_v, cfl_w /) )

$if($MPI)
  call mpi_allreduce(cfl, cfl_buf, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
  cfl = cfl_buf
$endif

return
end function get_max_cfl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function get_cfl_dt() result(dt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This functions determines the maximum allowable time step based on the CFL
! value specified in the param module
!
use types, only : rprec
use param, only : cfl, dx, dy, dz, nx, ny, nz
use sim_param, only : u,v,w

$if($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif

implicit none

real(rprec) :: dt

! dt inverse
real(rprec) :: dt_inv_u, dt_inv_v, dt_inv_w

$if($MPI)
real(rprec) :: dt_buf
$endif

! Avoid division by computing max dt^-1
dt_inv_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
dt_inv_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy 
dt_inv_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz

dt = cfl / maxval( (/ dt_inv_u, dt_inv_v, dt_inv_w /) )

$if($MPI)
  call mpi_allreduce(dt, dt_buf, 1, MPI_RPREC, MPI_MIN, MPI_COMM_WORLD, ierr)
  dt = dt_buf
$endif


return
end function get_cfl_dt

end module cfl_util
