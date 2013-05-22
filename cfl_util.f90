!!
!!  Copyright 2010,2011,2012 Johns Hopkins University
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
