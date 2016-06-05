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
use types, only : rprec
$if($MPI)
use mpi
$endif
save
$if ($MPI)
!include 'mpif.h'
$endif
private

public get_max_cfl, get_cfl_dt

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function get_max_cfl() result(cfl)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This function provides the value of the maximum CFL in the entire 
! domain
use param, only : dt, dx, dy, dz, nx, ny, nz, fourier, nxp
use sim_param, only : u, v, w
use sim_param, only : uF,vF,wF
$if($MPI)
use param, only : ierr, MPI_RPREC
$endif
implicit none
real(rprec) :: cfl,cfl_u,cfl_v,cfl_w

$if($MPI)
real(rprec) :: cfl_buf
$endif

if (fourier) then   !! remember dx = L_x / nxp (if fourier=true)
   cfl_u = maxval( abs(uF(1:nxp,1:ny,1:nz-1)) ) / dx
   cfl_v = maxval( abs(vF(1:nxp,1:ny,1:nz-1)) ) / dy
   cfl_w = maxval( abs(wF(1:nxp,1:ny,1:nz-1)) ) / dz
else
   cfl_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
   cfl_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
   cfl_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz
endif

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
use param, only : cfl, dx, dy, dz, nx, ny, nz, fourier, nxp
use sim_param, only : u,  v,  w
use sim_param, only : uF, vF, wF
$if($MPI)
use param, only : ierr, MPI_RPREC
$endif
implicit none
real(rprec) :: dt,dt_inv_u,dt_inv_v,dt_inv_w

$if($MPI)
real(rprec) :: dt_buf
$endif

! Avoid division by computing max dt^-1
if (fourier) then
   dt_inv_u = maxval( abs(uF(1:nxp,1:ny,1:nz-1)) ) / dx
   dt_inv_v = maxval( abs(vF(1:nxp,1:ny,1:nz-1)) ) / dy 
   dt_inv_w = maxval( abs(wF(1:nxp,1:ny,1:nz-1)) ) / dz
else
   dt_inv_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
   dt_inv_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy 
   dt_inv_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz
endif

dt = cfl / maxval( (/ dt_inv_u, dt_inv_v, dt_inv_w /) )

$if($MPI)
  call mpi_allreduce(dt, dt_buf, 1, MPI_RPREC, MPI_MIN, MPI_COMM_WORLD, ierr)
  dt = dt_buf
$endif

return
end function get_cfl_dt

end module cfl_util
