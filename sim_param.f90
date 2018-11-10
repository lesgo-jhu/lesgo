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

!*******************************************************************************
module sim_param
!*******************************************************************************
use types, only : rprec
use param, only : ld, nx, ny, nz, grid
use sim_var_3d_m
implicit none

save
public

logical :: sim_param_initialized = .false.

real(rprec), dimension(:,:,:), allocatable ::RHSx, RHSy, RHSz, RHSx_f, RHSy_f, &
    RHSz_f, dpdx, dpdy, dpdz, fx, fy, fz, fxa, fya, fza
real(rprec), target, dimension(:,:,:), allocatable :: p

real(rprec), dimension(:,:,:), pointer :: u, dudx, dudy, dudz
real(rprec), dimension(:,:,:), pointer :: v, dvdx, dvdy, dvdz
real(rprec), dimension(:,:,:), pointer :: w, dwdx, dwdy, dwdz
real(rprec), dimension(:,:,:), pointer :: txx, txy, tyy, txz, tyz, tzz
real(rprec), dimension(:,:,:), pointer :: divtx, divty, divtz

type(sim_var_3d_t) :: u_var, dudx_var, dudy_var, dudz_var, u_filt_var
type(sim_var_3d_t) :: v_var, dvdx_var, dvdy_var, dvdz_var, v_filt_var
type(sim_var_3d_t) :: w_var, dwdx_var, dwdy_var, dwdz_var, w_filt_var
type(sim_var_3d_t) :: txx_var, txy_var, tyy_var, txz_var, tyz_var, tzz_var
type(sim_var_3d_t) :: divtx_var, divty_var, divtz_var
type(sim_var_3d_t) :: dtxxdx, dtxydy, dtxzdz
type(sim_var_3d_t) :: dtxydx, dtyydy, dtyzdz
type(sim_var_3d_t) :: dtxzdx, dtyzdy, dtzzdz

contains

!*******************************************************************************
subroutine sim_param_init ()
!*******************************************************************************
!
! This subroutine initilizes all global arrays defined in the sim_param
! module. Here they are allocated and initialized to zero.
!

implicit none

u_var = sim_var_3d_t(grid, UV_GRID)
dudx_var = sim_var_3d_t(grid, UV_GRID)
dudy_var = sim_var_3d_t(grid, UV_GRID)
dudz_var = sim_var_3d_t(grid, W_GRID)
u_filt_var = sim_var_3d_t(grid, UV_GRID)
u => u_var%real
dudx => dudx_var%real
dudy => dudy_var%real
dudz => dudz_var%real

v_var = sim_var_3d_t(grid, UV_GRID)
dvdx_var = sim_var_3d_t(grid, UV_GRID)
dvdy_var = sim_var_3d_t(grid, UV_GRID)
dvdz_var = sim_var_3d_t(grid, W_GRID)
v_filt_var = sim_var_3d_t(grid, UV_GRID)
v => v_var%real
dvdx => dvdx_var%real
dvdy => dvdy_var%real
dvdz => dvdz_var%real

w_var = sim_var_3d_t(grid, W_GRID)
dwdx_var = sim_var_3d_t(grid, W_GRID)
dwdy_var = sim_var_3d_t(grid, W_GRID)
dwdz_var = sim_var_3d_t(grid, UV_GRID)
w_filt_var = sim_var_3d_t(grid, W_GRID)
w => w_var%real
dwdx => dwdx_var%real
dwdy => dwdy_var%real
dwdz => dwdz_var%real

txx_var = sim_var_3d_t(grid, UV_GRID)
txy_var = sim_var_3d_t(grid, UV_GRID)
tyy_var = sim_var_3d_t(grid, UV_GRID)
txz_var = sim_var_3d_t(grid, W_GRID)
tyz_var = sim_var_3d_t(grid, W_GRID)
tzz_var = sim_var_3d_t(grid, UV_GRID)
txx => txx_var%real
txy => txy_var%real
tyy => tyy_var%real
txz => txz_var%real
tyz => tyz_var%real
tzz => tzz_var%real

dtxxdx = sim_var_3d_t(grid, UV_GRID)
dtxydy = sim_var_3d_t(grid, UV_GRID)
dtxzdz = sim_var_3d_t(grid, UV_GRID)
dtxydx = sim_var_3d_t(grid, UV_GRID)
dtyydy = sim_var_3d_t(grid, UV_GRID)
dtyzdz = sim_var_3d_t(grid, UV_GRID)
dtxzdx = sim_var_3d_t(grid, W_GRID)
dtyzdy = sim_var_3d_t(grid, W_GRID)
dtzzdz = sim_var_3d_t(grid, W_GRID)

divtx_var = sim_var_3d_t(grid, UV_GRID)
divty_var = sim_var_3d_t(grid, UV_GRID)
divtz_var = sim_var_3d_t(grid, W_GRID)
divtx => divtx_var%real
divty => divty_var%real
divtz => divtz_var%real

! allocate ( u(ld, ny, 0:nz) ); u = 0.0_rprec
! allocate ( v(ld, ny, 0:nz) ); v = 0.0_rprec
! allocate ( w(ld, ny, 0:nz) ); w = 0.0_rprec
! allocate( dudx(ld, ny, 0:nz) ); dudx = 0.0_rprec
! allocate( dudy(ld, ny, 0:nz) ); dudy = 0.0_rprec
! allocate( dudz(ld, ny, 0:nz) ); dudz = 0.0_rprec
! allocate( dvdx(ld, ny, 0:nz) ); dvdx = 0.0_rprec
! allocate( dvdy(ld, ny, 0:nz) ); dvdy = 0.0_rprec
! allocate( dvdz(ld, ny, 0:nz) ); dvdz = 0.0_rprec
! allocate( dwdx(ld, ny, 0:nz) ); dwdx = 0.0_rprec
! allocate( dwdy(ld, ny, 0:nz) ); dwdy = 0.0_rprec
! allocate( dwdz(ld, ny, 0:nz) ); dwdz = 0.0_rprec
allocate( RHSx(ld, ny, 0:nz) ); RHSx = 0.0_rprec
allocate( RHSy(ld, ny, 0:nz) ); RHSy = 0.0_rprec
allocate( RHSz(ld, ny, 0:nz) ); RHSz = 0.0_rprec
allocate( RHSx_f(ld, ny, 0:nz) ); RHSx_f = 0.0_rprec
allocate( RHSy_f(ld, ny, 0:nz) ); RHSy_f = 0.0_rprec
allocate( RHSz_f(ld, ny, 0:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
! allocate ( txx(ld, ny, 0:nz) ); txx = 0.0_rprec
! allocate ( txy(ld, ny, 0:nz) ); txy = 0.0_rprec
! allocate ( tyy(ld, ny, 0:nz) ); tyy = 0.0_rprec
! allocate ( txz(ld, ny, 0:nz) ); txz = 0.0_rprec
! allocate ( tyz(ld, ny, 0:nz) ); tyz = 0.0_rprec
! allocate ( tzz(ld, ny, 0:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
! allocate ( divtx(ld, ny, 0:nz) ); divtx = 0.0_rprec
! allocate ( divty(ld, ny, 0:nz) ); divty = 0.0_rprec
! allocate ( divtz(ld, ny, 0:nz) ); divtz = 0.0_rprec

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
allocate ( fxa(ld, ny, 0:nz) ); fxa = 0.0_rprec
allocate ( fya(ld, ny, 0:nz) ); fya = 0.0_rprec
allocate ( fza(ld, ny, 0:nz) ); fza = 0.0_rprec
#endif

#if defined(PPLVLSET) || defined(PPATM)
allocate ( fx(ld, ny, nz) ); fx = 0.0_rprec
allocate ( fy(ld, ny, nz) ); fy = 0.0_rprec
allocate ( fz(ld, ny, nz) ); fz = 0.0_rprec
#endif

sim_param_initialized = .true.

end subroutine sim_param_init

end module sim_param
