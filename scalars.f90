!!
!!  Copyright (C) 2019  Johns Hopkins University
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
module scalars
!*******************************************************************************
! This module contains all of the subroutines associated with scalar transport
use types, only : rprec
use param, only : path
implicit none

save
private

public :: scalars_init, ic_scal, buoyancy_force, scalars_transport,            &
    scalars_checkpoint, obukhov, scalars_deriv

real(rprec), public, dimension(:,:,:), allocatable :: theta, dTdx, dTdy, dTdz, &
    RHS_T, RHS_Tf, u_big, v_big, w_big, dTdx_big, dTdy_big, dTdz_big, RHS_big, &
    pi_x, pi_y, pi_z, div_pi, temp_var

! SGS values
real(rprec), public :: Pr_sgs = 0.5
real(rprec), public, dimension(:,:,:), allocatable :: Kappa_t

! Monin-Obukhov BC
real(rprec), public, dimension(:,:), allocatable :: psi_m, phi_m, psi_h, phi_h,&
    L, tstar_lbc

! Gravitational acceleration
real(rprec), public :: g = 9.81_rprec
! Roughness length for scalars. typically zo/10
real(rprec), public :: zo_s = 0.00001_rprec
! Treat theta as passive_scalar (no buoyancy)
logical, public :: passive_scalar = .false.
! Whether to initialize theta field
logical, public :: inits = .true.
! Name of file for restarting
character(64) :: fname
! Reference temperature scale
real(rprec), public :: T_scale = 300._rprec

! Boundary conditions
! lbc: lower boundary condition
! ubc: upper boundary condition
!       0 - prescribed temperature, 1 - prescribed flux
integer, public :: lbc_scal = 0
real(rprec), public :: scal_bot = 300._rprec
real(rprec), public :: flux_bot = 0._rprec
real(rprec), dimension(:), allocatable, public :: ic_z, ic_theta
real(rprec), public :: ic_no_vel_noise_z
integer, public :: ic_nloc
real(rprec), public :: lapse_rate = 0._rprec
logical, public :: read_lbc_scal = .false.

! Interpolation of bottom boundary condition
real(rprec), dimension(:), allocatable :: t_interp, lbc_interp

contains

!*******************************************************************************
subroutine scalars_init
!*******************************************************************************
! This subroutine initializes the variables for the scalars module
use param, only : lbz, ld, ld_big, nx, ny, nz, ny2, u_star, z_i
use functions, only : count_lines
integer :: i, num_t, fid

! Allocate simulation variables
allocate ( theta(ld, ny, lbz:nz) ); theta = 0._rprec
allocate ( dTdx(ld, ny, lbz:nz) ); dTdx = 0._rprec
allocate ( dTdy(ld, ny, lbz:nz) ); dTdy = 0._rprec
allocate ( dTdz(ld, ny, lbz:nz) ); dTdz = 0._rprec
allocate ( RHS_T(ld, ny, lbz:nz) ); RHS_T = 0._rprec
allocate ( RHS_Tf(ld, ny, lbz:nz) ); RHS_Tf = 0._rprec
allocate ( u_big(ld_big, ny2, lbz:nz)); u_big = 0._rprec
allocate ( v_big(ld_big, ny2, lbz:nz)); v_big = 0._rprec
allocate ( w_big(ld_big, ny2, lbz:nz)); w_big = 0._rprec
allocate ( dTdx_big(ld_big, ny2, lbz:nz)); dTdx_big = 0._rprec
allocate ( dTdy_big(ld_big, ny2, lbz:nz)); dTdy_big = 0._rprec
allocate ( dTdz_big(ld_big, ny2, lbz:nz)); dTdz_big = 0._rprec
allocate ( RHS_big(ld_big, ny2, lbz:nz)); RHS_big = 0._rprec
allocate ( pi_x(ld, ny, lbz:nz) ); pi_x = 0._rprec
allocate ( pi_y(ld, ny, lbz:nz) ); pi_y = 0._rprec
allocate ( pi_z(ld, ny, lbz:nz) ); pi_z = 0._rprec
allocate ( div_pi(ld, ny, lbz:nz) ); div_pi = 0._rprec
allocate ( temp_var(ld, ny, lbz:nz) ); temp_var = 0._rprec
allocate ( Kappa_t(ld, ny, lbz:nz) ); Kappa_t = 0._rprec

! Obukhov values (defaults for passive scalars)
allocate ( psi_m(nx, ny) ); psi_m = 0._rprec
allocate ( phi_m(nx, ny) ); phi_m = 1._rprec
allocate ( psi_h(nx, ny) ); psi_h = 0._rprec
allocate ( phi_h(nx, ny) ); phi_h = 1._rprec
allocate ( L(nx, ny) ); L = 0._rprec
allocate ( tstar_lbc(nx, ny) ); tstar_lbc = 0._rprec

! Nondimensionalize variables
g = g*(z_i/(u_star**2))
flux_bot = flux_bot/u_star/T_scale
scal_bot = scal_bot/T_scale
lapse_rate = lapse_rate/T_scale*z_i
ic_theta = ic_theta/T_scale
ic_z = ic_z/z_i
ic_no_vel_noise_z = ic_no_vel_noise_z/z_i

! Read values from file
if (read_lbc_scal) then
    ! Count number of entries and allocate
    num_t = count_lines('lbc_scal.dat')
    allocate( t_interp(num_t) )
    allocate( lbc_interp(num_t) )

    ! Read entries
    open(newunit=fid, file='lbc_scal.dat', status='unknown', form='formatted', &
        position='rewind')
    do i = 1, num_t
        read(fid,*) t_interp(i), lbc_interp(i)
    end do
end if

end subroutine scalars_init

!*******************************************************************************
subroutine ic_scal(interp_flag)
!*******************************************************************************
! Set initial profile for scalar
use param, only : coord
use string_util
use grid_m
logical :: interp_flag

fname = path // 'scal.out'
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif
inquire (file=fname, exist=inits)

if (inits .and. .not.interp_flag) then
    inits = .true.
else
    inits = .false.
end if


if (inits) then
    write(*,*)  "--> Reading initial scalar field from file"
    call ic_scal_file
elseif (interp_flag) then
    write(*,*)  "--> Interpolating initial scalar field from file"
    call ic_scal_interp
else
    write(*,*)  "--> Creating initial boundary layer scalar field with LES BCs"
    call ic_scal_les
end if

end subroutine ic_scal

!*******************************************************************************
subroutine ic_scal_file
!*******************************************************************************
! Read initial profile for scalar from file
use param, only : nx, nz, read_endian
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP
use grid_m

open(12, file=fname, form='unformatted', convert=read_endian)
read(12) theta(:, :, 1:nz), RHS_T(:, :, 1:nz), psi_m(1:nx, :)
close(12)

#ifdef PPMPI
call mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(RHS_T, 0, MPI_SYNC_DOWNUP)
#endif

end subroutine ic_scal_file

!*******************************************************************************
subroutine ic_scal_les
!*******************************************************************************
use param, only : lbz, nz
use grid_m
use functions, only : linear_interp

integer :: i

do i = lbz, nz
    if (grid%z(i) < ic_z(ic_nloc)) then
        theta(:,:,i) = linear_interp(ic_z, ic_theta, grid%z(i))
    else
        theta(:,:,i) = ic_theta(ic_nloc) + lapse_rate*(grid%z(i) - ic_z(ic_nloc))
    end if
end do

end subroutine ic_scal_les

!*******************************************************************************
subroutine ic_scal_interp()
!*******************************************************************************
! This subroutine reads the initial conditions from a checkpoint file and
! interpolates onto the current grid
!
use param, only : nx, ny, nz, lbz, read_endian, path
use grid_m
use functions
integer :: nproc_f, Nx_f, Ny_f, Nz_f
real(rprec) :: Lx_f, Ly_f, Lz_f
integer :: i, j, k, z1, z2, ld_f, lh_f, Nz_tot_f
real(rprec) :: dx_f, dy_f, dz_f
integer :: i1, i2, j1, j2, k1, k2
real(rprec) :: ax, ay, az, bx, by, bz, xx, yy
real(rprec), allocatable, dimension(:) :: x_f, y_f!, z_f, zw_f
real(rprec), allocatable, dimension(:,:,:) :: theta_f
character(64) :: ff
integer :: npr1, npr2, nproc_r, Nz_tot_r
real(rprec), allocatable, dimension(:) :: z_r, zw_r
! real(rprec), allocatable, dimension(:,:,:) :: u_r, v_r, w_r

! Read grid information from file
open(12, file=path // 'grid.out', form='unformatted', convert=read_endian)
read(12) nproc_f, Nx_f, Ny_f, Nz_f, Lx_f, Ly_f, Lz_f
close(12)

! Compute intermediate values
lh_f = nx_f/2+1
ld_f = 2*lh_f
Nz_tot_f = ( nz_f - 1 ) * nproc_f + 1
dx_f = Lx_f / nx_f
dy_f = Ly_f / ny_f
dz_f = Lz_f / (nz_tot_f - 1)

! Figure out which processors actually need to be read
npr1 = max(floor(grid%z(lbz)/dz_f/Nz_f), 0)
npr2 = min(ceiling(grid%z(nz)/dz_f/Nz_f), nproc_f-1)
nproc_r = npr2-npr1+1
Nz_tot_r = ( nz_f - 1 ) * nproc_r + 1
! write(*,*) coord, npr1, npr2, nproc_r, Nz_tot_r

! Create file grid
allocate( z_r(nz_tot_r), zw_r(nz_tot_r))
do i = 1, nz_tot_r
    zw_r(i) = (i - 1 + npr1*(nz_f-1)) * dz_f
    z_r(i) = zw_r(i) + 0.5*dz_f
end do

! Create file grid
allocate( x_f(nx_f), y_f(ny_f))
do i = 1, nx_f
    x_f(i) = (i-1) * Lx_f/(nx_f)
end do
do i = 1, ny_f
    y_f(i) = (i-1) * Ly_f/ny_f
end do

! Read velocities from file
allocate( theta_f(ld_f, ny_f, nz_tot_r) )

! Loop through all levels
do i = 1, nproc_r
    ! Set level bounds
    z1 = nz_tot_f / nproc_f * (i-1) + 1
    z2 = nz_tot_f / nproc_f * i  + 1

    ! Read from file
    write(ff,*) i+npr1-1
    ff = path // "scal.out.c"//trim(adjustl(ff))
    open(12, file=ff,  action='read', form='unformatted')
    read(12) theta_f(:, :, z1:z2)
    close(12)
end do

! Calculate velocities
do i = 1, nx
    xx = grid%x(i) - floor(grid%x(i)/Lx_f) * Lx_f
    i1 = binary_search(x_f, xx)
    i2 = mod(i1, nx_f) + 1
    bx = (xx - x_f(i1)) / dx_f
    ax = 1._rprec - bx
    do j = 1, ny
        yy = grid%y(j) - floor(grid%y(j)/Ly_f) * Ly_f
        j1 = binary_search(y_f, yy)
        j2 = mod(j1, ny_f) + 1
        by = (yy - y_f(j1)) / dy_f
        ay = 1._rprec - by
        do k = 1, nz
            ! for u and v
            k1 = binary_search(z_r, grid%z(k))
            k2 = k1 + 1
            if (k1 == nz_tot_r) then
                theta(i,j,k) = ax*ay*theta_f(i1,j1,k1)                         &
                             + bx*ay*theta_f(i2,j1,k1)                         &
                             + ax*by*theta_f(i1,j2,k1)                         &
                             + bx*by*theta_f(i2,j2,k1)
            else if (k1 == 0) then
                theta(i,j,k) = ax*ay*theta_f(i1,j1,k2)                         &
                             + bx*ay*theta_f(i2,j1,k2)                         &
                             + ax*by*theta_f(i1,j2,k2)                         &
                             + bx*by*theta_f(i2,j2,k2)
            else
                bz = (grid%z(k) - z_r(k1)) / dz_f
                az = 1._rprec - bz
                theta(i,j,k) = ax*ay*az*theta_f(i1,j1,k1)                      &
                             + bx*ay*az*theta_f(i2,j1,k1)                      &
                             + ax*by*az*theta_f(i1,j2,k1)                      &
                             + bx*by*az*theta_f(i2,j2,k1)                      &
                             + ax*ay*bz*theta_f(i1,j1,k2)                      &
                             + bx*ay*bz*theta_f(i2,j1,k2)                      &
                             + ax*by*bz*theta_f(i1,j2,k2)                      &
                             + bx*by*bz*theta_f(i2,j2,k2)
            end if

        end do
    end do
end do

end subroutine ic_scal_interp

!*******************************************************************************
subroutine scalars_checkpoint
!*******************************************************************************
use param, only : nx, nz, write_endian

!  Open scal.out (lun_default in io) for final output
open(11, file=fname, form='unformatted', convert=write_endian,                 &
    status='unknown', position='rewind')
write(11) theta(:, :, 1:nz), RHS_T(:, :, 1:nz), psi_m(1:nx, :)
close(11)

end subroutine scalars_checkpoint

!*******************************************************************************
subroutine scalars_deriv
!*******************************************************************************
use param, only : lbz, nz, coord, nproc
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP
use derivatives, only : filt_da, ddz_uv

! Calculate derivatives of theta
call filt_da(theta, dTdx, dTdy, lbz)
call ddz_uv(theta, dTdz, lbz)

#ifdef PPMPI
call mpi_sync_real_array(dTdz, 0, MPI_SYNC_DOWNUP)
#endif

! Top boundary condition
if (coord == nproc-1) dTdz(:,:,nz) = lapse_rate

end subroutine scalars_deriv

!*******************************************************************************
subroutine obukhov(u_avg)
!*******************************************************************************
use param, only : vonk, dz, zo, nx, ny, ld, u_star, lbz, total_time_dim
use sim_param, only : ustar_lbc
use coriolis, only : repeat_interval
use functions, only : linear_interp
use test_filtermodule

real(rprec), dimension(nx, ny), intent(in) :: u_avg
real(rprec), dimension(ld, ny) :: theta1

integer :: i, j
! Use previous ustar_lbc to compute stability correction
if (passive_scalar) then
    ustar_lbc = u_avg*vonk/log(0.5_rprec*dz/zo)
    return
end if

theta1 = theta(:,:,1)
call test_filter(theta1)

! Using previous time step's psi_m and psi_h to calculate obukhov length and
! stability functions
ustar_lbc = u_avg*vonk/(log(0.5_rprec*dz/zo) + psi_m)
tstar_lbc = (theta1(1:nx,:) - scal_bot)*vonk                                   &
    / (log(0.5_rprec*dz/zo_s) + psi_h)

L = ustar_lbc**2*theta1(1:nx,:)/(vonk*g*tstar_lbc)
do i = 1, nx
    do j = 1, ny
        call stability(L(i,j), zo_s, phi_m(i,j), phi_h(i,j), psi_m(i,j),       &
            psi_h(i,j))
    end do
end do

! Recompute ustar_lbc using new values
ustar_lbc = u_avg*vonk/(log(0.5_rprec*dz/zo) + psi_m)

! Get boundary condition if reading from file
if (read_lbc_scal) then
    if (lbc_scal == 0) then
        scal_bot = linear_interp(t_interp, lbc_interp,                         &
            mod(total_time_dim, repeat_interval))/T_scale
    else
        flux_bot = linear_interp(t_interp, lbc_interp,                         &
            mod(total_time_dim, repeat_interval))/T_scale/u_star
    end if
end if

! Calculate tstar_lbc based on boundary condition
if (lbc_scal == 0) then
    tstar_lbc = (theta1(1:nx,:) - scal_bot)*vonk                               &
        / (log(0.5_rprec*dz/zo_s) + psi_h)
else
    tstar_lbc = -flux_bot/ustar_lbc
end if

! Calculate temperature gradient and flux
dTdz(1:nx,:,1) = tstar_lbc/(vonk*dz*0.5_rprec)*phi_h
pi_z(1:nx,:,1) = -tstar_lbc*ustar_lbc

end subroutine obukhov

!*******************************************************************************
subroutine scalars_transport()
!*******************************************************************************
use param, only : lbz, nx, nz, nx2, ny2, nproc, coord, dt, tadv1, tadv2,       &
    jt_total, dt
use param, only : lbc_mom, ubc_mom, dz
use sim_param, only : u, v, w
use sgs_param, only : Nu_t
use derivatives, only : filt_da, ddx, ddy, ddz_uv, ddz_w
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP
use test_filtermodule
use fft
use messages, only : error

integer :: k, jz_min, jz_max
real(rprec) :: const

! We do not advance the ground nodes, so start at k=2.
! For the MPI case, the means that we start from jz=2
! for coord=0 and jz=1 otherwise.
#ifdef PPMPI
    if (coord == 0) then
        jz_min = 2
    else
        jz_min = 1
    end if
    if (coord == nproc-1) then
        jz_max = nz-2
    else
        jz_max = nz-1
    end if
#else
    jz_max = nz-2
    jz_min = 2
#endif

! Save previous timestep's RHS
RHS_Tf = RHS_T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advective term u_i d_i \theta computed using dealiasing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We could save memory and computional time by saving u_big, v_big, w_big from
! the convec subroutine.
! (put in sim_param, not as a saved variable in the subroutine)

! Set variables onto big domain for multiplication in physical space and
! dealiasing
call to_big(u, u_big)
call to_big(v, v_big)
call to_big(w, w_big)
call to_big(dTdx, dTdx_big)
call to_big(dTdy, dTdy_big)
call to_big(dTdz, dTdz_big)

! Normalization for FFTs
const=1._rprec/(nx2*ny2)

! Interior of domain
do k = jz_min, jz_max
    RHS_big(:,:,k) = const*(u_big(:,:,k)*dTdx_big(:,:,k)                       &
        + v_big(:,:,k)*dTdy_big(:,:,k)                                         &
        + 0.5_rprec*w_big(:,:,k+1)*dTdz_big(:,:,k+1)                           &
        + 0.5_rprec*w_big(:,:,k)*dTdz_big(:,:,k))
end do

! Bottom of domain
if (coord == 0) then
    RHS_big(:,:,1) = const*(u_big(:,:,1)*dTdx_big(:,:,1)                       &
        + v_big(:,:,1)*dTdy_big(:,:,1)                                         &
        + 0.5_rprec*w_big(:,:,2)*dTdz_big(:,:,2))
end if

! Top of domain
if (coord == nproc-1) then
    RHS_big(:,:,nz-1) = const*(u_big(:,:,nz-1)*dTdx_big(:,:,nz-1)              &
        + v_big(:,:,nz-1)*dTdy_big(:,:,nz-1)                                   &
        + 0.5_rprec*w_big(:,:,nz-1)*dTdz_big(:,:,nz-1))
end if

! Put back on the smaller grid
do k = 1, nz-1
    call dfftw_execute_dft_r2c(forw_big, RHS_big(:,:,k), RHS_big(:,:,k))
    call unpadd(RHS_T(:,:,k), RHS_big(:,:,k))
    call dfftw_execute_dft_c2r(back, RHS_T(:,:,k), RHS_T(:,:,k))
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subgrid stress
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate eddy diffusivity
Kappa_t = Nu_t/Pr_sgs

! Bottom boundary
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free: Kappa_t is stored on w-nodes
        case (0)
            pi_x(:,:,1) = -0.5_rprec*(Kappa_t(:,:,1) + Kappa_t(:,:,2))*dTdx(:,:,1)
            pi_y(:,:,1) = -0.5_rprec*(Kappa_t(:,:,1) + Kappa_t(:,:,2))*dTdy(:,:,1)

        ! Wall: Kappa_t is stored on uvp-nodes
        case (1:)
            pi_x(:,:,1) = -Kappa_t(:,:,1)*dTdx(:,:,1)
            pi_y(:,:,1) = -Kappa_t(:,:,1)*dTdy(:,:,1)

    end select
end if

! Top boundary
if (coord == nproc-1) then
    select case (ubc_mom)

      ! Stress free: Kappa_t is stored on w-nodes
      case (0)
          pi_x(:,:,nz-1) = -0.5_rprec*(Kappa_t(:,:,nz-1) + Kappa_t(:,:,nz))*dTdx(:,:,nz-1)
          pi_y(:,:,nz-1) = -0.5_rprec*(Kappa_t(:,:,nz-1) + Kappa_t(:,:,nz))*dTdy(:,:,nz-1)
          pi_z(:,:,nz) = -Kappa_t(:,:,nz-1)*dTdz(:,:,nz)

      ! Wall: Kappa_t is stored on uvp-nodes
      case (1:)
          pi_x(:,:,nz-1) = -Kappa_t(:,:,nz-1)*dTdx(:,:,nz-1)
          pi_y(:,:,nz-1) = -Kappa_t(:,:,nz-1)*dTdy(:,:,nz-1)
          pi_z(:,:,nz) = -Kappa_t(:,:,nz-1)*dTdz(:,:,nz)

      end select
end if

! Calculate rest of the domain: Kappa_t is on w nodes
do k= jz_min, jz_max
    pi_x(:,:,k) = -0.5_rprec*(Kappa_t(:,:,k) + Kappa_t(:,:,k+1))*dTdx(:,:,k)
    pi_y(:,:,k) = -0.5_rprec*(Kappa_t(:,:,k) + Kappa_t(:,:,k+1))*dTdy(:,:,k)
    pi_z(:,:,k) = -Kappa_t(:,:,k)*dTdz(:,:,k)
end do
pi_z(:,:,jz_max+1) = -Kappa_t(:,:,jz_max+1)*dTdz(:,:,jz_max+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Divergence of heat flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Store the derivatives in the stress values...can change if we need to output
! stuff
call ddx(pi_x, div_pi, lbz)
call ddy(pi_y, temp_var, lbz)
div_pi = div_pi + temp_var
call ddz_w(pi_z, temp_var, lbz)
div_pi = div_pi + temp_var

do k = 1, nz-1
    RHS_T(1:nx,:,k) = -RHS_T(1:nx,:,k) - div_pi(1:nx,:,k)
end do

! Euler integration check
if ((jt_total == 1) .and. (inits)) then
    RHS_Tf = RHS_T
end if

! Take a step
theta(1:nx,:,1:nz-1) = theta(1:nx,:,1:nz-1)                                    &
    + dt*(tadv1*RHS_T(1:nx,:,1:nz-1) + tadv2*RHS_Tf(1:nx,:,1:nz-1))

#ifdef PPMPI
call mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP)
#endif

! Use gradient at top to project temperature above domain
if (coord == nproc-1) then
    theta(:,:,nz) = theta(:,:,nz-1) + lapse_rate*dz
end if

end subroutine scalars_transport

!*******************************************************************************
subroutine to_big(a, a_big)
!*******************************************************************************
use fft
use param, only : lbz, nx, ny, nz

real(rprec), dimension(ld, ny, lbz:nz), intent(inout) ::  a
real(rprec), dimension(ld_big, ny2, lbz:nz), intent(inout) :: a_big

integer :: jz
real(rprec) :: const

! Set variables onto big domain for multiplication in physical space and
! dealiasing
const = 1._rprec/(nx*ny)
do jz = lbz, nz
    temp_var(:,:,jz) = const*a(:,:,jz)
    call dfftw_execute_dft_r2c(forw, temp_var(:,:,jz), temp_var(:,:,jz))
    call padd(a_big(:,:,jz), temp_var(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, a_big(:,:,jz), a_big(:,:,jz))
end do

end subroutine to_big

!*******************************************************************************
subroutine buoyancy_force
!*******************************************************************************
! This subroutine calculates the buoyancy term due to temperature to be added to
! the RHS of the vertical momentum equation.
use param, only : coord, nx, ny, nz
use sim_param, only :  RHSz

integer :: k, jz_min
real(rprec) :: theta_bar

! We do not advance the ground nodes, so start at k=2.
! For the MPI case, the means that we start from jz=2
! for coord=0 and jz=1 otherwise.
#ifdef PPMPI
   if (coord == 0) then
      jz_min = 2
   else
      jz_min = 1
   end if
#else
   jz_min = 2
#endif

! Add to RHSz
if ( .not.passive_scalar ) then
    do k = jz_min, nz-1
        theta_bar = sum(0.5_rprec*(theta(1:nx,:,k)+theta(1:nx,:,k-1)))/nx/ny
        RHSz(1:nx,:,k) = RHSz(1:nx,:,k) + g*(0.5_rprec*(theta(1:nx,:,k)+theta(1:nx,:,k-1)) - theta_bar)
    end do
end if

end subroutine buoyancy_force

end module scalars
