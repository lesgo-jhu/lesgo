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
integer, public :: sgs_model_scal
real(rprec), public :: Pr_sgs = 0.5
integer, public :: cs_count_scal, dyn_init_scal
real(rprec), public, dimension(:,:,:), allocatable :: Kappa_t
real(rprec), dimension(:), allocatable :: sigma_theta
real(rprec), dimension(:,:,:), allocatable :: I_LM_t, I_MM_t, I_QN_t, I_NN_t,  &
    s_Beta_t, s_Tn_all_t, Ds_opt2, ds2_clips
real(rprec), dimension(:,:), allocatable :: theta_bar, theta_hat, L1, L2, L3,  &
    M1, M2, M3, Q1, Q2, Q3, N1, N2, N3, dTdx_bar, dTdy_bar, dTdz_bar, dTdx_hat,&
    dTdy_hat, dTdz_hat, S_dTdx_bar, S_dTdy_bar, S_dTdz_bar, S_dTdx_hat,        &
    S_dTdy_hat, S_dTdz_hat, LM, MM, QN, NN, s_Tn, dumfac, epsi, Ds_opt2_2d,    &
    Ds_opt2_4d
logical :: I_LM_t_MM_init = .false.
logical :: I_QN_t_NN_init = .false.
logical :: inilag_scalar
real(rprec) :: lagran_dt_scalar = 0._rprec

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
allocate ( sigma_theta(lbz:nz) ); sigma_theta = 0._rprec

! Obukhov values (defaults for passive scalars)
allocate ( psi_m(nx, ny) ); psi_m = 0._rprec
allocate ( phi_m(nx, ny) ); phi_m = 1._rprec
allocate ( psi_h(nx, ny) ); psi_h = 0._rprec
allocate ( phi_h(nx, ny) ); phi_h = 1._rprec
allocate ( L(nx, ny) ); L = 0._rprec
allocate ( tstar_lbc(nx, ny) ); tstar_lbc = 0._rprec

!For Lagrangian dynamic scale dependent model - temperature
allocate ( I_LM_t(ld,ny,lbz:nz) ); I_LM_t = 0._rprec
allocate ( I_MM_t(ld,ny,lbz:nz) ); I_MM_t = 0._rprec
allocate ( I_QN_t(ld,ny,lbz:nz) ); I_QN_t = 0._rprec
allocate ( I_NN_t(ld,ny,lbz:nz) ); I_NN_t = 0._rprec
allocate ( s_Beta_t(ld,ny,lbz:nz) ); s_Beta_t = 0._rprec
allocate ( s_Tn_all_t(ld,ny,lbz:nz) ); s_Tn_all_t = 0._rprec
allocate ( Ds_opt2(ld,ny,lbz:nz) ); Ds_opt2 = 0._rprec
allocate ( ds2_clips(ld,ny,lbz:nz) ); ds2_clips = 0._rprec
allocate ( theta_bar(ld,ny) ); theta_bar = 0._rprec
allocate ( theta_hat(ld,ny) ); theta_hat = 0._rprec
allocate ( L1(ld,ny) ); L1 = 0._rprec
allocate ( L2(ld,ny) ); L2 = 0._rprec
allocate ( L3(ld,ny) ); L3 = 0._rprec
allocate ( M1(ld,ny) ); M1 = 0._rprec
allocate ( M2(ld,ny) ); M2 = 0._rprec
allocate ( M3(ld,ny) ); M3 = 0._rprec
allocate ( Q1(ld,ny) ); Q1 = 0._rprec
allocate ( Q2(ld,ny) ); Q2 = 0._rprec
allocate ( Q3(ld,ny) ); Q3 = 0._rprec
allocate ( N1(ld,ny) ); N1 = 0._rprec
allocate ( N2(ld,ny) ); N2 = 0._rprec
allocate ( N3(ld,ny) ); N3 = 0._rprec
allocate ( dTdx_bar(ld,ny) ); dTdx_bar = 0._rprec
allocate ( dTdy_bar(ld,ny) ); dTdy_bar = 0._rprec
allocate ( dTdz_bar(ld,ny) ); dTdz_bar = 0._rprec
allocate ( dTdx_hat(ld,ny) ); dTdx_hat = 0._rprec
allocate ( dTdy_hat(ld,ny) ); dTdy_hat = 0._rprec
allocate ( dTdz_hat(ld,ny) ); dTdz_hat = 0._rprec
allocate ( S_dTdx_bar(ld,ny) ); S_dTdx_bar = 0._rprec
allocate ( S_dTdy_bar(ld,ny) ); S_dTdy_bar = 0._rprec
allocate ( S_dTdz_bar(ld,ny) ); S_dTdz_bar = 0._rprec
allocate ( S_dTdx_hat(ld,ny) ); S_dTdx_hat = 0._rprec
allocate ( S_dTdy_hat(ld,ny) ); S_dTdy_hat = 0._rprec
allocate ( S_dTdz_hat(ld,ny) ); S_dTdz_hat = 0._rprec
allocate ( LM(ld,ny) ); LM = 0._rprec
allocate ( MM(ld,ny) ); MM = 0._rprec
allocate ( QN(ld,ny) ); QN = 0._rprec
allocate ( NN(ld,ny) ); NN = 0._rprec
allocate ( s_Tn(ld,ny) ); s_Tn = 0._rprec
allocate ( dumfac(ld,ny) ); dumfac = 0._rprec
allocate ( epsi(ld,ny) ); epsi = 0._rprec
allocate ( Ds_opt2_2d(ld,ny) ); Ds_opt2_2d = 0._rprec
allocate ( Ds_opt2_4d(ld,ny) ); Ds_opt2_4d = 0._rprec

! Nondimensionalize variables
g = g*(z_i/(u_star**2))
flux_bot = flux_bot/u_star/T_scale
scal_bot = scal_bot/T_scale
lapse_rate = lapse_rate/T_scale*z_i
ic_theta = ic_theta/T_scale
ic_z = ic_z/z_i

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
subroutine ic_scal
!*******************************************************************************
! Set initial profile for scalar
use param, only : coord
use string_util
use grid_m

fname = path // 'scal.out'
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif
inquire (file=fname, exist=inits)
inits = .not.inits

if (inits) then
    write(*,*)  "--> Creating initial boundary layer scalar field with LES BCs"
    call ic_scal_les
    inilag_scalar = .true.
else
    write(*,*)  "--> Reading initial scalar field from file"
    call ic_scal_file
    inilag_scalar = .false.
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
read(12) theta(:, :, 1:nz), RHS_T(:, :, 1:nz), psi_m(1:nx, :),                 &
    Ds_opt2(:,:,1:nz), I_LM_t(:,:,1:nz), I_MM_t(:,:,1:nz),                     &
    I_QN_t(:,:,1:nz), I_NN_t(:,:,1:nz)
close(12)

#ifdef PPMPI
call mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(RHS_T, 0, MPI_SYNC_DOWNUP)
#endif

end subroutine ic_scal_file

!*******************************************************************************
subroutine ic_scal_les
!*******************************************************************************
use param, only : nz, lbz
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
subroutine scalars_checkpoint
!*******************************************************************************
use param, only : nx, nz, write_endian

!  Open scal.out (lun_default in io) for final output
open(11, file=fname, form='unformatted', convert=write_endian,                 &
    status='unknown', position='rewind')
write (11) theta(:, :, 1:nz), RHS_T(:, :, 1:nz), psi_m(1:nx, :),               &
    Ds_opt2(:,:,1:nz), I_LM_t(:,:,1:nz), I_MM_t(:,:,1:nz),                     &
    I_QN_t(:,:,1:nz), I_NN_t(:,:,1:nz)
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
    jt_total, dt, use_cfl_dt, jt, initu
use param, only : lbc_mom, ubc_mom, dz
use sim_param, only : u, v, w
use sgs_param, only : Nu_t, delta, S,S11, S12, S13, S22, S23, S33
use derivatives, only : filt_da, ddx, ddy, ddz_uv, ddz_w
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP
use test_filtermodule
use fft
use messages, only : error

integer :: k, jz_min, jz_max, jx, jy, jz
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
if (sgs_model_scal == 1) then
    Kappa_t = Nu_t*Pr_sgs

else if (sgs_model_scal == 5) then
    if (use_cfl_dt) then
        if ( ( jt .GE. dyn_init_scal-cs_count_scal + 1 ) .OR.  initu ) then
            lagran_dt_scalar = lagran_dt_scalar + dt
        endif
    else
    lagran_dt_scalar = cs_count_scal*dt
    end if

    ! Update Sij, Cs every cs_count timesteps (specified in param)
    elseif ( ((jt >= dyn_init_scal).OR.(initu)) .and.                        &
        (mod(jt_total, cs_count_scal)==0) ) then

        if (jt == dyn_init_scal) then
            write(*,*) "running Lagrangian dynamic sgs_model for teperature"
        end if

    call scalars_lagrange_Sdep()

    !Calculate eddy diffusivity kappa
    !Stored on w-nodes for entire domain except
    !on uvp node for jz=1 and 'wall' BC
    do jz=1,nz
       do jy=1,ny
          do jx=1,nx
             S(jx,jy) = sqrt( 2._rprec*(S11(jx,jy,jz)**2 +           S22(jx,jy,jz)**2 +&
                                        S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +&
                                        S13(jx,jy,jz)**2 +           S23(jx,jy,jz)**2 )))
             kappa_t(jx,jy,jz)=S(jx,jy)*Ds_opt2(jx,jy,jz)*delta**2
          end do
       end do
    end do

else
    call error('scalar_transport', 'invalid sgs_model_scal')
end if

! Bottom boundary
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free: Kappa_t is stored on w-nodes
        case (0)
            pi_x(:,:,1) = -(Kappa_t(:,:,1) + Kappa_t(:,:,2))*dTdx(:,:,1)
            pi_y(:,:,1) = -(Kappa_t(:,:,1) + Kappa_t(:,:,2))*dTdy(:,:,1)

        ! Wall: Kappa_t is stored on uvp-nodes
        case (1:)
            pi_x(:,:,1) = -2*Kappa_t(:,:,1)*dTdx(:,:,1)
            pi_y(:,:,1) = -2*Kappa_t(:,:,1)*dTdy(:,:,1)

    end select
end if

! Top boundary
if (coord == nproc-1) then
    select case (ubc_mom)

      ! Stress free: Kappa_t is stored on w-nodes
      case (0)
          pi_x(:,:,nz-1) = -(Kappa_t(:,:,nz-1) + Kappa_t(:,:,nz))*dTdx(:,:,nz-1)
          pi_y(:,:,nz-1) = -(Kappa_t(:,:,nz-1) + Kappa_t(:,:,nz))*dTdy(:,:,nz-1)
          pi_z(:,:,nz) = -2*Kappa_t(:,:,nz-1)*dTdz(:,:,nz)

      ! Wall: Kappa_t is stored on uvp-nodes
      case (1:)
          pi_x(:,:,nz-1) = -2*Kappa_t(:,:,nz-1)*dTdx(:,:,nz-1)
          pi_y(:,:,nz-1) = -2*Kappa_t(:,:,nz-1)*dTdy(:,:,nz-1)
          pi_z(:,:,nz) = -2*Kappa_t(:,:,nz-1)*dTdz(:,:,nz)

      end select
end if

! Calculate rest of the domain: Kappa_t is on w nodes
do k= jz_min, jz_max
    pi_x(:,:,k) = -(Kappa_t(:,:,k) + Kappa_t(:,:,k+1))*dTdx(:,:,k)
    pi_y(:,:,k) = -(Kappa_t(:,:,k) + Kappa_t(:,:,k+1))*dTdy(:,:,k)
    pi_z(:,:,k) = -2*Kappa_t(:,:,k)*dTdz(:,:,k)
end do
pi_z(:,:,jz_max+1) = -2*Kappa_t(:,:,jz_max+1)*dTdz(:,:,jz_max+1)

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

!*******************************************************************************
subroutine scalars_lagrange_Sdep()
!*******************************************************************************
! Lagrangian scale-dependent dynamic model to calculate the
! equivalent of the Smagorinsky coefficient for the eddy diffusivity
! for scalarsthis is done layer-by-layer to save memory.
! note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary.  mem requirement ~ n^2, not n^3
use param, only : ld, nx, ny, nz, coord, jt, nproc, use_cfl_dt, dt
use sgs_param, only : opftime, delta, S11, S12, S13, S22, S23, S33, delta, S,  &
    u_bar, v_bar, w_bar, u_hat, v_hat, w_hat, S_hat, S11_hat, S12_hat, S13_hat,&
    S22_hat, S23_hat, S33_hat, S_bar, S11_bar, S12_bar, S13_bar, S22_bar,      &
    S23_bar, S33_bar
use sim_param, only : u, v, w
use test_filtermodule
#ifdef PPMPI
use mpi_defs, only:mpi_sync_real_array,MPI_SYNC_DOWNUP
#endif

integer :: jx, jy, jz, i, j

real(rprec):: tf1,tf2,tf1_2,tf2_2 ! Size of the second test filter
real(rprec) :: fractus
real(rprec) :: Betaclip  !--scalar to save mem., otherwise (ld,ny,nz)
real(rprec) :: const
real(rprec) :: opftdelta,powcoeff

real(rprec), parameter :: zero = 1.e-24_rprec ! zero = infimum(0)
integer :: count_ds2_all, count_ds2_clip

! Set coefficients
opftdelta = opftime*delta
powcoeff = -1._rprec/4._rprec
fractus= 1._rprec/real(ny*nx,kind=rprec)
const = (delta**2) !Note difference from const in eddy viscosity model
tf1=2._rprec
tf2=4._rprec
tf1_2=tf1**2
tf2_2=tf2**2

!Get standard deviation of scalar quantity for the calculation of Tn
call get_scalar_std(theta, sigma_theta)

!if (coord==6) write(*,*) 'A. jt, coord, sigma_theta(nz-1)', jt, coord, sigma_theta(nz-1)
!if (coord==6) write(*,*) 'A. jt, coord, sigma_theta(nz)', jt, coord, sigma_theta(nz)
!if (coord==7) write(*,*) 'A. jt, coord, sigma_theta(0)', jt, coord, sigma_theta(0)
!if (coord==7) write(*,*) 'A. jt, coord, sigma_theta(1)', jt, coord, sigma_theta(1)

! "Rearrange" I_* (running averages) so that their new positions (i,j,k)
!  correspond to the current (i,j,k) particle
call scalars_interpolag_Sdep()
!
! For each horizontal level, calculate L_i(:,:), Q_i(:,:), M_i(:,:), and N_i(:,:).
! Then update the running averages, I_*(:,:,jz), which are used to
! calculate Ds_opt2(:,:,jz).
do jz = 1, nz
   count_ds2_all = 0
   count_ds2_clip = 0

   ! Calculate L_i
   !Interp u,v,w,scalar onto w-nodes and store result as u_bar,v_bar,w_bar
   !(except for very first level which should be on uvp-nodes
   if ( (coord == 0) .and. (jz == 1) ) then
      u_bar(:,:) = u(:,:,1)
      v_bar(:,:) = v(:,:,1)
      w_bar(:,:) = 0.25_rprec*w(:,:,2)
      theta_bar(:,:) = theta(:,:,1)
   else !w-nodes
      u_bar(:,:) = 0.5_rprec*( u(:,:,jz) + u(:,:,jz-1) )
      v_bar(:,:) = 0.5_rprec*( v(:,:,jz) + v(:,:,jz-1) )
      w_bar(:,:) = w(:,:,jz)
      theta_bar(:,:) = 0.5_rprec*( theta(:,:,jz) + theta(:,:,jz-1) )
   end if
   u_hat = u_bar
   v_hat = v_bar
   w_hat = w_bar
   theta_hat = theta_bar

   !First term before filtering
   L1 = u_bar*theta_bar
   L2 = v_bar*theta_bar
   L3 = w_bar*theta_bar

   Q1 = L1
   Q2 = L2
   Q3 = L3

   !Filter first term and add the second term to get the final value
   call test_filter ( u_bar )   ! in-place filtering
   call test_filter ( v_bar )
   call test_filter ( w_bar )
   call test_filter ( theta_bar )
   call test_filter ( L1 )
   L1 = L1 - u_bar*theta_bar
   call test_filter ( L2 )
   L2 = L2 - v_bar*theta_bar
   call test_filter ( L3 )
   L3 = L3 - w_bar*theta_bar

   call test_test_filter ( u_hat )
   call test_test_filter ( v_hat )
   call test_test_filter ( w_hat )
   call test_test_filter ( theta_hat )
   call test_test_filter ( Q1 )
   Q1 = Q1 - u_hat*theta_hat
   call test_test_filter ( Q2 )
   Q2 = Q2 - v_hat*theta_hat
   call test_test_filter ( Q3 )
   Q3 = Q3 - w_hat*theta_hat
   !
   !calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +&
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

   !already on w-nodes
   S11_bar(:,:) = S11(:,:,jz)
   S12_bar(:,:) = S12(:,:,jz)
   S13_bar(:,:) = S13(:,:,jz)
   S22_bar(:,:) = S22(:,:,jz)
   S23_bar(:,:) = S23(:,:,jz)
   S33_bar(:,:) = S33(:,:,jz)

   S11_hat = S11_bar
   S12_hat = S12_bar
   S13_hat = S13_bar
   S22_hat = S22_bar
   S23_hat = S23_bar
   S33_hat = S33_bar

   call test_filter ( S11_bar )
   call test_filter ( S12_bar )
   call test_filter ( S13_bar )
   call test_filter ( S22_bar )
   call test_filter ( S23_bar )
   call test_filter ( S33_bar )

   call test_test_filter ( S11_hat )
   call test_test_filter ( S12_hat )
   call test_test_filter ( S13_hat )
   call test_test_filter ( S22_hat )
   call test_test_filter ( S23_hat )
   call test_test_filter ( S33_hat )

   !Calculate |S_bar| (the test_filtered Sij)
   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   !Calculate |S_hat| (the test_test_filtered Sij)
   S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
        2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

   !Filter the temperature gradients
   if ( (coord==0).and.(jz==1) ) then !store on uvp-nodes
      dTdx_bar(:,:) = dTdx(:,:,1)
      dTdy_bar(:,:) = dTdy(:,:,1)
      dTdz_bar(:,:) = dTdz(:,:,2) !Try this
   else
      dTdx_bar(:,:) = 0.5_rprec*( dTdx(:,:,jz) + dTdx(:,:,jz-1) )
      dTdy_bar(:,:) = 0.5_rprec*( dTdy(:,:,jz) + dTdy(:,:,jz-1) )
      dTdz_bar(:,:) = dTdz(:,:,jz)
   end if

   dTdx_hat = dTdx_bar
   dTdy_hat = dTdy_bar
   dTdz_hat = dTdz_bar

   !Calculate |S|ds_i
   S_dTdx_bar = S * dTdx_bar
   S_dTdy_bar = S * dTdy_bar
   S_dTdz_bar = S * dTdz_bar

   S_dTdx_hat = S_dTdx_bar
   S_dTdy_hat = S_dTdy_bar
   S_dTdz_hat = S_dTdz_bar

   !Test filter temperature gradients
   call test_filter ( dTdx_bar )
   call test_filter ( dTdy_bar )
   call test_filter ( dTdz_bar )

   call test_test_filter ( dTdx_hat )
   call test_test_filter ( dTdy_hat )
   call test_test_filter ( dTdz_hat )

   call test_filter ( S_dTdx_bar )
   call test_filter ( S_dTdy_bar )
   call test_filter ( S_dTdz_bar )

   call test_test_filter ( S_dTdx_hat )
   call test_test_filter ( S_dTdy_hat )
   call test_test_filter ( S_dTdz_hat )

   !Calculate M_i and N_i
   M1 = const*(S_dTdx_bar - tf1_2*S_bar*dTdx_bar)
   M2 = const*(S_dTdy_bar - tf1_2*S_bar*dTdy_bar)
   M3 = const*(S_dTdz_bar - tf1_2*S_bar*dTdz_bar)

   N1 = const*(S_dTdx_hat - tf2_2*S_hat*dTdx_hat)
   N2 = const*(S_dTdy_hat - tf2_2*S_hat*dTdy_hat)
   N3 = const*(S_dTdz_hat - tf2_2*S_hat*dTdz_hat)
   !
   !Calculate L_i*M_i, M_i*M_i, Q_i*N_i, and N_i*N_i for each
   !point in the plane
   LM = L1*M1 + L2*M2 + L3*M3
   MM = M1*M1 + M2*M2 + M3*M3
   QN = Q1*N1 + Q2*N2 + Q3*N3
   NN = N1*N1 + N2*N2 + N3*N3

   !Initialize sgs quantities
   if (inilag_scalar) then
      if ((.not.I_LM_t_MM_init).and.(jt==cs_count_scal.or.jt==dyn_init_scal)) then
         !print*,'I_LM_t and I_MM_t initialized'
         I_MM_t (:,:,jz) = MM
         I_LM_t (:,:,jz) = 0.03_rprec * MM
         I_MM_t (ld-1:ld,:,jz) = 1._rprec
         I_LM_t (ld-1:ld,:,jz) = 1._rprec

         if (jz==nz) I_LM_t_MM_init=.true.
      end if
   end if

   !Update running averages (I_LM_t, I_MM_t)
   !Determine averaging timescale (for 2-delta filter)
   s_Tn = max( I_LM_t(:,:,jz)*I_MM_t(:,:,jz),zero )
   !s_Tn = opftdelta*(s_Tn**powcoeff)
   s_Tn = opftdelta*sigma_theta(jz)*(s_Tn**powcoeff)
   !Clip if necessary
   s_Tn(:,:) = max( zero, s_Tn(:,:) )

   !Calculate new running average
   dumfac = lagran_dt_scalar/s_Tn
   epsi = dumfac / (1._rprec+dumfac)

   I_LM_t(:,:,jz) = epsi*LM + (1._rprec-epsi)*I_LM_t(:,:,jz)
   I_MM_t(:,:,jz) = epsi*MM + (1._rprec-epsi)*I_MM_t(:,:,jz)
   !Clip if necessary
   I_LM_t(:,:,jz) = max( zero, I_LM_t(:,:,jz) )

   !Calculate Ds_opt2 (for 2-delta filter)
   ! Add +zero in denominator to avoid division by identically zero
   Ds_opt2_2d(:,:) = I_LM_t(:,:,jz)/(I_MM_t(:,:,jz)+zero)
   Ds_opt2_2d(ld,:) = zero
   Ds_opt2_2d(ld-1,:) = zero
   !Clip if necessary
   Ds_opt2_2d(:,:) = max( zero, Ds_opt2_2d(:,:) )

   !Initialize sgs quantities
   if (inilag_scalar) then
      if ((.not.I_QN_t_NN_init).and.(jt==cs_count_scal.or.jt==dyn_init_scal)) then
         !print*, 'I_NN_t and I_QN_t initialized'
         I_NN_t(:,:,jz) = NN
         I_QN_t(:,:,jz) = 0.03_rprec*NN
         I_NN_t(ld-1:ld,:,jz) = 1._rprec
         I_QN_t(ld-1:ld,:,jz) = 1._rprec

         if (jz==nz) I_QN_t_NN_init=.true.
      end if
   end if

   !Update running averages
   s_Tn = max( I_QN_t(:,:,jz)*I_NN_t(:,:,jz), zero )
   !s_Tn = opftdelta*(s_Tn**powcoeff)
   s_Tn = opftdelta*sigma_theta(jz)*(s_Tn**powcoeff)
   !Clip, if necessary
   s_Tn(:,:) = max( zero, s_Tn(:,:) )

   !Calculate new running average
   dumfac = lagran_dt_scalar/s_Tn
   epsi = dumfac / (1._rprec+dumfac)

   I_QN_t(:,:,jz) = epsi*QN + (1._rprec-epsi)*I_QN_t(:,:,jz)
   I_NN_t(:,:,jz) = epsi*NN + (1._rprec-epsi)*I_NN_t(:,:,jz)
   !Clip if necessary
   I_QN_t(:,:,jz) = max( zero, I_QN_t(:,:,jz) )

   !Calculate Ds_opt2 (for 4-delta filter)
   !Add +zero in denominator to avoid division by identically zero
   Ds_opt2_4d(:,:) = I_QN_t(:,:,jz)/(I_NN_t(:,:,jz) + zero)
   Ds_opt2_4d(ld,:) = zero
   Ds_opt2_4d(ld-1,:) = zero
   !Clip if necessary
   Ds_opt2_4d(:,:) = max( zero, Ds_opt2_4d(:,:) )

   s_Beta_t(:,:,jz) = &
       (Ds_opt2_4d(:,:)/Ds_opt2_2d(:,:))**(log(tf1)/(log(tf2)-log(tf1)))

   !--MPI
#ifdef PPMPI
   if ((coord==nproc-1).and.(jz==nz)) then
       s_Beta_t(:,:,jz) = 1._rprec
   end if
#else
   if (jz==nz) then
       s_Beta_t(:,:,jz) = 1._rprec
   end if
#endif

   !Clip s_Beta_t and
   do jy=1,ny
      do jx=1,ld
         Betaclip = max(s_Beta_t(jx,jy,jz),1._rprec/8._rprec)
         Ds_opt2(jx,jy,jz) = Ds_opt2_2d(jx,jy)/Betaclip
      end do
   end do
   Ds_opt2(ld,:,jz) = zero
   Ds_opt2(ld-1,:,jz) = zero


   !Count how often Ds is clipped
   !if (coord == 0) then
   do i=1,nx
      do j=1,ny
         if (real(Ds_opt2(i,j,jz)).lt.real(zero, kind=rprec)) then
            count_ds2_clip = count_ds2_clip + 1
            ds2_clips(i,j,jz) = 1._rprec
         end if
         count_ds2_all = count_ds2_all + 1
       end do
       !print*, 'jt, coord, count_ds2_clip, sum(ds2_clips(:,:,jz))', jt, coord, count_ds2_clip, sum(ds2_clips(:,:,jz))
   end do
   ! end if

   !Clip if necessary
   Ds_opt2(:,:,jz) = max( zero, Ds_opt2(:,:,jz) )

   !Save s_Tn to 3D array for use with tavg_scalar_sgs
   s_Tn_all_t(:,:,jz) = s_Tn(:,:)

end do

!Share new data between overlapping nodes
#ifdef PPMPI
   call mpi_sync_real_array( I_LM_t, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( I_MM_t, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( I_QN_t, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( I_NN_t, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( s_Tn_all_t, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( ds2_clips, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( s_Beta_t, 0, MPI_SYNC_DOWNUP )
#endif

!Reset variable for use during next set of cs_count_scal timesteps
if (use_cfl_dt) lagran_dt_scalar = 0._rprec

end subroutine scalars_lagrange_Sdep

!*******************************************************************************
subroutine get_scalar_std(theta,sigma_theta)
!*******************************************************************************
! use types, only: rprec
use param, only: ld, nx, ny, nz, lbz, coord, jt
use functions, only: interp_to_w_grid
#ifdef PPMPI
use mpi_defs, only: mpi_sync_real_array, MPI_SYNC_DOWNUP
#endif
!
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: theta
real(rprec), dimension(lbz:nz), intent(out) :: sigma_theta
real(rprec), dimension(lbz:nz) :: theta_avg_zplane, temp
real(rprec), dimension(ld,ny,lbz:nz) :: theta_w
real(rprec):: favg
integer :: jx, jy, jz

theta_w = interp_to_w_grid(theta,lbz)

favg = real(nx*ny,kind=rprec)


do jz=lbz,nz
   theta_avg_zplane(jz) = 0._rprec
   do jy=1,ny
      do jx=1,nx
         theta_avg_zplane(jz) = theta_avg_zplane(jz) + theta_w(jx,jy,jz)
      end do
   end do
   theta_avg_zplane(jz) = theta_avg_zplane(jz)/favg
end do

do jz=lbz,nz
   temp(jz) = 0.0_rprec
   do jy=1,ny
      do jx=1,nx
         temp(jz) = temp(jz) + (theta_w(jx,jy,jz) - theta_avg_zplane(jz))**2._rprec
      end do
   end do
   sigma_theta(jz) = sqrt(temp(jz)/favg)
end do

end subroutine get_scalar_std

!*******************************************************************************
subroutine scalars_interpolag_Sdep()
!*******************************************************************************
! This subroutine takes the arrays I_{LM,MM,QN,NN} from the previous
!   timestep and essentially moves the values around to follow the
!   corresponding particles. The (x, y, z) value at the current
!   timestep will be the (x-u*dt, y-v*dt, z-w*dt) value at the
!   previous timestep.  Since particle motion does not conform to
!   the grid, an interpolation will be required.  Variables should
!   be on the w-grid.

! This subroutine assumes that dt and cs_count_scal are chosen such that
!   the Lagrangian CFL in the z-direction will never exceed 1.  If the
!   Lag. CFL in the x-direction is less than one this should generally
!   be satisfied.

use param
use sim_param, only : u, v, w
use grid_m, only : grid
use functions, only : trilinear_interp
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#endif
use cfl_util, only : get_max_cfl

real(rprec), dimension(3) :: xyz_past

real(rprec), dimension(ld,ny,lbz:nz) :: tempI_LM_t, tempI_MM_t, tempI_QN_t, tempI_NN_t
integer :: i, j, k, kmin

real (rprec) :: lcfl

real(rprec), pointer, dimension(:) :: x,y,z


nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

! Perform (backwards) Lagrangian interpolation
! I_* arrays should be synced at this point (for MPI)

! Create dummy arrays so information will not be overwritten during interpolation
tempI_LM_t = I_LM_t
tempI_MM_t = I_MM_t
tempI_QN_t = I_QN_t
tempI_NN_t = I_NN_t

! Loop over domain (within proc): for each, calc xyz_past then trilinear_interp
! Variables x,y,z, F_LM, F_MM, F_QN, F_NN, etc are on w-grid
! Interpolation out of top/bottom of domain is not permitted.
! Note: x,y,z values are only good for k=1:nz-1 within each proc
if ( coord.eq.0 ) then
    kmin = 2
    ! At the bottom-most level (at the wall) the velocities are zero.
    ! Since there is no movement the values of F_LM, F_MM, etc should
    !   not change and no interpolation is necessary.
else
    kmin = 1
endif

! Intermediate levels
do k=kmin,nz-1
do j=1,ny
do i=1,nx
    ! Determine position at previous timestep (u,v interp to w-grid)
    xyz_past(1) = x(i) - 0.5_rprec*(u(i,j,k-1)+u(i,j,k))*lagran_dt_scalar
    xyz_past(2) = y(j) - 0.5_rprec*(v(i,j,k-1)+v(i,j,k))*lagran_dt_scalar
    xyz_past(3) = z(k) - w(i,j,k)*lagran_dt_scalar

    ! Interpolate
    I_LM_t(i,j,k) = trilinear_interp(tempI_LM_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
    I_MM_t(i,j,k) = trilinear_interp(tempI_MM_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
    I_QN_t(i,j,k) = trilinear_interp(tempI_QN_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
    I_NN_t(i,j,k) = trilinear_interp(tempI_NN_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
enddo
enddo
enddo

! Top-most level should not allow negative w
#ifdef PPMPI
if (coord.eq.nproc-1) then
#endif
    k = nz
    do j=1,ny
    do i=1,nx
        ! Determine position at previous timestep (u,v interp to w-grid)
        xyz_past(1) = x(i) - 0.5_rprec*(u(i,j,k-1)+u(i,j,k))*lagran_dt_scalar
        xyz_past(2) = y(j) - 0.5_rprec*(v(i,j,k-1)+v(i,j,k))*lagran_dt_scalar
        xyz_past(3) = z(k) - max(0.0_rprec,w(i,j,k))*lagran_dt_scalar

        ! Interpolate
        I_LM_t(i,j,k) = trilinear_interp(tempI_LM_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
        I_MM_t(i,j,k) = trilinear_interp(tempI_MM_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
        I_QN_t(i,j,k) = trilinear_interp(tempI_QN_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
        I_NN_t(i,j,k) = trilinear_interp(tempI_NN_t(1:nx,1:ny,lbz:nz),lbz,xyz_past)
    enddo
    enddo
#ifdef PPMPI
endif
#endif

 ! Share new data between overlapping nodes
#ifdef PPMPI
call mpi_sync_real_array( I_LM_t, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( I_MM_t, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( I_QN_t, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( I_NN_t, 0, MPI_SYNC_DOWNUP )
#endif

! Compute the Lagrangian CFL number and print to screen
!   Note: this is only in the x-direction... not good for complex geometry cases
if (mod (jt_total, lag_cfl_count) .eq. 0) then
    lcfl = get_max_cfl()
    lcfl = lcfl*lagran_dt_scalar/dt
#ifdef PPMPI
    if(coord.eq.0) print*, 'Lagrangian CFL condition= ', lcfl
#else
    print*, 'Lagrangian CFL condition= ', lcfl
#endif
endif

nullify(x,y,z)

end subroutine scalars_interpolag_Sdep

end module scalars
