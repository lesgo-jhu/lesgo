!!
!!  Copyright (C) 2012-2013  Johns Hopkins University
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

module tcm
use types, only : rprec
use param, only : coord
use iso_c_binding
implicit none

! Values read from input
logical     :: receding_horizon
real(rprec) :: t_advance
real(rprec) :: tcm_u            ! Freestream velocity (m/s)
!real(rprec) :: tcm_k            ! Wake expansion coefficient
real(rprec), dimension(:), allocatable :: tcm_k                 !wake expansion coefficient
real(rprec) :: tcm_delta        ! Forcing width (m)
integer     :: tcm_Nx           ! Number of cells
real(rprec) :: tcm_T            ! duration of simulation (s)
real(rprec) :: tcm_cfl          ! CFL condition
real(rprec) :: tcm_alpha
real(rprec) :: tcm_gamma
real(rprec) :: tcm_eta
integer     :: tcm_maxIter
real(rprec) :: steady_state_power
real(rprec) :: tcm_t_eval
logical     :: scale_Pref

! Other derived values
real(rprec), dimension(:), allocatable :: s                 ! streamwise turbine locations

! Intermediate arrays stored in lesgo
integer :: num_t
real(rprec), dimension(:), allocatable :: vs                ! Optimization states
real(rprec), dimension(:), allocatable :: u_hat             ! estimated local velocities
real(rprec), dimension(:), allocatable :: t_carray          ! time array for tcm
real(rprec), dimension(:), allocatable :: Pref_carray       ! power reference for tcm
real(rprec), dimension(:), allocatable :: Ct_prime_carray 
real(rprec), dimension(:), allocatable :: u_meas

interface
    subroutine initialize_tcm (t, Pref, vs, u_hat, s, N, Nt, Nx, u, D, k, delta, cfl, Ctp_ref, alpha, gamma, eta, maxIter, P_ss) bind(c, name='initialize_tcm')
        import :: c_double, c_int
        real(c_double), intent(in)          :: t(*), s(*), k(*)
        real(c_double), intent(out)         :: Pref(*), vs(*), u_hat(*)
        real(c_double), intent(out)         :: P_ss
        integer(c_int), intent(in), value   :: N, Nt, Nx, maxIter
        real(c_double), intent(in), value   :: u, D, delta, cfl, Ctp_ref, alpha, gamma, eta
    end subroutine initialize_tcm
    subroutine run_tcm_wrapped (t, Ctp, Pref, vs, u_hat, s, N, Nt, Nx, u, TT, D, k, delta, cfl, Ctp_ref, alpha, gamma, eta, maxIter, t_eval, e, print_result, t_write) bind(c, name='run_model')
        import :: c_double, c_int, c_bool
        real(c_double), intent(in)          :: t(*), s(*), Pref(*), e(*), k(*)
        real(c_double), intent(out)         :: Ctp(*), vs(*), u_hat(*)
        integer(c_int), intent(in), value   :: N, Nt, Nx, maxIter
        real(c_double), intent(in), value   :: u, TT, D, delta, cfl, Ctp_ref, alpha, gamma, eta, t_eval, t_write
        logical(c_bool), intent(in),value   :: print_result
    end subroutine run_tcm_wrapped
end interface

contains

!**********************************************************************
subroutine run_tcm(first_eval_i)
!**********************************************************************
use types, only : rprec
use param, only : z_i, u_star, total_time_dim, dx, cfl
!use open_file_fid_mod
!use turbines
use turbines_base, only : Pref_list, Pref_t_list, Ct_prime_list, Ct_prime_t_list, read_values_from_file, &
    get_number_of_lines, interpolate, interpolate_vec, num_x, num_y, control, dia_all, Ct_prime_all
use stat_defs, only : wind_farm
implicit none
logical, optional, intent(in) :: first_eval_i

integer :: i, k
real(rprec), dimension(num_x) :: u_hat_prev
real(rprec), dimension(num_t) :: t_carray_dummy
logical :: first_eval
logical :: print_result

if (coord == 0) then
    print_result = .true.
else
    print_result = .false.
endif

if (present(first_eval_i)) then
    first_eval = first_eval_i
else
    first_eval = .false.
endif

u_hat_prev = u_hat

! Interpolate reference signal onto t_carray
!do k = 1, num_t
!    Pref_carray(k) = interpolate(Pref_t_list(1,:), Pref_list(1,:), t_carray(k))
!enddo
Pref_carray = interpolate_vec(Pref_t_list(1,:) - tcm_t_eval, Pref_list(1,:), t_carray)
Pref_carray = steady_state_power * Pref_carray

! This always needs to be interpolated onto t = [0,T]
if (first_eval) then
    t_carray_dummy = t_carray
else
    t_carray_dummy = t_carray - t_advance
endif
do i = 1, num_x
    Ct_prime_carray((i-1)*num_t+1:i*num_t) = interpolate_vec(t_carray_dummy,Ct_prime_carray((i-1)*num_t+1:i*num_t), t_carray)
enddo

! Call tcm
write(*,*) "u_hat = ", u_hat
write(*,*) "u_meas = ", u_meas
!u_meas(:) = 0
call run_tcm_wrapped(t_carray, Ct_prime_carray, Pref_carray, vs, u_hat, s, num_x, num_t, tcm_Nx, &
    tcm_u, tcm_T, dia_all*z_i, tcm_k, tcm_delta, tcm_cfl, Ct_prime_all, tcm_alpha, tcm_gamma, tcm_eta, &
    tcm_maxIter, t_advance, u_meas - u_hat, print_result, tcm_t_eval)

! Apply to Ct_prime_list
! This should be in the interval t = [tcm_t_eval,tcm_t_eval + T]
do i = 1, num_x
    Ct_prime_list(i,:) = Ct_prime_carray((i-1)*num_t+1:i*num_t)
    Ct_prime_t_list(i,:) = t_carray + tcm_t_eval
enddo

! Correct u_hat to not include error correction
u_hat = u_hat - (u_meas - u_hat_prev)

! Set next evaluation time
tcm_t_eval = tcm_t_eval + t_advance

end subroutine run_tcm

!**********************************************************************
subroutine tcm_init()
!**********************************************************************
use types, only : rprec
use param, only : z_i, u_star, total_time_dim, dx, cfl
!use open_file_fid_mod
!use turbines
use turbines_base, only : Pref_list, Pref_t_list, Ct_prime_list, Ct_prime_t_list, read_values_from_file, &
    get_number_of_lines, interpolate, interpolate_vec, num_x, num_y, control, dia_all, Ct_prime_all
use stat_defs, only : wind_farm
use messages
implicit none

integer :: i, j, k
real(rprec) :: dt_dummy
real(rprec), dimension(:), allocatable :: Pref_dummy

! Exit if tcm is not used
if (control < 6) return

tcm_t_eval = total_time_dim

! Get a rough estimate of the times needed for the controller.
! These values will be interpolated, so one value per several time steps is ok.
dt_dummy = cfl * dx / 18. / u_star
num_t = ceiling((tcm_T - total_time_dim)/dt_dummy) + 1
allocate(t_carray(num_t));
allocate(Pref_carray(num_t));
allocate(Pref_dummy(num_t));
t_carray = dt_dummy * [(i, i = 0, num_t-1)]

! Interpolate reference signal onto t_carray
Pref_carray = interpolate_vec(Pref_t_list(1,:) - tcm_t_eval, Pref_list(1,:), t_carray)
Pref_dummy = Pref_carray

! Find streamwise locations of turbines
allocate(s(num_x))
do k = 1, num_x
    s(k) = wind_farm%turbine(1 + (k-1)*(num_y))%xloc * z_i
enddo

! Allocate some variables
allocate(vs(num_x*(tcm_Nx+1)))
allocate(u_hat(num_x))
allocate(u_meas(num_x))
allocate(Ct_prime_t_list(num_x,num_t))
allocate(Ct_prime_list(num_x,num_t))

! Get Ct_prime list from control model
allocate(Ct_prime_carray(num_t * num_x))
Ct_prime_carray(:) = Ct_prime_all
if (scale_Pref) then
    call initialize_tcm(t_carray, Pref_carray, vs, u_hat, s, num_x, num_t, tcm_Nx, tcm_u, dia_all*z_i, &
        tcm_k, tcm_delta, tcm_cfl, Ct_prime_all, tcm_alpha, tcm_gamma, tcm_eta, tcm_maxIter, steady_state_power)
else
    call initialize_tcm(t_carray, Pref_dummy, vs, u_hat, s, num_x, num_t, tcm_Nx, tcm_u, dia_all*z_i, &
        tcm_k, tcm_delta, tcm_cfl, Ct_prime_all, tcm_alpha, tcm_gamma, tcm_eta, tcm_maxIter, steady_state_power)
    steady_state_power = 1.0
endif
!write(*,*) Pref_carray
u_meas = u_hat
call run_tcm(.true.)

end subroutine tcm_init

end module
