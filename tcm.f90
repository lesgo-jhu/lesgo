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
use param
use iso_c_binding
implicit none

! Values read from input
real(rprec) :: tcm_u            ! Freestream velocity (m/s)
real(rprec) :: tcm_k            ! Wake expansion coefficient
real(rprec) :: tcm_delta        ! Forcing width (m)
integer     :: tcm_Nx           ! Number of cells
real(rprec) :: tcm_T            ! duration of simulation (s)
real(rprec) :: tcm_cfl          ! CFL condition
real(rprec) :: tcm_alpha
real(rprec) :: tcm_gamma
real(rprec) :: tcm_eta
integer :: tcm_maxIter

! Other derived values
real(rprec), dimension(:), allocatable :: s                 ! streamwise turbine locations

! Optimization values
real(rprec), dimension(:), allocatable :: vs                ! Optimization states
real(rprec), dimension(:), allocatable :: u_hat             ! estimated local velocities

interface
    subroutine run_tcm (t, Ctp, Pref, vs, u_hat, s, N, Nt, Nx, u, D, k, delta, cfl, Ctp_ref, alpha, gamma, eta, maxIter) bind(c, name='run_model')
        import :: c_double, c_int
        real(c_double), intent(in)          :: t(*), s(*)
        real(c_double), intent(out)         :: Ctp(*), Pref(*), vs(*), u_hat(*)
        integer(c_int), intent(in), value   :: N, Nt, Nx, maxIter
        real(c_double), intent(in), value   :: u, D, k, delta, cfl, Ctp_ref, alpha, gamma, eta
    end subroutine run_tcm
end interface

contains
!**********************************************************************
subroutine tcm_init()
!**********************************************************************
use types, only : rprec
use param, only : z_i, u_star, total_time_dim, dx, cfl
!use open_file_fid_mod
!use turbines
use turbines_base, only : Pref_list, Pref_t_list, Ct_prime_list, Ct_prime_t_list, read_values_from_file, get_number_of_lines, interpolate, num_x, num_y, control, dia_all, Ct_prime_all
use stat_defs, only : wind_farm
use messages
implicit none

integer :: i, j, k
real(rprec), dimension(:), allocatable :: t_carray, Pref_carray, Ct_prime_carray
real(rprec) :: dt_dummy
integer :: num_t

! Exit if tcm is not used
if (control /= 6) return

! Get a rough estimate of the times needed for the controller.
! These values will be interpolated, so one value per several time steps is ok.
dt_dummy = cfl * dx / 18. / u_star
num_t = ceiling((tcm_T - total_time_dim)/dt_dummy) + 1
allocate(t_carray(num_t));
allocate(Pref_carray(num_t));
t_carray = dt_dummy * [(i, i = 0, num_t-1)]

! Interpolate reference signal onto t_carray
do k = 1, num_t
    Pref_carray(k) = interpolate(Pref_t_list(1,:), Pref_list(1,:), t_carray(k))
end do

! Find streamwise locations of turbines
allocate(s(num_x))
do k = 1, num_x
    s(k) = wind_farm%turbine(1 + (k-1)*(num_y))%xloc * z_i
enddo

! Allocate some variables
allocate(vs(num_x*tcm_Nx))
allocate(u_hat(num_x))

! Get Ct_prime list from control model
allocate(Ct_prime_carray(num_t * num_x))
call run_tcm(t_carray, Ct_prime_carray, Pref_carray, vs, u_hat, s, num_x, num_t, tcm_Nx, tcm_u, dia_all*z_i, tcm_k, tcm_delta, tcm_cfl, Ct_prime_all, tcm_alpha, tcm_gamma, tcm_eta, tcm_maxIter)

! Apply to Ct_prime_list
allocate(Ct_prime_t_list(num_x,num_t))
allocate(Ct_prime_list(num_x,num_t))
do i = 1, num_x
    Ct_prime_list(i,:) = Ct_prime_carray((i-1)*num_t+1:i*num_t)
    Ct_prime_t_list(i,:) = t_carray        
enddo

end subroutine tcm_init

end module
