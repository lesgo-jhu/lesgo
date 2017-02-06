! rosenbrock(x, f, g)

program test

use types, only : rprec
use wake_model
use wake_model_adjoint
! use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
! use rh_control
use open_file_fid_mod
use functions, only : linear_interp
use cubic_spline
use lbfgsb_class
use minimize
use util, only : rosenbrock
implicit none

!
! ! common variables
! integer :: i, j
! real(rprec) :: cfl, dt
!
! ! wake model variables
! type(wake_model_t) :: wm
! type(wake_model_adjoint_t) :: wma
! real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
! real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
! integer :: N, Nx, Nt
!
! real(rprec) :: Pref

type(minimize_t) mini
type(lbfgsb) :: l
real(rprec), dimension(2) :: x, ox

x = 111._rprec
mini = minimize_t(rosenbrock)
l = lbfgsb(mini)
call l%minimize(x, ox)
write(*,*) ox
!
! ! initialize wake model
! cfl = 0.99_rprec
! Dia = 126._rprec
! Delta = 0.5_rprec * Dia
! rho = 1.225_rprec
! inertia = 4.0469e+07_rprec
! torque_gain = 2.1648e6
! U_infty = 9._rprec
! N = 7
! Nx = 256
! Nt = 5!*Nx
! allocate(s(N))
! allocate(k(N))
! allocate(beta(N))
! allocate(gen_torque(N))
! k = 0.05_rprec
! beta = 0._rprec
! do i = 1, N
!     s(i) = 7._rprec * Dia * i
! end do
! call generate_splines
! wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
!     wm_Ct_prime_spline, wm_Cp_prime_spline)
! wma = wake_model_adjoint_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,        &
!     wm_Ct_prime_spline, wm_Cp_prime_spline)
!
! ! integrate the wake model forward in time at least 2 flow through times
! dt = cfl * wm%dx / wm%U_infty
! do i = 1, Nt
!     do j = 1, wm%N
!         gen_torque = torque_gain * wm%omega**2 / wm%TIME**2 / wm%TORQUE
!     end do
!     call wm%advance(beta, gen_torque, dt)
!     call wm%adjoint_values(Pref, fstar(i,:,:), Adu(i,:), Aw(i,:), Bj(i,:), Bdu(i,:), Bw(i,:))
! end do
!
! ! integrate adjoints backwards in time
! do i = 1, Nt
!     call wma%retract(fstar(i,:,:), Adu(i,:), Aw(i,:), Bj(i,:), Bdu(i,:), Bw(i,:), dt)
!     write(*,*) wma%uhat_star
! end do

end program test
