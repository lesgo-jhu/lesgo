
module test_mod

contains
!*******************************************************************************
subroutine P_minimize(x, f, g)
!*******************************************************************************
use types, only : rprec
use turbines, only : wm_Ct_prime_spline, wm_Cp_prime_spline
implicit none
real(rprec), dimension(:), intent(in) :: x
real(rprec), intent(inout) :: f
real(rprec), dimension(:), intent(inout) :: g
real(rprec) :: lp, b, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
real(rprec) :: a, u

b = x(1)
lp = x(2)

call wm_Ct_prime_spline%interp(b, lp, Ctp, dCtp_dbeta, dCtp_dlambda)
call wm_Cp_prime_spline%interp(b, lp, Cpp, dCpp_dbeta, dCpp_dlambda)

a = Ctp/(4._rprec + Ctp)
u = (1 - a)
f = -Cpp * u**3

g(1) = -dCpp_dbeta * u**3                                                      &
    + 3._rprec * Cpp * u**2 * 4._rprec / (4._rprec + Ctp)**2 * dCtp_dbeta
g(2) = -dCpp_dlambda * u**3                                                    &
    + 3._rprec * Cpp * u**2 * 4._rprec / (4._rprec + Ctp)**2 * dCtp_dlambda

end subroutine P_minimize

end module test_mod

!*******************************************************************************
program test
!*******************************************************************************
use types, only : rprec
use wake_model
use test_mod
use minimize
use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
use open_file_fid_mod
use functions, only : linear_interp
use cubic_spline
use turbines_mpc

use lbfgsb
use conjugate_gradient
! use minimize
! use util, only : rosenbrock
implicit none

! common variables
integer :: i, j
real(rprec) :: cfl, dt

! wake model variables
type(wake_model_t) :: wm, wmi
real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
integer :: N, Nx, Nt

type(turbines_mpc_t) :: controller
type(conjugate_gradient_t) :: m
! type(lbfgsb_t) :: m

real(rprec), dimension(:), allocatable :: Pref
real(rprec), dimension(:), allocatable :: time
! real(rprec), dimension(:), allocatable :: cx, g
! real(rprec) :: f
!
allocate(time(2))
allocate(Pref(2))
time(1) = 0._rprec
time(2) = 600._rprec
!
! type(minimize_t) mini
! type(lbfgsb_t) :: l
! real(rprec), dimension(2) :: x, ox
! real(rprec) :: lb, ub
! lb = -10000000._rprec
! ub = 0.5_rprec
!
! x = 111._rprec
! mini = minimize_t(rosenbrock)
! l = lbfgsb_t(mini)!, 1000, lb, ub)
! call l%minimize(x, ox)
! write(*,*) ox

! initialize wake model
cfl = 0.01_rprec
Dia = 126._rprec
Delta = 0.5_rprec * Dia
rho = 1.225_rprec
inertia = 4.0469e+07_rprec
torque_gain = 2.1648e6
U_infty = 9._rprec
N = 1
Nx = 256
Nt = 100*Nx
allocate(s(N))
allocate(k(N))
allocate(beta(N))
allocate(gen_torque(N))

k = 0.05_rprec
beta = 0._rprec
do i = 1, N
    s(i) = 7._rprec * Dia * i
end do
call generate_splines
wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
    wm_Ct_prime_spline, wm_Cp_prime_spline)

! integrate the wake model forward in time to get reference power
dt = cfl * wm%dx / wm%U_infty
do i = 1, Nt
    gen_torque = torque_gain * wm%omega**2
    call wm%advance(beta, gen_torque, dt)
end do
Pref = 0.9*sum(wm%Phat)
! Reset wake model and create controller
! wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
!     wm_Ct_prime_spline, wm_Cp_prime_spline)
controller = turbines_mpc_t(wm, 0._rprec, time(2), 0.99_rprec, time, Pref)
controller%beta = 0._rprec
do i = 1, N
    controller%gen_torque(i,:) = gen_torque(i)
end do
! write(*,*) gen_torque
! write(*,*) controller%gen_torques
call controller%makeDimensionless
call controller%run()
write(*,*) "dJdb:", controller%grad_beta
write(*,*) "dJdT:", controller%grad_gen_torque
!
m = conjugate_gradient_t(controller, 2)
! m = lbfgsb_t(controller, 300)
call m%minimize( controller%get_control_vector() )
call controller%run()
call controller%finite_difference_gradient()
write(*,*) "dJdb:", controller%grad_beta
write(*,*) "dJdT:", controller%grad_gen_torque
write(*,*) "fddJdb:", controller%fdgrad_beta
write(*,*) "fddJdT:", controller%fdgrad_gen_torque
! write(*,*) "Pfarm:", controller%Pfarm
! write(*,*) "Pref:", controller%Pref
! write(*,*) "beta:", controller%beta
!
! ! Run in the wake model
! call controller%makeDimensional
! do i = 1, controller%Nt
!     call wm%advance(controller%beta(:,i), controller%gen_torque(:,i), controller%dt)
!     write(*,*) wm%beta
!     write(*,*) wm%gen_torque
! end do


! write(*,*) controller%gen_torque
! allocate(g(size(cx)))
! call controller%eval(cx,f,g)
!
!
!
! real(rprec) :: Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! type(lbfgsb_t) :: lb
! type(minimize_t) :: Pm
! real(rprec), dimension(2) :: x
!
! call generate_splines()
!
! write(*,*) "beta, lambda_prime, Ctp, Cpp, dCtp_dbeta, Ctp_dlambda, dCpp_dbeta, Cpp_dlambda"
! write(*,*) "This should be near the peak power point"
! call wm_Ct_prime_spline%interp(0._rprec, 11._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(0._rprec, 11._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 0._rprec, 11._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! write(*,*) "other points"
! call wm_Ct_prime_spline%interp(-5._rprec, 5._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-5._rprec, 5._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -5._rprec, 5._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(5._rprec, 5._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(5._rprec, 5._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 5._rprec, 5._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(-5._rprec, 15._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-5._rprec, 15._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -5._rprec, 15._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(5._rprec, 15._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(5._rprec, 15._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 5._rprec, 15._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! write(*,*) "Everything should be 0:"
! call wm_Ct_prime_spline%interp(100._rprec, -1._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(100._rprec, -1._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 100._rprec, -1._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(0._rprec, -1._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(0._rprec, -1._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 0._rprec, -1._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(-100._rprec, -1._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-100._rprec, -1._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -100._rprec, -1._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(-100._rprec, 11._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-100._rprec, 11._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -100._rprec, 11._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(-100._rprec, 100._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-100._rprec, 100._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -100._rprec, 100._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! write(*,*) "Cp' = 0, Ct' is near zero. All derivatives should be zero, except maybe dCt'/dbeta"
! call wm_Ct_prime_spline%interp(-20._rprec, 100._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-20._rprec, 100._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -20._rprec, 100._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! write(*,*) "Cp' = 0, Ct' has a value, and DCt'/dlambda is positive."
! call wm_Ct_prime_spline%interp(100._rprec, 5._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(100._rprec, 5._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 100._rprec, 5._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! write(*,*) "Cp' = 0, Ct' has a value, and DCt'/dbeta is positive."
! call wm_Ct_prime_spline%interp(-5._rprec, 100._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-5._rprec, 100._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -5._rprec, 100._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! write(*,*) "Cp' = 0, Ct'=4, and all derivatives are 0:"
! call wm_Ct_prime_spline%interp(100._rprec, 20._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(100._rprec, 20._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 100._rprec, 20._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(100._rprec, 100._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(100._rprec, 100._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 100._rprec, 100._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(30._rprec, 100._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(30._rprec, 100._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 30._rprec, 100._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
!
! Pm = minimize_t(P_minimize)
! lb = lbfgsb_t(Pm)
!
! x(1) = 0._rprec
! x(2) = 17._rprec
! call lb%minimize(x,x)
! write(*,*) x

! real(rprec), dimension(:), allocatable :: x, v, xi, vi
! call wm_Ct_prime_spline%interp(-100._rprec, 100._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-100._rprec, 100._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) -100._rprec, 100._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! integer :: N, i, Ni
! real(rprec) :: dx, pi, L, dxi, x0
! type(cubic_spline_t) :: sp
! real(rprec) :: Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
!
! call generate_splines()
!
! call wm_Ct_prime_spline%interp(100._rprec, 1000._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(100._rprec, 1000._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) Ctp, dCtp_dbeta, dCtp_dlambda
! write(*,*) Cpp, dCpp_dbeta, dCpp_dlambda
! !
! N = 50
! allocate( x(N) )
! allocate( v(N) )
!
! pi = 4._rprec*atan(1._rprec)
! L = 7._rprec*pi/8._rprec
! x0 = -7._rprec*pi/16._rprec
! dx = L/(N-1)
! do i = 1, N
!     x(i) = tan(dx * (i-1) + x0)
!     v(i) = sin(x(i))
! end do
!
! write(*,*) "x = [", x, "];"
! write(*,*) "v = [", v, "];"
!
! sp = cubic_spline_t(x, v, 'not-a-knot', 'not-a-knot')!, cos(x(1)), cos(x(N)))
!
! Ni = 256
! allocate( xi(Ni) )
! allocate( vi(Ni) )
!
! L = x(N) - x(1)
! dxi = 2*L/(Ni-1)
! do i = 1, Ni
!     xi(i) = dxi * (i-1) - L/2 + x(1)
! end do
!
! call sp%interp(xi, vi)
!
! write(*,*) "xi = [", xi, "];"
! write(*,*) "vi = [", vi, "];"
!
! deallocate(x)
! deallocate(v)
! deallocate(xi)
! deallocate(vi)


end program test
