program test

use types, only : rprec
use wake_model
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
! 
! ! common variables
! integer :: i, j
! real(rprec) :: cfl, dt
! 
! ! wake model variables
! type(wake_model_t) :: wm, wmi
! real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
! real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
! integer :: N, Nx, Nt
! 
! type(turbines_mpc_t) :: controller
! ! type(conjugate_gradient_t) :: m
! type(lbfgsb_t) :: m
! 
! real(rprec), dimension(:), allocatable :: Pref
! real(rprec), dimension(:), allocatable :: time
! real(rprec), dimension(:), allocatable :: cx, g
! real(rprec) :: f
! 
! allocate(time(2))
! allocate(Pref(2))
! time(1) = 0._rprec
! time(2) = 600._rprec
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
! 
! ! initialize wake model
! cfl = 0.01_rprec
! Dia = 126._rprec
! Delta = 0.5_rprec * Dia
! rho = 1.225_rprec
! inertia = 4.0469e+07_rprec
! torque_gain = 2.1648e6
! U_infty = 9._rprec
! N = 7
! Nx = 256
! Nt = 100*Nx
! allocate(s(N))
! allocate(k(N))
! allocate(beta(N))
! allocate(gen_torque(N))
! 
! k = 0.05_rprec
! beta = 0._rprec
! beta(N) = 5._rprec
! do i = 1, N
!     s(i) = 7._rprec * Dia * i
! end do
! call generate_splines
! wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
!     wm_Ct_prime_spline, wm_Cp_prime_spline)
! 
! ! integrate the wake model forward in time to get reference power
! dt = cfl * wm%dx / wm%U_infty
! do i = 1, Nt
!     gen_torque = torque_gain * wm%omega**2
!     call wm%advance(beta, gen_torque, dt)
! end do
! Pref = 1.05*sum(wm%Phat)
! 
! ! Reset wake model and create controller
! ! wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
! !     wm_Ct_prime_spline, wm_Cp_prime_spline)
! controller = turbines_mpc_t(wm, 0._rprec, time(2), 0.99_rprec, time, Pref)
! controller%beta = 0._rprec
! controller%beta(N,:) = 5._rprec
! do i = 1, N
!     controller%gen_torque(i,:) = gen_torque(i)
! end do
! ! write(*,*) gen_torque
! ! write(*,*) controller%gen_torque
! ! call controller%makeDimensionless
! call controller%run()
! ! 
! ! m = conjugate_gradient_t(controller, 1000)
! m = lbfgsb_t(controller, 300)
! call m%minimize( controller%get_control_vector() )
! ! 
! ! ! Run in the wake model
! ! call controller%makeDimensional
! ! do i = 1, controller%Nt
! !     call wm%advance(controller%beta(:,i), controller%gen_torque(:,i), controller%dt)
! !     write(*,*) wm%beta
! !     write(*,*) wm%gen_torque
! ! end do
! 
! 
! ! write(*,*) controller%gen_torque
! ! allocate(g(size(cx)))
! ! call controller%eval(cx,f,g)
! 
! 
! 
! real(rprec) :: Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! 
! call generate_splines()
! 
! write(*,*) "beta, lambda_prime, Ctp, Cpp, dCtp_dbeta, Ctp_dlambda, dCpp_dbeta, Cpp_dlambda"
! call wm_Ct_prime_spline%interp(0._rprec, 11._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(0._rprec, 11._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 0._rprec, 11._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(-51._rprec, 11._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(-51._rprec, 11._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 25._rprec, 11._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda
! call wm_Ct_prime_spline%interp(25._rprec, 11._rprec, Ctp, dCtp_dbeta, dCtp_dlambda)
! call wm_Cp_prime_spline%interp(25._rprec, 11._rprec, Cpp, dCpp_dbeta, dCpp_dlambda)
! write(*,*) 25._rprec, 11._rprec, Ctp, Cpp, dCtp_dbeta, dCtp_dlambda, dCpp_dbeta, dCpp_dlambda

real(rprec), dimension(:), allocatable :: x, v, xi, vi
integer :: N, i, Ni
real(rprec) :: dx, pi, L, dxi, x0
type(cubic_spline_t) :: sp

N = 50
allocate( x(N) )
allocate( v(N) )

pi = 4._rprec*atan(1._rprec)
L = 7._rprec*pi/8._rprec
x0 = -7._rprec*pi/16._rprec
dx = L/(N-1)
do i = 1, N
    x(i) = tan(dx * (i-1) + x0)
    v(i) = sin(x(i))
end do

write(*,*) "x = [", x, "];"
write(*,*) "v = [", v, "];"

sp = cubic_spline_t(x, v, 'not-a-knot', 'not-a-knot')!, cos(x(1)), cos(x(N)))

Ni = 256
allocate( xi(Ni) )
allocate( vi(Ni) )

L = x(N) - x(1)
dxi = 2*L/(Ni-1)
do i = 1, Ni
    xi(i) = dxi * (i-1) - L/2 + x(1)
end do

call sp%interp(xi, vi)

write(*,*) "xi = [", xi, "];"
write(*,*) "vi = [", vi, "];"

deallocate(x)
deallocate(v)
deallocate(xi)
deallocate(vi)


end program test
