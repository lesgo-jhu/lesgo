
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
integer :: i, j, ii
real(rprec) :: cfl, dt

! wake model variables
type(wake_model_t) :: wm, wmi
real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
integer :: N, Nx, Nt, Nskip

type(turbines_mpc_t) :: controller
! type(conjugate_gradient_t) :: m
type(lbfgsb_t) :: m

real(rprec), dimension(:), allocatable :: Pref
real(rprec), dimension(:), allocatable :: time
real(rprec), dimension(:,:), allocatable :: beta_c, alpha_c, gen_torque_c
integer, parameter :: omega_fid=1, beta_fid=2, gen_torque_fid=3, uhat_fid=4
integer, parameter :: Ctp_fid=5, Cpp_fid=60, Pref_fid=7, Pfarm_fid=8, alpha_fid=9

! initialize wake model
cfl = 0.01_rprec
Dia = 126._rprec
Delta = 0.5_rprec * Dia
rho = 1.225_rprec
inertia = 4.0469e+07_rprec
torque_gain = 2.1648e6
U_infty = 9._rprec
N = 7
Nx = 128
Nt = N**2*Nx
allocate(s(N))
allocate(k(N))
allocate(beta(N))
allocate(gen_torque(N))

k = 0.05_rprec
beta = 3._rprec
do i = 1, N
    s(i) = 7._rprec * Dia * i
end do
call generate_splines
wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
    wm_Ct_prime_spline, wm_Cp_prime_spline)

! integrate the wake model forward in time to get reference power
dt = cfl * wm%dx / wm%U_infty
do i = 1, Nt
    gen_torque = 1.9*torque_gain * wm%omega**2
    call wm%advance(beta, gen_torque, dt)
end do

! Create power reference
allocate(time(2))
allocate(Pref(2))
time(1) = 0._rprec
time(2) = 2._rprec*wm%x(wm%Nx)/wm%U_infty
Pref = 1.1*sum(wm%Phat)

! Create controller
controller = turbines_mpc_t(wm, 0._rprec, time(2), 0.99_rprec, time, Pref)
controller%beta(:,2:) = 3._rprec
controller%alpha(:,2:) = 0._rprec
call controller%makeDimensionless
call controller%run()

! Do the initial optimization
! m = conjugate_gradient_t(controller, 500)
m = lbfgsb_t(controller, 10)
call m%minimize( controller%get_control_vector() )
call controller%run()

! Allocate control vectors
allocate(beta_c(controller%N, controller%Nt))
allocate(alpha_c(controller%N, controller%Nt))
allocate(gen_torque_c(controller%N, controller%Nt))

! open files
open(omega_fid,file='omega.dat')
open(beta_fid,file='beta.dat')
open(gen_torque_fid,file='gen_torque.dat')
open(uhat_fid,file='uhat.dat')
open(Ctp_fid,file='Ctp.dat')
open(Cpp_fid,file='Cpp.dat')
open(Pref_fid,file='Pref.dat')
open(Pfarm_fid,file='Pfarm.dat')
open(alpha_fid,file='alpha.dat')

write(omega_fid,*) wm%omega
write(beta_fid,*) wm%beta
write(gen_torque_fid,*) wm%gen_torque
write(uhat_fid,*) wm%uhat
write(Ctp_fid,*) wm%Ctp
write(Cpp_fid,*) wm%Cpp
write(Pref_fid,*) Pref(1)
write(Pfarm_fid,*) sum(wm%Phat)
write(alpha_fid,*) wm%Phat/wm%Paero - 1._rprec

Nskip = 10
do j = 1, 20
    ! Copy over control vectors
    call controller%MakeDimensional
    beta_c = controller%beta
    alpha_c = controller%alpha
    gen_torque_c = controller%gen_torque

    ! Advance wake model
    do i = 2, Nskip+1
        call wm%advance(beta_c(:,i), gen_torque_c(:,i), controller%dt)
        write(omega_fid,*) wm%omega
        write(beta_fid,*) wm%beta
        write(gen_torque_fid,*) wm%gen_torque
        write(uhat_fid,*) wm%uhat
        write(Ctp_fid,*) wm%Ctp
        write(Cpp_fid,*) wm%Cpp
        write(Pref_fid,*) controller%Pref(i)
        write(Pfarm_fid,*) sum(wm%Phat)
        write(alpha_fid,*) wm%Phat/wm%Paero - 1._rprec
        write(*,*) "Percent error:",                                           &
            (sum(wm%Phat) - controller%Pref(i))/controller%Pref(i)*100._rprec
    end do

    ! create controller
    controller = turbines_mpc_t(wm, 0._rprec, time(2), 0.99_rprec, time, Pref)
    controller%beta(:,:controller%Nt-Nskip) = beta_c(:,Nskip+1:)
    controller%alpha(:,:controller%Nt-Nskip) = alpha_c(:,Nskip+1:)
    do ii = controller%Nt-Nskip, controller%Nt
        controller%beta(:,ii) = beta_c(:,controller%Nt)
        controller%alpha(:,ii) = 0._rprec
    end do
    call controller%makeDimensionless
    call controller%run()

    ! minimize
    ! m = conjugate_gradient_t(controller, 500)
    m = lbfgsb_t(controller, 10)
    call m%minimize( controller%get_control_vector() )
    call controller%run()

end do

! close files
open(omega_fid)
open(beta_fid)
open(gen_torque_fid)
open(uhat_fid)
open(Ctp_fid)
open(Cpp_fid)
open(Pref_fid)
open(Pfarm_fid)
open(alpha_fid)

write(*,*) controller%dt
! write(*,*) "beta:", controller%beta
! write(*,*) "alpha", controller%alpha
! write(*,*) "Pref", controller%Pref
! write(*,*) "Pfarm", controller%Pfarm
! call controller%finite_difference_gradient()
! open(1,file='out.dat')
! write(1,*) controller%grad_beta
! write(1,*) controller%grad_alpha
! write(1,*) controller%fdgrad_beta
! write(1,*) controller%fdgrad_alpha
! close(1)
! write(*,*) controller%dt
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
