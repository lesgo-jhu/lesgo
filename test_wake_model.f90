program test

use types, only : rprec
use wake_model
use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
use open_file_fid_mod
use functions, only : linear_interp
use cubic_spline
use turbines_mpc

! use lbfgsb
! use minimize
! use util, only : rosenbrock
implicit none

! common variables
integer :: i
real(rprec) :: cfl, dt

! wake model variables
type(wake_model_t) :: wm
! type(wake_model_adjoint_t) :: wma
real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
! real(rprec), dimension(:,:,:), allocatable :: fstar
! real(rprec), dimension(:,:), allocatable :: Adu, Aw, Bj, Bdu, Bw
! allocate(fstar(Nt,N,Nx))
! allocate(Adu(Nt,N))
! allocate(Aw(Nt,N))
! allocate(Bj(Nt,N))
! allocate(Bdu(Nt,N))
! allocate(Bw(Nt,N))
integer :: N, Nx, Nt

type(turbines_mpc_t) :: controller

real(rprec), dimension(:), allocatable :: Pref
real(rprec), dimension(:), allocatable :: time

allocate(time(2))
allocate(Pref(2))
time(1) = 0._rprec
time(2) = 60._rprec
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
cfl = 0.99_rprec
Dia = 126._rprec
Delta = 0.5_rprec * Dia
rho = 1.225_rprec
inertia = 4.0469e+07_rprec
torque_gain = 2.1648e6
U_infty = 9._rprec
N = 7
Nx = 256
Nt = 5*Nx
allocate(s(N))
allocate(k(N))
allocate(beta(N))
allocate(gen_torque(N))

k = 0.05_rprec
beta = -1._rprec
do i = 1, N
    s(i) = 7._rprec * Dia * i
end do
call generate_splines
wm = wake_model_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,                 &
    wm_Ct_prime_spline, wm_Cp_prime_spline)

! integrate the wake model forward in time
dt = cfl * wm%dx / wm%U_infty
do i = 1, Nt
    gen_torque = torque_gain * wm%omega**2
    call wm%advance(beta, gen_torque, dt)
end do

! Create controller, set input values, and make dimensionless
Pref = 0.9*sum(wm%Phat)
controller = turbines_mpc_t(wm, 0._rprec, time(2), 0.99_rprec, time, Pref)
controller%beta = -1._rprec
do i = 1,controller%Nt
    controller%gen_torque(:,i) = gen_torque
end do
call controller%makeDimensionless

! run with adjoints and compare to finite difference gradient
call controller%run()
call controller%finite_difference_gradient
write(*,*) "grad_beta:", controller%grad_beta
write(*,*) "fdgrad_beta:", controller%fdgrad_beta
write(*,*) "grad_gen_torque:", controller%grad_gen_torque
write(*,*) "fdgrad_gen_torque:", controller%fdgrad_gen_torque
! write(*,*) controller%fdgrad_gen_torque

end program test
