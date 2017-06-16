!*******************************************************************************
program mpc
!*******************************************************************************
use types, only : rprec
use wake_model
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
real(rprec) :: tt, T
real(rprec) :: Ca, Cb
integer, parameter :: omega_fid=1, beta_fid=2, gen_torque_fid=3, uhat_fid=4
integer, parameter :: Ctp_fid=5, Cpp_fid=60, Pref_fid=7, Pfarm_fid=8
integer, parameter :: alpha_fid=9, u_fid=61

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

! Create power reference
allocate(time(3))
allocate(Pref(3))
time(1) = 0._rprec
time(2) = 300._rprec
time(3) = 600._rprec
time(4) = 2000._rprec
Pref(1) = sum(wm%Phat)
Pref(2) = sum(wm%Phat)
Pref(3) = 0.75*sum(wm%Phat)
Pref(4) = 0.75*sum(wm%Phat)

write(*,*) time
write(*,*) Pref

! Create controller
tt = 0._rprec
T = 2._rprec*wm%x(wm%Nx)/wm%U_infty
controller = turbines_mpc_t(wm, 0._rprec, T, 0.99_rprec, time, Pref)
controller%beta(:,2:) = 0._rprec
controller%alpha(:,2:) = 0._rprec
call controller%makeDimensionless
call controller%rescale_gradient
Ca = controller%Ca
Cb = controller%Cb

! Do the initial optimization
! m = conjugate_gradient_t(controller, 500)
m = lbfgsb_t(controller, 10)
call m%minimize( controller%get_control_vector() )
call controller%run()
!write(*,*) "beta"
!write(*,*) controller%grad_beta
!write(*,*) "alpha"
!write(*,*) controller%grad_alpha

! Allocate control vectors
allocate(beta_c(controller%N, controller%Nt))
allocate(alpha_c(controller%N, controller%Nt))
allocate(gen_torque_c(controller%N, controller%Nt))

! open files
call system ( "mkdir -p output-mpc" )
open(omega_fid,file='output-mpc/omega.dat')
open(beta_fid,file='output-mpc/beta.dat')
open(gen_torque_fid,file='output-mpc/gen_torque.dat')
open(uhat_fid,file='output-mpc/uhat.dat')
open(Ctp_fid,file='output-mpc/Ctp.dat')
open(Cpp_fid,file='output-mpc/Cpp.dat')
open(Pref_fid,file='output-mpc/Pref.dat')
open(Pfarm_fid,file='output-mpc/Pfarm.dat')
open(alpha_fid,file='output-mpc/alpha.dat')
open(u_fid,file='output-mpc/u.dat')

write(omega_fid,*) wm%omega
write(beta_fid,*) wm%beta
write(gen_torque_fid,*) wm%gen_torque
write(uhat_fid,*) wm%uhat
write(Ctp_fid,*) wm%Ctp
write(Cpp_fid,*) wm%Cpp
write(Pref_fid,*) Pref(1)
write(Pfarm_fid,*) sum(wm%Phat)
write(alpha_fid,*) wm%Phat/wm%Paero - 1._rprec
write(u_fid,*) wm%u

Nskip = 4
do j = 1,100
    ! Copy over control vectors
    call controller%MakeDimensional
    beta_c = controller%beta
    alpha_c = controller%alpha
    gen_torque_c = controller%gen_torque

    ! Advance wake model
    do i = 2, Nskip+1
        tt = tt + controller%dt
        write(*,*) tt
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
        write(u_fid,*) wm%u
        write(*,*) "Percent error:",                                           &
            (sum(wm%Phat) - controller%Pref(i))/controller%Pref(i)*100._rprec
    end do

    ! create controller
    controller = turbines_mpc_t(wm, 0._rprec, T, 0.99_rprec, time-tt, Pref)
    controller%beta(:,:controller%Nt-Nskip) = beta_c(:,Nskip+1:)
    controller%alpha(:,:controller%Nt-Nskip) = alpha_c(:,Nskip+1:)
    do ii = controller%Nt-Nskip, controller%Nt
        controller%beta(:,ii) = beta_c(:,controller%Nt)
        controller%alpha(:,ii) = 0._rprec
    end do
    call controller%makeDimensionless
    call controller%rescale_gradient(Ca, Cb)

    ! minimize
    ! m = conjugate_gradient_t(controller, 500)
    m = lbfgsb_t(controller, 10)
    call m%minimize( controller%get_control_vector() )
    call controller%run()

end do

! close files
close(omega_fid)
close(beta_fid)
close(gen_torque_fid)
close(uhat_fid)
close(Ctp_fid)
close(Cpp_fid)
close(Pref_fid)
close(Pfarm_fid)
close(alpha_fid)
close(u_fid)
call controller%makeDimensional()
write(*,*) controller%dt

end program mpc
