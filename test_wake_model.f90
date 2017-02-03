program test

use types, only : rprec
use wake_model
use wake_model_adjoint
use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
! use rh_control
use open_file_fid_mod
use functions, only : linear_interp
use cubic_spline
implicit none

! common variables
integer :: i, j
real(rprec) :: cfl, dt

! wake model variables
type(wake_model_t) :: wm
type(wake_model_adjoint_t) :: wma
real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
integer :: N, Nx, Nt

real(rprec) :: Pref

! adjoint variables
real(rprec), dimension(:,:,:), allocatable :: fstar
real(rprec), dimension(:,:), allocatable :: Adu, Aw, Bj, Bdu, Bw

!
! ! minimizer
! type(MinimizedFarm) :: mf
! real(rprec) :: t0, T, tau
! real(rprec), dimension(:), allocatable :: time, Pref
! real(rprec), dimension(:,:), allocatable :: phi_in

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
Nt = 5!*Nx
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
wma = wake_model_adjoint_t(s, U_infty, Delta, k, Dia, rho, inertia, Nx,        &
    wm_Ct_prime_spline, wm_Cp_prime_spline)

! adjoint variables
allocate(fstar(Nt, N, Nx))
allocate(Adu(Nt, N))
allocate(Aw(Nt, N))
allocate(Bj(Nt, N))
allocate(Bdu(Nt, N))
allocate(Bw(Nt, N))
fstar = 0._rprec
Adu = 0._rprec
Aw = 0._rprec
Bj = 0._rprec
Bdu = 0._rprec
Bw = 0._rprec

call wm%makeDimensionless

! integrate the wake model forward in time at least 2 flow through times
dt = cfl * wm%dx / wm%U_infty
do i = 1, Nt
    do j = 1, wm%N
        gen_torque = torque_gain * wm%omega**2 / wm%TIME**2 / wm%TORQUE
    end do
    call wm%advance(beta, gen_torque, dt)
    call wm%adjoint_values(Pref, fstar(i,:,:), Adu(i,:), Aw(i,:), Bj(i,:), Bdu(i,:), Bw(i,:))
end do

! integrate adjoints backwards in time
do i = 1, Nt
    call wma%retract(fstar(i,:,:), Adu(i,:), Aw(i,:), Bj(i,:), Bdu(i,:), Bw(i,:), dt)
    write(*,*) wma%uhat_star
end do


!
! write(*,*) "Adu = ", Adu(Nt,:)
! write(*,*) "Aw = ", Aw(Nt,:)
! write(*,*) "Bj = ", Bj(Nt,:)
! write(*,*) "Bdu = ", Bdu(Nt,:)
! write(*,*) "Bw = ", Bw(Nt,:)
!
! deallocate(fstar)
! deallocate(Adu)
! deallocate(Aw)
! deallocate(Bj)
! deallocate(Bdu)
! deallocate(Bw)
!
! write(*,*) "int(G) = ", sum(wm%G,2)*wm%dx
!
! ! call wm%makeDimensional
! write(*,*) "uhat = ", wm%uhat
! write(*,*) "gen_torque = ", wm%gen_torque
! write(*,*) "Paero = ", wm%Paero
! write(*,*) "omega = ", wm%omega
! write(*,*) "J = ", wm%inertia
!
! write(*,*) "TIME = ", wm%TIME
! write(*,*) "LENGTH = ", wm%LENGTH
! write(*,*) "VELOCITY = ", wm%VELOCITY
! write(*,*) "MASS = ", wm%MASS
! write(*,*) "TORQUE = ", wm%TORQUE
! write(*,*) "POWER = ", wm%POWER

!
! ! create minimizer
! t0 = 0
! T = 5._rprec * 60._rprec
! allocate(time(2))
! allocate(Pref(2))
! time(1) = 0
! time(2) = T
! Pref = 0.2_rprec * 1.33_rprec * U_infty**3 * N
! tau = 120._rprec
! mf = MinimizedFarm(wm, t0, T, cfl, time, Pref, tau)
!
! allocate(phi_in(N, 2))
! phi_in = 1.33_rprec
! call mf%run(time, phi_in)
!
! call mf%finiteDifferenceGradient
! write(*,*) mf%fdgrad - mf%grad

end program test
