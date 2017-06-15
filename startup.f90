!*******************************************************************************
program wake_model_startup
!*******************************************************************************
use types, only : rprec
use wake_model
use minimize
use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
use open_file_fid_mod
use functions, only : linear_interp
use cubic_spline
use turbines_mpc

implicit none

! common variables
integer :: i, j, ii
real(rprec) :: cfl, dt

! wake model variables
type(wake_model_t) :: wm, wmi
real(rprec), dimension(:), allocatable :: s, k, beta, gen_torque
real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
integer :: N, Nx, Nt, Nskip

! Output files
integer, parameter :: omega_fid=1, beta_fid=2, gen_torque_fid=3, uhat_fid=4
integer, parameter :: Ctp_fid=5, Cpp_fid=60, u_fid=61

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
Nt = (N+1)**2*Nx
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
    
! open files
call system ( "mkdir -p output-startup" )
open(omega_fid,file='output-startup/omega.dat')
open(beta_fid,file='output-startup/beta.dat')
open(gen_torque_fid,file='output-startup/gen_torque.dat')
open(uhat_fid,file='output-startup/uhat.dat')
open(Ctp_fid,file='output-startup/Ctp.dat')
open(Cpp_fid,file='output-startup/Cpp.dat')
open(u_fid,file='output-startup/u.dat')

! integrate the wake model forward
dt = cfl * wm%dx / wm%U_infty
write(*,*) dt
do i = 1, Nt
    gen_torque = torque_gain * wm%omega**2
    call wm%advance(beta, gen_torque, dt)
    write(omega_fid,*) wm%omega
    write(beta_fid,*) wm%beta
    write(gen_torque_fid,*) wm%gen_torque
    write(uhat_fid,*) wm%uhat
    write(Ctp_fid,*) wm%Ctp
    write(Cpp_fid,*) wm%Cpp
    write(u_fid,*) wm%u
end do

! close files
close(omega_fid)
close(beta_fid)
close(gen_torque_fid)
close(uhat_fid)
close(Ctp_fid)
close(Cpp_fid)
close(u_fid)

end program wake_model_startup