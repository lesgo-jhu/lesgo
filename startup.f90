!*******************************************************************************
program wake_model_startup
!*******************************************************************************
use types, only : rprec
use wake_model
use minimize_m
use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
use open_file_fid_mod
use functions, only : linear_interp
use cubic_spline

implicit none

! common variables
integer :: i, j
real(rprec) :: cfl, dt

! wake model variables
type(wake_model_t) :: wm
real(rprec), dimension(:), allocatable :: sx, sy, k, beta, gen_torque
real(rprec) :: U_infty, Delta, Dia, rho, inertia, torque_gain
integer :: N, Nx, Ny, Nt, Numx, Numy

! Output files
integer, parameter :: omega_fid=1, beta_fid=2, gen_torque_fid=3, uhat_fid=4
integer, parameter :: Ctp_fid=5, Cpp_fid=60, u_fid=61, x_fid=62, y_fid=63

! initialize wake model
cfl = 0.1_rprec
Dia = 126._rprec
Delta = 0.25_rprec * Dia
rho = 1.225_rprec
inertia = 4.0469e+07_rprec
torque_gain = 2.1648e6
U_infty = 9._rprec
Numx = 4
Numy = 2
N = Numx*Numy
Nx = 128
Ny = 96
Nt = 8*Nx
! Nt = 2
allocate(sx(N))
allocate(sy(N))
allocate(k(N))
allocate(beta(N))
allocate(gen_torque(N))
k = 0.05_rprec
beta = -2._rprec
do i = 1, Numx
    do j = 1, Numy
        sx((i-1)*2+j) = 7._rprec * Dia * i
        sy((i-1)*2+j) = 5._rprec * Dia * j! + 1.5_rprec*Dia*i
    end do
end do
call generate_splines
wm = wake_model_t(sx, sy, U_infty, Delta, k, Dia, rho, inertia, Nx, Ny,        &
    wm_Ct_prime_spline, wm_Cp_prime_spline, torque_gain)

! open files
call system ( "mkdir -p output-startup" )
open(omega_fid,file='output-startup/omega.dat')
open(beta_fid,file='output-startup/beta.dat')
open(gen_torque_fid,file='output-startup/gen_torque.dat')
open(uhat_fid,file='output-startup/uhat.dat')
open(Ctp_fid,file='output-startup/Ctp.dat')
open(Cpp_fid,file='output-startup/Cpp.dat')
open(u_fid,file='output-startup/u.dat')
open(x_fid,file='output-startup/x.dat')
open(y_fid,file='output-startup/y.dat')

! integrate the wake model forward
dt = cfl * wm%dx / wm%U_infty
do i = 1, Nt
    gen_torque = torque_gain * wm%omega**2
    call wm%advance(beta, gen_torque, dt)
end do

! integrate the wake model forward
dt = cfl * wm%dx / wm%U_infty
beta = 0._rprec
do i = 1, Nt
    gen_torque = torque_gain * wm%omega**2
    call wm%advance(beta, gen_torque, dt)
    write(omega_fid,*) wm%omega
    write(beta_fid,*) wm%beta
    write(gen_torque_fid,*) wm%gen_torque
    write(uhat_fid,*) wm%uhat
    write(Ctp_fid,*) wm%Ctp
    write(Cpp_fid,*) wm%Cpp
end do

write(u_fid,*) wm%u
write(x_fid,*) wm%x
write(y_fid,*) wm%y

! close files
close(omega_fid)
close(beta_fid)
close(gen_torque_fid)
close(uhat_fid)
close(Ctp_fid)
close(Cpp_fid)
close(u_fid)
close(x_fid)
close(y_fid)

end program wake_model_startup
