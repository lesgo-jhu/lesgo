module get_wm_m
implicit none

contains

!*******************************************************************************
subroutine read_Pref(time, Pref)
!*******************************************************************************
use turbines, only : count_lines
use param, only : path
use types, only : rprec
implicit none
integer :: N, fid, i
real(rprec), dimension(:), allocatable, intent(out) :: time, Pref
character(*), parameter :: Pref_dat = path // 'input_turbines/Pref.dat'

! Create power reference
! Count number of entries and allocate
N = count_lines(Pref_dat)
allocate( time(N) )
allocate( Pref(N) )

! Read values from file
open(newunit=fid, file=Pref_dat, position='rewind', form='formatted')
do i = 1, N
    read(fid,*) time(i), Pref(i)
end do
close(fid)

end subroutine read_Pref

!*******************************************************************************
function get_wm(U_infty) result(wm)
!*******************************************************************************
use turbines, only : turbines_init, num_x, num_y, inertia_all, dia_all
use turbines, only : generate_splines, wm_Ct_prime_spline, wm_Cp_prime_spline
use turbines, only : torque_gain
use stat_defs, only : wind_farm
use param, only : z_i, nx, ny
use types, only : rprec
use input_util
use grid_m
use cubic_spline
use wake_model
implicit none

real(rprec), intent(in) :: U_infty
type(wake_model_t) :: wm
real(rprec) :: rho
real(rprec), dimension(:), allocatable :: k

! Call read lesgo's input file
call read_input_conf

! Initialize uv grid (calculate x,y,z vectors)
call grid%build()

! Initialize the turbines function
call turbines_init

rho = 1.225

allocate(k(2))
k = 0.05_rprec

wm = wake_model_t(wind_farm%turbine(1:2)%xloc*z_i, wind_farm%turbine(1:2)%yloc*z_i,&
    U_infty, 0.5*dia_all*z_i, k, dia_all*z_i, rho, inertia_all, 25, 25,      &
    wm_Ct_prime_spline, wm_Cp_prime_spline, torque_gain)

end function get_wm

end module get_wm_m

!*******************************************************************************
program mpc
!*******************************************************************************
use types, only : rprec
use wake_model
use minimize_m
use turbines, only : torque_gain
use open_file_fid_mod
use functions, only : linear_interp
use turbines_mpc
use get_wm_m

use lbfgsb
implicit none

! common variables
integer :: i, j, ii
real(rprec) :: cfl, dt

! wake model variables
real(rprec), dimension(:), allocatable :: beta, gen_torque
real(rprec) :: U_infty
integer :: Nt, Nskip

type(wake_model_t) :: wm
type(turbines_mpc_t) :: controller
! type(conjugate_gradient_t) :: m
type(lbfgsb_t) :: m

real(rprec), dimension(:), allocatable :: Pref
real(rprec), dimension(:), allocatable :: time
real(rprec), dimension(:,:), allocatable :: beta_c, torque_gain_c
real(rprec) :: tt, T
real(rprec) :: Ca, Cb
integer :: omega_fid, beta_fid, torque_gain_fid, uhat_fid, Ctp_fid, Cpp_fid
integer :: Pref_fid, Pfarm_fid, u_fid

! Create the wake model
U_infty = 9._rprec
wm = get_wm(U_infty)
allocate(beta(wm%N))
allocate(gen_torque(wm%N))
cfl = 0.1_rprec
Nt = floor(wm%Nx/cfl)

! integrate the wake model forward in time to get reference power
dt = cfl * wm%dx / wm%U_infty
do i = 1, Nt
    call wm%advance(dt)
end do

! Read the reference signal
call read_Pref(time, Pref)
Pref = Pref * sum(wm%Phat)

write(*,*) Pref

! Create controller
tt = 0._rprec
T = 2._rprec*wm%x(wm%Nx)/wm%U_infty
controller = turbines_mpc_t(wm, 0._rprec, T, 0.99_rprec, time, Pref, 0.723_rprec,  1.267_rprec)
controller%beta = 0._rprec
controller%torque_gain = wm%torque_gain(1)
call controller%makeDimensionless
write(*,*) controller%torque_gain
call controller%run
call controller%rescale_gradient
Ca = controller%Ca
Cb = controller%Cb
!
write(*,*) "Ca, Cb: ", controller%Ca, controller%Cb

! Do the initial optimization
m = lbfgsb_t(controller, 10, controller%get_lower_bound(), controller%get_upper_bound() )
call controller%run()
call m%minimize( controller%get_control_vector() )
call controller%run()
! call controller%rescale_gradient
! Ca = controller%Ca
! Cb = controller%Cb

! Allocate control vectors
allocate(beta_c(controller%N, controller%Nt))
allocate(torque_gain_c(controller%N, controller%Nt))

! open files
call system ( "mkdir -p output-mpc" )
open(newunit=omega_fid,file='output-mpc/omega.dat')
open(newunit=beta_fid,file='output-mpc/beta.dat')
open(newunit=torque_gain_fid,file='output-mpc/torque_gain.dat')
open(newunit=uhat_fid,file='output-mpc/uhat.dat')
open(newunit=Ctp_fid,file='output-mpc/Ctp.dat')
open(newunit=Cpp_fid,file='output-mpc/Cpp.dat')
open(newunit=Pref_fid,file='output-mpc/Pref.dat')
open(newunit=Pfarm_fid,file='output-mpc/Pfarm.dat')
open(newunit=u_fid,file='output-mpc/u.dat')

write(omega_fid,*) wm%omega
write(beta_fid,*) wm%beta
write(torque_gain_fid,*) wm%torque_gain
write(uhat_fid,*) wm%uhat
write(Ctp_fid,*) wm%Ctp
write(Cpp_fid,*) wm%Cpp
write(Pref_fid,*) Pref(1)
write(Pfarm_fid,*) sum(wm%Phat)
write(u_fid,*) wm%u

write(*,*) wm%omega

Nskip = 1
call controller%MakeDimensional
write(*,*) "dt = ", controller%dt

do j = 1, 3*ceiling( time(size(time)) / controller%dt ) / Nskip
    ! Copy over control vectors
    call controller%MakeDimensional
    beta_c = controller%beta
    torque_gain_c = controller%torque_gain
    ! Advance wake model
    do i = 2, Nskip+1
        tt = tt + controller%dt
        write(*,*) tt
        call wm%advance(controller%dt, beta_c(:,i), torque_gain_c(:,i))
        write(omega_fid,*) wm%omega
        write(beta_fid,*) wm%beta
        write(torque_gain_fid,*) wm%torque_gain
        write(uhat_fid,*) wm%uhat
        write(Ctp_fid,*) wm%Ctp
        write(Cpp_fid,*) wm%Cpp
        write(Pref_fid,*) controller%Pref(i)
        write(Pfarm_fid,*) sum(wm%Phat)
        write(u_fid,*) wm%u
        write(*,*) "Percent error:",                                           &
            (sum(wm%Phat) - controller%Pref(i))/controller%Pref(i)*100._rprec
    end do

    ! create controller
    controller = turbines_mpc_t(wm, 0._rprec, T, 0.99_rprec, time-tt, Pref, 0.723_rprec,  1.267_rprec)
    controller%beta(:,:controller%Nt-Nskip) = beta_c(:,Nskip+1:)
    controller%torque_gain(:,:controller%Nt-Nskip) = torque_gain_c(:,Nskip+1:)
    do ii = controller%Nt-Nskip, controller%Nt
        controller%beta(:,ii) = 0._rprec
        controller%torque_gain(:,ii) = torque_gain
    end do
    call controller%makeDimensionless
    !call controller%rescale_gradient
    controller%Ca = Ca
    controller%Cb = Cb

    ! minimize
!     m = conjugate_gradient_t(controller, 500)
    m = lbfgsb_t(controller, 10, controller%get_lower_bound(),                 &
        controller%get_upper_bound() )
    call m%minimize( controller%get_control_vector() )
    call controller%run()

end do

! close files
close(omega_fid)
close(beta_fid)
close(torque_gain_fid)
close(uhat_fid)
close(Ctp_fid)
close(Cpp_fid)
close(Pref_fid)
close(Pfarm_fid)
close(u_fid)

end program mpc
