program test
use types, only : rprec
use wake_model_class
use rh_control

! common variables
integer :: i
real(rprec) :: cfl, dt

! wake model variables
type(WakeModel) :: wm
real(rprec), dimension(:), allocatable :: s, k, Ctp
real(rprec) :: U_infty, Delta, Dia
integer :: N, Nx

! minimizer
type(MinimizedFarm) :: mf
real(rprec) :: t0, T, tau
real(rprec), dimension(:), allocatable :: time, Pref
real(rprec), dimension(:,:), allocatable :: phi_in

! initialize wake model
cfl = 0.99_rprec
Dia = 100._rprec
Delta = 0.5_rprec * Dia
U_infty = 9._rprec
N = 7
Nx = 64
allocate(s(N))
allocate(k(N))
allocate(Ctp(N))
k = 0.05_rprec
Ctp = 1.33_rprec
do i = 1, N
    s(i) = 7._rprec * Dia * i
end do
wm = WakeModel(s, U_infty, Delta, k, Dia, Nx)

! integrate the wake model forward in time at least 2 flow through times
dt = cfl * wm%dx / U_infty
do i = 1, 2*wm%N
    call wm%advance(Ctp, dt)
end do

! create minimizer
t0 = 0
T = 5._rprec * 60._rprec
allocate(time(2))
allocate(Pref(2))
time(1) = 0
time(2) = T
Pref = 0.2_rprec * 1.33_rprec * U_infty**3 * N
tau = 120._rprec
mf = MinimizedFarm(wm, t0, T, cfl, time, Pref, tau)

allocate(phi_in(N, 2))
phi_in = 1.33_rprec
call mf%run(time, phi_in)

call mf%finiteDifferenceGradient
write(*,*) mf%fdgrad - mf%grad

end program test
