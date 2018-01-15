program test
use types, only : rprec
use wake_model_class
!use rh_control

use open_file_fid_mod
use param, only: CHAR_BUFF_LENGTH
use string_util, only : string_splice

! common variables
integer :: i, j, loc, scl, m
real(rprec) :: cfl, dt

! wake model variables
type(WakeModel) :: wm
real(rprec), dimension(:), allocatable      :: k, Ctp1, Ctp2, dval
real(rprec), dimension(:,:), allocatable    :: s, Ctp, u_ev, gridx
real(rprec), dimension(:,:,:), allocatable  :: u_tot
real(rprec)                                 :: U_infty, Delta, Dia, diff, ndiff
integer                             :: N, Nx, Ny, fid, xspace, row, col
character(CHAR_BUFF_LENGTH)                 :: name_ctp, name_uev
character(4)                                :: save_file

!! minimizer
!type(MinimizedFarm) :: mf
!real(rprec) :: t0, T, tau
!real(rprec), dimension(:), allocatable :: time, Pref
!real(rprec), dimension(:,:), allocatable :: phi_in

! initialize wake model
cfl = 0.5_rprec
Dia = 126._rprec
Delta = 0.5_rprec * Dia
U_infty = 8._rprec

N = 84
row = 7
col = 12
Ny = 210
Nx = Ny
scl = 13
allocate( s(N, 2) )
allocate( k(N) )
allocate( Ctp1(N) )
allocate( Ctp2(N) )
allocate( dval(Nx) )

allocate( Ctp( scl*N, N) )
allocate( u_ev( scl*N, N) )
allocate( u_tot( scl*N, Nx, Ny) )
allocate( gridx( Nx, Ny) )

k = 0.05_rprec
Ctp1 = 0.0_rprec
Ctp2 = 1.33_rprec
!s(:,1) = [900_rprec, 900_rprec, 900_rprec, 900_rprec, 900_rprec,    &
!         1620_rprec, 1620_rprec, 1620_rprec, 1620_rprec, 1620_rprec, &
!         1260_rprec, 1260_rprec, 1260_rprec, 1260_rprec, 1260_rprec, &
!         1980_rprec, 1980_rprec, 1980_rprec, 1980_rprec, 1980_rprec, &
!         2500_rprec, 2800_rprec, 3200_rprec, 2600_rprec, 3500_rprec, &
!         3300_rprec, 3700_rprec, 4000_rprec, 4000_rprec, 4400_rprec, &
!         4400_rprec, 4700_rprec, 5200_rprec, 4700_rprec]
loc = 0
do i = 1, row
   do j = 1, col 
      s(j+loc,1) = 7._rprec * Dia * i
      s(j+loc,2) = 5._rprec * Dia * j
   enddo
   loc =  12 * i
enddo

!s(row+1:2*row,1) = s(1:row,1)
!s(:,2) = [2100_rprec, 2780_rprec, 3460_rprec, 4140_rprec, 4820_rprec, &
!         2100_rprec, 2780_rprec, 3460_rprec, 4140_rprec, 4820_rprec,  &
!         2440_rprec, 3120_rprec, 3800_rprec, 4480_rprec, 5160_rprec,  &
!         2440_rprec, 3120_rprec, 3800_rprec, 4480_rprec, 5160_rprec,  &
!         2340_rprec, 3360_rprec, 2900_rprec, 4140_rprec, 4400_rprec,  &
!         3560_rprec, 4900_rprec, 2280_rprec, 3100_rprec, 2600_rprec,  &
!         3900_rprec, 4300_rprec, 3380_rprec, 5000_rprec]
!s(:,2) = s(:,2) - 600
!s(row+1:2*row,2) = 2000._rprec
!s(1:row,2) = 1000._rprec

wm = WakeModel(s, U_infty, Delta, k, Dia, Nx, Ny)

!  integrate the wake model forward in time at least 2 flow through times
!  (start-up case: Ctp changes)
dt = cfl * wm%dx / U_infty
name_ctp = 'Ctp1d_test.dat'
name_uev = 'uhat2d_test.dat'
!
do i = 1, 2*wm%N
    call wm%advance(Ctp1, dt)
!    print *, Ctp1
    Ctp(i,1:N) = Ctp1
    u_ev(i,1:N) = wm%uhat
    u_tot(i, :, :) = wm%u
!    write(save_file,'(i3)') i
!    open(unit=20, file='vel_movie'//trim(save_file)//'.dat', form='formatted', action='write')
!         do m = 1,wm%Nx
!             write(20, "(210(2x,F6.4))") (u_tot(i,m,j), j = 1,Ny)
!         end do
!    close(20)
    print *, Ctp(i,1)
end do
do i = 1, (scl-2)*wm%N
    call wm%advance(Ctp2, dt)
    Ctp(i+2*wm%N,1:N) = Ctp2
    u_ev(i+2*wm%N,1:N) = wm%uhat
    u_tot(i+2*wm%N,:, :) = wm%u
!    write(save_file,'(i3)') i+2*wm%N
!    open(unit=20, file='vel_movie'//trim(save_file)//'.dat', form='formatted', action='write')
!         do m = 1,wm%Nx
!             write(20, "(210(2x,F6.4))") (u_tot(i+2*wm%N,m,j), j = 1,Ny)
!         end do
!    close(20)
    print *, Ctp(i+2*wm%N,1)
end do


open(unit=8, file='vel_field.dat', form = 'formatted', action='write')
  do i = 1,wm%Nx
    write(8, "(210(2x,F6.4))") (u_tot(scl*N,i,j), j = 1,Ny)
  end do
close(8)

open(unit=3, file=name_uev, form='formatted', action='write')
  do i = 1, scl*wm%N
     write(3,"(84(2x,F8.5))") (u_ev(i,j), j = 1,84)
  end do
close(3)

print *, cfl, dt, scl*wm%N / dt
print *, wm%dy, wm%dx
end program test
