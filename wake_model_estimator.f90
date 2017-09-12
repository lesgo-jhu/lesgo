!!
!!  Copyright (C) 2016  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module wake_model_estimator
!*******************************************************************************
use types, only : rprec
use util,  only : logistic, softplus, gaussian
use wake_model
use messages
use bi_pchip
use param, only : pi
implicit none

private
public wake_model_estimator_t

type :: wake_model_estimator_t
    type(wake_model_t), dimension(:), allocatable :: ensemble
    type(wake_model_t) :: wm
    integer :: Ne, Nm, Ns
    real(rprec), dimension(:,:), allocatable :: A, Aprime, Ahat, Ahatprime
    real(rprec), dimension(:,:), allocatable :: E, D, Dprime
    real(rprec), dimension(:), allocatable :: Abar, Ahatbar
    real(rprec) :: sigma_du, sigma_k, sigma_uhat
    real(rprec) :: tau
contains
    procedure, public :: initialize_val
    procedure, private :: initialize_file
    procedure, public :: write_to_file
    procedure, public :: calc_U_infty
    procedure, public :: advance
    procedure, private :: advance_ensemble_val
    procedure, private :: advance_ensemble_noval
    generic, public :: advance_ensemble => advance_ensemble_val,               &
        advance_ensemble_noval
    procedure, public :: generate_initial_ensemble
end type wake_model_estimator_t

interface wake_model_estimator_t
    module procedure :: constructor_val
    module procedure :: constructor_file
end interface wake_model_estimator_t

contains

!*******************************************************************************
function constructor_val(i_Ne, i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia,     &
    i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline, i_torque_gain,   &
    i_sigma_du, i_sigma_k, i_sigma_uhat, i_tau) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_estimator_t) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_sx, i_sy, i_k
integer, intent(in) :: i_Ne, i_Nx, i_Ny
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec), intent(in) :: i_torque_gain
real(rprec), intent(in) :: i_sigma_du, i_sigma_k, i_sigma_uhat
real(rprec), intent(in) :: i_tau

call this%initialize_val(i_Ne, i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia,     &
    i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline, i_torque_gain,   &
    i_sigma_du, i_sigma_k, i_sigma_uhat, i_tau)

end function constructor_val

!*******************************************************************************
function constructor_file(fpath, i_Ctp_spline, i_Cpp_spline, i_sigma_du,       &
    i_sigma_k, i_sigma_uhat, i_tau) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_estimator_t) :: this
character(*), intent(in) :: fpath
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec), intent(in) :: i_sigma_du, i_sigma_k, i_sigma_uhat
real(rprec), intent(in) :: i_tau

call this%initialize_file(fpath, i_Ctp_spline, i_Cpp_spline, i_sigma_du,       &
    i_sigma_k, i_sigma_uhat, i_tau)

end function constructor_file

!*******************************************************************************
subroutine initialize_val(this, i_Ne, i_sx, i_sy, i_U_infty, i_Delta, i_k,     &
    i_Dia, i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline,           &
    i_torque_gain, i_sigma_du, i_sigma_k, i_sigma_uhat, i_tau)
!*******************************************************************************
use grid_m
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_sx, i_sy, i_k
integer, intent(in) :: i_Ne, i_Nx, i_Ny
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec), intent(in) :: i_torque_gain
real(rprec), intent(in) :: i_sigma_du, i_sigma_k, i_sigma_uhat
real(rprec), intent(in) :: i_tau
integer :: i

! Set ensemble parameters
this%sigma_du = i_sigma_du
this%sigma_k = i_sigma_k
this%sigma_uhat = i_sigma_uhat
this%tau = i_tau

! Create wake model
this%wm = wake_model_t(i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia, i_rho,      &
    i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline, i_torque_gain)

! Create ensemble
this%Ne = i_Ne
allocate( this%ensemble(this%Ne) )
do i = 1, this%Ne
    this%ensemble(i) = wake_model_t(i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia,&
        i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline, i_torque_gain)
end do

! Allocate filter matrices
this%Nm = this%wm%N                             ! Number of measurements
this%Ns = (this%wm%Nx + 1) * this%wm%N          ! Number of states
allocate( this%Abar(this%Ns) )
allocate( this%Ahatbar(this%Nm) )
allocate( this%A(this%Ns, this%Ne) )
allocate( this%Aprime(this%Ns, this%Ne) )
allocate( this%Ahat(this%Nm, this%Ne) )
allocate( this%Ahatprime(this%Nm, this%Ne) )
allocate( this%E(this%Nm, this%Ne) )
allocate( this%D(this%Nm, this%Ne) )
allocate( this%Dprime(this%Nm, this%Ne) )

end subroutine initialize_val

!*******************************************************************************
subroutine initialize_file(this, fpath, i_Ctp_spline, i_Cpp_spline, i_sigma_du,&
    i_sigma_k, i_sigma_uhat, i_tau)
!*******************************************************************************
! Constructor for wake model with values given
use param, only : CHAR_BUFF_LENGTH
use string_util, only : string_splice
use open_file_fid_mod
implicit none
class(wake_model_estimator_t), intent(inout) :: this
character(*), intent(in) :: fpath
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec), intent(in) :: i_sigma_du, i_sigma_k, i_sigma_uhat
real(rprec), intent(in) :: i_tau
integer :: i, fid
character(CHAR_BUFF_LENGTH) :: fstring

! set std deviations for noise and filter time constant
this%sigma_du = i_sigma_du
this%sigma_k = i_sigma_k
this%sigma_uhat = i_sigma_uhat
this%tau = i_tau

! Read current estimator matrices
fid = open_file_fid(fpath // '/wm_est.dat', 'rewind', 'unformatted')
read(fid) this%Ne, this%Nm, this%Ns

! allocate matrices
allocate( this%Abar(this%Ns) )
allocate( this%Ahatbar(this%Nm) )
allocate( this%A(this%Ns, this%Ne) )
allocate( this%Aprime(this%Ns, this%Ne) )
allocate( this%Ahat(this%Nm, this%Ne) )
allocate( this%Ahatprime(this%Nm, this%Ne) )
allocate( this%E(this%Nm, this%Ne) )
allocate( this%D(this%Nm, this%Ne) )
allocate( this%Dprime(this%Nm, this%Ne) )

! read matrices
read(fid) this%A
read(fid) this%Aprime
read(fid) this%Ahat
read(fid) this%Ahatprime
read(fid) this%E
read(fid) this%D
read(fid) this%Dprime
read(fid) this%Abar
read(fid) this%Ahatbar
close(fid)

! write the current estimate
this%wm = wake_model_t(fpath // '/wm.dat', i_Ctp_spline, i_Cpp_spline)

! write the ensemble
allocate( this%ensemble(this%Ne) )
do i = 1, this%Ne
    call string_splice(fstring, fpath // '/ensemble_', i, '.dat')
    this%ensemble(i) = wake_model_t(fstring, i_Ctp_spline, i_Cpp_spline)
end do

end subroutine initialize_file

!*******************************************************************************
subroutine write_to_file(this, fpath)
!*******************************************************************************
use string_util, only : string_splice
use param, only : CHAR_BUFF_LENGTH
use open_file_fid_mod
implicit none
class(wake_model_estimator_t), intent(inout) :: this
character(*), intent(in) :: fpath
character(CHAR_BUFF_LENGTH) :: fstring
integer :: i, fid

! Create the folder if necessary
call system('mkdir -vp ' // fpath)

! write the current estimator matrices
fid = open_file_fid(fpath // '/wm_est.dat', 'rewind', 'unformatted')
write(fid) this%Ne, this%Nm, this%Ns
write(fid) this%A
write(fid) this%Aprime
write(fid) this%Ahat
write(fid) this%Ahatprime
write(fid) this%E
write(fid) this%D
write(fid) this%Dprime
write(fid) this%Abar
write(fid) this%Ahatbar
close(fid)

! write the current estimate
call this%wm%write_to_file(fpath // '/wm.dat')

! write the ensemble
do i = 1, this%Ne
    call string_splice( fstring, fpath // '/ensemble_', i, '.dat' )
    call this%ensemble(i)%write_to_file(fstring)
end do

end subroutine write_to_file

!*******************************************************************************
subroutine calc_U_infty(this, um, alpha)
!*******************************************************************************
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: um
real(rprec), intent(in) :: alpha
real(rprec) :: Uinftyi
integer :: i, N, N_unwaked

N = this%wm%N

! Check input size
if (size(um) /= N) then
    call error('wake_model_t.generate_initial_ensemble','um must be size N')
end if
if (alpha > 1._rprec .or. alpha < 0._rprec) then
    call error('wake_model_t.generate_initial_ensemble',                       &
        'Required: 0 <= alpha <=1')
end if

! Make dimensional and calculate the new Uinfty
! Do not change scalings for the ensemble members, because they will not
! actually be made non-dimensional with a controller
Uinftyi = 0._rprec
N_unwaked = this%wm%N - this%wm%Nwaked
do i = 1, N
    if (.not.this%wm%waked(i)) then
        Uinftyi = Uinftyi + (4._rprec + this%wm%Ctp(i))/4._rprec*um(i)/N_unwaked
    end if
end do
this%wm%U_infty = alpha * Uinftyi + (1._rprec - alpha) * this%wm%U_infty
do i = 1, this%Ne
    this%ensemble(i)%U_infty = this%wm%U_infty
end do
this%wm%VELOCITY = this%wm%U_infty
this%wm%TIME  = this%wm%LENGTH / this%wm%VELOCITY
this%wm%TORQUE = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**2
this%wm%POWER = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**3

end subroutine calc_U_infty

!*******************************************************************************
subroutine generate_initial_ensemble(this)
!*******************************************************************************
use util, only : random_normal
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), parameter :: cfl = 0.2
real(rprec) :: dt, FTT
integer:: ii, i, j, N, Nx, jstart, jend
real(rprec), dimension(:), allocatable :: beta

write(*,*) 'Generating initial wake model filter ensemble...'

! Initialize random number generator
call init_random_seed

! Calculate safe dt
dt = cfl * this%wm%dx / this%wm%U_infty
FTT = this%wm%x(this%wm%Nx) / this%wm%U_Infty

! set integer
N = this%wm%N
Nx = this%wm%Nx

! create dummy beta array
allocate( beta(N) )
beta = 0._rprec

! Do 1 FFT of k's
do i = 1, this%Ne
    do ii = 1, floor(FTT / dt)
        do j = 1, N
            this%ensemble(i)%k(j) = max(this%ensemble(i)%k(j)                  &
                + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
        end do
    end do
    call this%ensemble(i)%compute_wake_expansion
end do

! Do at least 1 FTT of simulations to get good ensemble statistics
do ii = 1, floor(FTT / dt)
    write(*,*) "Advancing ensemble at time step", ii
    ! Advance step for objects
    ! always safeguard against negative k's, du's, and omega's
    do i = 1, this%Ne
        do j = 1, N
      !      this%ensemble(i)%k(j) = max(this%ensemble(i)%k(j)                      &
      !          + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
            this%ensemble(i)%du(j,:) = max(this%ensemble(i)%du(j,:)            &
                + sqrt(dt) * this%sigma_du * random_normal()                   &
                * this%ensemble(i)%G(j,:) * sqrt(2*pi)                         &
                * this%ensemble(i)%Delta, 0._rprec)
        end do
        call this%ensemble(i)%advance(dt)
    end do
    call this%wm%advance(dt)
end do

! Place ensemble into a matrix with each member in a column
this%Abar = 0._rprec
this%Ahatbar = 0._rprec
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend = j*Nx
        this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
    end do
    this%A((N*Nx+1):,i) = this%ensemble(i)%k(:)
    this%Ahat(:,i) = this%ensemble(i)%uhat
    this%Abar = this%Abar + this%A(:,i) / this%Ne
    this%Ahatbar = this%Ahatbar + this%Ahat(:,i) / this%Ne
end do
do j = 1, this%Ne
    this%Aprime(:,j) = this%A(:,j) - this%Abar
    this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
end do

write(*,*) "Done generating initial ensemble"

end subroutine generate_initial_ensemble

!*******************************************************************************
subroutine advance(this, um, omegam, beta, gen_torque, dt)
!*******************************************************************************
use util, only : random_normal, inverse
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), intent(in) :: dt
real(rprec), dimension(:), intent(in) :: um, omegam, beta, gen_torque
real(rprec) :: alpha
integer :: i, j
integer :: N, Nx
integer :: jstart, jend

N = this%wm%N
Nx = this%wm%Nx

! Check size of inputs
if (size(um) /= N) then
    call error('wake_model_t.advance','um must be size N')
end if
if (size(omegam) /= N) then
    call error('wake_model_t.advance','omegam must be size N')
end if
if (size(beta) /= N) then
    call error('wake_model_t.advance','beta must be size N')
end if
if (size(gen_torque) /= N) then
    call error('wake_model_t.advance','gen_torque must be size N')
end if

! Calculate noisy measurements
this%E = 0._rprec
do i = 1, N
    do j = 1, this%Ne
        this%E(i, j) = this%sigma_uhat * random_normal()
        this%D(i, j) = um(i) + this%E(i, j)
    end do
end do
this%Dprime = this%D - this%Ahat

! Update Anew = A + A'*Ahat'^T * (Ahat'*Ahat'^T + E*E^T)^-1 * D'
! Since the dimension is small, we don't bother doing the SVD. If the matrix
! becomes singular, then this should be considered as in section 4.3.2 of
! Everson(2003)
this%A = this%A + matmul( matmul(this%Aprime, transpose(this%Ahatprime)),      &
    matmul(inverse(matmul(this%Ahatprime, transpose(this%Ahatprime)) +         &
    matmul(this%E, transpose(this%E))), this%Dprime))

! Protect states from being negative
do i = 1, this%Ne
    do j = 1, this%Ns
        this%A(j,i) = max(this%A(j,i), 0._rprec)
    end do
end do

! Compute mean
this%Abar = 0._rprec
do i = 1, this%Ne
    this%Abar = this%Abar + this%A(:,i) / this%Ne;
end do

! Filter the freestream velocity based on unwaked turbines
alpha = dt / (this%tau + dt)
call this%calc_U_infty(um, alpha)

! Place omega measurement into the wake model and each ensemble member
this%wm%omega = omegam
do i = 1, this%Ne
    this%ensemble(i)%omega = omegam
end do

! Fill into objects
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend = j*Nx
        this%ensemble(i)%du(j,:)  = max(this%A(jstart:jend,i), 0._rprec)
        this%ensemble(i)%k(j) = max(this%A(Nx*N+j,i), 0._rprec)
    end do
    this%ensemble(i)%U_infty  = this%wm%U_infty
    this%ensemble(i)%VELOCITY = this%wm%U_infty
    this%ensemble(i)%TIME  = this%wm%LENGTH / this%wm%VELOCITY
    this%ensemble(i)%TORQUE = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**2
    this%ensemble(i)%POWER = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**3
end do
do j = 1, N
    jstart = (j-1)*Nx+1
    jend = j*Nx
    this%wm%du(j,:) = max(this%Abar(jstart:jend), 0._rprec)
end do
this%wm%k(:) = max(this%Abar((N*Nx+1):), 0._rprec)

! Fix previous timestep outputs based on corrected states
!call this%advance_ensemble(beta, gen_torque, 0._rprec)
! Advance ensemble and mean estimate
call this%advance_ensemble(beta, gen_torque, dt)

! Place ensemble into a matrix with each member in a column
this%Abar = 0
this%Ahatbar = 0
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend = j*Nx
        this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
    end do
    this%A((N*Nx+1):,i) = this%ensemble(i)%k(:)
    this%Ahat(:,i) = this%ensemble(i)%uhat
    this%Abar = this%Abar + this%A(:,i) / this%Ne
    this%Ahatbar = this%Ahatbar + this%Ahat(:,i) / this%Ne
end do
do j = 1, this%Ne
    this%Aprime(:,j) = this%A(:,j) - this%Abar
    this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
end do

end subroutine advance

!*******************************************************************************
subroutine advance_ensemble_val(this, beta, gen_torque, dt)
!*******************************************************************************
use param, only : pi
use util, only  : random_normal
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: beta, gen_torque
real(rprec), intent(in) :: dt
integer :: i, j, N

! Advance step for objects
! always safeguard against negative k's, du's, and omega's
N = this%wm%N
do i = 1, this%Ne
    do j = 1, N
        this%ensemble(i)%k(j) = max(this%ensemble(i)%k(j)                      &
            + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
        this%ensemble(i)%du(j,:) = max(this%ensemble(i)%du(j,:)                &
            + sqrt(dt) * this%sigma_du * random_normal()                       &
            * this%ensemble(i)%G(j,:) * sqrt(2*pi) * this%ensemble(i)%Delta,   &
            0._rprec)
!        this%ensemble(i)%omega(j) = max(this%ensemble(i)%omega(j)              &
!            + sqrt(dt) * 10 * this%sigma_omega * random_normal(), 0._rprec)
    end do
    call this%ensemble(i)%compute_wake_expansion
    call this%ensemble(i)%advance(beta, gen_torque, dt)
end do
call this%wm%compute_wake_expansion
call this%wm%advance(beta, gen_torque, dt)

end subroutine advance_ensemble_val

!*******************************************************************************
subroutine advance_ensemble_noval(this, dt)
!*******************************************************************************
use param, only : pi
use util, only  : random_normal
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), intent(in) :: dt
integer :: i, j, N

! Advance step for objects
! always safeguard against negative k's, du's, and omega's
N = this%wm%N
do i = 1, this%Ne
    do j = 1, N
        this%ensemble(i)%k(j) = max(this%ensemble(i)%k(j)                      &
            + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
        this%ensemble(i)%du(j,:) = max(this%ensemble(i)%du(j,:)                &
            + sqrt(dt) * this%sigma_du * random_normal()                       &
            * this%ensemble(i)%G(j,:) * sqrt(2*pi) * this%ensemble(i)%Delta,   &
            0._rprec)
!        this%ensemble(i)%omega(j) = max(this%ensemble(i)%omega(j)              &
!            + sqrt(dt) * 10 * this%sigma_omega * random_normal(), 0._rprec)
    end do
    call this%ensemble(i)%advance(dt)
end do
call this%wm%advance(dt)

end subroutine advance_ensemble_noval

end module wake_model_estimator
