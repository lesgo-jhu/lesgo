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
    real(rprec) :: sigma_du, sigma_k, sigma_omega, sigma_uhat
    real(rprec) :: tau
contains
    procedure, public :: initialize_val
    procedure, public :: advance
    procedure, public :: advance_ensemble
    procedure, public :: generate_initial_ensemble
end type wake_model_estimator_t

interface wake_model_estimator_t
    module procedure :: constructor_val
!     module procedure :: constructor_file
end interface wake_model_estimator_t

contains

!*******************************************************************************
function constructor_val(i_Ne, i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia,     &
    i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline,                  &
    i_sigma_du, i_sigma_k, i_sigma_omega, i_sigma_uhat, i_tau) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_estimator_t) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_sx, i_sy, i_k
integer, intent(in) :: i_Ne, i_Nx, i_Ny
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec), intent(in) :: i_sigma_du, i_sigma_k, i_sigma_omega, i_sigma_uhat
real(rprec), intent(in) :: i_tau

call this%initialize_val(i_Ne, i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia,     &
    i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline,                  &
    i_sigma_du, i_sigma_k, i_sigma_omega, i_sigma_uhat, i_tau)
    
end function constructor_val

!*******************************************************************************
subroutine initialize_val(this, i_Ne, i_sx, i_sy, i_U_infty, i_Delta, i_k,     &
    i_Dia, i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline,           &
    i_sigma_du, i_sigma_k, i_sigma_omega, i_sigma_uhat, i_tau)
!*******************************************************************************
use param, only : coord
use grid_m
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_sx, i_sy, i_k
integer, intent(in) :: i_Ne, i_Nx, i_Ny
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec), intent(in) :: i_sigma_du, i_sigma_k, i_sigma_omega, i_sigma_uhat
real(rprec), intent(in) :: i_tau
integer :: i

! Set ensemble parameters
this%sigma_du = i_sigma_du
this%sigma_k = i_sigma_k
this%sigma_omega = i_sigma_omega
this%sigma_uhat = i_sigma_uhat
this%tau = i_tau

! Create wake model
this%wm = wake_model_t(i_sx, i_sy, i_U_infty, i_Delta, i_k,                    &
    i_Dia, i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline)
    
! Create ensemble
this%Ne = i_Ne
allocate( this%ensemble(this%Ne) )
do i = 1, this%Ne
    this%ensemble(i) = wake_model_t(i_sx, i_sy, i_U_infty, i_Delta, i_k,       &
        i_Dia, i_rho, i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline)
end do

! Allocate filter matrices
this%Nm = 2 * this%wm%N                         ! Number of measurements
this%Ns = (this%wm%N + 2) * this%wm%Nx - 1      ! Number of states
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
subroutine generate_initial_ensemble(this, torque_gain)
!*******************************************************************************
use util, only : random_normal
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), intent(in) :: torque_gain
real(rprec), parameter :: cfl = 0.99
real(rprec) :: dt, FTT
integer:: i, j, N, Nx, jstart, jend
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

! Do at least 1 FTT of simulations to get good ensemble statistics  
do i = 1, floor(FTT / dt)
    call this%advance_ensemble(beta, torque_gain*this%wm%omega**2, dt)
end do

! call advance_ensemble(this, beta, torque_gain*this%wm%omega**2, dt)
! 
! ! Place ensemble into a matrix with each member in a column
! this%Abar = 0
! this%Ahatbar = 0
! do i = 1, this%Ne
!     do j = 1, N
!         jstart = (j-1)*Nx+1
!         jend   = j*Nx
!         this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
!     end do
!     this%A(Nx*N+1:,i) = this%ensemble(i)%k(1:N-1)
!     this%Ahat(:,i)    = this%ensemble(i)%uhat
!     this%Abar         = this%Abar + this%A(:,i) / this%Ne
!     this%Ahatbar      = this%Ahatbar + this%Ahat(:,i) / this%Ne
! end do
! do j = 1, this%Ne
!     this%Aprime(:,j) = this%A(:,j) - this%Abar
!     this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
! end do
    
end subroutine generate_initial_ensemble

!*******************************************************************************
subroutine advance(this, um, omegam, beta, gen_torque, dt)
!*******************************************************************************
use sim_param, only : u
use util, only : random_normal, inverse
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), intent(in) :: dt
real(rprec), dimension(:), intent(in) :: um, omegam, beta, gen_torque
real(rprec) :: Uinftyi, alpha
integer :: i, j, N_unwaked
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
do i = N+1, 2*N
    do j = 1, this%Ne
        this%E(i, j) = this%sigma_omega * random_normal()
        this%D(i, j) = omegam(i) + this%E(i, j)
    end do
end do
this%Dprime = this%D - this%Ahat

! Update Anew = A + A'*Ahat'^T * (Ahat'*Ahat'^T + E*E^T)^-1 * D'
! Since the dimension is small, we don't bother doing the SVD. If the matrix becomes
! singular, then this should be considered as in section 4.3.2 of Everson(2003)
this%A = this%A + matmul( matmul(this%Aprime, transpose(this%Ahatprime)),      &
    matmul(inverse(matmul(this%Ahatprime, transpose(this%Ahatprime)) +         &
    matmul(this%E, transpose(this%E))), this%Dprime))
    
! Compute mean
this%Abar = 0._rprec
do i = 1, this%Ne
    this%Abar = this%Abar + this%A(:,i) / this%Ne;
end do

! Filter the freestream velocity based on unwaked turbines
Uinftyi = 0._rprec
N_unwaked = this%wm%N - this%wm%Nwaked
alpha = dt / (this%tau + dt)
do i = 1, this%wm%N
    if (.not.this%wm%waked(i)) then
        Uinftyi = Uinftyi + (4._rprec + this%wm%Ctp(i))/4._rprec*um(i)/N_unwaked
    end if
end do
this%wm%U_infty = alpha * Uinftyi + (1 - alpha) * this%wm%U_infty
this%wm%VELOCITY = this%wm%U_infty
this%wm%TIME  = this%wm%LENGTH / this%wm%VELOCITY
this%wm%TORQUE = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**2
this%wm%POWER = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**3
write(*,*) this%wm%U_infty

! Fill into objects
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend = j*Nx
        this%ensemble(i)%du(j,:)  = this%A(jstart:jend,i)
    end do
    this%ensemble(i)%U_infty  = this%wm%U_infty
    this%ensemble(i)%VELOCITY = this%wm%U_infty
    this%ensemble(i)%TIME  = this%wm%LENGTH / this%wm%VELOCITY
    this%ensemble(i)%TORQUE = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**2
    this%ensemble(i)%POWER = this%wm%MASS * this%wm%LENGTH**2 / this%wm%TIME**3
    this%ensemble(i)%omega(:) = this%A(Nx*N+1:(Nx+1)*N,i)
    this%ensemble(i)%k(1:N-1) = this%A((Nx+1)*N+1:,i)
    this%ensemble(i)%k(N) = this%ensemble(i)%k(N-1)
end do
do j = 1, N
    jstart = (j-1)*Nx+1
    jend = j*Nx
    this%wm%du(j,:) = this%Abar(jstart:jend)
end do
this%wm%omega(:) = this%Abar(Nx*N+1:(Nx+1)*N)
this%wm%k(1:N-1) = this%Abar((Nx+1)*N+1:)
this%wm%k(N) = this%wm%k(N-1)

! Advance ensemble and mean estimate
call this%advance_ensemble(beta, gen_torque, dt)

! Place ensemble into a matrix with each member in a column
this%Abar = 0
this%Ahatbar = 0
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend   = j*Nx
        this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
    end do
    this%A(Nx*N+1:Nx*(N+1),i) = this%ensemble(i)%omega(:)
    this%A(Nx*(N+1)+1:,i) = this%ensemble(i)%k(1:N-1)
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
subroutine advance_ensemble(this, beta, gen_torque, dt)
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
    do j = 1,N
        this%ensemble(i)%k(j) = max(this%ensemble(i)%k(j)                      &
            + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
        this%ensemble(i)%du(j,:) = max(this%ensemble(i)%du(j,:)                &
            + sqrt(dt) * this%sigma_du * random_normal()                       &
            * this%ensemble(i)%G(j,:) * sqrt(2*pi) * this%ensemble(i)%Delta,   &
            0._rprec)
        this%ensemble(i)%omega(j) = max(this%ensemble(i)%omega(j)              &
            + sqrt(dt) * this%sigma_omega * random_normal(), 0._rprec)
                              
    end do
    this%ensemble(i)%k(N) = this%ensemble(i)%k(N-1)
    call this%ensemble(i)%compute_wake_expansion
    call this%ensemble(i)%advance(beta, gen_torque, dt)
end do
call this%wm%compute_wake_expansion
call this%wm%advance(beta, gen_torque, dt)
    
end subroutine advance_ensemble


end module wake_model_estimator
