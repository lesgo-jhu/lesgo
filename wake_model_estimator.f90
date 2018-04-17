!!
!!  Copyright (C) 2018  Johns Hopkins University
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
use messages
use wake_model
implicit none

private
public wake_model_estimator_t

type :: wake_model_estimator_t
    type(wake_model_t), dimension(:), allocatable :: ensemble
    type(wake_model_t)                            :: wm
    real(rprec)                                :: sigma_du, sigma_k, sigma_Phat
    integer                                    :: Ne, Nm, Ns
    real(rprec), dimension(:,:), allocatable   :: A, Aprime, Ahat, Ahatprime, E, D, Dprime
    real(rprec), dimension(:), allocatable     :: Abar, Ahatbar
    real(rprec)                                :: tau_U_infty = 300
contains
    procedure, private :: initialize_val
    procedure, private :: initialize_file
    procedure, public  :: write_to_file
    procedure, public  :: generateInitialEnsemble
    procedure, public  :: advance
    procedure, private :: advanceEnsemble
end type wake_model_estimator_t

interface wake_model_estimator_t
    module procedure :: constructor_val
    module procedure :: constructor_file
end interface wake_model_estimator_t

contains

! Constructor for wake model
!*******************************************************************************
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ne,                &
                         i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau) result(this)
!*******************************************************************************
implicit none
type(wake_model_estimator_t)              :: this
real(rprec), intent(in)               :: i_U_infty, i_Delta, i_Dia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in)                   :: i_Nx
integer, intent(in)                   :: i_Ne
real(rprec), intent(in)               :: i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau

call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ne,            &
                         i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)

end function constructor_val

! Constructor for wake model that reads from file
!*******************************************************************************
function constructor_file(fpath, i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau) result(this)
!*******************************************************************************
use param, only : CHAR_BUFF_LENGTH
implicit none

type(wake_model_estimator_t) :: this
character(*), intent(in) :: fpath
real(rprec), intent(in)  :: i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau

call this%initialize_file(fpath, i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)

end function constructor_file

!*******************************************************************************
subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ne,         &
                          i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)
!*******************************************************************************
implicit none
class(wake_model_estimator_t), intent(inout)    :: this
real(rprec), intent(in)                     :: i_sigma_du, i_U_infty, i_Delta, i_Dia
real(rprec), dimension(:), intent(in)       :: i_s, i_k
integer, intent(in)                         :: i_Nx
integer, intent(in)                         :: i_Ne
real(rprec), intent(in)                     :: i_sigma_k, i_sigma_Phat, i_tau
integer                                     :: i

! Set std deviations for noise
this%sigma_du   = i_sigma_du
this%sigma_k    = i_sigma_k
this%sigma_Phat = i_sigma_Phat

! Filter time for U_infty
this%tau_U_infty = i_tau

! Create ensemble members
this%Ne = i_Ne
allocate( this%ensemble(this%Ne) )
do i = 1, this%Ne
    this%ensemble(i) = wake_model_t(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)
end do

! Create wake model estimate
this%wm = wake_model_t(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)

! Allocate filter matrices
this%Nm = this%wm%N                              ! Number of measurements
this%Ns = this%wm%N * this%wm%Nx + this%wm%N - 1 ! Number of states
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
subroutine initialize_file(this, fpath, i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)
!*******************************************************************************
use param, only : CHAR_BUFF_LENGTH
use string_util, only : string_splice
implicit none

class(wake_model_estimator_t)   :: this
character(*), intent(in)    :: fpath
character(CHAR_BUFF_LENGTH) :: fstring
integer                     :: i, fid
real(rprec), intent(in)     :: i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau

! Set std deviations for noise
this%sigma_du   = i_sigma_du
this%sigma_k    = i_sigma_k
this%sigma_Phat = i_sigma_Phat

! Filter time for U_infty
this%tau_U_infty = i_tau

fstring = fpath // '/wm_est.dat'
open(newunit=fid, file=fstring, status='unknown', form='unformatted',        &
    position='rewind')
read(fid) this%Ne, this%Nm, this%Ns

allocate( this%Abar(this%Ns) )
allocate( this%Ahatbar(this%Nm) )
allocate( this%A(this%Ns, this%Ne) )
allocate( this%Aprime(this%Ns, this%Ne) )
allocate( this%Ahat(this%Nm, this%Ne) )
allocate( this%Ahatprime(this%Nm, this%Ne) )
allocate( this%E(this%Nm, this%Ne) )
allocate( this%D(this%Nm, this%Ne) )
allocate( this%Dprime(this%Nm, this%Ne) )

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

fstring = fpath // '/wm.dat'
this%wm = wake_model_t(fstring)

allocate( this%ensemble(this%Ne) )
do i = 1, this%Ne
    call string_splice( fstring, fpath // '/ensemble_', i, '.dat' )
    this%ensemble(i) = wake_model_t(fstring)
end do

end subroutine initialize_file

!*******************************************************************************
subroutine write_to_file(this, fpath)
!*******************************************************************************
use param, only : CHAR_BUFF_LENGTH
use string_util, only : string_splice
implicit none

class(wake_model_estimator_t)   :: this
character(*), intent(in)    :: fpath
character(CHAR_BUFF_LENGTH) :: fstring
integer                     :: i, fid

call system('mkdir -vp ' // fpath)
fstring = fpath // '/wm_est.dat'
open(newunit=fid, file=fstring, status='unknown', form='unformatted',      &
    position='rewind')
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

fstring = fpath // '/wm.dat'
call this%wm%write_to_file(fstring)

do i = 1, this%Ne
    call string_splice( fstring, fpath // '/ensemble_', i, '.dat' )
    call this%ensemble(i)%write_to_file(fstring)
end do

end subroutine write_to_file

!*******************************************************************************
subroutine generateInitialEnsemble(this, Ctp)
!*******************************************************************************
use util, only : random_normal
implicit none
class(wake_model_estimator_t), intent(inout)    :: this
real(rprec), dimension(:)                   :: Ctp
real(rprec)                                 :: dt
real(rprec), parameter                      :: cfl = 0.99
real(rprec)                                 :: FTT
integer                                     :: i, j, N, Nx, jstart, jend

if (size(Ctp) /= this%ensemble(1)%N) then
    call error('wake_model_estimator_t.generateInitialEnsemble','Ctp must be of size N')
end if

write(*,*) 'Generating initial wake model filter ensemble...'

! Initialize random number generator
call init_random_seed

! Calculate safe dt.
dt          = cfl * this%wm%dx / this%wm%U_infty
FTT         = this%wm%x(this%wm%Nx) / this%wm%U_Infty

! Do at least 1 FTT of simulations to get good ensemble statistics
N = this%wm%N
do i = 1, floor(FTT / dt)
    call this%advanceEnsemble(Ctp, dt)
end do

! Place ensemble into a matrix with each member in a column
this%Abar    = 0
this%Ahatbar = 0
N  = this%wm%N
Nx = this%wm%Nx
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend   = j*Nx
        this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
    end do
    this%A(Nx*N+1:,i) = this%ensemble(i)%k(1:N-1)
    this%Ahat(:,i)    = this%ensemble(i)%Phat
    this%Abar         = this%Abar + this%A(:,i) / this%Ne
    this%Ahatbar      = this%Ahatbar + this%Ahat(:,i) / this%Ne
end do
do j = 1, this%Ne
    this%Aprime(:,j) = this%A(:,j) - this%Abar
    this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
end do

end subroutine generateInitialEnsemble

!*******************************************************************************
subroutine advance(this, dt, Pm, Ctp)
!*******************************************************************************
use util, only : random_normal, inverse
implicit none
class(wake_model_estimator_t), intent(inout)    :: this
real(rprec), intent(in)                     :: dt
real(rprec), dimension(:), intent(in)       :: Pm, Ctp
integer                                     :: i, j, N, Nx, jstart, jend
real(rprec)                                 :: alpha, Uinftyi, r1_if

N = this%wm%N
Nx = this%wm%Nx

if (size(Pm) /= N ) then
    call error('wake_model_estimator_t.advance','Pm must be of size N')
end if
if (size(Ctp) /= N ) then
    call error('wake_model_estimator_t.advance','Ctp must be of size N')
end if

! Calculate noisy measurements
this%E = 0._rprec
do i = 1, N
    do j = 1, this%Ne
        this%E(i, j) = this%sigma_Phat * random_normal()
        this%D(i, j) = Pm(i) + this%E(i, j)
    end do
end do
this%Dprime = this%D - this%Ahat

! Update Anew = A + A'*Ahat'^T * (Ahat'*Ahat'^T + E*E^T)^-1 * D'
! Since the dimension is small, we don't bother doing the SVD. If the matrix becomes
! singular, then this should be considered as in section 4.3.2 of Everson(2003)
this%A = this%A + matmul( matmul(this%Aprime, transpose(this%Ahatprime)),            &
matmul(inverse(matmul(this%Ahatprime, transpose(this%Ahatprime)) +                   &
matmul(this%E, transpose(this%E))), this%Dprime))

! Compute mean
this%Abar = 0._rprec
do i = 1, this%Ne
    this%Abar = this%Abar + this%A(:,i) / this%Ne;
end do

! Filter U_infty
alpha = dt / (this%tau_U_infty + dt)
r1_if = this%wm%Ctp(1) / ( 4._rprec + this%wm%Ctp(1) )
Uinftyi = ( Pm(1) / this%wm%Ctp(1) )**(1._rprec/3._rprec) / (1._rprec - r1_if)
this%wm%U_infty = alpha * Uinftyi + (1 - alpha) * this%wm%U_infty
this%wm%VELOCITY = this%wm%U_infty
this%wm%TIME  = this%wm%LENGTH / this%wm%VELOCITY
this%wm%FORCE = this%wm%VELOCITY / this%wm%TIME

! Fill into objects
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend   = j*Nx
        this%ensemble(i)%du(j,:)  = this%A(jstart:jend,i)
    end do
    this%ensemble(i)%U_infty  = this%wm%U_infty
    this%ensemble(i)%VELOCITY = this%wm%VELOCITY
    this%ensemble(i)%TIME  = this%wm%TIME
    this%ensemble(i)%FORCE = this%wm%FORCE
    this%ensemble(i)%k(1:N-1) = this%A(Nx*N+1:,i)
    this%ensemble(i)%k(N)     = this%ensemble(i)%k(N-1)
end do
do j = 1, N
    jstart = (j-1)*Nx+1
    jend   = j*Nx
    this%wm%du(j,:)  = this%Abar(jstart:jend)
end do
this%wm%k(1:N-1) = this%Abar(Nx*N+1:)
this%wm%k(N)     = this%wm%k(N-1)

! Advance ensemble and mean estimate
call this%advanceEnsemble(Ctp, dt)

! Place ensemble into a matrix with each member in a column
this%Abar    = 0
this%Ahatbar = 0
N  = this%wm%N
Nx = this%wm%Nx
do i = 1, this%Ne
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend   = j*Nx
        this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
    end do
    this%A(Nx*N+1:,i) = this%ensemble(i)%k(1:N-1)
    this%Ahat(:,i)    = this%ensemble(i)%Phat
    this%Abar         = this%Abar + this%A(:,i) / this%Ne
    this%Ahatbar      = this%Ahatbar + this%Ahat(:,i) / this%Ne
end do
do j = 1, this%Ne
    this%Aprime(:,j) = this%A(:,j) - this%Abar
    this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
end do

end subroutine

!*******************************************************************************
subroutine advanceEnsemble(this, Ctp, dt)
!*******************************************************************************
use param, only : pi
use util, only  : random_normal
implicit none
class(wake_model_estimator_t), intent(inout) :: this
real(rprec), dimension(:), intent(in)    :: Ctp
real(rprec), intent(in)                  :: dt
integer                                  :: i, j, N

! Advance step for objects
N = this%wm%N
do i = 1, this%Ne
    do j = 1,N
        this%ensemble(i)%k(j) = this%ensemble(i)%k(j)                                &
                              + sqrt(dt) * this%sigma_k * random_normal()
        this%ensemble(i)%du(j,:) = this%ensemble(i)%du(j,:)                          &
            + sqrt(dt) * this%sigma_du * random_normal()                             &
            * this%ensemble(i)%G(j,:) * sqrt(2*pi) * this%ensemble(i)%Delta
    end do
    this%ensemble(i)%k(N) = this%ensemble(i)%k(N-1)
    call this%ensemble(i)%computeWakeExpansionFunctions
    call this%ensemble(i)%advance(Ctp, dt)
end do
call this%wm%computeWakeExpansionFunctions
call this%wm%advance(Ctp, dt)

end subroutine advanceEnsemble

end module wake_model_estimator
