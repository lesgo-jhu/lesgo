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
module wake_model
!*******************************************************************************
use types, only : rprec
use util,  only : logistic, softplus, gaussian
use wake_model_base
use messages
use bicubic_spline
use param, only : pi
implicit none

private
public wake_model_t

type, extends(wake_model_base_t) :: wake_model_t
    ! velocity deficit (turbine, space)
    real(rprec), dimension(:,:), allocatable :: du
    ! superimposed velocity (space)
    real(rprec), dimension(:), allocatable :: u
    ! estimated local turbine velocity (turbine)
    real(rprec), dimension(:), allocatable :: uhat
    ! estimated generator power (turbine)
    real(rprec), dimension(:), allocatable :: Phat
    ! estimated aerodynamic power (turbine)
    real(rprec), dimension(:), allocatable :: Paero
    ! local thrust coefficient (turbine)
    real(rprec), dimension(:), allocatable :: Ctp
    ! local thrust coefficient (turbine)
    real(rprec), dimension(:), allocatable :: Cpp
    ! rotational speed (turbine)
    real(rprec), dimension(:), allocatable :: omega
    ! blade pitch angle
    real(rprec), dimension(:), allocatable :: beta
    ! local tip speed ratio
    real(rprec), dimension(:), allocatable :: lambda_prime
    ! generator torque
    real(rprec), dimension(:), allocatable :: gen_torque
contains
    procedure, public :: initialize_val
!     procedure, private :: initialize_file
!     procedure, public :: print
!     procedure, public :: write_to_file
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: advance
    procedure, private :: rhs
    procedure, public :: adjoint_values
end type wake_model_t

interface wake_model_t
    module procedure :: constructor_val
!     module procedure :: constructor_file
end interface wake_model_t

contains

!*******************************************************************************
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,           &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_t) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in) :: i_Nx
type(bicubic_spline_t), intent(in) :: i_Ctp_spline, i_Cpp_spline

call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,           &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline)
end function constructor_val
!
! !*******************************************************************************
! function constructor_file(fstring) result(this)
! !*******************************************************************************
! ! Constructor for wake model that reads from file
! use open_file_fid_mod
! use param, only : CHAR_BUFF_LENGTH
! implicit none
!
! type(wake_model_t) :: this
! character(*), intent(in) :: fstring
!
! call this%initialize_file(fstring)
!
! end function constructor_file

!*******************************************************************************
subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,    &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline)
!*******************************************************************************
implicit none
class(wake_model_t), intent(inout) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in) :: i_Nx
type(bicubic_spline_t), intent(in) :: i_Ctp_spline, i_Cpp_spline

! Call base class initializer
call this%wake_model_base_t%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia,&
    i_rho, i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline)

allocate( this%du(this%N, this%Nx) )
allocate( this%u(this%Nx) )
allocate( this%uhat(this%N) )
allocate( this%Phat(this%N) )
allocate( this%Paero(this%N) )
allocate( this%Ctp(this%N)  )
allocate( this%Cpp(this%N)  )
allocate( this%omega(this%N)  )
allocate( this%beta(this%N)  )
allocate( this%lambda_prime(this%N)  )
allocate( this%gen_torque(this%N)  )

this%du = 0._rprec
this%u = this%U_infty
this%uhat = this%U_infty
this%Phat = 0._rprec
this%Paero = 0._rprec
this%Ctp = 0._rprec
this%Cpp = 0._rprec
this%omega = 1._rprec
this%beta = 0._rprec
this%lambda_prime = 0.5_rprec * this%Dia / this%U_infty
this%gen_torque = 0._rprec

end subroutine initialize_val
!
! !*******************************************************************************
! subroutine initialize_file(this, fstring)
! !*******************************************************************************
! use open_file_fid_mod
! use param, only : CHAR_BUFF_LENGTH
! implicit none
!
! class(wake_model_t), intent(inout) :: this
! character(*), intent(in) :: fstring
! integer :: i, fid
!
! !  Open vel.out (lun_default in io) for final output
! fid = open_file_fid(fstring, 'rewind', 'unformatted')
! read(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                      &
!           this%U_infty, this%isDimensionless
! read(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
!
! allocate( this%s(this%N)    )
! allocate( this%k(this%N)    )
! allocate( this%uhat(this%N) )
! allocate( this%Phat(this%N) )
! allocate( this%Ctp(this%N)  )
! allocate( this%x(this%Nx) )
! allocate( this%u(this%Nx) )
! allocate( this%G(this%N,  this%Nx) )
! allocate( this%d(this%N,  this%Nx) )
! allocate( this%dp(this%N, this%Nx) )
! allocate( this%w(this%N,  this%Nx) )
! allocate( this%f(this%N, this%Nx) )
! allocate( this%du(this%N, this%Nx) )
!
! read(fid) this%s
! read(fid) this%k
! read(fid) this%x
! read(fid) this%du
! read(fid) this%u
! read(fid) this%uhat
! read(fid) this%Phat
! read(fid) this%Ctp
! close(fid)
!
! do i = 1, this%N
!     this%G(i,:) = gaussian(this%x, this%s(i), this%Delta)
!     this%G(i,:) = this%G(i,:) / sum(this%G(i,:)) / this%dx
! end do
! call this%computeWakeExpansionFunctions
!
! end subroutine initialize_file

!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(wake_model_t), intent(inout) :: this

if (.not.this%isDimensionless) then
    call this%wake_model_base_t%makeDimensionless
    ! units V
    this%du = this%du / this%VELOCITY
    this%u = this%u / this%VELOCITY
    this%uhat = this%uhat / this%VELOCITY
    ! units power
    this%Phat = this%Phat / this%POWER
    this%Paero = this%Paero / this%POWER
    ! units T^-1
    this%omega = this%omega * this%TIME
    ! units TORQUE
    this%gen_torque = this%gen_torque / this%TORQUE
end if

end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(wake_model_t), intent(inout) :: this

if (this%isDimensionless) then
    call this%wake_model_base_t%makeDimensional
    ! units V
    this%du = this%du * this%VELOCITY
    this%u = this%u * this%VELOCITY
    this%uhat = this%uhat * this%VELOCITY
    ! units power
    this%Phat = this%Phat * this%POWER
    this%Paero = this%Paero * this%POWER
    ! units T^-1
    this%omega = this%omega / this%TIME
    ! units TORQUE
    this%gen_torque = this%gen_torque * this%TORQUE
end if

end subroutine makeDimensional
!
! !*******************************************************************************
! subroutine write_to_file(this, fstring)
! !*******************************************************************************
! ! Writes object to file
! use open_file_fid_mod
! use param, only : CHAR_BUFF_LENGTH
! implicit none
! class(wake_model_t), intent(in) :: this
! character(CHAR_BUFF_LENGTH), intent(in) :: fstring
! integer :: fid
! !
! ! !  Open vel.out (lun_default in io) for final output
! ! fid = open_file_fid(fstring, 'rewind', 'unformatted')
! ! write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                     &
! !            this%U_infty, this%isDimensionless
! ! write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
! ! write(fid) this%s
! ! write(fid) this%k
! ! write(fid) this%x
! ! write(fid) this%du
! ! write(fid) this%u
! ! write(fid) this%uhat
! ! write(fid) this%Phat
! ! write(fid) this%Ctp
! ! close(fid)
!
! end subroutine write_to_file
!
! !*******************************************************************************
! subroutine print(this)
! !*******************************************************************************
! ! Prints all variables of the class to standard output
! implicit none
! class(wake_model_t), intent(in) :: this
! integer :: i
!
! write(*,*) ' U_infty = ', this%U_infty
! write(*,*) ' Delta   = ', this%Delta
! write(*,*) ' Dia     = ', this%Dia
! write(*,*) ' Nx      = ', this%Nx
! write(*,*) ' x       = ', this%x
! write(*,*) ' dx      = ', this%dx
! write(*,*) ' isDimensionless = ', this%isDimensionless
! write(*,*) ' LENGTH  = ', this%LENGTH
! write(*,*) ' VELOCITY = ', this%VELOCITY
! write(*,*) ' TIME    = ', this%TIME
! ! write(*,*) ' FORCE   = ', this%FORCE
! write(*,*) ' u       = ', this%u
! do i = 1, this%N
!     write(*,*) ' Wake', i,':'
!     write(*,*) '  uhat = ', this%uhat(i)
!     write(*,*) '  s = ', this%s(i)
!     write(*,*) '  k = ', this%k(i)
!     write(*,*) '  G = ', this%G(i,:)
!     write(*,*) '  d = ', this%d(i,:)
!     write(*,*) '  dp = ', this%dp(i,:)
!     write(*,*) '  w = ', this%w(i,:)
! end do
!
! end subroutine print

!*******************************************************************************
subroutine advance(this, beta, gen_torque, dt)
!*******************************************************************************
! Note: every input value is for time step n. Values at time step n-1 were saved
! during the previous call the advance and are used before being reassigned.
implicit none
class(wake_model_t), intent(inout) :: this
real(rprec), intent(in) :: dt
real(rprec), dimension(:), intent(in) :: beta, gen_torque
integer :: i
real(rprec), dimension(:), allocatable :: du_superimposed

if (size(beta) /= this%N .or. size(gen_torque) /= this%N) then
    call error('wake_model_t.advance','beta and gen_torque must be size N')
end if

! Compute new wake deficit and superimpose wakes
allocate(du_superimposed(this%Nx))
du_superimposed = 0._rprec
do i = 1, this%N
    this%du(i,:) = this%du(i,:) +  dt * this%rhs(this%du(i,:),                 &
        this%f(i,:) * this%Ctp(i) / (4._rprec + this%Ctp(i)), i)
    du_superimposed = du_superimposed + this%du(i,:)**2
end do
du_superimposed = sqrt(du_superimposed)

! Calculate new rotational speed
do i = 1, this%N
    this%omega(i) = this%omega(i) + dt * (this%Paero(i) / this%omega(i)       &
        - this%gen_torque(i)) / this%inertia
end do

! Find the velocity field
this%u = this%U_infty - du_superimposed

! Find estimated velocities, coefficients, and power
do i = 1, this%N
    this%beta(i) = beta(i)
    this%gen_torque(i) = gen_torque(i)
    this%uhat(i) = sum(this%G(i,:) * this%u * this%dx)
    this%lambda_prime(i) = 0.5_rprec * this%omega(i) * this%Dia / this%uhat(i)
    call this%Ctp_spline%interp(this%beta(i), this%lambda_prime(i), this%Ctp(i))
    call this%Cpp_spline%interp(this%beta(i), this%lambda_prime(i), this%Cpp(i))
    this%Paero(i) = this%rho * pi * this%Dia**2 * this%Cpp(i)                  &
                    * this%uhat(i)**3 / 8._rprec
    this%Phat(i) = this%gen_torque(i) * this%omega(i)
end do

end subroutine advance

!*******************************************************************************
function rhs(this, du, f, i) result(ddudt)
!*******************************************************************************
! Evaluates RHS of wake equation
use util, only : ddx_upwind1
implicit none
class(wake_model_t), intent(in) :: this
real(rprec), dimension(:), intent(in) :: f, du
integer, intent(in) :: i
real(rprec), dimension(:), allocatable :: ddudt, ddudx

allocate(ddudt(this%Nx))
allocate(ddudx(this%Nx))

ddudx = ddx_upwind1(du, this%dx)
ddudt = -this%U_infty * ddudx - this%w(i,:) * du + f

end function rhs

!*******************************************************************************
subroutine adjoint_values(this, Pref, fstar, Uw, Udu, Wj, Ww, Wdu, Bw, Bdu)
!*******************************************************************************

implicit none
class(wake_model_t), intent(in) :: this
real(rprec), intent(in) :: Pref
real(rprec), dimension(:,:), intent(out) :: fstar
real(rprec), dimension(:), intent(out) :: Uw, Udu, Wj, Ww, Wdu, Bw, Bdu

real(rprec), dimension(:), allocatable :: du_super, dCt_dbeta, dCt_dlambda
real(rprec), dimension(:), allocatable :: dCp_dbeta, dCp_dlambda
real(rprec) :: dummy
integer :: n, i
allocate(du_super(this%Nx))
allocate(dCt_dbeta(this%N))
allocate(dCt_dlambda(this%N))
allocate(dCp_dbeta(this%N))
allocate(dCp_dlambda(this%N))

! Calculate fstar and derivatives of Ct and Cp
du_super = this%U_infty - this%u
do n = 1, this%N
    fstar(n, :) = - this%du(n,:) / du_super
    do i = 1, this%Nx
        if ( du_super(i) <= 1E-10 )   fstar(n, i) = 0._rprec
    end do
    call this%Ctp_spline%interp(this%beta(n), this%lambda_prime(n), dummy,     &
                           dCt_dbeta(n), dCt_dlambda(n))
    call this%Cpp_spline%interp(this%beta(n), this%lambda_prime(n), dummy,     &
                           dCp_dbeta(n), dCp_dlambda(n))
end do

! Everything else can be calculated at once
Uw = 3._rprec * this%Paero / this%uhat / this%omega / this%inertia             &
    - this%Paero * dCp_dlambda * 0.5_rprec*this%Dia / this%uhat**2 /this%Cpp / this%inertia
Udu = - 4._rprec / (4._rprec + this%Ctp)**2 * dCt_dlambda * this%omega         &
    * 0.5_rprec * this%Dia / this%uhat**2
Wj = - 2._rprec * (sum(this%Phat) - Pref) * this%Phat / this%omega
Ww = -this%Paero / this%omega**2 / this%inertia                                &
    + this%Paero * dCp_dlambda * 0.5_rprec*this%Dia / this%uhat / this%omega /this%Cpp / this%inertia
Wdu = 4._rprec / (4._rprec + this%Ctp)**2 * dCt_dlambda * 0.5_rprec * this%Dia / this%uhat
Bw = - dCp_dbeta * this%Paero / this%omega / this%Cpp / this%inertia
Bdu = -8._rprec * this%U_infty**2 / (4._rprec + this%Ctp)**2 * dCt_dbeta
do n = 1, this%N
    if (this%Cpp(n) <= 1E-10) then
        Uw = 0._rprec
        Ww = 0._rprec
        Bw = 0._rprec
    end if
end do

deallocate(du_super)
deallocate(dCt_dbeta)
deallocate(dCt_dlambda)
deallocate(dCp_dbeta)
deallocate(dCp_dlambda)

end subroutine adjoint_values

end module wake_model
