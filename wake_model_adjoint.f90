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
module wake_model_adjoint
!*******************************************************************************
use types, only : rprec
use util,  only : logistic, softplus, gaussian
use wake_model_base
use messages
use bicubic_spline
use param, only : pi
implicit none

private
public wake_model_adjoint_t

type, extends(wake_model_base_t) :: wake_model_adjoint_t
    ! adjoint velocity deficit (turbine,space)
    real(rprec), dimension(:,:), allocatable :: du_star
    ! adjoint superimposed velocity (space)
    real(rprec), dimension(:), allocatable :: u_star
    ! adjoint estimated local turbine velocity (turbine)
    real(rprec), dimension(:), allocatable :: uhat_star
    ! adjoint rotational speed (turbine)
    real(rprec), dimension(:), allocatable :: omega_star
contains
    procedure, public :: initialize_val
    ! procedure, private :: initialize_file
!     procedure, public :: print
!     procedure, public :: write_to_file
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: retract
    procedure, private :: rhs
end type wake_model_adjoint_t

interface wake_model_adjoint_t
    module procedure :: constructor_val
    ! module procedure :: constructor_file
end interface wake_model_adjoint_t

contains

!*******************************************************************************
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,           &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_adjoint_t) :: this
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
! type(wake_model_adjoint_t) :: this
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
class(wake_model_adjoint_t), intent(inout) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in) :: i_Nx
type(bicubic_spline_t), intent(in) :: i_Ctp_spline, i_Cpp_spline

! Call base class initializer
call this%wake_model_base_t%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia,&
    i_rho, i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline)

allocate( this%du_star(this%N, this%Nx) )
allocate( this%u_star(this%Nx) )
allocate( this%uhat_star(this%N) )
allocate( this%omega_star(this%N)  )
!
this%du_star = 0._rprec
this%u_star = 0._rprec
this%uhat_star = 0._rprec
this%omega_star = 0._rprec

end subroutine initialize_val
!
! !*******************************************************************************
! subroutine initialize_file(this, fstring)
! !*******************************************************************************
! use open_file_fid_mod
! use param, only : CHAR_BUFF_LENGTH
! implicit none
!
! class(wake_model_adjoint_t), intent(inout) :: this
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
! allocate( this%fp(this%N, this%Nx) )
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
!
!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(wake_model_adjoint_t), intent(inout) :: this

if (.not.this%isDimensionless) then
    call this%wake_model_base_t%makeDimensionless
    ! units V
    this%du_star = this%du_star / this%VELOCITY
    this%u_star = this%u_star / this%VELOCITY
    this%uhat_star = this%uhat_star / this%VELOCITY
    ! units T^(-1)
    this%omega_star = this%omega_star * this%TIME
end if

end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(wake_model_adjoint_t), intent(inout) :: this

if (this%isDimensionless) then
    call this%wake_model_base_t%makeDimensional
    ! units V
    this%du_star = this%du_star * this%VELOCITY
    this%u_star = this%u_star * this%VELOCITY
    this%uhat_star = this%uhat_star * this%VELOCITY
    ! units T^(-1)
    this%omega_star = this%omega_star / this%TIME
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
! class(wake_model_adjoint_t), intent(in) :: this
! character(CHAR_BUFF_LENGTH), intent(in) :: fstring
! integer :: fid
!
! !  Open vel.out (lun_default in io) for final output
! fid = open_file_fid(fstring, 'rewind', 'unformatted')
! write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                     &
!            this%U_infty, this%isDimensionless
! write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
! write(fid) this%s
! write(fid) this%k
! write(fid) this%x
! write(fid) this%du
! write(fid) this%u
! write(fid) this%uhat
! write(fid) this%Phat
! write(fid) this%Ctp
! close(fid)
!
! end subroutine write_to_file
!
! !*******************************************************************************
! subroutine print(this)
! !*******************************************************************************
! ! Prints all variables of the class to standard output
! implicit none
! class(wake_model_adjoint_t), intent(in) :: this
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
! write(*,*) ' FORCE   = ', this%FORCE
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
!
! !*******************************************************************************
! subroutine advance(this, beta, gen_torque, dt)
! !*******************************************************************************
! implicit none
! class(wake_model_adjoint_t), intent(inout) :: this
! real(rprec), intent(in) :: dt
! real(rprec), dimension(:), intent(in) :: beta, gen_torque
! integer :: i
! real(rprec), dimension(:), allocatable :: du_superimposed
! real(rprec) :: lambda_prime, Paero
!
! if (size(beta) /= this%N .or. size(gen_torque) /= this%N) then
!     call error('wake_model_adjoint_t.advance','beta and gen_torque must be size N')
! end if
!
! ! Compute new wake deficit and superimpose wakes
! allocate(du_superimposed(this%Nx))
! du_superimposed = 0.0
! do i = 1, this%N
!     this%du(i,:) = this%du(i,:) +  dt * this%rhs(this%du(i,:),                 &
!         this%fp(i,:) * this%Ctp(i) / (4.0 + this%Ctp(i)), i)
!     du_superimposed = du_superimposed + this%du(i,:)**2
! end do
! du_superimposed = sqrt(du_superimposed)
!
! ! Calculate new rotational speed
! do i = 1, this%N
!     Paero = this%rho * pi * this%Dia**2 * this%Cpp(i) * this%uhat(i)**3 / 8.d0
!     this%omega(i) = this%omega(i) + dt * (Paero / this%omega(i)         &
!         - this%gen_torque(i)) / this%inertia
!
! end do
!
! ! Find the velocity field
! this%u = this%U_infty - du_superimposed
!
! ! Find estimated velocities
! do i = 1, this%N
!     this%beta(i) = beta(i)
!     this%gen_torque(i) = gen_torque(i)
!     this%uhat(i) = sum(this%G(i,:) * this%u * this%dx)
!     lambda_prime = 0.5_rprec * this%omega(i) * this%Dia / this%uhat(i)
!     call this%Ctp_spline%interp(this%beta(i), lambda_prime, this%Ctp(i))
!     call this%Cpp_spline%interp(this%beta(i), lambda_prime, this%Cpp(i))
!     this%Phat(i) = this%gen_torque(i) * this%omega(i)
! end do
!
! end subroutine advance

!*******************************************************************************
subroutine retract(this, fstar, Adu, Aw, Bj, Bdu, Bw, dt)!, g)
    !***************************************************************************
    implicit none
    class(wake_model_adjoint_t), intent(inout) :: this
    real(rprec), dimension(:,:), intent(in) :: fstar
    real(rprec), dimension(:), intent(in) :: Adu, Aw, Bj, Bdu, Bw
    real(rprec), intent(in) :: dt
    real(rprec) :: fdustar
    integer :: i

    ! adjoint of velocity field is a sum. Set to 0
    this%u_star = 0._rprec

    ! compute adjoints of everything except velocity field
    do i = 1, this%N
        ! evaluate \int_0^L f_n(x) (du_star)_n(x,t) \, dx
        fdustar = sum(this%f(i,:)*this%du_star(i,:)) * this%dx
        ! forward differencing of adjoint rotational speed equation
        this%omega_star(i) = this%omega_star(i) + dt * ( Bj(i) + Bdu(i)*fdustar&
                             + Bw(i) * this%omega_star(i) )
        ! adjoint of estimated velocity
        this%uhat_star(i) = Adu(i)*fdustar + Aw(i)*this%omega_star(i)
        ! superimpose adjoints of estimated velocities
        this%u_star = this%u_star + this%G(i,:)*this%uhat_star(i)
    enddo

    ! compute velocity field adjoints
    do i = 1, this%N
        this%du_star(i,:) = this%du_star(i,:) - dt                             &
            * this%rhs(this%du_star(i,:), fstar(i,:)*this%u_star, i)
    enddo

end subroutine retract

!*******************************************************************************
function rhs(this, du, f, i) result(ddudt)
!*******************************************************************************
! Evaluates RHS of wake equation
use util, only : ddx_downwind1
implicit none
class(wake_model_adjoint_t), intent(in) :: this
real(rprec), dimension(:), intent(in) :: f, du
integer, intent(in) :: i
real(rprec), dimension(:), allocatable :: ddudt, ddudx

allocate(ddudt(this%Nx))
allocate(ddudx(this%Nx))

ddudx = ddx_downwind1(du, this%dx)
ddudt = -this%U_infty * ddudx - this%w(i,:) * du - f

end function rhs

end module wake_model_adjoint
