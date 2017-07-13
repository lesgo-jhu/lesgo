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
module wake_model_base
!*******************************************************************************
use types, only : rprec
use util, only : logistic, softplus, gaussian
use messages
use bi_pchip
use param, only : pi
implicit none

private
public wake_model_base_t

type wake_model_base_t
    ! streamwise turbine location
    real(rprec), dimension(:), allocatable :: sx
    ! spanwise turbine location
    real(rprec), dimension(:), allocatable :: sy
    ! wake expansion coefficient
    real(rprec), dimension(:), allocatable :: k
    ! streamwise coordinate
    real(rprec), dimension(:), allocatable :: x
    ! spanwise coordinate
    real(rprec), dimension(:), allocatable :: y
    ! Gaussian forcing function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: G
    integer, dimension(:), allocatable :: Gstart, Gend
    ! Indicator functions for turbine wakes (turbine, space)
    integer, dimension(:,:), allocatable :: Istart, Iend
    real(rprec), dimension(:,:), allocatable :: Isum
    ! dimensionless wake diameter (turbine, space)
    real(rprec), dimension(:,:), allocatable :: d
    ! d/dx of d (turbine, space)
    real(rprec), dimension(:,:), allocatable :: dp
    ! wake expansion function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: w
    ! forcing prefactor (turbine, space)
    real(rprec), dimension(:,:), allocatable :: f
    ! inlet velocity
    real(rprec) :: U_infty = 0
    ! Gaussian forcing width
    real(rprec) :: Delta = 0
    ! rotor diameter
    real(rprec) :: Dia = 0
    ! air density
    real(rprec) :: rho = 0
    ! Rotor inertia
    real(rprec) :: inertia = 0
    ! grid spacing
    real(rprec) :: dx = 0
    real(rprec) :: dy = 0
    ! Number of points
    integer :: Nx = 0
    integer :: Ny = 0
    ! Number of turbines
    integer :: N = 0
    ! Specifies whether the wake model is in a dimensionless state
    logical :: isDimensionless = .false.
    ! Dimensional scales
    real(rprec) :: LENGTH = 0, VELOCITY = 0, TIME = 0, MASS = 0, TORQUE = 0, POWER = 0
    ! Splines for Ctp and Cpp
    type(bi_pchip_t) Ctp_spline, Cpp_spline
    ! Determine whether a turbine is waked or not
    logical, dimension(:), allocatable :: waked
    integer :: Nwaked
contains
    procedure, public :: initialize_val
    ! procedure, public :: print
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: compute_wake_expansion
end type wake_model_base_t

interface wake_model_base_t
    module procedure :: constructor_val
end interface wake_model_base_t

contains

!*******************************************************************************
function constructor_val(i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia, i_rho,    &
    i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_base_t) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_sx, i_sy, i_k
integer, intent(in) :: i_Nx, i_Ny
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline

call this%initialize_val(i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia, i_rho,    &
    i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline)

end function constructor_val

!*******************************************************************************
subroutine initialize_val(this, i_sx, i_sy, i_U_infty, i_Delta, i_k, i_Dia,    &
    i_rho,i_inertia, i_Nx, i_Ny, i_Ctp_spline, i_Cpp_spline)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_sx, i_sy, i_k
integer, intent(in) :: i_Nx, i_Ny
integer :: i, j, ii
type(bi_pchip_t), intent(in) :: i_Ctp_spline, i_Cpp_spline
real(rprec) :: xstart, xend, ystart, yend
integer, dimension(1) :: temp_int
real(rprec) :: R

! Allocate based on number of turbines
this%N = size(i_sx)
if ( size(i_k) /= this%N .or. size(i_sy) /= this%N) then
    call error('wake_model_base_t.initialize',                                 &
        'sx, sy, and k must be the same size')
end if
allocate( this%sx(this%N) )
allocate( this%sy(this%N) )
allocate( this%k(this%N) )
allocate( this%Gstart(this%N) )
allocate( this%Gend(this%N) )

! Assign splines
this%Ctp_spline = i_Ctp_spline
this%Cpp_spline = i_Cpp_spline

! Allocate based on number of gridpoints
this%Nx = i_Nx
this%Ny = i_Ny
allocate( this%x(this%Nx) )
allocate( this%y(this%Ny) )
allocate( this%G(this%N,  this%Nx) )
allocate( this%d(this%N,  this%Nx) )
allocate( this%dp(this%N, this%Nx) )
allocate( this%w(this%N,  this%Nx) )
allocate( this%f(this%N, this%Nx) )
allocate( this%Istart(this%N, this%Nx) )
allocate( this%Iend(this%N, this%Nx) )
allocate( this%Isum(this%N, this%Nx) )
allocate( this%waked(this%N) ) 

! Assign input arguments
this%sx = i_sx
this%sy = i_sy
this%U_infty = i_U_infty
this%Delta = i_Delta
this%k = i_k
this%Dia = i_Dia
this%rho = i_rho
this%inertia = i_inertia

! Normalization constants
this%VELOCITY = i_U_infty
this%LENGTH = i_Dia
this%TIME = this%LENGTH / this%VELOCITY
this%MASS = this%rho * this%LENGTH**3
this%TORQUE = this%MASS * this%LENGTH**2 / this%TIME**2
this%POWER = this%MASS * this%LENGTH**2 / this%TIME**3

! Streamwise coordinate
xstart = minval(this%sx) - 2.5_rprec*this%Dia
xend = maxval(this%sx) + 7.5_rprec*this%Dia
this%dx = ( xend-xstart ) / this%Nx
this%x(1) = xstart + this%dx/2
do i = 2, this%Nx
    this%x(i) = this%x(i-1) + this%dx
end do

! Spanwise coordinate
ystart = minval(this%sy) - 5._rprec*this%Dia
yend = maxval(this%sy) + 5._rprec*this%Dia
this%dy = ( yend-ystart ) / this%Ny
this%y(1) = ystart + this%dy/2
do i = 2, this%Ny
    this%y(i) = this%y(i-1) + this%dy
end do

! Gaussian functions
this%G = 0._rprec
do i = 1, this%N
    temp_int = minloc(abs(this%x - this%sx(i) + 4*this%Delta))
    this%Gstart(i) = temp_int(1)
    temp_int = minloc(abs(this%x -  this%sx(i) - 4*this%Delta))
    this%Gend(i) = temp_int(1)
    this%G(i,this%Gstart(i):this%Gend(i)) =                                 &
        gaussian(this%x(this%Gstart(i):this%Gend(i)), this%sx(i), this%Delta)
    this%G(i,:) = this%G(i,:) / sum(this%G(i,:)) / this%dx
end do

call this%compute_wake_expansion

end subroutine initialize_val

!*******************************************************************************
subroutine compute_wake_expansion(this)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this
integer :: i, ii, j
integer, dimension(1) :: temp_int
real(rprec) :: R

! Compute wake expansion functions
do i = 1, this%N
    this%d(i,:) = 1.0 + this%k(i) * softplus(2.0 *                             &
        (this%sx(i) + 2.0 * this%Delta)/this%Dia, 2.0*this%x/this%Dia)
    this%dp(i,:) = 2.0 * this%k(i) * logistic(2.0 *                            &
        (this%sx(i) + 2.0 * this%Delta)/this%Dia, 2.0*this%x/this%Dia)/this%Dia
    this%w(i,:) = 2.0 * this%U_infty * this%dp(i,:) / this%d(i,:)
    this%f(i,:) = 2.0 * this%U_infty**2 * this%G(i,:) / ( this%d(i,:)**2 )
end do

! Indicator functions for wakes
do i = 1, this%N
    do ii = 1, this%Nx
        temp_int = minloc(abs(this%y - this%sy(i) + 0.5_rprec*this%Dia*this%d(i,ii)))
        this%Istart(i,ii) = temp_int(1)
        temp_int = minloc(abs(this%y - this%sy(i) - 0.5_rprec*this%Dia*this%d(i,ii)))
        this%Iend(i,ii) = temp_int(1)
        this%Isum(i,ii) = (this%Iend(i,ii) - this%Istart(i,ii) + 1) * this%dy
    end do
end do

! Determine whether a turbine is waked or not
this%waked = .false.
this%Nwaked = 0
do i = 1, this%N
    do j = 1, this%N
        if (this%waked(i)) cycle
        R = 0.5 * this%Dia * (1._rprec + this%k(j) * softplus(2._rprec *       &
        (this%sx(j) + 2.0 * this%Delta)/this%Dia, 2.0 * this%sx(i)/this%Dia) )
        if (this%sx(i)>this%sx(j) .and. (R-abs(this%sy(i)-this%sy(j)))>0) then
            this%waked(i) = .true.
            this%Nwaked = this%Nwaked + 1
            cycle
        end if
    end do
end do

end subroutine compute_wake_expansion

!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this

if (.not.this%isDimensionless) then
    this%isDimensionless = .true.
    ! units L
    this%sx = this%sx / this%LENGTH
    this%sy = this%sy / this%LENGTH
    this%x = this%x / this%LENGTH
    this%y = this%y / this%LENGTH
    this%Delta = this%Delta / this%LENGTH
    this%Dia = this%Dia / this%LENGTH
    this%dx = this%dx / this%LENGTH
    this%dy = this%dy / this%LENGTH
    this%Isum = this%Isum / this%LENGTH
    ! units L^-1
    this%G = this%G * this%LENGTH
    this%dp = this%dp * this%LENGTH
    ! units T^-1
    this%w = this%w * this%TIME
    ! units V*T^-2
    this%f = this%f / this%VELOCITY * this%TIME
    ! units V
    this%U_infty = this%U_infty / this%VELOCITY
    ! units M*L^-3
    this%rho = this%rho / this%MASS * this%LENGTH**3
    ! units TORQUE*T^2
    this%inertia = this%inertia / ( this%TORQUE * this%TIME**2 )
end if
end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this

if (this%isDimensionless) then
    ! units L
    this%sx = this%sx * this%LENGTH
    this%sy = this%sy * this%LENGTH
    this%x = this%x * this%LENGTH
    this%y = this%y * this%LENGTH
    this%Delta = this%Delta * this%LENGTH
    this%Dia = this%Dia * this%LENGTH
    this%dx = this%dx * this%LENGTH
    this%dy = this%dy * this%LENGTH
    this%Isum = this%Isum * this%LENGTH
    ! units L^-1
    this%G = this%G / this%LENGTH
    ! units T^-1
    this%w = this%w / this%TIME
    ! units V*T^-1
    this%f = this%f * this%VELOCITY / this%TIME
    ! units V
    this%U_infty = this%U_infty * this%VELOCITY
    ! units M*L^-3
    this%rho = this%rho * this%MASS / this%LENGTH**3
    ! units TORQUE*T^2
    this%inertia = this%inertia * this%TORQUE * this%TIME**2
end if
end subroutine makeDimensional

! ********************************************************************************
! subroutine write_to_file(this, fstring)
! ********************************************************************************
! Writes object to file
! use open_file_fid_mod
! use param, only : CHAR_BUFF_LENGTH
! implicit none
! class(wake_model_base_t), intent(in) :: this
! character(CHAR_BUFF_LENGTH), intent(in) :: fstring
! integer :: i, fid
!
! !  Open vel.out (lun_default in io) for final output
! fid = open_file_fid(fstring, 'rewind', 'unformatted')
! write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                     &
!            this%U_infty, this%isDimensionless
! write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
! write(fid) this%s
! write(fid) this%k
! write(fid) this%x
!
!
! close(fid)
!
! end subroutine write_to_file
!
! !*******************************************************************************
! subroutine print(this)
! !*******************************************************************************
! ! Prints all variables of the class to standard output
! implicit none
! class(wake_model_base_t), intent(in) :: this
! integer :: i
!
! ! write(*,*) ' U_infty  = ', this%U_infty
! ! write(*,*) ' Delta    = ', this%Delta
! ! write(*,*) ' Dia      = ', this%Dia
! ! write(*,*) ' Nx       = ', this%Nx
! ! write(*,*) ' x        = ', this%x
! ! write(*,*) ' dx       = ', this%dx
! ! write(*,*) ' isDimensionless = ', this%isDimensionless
! ! write(*,*) ' LENGTH   = ', this%LENGTH
! ! write(*,*) ' VELOCITY = ', this%VELOCITY
! ! write(*,*) ' TIME     = ', this%TIME
! ! ! write(*,*) ' FORCE    = ', this%FORCE
! ! do i = 1, this%N
! !     write(*,*) ' Wake', i,':'
! !     write(*,*) '  s = ', this%sx(i)
! !     write(*,*) '  k = ', this%k(i)
! !     write(*,*) '  G = ', this%G(i,:)
! !     write(*,*) '  d = ', this%d(i,:)
! !     write(*,*) '  dp = ', this%dp(i,:)
! !     write(*,*) '  w = ', this%w(i,:)
! ! end do
! end subroutine print

end module wake_model_base
