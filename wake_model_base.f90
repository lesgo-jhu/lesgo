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
use bicubic_spline
implicit none

private
public wake_model_base_t

type wake_model_base_t
    ! turbine location
    real(rprec), dimension(:), allocatable :: s
    ! wake expansion coefficient
    real(rprec), dimension(:), allocatable :: k
    ! streamwise coordinate
    real(rprec), dimension(:), allocatable :: x
    ! Gaussian forcing function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: G   
    ! dimensionless wake diameter (turbine, space)
    real(rprec), dimension(:,:), allocatable :: d
    ! d/dx of d (turbine, space)
    real(rprec), dimension(:,:), allocatable :: dp
    ! wake expansion function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: w
    ! forcing prefactor f = fp*Ctp/(4+Ctp)  (turbine, space)
    real(rprec), dimension(:,:), allocatable :: fp
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
    ! rotor diameter
    real(rprec) :: dx = 0
    ! Number of streamwise points
    integer :: Nx = 0
    ! Number of turbines
    integer :: N = 0
    ! Specifies whether the wake model is in a dimensionless state
    logical :: isDimensionless = .false.
    ! Dimensional scales
    real(rprec) :: LENGTH = 0, VELOCITY = 0, TIME = 0, FORCE = 0
    ! Splines for Ctp and Cpp
    type(bicubic_spline_t) Ctp_spline, Cpp_spline
contains
    procedure, public :: initialize_val
    procedure, public :: print
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: computeWakeExpansionFunctions
end type wake_model_base_t

interface wake_model_base_t
    module procedure :: constructor_val
end interface wake_model_base_t

contains

!*******************************************************************************
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,           &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_base_t) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in) :: i_Nx
type(bicubic_spline_t), intent(in) :: i_Ctp_spline, i_Cpp_spline

call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,           &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline)

end function constructor_val

!*******************************************************************************
subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_rho,    &
    i_inertia, i_Nx, i_Ctp_spline, i_Cpp_spline)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this
real(rprec), intent(in) :: i_U_infty, i_Delta, i_Dia, i_rho, i_inertia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in) :: i_Nx
integer :: i
type(bicubic_spline_t), intent(in) :: i_Ctp_spline, i_Cpp_spline

! Allocate based on number of turbines
this%N = size(i_s)
if ( size(i_k) /= this%N ) then
    call error('wake_model_base_t.initialize','s and k must be the same size')
end if
allocate( this%s(this%N) )
allocate( this%k(this%N) )

! Assign splines
this%Ctp_spline = i_Ctp_spline
this%Cpp_spline = i_Cpp_spline

! Allocate based on number of gridpoints
this%Nx = i_Nx
allocate( this%x(this%Nx) )
allocate( this%G(this%N,  this%Nx) )
allocate( this%d(this%N,  this%Nx) )
allocate( this%dp(this%N, this%Nx) )
allocate( this%w(this%N,  this%Nx) )
allocate( this%fp(this%N, this%Nx) )

! Assign input arguments
this%s   = i_s
this%U_infty = i_U_infty
this%Delta = i_Delta
this%k   = i_k
this%Dia = i_Dia
this%rho = i_rho
this%inertia = i_inertia

! Normalization constants
this%VELOCITY = i_U_infty
this%LENGTH = i_Dia
this%TIME = this%LENGTH / this%VELOCITY
this%FORCE = this%VELOCITY / this%TIME

! Calculate other variables
this%dx = ( this%s(1) + this%s(this%N) ) / this%Nx
this%x(1) = this%dx/2
do i = 2, this%Nx
    this%x(i) = this%x(i-1) + this%dx
end do
do i = 1, this%N
    this%G(i,:) = gaussian(this%x, this%s(i), this%Delta)
    this%G(i,:) = this%G(i,:) / sum(this%G(i,:)) / this%dx
end do

call this%computeWakeExpansionFunctions
    
end subroutine initialize_val

!*******************************************************************************
subroutine computeWakeExpansionFunctions(this)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this
integer :: i

do i = 1, this%N
    this%d(i,:) = 1.0 + this%k(i) * softplus(2.0 *                            &
        (this%s(i) + 2.0 * this%Delta)/this%Dia, 2.0*this%x/this%Dia)
    this%dp(i,:) = 2.0 * this%k(i) * logistic(2.0 *                            &
        (this%s(i) + 2.0 * this%Delta)/this%Dia, 2.0*this%x/this%Dia)/this%Dia
    this%w(i,:) = 2.0 * this%U_infty * this%dp(i,:) / this%d(i,:)
    this%fp(i,:) = 2.0 * this%U_infty**2 * this%G(i,:) / ( this%d(i,:)**2 )   
end do

end subroutine computeWakeExpansionFunctions

!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this

if (.not.this%isDimensionless) then
    this%isDimensionless = .true.
    this%s = this%s / this%LENGTH
    this%U_infty = this%U_infty / this%VELOCITY
    this%Delta = this%Delta / this%LENGTH
    this%Dia = this%Dia / this%LENGTH
    this%dx = this%dx / this%LENGTH
    this%x = this%x / this%LENGTH
    this%G = this%G * this%LENGTH         ! G has units 1/length
    this%dp = this%dp * this%LENGTH        ! dp has units 1/length
    this%w = this%w * this%TIME           ! w has units 1/time
    this%fp = this%fp / this%FORCE
end if
end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(wake_model_base_t), intent(inout) :: this

if (this%isDimensionless) then
    this%isDimensionless = .false.
    this%s = this%s * this%LENGTH
    this%U_infty = this%U_infty * this%VELOCITY
    this%Delta = this%Delta * this%LENGTH
    this%Dia = this%Dia * this%LENGTH
    this%dx = this%dx * this%LENGTH
    this%x = this%x * this%LENGTH
    this%G = this%G / this%LENGTH         ! G has units 1/length
    this%dp = this%dp / this%LENGTH        ! dp has units 1/length
    this%w = this%w / this%TIME           ! w has units 1/time   
    this%fp = this%fp * this%FORCE  
end if
end subroutine makeDimensional


!*******************************************************************************
! subroutine write_to_file(this, fstring)
!*******************************************************************************
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
! write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,               &
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


!*******************************************************************************
subroutine print(this)
!*******************************************************************************
! Prints all variables of the class to standard output
implicit none
class(wake_model_base_t), intent(in) :: this
integer :: i

write(*,*) ' U_infty  = ', this%U_infty
write(*,*) ' Delta    = ', this%Delta
write(*,*) ' Dia      = ', this%Dia
write(*,*) ' Nx       = ', this%Nx
write(*,*) ' x        = ', this%x
write(*,*) ' dx       = ', this%dx
write(*,*) ' isDimensionless = ', this%isDimensionless
write(*,*) ' LENGTH   = ', this%LENGTH
write(*,*) ' VELOCITY = ', this%VELOCITY
write(*,*) ' TIME     = ', this%TIME
write(*,*) ' FORCE    = ', this%FORCE
do i = 1, this%N
    write(*,*) ' Wake', i,':'
    write(*,*) '  s = ', this%s(i)
    write(*,*) '  k = ', this%k(i)
    write(*,*) '  G = ', this%G(i,:)
    write(*,*) '  d = ', this%d(i,:)
    write(*,*) '  dp = ', this%dp(i,:)
    write(*,*) '  w = ', this%w(i,:)
end do
end subroutine print

end module wake_model_base
