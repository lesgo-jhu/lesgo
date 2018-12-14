!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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
module turbine_indicator
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lh

private
public :: turb_ind_func_t

! Indicator function calculator
type turb_ind_func_t
    real(rprec), dimension(:), allocatable :: r
    real(rprec), dimension(:), allocatable :: R2
    real(rprec) :: sqrt6overdelta1, ell, delta1, delta2, thk
    real(rprec) :: M
contains
    procedure, public :: init
    procedure, public :: val
end type turb_ind_func_t

contains

!*******************************************************************************
function val(this, r, x) result(Rval)
!*******************************************************************************
use functions, only : linear_interp
implicit none
class(turb_ind_func_t), intent(in) :: this
real(rprec), intent(in) :: r, x
real(rprec) :: R1, R2, Rval, rr, xx

! Make dimensionless to lookup value
rr = r / this%ell
xx = x / this%ell

! Calculate as a dimensionless value, then return as dimensional
R2 = linear_interp(this%r, this%R2, rr)
R1 = 0.5_rprec / this%thk * ( erf(this%sqrt6overdelta1*(xx + 0.5*this%thk))    &
    - erf(this%sqrt6overdelta1*(xx - 0.5*this%thk)) )
Rval = R1 * R2 / this%ell**3

end function val

!*******************************************************************************
subroutine init(this, delta1, delta2, thk, dia)
!*******************************************************************************
use param, only : write_endian, path, pi
use functions, only : bilinear_interp
implicit none
include'fftw3.f'

class(turb_ind_func_t), intent(inout) :: this
real(rprec), intent(in) :: delta1, delta2, thk, dia

real(rprec) :: L, d, dr
integer :: i, j, N
real(rprec), dimension(:), allocatable :: yz
real(rprec), dimension(:,:), allocatable :: g, f, h

integer*8 plan
complex(rprec), dimension(:,:), allocatable :: ghat, fhat, hhat

! Units for all values
this%ell = 0.5_rprec * dia
this%delta1 = delta1 / this%ell
this%delta2 = delta2 / this%ell
this%thk = thk / this%ell
this%sqrt6overdelta1 = sqrt(6._rprec) / this%delta1

! Size of radial domain
L = 1 + 10*this%delta2/sqrt(12._rprec)
N = 2 * 1024
d = 2 * L / N

! Create yz grid
allocate(yz(N))
do i = 1, N
    yz(i) = -L + d * (i - 0.5_rprec)
end do

! Create Circle and Gaussian
allocate(h(N, N))
allocate(g(N, N))
do i = 1, N
    do j = 1, N
        g(i,j) = 6._rprec/(pi*this%delta2**2)                                  &
            * exp(-6*(yz(i)**2+yz(j)**2)/this%delta2**2)
        if (sqrt(yz(i)**2 + yz(j)**2) < 1) then
            h(i,j) = 1.0_rprec / pi
        else
            h(i,j) = 0.0
        end if
    end do
end do

! Do the convolution f = g*h in fourier space
allocate(ghat(N/2+1, N))
allocate(hhat(N/2+1, N))
allocate(fhat(N/2+1, N))
allocate(f(N, N))

call dfftw_plan_dft_r2c_2d(plan, N, N, g, ghat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, g, ghat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, h, hhat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, h, hhat)
call dfftw_destroy_plan(plan)

fhat = ghat * hhat * d**2!* d**2!*d**2!* d**2! * (d)**2

call dfftw_plan_dft_c2r_2d(plan, N, N, fhat, f, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fhat, f)
call dfftw_destroy_plan(plan)

! Place into the lookup table
allocate(this%r(N/2))
allocate(this%R2(N/2))
this%r = yz(1:N/2) + L
this%R2 = f(1,1:N/2) / N**2

! Normalize to integrate to unity
dr = this%r(2) - this%r(1)
this%R2 = this%R2/sum(2*pi*this%r*this%R2*dr)

! Calculate correction factor
this%M = sum(2*pi*this%r*this%R2**2*dr)*pi

end subroutine init

end module turbine_indicator
