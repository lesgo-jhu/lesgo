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
    real(rprec) :: sqrt6overdelta, t_half
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
real(rprec) :: R1, R2, Rval

R2 = linear_interp(this%r, this%R2, r)
R1 = 0.5_rprec * ( erf(this%sqrt6overdelta*(x + this%t_half)) -                &
    erf(this%sqrt6overdelta*(x - this%t_half)) )
Rval = R1 * R2

end function val

!*******************************************************************************
subroutine init(this, delta2, thk, dia)
!*******************************************************************************
use param, only : write_endian, path, pi
use functions, only : bilinear_interp
implicit none
include'fftw3.f'

class(turb_ind_func_t), intent(inout) :: this
real(rprec), intent(in) :: delta2, thk, dia

real(rprec) :: L, d, R
integer :: N
integer, dimension(:), allocatable :: ind
real(rprec), dimension(:), allocatable :: yz
real(rprec), dimension(:,:), allocatable :: g, f, h
real(rprec), dimension(:), allocatable :: xi
real(rprec) :: dr, Lr
integer :: i, j

integer*8 plan
complex(rprec), dimension(:,:), allocatable :: ghat, fhat, hhat

L = 3._rprec * dia
N = 4*ceiling(2._rprec*L / sqrt(delta2))
d = L / N
R = 0.5 * dia

allocate(yz(N))
allocate(ind(N))
allocate(g(N, N))
allocate(h(N, N))
allocate(f(N, N))
allocate(ghat(N/2+1, N))
allocate(hhat(N/2+1, N))
allocate(fhat(N/2+1, N))

! Calculate constants
this%t_half = 0.5 * thk
this%sqrt6overdelta = sqrt(6._rprec) / sqrt(delta2)

! Calculate yz and indices to sort the result
do i = 1, N/2
    yz(i) = d*(i-0.5)
    ind(i) = N/2+i
end do
do i = N/2+1, N
    yz(i) = -L + d*(i-0.5)
    ind(i) = i-N/2
end do

! Calculate g and f
do j = 1, N
    do i = 1, N
        g(i,j) = exp(-6*(yz(i)**2+yz(j)**2)/delta2)
        if (sqrt(yz(i)**2 + yz(j)**2) < R) then
            h(i,j) = 1.0
        else
            h(i,j) = 0.0
        end if
    end do
end do

! Do the convolution f = g*h in fourier space
call dfftw_plan_dft_r2c_2d(plan, N, N, g, ghat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, g, ghat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, h, hhat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, h, hhat)
call dfftw_destroy_plan(plan)

fhat = ghat*hhat

! Compute the inverse fft of fhat
call dfftw_plan_dft_c2r_2d(plan, N, N, fhat, f, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fhat, f)
call dfftw_destroy_plan(plan)

! Normalize
f = f / N**2 * d**2

! Sort the results
f = f(ind,ind)
yz = yz(ind);

! Interpolate onto the lookup table
allocate(xi(N))
if (allocated(this%r) ) then
    deallocate(this%r)
end if
allocate( this%r(N) )
allocate( this%R2(N) )

Lr = R + 4 * sqrt(delta2)
dr = Lr / (N - 1)
do i = 1,N
    this%r(i) = (i-1)*dr
    xi(i) = 0
end do
this%R2 = max(6._rprec / pi / delta2 * bilinear_interp(yz, yz, f, xi, this%r), &
    0._rprec)

end subroutine init

end module turbine_indicator
