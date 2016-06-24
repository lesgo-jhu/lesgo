!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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

!**********************************************************************
! fftw 3.X version
!**********************************************************************
module fft
use types,only:rprec
use param,only:ld,lh,ny,ld_big, ny2, spectra_calc
use iso_c_binding
implicit none
include 'fftw3.f'
save

public :: padd, unpadd, init_fft

public :: kx, ky, k2, forw_spectra
public :: forw, back, forw_big, back_big

real(rprec), allocatable, dimension(:,:) :: kx, ky, k2
integer*8::forw_spectra
integer*8::forw,back,forw_big,back_big

real (rprec), dimension (:, :), allocatable :: data, data_big

contains

!**********************************************************************
subroutine padd (u_big,u)
!**********************************************************************
! puts arrays into larger, zero-padded arrays 
! automatically zeroes the oddballs
use types,only:rprec
use param,only:ld,ld_big,nx,ny,ny2
implicit none

!  u and u_big are interleaved as complex arrays
real(rprec), dimension(ld,ny), intent(in) :: u
real(rprec), dimension(ld_big,ny2), intent(out) :: u_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

! make sure the big array is zeroed!
u_big(:,:) = 0._rprec

! note: split access in an attempt to maintain locality
u_big(:nx,:ny_h) = u(:nx,:ny_h)

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

u_big(:nx,j_big_s:ny2) = u(:nx,j_s:ny)

return
end subroutine padd

!**********************************************************************
subroutine unpadd(cc,cc_big)
!**********************************************************************
use types,only:rprec
use param,only:ld,nx,ny,ny2,ld_big
implicit none

!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( ld, ny ) :: cc
real(rprec), dimension( ld_big, ny2 ) :: cc_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

cc(:nx,:ny_h) = cc_big(:nx,:ny_h)

! oddballs
cc(ld-1:ld,:) = 0._rprec
cc(:,ny_h+1) = 0._rprec

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2
cc(:nx,j_s:ny) = cc_big(:nx,j_big_s:ny2)

end subroutine unpadd

!**********************************************************************
subroutine init_fft()
!**********************************************************************
use param,only:nx,ny,nx2,ny2
implicit none

! Allocate temporary arrays for creating the FFTW plans
allocate( data(ld, ny) )
allocate( data_big(ld_big, ny2) )

! Create the forward and backward plans for the unpadded and padded
! domains. Notice we are using FFTW_UNALIGNED since the arrays used will not be
! guaranteed to be memory aligned. 
call dfftw_plan_dft_r2c_2d(forw    ,nx ,ny ,data    ,data    ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back    ,nx ,ny ,data    ,data    ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_2d(forw_big,nx2,ny2,data_big,data_big,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back_big,nx2,ny2,data_big,data_big,FFTW_PATIENT,FFTW_UNALIGNED)

deallocate(data)
deallocate(data_big)

call init_wavenumber()
end subroutine init_fft

!**********************************************************************
subroutine init_wavenumber()
!**********************************************************************
use param,only:lh,ny,L_x,L_y,pi
implicit none
integer :: jx,jy

! Allocate wavenumbers
allocate( kx(lh,ny), ky(lh,ny), k2(lh,ny) )

do jx=1,lh-1
    kx(jx,:) = real(jx-1,kind=rprec)
end do

do jy=1,ny
    ky(:,jy) = real(modulo(jy - 1 + ny/2,ny) - ny/2,kind=rprec)
end do

! Nyquist: makes doing derivatives easier
kx(lh,:)=0._rprec
ky(lh,:)=0._rprec
kx(:,ny/2+1)=0._rprec
ky(:,ny/2+1)=0._rprec

! for the aspect ratio change
kx=2._rprec*pi/L_x*kx
ky=2._rprec*pi/L_y*ky 

! magnitude squared: will have 0's around the edge
k2 = kx*kx + ky*ky

end subroutine init_wavenumber

end module fft
