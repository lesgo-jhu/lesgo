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
! fftw 2.1.X version
!**********************************************************************
module fft
use types,only:rprec
use param,only:ld,lh,ny,ld_big, ny2, spectra_calc
use iso_c_binding
implicit none
include'fftw3.f'
save

public :: kx, ky, k2, init_fft, forw_spectra
public :: forw, back, forw_big, back_big

real(rprec), allocatable, dimension(:,:) :: kx, ky, k2
integer*8::forw_spectra
integer*8::forw,back,forw_big,back_big

real (rprec), dimension (:, :), allocatable :: data, data_big

contains

!**********************************************************************
subroutine init_fft()
!**********************************************************************
use param,only:nx,ny,nx2,ny2
implicit none

! Allocate temporary arrays for creating the FFTW plans
allocate( data(ld, ny) )
allocate( data_big(ld_big, ny2) )

! Initialize and implement with threads
!call dfftw_init_threads(iret)
!call dfftw_plan_with_nthreads(2)

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
use param,only:lh,nx,ny,L_x,L_y,pi
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
