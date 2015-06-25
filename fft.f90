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
use param,only:ld,lh,ny,ld_big, ny2, spectra_calc, kx_num, kx_vec
$if ($FFTW3)
!use, intrinsic :: iso_c_binding
use iso_c_binding
$endif
implicit none
$if ($FFTW3)
include'fftw3.f'
!include 'fftw3-mpi.f'
$endif
save

public :: kx, ky, k2, init_fft, forw_spectra
public :: expp, expn    !! jb
!!public :: ky_vec                                           !!jb
public :: forw, back, forw_big, back_big
public :: forw_1d, back_1d, forw_complex, forw_y, back_y   !!jb

real(rprec), allocatable, dimension(:,:) :: kx, ky, k2
complex(rprec), allocatable, dimension(:,:) :: expp, expn      !!jb
!!real(rprec), allocatable, dimension(:) :: ky_vec            !!jb
integer, allocatable, dimension(:) :: kx_veci            !!jb
integer*8::forw_spectra
integer*8::forw,back,forw_big,back_big
integer*8::forw_1d, back_1d, forw_complex, forw_y, back_y   !!jb

$if ($FFTW3)
real (rprec), dimension (:, :), allocatable :: data, data_big
real (rprec), dimension (:), allocatable :: data_x, data_in   !!jb
complex (rprec), dimension (:), allocatable :: data_out       !!jb
real (rprec), dimension (:), allocatable :: data_in_y         !!jb
complex (rprec), dimension (:), allocatable :: data_out_y     !!jb
$else

!public
private

public ::  FFTW_FORWARD, FFTW_BACKWARD,&
     FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE,FFTW_MEASURE,&
     FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
public :: fftwNull_p

! plans

! fftw 2.1.3 stuff
integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
integer, parameter :: FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1
ffinteger, parameter :: FFTW_ESTIMATE=0,FFTW_MEASURE=1
integer, parameter :: FFTW_OUT_OF_PLACE=0
integer, parameter :: FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16
integer, parameter :: FFTW_THREADSAFE=128
integer, parameter :: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
integer, parameter :: FFTW_SCRAMBLED_INPUT=8192
integer, parameter :: FFTW_SCRAMBLED_OUTPUT=16384

! Null pointer for fftw2 dummy argument
integer(2), pointer :: fftwNull_p
$endif
contains

!**********************************************************************
subroutine init_fft()
!**********************************************************************
use param,only:nx,ny,nx2,ny2
implicit none

$if ($FFTW3)
! Allocate temporary arrays for creating the FFTW plans
allocate( data(ld, ny) )
allocate( data_big(ld_big, ny2) )
allocate( data_x(ld) )             !!jb
allocate( data_in(nx) )            !!jb
allocate( data_out(nx/2+1) )       !!jb
allocate( data_in_y(ny) )          !!jb
allocate( data_out_y(ny/2+1) )     !!jb

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

!!jb >>>
call dfftw_plan_dft_r2c_1d(forw_1d ,nx   ,data_x  ,data_x  ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(back_1d ,nx   ,data_x  ,data_x  ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_complex ,nx   ,data_in  ,data_out ,FFTW_PATIENT,FFTW_UNALIGNED)

call dfftw_plan_dft_r2c_1d(forw_y, ny, data_in_y, data_out_y, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(back_y, ny, data_out_y, data_in_y, FFTW_PATIENT, FFTW_UNALIGNED)
!! <<<<

deallocate(data)
deallocate(data_big)
deallocate(data_x)
deallocate(data_in)
deallocate(data_out)
deallocate(data_in_y)
deallocate(data_out_y)

$else
call rfftw2d_f77_create_plan(forw,nx,ny,FFTW_REAL_TO_COMPLEX,&
     FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call rfftw2d_f77_create_plan(back,nx,ny,FFTW_COMPLEX_TO_REAL,&
     FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call rfftw2d_f77_create_plan(forw_big,nx2,ny2,&
     FFTW_REAL_TO_COMPLEX,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call rfftw2d_f77_create_plan(back_big,nx2,ny2,&
     FFTW_COMPLEX_TO_REAL,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)

if(spectra_calc) then
  call rfftw_f77_create_plan(forw_spectra, Nx, FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE)
endif

$endif

call init_wavenumber()
end subroutine init_fft

!**********************************************************************
subroutine init_wavenumber()
!**********************************************************************
use param,only:lh,nx,ny,L_x,L_y,pi,kx_limit,kx_allow,kx_dft,coord
implicit none
integer :: jx,jy,ii
complex(rprec) :: c1   !!jb
complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)   !!jb

! Allocate wavenumbers
allocate( kx(lh,ny), ky(lh,ny), k2(lh,ny) )

if (kx_dft) then
  allocate( expp(kx_num, nx), expn(kx_num, nx) )
  !allocate( ky_vec( ny ) )
  allocate( kx_veci ( kx_num ) )    !!jb

  kx_veci = int( kx_vec ) + 1    !!jb
  kx_vec = kx_vec * 2._rprec * pi / L_x     !!jb , aspect ratio change

  c1 = L_x / real(nx,rprec) * imag   !!jb
  do jx=1,kx_num
  do ii=1,nx
  expp(jx,ii) = exp( c1 * kx_vec(jx) * real(ii-1,rprec) )
  expn(jx,ii) = exp(-c1 * kx_vec(jx) * real(ii-1,rprec) )
  enddo
  enddo
endif

do jx=1,lh-1
   kx(jx,:) = real(jx-1,kind=rprec)
end do

!do jx=1,kx_num
!   kx(jx,:) = kx_vec(jx)
!end do

if (kx_limit) then     !!jb
   kx(2,:) = real(kx_allow,kind=rprec)
endif

do jy=1,ny
   ky(:,jy) = real(modulo(jy - 1 + ny/2,ny) - ny/2,kind=rprec)
end do

!do jy=1,ny/2      !!jb
!   ky_vec(jy) = real( jy-1 , kind=rprec )
!   ky_vec(ny - jy + 2) =  -ky_vec(jy)
!end do
!ky_vec(ny/2+1) = real(ny/2, rprec)

! Nyquist: makes doing derivatives easier
      kx(lh,:)=0._rprec
      ky(lh,:)=0._rprec
      kx(:,ny/2+1)=0._rprec
      ky(:,ny/2+1)=0._rprec

! for the aspect ratio change
      kx=2._rprec*pi/L_x*kx
      ky=2._rprec*pi/L_y*ky 
      !ky_vec=2._rprec*pi/L_y*ky_vec     !!jb 

! magnitude squared: will have 0's around the edge
      k2 = kx*kx + ky*ky

if (coord == 0) then
   write (*,*) 'jx, jy, kx, ky, k2 =========================== '
   do jx=1,lh
      do jy=1,ny
         write(*,*) jx,jy,kx(jx,jy),ky(jx,jy),k2(jx,jy)
      enddo
   enddo
endif


end subroutine init_wavenumber

end module fft
