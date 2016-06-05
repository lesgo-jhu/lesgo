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
use types, only: rprec
use param, only: ld,lh,ny,ld_big, ny2, spectra_calc, kx_num, kx_vec
use param, only: fourier, kxs_in, kxs
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
public :: expp, expn, expp_big, expn_big    !! jb
public :: expp_n, expn_n, expp_big_n, expn_big_n    !! jb
!!public :: ky_vec                                           !!jb
public :: forw, back, forw_big, back_big
public :: forw_1d, back_1d, forw_y, back_y   !!jb
public :: forw_complex, forw_complex_fourier  !!jb
public :: forw_y_big, back_y_big   !!jb
public :: forw_span_spectra   !! jb
public :: forw_span_spectra2  !! jb
public :: ccomp_forw, rcomp_forw        !! jb
public :: ccomp_back, rcomp_back        !! jb
public :: ycomp_back, ycomp_forw        !! jb
public :: ycomp_back_big, ycomp_forw_big        !! jb
public :: forw_fourier, back_fourier     !!jb

real(rprec), allocatable, dimension(:,:) :: kx, ky, k2
complex(rprec), allocatable, dimension(:,:) :: expp, expn      !!jb
complex(rprec), allocatable, dimension(:,:) :: expp_big, expn_big      !!jb
complex(rprec), allocatable, dimension(:,:) :: expp_n, expn_n      !!jb
complex(rprec), allocatable, dimension(:,:) :: expp_big_n, expn_big_n      !!jb
!!real(rprec), allocatable, dimension(:) :: ky_vec            !!jb
integer, allocatable, dimension(:) :: kx_veci            !!jb
integer, allocatable, dimension(:) :: kx_veci_n            !!jb
integer, allocatable, dimension(:) :: kx_veci_n_big            !!jb
integer*8::forw_spectra, forw_span_spectra, forw_span_spectra2
integer*8::forw,back,forw_big,back_big
integer*8::forw_1d, back_1d, forw_y, back_y   !!jb
integer*8::forw_complex, forw_complex_fourier  !!jb
integer*8::forw_y_big, back_y_big   !!jb
integer*8::ccomp_forw, rcomp_forw    !!jb
integer*8::ccomp_back, rcomp_back    !!jb
integer*8::ycomp_back, ycomp_forw    !!jb
integer*8::ycomp_back_big, ycomp_forw_big    !!jb
integer*8::forw_fourier, back_fourier

$if ($FFTW3)
real (rprec), dimension (:, :), allocatable :: data, data_big
real (rprec), dimension (:, :), allocatable :: data_fourier  !!jb
real (rprec), dimension (:), allocatable :: data_x, data_in   !!jb
real (rprec), dimension (:), allocatable :: data_in_fourier   !!jb
complex (rprec), dimension (:), allocatable :: data_out       !!jb
complex (rprec), dimension (:), allocatable :: data_out_fourier       !!jb
real (rprec), dimension (:), allocatable :: data_in_y         !!jb
complex (rprec), dimension (:), allocatable :: data_out_y     !!jb
!real (rprec), dimension (:), allocatable :: data_out_y     !!jb
real (rprec), dimension (:), allocatable :: data_in_y_big         !!jb
complex (rprec), dimension (:), allocatable :: data_out_y_big     !!jb
real (rprec), dimension(:), allocatable :: span_in
complex (rprec), dimension(:), allocatable :: span_out
real (rprec), dimension(:), allocatable :: span_in2
real (rprec), dimension(:), allocatable :: span_out2
real (rprec), dimension(:), allocatable :: rcomp_data
complex (rprec), dimension(:), allocatable :: ccomp_data
complex (rprec), dimension(:), allocatable :: ycomp_data
complex (rprec), dimension(:), allocatable :: ycomp_data_big
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
use param,only:nx,ny,nx2,ny2,span_spectra_calc,nxp
implicit none

$if ($FFTW3)
! Allocate temporary arrays for creating the FFTW plans
allocate( data(ld, ny) )
allocate( data_big(ld_big, ny2) )
allocate( data_x(ld) )             !!jb
allocate( data_in(nx) )            !!jb
allocate( data_out(nx/2+1) )       !!jb
allocate( data_in_fourier(nxp) )            !!jb
allocate( data_out_fourier(nxp/2+1) )       !!jb
allocate( data_in_y(ny) )          !!jb
allocate( data_out_y(ny/2+1) )     !!jb
!allocate( data_out_y(ny) )     !!jb
allocate( data_in_y_big(ny2) )          !!jb
allocate( data_out_y_big(ny2/2+1) )     !!jb
allocate( span_in(ny) )
allocate( span_out(ny/2+1) )

allocate( span_in2 ( nx ) )   !ny-2
allocate( span_out2 ( ld ) )    !ny

allocate( ccomp_data(ny) )    !! complex   !!jb
allocate( rcomp_data(2*ny) )    !! complex stored as real   !!jb
allocate( ycomp_data(ny) )    !! complex stored as real   !!jb
allocate( ycomp_data_big(ny2) )    !! complex stored as real   !!jb

allocate( data_fourier(nxp+2, ny) )

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
call dfftw_plan_dft_r2c_2d(forw_fourier    ,nxp ,ny ,data_fourier    ,data_fourier    ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back_fourier    ,nxp ,ny ,data_fourier    ,data_fourier    ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_1d ,nx   ,data_x  ,data_x  ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(back_1d ,nx   ,data_x  ,data_x  ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_complex, nx, data_in, data_out, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_complex_fourier, nxp, data_in_fourier, data_out_fourier, FFTW_PATIENT, FFTW_UNALIGNED)

call dfftw_plan_dft_r2c_1d(forw_y, ny, data_in_y, data_out_y, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(back_y, ny, data_out_y, data_in_y, FFTW_PATIENT, FFTW_UNALIGNED)

call dfftw_plan_dft_r2c_1d(forw_y_big, ny2, data_in_y_big, data_out_y_big, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_1d(back_y_big, ny2, data_out_y_big, data_in_y_big, FFTW_PATIENT, FFTW_UNALIGNED)

!if (span_spectra_calc) then
call dfftw_plan_dft_r2c_1d(forw_span_spectra , ny, span_in, span_out, FFTW_PATIENT, FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_span_spectra2, nx, span_in2, span_out2, FFTW_PATIENT, FFTW_UNALIGNED) !ny-2
call dfftw_plan_dft_1d(ccomp_forw, ny, ccomp_data, ccomp_data, FFTW_FORWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
call dfftw_plan_dft_1d(rcomp_forw, ny, rcomp_data, rcomp_data, FFTW_FORWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
call dfftw_plan_dft_1d(ccomp_back, ny, ccomp_data, ccomp_data, FFTW_BACKWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
call dfftw_plan_dft_1d(rcomp_back, ny, rcomp_data, rcomp_data, FFTW_BACKWARD, FFTW_PATIENT, FFTW_UNALIGNED) 

call dfftw_plan_dft_1d(ycomp_forw, ny, ycomp_data, ycomp_data, FFTW_FORWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
call dfftw_plan_dft_1d(ycomp_back, ny, ycomp_data, ycomp_data, FFTW_BACKWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
call dfftw_plan_dft_1d(ycomp_forw_big, ny2, ycomp_data_big, ycomp_data_big, FFTW_FORWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
call dfftw_plan_dft_1d(ycomp_back_big, ny2, ycomp_data_big, ycomp_data_big, FFTW_BACKWARD, FFTW_PATIENT, FFTW_UNALIGNED) 
!endif

!! <<<<

deallocate(data)
deallocate(data_big)
deallocate(data_fourier)
deallocate(data_x)
deallocate(data_in)
deallocate(data_out)
deallocate(data_in_fourier)
deallocate(data_out_fourier)
deallocate(data_in_y)
deallocate(data_out_y)
deallocate(data_in_y_big)
deallocate(data_out_y_big)
deallocate(span_in)
deallocate(span_out)
deallocate(span_in2)
deallocate(span_out2)
deallocate(ccomp_data)
deallocate(rcomp_data)
deallocate(ycomp_data)
deallocate(ycomp_data_big)

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
use param, only: lh, nx, ny, L_x, L_y, pi, coord, nx2
use param, only: kxs_in, kxs, kx_limit, kx_allow, kx_dft
implicit none
integer :: jx,jy,ii
complex(rprec) :: c1, c1_big   !!jb
complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)   !!jb

! Allocate wavenumbers
!!if (kx_dft) then             !!jb
!!  allocate( kx(kx_num,ny), ky(kx_num,ny), k2(kx_num,ny) )
!!else
  allocate( kx(lh,ny), ky(lh,ny), k2(lh,ny) )
!!endif
  
  allocate( kxs ( kx_num ) )    !!jb
  allocate( kx_veci ( kx_num ) )    !!jb
  allocate( kx_veci_n ( kx_num + kx_num - 2 ) )    !!jb
  allocate( kx_veci_n_big ( kx_num + kx_num - 2 ) )    !!jb
  kx_veci = int( kx_vec ) + 1    !!jb
  kx_veci_n(1:kx_num) = kx_veci(1:kx_num)
  kx_veci_n_big(1:kx_num) = kx_veci(1:kx_num)
  kx_veci_n(kx_num+1 : kx_num+kx_num-2) = nx+2 - kx_veci(2:kx_num-1)
  kx_veci_n_big(kx_num+1 : kx_num+kx_num-2) = nx2+2 - kx_veci(2:kx_num-1)

if (kx_dft) then
  allocate( expp(kx_num, nx), expn(kx_num, nx) )
  allocate( expp_big(kx_num, nx2), expn_big(kx_num, nx2) )
  allocate( expp_n(nx, nx), expn_n(nx, nx) )
  allocate( expp_big_n(nx2, nx2), expn_big_n(nx2, nx2) )
  !allocate( ky_vec( ny ) )

  kx_vec = kx_vec * 2._rprec * pi / real(L_x,rprec)     !!jb , aspect ratio change
  kxs(:) = int( kxs_in(:) )

  c1     = real(L_x, rprec) / real(nx,rprec) * imag   !!jb
  c1_big = real(L_x,rprec) / real(nx2,rprec) * imag   !!jb
  do jx=1,kx_num
  do ii=1,nx
  expp(jx,ii) = exp( c1 * kx_vec(jx) * real(ii-1,rprec) )
  expn(jx,ii) = exp(-c1 * kx_vec(jx) * real(ii-1,rprec) )
  expp_big(jx,ii) = exp( c1_big * kx_vec(jx) * real(ii-1,rprec) )
  expn_big(jx,ii) = exp(-c1_big * kx_vec(jx) * real(ii-1,rprec) )
  enddo
  enddo
  do jx=1,kx_num
  do ii=1,nx2
  expp_big(jx,ii) = exp( c1_big * kx_vec(jx) * real(ii-1,rprec) )
  expn_big(jx,ii) = exp(-c1_big * kx_vec(jx) * real(ii-1,rprec) )
  enddo
  enddo

  do jx=1,nx
  do ii=1,nx
  expp_n(jx,ii) = exp( c1 * real(jx-1,rprec) * real(ii-1,rprec) * 2._rprec * pi / real(L_x,rprec) )
  expn_n(jx,ii) = exp(-c1 * real(jx-1,rprec) * real(ii-1,rprec) * 2._rprec * pi / real(L_x,rprec) )
  enddo
  enddo
  do jx=1,nx2
  do ii=1,nx2
  expp_big_n(jx,ii) = exp( c1_big * real(jx-1,rprec) * real(ii-1,rprec) * 2._rprec * pi / real(L_x,rprec) )
  expn_big_n(jx,ii) = exp(-c1_big * real(jx-1,rprec) * real(ii-1,rprec) * 2._rprec * pi / real(L_x,rprec) )
  enddo
  enddo

endif

do jx=1,lh-1
   kx(jx,:) = real(jx-1,kind=rprec)
end do

!!$kx(1,:) = 0._rprec
!!$kx(2,:) = 2._rprec
!!$kx(3,:) = 3._rprec
!!$kx(4,:) = 7._rprec
!!$kx(5,:) = 0._rprec


!!$if (kx_dft) then   !!jb
!!$ do jx=1,kx_num
!!$    kx(jx,:) = kx_vec(jx)
!!$ end do
!!$endif

if (kx_limit) then     !!jb
   kx(:,:) = real(kx_allow,kind=rprec) * kx(:,:)
endif

!!$if (kx_limit) then     !!jb
!!$   kx(2,:) = real(kx_allow,kind=rprec)
!!$endif

if (fourier) then   !!jb
 do jx=1,kx_num
    kx(jx,:) = kxs_in(jx)
 end do
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
      kx=2._rprec*real(pi/L_x*kx,rprec)
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
