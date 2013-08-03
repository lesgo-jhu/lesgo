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
use param,only:lh,ny,spectra_calc
$if ($FFTW3)
!use, intrinsic :: iso_c_binding
use iso_c_binding
$endif
implicit none
$if ($FFTW3)
!include'fftw3.f'
!include'fftw3-mpi.f03'
!include'fftw3.f'
include'fftw3.f03'
!include'fftw3l.f03'
!include 'fftw3-mpi.f03'
$endif
save

public :: kx, ky, k2, init_fft, forw_spectra
real(rprec), allocatable, dimension(:,:) :: kx, ky, k2
integer*8::forw_spectra

!RICHARD:FFTW3
$if ($FFTW3)
public :: plan_forward,plan_backward,plan_forward_big,plan_backward_big
public :: in,inbig,out,outbig
public :: in2,inbig2,out2,outbig2
public :: cdata,cdata2

real (rprec), pointer :: in(:,:)
real (rprec), pointer :: inbig(:,:)
real (rprec), pointer :: out(:,:)
real (rprec), pointer :: outbig(:,:)

real (rprec), dimension (:, :), allocatable :: in2
real (rprec), dimension (:, :), allocatable :: inbig2
real (rprec), dimension (:, :), allocatable :: out2
real (rprec), dimension (:, :), allocatable :: outbig2

integer*8::iret,plan_forward,plan_backward,plan_forward_big,plan_backward_big
type(C_PTR) :: cdata
type(C_PTR) :: cdata2
$else

!public
private
public :: forw, back, forw_big, back_big
public ::  FFTW_FORWARD, FFTW_BACKWARD,&
     FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE,FFTW_MEASURE,&
     FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
public :: fftwNull_p

! plans
integer*8::forw,back,forw_big,back_big
! fftw 2.1.3 stuff
integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
integer, parameter :: FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1
integer, parameter :: FFTW_ESTIMATE=0,FFTW_MEASURE=1
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

!RICHARD FFTW3
$if ($FFTW3)
logical, save :: arrays_allocatedfftw3 = .false. 
if( .not. arrays_allocatedfftw3 ) then
   allocate(in2     (nx  , ny ))  
   allocate(inbig2  (nx2 , ny2))
   allocate(out2    (nx+2, ny ))
   allocate(outbig2 (nx2+2, ny2))
   arrays_allocatedfftw3 = .true. 
endif

!call fftw_mpi_init()

cdata = fftw_alloc_real(int(nx*ny, C_SIZE_T))
cdata2= fftw_alloc_real(int((nx+2)*ny, C_SIZE_T))
call c_f_pointer(cdata, in, [nx,ny])
call c_f_pointer(cdata2,out, [nx+2,ny])
$endif

$if ($FFTW3) 
call dfftw_init_threads(iret)
call dfftw_plan_with_nthreads(2)

call dfftw_plan_dft_r2c_2d(plan_forward     ,nx ,ny ,in2     ,out2   ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(plan_backward    ,nx ,ny ,out2    ,in2    ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_2d(plan_forward_big ,nx2,ny2,inbig2  ,outbig2,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(plan_backward_big,nx2,ny2,outbig2 ,inbig2 ,FFTW_PATIENT,FFTW_UNALIGNED)
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
