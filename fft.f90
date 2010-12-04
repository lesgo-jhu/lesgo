!**********************************************************************
! fftw 2.1.X version
!**********************************************************************
module fft
use types,only:rprec
use param,only:lh,ny
implicit none

save

!public
private
public :: kx,ky,k2,eye,forw,back,forw_big,back_big,init_fft,init_wavenumber
public ::  FFTW_FORWARD, FFTW_BACKWARD,&
     FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE,FFTW_MEASURE,&
     FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM

! plans
integer*8::forw,back,forw_big,back_big
real(kind=rprec),dimension(lh,ny) :: kx,ky,k2
complex(kind=rprec), parameter :: eye = (0._rprec,1._rprec)
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

contains

!**********************************************************************
subroutine init_fft()
!**********************************************************************
use param,only:nx,ny,nx2,ny2
implicit none
! formulate the fft plans--may want to use FFTW_USE_WISDOM
call rfftw2d_f77_create_plan(forw,nx,ny,FFTW_REAL_TO_COMPLEX,&
     FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call rfftw2d_f77_create_plan(back,nx,ny,FFTW_COMPLEX_TO_REAL,&
     FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call rfftw2d_f77_create_plan(forw_big,nx2,ny2,&
     FFTW_REAL_TO_COMPLEX,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call rfftw2d_f77_create_plan(back_big,nx2,ny2,&
     FFTW_COMPLEX_TO_REAL,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_THREADSAFE)
call init_wavenumber()
end subroutine init_fft

!**********************************************************************
subroutine init_wavenumber()
!**********************************************************************
use param,only:lh,nx,ny,L_x,L_y,pi
implicit none
integer :: jx,jy

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
