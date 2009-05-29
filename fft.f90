!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fftw 2.1.X version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fft
use types,only:rprec
use param2,only:lh,ny
implicit none

save

!public
!private
!public :: kx,ky,k2,eye,forw,back,forw_big,back_big
!public ::  FFTW_FORWARD, FFTW_BACKWARD,&
!     FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE,FFTW_MEASURE,&
!     FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM,FFTW_THREADSAFE
! plans
integer*8::forw,back,forw_big,back_big
real(kind=rprec), allocatable, dimension(:,:)::kx,ky,k2
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

end module fft




