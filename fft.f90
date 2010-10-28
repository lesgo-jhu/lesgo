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
public :: kx,ky,k2,eye,forw,back,forw_big,back_big,init_fft
public ::  FFTW_FORWARD, FFTW_BACKWARD,&
     FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE,FFTW_MEASURE,&
     FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
! plans
integer*8::forw,back,forw_big,back_big
real(kind=rprec),dimension(lh,ny)::kx,ky,k2
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

!**********************************************************************
subroutine emul_complex_mult_real_complex_2D( a, a_c, nx_r, nx_c, ny, b )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type.
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (complex,size(nx_c,ny))- input complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension( nx_r, ny ), intent(in) :: a
complex(rprec), dimension( nx_c, ny), intent(in) :: a_c

real(rprec), dimension(nx_r, ny), intent(out) :: b

!  Cached variables
real(rprec) ::  a_r, a_i, a_c_r, a_c_i

integer :: i,j,ir,ii

!  Emulate complex multiplication
!  Using outer loop to get contingious memory access
do j=1, ny
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)
    a_c_r = real(a_c(i,j),kind=rprec)
    $if($DBLPREC)
    a_c_i = dimag(a_c(i,j))
    $else
    a_c_i = aimag(a_c(i,j))
    $endif
    
    !  Perform multiplication
    b(ir,j) = a_r * a_c_r - a_i * a_c_i
    b(ii,j) = a_r * a_c_i + a_i * a_c_r

  enddo
enddo

return

end subroutine emul_complex_mult_real_complex_2D

!**********************************************************************
subroutine emul_complex_mult_inplace_real_complex_2D( a, a_c, nx_r, nx_c, ny )
!**********************************************************************
!
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. The resulting
!  multiplication is returned as array a
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input/output real array
!    a_c (complex,size(nx_c,ny))- input complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension( nx_r, ny ), intent(inout) :: a
complex(rprec), dimension( nx_c, ny), intent(in) :: a_c

!  Cached variables
real(rprec) ::  a_r, a_i, a_c_r, a_c_i

integer :: i,j,ir,ii

!  Emulate complex multiplication
!  Using outer loop to get contingious memory access
do j=1, ny
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)

    a_c_r = real(a_c(i,j),kind=rprec)
    $if($DBLPREC)
    a_c_i = dimag(a_c(i,j))
    $else
    a_c_i = aimag(a_c(i,j))
    $endif
    
    !  Perform multiplication
    a(ir,j) = a_r * a_c_r - a_i * a_c_i
    a(ii,j) = a_r * a_c_i + a_i * a_c_r

  enddo
enddo

return

end subroutine emul_complex_mult_inplace_real_complex_2D

!**********************************************************************
subroutine emul_complex_mult_real_complex_imag_2D( a, a_c, nx_r, nx_c, ny, b )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0)
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (real,size(nx_c,ny))   - input imaginary part of complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), target, dimension(nx_r, ny), intent(in) :: a
real(rprec), target, dimension(nx_c, ny), intent(in) :: a_c

real(rprec), target, dimension(nx_r, ny), intent(out) :: b

!  Cached variables
real(rprec) ::  a_c_i

integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contingious memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_c_i = a_c(i,j)

    !  Perform multiplication
    b(ir,j) = - a(ii,j) * a_c_i
    b(ii,j) =  a(ir,j) * a_c_i

  enddo
enddo

return

end subroutine emul_complex_mult_real_complex_imag_2D

!**********************************************************************
subroutine emul_complex_mult_inplace_real_complex_imag_2D( a, a_c, nx_r, nx_c, ny )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0)
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (real,size(nx_c,ny))   - input imaginary part of complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), target, dimension(nx_r, ny), intent(inout) :: a
real(rprec), target, dimension(nx_c, ny), intent(in) :: a_c

!  Cached variables
real(rprec) :: a_r, a_i, a_c_i

integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contingious memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)
    a_c_i = a_c(i,j)

    !  Perform multiplication
    a(ir,j) = - a_i * a_c_i
    a(ii,j) =  a_r * a_c_i

  enddo
enddo

return

end subroutine emul_complex_mult_inplace_real_complex_imag_2D

end module fft
