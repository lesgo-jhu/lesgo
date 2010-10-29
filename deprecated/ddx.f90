!**********************************************************************
subroutine ddx(f,dfdx)              
!**********************************************************************
!
!  This subroutine computes the partial derivative of f with respect to
!  x using spectral decomposition.
!  

use types,only:rprec
use param,only:ld,lh,nx,ny,nz,dz
use fft

implicit none

integer::jz

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

!  Original complex versions
!complex(kind=rprec),dimension(lh,ny,$lbz:nz),intent(in)::f_c
!complex(kind=rprec),dimension(lh,ny,$lbz:nz),intent(inout)::dfdx_c  

real(rprec), dimension(ld,ny,$lbz:nz), intent(in) :: f
real(rprec), dimension(ld,ny,$lbz:nz), intent(out) :: dfdx

real(kind=rprec), parameter ::const = 1._rprec/(nx*ny)

! Loop through horizontal slices
do jz=$lbz,nz

  !  Use dfdx to hold f; since we are doing IN_PLACE FFT's this is required
  dfdx(:,:,jz)=const*f(:,:,jz)
  call rfftwnd_f77_one_real_to_complex(forw,dfdx(:,:,jz),null())       

  ! Zero padded region and Nyquist frequency
  !  dfdx_c(lh,:,jz)=0._rprec ! Complex version
  !  dfdx_c(:,ny/2+1,jz)=0._rprec !  Complex version
  dfdx(ld-1:ld,:,jz) = 0._rprec
  dfdx(:,ny/2+1,jz) = 0._rprec

  ! Compute coefficients for pseudospectral derivative calculation
  !  dfdx_c(:,:,jz)=eye*kx(:,:)*dfdx_c(:,:,jz) ! Complex version

  !  Use complex emulation of dfdx to perform complex multiplication
  !  call emul_complex_mult_real_complex_2D( dfdx(:,:,jz), eye*kx, ld, lh, ny, dfdx(:,:,jz) )
  
  !  Optimized version for real(eye*kx)=0; only passing imaginary part of eye*kx
  call emul_complex_mult_inplace_real_complex_imag_2D( dfdx(:, :, jz), kx, ld, lh, ny )

  ! Perform inverse transform to get pseudospectral derivative
  call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),null())

enddo

return

end subroutine ddx
