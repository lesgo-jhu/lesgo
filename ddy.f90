!**********************************************************************
subroutine ddy(f,dfdy)              
!**********************************************************************
!
!  This subroutine computes the partial derivative of f with respect to
!  y using spectral decomposition.
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
real(rprec), dimension(ld,ny,$lbz:nz), intent(out) :: dfdy

real(kind=rprec), parameter ::const = 1._rprec/(nx*ny)

! Loop through horizontal slices
do jz=$lbz,nz    

  !  Use dfdy to hold f; since we are doing IN_PLACE FFT's this is required
  dfdy(:,:,jz)=const*f(:,:,jz)
  call rfftwnd_f77_one_real_to_complex(forw,dfdy(:,:,jz),null())     

  ! Zero padded region and Nyquist frequency
  !  dfdy_c(lh,:,jz)=0._rprec ! Complex version
  !  dfdy_c(:,ny/2+1,jz)=0._rprec !  Complex version
  dfdy(ld-1:ld,:,jz) = 0._rprec
  dfdy(:,ny/2+1,jz) = 0._rprec   

  ! Compute coefficients for pseudospectral derivative calculation
  !  dfdy_c(:,:,jz)=eye*ky(:,:)*dfdy_c(:,:,jz) ! Complex version

  !  Use complex emulation of dfdy to perform complex multiplication
  !  call emul_complex_mult_real_complex_2D( dfdy(:,:,jz), eye*ky, ld, lh, ny, dfdy(:,:,jz) )
  
  !  Optimized version for real(eye*ky)=0; only passing imaginary part of eye*ky
  call emul_complex_mult_real_complex_imag_2D( dfdy(:, :, jz), ky, ld, lh, ny, dfdy(:,:,jz) )

  ! Perform inverse transform to get pseudospectral derivative
  call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),null())   

end do

return

end subroutine ddy
