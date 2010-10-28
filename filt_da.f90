!**********************************************************************
subroutine filt_da(f,dfdx,dfdy)
!**********************************************************************
!
! kills the oddball components in f, and calculates in horizontal derivatives
!--supplies results for jz=$lbz:nz
!--MPI: all these arrays are 0:nz
!--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
!  except on bottom process (0 level set to BOGUS, starts at 1)
!
use types,only:rprec
use param,only:ld,lh,nx,ny,nz
use fft
implicit none
integer::jz
! only need complex treatment
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rprec), dimension(ld, ny, $lbz:nz), intent(inout) :: f
real(rprec), dimension(ld, ny, $lbz:nz), intent(out) :: dfdx, dfdy
!complex(rprec), dimension (lh, ny, $lbz:nz) :: f_c, dfdx_c, dfdy_c

real(kind=rprec), parameter ::const = 1._rprec/(nx*ny)

! loop through horizontal slices
do jz=$lbz,nz
  
  f(:,:,jz)=const*f(:,:,jz)   !normalize

  call rfftwnd_f77_one_real_to_complex(forw,f(:,:,jz),null())

  ! what exactly is the filter doing here? in other words, why is routine
  ! called filt_da? not filtering anything

  ! Zero padded region and Nyquist frequency
  !  f_c(lh,:,jz)=0._rprec ! Complex version
  !  f_c(:,ny/2+1,jz)=0._rprec !  Complex version
  f(ld-1:ld,:,jz) = 0._rprec
  f(:,ny/2+1,jz) = 0._rprec   

  !  Compute in-plane derivatives
  !  dfdy_c(:,:,jz)=eye*ky(:,:)*f_c(:,:,jz) !  complex version
  !  dfdx_c(:,:,jz)=eye*kx(:,:)*f_c(:,:,jz) !  complex version
  call emul_complex_mult_real_complex_imag_2D( f(:,:,jz), kx, ld, lh, ny, dfdx(:,:,jz) )
  call emul_complex_mult_real_complex_imag_2D( f(:,:,jz), ky, ld, lh, ny, dfdy(:,:,jz) )

  ! the oddballs for derivatives should already be dead, since they are for f
  ! inverse transform 
  call rfftwnd_f77_one_complex_to_real(back,f(:,:,jz),null())
  call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),null())
  call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),null())

end do

return
end subroutine filt_da
