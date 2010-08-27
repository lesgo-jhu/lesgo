! kills the oddball components in f, and calculates in horizontal derivatives
!--supplies results for jz=$lbz:nz
  !--MPI: all these arrays are 0:nz
  !--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
  !  except on bottom process (0 level set to BOGUS, starts at 1)
subroutine filt_da(f_c,dfdx_c,dfdy_c)
use types,only:rprec
use param,only:lh,nx,ny,nz
use fft
implicit none
integer::jz
! only need complex treatment
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
complex(rprec), dimension (lh, ny, $lbz:nz) :: f_c, dfdx_c, dfdy_c
real(kind=rprec)::const, ignore_me
! loop through horizontal slices
const=1._rprec/(nx*ny)
do jz=$lbz,nz
   f_c(:,:,jz)=const*f_c(:,:,jz)   !normalize
   call rfftwnd_f77_one_real_to_complex(forw,f_c(:,:,jz),ignore_me)
! what exactly is the filter doing here? in other words, why is routine
! called filt_da? not filtering anything
! kill the oddballs
   f_c(lh,:,jz)=0._rprec
   f_c(:,ny/2+1,jz)=0._rprec
! in-plane derivatives
   dfdy_c(:,:,jz)=eye*ky(:,:)*f_c(:,:,jz)
   dfdx_c(:,:,jz)=eye*kx(:,:)*f_c(:,:,jz)
! the oddballs for derivatives should already be dead, since they are for f
! inverse transform 
   call rfftwnd_f77_one_complex_to_real(back,f_c(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back,dfdx_c(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back,dfdy_c(:,:,jz),ignore_me)
end do
end subroutine filt_da
