subroutine ddx(dfdx, f)              
! f remains untouched
use types,only:rprec
use param,only:lh,nx,ny,nz,dz
use fft
implicit none
integer::jz
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
complex(kind=rprec),dimension(lh,ny,$lbz:nz),intent(in)::f
complex(kind=rprec),dimension(lh,ny,$lbz:nz),intent(inout)::dfdx
real(kind=rprec)::const, ignore_me

const=1._rprec/(nx*ny)
! Loop through horizontal slices
do jz=$lbz,nz
! temporary storage in dfdx--even though they're complex arrays!
   dfdx(:,:,jz)=const*f(:,:,jz)
   call rfftwnd_f77_one_real_to_complex(forw,dfdx(:,:,jz),ignore_me)       
! oddballs
   dfdx(lh,:,jz)=0._rprec
   dfdx(:,ny/2+1,jz)=0._rprec
! compute coefficients for pseudospectral derivative calculation
   dfdx(:,:,jz)=eye*kx(:,:)*dfdx(:,:,jz)
! inverse transform to get pseudospectral derivative
   call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),ignore_me)
end do

return
end subroutine ddx
