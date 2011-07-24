!**********************************************************************
subroutine ddy (dfdy, f)              
!**********************************************************************
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
complex(kind=rprec),dimension(lh,ny,$lbz:nz),intent(inout)::dfdy
real(kind=rprec)::const, ignore_me

const=1._rprec/(nx*ny)   
! Loop through horizontal slices
do jz=$lbz,nz
! temporary storage in dfdx
   dfdy(:,:,jz)=const*f(:,:,jz)
   call rfftwnd_f77_one_real_to_complex(forw,dfdy(:,:,jz),ignore_me)       
! oddballs
   dfdy(lh,:,jz)=0._rprec
   dfdy(:,ny/2+1,jz)=0._rprec
! compute coefficients for pseudospectral derivative calculation.
   dfdy(:,:,jz)=eye*ky(:,:)*dfdy(:,:,jz)
! inverse transform to get pseudospectral derivative.
   call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),ignore_me)
end do
end subroutine ddy
