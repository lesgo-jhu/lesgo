subroutine ddxy (dfdx_c, dfdy_c, f_c)              
use types,only:rprec
use param2,only:lh,nx,ny,nz,dz
use fft
implicit none
integer::jz
! only need complex treatment
complex(kind=rprec),dimension(lh,ny,nz),intent(in)::f_c
complex(kind=rprec),dimension(lh,ny,nz),intent(inout)::dfdx_c,dfdy_c
real(kind=rprec)::const, ignore_me

const=1._rprec/(nx*ny)
!...Loop through horizontal slices
do jz=1,nz
! temporay storage in dfdx_c, this was don't mess up f_c
   dfdx_c(:,:,jz)=const*f_c(:,:,jz)   !normalize
   call rfftwnd_f77_one_real_to_complex(forw,dfdx_c(:,:,jz),ignore_me)
   dfdx_c(lh,:,jz)=0._rprec
   dfdx_c(:,ny/2+1,jz)=0._rprec
! derivatives: must to y's first here, because we're using dfdx as storage
   dfdy_c(:,:,jz)=eye*ky(:,:)*dfdx_c(:,:,jz)
   dfdx_c(:,:,jz)=eye*kx(:,:)*dfdx_c(:,:,jz)
! the oddballs for derivatives should already be dead, since they are for f

! inverse transform 
   call rfftwnd_f77_one_complex_to_real(back,dfdx_c(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back,dfdy_c(:,:,jz),ignore_me)
end do
end subroutine ddxy
