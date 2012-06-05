subroutine dealias1 (u,u_big)
! puts an array into the 'big' size
! doesn't dealias anything
! note: this sort of trashes u
use types,only:rprec
use param,only:ld,ld_big,nx,nx2,ny,nz,ny2,dx,dy
use fft
implicit none
integer::jz
real(kind=rprec),dimension(ld,ny,nz),intent(in):: u
real(kind=rprec),dimension(ld_big,ny2,nz),intent(inout)::u_big
real(kind=rprec)::const
! ahh get rid of this
real(kind=rprec),dimension(ld,ny,nz)::temp    
! be careful using u after calling this subroutine!
const=1._rprec/(nx*ny)
temp=const*u

do jz=1,nz
! still need to normalize
  $if ($FFTW3)
  in2(1:nx,1:ny)=temp(1:nx,1:ny,jz)
  call dfftw_execute_dft_r2c(plan_forward,in2(1:nx,1:ny),temp(1:nx+2,1:ny,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw,temp(:,:,jz),fftwNull_p)
  $endif
! check padd syntax
   call padd(u_big(:,:,jz),temp(:,:,jz))
   $if ($FFTW3)
   inbig2(1:nx2,1:ny2)=u_big(1:nx2,1:ny2,jz)
   call dfftw_execute_dft_c2r(plan_backward_big,inbig2(1:nx2+2,1:ny2),   u_big(1:nx2,1:ny2,jz))   
   $else
   call rfftwnd_f77_one_complex_to_real(back_big,u_big(:,:,jz),fftwNull_p)
   $endif
end do
end subroutine dealias1
