subroutine dealias2(u,u_big)
! puts back into small array
use types,only:rprec
use param,only:ld,ld_big,nx,ny,nz,nx2,ny2,dx,dy
use fft
implicit none
integer::jz
real(kind=rprec),dimension(ld,ny,nz),intent(in):: u
real(kind=rprec),dimension(ld_big,ny2,nz),intent(inout)::u_big
real(kind=rprec)::const

! normalize
const=1._rprec/(nx2*ny2)
u_big=const*u_big
! Loop through horizontal slices
do jz=1,nz
! perform forward FFT
  $if ($FFTW3)
  inbig2(1:nx2,1:ny2)=u_big(1:nx2,1:ny2,jz)
  call dfftw_execute_dft_r2c(plan_forward,inbig2(1:nx2,1:ny2),u_big(1:nx2+2,1:ny2,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw_big,u_big(:,:,jz),fftwNull_p)
  $endif
  call unpadd(u(:,:,jz),u_big(:,:,jz))
! Back to physical space
   $if ($FFTW3)
   in2(1:nx,1:ny)=u(1:nx,1:ny,jz)
   call dfftw_execute_dft_c2r(plan_backward_big,in2(1:nx+2,1:ny),   u(1:nx,1:ny,jz))     
   $else
   call rfftwnd_f77_one_complex_to_real(back,u(:,:,jz),fftwNull_p)
   $endif
end do

! sc: do we need this?
!.....Making filter circular !!!!!!!!!!!
!	call filt_da(uu)
!.......................................
end subroutine dealias2
