subroutine dealias2(u,u_big,lbz)
! puts back into small array
use types,only:rprec
use param,only:ld,ld_big,nx,ny,nz,nx2,ny2,dx,dy
use fft
implicit none
integer::jz

integer, intent(in) :: lbz
real(kind=rprec),dimension(ld,ny,lbz:nz),intent(in):: u
real(kind=rprec),dimension(ld_big,ny2,lbz:nz),intent(inout)::u_big
real(kind=rprec)::ignore_me,const

! normalize
const=1._rprec/(nx2*ny2)
u_big=const*u_big
! Loop through horizontal slices

do jz=lbz,nz
! perform forward FFT
   call rfftwnd_f77_one_real_to_complex(forw_big,u_big(:,:,jz),ignore_me)    
   call unpadd(u(:,:,jz),u_big(:,:,jz))
! Back to physical space
   call rfftwnd_f77_one_complex_to_real(back,u(:,:,jz),ignore_me)
end do

! sc: do we need this?
!.....Making filter circular !!!!!!!!!!!
!	call filt_da(uu)
!.......................................
end subroutine dealias2
