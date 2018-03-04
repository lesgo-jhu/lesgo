subroutine dealias1 (u,u_big,lbz)
! puts an array into the 'big' size
! doesn't dealias anything
! note: this sort of trashes u
use types,only:rprec
use param,only:ld,ld_big,nx,ny,nz,ny2
use fft

implicit none

integer::jz

integer, intent(in) :: lbz
real(rprec), dimension(ld,ny,lbz:nz),intent(in) :: u
real(rprec),dimension(ld_big,ny2,lbz:nz),intent(inout) :: u_big
real(rprec)::ignore_me,const
! ahh get rid of this
real(rprec),allocatable,dimension(:,:,:) :: temp    


allocate(temp(ld,ny,lbz:nz))


! be careful using u after calling this subroutine!
const=1._rprec/(nx*ny)
temp=const*u

do jz=lbz,nz
   temp(:,:,jz) = const*u(:,:,jz)
! still need to normalize
   call rfftwnd_f77_one_real_to_complex(forw,temp(:,:,jz),ignore_me)
!   check padd syntax
   call padd(u_big(:,:,jz),temp(:,:,jz))
   call rfftwnd_f77_one_complex_to_real(back_big,u_big(:,:,jz),ignore_me)
end do

end subroutine dealias1
