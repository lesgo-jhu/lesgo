subroutine dealias1 (u,u_big)
! puts an array into the 'big' size
! doesn't dealias anything
! note: this sort of trashes u
use types,only:rprec
use param,only:ld,ld_big,nx,ny,nz,ny2,dx,dy

$if($CUDA)
use cudafor
use cuda_fft
use cuda_padd_mod
$else
use fft
$endif

implicit none
integer::jz
real(kind=rprec),dimension(ld,ny,nz),intent(in):: u
real(kind=rprec),dimension(ld_big,ny2,nz),intent(inout)::u_big
real(kind=rprec)::ignore_me,const

$if($CUDA)
real(rprec), device, allocatable(:,:) :: u_dev, u_big_dev
$else
! ahh get rid of this
real(kind=rprec),dimension(ld,ny,nz)::temp    
$endif
! be careful using u after calling this subroutine!
const=1._rprec/(nx*ny)

$if($CUDA)

$else
temp=const*u
$endif

do jz=1,nz

  $if($CUDA) 
  !  Copy data to device
  u_dev = const*u(:,:,jz)
  !  Perform fft 
  call cufftExecD2Z(cuda_forw,u_dev,u_dev)  

  ! padd u
  call cuda_padd_zero<<< dimGrid_big, dimBlock >>>( u_big_dev )
  call cuda_padd<<< dimGrid, dimBlock >>>( u_dev, u_big_dev )

  call cufftExecZ2D(cuda_back_big,u_big_dev,u_big_dev) 

  !  Copy back to host
  u_big(:,:,jz) = u_big_dev  

  $else

  ! still need to normalize
  call rfftwnd_f77_one_real_to_complex(forw,temp(:,:,jz),ignore_me)
  ! check padd syntax
  call padd(u_big(:,:,jz),temp(:,:,jz))
  call rfftwnd_f77_one_complex_to_real(back_big,u_big(:,:,jz),ignore_me)

  $endif

end do

$if($CUDA)
deallocate(u_dev, u_big_dev)
$endif

return
end subroutine dealias1
