subroutine dealias2(u,u_big)
! puts back into small array
use types,only:rprec
use param,only:ld,ld_big,nx,ny,nz,nx2,ny2,dx,dy
$if($CUDA)
use cudafor
use cuda_defs
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
real(rprec), device, allocatable, dimension(:,:) :: u_dev, u_big_dev
$endif

$if($CUDA)
allocate(u_dev(ld,ny), u_big_dev(ld_big,ny2))
$endif

! normalize
const=1._rprec/(nx2*ny2)
u_big=const*u_big
! Loop through horizontal slices
do jz=1,nz

  $if($CUDA)
  !  Copy data to device
  u_big_dev = u_big(:,:,jz)

  ! perform forward FFT
  call cufftExecD2Z_2D(cuda_forw_big,u_big_dev,u_big_dev) 

  call cuda_unpadd_zero<<< dimGrid, dimBlock >>>( u_dev )
  call cuda_unpadd <<< dimGrid, dimBlock >>>( u_big_dev, u_dev )

  ! Back to physical space
   call cufftExecZ2D_2D( cuda_back, u_dev, u_dev )

  !  Copy data to host
  u(:,:,jz) = u_dev

  $else

  ! perform forward FFT
  call rfftwnd_f77_one_real_to_complex(forw_big,u_big(:,:,jz),ignore_me)    
  call unpadd(u(:,:,jz),u_big(:,:,jz))
  ! Back to physical space
  call rfftwnd_f77_one_complex_to_real(back,u(:,:,jz),ignore_me)

  $endif

end do

! sc: do we need this?
!.....Making filter circular !!!!!!!!!!!
!	call filt_da(uu)
!.......................................

$if($CUDA)
deallocate(u_dev, u_big_dev)
$endif

return
end subroutine dealias2
