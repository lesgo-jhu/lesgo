subroutine padd (u_big,u)
! puts arrays into larger, zero-padded arrays 
! automatically zeroes the oddballs
use types,only:rprec
use param,only:ld,ld_big,nx,ny,ny2
implicit none
integer::jx,jy
integer :: ir, ii
integer :: jy_off, jy_off_big

! note we're calling with 2D arrays
!complex(kind=rprec),dimension(lh,ny),intent(in)::u
!complex(kind=rprec),dimension(lh_big,ny2),intent(out)::u_big
real(kind=rprec), dimension(ld,ny), intent(in) :: u
real(kind=rprec), dimension(ld_big,ny2), intent(out) :: u_big

! make sure the big array is zeroed!
u_big=0._rprec
! note: the loops are split in an attempt to maintain locality
! test it                   !
do jy=1,ny/2
  do jx=1,nx/2

    ii = 2*jx
    ir = ii - 1

    ! skip the Nyquist frequency since it should be zero anyway
    !u_big(jx,jy)=u(jx,jy)
    u_big(ir:ii,jy) = u(ir:ii,jy)

  end do
end do

do jy=1,ny/2-1

  !  Cache index
  jy_off_big = jy+ny2-ny/2+1
  jy_off     = jy+ny/2+1

  do jx=1,nx/2
  
    ii = 2*jx
    ir = ii - 1
  
    !u_big(jx,jy+ny2-ny/2+1)=u(jx,jy+ny/2+1)
    u_big( ir:ii, jy_off_big ) = u( ir:ii, jy_off )

  end do
end do

return
end subroutine padd
