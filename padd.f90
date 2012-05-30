subroutine padd (u_big,u)
! puts arrays into larger, zero-padded arrays 
! automatically zeroes the oddballs
use types,only:rprec
use param,only:ld,ld_big,nx,ny,ny2
implicit none
integer::jx,jy
!  u and u_big are interleaved as complex arrays
real(kind=rprec), dimension(ld,ny), intent(in) :: u
real(kind=rprec), dimension(ld_big,ny2), intent(out) :: u_big

! make sure the big array is zeroed!
do jy=1,ny2
  do jx=1,ld_big
    u_big(jx,jy) = 0._rprec
  enddo
enddo

! note: the loops are split in an attempt to maintain locality
do jy=1,ny/2
  do jx=1,nx
    u_big(jx,jy) = u(jx,jy)
  end do
end do

do jy=1,ny/2-1
  do jx=1,nx
    u_big(jx,jy+ny2-ny/2+1) = u(jx,jy+ny/2+1)
  end do
end do

return
end subroutine padd