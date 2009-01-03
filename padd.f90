subroutine padd (u_big,u)
! puts arrays into larger, zero-padded arrays 
! automatically zeroes the oddballs
use types,only:rprec
use param,only:lh,lh_big,nx,ny,ny2
implicit none
integer::jx,jy
! note we're calling with 2D arrays
complex(kind=rprec),dimension(lh,ny),intent(in)::u
complex(kind=rprec),dimension(lh_big,ny2),intent(out)::u_big
! make sure the big array is zeroed!
u_big=0._rprec
! note: the loops are split in an attempt to maintain locality
! test it                   !
do jy=1,ny/2
do jx=1,nx/2
! skip the Nyquist frequency since it should be zero anyway
   u_big(jx,jy)=u(jx,jy)
end do
end do
do jy=1,ny/2-1
do jx=1,nx/2
   u_big(jx,jy+ny2-ny/2+1)=u(jx,jy+ny/2+1)
end do
end do
end subroutine padd
