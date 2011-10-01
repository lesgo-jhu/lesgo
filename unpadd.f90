subroutine unpadd(cc,cc_big)
use types,only:rprec
use param,only:ld,nx,ny,ny2,nz,ld_big,lh
implicit none
integer::jx,jy
! note using 2d definitions here!
!complex(kind=rprec),dimension(lh_big,ny2)::cc_big
!complex(kind=rprec),dimension(lh,ny)::cc

!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( ld, ny ) :: cc
real(rprec), dimension( ld_big, ny2 ) :: cc_big

integer :: ir, ii, jy_off, jy_off_big

do jy=1,ny/2
do jx=1,nx/2

   ii = 2*jx
   ir = ii - 1
   !cc(jx,jy)=cc_big(jx,jy)
   cc(ir:ii,jy) = cc_big(ir:ii,jy)

end do
end do

! oddballs
!cc(nx/2+1,:)=0._rprec ! equivalent to cc(lh,:)=0._rprec
!cc(:,ny/2+1)=0._rprec
cc(ld-1:ld,:) = 0._rprec
cc(:,ny/2+1) = 0._rprec

do jy=1,ny/2-1

  !  Cache index
  jy_off_big = jy+ny2-ny/2+1
  jy_off     = jy+ny/2+1

  do jx=1,nx/2

    ii = 2*jx
    ir = ii - 1

    !cc(jx,jy+ny/2+1)=cc_big(jx,jy+ny2-ny/2+1)
    cc(ir:ii,jy_off) = cc_big(ir:ii,jy_off_big)

  end do
end do
end subroutine unpadd
