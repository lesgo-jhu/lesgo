subroutine unpadd(cc,cc_big)
use types,only:rprec
use param,only:ld,nx,ny,ny2,nz,ld_big,lh
implicit none
integer::jx,jy
!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( ld, ny ) :: cc
real(rprec), dimension( ld_big, ny2 ) :: cc_big

do jy=1,ny/2
  do jx=1,nx
   cc(jx,jy) = cc_big(jx,jy)
  end do
end do

! oddballs
cc(ld-1:ld,:) = 0._rprec
cc(:,ny/2+1) = 0._rprec

do jy=1,ny/2-1
  do jx=1,nx
    cc(jx,jy+ny/2+1) = cc_big(jx,jy+ny2-ny/2+1)
  end do
end do

end subroutine unpadd
