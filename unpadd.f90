subroutine unpadd(cc,cc_big)
use types,only:rprec
use param,only:ld,nx,ny,ny2,nz,lh_big,lh
implicit none
integer::jx,jy
! note using 2d definitions here!
complex(kind=rprec),dimension(lh_big,ny2)::cc_big
complex(kind=rprec),dimension(lh,ny)::cc
do jy=1,ny/2
do jx=1,nx/2
   cc(jx,jy)=cc_big(jx,jy)
end do
end do
! oddballs
cc(nx/2+1,:)=0._rprec
cc(:,ny/2+1)=0._rprec
do jy=1,ny/2-1
do jx=1,nx/2
   cc(jx,jy+ny/2+1)=cc_big(jx,jy+ny2-ny/2+1)
end do
end do
end subroutine unpadd
