!**********************************************************************
subroutine ddz_uv(f, dfdz)
!**********************************************************************
!
! first derivative in z direction for boundary layer (2nd order numerics)
!--provides dfdz(:, :, 2:nz), value at jz=1 is not touched
!--MPI: provides dfdz(:, :, 1:nz) where f is on uvp-node, except at
!  bottom process it only supplies 2:nz
!
use types,only:rprec
use param,only:ld,nx,ny,nz,dz, USE_MPI, coord, nproc, BOGUS
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real (rprec), dimension (ld, ny, $lbz:nz), intent (in) :: f
real (rprec), dimension(ld, ny, $lbz:nz), intent (inout) :: dfdz

integer::jx,jy,jz
real (rprec) :: const

const=1._rprec/dz

$if ($MPI)

  dfdz(:, :, 0) = BOGUS

  if (coord > 0) then
    !--ghost node information is available here
    !--watch the z-dimensioning!
    !--if coord == 0, dudz(1) will be set in wallstress
    do jy=1,ny
    do jx=1,nx    
       dfdz(jx,jy,1)=const*(f(jx,jy,1)-f(jx,jy,0))
    end do
    end do
  end if

$endif

do jz=2,nz-1
do jy=1,ny
do jx=1,nx    
   dfdz(jx,jy,jz)=const*(f(jx,jy,jz)-f(jx,jy,jz-1))
end do
end do
end do

!--should integrate this into the above loop, explicit zeroing is not
!  needed, since dudz, dvdz are forced to zero by copying the u, v fields
!--also, what happens when called with tzz? 
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  dfdz(:,:,nz)=0._rprec  !--do not need to do this...
else
  do jy=1,ny
  do jx=1,nx
     dfdz(jx,jy,nz)=const*(f(jx,jy,nz)-f(jx,jy,nz-1))
  end do
  end do
end if

end subroutine ddz_uv
