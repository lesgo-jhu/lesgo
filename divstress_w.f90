!--provides divt 1:nz
!--nz may not be used anyway (BC is used instead)
!--MPI: provides 1:nz-1, except at top 1:nz
subroutine divstress_w(divt, tx, ty, tz)
use types,only:rprec
use param,only:ld,nx,ny,nz, USE_MPI, nproc, coord, BOGUS, lbz
use derivatives, only : ddx, ddy, ddz_uv

implicit none

real(kind=rprec),dimension(ld,ny,lbz:nz),intent(out)::divt
real (rprec), dimension (ld, ny, lbz:nz), intent (in) :: tx, ty, tz
real(kind=rprec),dimension(ld,ny,lbz:nz)::dtxdx,dtydy, dtzdz
integer::jx,jy,jz

$if ($VERBOSE)
write (*, *) 'started divstress_w'
$endif

! compute stress gradients      
!--tx 1:nz => dtxdx 1:nz
call ddx(tx, dtxdx, lbz)
$if ($MPI)
  dtxdx(:, :, 0) = BOGUS
$endif

!--ty 1:nz => dtydy 1:nz
call ddy(ty, dtydy, lbz)
$if ($MPI)
  dtydy(:, :, 0) = BOGUS
$endif

!--tz 0:nz-1 (special case) => dtzdz 1:nz-1 (default), 2:nz-1 (bottom),
!                                    1:nz (top)
call ddz_uv(tz, dtzdz, lbz)
$if ($MPI)
  dtzdz(:, :, 0) = BOGUS
$endif

$if ($MPI)
  divt(:, :, 0) = BOGUS
$endif

if (coord == 0) then
  ! at wall we have to assume that dz(tzz)=0.0.  Any better ideas?
  do jy=1,ny
  do jx=1,nx
  ! in old version, it looks like some people tried some stuff with dwdz here
  ! but they were zeroed out, so they were left out of this version
     divt(jx,jy,1)=dtxdx(jx,jy,1)+dtydy(jx,jy,1)
  end do
  end do
else
  do jy=1,ny
  do jx=1,nx              
     divt(jx,jy,1)=dtxdx(jx,jy,1)+dtydy(jx,jy,1)+dtzdz(jx,jy,1)
  end do
  end do
end if

do jz=2,nz-1
do jy=1,ny
do jx=1,nx              
   divt(jx,jy,jz)=dtxdx(jx,jy,jz)+dtydy(jx,jy,jz)+dtzdz(jx,jy,jz)
end do
end do
end do

!--set ld-1, ld to 0 (could maybe do BOGUS)
divt(ld-1:ld, :, 1:nz-1) = 0._rprec

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  do jy=1,ny
  do jx=1,nx              
     divt(jx,jy,nz)=dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)+dtzdz(jx,jy,nz)
  end do
  end do

  divt(ld-1:ld, :, nz) = 0._rprec
else
  divt(:, :, nz) = BOGUS
end if

$if ($VERBOSE)
write (*, *) 'finished divstress_w'
$endif

end subroutine divstress_w
