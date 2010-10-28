!**********************************************************************
subroutine ddz_w(f, dfdz)
!**********************************************************************
!
!  First deriv in z direction for boundary layer (2nd order numerics)
!  F is on w nodes and dFdz is on uvp nodes
!    -23 January 1996
!--MPI: provides 0:nz-1, except top process has 0:nz, and bottom process
!  has 1:nz-1
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
real(kind=rprec),dimension(ld, ny, $lbz:nz), intent (out) :: dfdz

real(kind=rprec)::const
integer::jx,jy,jz

const=1._rprec/dz
do jz=$lbz,nz-1
do jy=1,ny
do jx=1,nx
   dfdz(jx,jy,jz)=const*(f(jx,jy,jz+1)-f(jx,jy,jz))
end do
end do
end do

if (USE_MPI .and. coord == 0) then
  !--bottom process cannot calculate dfdz(jz=0)
  dfdz(:, :, $lbz) = BOGUS
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  !--this is not needed
  ! Stress free lid (ubc)
  !if (ubc==0) then stop

  !else
  !--do not need to do this
  dfdz(:,:,nz)=0._rprec !dfdz(:,:,Nz-1) ! c? any better ideas for sponge?
  !end if
else
  !--does not supply dfdz(:, :, nz)
  dfdz(:, :, nz) = BOGUS
end if

end subroutine ddz_w
