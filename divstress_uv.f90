!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

subroutine divstress_uv (divtx, divty, txx, txy, txz, tyy, tyz)
use types,only:rprec
use param,only:ld,ny,nz,BOGUS,lbz,fourier,coord,nx
use derivatives, only : ddx, ddy, ddz_w, ddxy
use derivatives, only : ddx_kxspace, ddy_kxspace, wave2phys   !!jb

implicit none

real(kind=rprec),dimension(ld,ny,lbz:nz),intent(out)::divtx,divty
real (rprec), dimension (ld, ny, lbz:nz), intent (in) :: txx, txy, txz, tyy, tyz
! sc: we should be able to save some memory here!
! do this by writing a divergence subroutine--then do not store derivs 
real(kind=rprec),dimension(ld,ny,lbz:nz)::dtxdx,dtydy, dtzdz
real(kind=rprec),dimension(ld,ny,lbz:nz)::dtxdx2,dtydy2, dtzdz2

integer::jx,jy

$if ($DEBUG)
logical, parameter :: DEBUG = .true.
$endif

$if ($VERBOSE)
write (*, *) 'started divstress_uv'
$endif

! compute stress gradients      
!--MPI: tx 1:nz-1 => dtxdx 1:nz-1
if (.not. fourier) then
  call ddx(txx, dtxdx, lbz)  !--really should replace with ddxy (save an fft)
else
  call ddx_kxspace(txx, dtxdx, lbz)
endif

!--MPI: ty 1:nz-1 => dtdy 1:nz-1
if (fourier) then
  call ddy_kxspace(txy, dtydy, lbz)        !!jb
endif

!--MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(txz, dtzdz, lbz)

! compute stress gradients      
!--MPI: tx 1:nz-1 => dtxdx 1:nz-1
if (fourier) then
  call ddx_kxspace(txy, dtxdx2, lbz)     !!jb
endif

!--MPI: ty 1:nz-1 => dtdy 1:nz-1
if (.not. fourier) then
  call ddy(tyy, dtydy2, lbz)
else
  call ddy_kxspace(tyy, dtydy2, lbz)     !!jb
endif

!--MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(tyz, dtzdz2, lbz)

if (.not. fourier) then
  call ddxy(txy , dtxdx2, dtydy, lbz)     !!jb         
endif

!!$if (coord == 0) then
!!$   if (fourier) then
!!$      print*, 'heree5'
!!$      call wave2phys(dtxdx)
!!$      call wave2phys(dtydy)
!!$      call wave2phys(dtzdz)
!!$      call wave2phys(dtxdx2)
!!$      call wave2phys(dtydy2)
!!$      call wave2phys(dtzdz2)
!!$   endif
!!$   print*, 'qwe 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
!!$   do jx=1,nx
!!$      do jy=1,ny
!!$         write(*,*) jx,jy,dtxdx(jx,jy,1),dtydy(jx,jy,1),dtzdz(jx,jy,1)
!!$         !!write(*,*) dtzdz(jx,jy,0:3)
!!$      enddo
!!$   enddo
!!$   print*, ' qwe 2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
!!$   do jx=1,nx
!!$      do jy=1,ny
!!$         write(*,*) jx,jy,dtxdx2(jx,jy,1),dtydy2(jx,jy,1),dtzdz2(jx,jy,1)
!!$      enddo
!!$   enddo
!!$endif

!--MPI following comment only true at bottom process
! the following gives bad results...but it seems like i the
! issue should be taken care of somewhere
! need to correct wall level, since tz will be on uv-node there
!      dtzdz(:,:,1) = (tz(:,:,2)-tz(:,:,1))/(0.5*dz)

!--only 1:nz-1 are valid
divtx(:, :, 1:nz-1) = dtxdx(:, :, 1:nz-1) + dtydy(:, :, 1:nz-1) +  &
                     dtzdz(:, :, 1:nz-1)

!--Set ld-1, ld to 0 (or could do BOGUS)
divtx(ld-1:ld, :, 1:nz-1) = 0._rprec

$if ($SAFETYMODE)
$if ($MPI)
  divtx(:, :, 0) = BOGUS
$endif
divtx(:, :, nz) = BOGUS
$endif

!--only 1:nz-1 are valid
divty(:, :, 1:nz-1) = dtxdx2(:, :, 1:nz-1) + dtydy2(:, :, 1:nz-1) +  &
                     dtzdz2(:, :, 1:nz-1)

!--Set ld-1, ld to 0 (or could do BOGUS)
divty(ld-1:ld, :, 1:nz-1) = 0._rprec

$if ($SAFETYMODE)
$if ($MPI)
  divty(:, :, 0) = BOGUS
$endif
divty(:, :, nz) = BOGUS
$endif


$if ($VERBOSE)
write (*, *) 'finished divstress_uv'
$endif

end subroutine divstress_uv

