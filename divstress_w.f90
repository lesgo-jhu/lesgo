!!
!!  Copyright 2009,2010,2011,2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

!--provides divt 1:nz
!--nz may not be used anyway (BC is used instead)
!--MPI: provides 1:nz-1, except at top 1:nz
subroutine divstress_w(divt, tx, ty, tz)
use types,only:rprec
use param,only:ld,nx,ny,nz, nproc, coord, BOGUS, lbz
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
$if ($SAFETYMODE)
$if ($MPI)
  dtxdx(:, :, 0) = BOGUS
$endif
$endif


!--ty 1:nz => dtydy 1:nz
call ddy(ty, dtydy, lbz)
$if ($SAFETYMODE)
$if ($MPI)
  dtydy(:, :, 0) = BOGUS
$endif
$endif

!--tz 0:nz-1 (special case) => dtzdz 1:nz-1 (default), 2:nz-1 (bottom),
!                                    1:nz (top)
call ddz_uv(tz, dtzdz, lbz)
$if ($SAFETYMODE)
$if ($MPI)
  dtzdz(:, :, 0) = BOGUS
$endif
$endif

$if ($SAFETYMODE)
$if ($MPI)
  divt(:, :, 0) = BOGUS
$endif
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

$if ($MPI) 
  if (coord == nproc-1) then
    do jy=1,ny
    do jx=1,nx              
       divt(jx,jy,nz)=dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)+dtzdz(jx,jy,nz)
    end do
    end do
    divt(ld-1:ld, :, nz) = 0._rprec
  else
$if ($SAFETYMODE)
    divt(:, :, nz) = BOGUS
$endif    
  endif
$else
  do jy=1,ny
  do jx=1,nx              
     divt(jx,jy,nz)=dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)+dtzdz(jx,jy,nz)
  end do
  end do
  divt(ld-1:ld, :, nz) = 0._rprec
$endif

$if ($VERBOSE)
write (*, *) 'finished divstress_w'
$endif

end subroutine divstress_w
