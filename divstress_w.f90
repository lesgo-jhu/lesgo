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

!--provides divt 1:nz
!--nz may not be used anyway (BC is used instead)
!--MPI: provides 1:nz-1, except at top 1:nz
subroutine divstress_w(divt, tx, ty, tz)
use types, only : rprec
use param, only : ld, nx, ny, nz, coord, BOGUS, lbz
use param, only : ubc_mom   !! channel
use derivatives, only : ddx, ddy, ddz_uv
#ifdef PPMPI
use param,only:nproc,coord
#endif
implicit none

real(kind=rprec),dimension(ld,ny,lbz:nz),intent(out)::divt
real (rprec), dimension (ld, ny, lbz:nz), intent (in) :: tx, ty, tz
real(kind=rprec),dimension(ld,ny,lbz:nz)::dtxdx,dtydy, dtzdz
integer::jx,jy,jz

#ifdef PPVERBOSE
write (*, *) 'started divstress_w'
#endif

! compute stress gradients      
!--tx 1:nz => dtxdx 1:nz
call ddx(tx, dtxdx, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
  dtxdx(:, :, 0) = BOGUS
#endif
#endif


!--ty 1:nz => dtydy 1:nz
call ddy(ty, dtydy, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
  dtydy(:, :, 0) = BOGUS
#endif
#endif

!--tz 0:nz-1 (special case) => dtzdz 1:nz-1 (default), 2:nz-1 (bottom),
!                                    1:nz (top)
call ddz_uv(tz, dtzdz, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
  dtzdz(:, :, 0) = BOGUS
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef PPMPI
  divt(:, :, 0) = BOGUS
#endif
#endif

! TODO does this give the right behavior for slip wall at bottom
! or is it necessary to condition on lbc_mom like done below?
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

!!channel
! TODO not sure why this is necessary, since it is not necessary for lbc_mom
if (ubc_mom > 0 .and. coord == nproc-1) then  !--jb ( ubc_mom > 1 means wall )
  ! at wall we have to assume that dz(tzz)=0.0.  Any better ideas?
  do jy=1,ny
  do jx=1,nx
  ! in old version, it looks like some people tried some stuff with dwdz here
  ! but they were zeroed out, so they were left out of this version
     divt(jx,jy,nz)=dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)
  end do
  end do
else
  do jy=1,ny
  do jx=1,nx              
     divt(jx,jy,nz)=dtxdx(jx,jy,nz) + dtydy(jx,jy,nz) + dtzdz(jx,jy,nz)
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

! TODO: is this intentional? the if-logic is overlapping
if ( ubc_mom .ne. 1 ) then
#ifdef PPMPI 
  if (coord == nproc-1) then
    do jy=1,ny
    do jx=1,nx              
       divt(jx,jy,nz)=dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)+dtzdz(jx,jy,nz)
    end do
    end do
    divt(ld-1:ld, :, nz) = 0._rprec
  else
#ifdef PPSAFETYMODE
    divt(:, :, nz) = BOGUS
#endif    
  endif
#else
  do jy=1,ny
  do jx=1,nx              
     divt(jx,jy,nz)=dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)+dtzdz(jx,jy,nz)
  end do
  end do
  divt(ld-1:ld, :, nz) = 0._rprec
#endif
endif

#ifdef PPVERBOSE
write (*, *) 'finished divstress_w'
#endif

end subroutine divstress_w
