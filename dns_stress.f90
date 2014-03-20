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

subroutine dns_stress()   !txx,txy,txz,tyy,tyz,tzz)   !--jb
  ! using the 'sgs' sign convention for stress, so there is a - sign
  use types,only:rprec
  use param,only:ld,nx,ny,nz,z_i,u_star,nu_molec,coord,nproc,BOGUS,channel_bc
  use sim_param,only:dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  use sim_param,only:txx,txy,txz,tyy,tyz,tzz   !--jb

  $if ($MPI)
  use mpi_defs, only: mpi_sync_real_array, MPI_SYNC_DOWN
  $endif

  implicit none
  !!real(rprec),dimension(ld,ny,nz),intent(out)::txx,txy,txz,tyy, tyz,tzz   !--jb
  real(rprec)::S11,S12,S13,S22,S23,S33
  real(rprec)::nu
  integer::jx,jy,jz,jz_min,jz_max   !--jb

  ! non-dimensional molecular viscosity
  nu=nu_molec/(z_i*u_star)

  $if ($MPI)    !--jb, copied from sgs_stag_util.f90
    ! dudz calculated for 0:nz-1 (on w-nodes) except bottom process
    ! (only 1:nz-1) exchange information between processors to set
    ! values at nz from jz=1 above to jz=nz below
    call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
  $endif

  ! uvp-nodes
  do jz=1, nz-1
     do jy=1,ny
        do jx=1,nx
           S11=dudx(jx,jy,jz)
           S12=0.5_rprec*(dudy(jx,jy,jz)+dvdx(jx,jy,jz))
           S22=dvdy(jx,jy,jz)
           S33=dwdz(jx,jy,jz)
           txx(jx,jy,jz)=-2._rprec*nu*S11
           txy(jx,jy,jz)=-2._rprec*nu*S12
           tyy(jx,jy,jz)=-2._rprec*nu*S22
           tzz(jx,jy,jz)=-2._rprec*nu*S33
        end do
     end do
  end do

!!$if (channel_bc) then
!!$   $if ($MPI)
!!$   if (coord == nproc-1) then
!!$   $endif
!!$     ! uvp-nodes
!!$     do jy=1,ny
!!$     do jx=1,nx
!!$        S11=dudx(jx,jy,nz)
!!$        S12=0.5_rprec*(dudy(jx,jy,nz)+dvdx(jx,jy,nz))
!!$        S22=dvdy(jx,jy,nz)
!!$        S33=dwdz(jx,jy,nz)
!!$        txx(jx,jy,nz)=-2._rprec*nu*S11
!!$        txy(jx,jy,nz)=-2._rprec*nu*S12
!!$        tyy(jx,jy,nz)=-2._rprec*nu*S22
!!$        tzz(jx,jy,nz)=-2._rprec*nu*S33
!!$     end do
!!$     end do
!!$   $if ($MPI)
!!$   endif
!!$   $endif
!!$else
  !--if values not needed, set to bogus value (easier to catch errors)
  $if ($MPI)
  if (coord == nproc-1) then
  $endif
     ! top values of txx, txy, tyy, tzz not needed for stress free bc's
     !--if values not needed, set to bogus value (easier to catch errors)
     txx(:,:,nz) = BOGUS
     txy(:,:,nz) = BOGUS
     tyy(:,:,nz) = BOGUS
     tzz(:,:,nz) = BOGUS
  $if ($MPI)
  endif
  $endif

!!$endif

  ! w-nodes
  if (coord == 0) then
     ! leave the wall level alone: taken care of with wall stress
     !--assume here that wall stress has already been set (for MPI)
     jz_min = 2
  else
     jz_min = 1
  end if

  if (coord == nproc-1) then  !--jb
     jz_max = nz-1
  else
     jz_max = nz
  endif

  do jz = jz_min, jz_max    !nz-1
     do jy = 1, ny
        do jx = 1, nx
           S13 = 0.5_rprec*(dudz(jx,jy,jz)+dwdx(jx,jy,jz))
           S23 = 0.5_rprec*(dvdz(jx,jy,jz)+dwdy(jx,jy,jz))
           txz(jx,jy,jz) = -2._rprec*nu*S13
           tyz(jx,jy,jz) = -2._rprec*nu*S23
        end do
     end do
  end do

  $if ($MPI) 
  if (coord == nproc-1 .and. .not. channel_bc) then   ! stress-free lid
     txz(:,:,nz)=0._rprec
     tyz(:,:,nz)=0._rprec
  endif
!!$  else
!!$     !--nz here saves communication in MPI version: can only do this since
!!$     !  dudz, dwdx, dvdz, dwdy are available at nz (not true w/ all components)
!!$     do jy = 1, ny
!!$        do jx = 1, nx
!!$           S13 = 0.5_rprec*(dudz(jx,jy,nz)+dwdx(jx,jy,nz))
!!$           S23 = 0.5_rprec*(dvdz(jx,jy,nz)+dwdy(jx,jy,nz))
!!$           txz(jx,jy,nz) = -2._rprec*nu*S13
!!$           tyz(jx,jy,nz) = -2._rprec*nu*S23
!!$        enddo
!!$     enddo
!!$  endif
  $else
     if (.not. channel_bc) then    ! stress-free lid
        txz(:,:,nz)=0._rprec
        tyz(:,:,nz)=0._rprec
     endif
  $endif

$if ($MPI)    !--jb
  call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
  call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
$endif

end subroutine dns_stress
