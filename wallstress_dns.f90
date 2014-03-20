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

subroutine wallstress_dns ()
use types,only:rprec
use param,only:nx,ny,nz,nu_molec,z_i,u_star,dz,lbc_mom
use param, only:coord,nproc,utop,ubot
use sim_param,only:u,v,dudz,dvdz,txz,tyz
implicit none
integer::jx,jy

select case (lbc_mom)

  case (0) ! Stress free

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case (1) ! Wall
    if (coord == 0) then     !--bottom wall
       do jy=1,ny
       do jx=1,nx
          !! one-sided difference (1st order approximation)
          dudz(jx,jy,1)=( u(jx,jy,1) - ubot )/(0.5_rprec*dz)
          dvdz(jx,jy,1)=v(jx,jy,1)/(0.5_rprec*dz)
          txz(jx,jy,1)=-nu_molec/(z_i*u_star)*dudz(jx,jy,1)
          tyz(jx,jy,1)=-nu_molec/(z_i*u_star)*dvdz(jx,jy,1)
       end do
       end do
    endif
    !! if no MPI, then only one processor, so it must apply both top
    !! and bottom BCs
    if (coord == nproc - 1) then   !--top wall
       do jy=1,ny
       do jx=1,nx
          !! one-sided difference (1st order approximation)
          dudz(jx,jy,nz)=( utop - u(jx,jy,nz-1) )/(0.5_rprec*dz)
          dvdz(jx,jy,nz)=-v(jx,jy,nz-1)/(0.5_rprec*dz)
          txz(jx,jy,nz)=-nu_molec/(z_i*u_star)*dudz(jx,jy,nz)
          tyz(jx,jy,nz)=-nu_molec/(z_i*u_star)*dvdz(jx,jy,nz)
       end do
       end do
    endif

  case default

    write (*, *) 'invalid lbc_mom'

end select

end subroutine wallstress_dns
