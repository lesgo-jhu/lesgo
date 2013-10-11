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
use param,only:nx,ny,nu_molec,z_i,u_star,dz,lbc_mom
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
    do jy=1,ny
    do jx=1,nx
       txz(jx,jy,1)=-nu_molec/(z_i*u_star)*u(jx,jy,1)/(0.5_rprec*dz)
       tyz(jx,jy,1)=-nu_molec/(z_i*u_star)*v(jx,jy,1)/(0.5_rprec*dz)
       dudz(jx,jy,1)=u(jx,jy,1)/(0.5_rprec*dz)
       dvdz(jx,jy,1)=v(jx,jy,1)/(0.5_rprec*dz)
    end do
    end do

  case default

    write (*, *) 'invalid lbc_mom'

end select

end subroutine wallstress_dns
