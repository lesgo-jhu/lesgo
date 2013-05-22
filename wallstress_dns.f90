!!
!!  Copyright 2009,2011,2012 Johns Hopkins University
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

subroutine wallstress_dns ()
use types,only:rprec
use param,only:ld,nx,ny,nz,nu_molec,z_i,u_star,dz,lbc_mom
use sim_param,only:u,v,dudz,dvdz, txz, tyz
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
