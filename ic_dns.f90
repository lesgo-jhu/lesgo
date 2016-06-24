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

subroutine ic_dns()
use types,only:rprec
use param
use sim_param,only:u,v,w
implicit none
real(rprec),dimension(nz)::ubar
real(rprec)::rms,temp
integer::jx,jy,jz,z
real(rprec) :: dummy_rand

if (inflow) then

  ! uniform flow case:
  ubar = inflow_velocity

else

  ! calculate height of first uvp point in wall units
  ! lets do a laminar case (?)
  do jz=1,nz
  
     z=int((real(jz)-.5_rprec)*dz) ! non-dimensional
     ubar(jz)=(u_star*z_i/nu_molec)*z*(1._rprec-.5_rprec*z) ! non-dimensional
  !         ubar(jz)=0.
  end do  
end if

call init_random_seed

do jz=1,nz
  print *,'jz, ubar:',jz,ubar(jz)
end do
! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms=0.2_rprec
do jz=1,nz
  do jy=1,ny
     do jx=1,nx
       call random_number(dummy_rand)
       u(jx,jy,jz)=ubar(jz)+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
       call random_number(dummy_rand)
       v(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
       call random_number(dummy_rand)
       w(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
    end do
  end do
end do

! make sure w-mean is 0
temp=0._rprec
do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         temp=temp+w(jx,jy,jz)
      end do
   end do
end do
temp=temp/(nx*ny*nz)

do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         w(jx,jy,jz)=w(jx,jy,jz)-temp
      end do
   end do
end do
      
w(:,:,1)=0._rprec
w(:,:,nz)=0._rprec
u(:,:,nz)=u(:,:,nz-1)
v(:,:,nz)=v(:,:,nz-1)
end subroutine ic_dns
