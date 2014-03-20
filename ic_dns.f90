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
real(rprec)::rms,temp,ran3,z,vpert,wpert
integer::jx,jy,jz,seed

if (inflow) then

  ! uniform flow case:
  ubar = inflow_velocity

else

  seed=-112

  !! calculate height of first uvp point in wall units
  !! lets do a laminar case (?)
  if (ic_couette) then          !! linear laminar profile (couette)
     do jz=1,nz
        $if ($MPI)    !--jb
        z=(coord*(nz-1) + real(jz) - 0.5_rprec) * dz ! non-dimensional
        $else
        z = (real(jz) - 0.5_rprec) * dz ! non-dimensional
        $endif
        ubar(jz)= (utop-ubot)/L_z * z + ubot ! non-dimensional
     end do
  else
     do jz=1,nz             !! parabolic laminar profile (channel)
        $if ($MPI)    !--jb
        z=(coord*(nz-1) + real(jz) - 0.5_rprec) * dz ! non-dimensional
        $else
        z = (real(jz) - 0.5_rprec) * dz ! non-dimensional
        $endif
        ubar(jz)=(u_star*z_i/nu_molec) * z * (1._rprec - 0.5_rprec*z) ! non-dimensional
        !         ubar(jz)=0.
     end do
  endif
 
end if

if (channel_bc .and. .not. ic_couette) then
   vpert = 5.0_rprec    !! works for channel flow
   wpert = 2.0_rprec
else
   vpert = 0.0_rprec    !! works for couette flow
   wpert = 0.0_rprec
endif

! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms=0.2_rprec
do jz=1,nz
  do jy=1,ny
     do jx=1,nx
       u(jx,jy,jz)=ubar(jz)     +(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
       v(jx,jy,jz)=vpert*u_star +(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
       w(jx,jy,jz)=wpert*u_star +(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
    end do 
  end do 
end do

do jz=1,nz
  print *,'coord, jz, u, v, w:', coord, jz, u(1,1,jz), v(1,1,jz), w(1,1,jz)
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
u(:,:,nz)=u(:,:,nz-1)  ! shouldn't make a difference, but probably good to remove this line and next  
v(:,:,nz)=v(:,:,nz-1)    
end subroutine ic_dns
