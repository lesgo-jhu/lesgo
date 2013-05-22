!!
!!  Copyright 2009,2011 Johns Hopkins University
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

subroutine ic_dns()
use types,only:rprec
use param
use sim_param,only:u,v,w
implicit none
real(kind=rprec),dimension(nz)::ubar
real(kind=rprec)::rms,temp,ran3
integer::jx,jy,jz,seed,z

if (inflow) then

  ! uniform flow case:
  ubar = inflow_velocity

else

  seed=-112

  ! calculate height of first uvp point in wall units
  ! lets do a laminar case (?)
  do jz=1,nz
  
     z=(real(jz)-.5_rprec)*dz ! non-dimensional
     ubar(jz)=(u_star*z_i/nu_molec)*z*(1._rprec-.5_rprec*z) ! non-dimensional
  !         ubar(jz)=0.
  end do  
end if

do jz=1,nz
  print *,'jz, ubar:',jz,ubar(jz)
end do
! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms=0.2_rprec
do jz=1,nz
  do jy=1,ny
     do jx=1,nx
       u(jx,jy,jz)=ubar(jz)+(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
       v(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
       w(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
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
