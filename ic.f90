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

subroutine ic()
  use types,only:rprec
  use param
  use sim_param,only:u,v,w
  use messages, only : error
#ifdef PPTURBINES
  use turbines, only: turbine_vel_init
#endif 
  
  implicit none
  
#ifdef PPVERBOSE
  character(*), parameter :: sub_name = 'ic'
#endif

#ifdef PPTURBINES
    real(rprec) :: zo_turbines
#endif    
  
  integer::jx,jy,jz,seed
  integer :: jz_abs

  real(kind=rprec),dimension(nz)::ubar
  real(kind=rprec)::rms, noise, arg, arg2
  real(kind=rprec)::z,w_star

#ifdef PPTURBINES
  zo_turbines = 0._rprec
#endif

#ifdef PPCPS  

  call boundary_layer_ic()

#else
  
  if ( inflow ) then  !--no turbulence
     call uniform_ic()
  else
     call boundary_layer_ic()
  end if

#endif

  !VK Display the mean vertical profiles of the initialized variables on the
  !screen
  do jz=1,nz
#ifdef PPMPI
     z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#else
     z = (jz - 0.5_rprec) * dz
#endif
     write(6,7780) jz,z,sum(u(1:nx,:,jz))/float(nx*ny),sum(v(1:nx,:,jz))/&
          float(nx*ny),sum(w(1:nx,:,jz))/float(nx*ny)
  end do
7780 format('jz,z,ubar,vbar,wbar:',(1x,I3,4(1x,F9.4)))

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine uniform_ic()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

u = inflow_velocity 
!v = 0.05_rprec * inflow_velocity
v = 0._rprec
w = 0._rprec

return
end subroutine uniform_ic

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine boundary_layer_ic()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Log profile that is modified to flatten at z=z_i
! This is a hot mess (JSG 20111221)
implicit none

interface
   function ran3(idum)
     integer(4) idum
     real(8) :: ran3
   end function ran3
end interface

!w_star=(9.81_rprec/T_init*wt_s*z_i)**(1._rprec/3._rprec)
w_star = u_star


!      T_star=wt_s/w_star
!      q_star=T_star

!#ifdef PPMPI
!  seed = -80 - coord
!#else
!  seed=-80
!#endif

if( coord == 0 ) write(*,*) '------> Creating modified log profile for IC'
do jz=1,nz

#ifdef PPMPI
   z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#else
   z=(jz-.5_rprec)*dz
#endif

   ! IC in equilibrium with rough surface (rough dominates in effective zo)
   arg2=z/zo
   arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

#ifdef PPLVLSET
   ! Kludge to adjust magnitude of velocity profile
   ! Not critical - may delete
   arg = 0.357*arg
#endif

#ifdef PPTURBINES
   call turbine_vel_init (zo_turbines)
   arg2=z/zo_turbines
   arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z          
#endif        

   !ubar(jz)=arg

   ! Added by VK for making the u less than 1...need to change this
   ! initialization routine
   if (coriolis_forcing) then
      ubar(jz)=arg/30._rprec
   else
      ubar(jz)=arg

   end if

   if ((coriolis_forcing).and.(z.gt.(.5_rprec))) ubar(jz)=ug

end do

rms = 3._rprec
do jz=1,nz
#ifdef PPMPI
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i    !dimensions in meters
#else
   jz_abs = jz
   z = (jz-.5_rprec) * dz * z_i                        !dimensions in meters
#endif
   seed = -80 - jz_abs  !--trying to make consistent init for MPI
   do jy=1,ny
      do jx=1,nx
         !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
         !...Taking std dev of vel as 1 at all heights
         if (z.le.z_i) then
            noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
            u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
            noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
            v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star !noise
            noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
            w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star
! Tony ATM
!~ u(jx,jy,jz)=1._rprec
!~ v(jx,jy,jz)=0._rprec
!~ w(jx,jy,jz)=0._rprec

         else
! Tony ATM
!~ u(jx,jy,jz)=1._rprec
!~ v(jx,jy,jz)=0._rprec
!~ w(jx,jy,jz)=0._rprec
            noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
            u(jx,jy,jz)=noise*w_star/u_star*.01_rprec+ubar(jz)
            noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
            v(jx,jy,jz)=noise*w_star/u_star*.01_rprec
            noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
            w(jx,jy,jz)=noise*w_star/u_star*.01_rprec
         end if
      end do
   end do
end do

!...BC for W
if (coord == 0) then
   w(1:nx, 1:ny, 1) = 0._rprec
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif    
   w(1:nx, 1:ny, nz) = 0._rprec
   u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
   v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
#ifdef PPMPI
end if
#endif

return
end subroutine boundary_layer_ic

end subroutine ic
