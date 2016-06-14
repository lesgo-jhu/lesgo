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

! For use with staggered grid LES
! JDA, 23 Jan 96
!--provides txz, tyz (w-nodes) and dudz, dvdz (w-nodes) at jz=1

!xiang: when this function is called, it already ensures coord=0
subroutine wallstress ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,lbc_mom,zo
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use test_filtermodule
use iwmles, only : iwm_on, iwm_wallstress !xiang for the option of integral wall model
implicit none
integer::jx,jy
real(rprec),dimension(nx,ny)::denom,u_avg,ustar
real(rprec),dimension(ld,ny)::u1,v1
! No need to define phi_m or psi_m as a matrix as only a constant value is used
real(rprec)::const,phi_m,psi_m

select case (lbc_mom)
  case (0) ! Stress free
    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec
  case (1) ! Wall
   if(iwm_on == 0)then !if not using integral wall model...
    ! See John D. Albertson's dissertation, eqns (2.46)-(2.52)
    ! For dudz and dvdz at wall, we should use derivwall.f90
    ! Also, see:
    ! E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent Lagrangian dynamic model
    !   for large eddy simulation of complex turbulent flows" (2005) -- Appendix    
    
    !TS Remove the following line when obukhov.f is used
    psi_m=0._rprec
    phi_m=1._rprec

    u1=u(:,:,1)
    v1=v(:,:,1)
    call test_filter ( u1 )
    call test_filter ( v1 )
    denom=log(0.5_rprec*dz/zo)-psi_m
    u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
    ustar=u_avg*vonk/denom

    do jy=1,ny
    do jx=1,nx
       const=-(ustar(jx,jy)**2) /u_avg(jx,jy)
       txz(jx,jy,1)=const *u1(jx,jy)
       tyz(jx,jy,1)=const *v1(jx,jy)
    !TS REMOVE derivwall.f90 and add it here
    !this is as in Moeng 84
       dudz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,1)/u_avg(jx,jy)&
    !TS ADD for non-neutral case
           *phi_m
       dvdz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
    !TS ADD for non-neutral case
           *phi_m
       dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
       dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
    end do
    end do
   else
	call iwm_wallstress() !xiang: calculate stress using the integral wall model...
   endif
end select
end subroutine wallstress