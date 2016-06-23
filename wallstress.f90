!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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
!**********************************************************************
subroutine wallstress
!**********************************************************************
! 
use types, only : rprec
use param, only : lbc_mom
use messages, only : error
use iwmles, only : iwm_wallstress
use sim_param, only : txz, tyz, dudz, dvdz
implicit none
character(*), parameter :: sub_name = 'wallstress'

select case (lbc_mom)
    case (0)                        ! Stress free
        call ws_free

    case (1)                        ! DNS wall
        call ws_dns

    case (2)                        ! Equilibrium wall model
        call ws_equilibrium

    case (3)                        ! Integral wall model
        call iwm_wallstress()
    
    case default
        call error (sub_name, 'invalid lbc_mom')
        
end select

contains

!**********************************************************************
subroutine ws_free
!**********************************************************************
! 
implicit none

txz(:, :, 1) = 0._rprec
tyz(:, :, 1) = 0._rprec
dudz(:, :, 1) = 0._rprec
dvdz(:, :, 1) = 0._rprec

end subroutine ws_free

!**********************************************************************
subroutine ws_dns
!**********************************************************************
! 
use param, only : ld, nx, ny, nz, nu_molec, z_i, u_star, dz
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
       txz(i,j,1) = -nu_molec/(z_i*u_star)*u(i,j,1)/(0.5_rprec*dz)
       tyz(i,j,1) = -nu_molec/(z_i*u_star)*v(i,j,1)/(0.5_rprec*dz)
       dudz(i,j,1) = u(i,j,1)/(0.5_rprec*dz)
       dvdz(i,j,1) = v(i,j,1)/(0.5_rprec*dz)
    end do
end do

end subroutine ws_dns

!**********************************************************************
subroutine ws_equilibrium
!**********************************************************************
! 
! See John D. Albertson's dissertation, eqns (2.46)-(2.52)
! For dudz and dvdz at wall, we should use derivwall.f90
! Also, see:
! E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent Lagrangian dynamic model
!   for large eddy simulation of complex turbulent flows" (2005) -- Appendix    
! 
use param,only:dz,ld,lh,nx,ny,nz,vonk,zo
use messages, only : error
use sim_param,only:u,v
use test_filtermodule
implicit none
integer::i,j
real(rprec),dimension(nx,ny)::denom,u_avg,ustar
real(rprec),dimension(ld,ny)::u1,v1
! No need to define phi_m or psi_m as a matrix as only a constant value is used
real(rprec)::const,phi_m,psi_m

psi_m=0._rprec
phi_m=1._rprec

u1=u(:,:,1)
v1=v(:,:,1)
call test_filter ( u1 )
call test_filter ( v1 )
denom=log(0.5_rprec*dz/zo)-psi_m
u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar=u_avg*vonk/denom

do j=1,ny
do i=1,nx
   const=-(ustar(i,j)**2) /u_avg(i,j)
   txz(i,j,1)=const *u1(i,j)
   tyz(i,j,1)=const *v1(i,j)
!TS REMOVE derivwall.f90 and add it here
!this is as in Moeng 84
   dudz(i,j,1)=ustar(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)&
!TS ADD for non-neutral case
       *phi_m
   dvdz(i,j,1)=ustar(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)&
!TS ADD for non-neutral case
       *phi_m
   dudz(i,j,1)=merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
   dvdz(i,j,1)=merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
end do
end do

end subroutine ws_equilibrium

end subroutine wallstress