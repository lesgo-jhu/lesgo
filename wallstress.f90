!!
!!  Copyright 2009,2010,2011,2012 Johns Hopkins University
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

! For use with staggered grid LES
! JDA, 23 Jan 96
!--provides txz, tyz (w-nodes) and dudz, dvdz (w-nodes) at jz=1
subroutine wallstress ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,lbc_mom,zo
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use test_filtermodule
implicit none
integer::jx,jy
real(kind=rprec),dimension(nx,ny)::ustar,u_avg,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec)::const
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m

select case (lbc_mom)

  case (0) ! Stress free

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case (1) ! Wall
    ! See John D. Albertson's dissertation, eqns (2.46)-(2.52)
    ! For dudz and dvdz at wall, we should use derivwall.f90
    
    ! Also, see:
    ! E. Bou-Zeid, C. Meneveau & M.B. Parlange, “A scale-dependent Lagrangian dynamic model
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
           *phi_m(jx,jy)
       dvdz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
    !TS ADD for non-neutral case
           *phi_m(jx,jy)
       dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
       dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
    end do
    end do

  case default

    write (*, *) 'invalid lbc_mom'
    stop

end select

end subroutine wallstress
