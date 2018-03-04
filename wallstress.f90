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
subroutine wallstress ()
use types,only:rprec
use param,only:dz,ld,lh,nx,ny,nz,vonk,lbc_mom,zo,u_b, jt_total, &
               coord, MPI_RPREC, comm, ierr, domain_nstart, &
               domain_nend, domain_nskip
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use sim_param, only: ustar, ustar_inst_avg
use test_filtermodule
use ocean_base, only: ocean_flag
use mpi_defs

implicit none
integer::jx,jy
!real(kind=rprec),dimension(nx,ny)::ustar,u_avg,denom
real(kind=rprec),dimension(nx,ny)::u_avg,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec),dimension(ld,ny)::u_r1,v_r1
real(kind=rprec)::const
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m

! Ocean simulation -- Eshwan
!if (ocean_flag) then
!   txz(:,:,1) = 1._rprec ! Everything has been normalized by u_star
!   tyz(:,:,1) = 0._rprec 
!   dudz(:,:,1) = dudz(:,:,2)
!   dvdz(:,:,1) = 0._rprec
!endif
! End of section for ocean simulation


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

  case (2) !Prescribed stress - Eshwan
   txz(:,:,1) = 1.0_rprec ! Everything has been normalized by u_star
   tyz(:,:,1) = 0._rprec 
   dudz(:,:,1) = dudz(:,:,2)
   dvdz(:,:,1) = 0._rprec

  case (3) !Prescribed velocity - Eshwan
   
    psi_m=0._rprec
    phi_m=1._rprec

    u1=u(:,:,1) 
    v1=v(:,:,1)
    u_r1 = u1 - u_b
    v_r1 = v1 

    call test_filter ( u_r1 )
    call test_filter ( v_r1 )
    denom=log(0.5_rprec*dz/zo)
    u_avg=sqrt(u_r1(1:nx,1:ny)**2+v_r1(1:nx,1:ny)**2)
    ustar=u_avg*vonk/denom

    do jy=1,ny
    do jx=1,nx
       const=-(ustar(jx,jy)**2) /u_avg(jx,jy)

       txz(jx,jy,1)=const *u_r1(jx,jy)
       tyz(jx,jy,1)=const *v_r1(jx,jy)

       dudz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*u_r1(jx,jy)/u_avg(jx,jy)
       dvdz(jx,jy,1)=ustar(jx,jy)/(0.5_rprec*dz*vonK)*v_r1(jx,jy)/u_avg(jx,jy)
       dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u_r1(jx,jy).eq.0._rprec)
       dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v_r1(jx,jy).eq.0._rprec)
    end do
    end do

    !Calculate instantaneous average of ustar
    ustar_inst_avg = 0.0_rprec 

    if (coord==0) then
       do jy=1,ny
          do jx=1,nx
             ustar_inst_avg = ustar_inst_avg + ustar(jx,jy)
          end do
       end do
       ustar_inst_avg = ustar_inst_avg / real(nx*ny)
    end if

    !call mpi_bcast(ustar_inst_avg, 1, MPI_RPREC, 0, comm, ierr)

    !if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and. ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
    if (mod(jt_total,100)==0) then
    if (coord==0) write(*,*) 'Instantaneous average friction velocity', jt_total, ustar_inst_avg
    end if    

  case default

    write (*, *) 'invalid lbc_mom'
    stop

end select

end subroutine wallstress
