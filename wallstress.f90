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
use param,only:dz,nx,ny,vonk,lbc_mom,zo
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use derivatives, only: dft_direct_back_2d_n_yonlyC_big, dft_direct_forw_2d_n_yonlyC_big
use derivatives, only: convolve
use test_filtermodule
implicit none
integer::jx,jy
!!real(rprec),dimension(nx,ny) :: denom, ustar, u_avg   !!jb
real(rprec),dimension(ld,ny) :: ustar, u_avg   !!jb
real(rprec),dimension(ld,ny) :: u1, v1
real(rprec),dimension(ld_big,ny2) :: u1_big, v1_big, s1_big, txz_big, tyz_big
! No need to define phi_m or psi_m as a matrix as only a constant value is used
real(rprec)::const,phi_m,psi_m
real(rprec)::u_ss, denom

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
    ! E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent Lagrangian dynamic model
    !   for large eddy simulation of complex turbulent flows" (2005) -- Appendix    
    
    !TS Remove the following line when obukhov.f is used
    psi_m=0._rprec
    phi_m=1._rprec

    u1=u(:,:,1)
    v1=v(:,:,1)
    !call test_filter ( u1 )  !!jb - turning off for now (fourier)
    !call test_filter ( v1 )  !!jb  !here
    denom=log(0.5_rprec*dz/zo)-psi_m

    if (fourier) then
       u1_big = 0._rprec
       v1_big = 0._rprec
       s1_big = 0._rprec
       call padd(u1_big(:,:), u1(:,:))
       call padd(v1_big(:,:), v1(:,:))
       call dft_direct_back_2d_n_yonlyC_big(u1_big(:,:))
       call dft_direct_back_2d_n_yonlyC_big(v1_big(:,:))
       s1_big = convolve(u1_big(:,:),u1_big(:,:)) + convolve(v1_big(:,:),v1_big(:,:))
       call dft_direct_forw_2d_n_yonlyC_big(s1_big(:,:))
       call unpadd(u_avg(:,:), s1_big(:,:))
       u_avg(1,1) = sqrt( u_avg(1,1) )
       u_avg(2:ld,1:ny) = 0._rprec
       u_avg(1,2:ny) = 0._rprec
    else
       u_ss = 0._rprec
       do jy=1,ny
       do jx=1,nx
          u_ss = u_ss + u1(jx,jy)**2 + v1(jx,jy)**2
       enddo
       enddo
       u_ss = sqrt( u_ss / real(nx*ny,rprec) )
       u_avg(:,:) = u_ss
    endif

!!$    if (.not. fourier) then
!!$       u_avg = u_avg / real(nx*ny,rprec)
!!$       call dfftw_execute_dft_r2c(forw, u_avg(:,:), u_avg(:,:))
!!$    endif
!!$    print*, 'wallo: >>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$    do jx=1,nx
!!$    do jy=1,ny
!!$       write(*,*) jx, jy, u_avg(jx,jy)
!!$    enddo
!!$    enddo

    !!!u_avg(1:nx,1:ny) = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)  !here
    ustar=u_avg*vonk/denom

if (fourier) then
   const=-((vonk/denom)**2) * u_avg(1,1)
   s1_big = 0._rprec
   s1_big(1,1) = const
   call dft_direct_back_2d_n_yonlyC_big(s1_big(:,:))
   
   txz_big(:,:) = convolve(s1_big(:,:), u1_big(:,:) )
   tyz_big(:,:) = convolve(s1_big(:,:), v1_big(:,:) )
   call dft_direct_forw_2d_n_yonlyC_big(txz_big(:,:))
   call dft_direct_forw_2d_n_yonlyC_big(tyz_big(:,:))
   call unpadd(txz(:,:,1), txz_big(:,:))
   call unpadd(tyz(:,:,1), tyz_big(:,:))
   
   !! same as for non-fourier case, just simplified since u_avg and vonk cancel out
   dudz(:,:,1)=1._rprec/(0.5_rprec*dz*denom) * u(:,:,1) * phi_m
   dvdz(:,:,1)=1._rprec/(0.5_rprec*dz*denom) * v(:,:,1) * phi_m
   do jy=1,ny
   do jx=1,nx
   dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
   dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
   enddo
   enddo
else
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
 endif

$if($DEBUG)
  case default

    write (*, *) 'invalid lbc_mom'
    stop
$endif
end select
end subroutine wallstress
