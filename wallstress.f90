! For use with staggered grid LES
! JDA, 23 Jan 96
! zo is nondimensionalized, zo1 not!
!--provides txz, tyz, dudz, dvdz at jz=1
subroutine wallstress ()
use types,only:rprec
use param, only : vonk
use param2,only:dz,ld,lh,nx,ny,nz,lbc_mom
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use bottombc,only:zo,psi_m,phi_m
use test_filtermodule
implicit none
integer::jx,jy
real(kind=rprec),dimension(nx,ny)::ustar,u_avg,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec)::const

select case (lbc_mom)

  case ('wall')

    !TS Remove the following line when obukhov.f is used
    psi_m=0._rprec
    phi_m=1._rprec

    u1=u(:,:,1)
    v1=v(:,:,1)
    call test_filter(u1,G_test)
    call test_filter(v1,G_test)
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

  case ('stress free')

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case default

    write (*, *) 'invalid lbc_mom'
    stop

end select

end subroutine wallstress
