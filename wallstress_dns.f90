subroutine wallstress_dns ()
use types,only:rprec
use param,only:ld,nx,ny,nz,nu_molec,z_i,u_star,dz,lbc_mom
use sim_param,only:u,v,dudz,dvdz, txz, tyz
implicit none
integer::jx,jy

select case (lbc_mom)

  case (0) ! Stress free

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case (1) ! Wall
    do jy=1,ny
    do jx=1,nx
       txz(jx,jy,1)=-nu_molec/(z_i*u_star)*u(jx,jy,1)/(0.5_rprec*dz)
       tyz(jx,jy,1)=-nu_molec/(z_i*u_star)*v(jx,jy,1)/(0.5_rprec*dz)
       dudz(jx,jy,1)=u(jx,jy,1)/(0.5_rprec*dz)
       dvdz(jx,jy,1)=v(jx,jy,1)/(0.5_rprec*dz)
    end do
    end do

  case default

    write (*, *) 'invalid lbc_mom'

end select

end subroutine wallstress_dns
