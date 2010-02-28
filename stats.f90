!  This file contains (or will in the future) all subroutines
!  that compute statistics

!***************************************************************
subroutine tsum_compute()
!***************************************************************
!  This subroutine collects the stats for each flow 
!  variable quantity
use types, only : rprec
use stat_defs, only : tsum_t
use param, only : nx,ny,nz
use sim_param, only : u,v,w, dudz
use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
integer :: i,j,k
real(rprec) :: u_p, v_p, w_p, dudz_p

!  Make sure w has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)
call interp_to_uv_grid(dudz, dudz_uv, dudz_uv_tag)

do k=1,nz
  do j=1,ny
    do i=1,nx
!  Being cache friendly
      u_p = u(i,j,k)
      v_p = v(i,j,k)
!  Interpolate each w and dudz to uv grid
      w_p = w_uv(i,j,k)  
      dudz_p = dudz_uv(i,j,k)
      
      tsum_t%u(i,j,k)=tsum_t%u(i,j,k) + u_p
      tsum_t%v(i,j,k)=tsum_t%v(i,j,k) + v_p
      tsum_t%w(i,j,k)=tsum_t%w(i,j,k) + w_p
      tsum_t%u2(i,j,k)=tsum_t%u2(i,j,k) + u_p*u_p
      tsum_t%v2(i,j,k)=tsum_t%v2(i,j,k) + v_p*v_p
      tsum_t%w2(i,j,k)=tsum_t%w2(i,j,k) + w_p*w_p
      tsum_t%uw(i,j,k)=tsum_t%uw(i,j,k) + u_p*w_p
      tsum_t%vw(i,j,k)=tsum_t%vw(i,j,k) + v_p*w_p
      tsum_t%uv(i,j,k)=tsum_t%uv(i,j,k) + u_p*v_p
      tsum_t%dudz(i,j,k)=tsum_t%dudz(i,j,k) + dudz_p
    enddo
  enddo
enddo

return

end subroutine tsum_compute

!***************************************************************
subroutine tavg_compute()
!***************************************************************
!  This subroutine collects the stats for each flow 
!  variable quantity
use stat_defs, only : tavg_t
use param, only : nx,ny,nz
use sim_param, only : u,v,w, dudz
use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
implicit none
integer :: i,j,k, navg
double precision :: w_interp, dudz_interp, fa

!  Make sure w has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)
call interp_to_uv_grid(dudz, dudz_uv, dudz_uv_tag)

!  Initialize w_interp and dudz_interp
w_interp = 0.
dudz_interp = 0. 

!  Compute number of times to average over
navg = tavg_t%nend - tavg_t%nstart + 1
!  Averaging factor
fa=1./dble(navg)

do k=1,nz
  do j=1,ny
    do i=1,nx
!  Interpolate each w and dudz to uv grid
      w_interp = w_uv(i,j,k)  
      dudz_interp = dudz_uv(i,j,k)
      
      tavg_t%u(i,j,k)=tavg_t%u(i,j,k) + fa*u(i,j,k)
      tavg_t%v(i,j,k)=tavg_t%v(i,j,k) + fa*v(i,j,k)
      tavg_t%w(i,j,k)=tavg_t%w(i,j,k) + fa*w_interp
      tavg_t%u2(i,j,k)=tavg_t%u2(i,j,k) + fa*u(i,j,k)*u(i,j,k)
      tavg_t%v2(i,j,k)=tavg_t%v2(i,j,k) + fa*v(i,j,k)*v(i,j,k)
      tavg_t%w2(i,j,k)=tavg_t%w2(i,j,k) + fa*w_interp*w_interp
      tavg_t%uw(i,j,k)=tavg_t%uw(i,j,k) + fa*u(i,j,k)*w_interp
      tavg_t%vw(i,j,k)=tavg_t%vw(i,j,k) + fa*v(i,j,k)*w_interp
      tavg_t%uv(i,j,k)=tavg_t%uv(i,j,k) + fa*u(i,j,k)*v(i,j,k)
      tavg_t%dudz(i,j,k)=tavg_t%dudz(i,j,k) + fa*dudz_interp
    enddo
  enddo
enddo

return

end subroutine tavg_compute

!***************************************************************
subroutine rs_compute()
!***************************************************************
!  This subroutine computes Reynolds stresses from tavg_t u,v,w
!  values
use stat_defs, only : rs_t, tavg_t
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z
use sim_param, only : u,v,w

implicit none

integer :: i,j,k

do k=1,nz
  do j=1,ny
    do i=1,nx
      rs_t%up2(i,j,k)=tavg_t%u2(i,j,k) - tavg_t%u(i,j,k)*tavg_t%u(i,j,k)
      rs_t%vp2(i,j,k)=tavg_t%v2(i,j,k) - tavg_t%v(i,j,k)*tavg_t%v(i,j,k)
	  rs_t%wp2(i,j,k)=tavg_t%w2(i,j,k) - tavg_t%w(i,j,k)*tavg_t%w(i,j,k)
	  rs_t%upwp(i,j,k)=tavg_t%uw(i,j,k) - tavg_t%u(i,j,k)*tavg_t%w(i,j,k)
	  rs_t%vpwp(i,j,k)=tavg_t%vw(i,j,k) - tavg_t%v(i,j,k)*tavg_t%w(i,j,k)
	  rs_t%upvp(i,j,k)=tavg_t%uv(i,j,k) - tavg_t%u(i,j,k)*tavg_t%v(i,j,k)
	enddo
  enddo
enddo
  
return
end subroutine rs_compute

! !**********************************************************************
! subroutine plane_avg_compute(jt)
! !**********************************************************************
! use functions, only : linear_interp, interp_to_uv_grid
! use grid_defs, only : z
! use param, only : dx, dy, dz, Nx, Ny, Nz,coord
! use stat_defs,only:yplane_t, zplane_t
! use sim_param, only : u, v, w
! 
! implicit none
! 
! integer,intent(in)::jt
! 
! integer :: i,j,k
! double precision :: ui, vi ,wi
! 
! !  Determine if y-plane averaging should be performed
! if(yplane_t%calc .and. jt >= yplane_t%nstart .and. jt <= yplane_t%nend) then
!   
! !  Compute average for each y-plane (jya)
!   do j=1,yplane_t%nloc
!     do k=1,Nz
!       do i=1,Nx
!         ui = linear_interp(u(i,yplane_t%istart(j),k), &
!           u(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
!         vi = linear_interp(v(i,yplane_t%istart(j),k), &
!           v(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
!         wi = linear_interp(interp_to_uv_grid('w',i,yplane_t%istart(j),k), &
!           interp_to_uv_grid('w',i,yplane_t%istart(j)+1,k), dy, &
!           yplane_t%ldiff(j))
!         yplane_t%ua(i,j,k) = yplane_t%ua(i,j,k) + yplane_t%fa*ui 
!         yplane_t%va(i,j,k) = yplane_t%va(i,j,k) + yplane_t%fa*vi
!         yplane_t%wa(i,j,k) = yplane_t%wa(i,j,k) + yplane_t%fa*wi
! 
!       enddo
!     enddo
!   enddo
! endif 
! 
! !  Determine if y-plane averaging should be performed
! if(zplane_t%calc .and. jt >= zplane_t%nstart .and. jt <= zplane_t%nend) then
! !  Compute average for each y-plane (jya)
! 
!   do k=1,zplane_t%nloc
! 
!   $if ($MPI)
!   if(zplane_t%coord(k) == coord) then
!   $endif
! 
!     do j=1,Ny
!       do i=1,Nx
! 
!         ui = linear_interp(u(i,j,zplane_t%istart(k)),u(i,j,zplane_t%istart(k)+1), &
!           dz, zplane_t%ldiff(k))
!         vi = linear_interp(v(i,j,zplane_t%istart(k)),v(i,j,zplane_t%istart(k)+1), &
!           dz, zplane_t%ldiff(k))
!         wi = linear_interp(interp_to_uv_grid('w',i,j,zplane_t%istart(k)), &
!           interp_to_uv_grid('w',i,j,zplane_t%istart(k)+1), &
!           dz, zplane_t%ldiff(k))
! 
!         zplane_t%ua(i,j,k) = zplane_t%ua(i,j,k) + zplane_t%fa*ui
!         zplane_t%va(i,j,k) = zplane_t%va(i,j,k) + zplane_t%fa*vi
!         zplane_t%wa(i,j,k) = zplane_t%wa(i,j,k) + zplane_t%fa*wi
! 
!       enddo
!     enddo
! 
!     $if ($MPI)
!     endif
!     $endif
! 
!   enddo
! 
! endif
! 
! return
! end subroutine plane_avg_compute
