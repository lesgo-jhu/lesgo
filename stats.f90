!  This file contains (or will in the future) all subroutines
!  that compute statistics

!***************************************************************
subroutine tavg_compute()
!***************************************************************
!  This subroutine collects the stats for each flow 
!  variable quantity
use stat_defs, only : tavg_t, interp_to_uv_grid
use param, only : nx,ny,nz
use sim_param, only : u,v,w, dudz
integer :: i,j,k, navg
double precision :: w_interp, dudz_interp, fa

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
      w_interp = interp_to_uv_grid('w',i,j,k)  
	  dudz_interp = interp_to_uv_grid('dudz',i,j,k)
      
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

!**********************************************************************
subroutine plane_avg_compute(jt)
!**********************************************************************
use param, only : dx, dy, dz, Nx, Ny, Nz
use stat_defs,only:yplane_t, zplane_t, interp_to_uv_grid
use sim_param, only : u, v, w

implicit none

integer,intent(in)::jt

integer :: i,j,k
double precision :: ui, vi ,wi

!  Determine if y-plane averaging should be performed
if(jt >= yplane_t%nstart .and. jt <= yplane_t%nend) then
  
!  Compute average for each y-plane (jya)
  do j=1,yplane_t%na
    do k=1,Nz
      do i=1,Nx
	    ui = (u(i,yplane_t%istart(j)+1,k) - &
          u(i,yplane_t%istart(j),k))*yplane_t%ldiff(j)/dy + &
          u(i,yplane_t%istart(j),k)
		vi = (v(i,yplane_t%istart(j)+1,k) - &
          v(i,yplane_t%istart(j),k))*yplane_t%ldiff(j)/dy + &
          v(i,yplane_t%istart(j),k)
		wi = (interp_to_uv_grid('w',i,yplane_t%istart(j)+1,k) - &
          interp_to_uv_grid('w',i,yplane_t%istart(j),k))*yplane_t%ldiff(j)/dy + &
		  interp_to_uv_grid('w',i,yplane_t%istart(j),k)
        yplane_t%ua(i,j,k) = yplane_t%ua(i,j,k) + yplane_t%fa*ui 
	    yplane_t%va(i,j,k) = yplane_t%va(i,j,k) + yplane_t%fa*vi
	    yplane_t%wa(i,j,k) = yplane_t%wa(i,j,k) + yplane_t%fa*wi
	  enddo
	enddo
  enddo
endif 

!  Determine if y-plane averaging should be performed
if(jt >= zplane_t%nstart .and. jt <= zplane_t%nend) then
!  Compute average for each y-plane (jya)
  do k=1,zplane_t%na
    do j=1,Ny
      do i=1,Nx
        ui = (u(i,j,zplane_t%istart(k)+1) - u(i,j,zplane_t%istart(k)))*zplane_t%ldiff(k)/dz + &
		  u(i,j,zplane_t%istart(k))
        vi = (v(i,j,zplane_t%istart(k)+1) - v(i,j,zplane_t%istart(k)+1))*zplane_t%ldiff(k)/dz + &
		  v(i,j,zplane_t%istart(k)) 
        wi = (interp_to_uv_grid('w',i,j,zplane_t%istart(k)+1) - &
		  interp_to_uv_grid('w',i,j,zplane_t%istart(k)))*zplane_t%ldiff(k)/dz + &
          interp_to_uv_grid('w',i,j,zplane_t%istart(k))
        zplane_t%ua(i,j,k) = zplane_t%ua(i,j,k) + zplane_t%fa*ui
        zplane_t%va(i,j,k) = zplane_t%va(i,j,k) + zplane_t%fa*vi
        zplane_t%wa(i,j,k) = zplane_t%wa(i,j,k) + zplane_t%fa*wi
      enddo
    enddo
  enddo
endif 

return
end subroutine plane_avg_compute
