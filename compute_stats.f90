!***************************************************************
subroutine compute_stats()
!***************************************************************
!  This subroutine computes statistical values at the end
!  of the simulation for each i,j,k location
use stat_defs 
use param, only : nx,ny,nz,dx,dy,dz,z_i,L_x,L_y,L_z
use sim_param, only : u,v,w
use io, only : rs_write, tavg_write
implicit none

integer :: i,j,k
double precision :: u_avg

!  Check if Reynolds stresses are to be computed
if(rs_t%calc) then
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
  
  call rs_write()
endif

!  Check if average quantities are to be recorded
if(tavg_t%calc) call tavg_write()

    
return
end subroutine compute_stats

