!***************************************************************
subroutine compute_stats()
!***************************************************************
!  This subroutine computes statistical values at the end
!  of the simulation for each i,j,k location
use stat_defs 
use param, only : nx,ny,nz,dx,dy,dz,z_i,L_x,L_y,L_z
use sim_param, only : u,v,w
implicit none

integer :: i,j,k
double precision :: x(0:nx), y(0:ny), z(0:nz)
double precision :: dnx, dny, dnz, u_aver

!  Convert integer quantities to double precision
dnx=dble(nx)
dny=dble(ny)
dnz=dble(nz)

!  Compute x
x(0)=-dx
do i=1,nx
  x(i)=x(i-1)+dx
enddo
!  Compute y
y(0)=-dy
do j=1,ny
  y(j)=y(j-1)+dy
enddo
!  Compute z
z(0)=-dz*z_i/2.
do k=1,nz
  z(k)=z(k-1)+dz*z_i
enddo

!  Non-dimensionalize by L_x, L_y, and L_z
x=x/L_x
y=y/L_y
z=z/L_z

!  Check if Reynolds stresses are to be computed
if(rs_t%calc) then
  do k=1,nz
    do j=1,ny
      do i=1,nx
  	    rs_t%up2(i,j,k)=taver_t%u2(i,j,k) - taver_t%u(i,j,k)*taver_t%u(i,j,k)
	    rs_t%vp2(i,j,k)=taver_t%v2(i,j,k) - taver_t%v(i,j,k)*taver_t%v(i,j,k)
	    rs_t%wp2(i,j,k)=taver_t%w2(i,j,k) - taver_t%w(i,j,k)*taver_t%w(i,j,k)
	    rs_t%upwp(i,j,k)=taver_t%uw(i,j,k) - taver_t%u(i,j,k)*taver_t%w(i,j,k)
	    rs_t%vpwp(i,j,k)=taver_t%vw(i,j,k) - taver_t%v(i,j,k)*taver_t%w(i,j,k)
	    rs_t%upvp(i,j,k)=taver_t%uv(i,j,k) - taver_t%u(i,j,k)*taver_t%v(i,j,k)
	  enddo
    enddo
  enddo

!  Add temporary output information
  open(unit = 7,file = 'output/aver_rs.out')
  write(7,*) 'variables= "z", "<up2>", "<vp2>", "<wp2>", "<upwp>", "<vpwp>", "<upvp>"'
  do k=1,nz
!  Write spatially averaged, temporally averaged quantities
    write(7,*) z(k), sum(rs_t%up2(:,:,k))/(dnx*dny), sum(rs_t%vp2(:,:,k))/(dnx*dny), &
      sum(rs_t%wp2(:,:,k))/(dnx*dny), sum(rs_t%upwp(:,:,k))/(dnx*dny), &
	  sum(rs_t%vpwp(:,:,k))/(dnx*dny), sum(rs_t%upvp(:,:,k))/(dnx*dny)
  enddo
  close(7)
endif

!  Check if average quantities are to be recorded
if(aver_calc) then
  open(unit = 7,file = 'output/aver_uvw.out')
  open(unit = 8,file = 'output/aver_dudz.out')
  write(7,*) 'variables= "z", "<u>/u*", "<v>/u*", "<w>/u*"'
  write(8,*) 'variables= "z", "<dudz>/u*"'
  do k=1,nz
	write(7,*) z(k), sum(taver_t%u(:,:,k))/(dnx*dny), sum(taver_t%v(:,:,k))/(dnx*dny), & 
	  sum(taver_t%w(:,:,k))/(dnx*dny)
	write(8,*) z(k), sum(taver_t%dudz(:,:,k))/(dnx*dny)
  enddo
  close(7)
  close(8)
endif

!if(final_calc) then
!  open(unit = 7,file = 'output/final.out', status='unknown',form='formatted', &
!        action='write',position='rewind')
!  write(7,*) 'variables= "x/Lx", "y/Ly", "z/H", "u/u*", "v/u*", "w/u*", "(u - <u>)/u*"'
!  write(7,"(1a,i4,1a,i4,1a,i4,1a,i4,1a,i4)") 'ZONE T="', 1,'", DATAPACKING=POINT, i=', &
!    nx,', j=', ny,', k=', nz
!  write(7,"(1a)") ''//adjustl('DT=(SINGLE SINGLE SINGLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
!  
!  u_aver = sum(u)/(dnx*dny*dnz)
!  write(*,*) 'u_aver =', u_aver
!  do k=1,nz
!    do j=1,ny
!	  do i=1,nx
!	    write(7,*) x(i), y(j), z(k), u(i,j,k), v(i,j,k), interp_to_uv_grid('w',i,j,k), u(i,j,k) - u_aver
!      enddo
!	enddo
!  enddo
  
!  close(7)
  
!endif
	    
return
end subroutine compute_stats

