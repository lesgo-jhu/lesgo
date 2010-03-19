module turbines
use types,only:rprec
use param, only: dx,lh,ny
use stat_defs, only:turbine_t

implicit none

save
private

public :: turbines_init, turbines_forcing

real(rprec) :: Ct									!thrust coefficient
integer :: flag										!indicates fist time step (where u_d_T = u_d)
real :: eps											!epsilon used for disk velocity time-averaging

real, pointer, dimension(:) :: u_d_T				!running time-average of mean disk velocity
integer, pointer, dimension(:) :: num_nodes			!number of nodes associated with each turbine
integer, pointer, dimension(:,:) ::nodes			!(i,j,k) of each included node
real, pointer, dimension (:,:) :: n_hat				!(nx,ny,nz) of unit normal for each turbine

real, dimension(lh,ny) :: G_turbines_xy				!for filtering forces (transfer function)
integer*8::forw_xy,back_xy							!FFT plans
real, dimension(lh,ny)::kxT,kyT,k2T					!wavenumbers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_init()
!PROBLEMS WITH MPI - MAKE SURE TURBINES ARE COMPLETELY WITHIN ONE PROC'S DOMAIN
use param, only: dt, L_x, L_y, L_z, u_star


!##############################  SET BY USER  ############################################

turbine_t%nloc = 2      		!number of locations 
allocate(turbine_t%dia(turbine_t%nloc))
allocate(turbine_t%thk(turbine_t%nloc))
allocate(turbine_t%height(turbine_t%nloc))
allocate(turbine_t%xloc(turbine_t%nloc))
allocate(turbine_t%yloc(turbine_t%nloc))
allocate(turbine_t%theta1(turbine_t%nloc))
allocate(turbine_t%theta2(turbine_t%nloc))
	turbine_t%dia(1) = 100./1000. 		!turbine diameter  [norm by z_i]  
	turbine_t%thk(1) = dx				!turbine disk thickness 
	turbine_t%height(1) = 100./1000.	!turbine height    [norm by z_i] 
	turbine_t%xloc(1) = L_x*1./4.   !x locations [dimensionless - usually in terms of L_x]
	turbine_t%yloc(1) = L_y/2.    	!y locations [dimensionless - usually in terms of L_y]
	turbine_t%theta1(1) = 0.			!angle CCW(from above) from -x direction [degrees]
	turbine_t%theta2(1) = 0.			!angle above the horizontal [degrees]	

	turbine_t%dia(2) = 100./1000. 
	turbine_t%thk(2) = dx	
	turbine_t%height(2) = 100./1000.		
	turbine_t%xloc(2) = L_x*3./4.
	turbine_t%yloc(2) = L_y/2.
	turbine_t%theta1(2) = 0.
	turbine_t%theta2(2) = 0.
turbine_t%ifilter_xy = 2			!Filter type: 1->cut off 2->Gaussian 3->Top-hat
turbine_t%ifilter_z = 4				!Filter type: 4->Truncated Gaussian		
turbine_t%alpha_xy = 1.5			!filter size is alpha*(grid spacing)
turbine_t%alpha_z = 1.5  			!h=horizontal (xy), v=vertical (xz,yz)

Ct = 1.33		!thrust coefficient
!#########################################################################################

allocate(u_d_T(turbine_t%nloc))
allocate(num_nodes(turbine_t%nloc))
allocate(nodes(100*turbine_t%nloc,3))
allocate(n_hat(turbine_t%nloc,3))

!locate applicable nodes
	call turbines_nodes()
!initialize FFT (create plans and initialize wavenumbers)
	call turbines_fft()	
!initialize filter	
	call turbines_filter_init()
!initialize variables
	u_d_T = 0.
	flag = 1
	eps = dt/(0.27*L_z/u_star) / (1 + dt/(0.27*L_z/u_star))
	
end subroutine turbines_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_forcing()
use param, only: nx,ny,nz,pi,dx,dy,dz,dt
use sim_param, only: u,v,w
use io, only: w_uv, w_uv_tag, interp_to_uv_grid

real, dimension(nx,ny,nz) :: f_x=0., f_y=0., f_z=0.   !only initialized during first call?
real, dimension(turbine_t%nloc) :: u_d, f_t, cf
integer :: i,j,k,m,l,i2,j2,k2,count_i

!calculate total disk-averaged velocity for each turbine (current,instantaneous)
	call interp_to_uv_grid(w, w_uv, w_uv_tag)
	count_i = 0		!index count - used for reading array "nodes"
	do m=1,turbine_t%nloc
		u_d(m) = 0.
		do l=1,num_nodes(m)
			i2 = nodes(count_i+l,1)
			j2 = nodes(count_i+l,2)
			k2 = nodes(count_i+l,3)			
			u_d(m) = u_d(m) + (n_hat(m,1)*u(i2,j2,k2) + &
			n_hat(m,2)*v(i2,j2,k2) + n_hat(m,3)*w_uv(i2,j2,k2))/num_nodes(m)
		enddo
		count_i = count_i + num_nodes(m)
	enddo
!add this current value to the "running average" (first order relaxation, T=0.27*H/u_star)
	if (flag) then
		u_d_T = u_d
		flag = 0
	else
		u_d_T = (1.-eps)*u_d_T + eps*u_d
	endif

!calculate total thrust force (per unit mass) for each turbine (includes volume/mass correction factor)
	do m=1,turbine_t%nloc
		cf(m) = pi/4.*turbine_t%dia(m)*turbine_t%dia(m)*turbine_t%thk(m) / (dx*dy*dz*num_nodes(m))
		f_t(m) = -0.5*Ct*u_d_T(m)*u_d_T(m)/turbine_t%thk(m) * cf(m)
	enddo
	
!apply forcing to each node
	count_i = 0		!index count
	do m=1,turbine_t%nloc
		do l=1,num_nodes(m)
			i2 = nodes(count_i+l,1)
			j2 = nodes(count_i+l,2)
			k2 = nodes(count_i+l,3)			
			f_x(i2,j2,k2) = f_t(m)*n_hat(m,1) 
			f_y(i2,j2,k2) = f_t(m)*n_hat(m,2) 
			f_z(i2,j2,k2) = f_t(m)*n_hat(m,3) 
		enddo
		count_i = count_i + num_nodes(m)
	enddo
		
!filter forcing (currently only in the x-y planes)
	do k=1,nz
		call turbines_filter_xy(f_x(:,:,k))		!filters in place so input is ruined
		call turbines_filter_xy(f_y(:,:,k)) 
		call turbines_filter_xy(f_z(:,:,k)) 
	enddo

!apply forcing to RHS

end subroutine turbines_forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_nodes()
!This subroutine locates nodes for each turbine and builds the arrays: n_hat, num_nodes, and nodes
use param, only:dx,dy,dz,pi	
use grid_defs, only:x,y,z
	
integer :: imax,jmax,kmax,count_i,count_n,m,icp,jcp,kcp,i,j,k
real :: R_t,rx,ry,rz,r,r_norm,r_disk

count_i = 1		!index count - used for writing to array "nodes"
count_n = 0		!used for counting nodes for each turbine

do m=1,turbine_t%nloc
	R_t = turbine_t%dia(m)/2.
	imax = R_t/dx + 1
	jmax = R_t/dy + 1
	kmax = R_t/dz + 1
!determine unit normal vector for each turbine	
	n_hat(m,1) = -cos(pi*turbine_t%theta1(m)/180.)*cos(pi*turbine_t%theta2(m)/180.)
	n_hat(m,2) = -sin(pi*turbine_t%theta1(m)/180.)*cos(pi*turbine_t%theta2(m)/180.)
	n_hat(m,3) = sin(pi*turbine_t%theta2(m)/180.)
!determine nearest (i,j,k) to turbine center
	icp = nint(turbine_t%xloc(m)/dx)
	jcp = nint(turbine_t%yloc(m)/dy)
	kcp = nint(turbine_t%height(m)/dz + 0.5)
!check neighboring grid points	
	do i=(icp-imax),(icp+imax)
		do j=(jcp-jmax),(jcp+jmax)
			do k=(kcp-kmax),(kcp+kmax)
			!vector from center point to this node is (rx,ry,rz) with length r
				rx = x(i) - turbine_t%xloc(m)
				ry = y(j) - turbine_t%yloc(m)
				rz = z(k) - turbine_t%height(m)
				r = sqrt(rx*rx + ry*ry + rz*rz)
			!length projected onto unit normal for this turbine
				r_norm = abs(rx*n_hat(m,1) + ry*n_hat(m,2) + rz*n_hat(m,3))
			!(remaining) length projected onto turbine disk
				r_disk = sqrt(r*r - r_norm*r_norm)
			!if r_disk<R_t and r_norm within thk/2 from turbine -- this node is part of the turbine
				if ( (r <= R_t) .AND. (r_norm <= turbine_t%thk(m)/2.) ) then
					count_n = count_n + 1
					nodes(count_i,1) = i
					nodes(count_i,2) = j
					nodes(count_i,3) = k
					count_i = count_i + 1
				endif
			enddo
		enddo
	enddo
	num_nodes(m) = count_n
	count_n = 0
enddo
		
end subroutine turbines_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_filter_init()
! G_turbines should take care of normalization and Nyquist
! note the normalization for FFT's is already in G! (see 1/(nx*ny))
! currently: only prepares G_turbines_xy (not "_z)
use param, only: nx,ny,dx,dy,pi

real(kind=rprec):: delta, kc2

	G_turbines_xy=1._rprec/(nx*ny)  		! normalization for the forward FFT
	delta = turbine_t%alpha_xy*sqrt(dx*dy)  ! "2d-delta", not full 3d one

! 1. Spectral cutoff filter
	if(turbine_t%ifilter_xy==1) then 
		kc2 = (pi/(delta))**2
		where (real(k2T, rprec) >= kc2) G_turbines_xy = 0._rprec
! 2. Gaussian filter
	else if(turbine_t%ifilter_xy==2) then 
		G_turbines_xy=exp(-(2._rprec*delta)**2*k2T/(4._rprec*6._rprec))*G_turbines_xy       
! 3. Top-hat (Box) filter
	else if(turbine_t%ifilter_xy==3) then 
		G_turbines_xy= (sin(kxT*delta/2._rprec)*sin(kyT*delta/2._rprec)+1E-8)/&
        (kxT*delta/2._rprec*kyT*delta/2._rprec+1E-8)*G_turbines_xy
	endif

! since our k2T has zero at Nyquist, we have to do this by hand
	G_turbines_xy(nx/2+1,:) = 0._rprec
	G_turbines_xy(:,ny/2+1) = 0._rprec

end subroutine turbines_filter_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_filter_xy(real_a)
use param,only:nx,ny

real, dimension(nx,ny), intent(inout) :: real_a
complex, dimension(nx/2+1,ny) :: complex_a

!perform FFT (physical to wavenumber space)
	call rfftwnd_f77_one_real_to_complex(forw_xy, real_a, complex_a)
!perform filter convolution
	complex_a = complex_a * G_turbines_xy
!perform FFT (wavenumber to physical space)
	call rfftwnd_f77_one_complex_to_real(back_xy, complex_a, real_a)
	
!normalization is part of G_turbines (check that it is done correctly)

end subroutine turbines_filter_xy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_fft()
!arrays should be in column-major order
use param,only:nx,ny,L_x,L_y,pi

integer, parameter :: FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1
integer, parameter :: FFTW_ESTIMATE=0,FFTW_MEASURE=1
integer, parameter :: FFTW_OUT_OF_PLACE=0
integer, parameter :: FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16
integer, parameter :: FFTW_THREADSAFE=128

integer :: jx,jy

!create my own plans for FFT (that don't filter in place)
	call rfftw2d_f77_create_plan(forw_xy,nx,ny,FFTW_REAL_TO_COMPLEX,&
     FFTW_MEASURE+FFTW_THREADSAFE)
	call rfftw2d_f77_create_plan(back_xy,nx,ny,FFTW_COMPLEX_TO_REAL,&
     FFTW_MEASURE+FFTW_THREADSAFE)

!also initialize wavenumber arrays - to be used during filtering?
	do jx=1,(nx/2+1)-1
		kxT(jx,:) = real(jx-1,kind=rprec)
	end do

	do jy=1,ny
	kyT(:,jy) = real(modulo(jy - 1 + ny/2,ny) - ny/2,kind=rprec)
	end do

	! Nyquist: makes doing derivatives easier
    kxT(nx/2+1,:)=0._rprec
    kyT(nx/2+1,:)=0._rprec
    kxT(:,ny/2+1)=0._rprec
    kyT(:,ny/2+1)=0._rprec
	  
	! for the aspect ratio change
    kxT=2._rprec*pi/L_x*kxT
    kyT=2._rprec*pi/L_y*kyT
	  
	! magnitude squared: will have 0's around the edge
    k2T = kxT*kxT + kyT*kyT

end subroutine turbines_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module turbines
