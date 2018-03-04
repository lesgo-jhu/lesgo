module scalars_util

implicit none

integer, parameter :: tag_counter = 200

contains

!************************************************************************
subroutine scalar_transport (jt,theta,buoyancy)
!************************************************************************

use types, only: rprec
use derivatives, only: filt_da, ddz_uv, ddz_w
use param, only: nx, ny, nz, ld, dt, u_star, ld_big, &
                 initu, nx2, ny2, lbz, zo, z_i,   &
                 MPI_RPREC, down, up, comm, status, ierr, coord
use sgs_param, only : Nu_t
use sim_param, only: u, v, w, dTdx, dTdy, dTdz, RHS_T, RHS_Tf
use scalars_base, only: S_FLAG,inv_strength, T_scale, passive_scalar
use ic_scalar, only: init_scalar
!use ocean_base, only: inst_mixed_layer_depth

implicit none

integer :: jt, i, j, k,  jx, jy, jz!, k_ml_global
!integer, intent(inout) :: k_ml_global
real(rprec) :: zo_s
real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: theta
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: buoyancy
real(rprec), dimension(lbz:nz) :: theta_avg_z

!Assume scalar roughness length = (1/10)*momentum roughness length
zo_s=zo/10._rprec


$if($MPI) 
   call mpi_sendrecv (theta(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+7, &
                      theta(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+7, &
                      comm, status, ierr)
   call mpi_sendrecv (theta(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+8, &
                      theta(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+8, &
                     comm, status, ierr)
$endif


$if($MPI) 
   call mpi_sendrecv (dTdz(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+1, &
                      dTdz(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+1, &
                      comm, status, ierr)
   call mpi_sendrecv (dTdz(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+2, &
                      dTdz(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+2, &
                     comm, status, ierr)
$endif

!Calculate derivatives of theta
call filt_da(theta,dTdx,dTdy,lbz)

call ddz_uv(theta,dTdz,lbz)

!Top boundary condition
!dTdz(:,:,Nz) = inv_strength/T_scale*z_i

$if($MPI) 
   call mpi_sendrecv (w(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+3, &
                      w(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+3, &
                      comm, status, ierr)
   call mpi_sendrecv (w(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+4, &
                      w(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+4, &
                     comm, status, ierr)
$endif

$if($MPI) 
   call mpi_sendrecv (Nu_t(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+5, &
                      Nu_t(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+5, &
                      comm, status, ierr)
   call mpi_sendrecv (Nu_t(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+6, &
                      Nu_t(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+6, &
                     comm, status, ierr)
$endif

if (S_FLAG) then
   RHS_Tf = RHS_T

   !call inst_mixed_layer_depth(jt,theta,k_ml_global)

   call RHS_scalar_solver(theta,dTdx,dTdy,dTdz,zo_s,jt,RHS_T)

   if (.not.passive_scalar) then
      call buoyancy_temperature(theta,buoyancy)
   
   else
      buoyancy = 0._rprec
   end if

   call step_scalar(theta,RHS_T,RHS_Tf,jt)
end if 


!Investigate logarithmic temperature profile
!print*,'************',jt,'*************'
!do k= 0, Nz
!   $if ($MPI)
!       k_abs =  k
!   $else
!       k_abs = k
!   $endif
!   theta_avg_z(k) = 0.0_rprec
!   do j = 1, Ny
!      do i = 1, Nx
!         theta_avg_z(k) = theta_avg_z(k) + theta(i,j,k)
!      end do 
!   end do
!   theta_avg_z(k) = theta_avg_z(k)/real(Nx*Ny)
   
   !if (coord == 7) then 
   !    print*,'coord, k, theta_avg_z(k), buoyancy(16,16,k)=', coord, k, theta_avg_z(k), buoyancy (16,16,k)
   !end if
!end do

end subroutine scalar_transport

!************************************************************************
subroutine RHS_scalar_solver (theta,dsdx,dsdy,dsdz,zo_s,jt,RHS)
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, ld_big, z_i, dz, vonk, &
                 initu, nx2, ny2, lbz, zo,              &
                 coord, nproc, dz, ml_depth
use sim_param, only: u, v, w, sflux
use sgs_param, only: Nu_t
use derivatives, only: ddx, ddy, ddz_w, ddz_uv
use test_filtermodule
use ocean_base, only : prandtl_profile

implicit none

integer :: i, j, k, n, jt, jz, jy, jx, jz_min

real(rprec), dimension(ld,ny,lbz:nz),intent(in) :: theta
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdx, dsdy, dsdz
real(rprec), dimension(ld,ny,lbz:nz) :: RHS
real(rprec), dimension(nz) :: z
!integer, intent(in) :: k_ml_global
integer :: k_ml_global, k_abs
real(rprec), allocatable, dimension(:,:,:) :: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m,RHS_m
real(rprec), allocatable, dimension(:,:,:) :: temp,dtemp
real(rprec), dimension(ld,ny)::u1,v1
!real(rprec), dimension(nx,ny) :: sflux
real(rprec), intent(in) :: zo_s
real(rprec) :: sflux_mean, ustar_mean
real(rprec), dimension(nz) :: Pr

allocate(u_m(ld_big,ny2,lbz:nz))
allocate(v_m(ld_big,ny2,lbz:nz))
allocate(w_m(ld_big,ny2,lbz:nz))
allocate(dsdx_m(ld_big,ny2,lbz:nz))
allocate(dsdy_m(ld_big,ny2,lbz:nz))
allocate(dsdz_m(ld_big,ny2,lbz:nz))
allocate(RHS_m(ld_big,ny2,lbz:nz))
allocate(temp(ld,ny,lbz:nz))
allocate(dtemp(ld,ny,lbz:nz))

call dealias1(u,u_m,lbz)
call dealias1(v,v_m,lbz)
call dealias1(w,w_m,lbz)

call dealias1(dsdx,dsdx_m,lbz)
call dealias1(dsdy,dsdy_m,lbz)
call dealias1(dsdz,dsdz_m,lbz)

$if($MPI)
   if (coord ==0) then
      jz_min = 2
   else 
      jz_min = 1
   end if
$else
   jz_min = 2
$endif

!Compute RHS of scalar transport equation
!Advection term

do k=jz_min, Nz-1
   do j=1,Ny2
      do i=1,Nx2
         RHS_m(i,j,k) = u_m(i,j,k)*dsdx_m(i,j,k) + v_m(i,j,k)*dsdy_m(i,j,k) &
                      + (w_m(i,j,k)*dsdz_m(i,j,k) + w_m(i,j,k+1)*dsdz_m(i,j,k+1))*0.5_rprec
      end do
   end do
end do

$if($MPI)
    if (coord ==0) then
       do j=1,Ny2
          do i=1,Nx2
             RHS_m(i,j,1)  = u_m(i,j,1)*dsdx_m(i,j,1) + v_m(i,j,1)*dsdy_m(i,j,1) &
                             + (0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
          end do
       end do
    end if
$else
    do j=1,Ny2
       do i=1,Nx2
          RHS_m(i,j,1)  = u_m(i,j,1)*dsdx_m(i,j,1) + v_m(i,j,1)*dsdy_m(i,j,1) &
                          + (0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
       end do
    end do
$endif

$if($MPI)
    if (coord == nproc-1) then
       do j=1,Ny2
          do i=1,Nx2
             RHS_m(i,j,Nz) = u_m(i,j,Nz)*dsdx_m(i,j,Nz) + v_m(i,j,Nz)*dsdy_m(i,j,Nz)
          end do
       end do
    end if
$else
    do j=1,Ny2
       do i=1,Nx2
          RHS_m(i,j,Nz) = u_m(i,j,Nz)*dsdx_m(i,j,Nz) + v_m(i,j,Nz)*dsdy_m(i,j,Nz)
       end do
    end do
$endif

call dealias2(RHS,RHS_m,lbz)

!SGS part of the RHS
!Find the vertical level of the end of the mixed layer
!call thermocline_level(k_ml)
!print*,'jt,k_ml',jt,k_ml


!call inst_mixed_layer_depth(jt,theta,k_ml_global)
call prandtl_profile(jt,theta,Pr)


$if($MPI)
   if (coord == 0) then
       do j=1,Ny
          do i=1,Nx
             temp(i,j,1)=(1._rprec/Pr(1))*Nu_t(i,j,1)*dsdx(i,j,1)
          end do
       end do
   end if
$else
    do j=1,Ny
       do i=1,Nx
          temp(i,j,1)=(1._rprec/Pr(1))*Nu_t(i,j,1)*dsdx(i,j,1)
       end do
    end do
$endif

do k=jz_min,Nz-1
!   k_abs = coord*(Nz-1) + k
!   if (k_abs.le.k_ml) then  
      do j=1,Ny
         do i=1,Nx
            temp(i,j,k) = (1._rprec/Pr(k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
         end do
      end do
!   else
!      do j=1,Ny
!         do i=1,Nx
!            temp(i,j,k) = (1._rprec/diff_lam)*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
!         end do
!      end do
!   end if
end do

call ddx(temp,dtemp,lbz)

$if($MPI)
   if (coord == 0) then
      do j=1,Ny
         do i=1,Nx
            RHS(i,j,1) = (-1._rprec*RHS(i,j,1)) + dtemp(i,j,1)
            temp(i,j,1) = (1._rprec/Pr(1))*Nu_t(i,j,1)*dsdy(i,j,1)
         end do
      end do
   end if
$else
    do j=1,Ny
       do i=1,Nx
          RHS(i,j,1) = (-1._rprec*RHS(i,j,1)) + dtemp(i,j,1)
          temp(i,j,1) = (1._rprec/Pr(1))*Nu_t(i,j,1)*dsdy(i,j,1)
       end do
    end do
$endif


do k=jz_min,Nz-1
   do j=1,Ny
      do i=1,Nx
         RHS(i,j,k) = (-1._rprec*RHS(i,j,k) + dtemp(i,j,k))
      end do 
   end do
end do

do k=jz_min,Nz-1
!   k_abs = coord*(Nz-1) + k
!   if (k_abs.le.k_ml) then  
      do j=1,Ny
         do i=1,Nx
            temp(i,j,k) = (1._rprec/Pr(k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
         end do 
      end do
!   else
!      do j=1,Ny
!         do i=1,Nx
!            temp(i,j,k) = (1._rprec/diff_lam)*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
!         end do 
!      end do
!   end if
end do

call ddy(temp,dtemp,lbz)

$if($MPI)
   if (coord == 0) then 
      do j=1,Ny
         do i=1,Nx
            RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
         end do
      end do
      ! Get surface flux from similarity relationship
      call surface_flux(theta)
      temp(:,:,1) = -1.0_rprec*sflux  
   end if
$else
    do j=1,Ny
       do i=1,Nx
          RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
       end do
    end do
    ! Get surface flux from similarity relationship
    call surface_flux(theta)
    temp(:,:,1) = -1.0_rprec*sflux  
$endif

do k=jz_min,Nz
   do j=1,Ny
      do i=1,Nx
         RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
      end do
   end do
end do


do k=jz_min,Nz
!   k_abs = coord*(Nz-1) + k
!   if (k_abs.le.k_ml) then  
      do j=1,Ny
         do i=1,Nx
            temp(i,j,k) = (1._rprec/Pr(k))*Nu_t(i,j,k)*dsdz(i,j,k)
         end do
      end do
!   else  
!      do j=1,Ny
!         do i=1,Nx
!            temp(i,j,k) = (1._rprec/diff_lam)*Nu_t(i,j,k)*dsdz(i,j,k)
!         end do
!      end do
!   end if
end do

call ddz_w(temp,dtemp,lbz)

do k=1,Nz-1
   do j=1,Ny
      do i=1,Nx
         RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
      end do
   end do
end do

end subroutine RHS_scalar_solver


!************************************************************************
subroutine step_scalar (scalar, RHS_pre, RHS_post, jt)
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, l_z, z_i, dz, lbz, dt, &
                 coord, nproc
 
implicit none

integer :: i, j, k, jt, jz

real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: scalar 
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: RHS_pre, RHS_post

do k=1,nz
   do j=1,ny
      do i=1,nx
         scalar(i,j,k) = scalar(i,j,k) + &
                         dt*(1.5_rprec*RHS_pre(i,j,k) - 0.5_rprec*RHS_post(i,j,k))
      end do
   end do
end do

$if($MPI)
    if (coord == nproc-1) then
    !For neurtral and passive scalars
    scalar(:,:,Nz-1) = scalar(:,:,Nz)
    end if
$else
    scalar(:,:,Nz-1) = scalar(:,:,Nz)
$endif

end subroutine step_scalar


!************************************************************************
subroutine surface_flux (theta)
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, dz, vonk, &
                 lbz, zo
use sim_param, only: u, v, ustar,sflux
use test_filtermodule

implicit none

integer :: i, j

real(rprec), dimension(ld,ny,lbz:nz),intent(in) :: theta
!real(rprec), dimension(nx,ny), intent(out) :: sflux
real(rprec),dimension(nx,ny)::u_avg
real(rprec),dimension(ld,ny)::u1,v1,theta1
real(rprec) :: denom, zo_s
!real(rprec) :: sflux_mean, ustar_mean

!Assume scalar roughness length = (1/10)*momentum roughness length
zo_s=zo/10._rprec

!!Wall model for dsdz
!Calculate friction velocity ustar
u1=u(:,:,1)
v1=v(:,:,1)
theta1=theta(:,:,1)
call test_filter ( u1 )
call test_filter ( v1 )
call test_filter ( theta1 ) 
denom=log(0.5_rprec*dz/zo)
u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar=u_avg*vonk/denom

do j=1,Ny
   do i=1,Nx
    !sflux(i,j) =  (330.0_rprec/300.0_rprec - theta1(i,j))/((1._rprec/vonk)*log(0.5_rprec*dz/zo_s)) * ustar(i,j)
    sflux(i,j) = 6.425e-7_rprec
   end do
end do

end subroutine surface_flux



!************************************************************************
subroutine buoyancy_temperature (theta,buoyancy)
!************************************************************************
!This subroutine calculates the buoyancy term due to temperature to be 
!added to the RHS of the vertical momentum equation. 

use types, only: rprec
use param, only: z_i, u_star, ld, nx, ny, nz, lbz, coord
use ocean_base, only: alpha_theta
use scalars_base, only: T_scale

implicit none

integer :: i,j,k,jz_min
real(rprec), parameter :: g=9.81_rprec
real(rprec), dimension(ld,ny,nz), intent(in) :: theta
real(rprec), dimension(ld,ny,nz), intent(out) :: buoyancy
real(rprec), dimension(nz) :: theta_bar
real(rprec) :: g_hat,above,below

!Non-dimensionalize gravity
g_hat = g*(z_i/(u_star**2))

! We do not advance the ground nodes, so start at k=2. 
! For the MPI case, the means that we start from jz=2
! for coord=0 and jz=1 otherwise. 
$if($MPI)
   if (coord ==0) then
      jz_min = 2
   else 
      jz_min = 1
   end if
$else
   jz_min = 2
$endif

!Ocean flag

do k=jz_min, nz-1
   do j=1,ny
      do i=1,nx
         above=theta(i,j,k) 
         below=theta(i,j,k-1) 
         buoyancy(i,j,k) = g_hat*(1._rprec - alpha_theta*((above+below)/2._rprec)*T_scale )
      end do
   end do
end do

end subroutine buoyancy_temperature

end module scalars_util
