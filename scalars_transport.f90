module scalars_transport

implicit none

integer, parameter :: tag_counter = 200

contains

!************************************************************************
subroutine temperature_transport ()
!************************************************************************
use types, only: rprec
use derivatives, only: filt_da, ddz_uv
use param, only: nx, ny, nz, ld, lbz, z_i, &
                 MPI_RPREC, down, up, comm, status, ierr, coord, & 
                 ssgs_init, nproc, dz, jt
use sgs_param, only : Nu_t, kappa_tt, Ds_opt2_t
use sim_param, only: u, v, w, theta, dTdx, dTdy, dTdz, RHS_T, RHS_Tf
use sim_param, only: t_flux, sigma_theta
use scalars_base, only: temperature_surface_bc
use scalars_base, only: dTdz_top 
use scalars_util, only: RHS_scalar_solver, step_scalar
use ocean_base, only: prandtl_profile
use sgs_param, only: s_Tn_all_t, s_Beta_t, I_LM_t, I_MM_t, I_QN_t, I_NN_t, ds2_clips_t
use sgs_param, only: I_LM_MM_init_t, I_QN_NN_init_t

implicit none

integer :: jx, jy, jz

!Calculate derivatives of theta
call filt_da(theta,dTdx,dTdy,lbz)

call ddz_uv(theta,dTdz,lbz)

$if($MPI)
   if (coord==nproc-1) then
      do jy=1,ny
         do jx=1,nx
            dTdz(jx,jy,nz) = (1._rprec/dz)*(theta(jx,jy,nz)-theta(jx,jy,nz-1))
         end do
      end do
   end if
$else
      do jy=1,ny
         do jx=1,nx
            dTdz(jx,jy,nz) = (1._rprec/dz)*(theta(jx,jy,nz)-theta(jx,jy,nz-1))
         end do
      end do
$endif

$if($MPI) 
   call mpi_sendrecv (dTdz(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+1, &
                      dTdz(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+1, &
                      comm, status, ierr)
   call mpi_sendrecv (dTdz(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+2, &
                      dTdz(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+2, &
                     comm, status, ierr)
$endif

RHS_Tf = RHS_T

call temperature_surface_bc()

!if (jt.lt.ssgs_init) then
   call prandtl_profile(kappa_tt)
!else
!   call scalars_sgs(theta,dTdx,dTdy,dTdz,Ds_opt2_t,kappa_tt,s_Tn_all_t,s_Beta_t,I_LM_t,I_MM_t,I_QN_t,I_NN_t,ds2_clips_t,I_LM_MM_init_t,I_QN_NN_init_t,sigma_theta)
!end if

call RHS_scalar_solver(theta,dTdx,dTdy,dTdz,RHS_T,t_flux,kappa_tt)

call step_scalar(theta,RHS_T,RHS_Tf,dTdz_top)

$if($MPI) 
   call mpi_sendrecv (theta(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+3, &
                      theta(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+3, &
                      comm, status, ierr)
   call mpi_sendrecv (theta(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+4, &
                      theta(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+4, &
                     comm, status, ierr)
$endif

!if (coord==6) write(*,*) 'jt coord, u(25,25,nz-1)', jt, coord, u(25,25,nz-1)
!if (coord==6) write(*,*) 'jt coord, u(25,25,nz)', jt, coord, u(25,25,nz)
!if (coord==7) write(*,*) 'jt coord, u(25,25,0)', jt, coord, u(25,25,0)
!if (coord==7) write(*,*) 'jt coord, u(25,25,1)', jt, coord, u(25,25,1)

end subroutine temperature_transport


!************************************************************************
subroutine salinity_transport ()
!************************************************************************
use types, only: rprec
use derivatives, only: filt_da, ddz_uv
use param, only: nx, ny, nz, ld, lbz, z_i, &
                 MPI_RPREC, down, up, comm, status, ierr, coord, & 
                 ssgs_init, nproc, dz, jt
use sgs_param, only : Nu_t, kappa_ts, Ds_opt2_s
use sim_param, only: u, v, w, sal, dSdx, dSdy, dSdz, RHS_S, RHS_Sf
use sim_param, only: sal_flux, sigma_sal
use scalars_base, only: salinity_surface_bc
use scalars_base, only: dSdz_top
use scalars_util, only: RHS_scalar_solver, step_scalar 
use ocean_base, only: schmidt_profile
use sgs_param, only: s_Tn_all_s, s_Beta_s, I_LM_s, I_MM_s, I_QN_s, I_NN_s, ds2_clips_s
use sgs_param, only: I_LM_MM_init_s, I_QN_NN_init_s

implicit none

integer :: jx, jy, jz

!Calculate derivatives of theta
call filt_da(sal,dSdx,dSdy,lbz)

call ddz_uv(sal,dSdz,lbz)

$if($MPI)
   if (coord==nproc-1) then
      do jy=1,ny
         do jx=1,nx
            dSdz(jx,jy,nz) = (1._rprec/dz)*(sal(jx,jy,nz)-sal(jx,jy,nz-1))
         end do
      end do
   end if
$else
      do jy=1,ny
         do jx=1,nx
            dSdz(jx,jy,nz) = (1._rprec/dz)*(sal(jx,jy,nz)-sal(jx,jy,nz-1))
         end do
      end do
$endif

$if($MPI) 
   call mpi_sendrecv (dSdz(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+5, &
                      dSdz(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+5, &
                      comm, status, ierr)
   call mpi_sendrecv (dSdz(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+6, &
                      dSdz(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+6, &
                      comm, status, ierr)
$endif

RHS_Sf = RHS_S

call salinity_surface_bc()

!if (jt.lt.ssgs_init) then
   call schmidt_profile(kappa_ts)
!else
!   call scalars_sgs(sal,dSdx,dSdy,dSdz,Ds_opt2_s,kappa_ts,s_Tn_all_s,s_Beta_s,I_LM_s,I_MM_s,I_QN_s,I_NN_s,ds2_clips_s,I_LM_MM_init_s,I_QN_NN_init_s,sigma_sal)
!end if

call RHS_scalar_solver(sal,dSdx,dSdy,dSdz,RHS_S,sal_flux,kappa_ts)

call step_scalar(sal,RHS_S,RHS_Sf,dSdz_top)

$if($MPI) 
   call mpi_sendrecv (sal(:,:,1),ld*ny, MPI_RPREC, down, tag_counter+7, &
                      sal(:,:,nz), ld*ny, MPI_RPREC, up, tag_counter+7, &
                      comm, status, ierr)
   call mpi_sendrecv (sal(:,:,nz-1),ld*ny, MPI_RPREC, up, tag_counter+8, &
                      sal(:,:,0), ld*ny, MPI_RPREC, down, tag_counter+8, &
                      comm, status, ierr)
$endif

!if (coord==6) write(*,*) '2.jt coord, u(25,25,nz-1)', jt, coord, u(25,25,nz-1)
!if (coord==6) write(*,*) '2.jt coord, u(25,25,nz)', jt, coord, u(25,25,nz)
!if (coord==7) write(*,*) '2.jt coord, u(25,25,0)', jt, coord, u(25,25,0)
!if (coord==7) write(*,*) '2.jt coord, u(25,25,1)', jt, coord, u(25,25,1)

end subroutine salinity_transport


!************************************************************************
subroutine buoyancy_force ()
!************************************************************************
!This subroutine calculates the buoyancy term due to temperature  and salinity
!to be added to the RHS of the vertical momentum equation. 

use types, only: rprec
use param, only: z_i, u_star, ld, nx, ny, nz, lbz, coord, jt, dz
use scalars_base, only : alpha_theta, beta_sal, T_scale, S_scale, &
                         passive_temperature, passive_salinity, &
                         THETA_FLAG, SAL_FLAG
use sim_param, only: theta, sal, buoyancy

implicit none

integer :: i,j,k,jz_min
real(rprec), parameter :: g=9.81_rprec
real(rprec) :: g_hat,t_above,t_below, s_above, s_below
real(rprec), dimension(1:nz) :: buoy_avg
real(rprec) :: z
integer :: jx,jy,jz,jz_abs

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

if ( THETA_FLAG.and.(.not.SAL_FLAG) ) then
   if ( .not.passive_temperature ) then
      do k=jz_min, nz-1
         do j=1,ny
            do i=1,nx
               t_above=theta(i,j,k) 
               t_below=theta(i,j,k-1) 
               buoyancy(i,j,k) = g_hat*(1._rprec - alpha_theta*((t_above+t_below)/2._rprec)*T_scale )
            end do
         end do
      end do
   end if
else if ( SAL_FLAG.and.(.not.THETA_FLAG) ) then
        if ( .not.passive_salinity ) then  
           do k=jz_min, nz-1
              do j=1,ny
                 do i=1,nx
                    s_above=sal(i,j,k) 
                    s_below=sal(i,j,k-1) 
                    buoyancy(i,j,k) = g_hat*(1._rprec + beta_sal*((s_above+s_below)/2._rprec)*S_scale )
                 end do
              end do
           end do
        end if
else if ( THETA_FLAG.and.SAL_FLAG ) then
        if ( (.not.passive_temperature).and.(passive_salinity) ) then
           do k=jz_min, nz-1
              do j=1,ny
                 do i=1,nx
                    t_above=theta(i,j,k) 
                    t_below=theta(i,j,k-1) 
                    buoyancy(i,j,k) = g_hat*(1._rprec - alpha_theta*((t_above+t_below)/2._rprec)*T_scale )
                 end do
              end do
           end do
        else if ( (passive_temperature).and.(.not.passive_salinity) ) then
                do k=jz_min, nz-1
                   do j=1,ny
                      do i=1,nx
                         s_above=sal(i,j,k) 
                         s_below=sal(i,j,k-1) 
                         buoyancy(i,j,k) = g_hat*(1._rprec + beta_sal*((s_above+s_below)/2._rprec)*S_scale )
                      end do
                   end do
                end do
        else if ( (.not.passive_temperature).and.(.not.passive_salinity) ) then
                do k=jz_min, nz-1
                   do j=1,ny
                      do i=1,nx
                         t_above=theta(i,j,k) 
                         t_below=theta(i,j,k-1) 
                         s_above=sal(i,j,k) 
                         s_below=sal(i,j,k-1) 
                         buoyancy(i,j,k) = g_hat*(1._rprec - alpha_theta*(((t_above+t_below)/2._rprec)*T_scale - T_scale) &
                                              + beta_sal*(((s_above+s_below)/2._rprec)*S_scale - S_scale))
                      end do
                   end do
                end do
        end if
end if

!do jz=1,nz
!   $if($MPI)
!   jz_abs = coord * (nz-1) + jz
!   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
!   $else
!   jz_abs = jz
!   z = (jz - 0.5_rprec) * dz * z_i
!   $endif
!   buoy_avg(jz) = 0.0_rprec
!   do jy=1,ny
!      do jx=1,nx
!         buoy_avg(jz) = buoy_avg(jz) + buoyancy(jx,jy,jz)
!      end do
!   end do
!buoy_avg(jz) = buoy_avg(jz)/real(Nx*Ny)
!write(*,*) jz_abs, z, buoy_avg(jz)
!end do

!if (coord==6) write(*,*) 'jt, coord, buoyancy(25,25,nz-1)', jt, coord, buoyancy(25,25,nz-1)
!if (coord==6) write(*,*) 'jt, coord, buoyancy(25,25,nz)', jt, coord, buoyancy(25,25,nz)
!if (coord==7) write(*,*) 'jt, coord, buoyancy(25,25,0)', jt, coord, buoyancy(25,25,0)
!if (coord==7) write(*,*) 'jt, coord, buoyancy(25,25,1)', jt, coord, buoyancy(25,25,1)

end subroutine buoyancy_force


!************************************************************************
subroutine buoyancy_force_gsw ()
!************************************************************************
!This subroutine calculates the buoyancy term using the gsw formulation. 
!Temperature and salinity must both be used. 

use types, only: rprec
use param, only: z_i, u_star, ld, nx, ny, nz, lbz, coord, jt, dz
use scalars_base, only : alpha_theta, beta_sal, T_scale, S_scale, &
                         passive_temperature, passive_salinity, &
                         THETA_FLAG, SAL_FLAG, rho_0, longitude, &
                         latitude
use sim_param, only: theta, sal, buoyancy, pressure_z

implicit none

integer :: i,j,k,jz_min
real(rprec), parameter :: g=9.81_rprec
real(rprec) :: g_hat,t_above,t_below, s_above, s_below
real(rprec) :: p_z_above, p_z_below
real(rprec) :: rho, psal, sa, t, p_z, p_ref,ct,pt
real(rprec) :: gsw_pot_rho_t_exact, gsw_ct_from_pt,&
               gsw_t_from_ct, gsw_sa_from_sp
real(rprec), dimension(1:nz) :: buoy_avg
real(rprec) :: z
integer :: jx,jy,jz,jz_abs

!Non-dimensionalize gravity
g_hat = g*(z_i/(u_star**2))

!Reference pressure
p_ref = 0._rprec

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

do k=jz_min, nz-1
   do j=1,ny
      do i=1,nx
         t_above=theta(i,j,k) 
         t_below=theta(i,j,k-1) 
         s_above=sal(i,j,k) 
         s_below=sal(i,j,k-1) 
         p_z_above=pressure_z(i,j,k)
         p_z_below=pressure_z(i,j,k-1)

         pt = ((t_above+t_below)/2._rprec)*T_scale -273.15_rprec
         psal = ((s_above+s_below)/2._rprec)*S_scale
         p_z = (p_z_above+p_z_below)/2._rprec
         sa = gsw_sa_from_sp(psal, p_z, longitude, latitude)
         ct = gsw_ct_from_pt(sa,pt)
         t = gsw_t_from_ct(sa,ct,p_z)
         rho = gsw_pot_rho_t_exact(sa,t,p_z,p_ref) 

         buoyancy(i,j,k) = (1._rprec-((rho/rho_0) - 1._rprec))*g_hat
         end do
      end do
end do

!do jz=1,nz
!   $if($MPI)
!   jz_abs = coord * (nz-1) + jz
!   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
!   $else
!   jz_abs = jz
!   z = (jz - 0.5_rprec) * dz * z_i
!   $endif
!   buoy_avg(jz) = 0.0_rprec
!   do jy=1,ny
!      do jx=1,nx
!         buoy_avg(jz) = buoy_avg(jz) + buoyancy(jx,jy,jz)
!      end do
!   end do
!buoy_avg(jz) = buoy_avg(jz)/real(Nx*Ny)
!write(*,*) jz_abs, z, buoy_avg(jz)
!end do

!if (coord==6) write(*,*) 'jt, coord, buoyancy(25,25,nz-1)', jt, coord, buoyancy(25,25,nz-1)
!if (coord==6) write(*,*) 'jt, coord, buoyancy(25,25,nz)', jt, coord, buoyancy(25,25,nz)
!if (coord==7) write(*,*) 'jt, coord, buoyancy(25,25,0)', jt, coord, buoyancy(25,25,0)
!if (coord==7) write(*,*) 'jt, coord, buoyancy(25,25,1)', jt, coord, buoyancy(25,25,1)


end subroutine buoyancy_force_gsw

end module scalars_transport
