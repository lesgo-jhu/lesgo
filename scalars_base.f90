module scalars_base

use types, only: rprec
$if ($MPI)
use mpi_defs
$endif

implicit none

real(rprec), parameter :: g=9.81_rprec
!Roughness length for scalars: (1/10)*momentum roughness length
real(rprec) :: zo_s = 0.00001_rprec
real(rprec) :: rho_0 = 1024._rprec !reference density

!THETA_FLAG=1 for tempearture, THETA_FLAG=0 for no temperature
logical, parameter :: THETA_FLAG = .true.
logical, parameter :: passive_temperature = .false.
real(rprec), parameter :: alpha_theta = -0.0001_rprec !thermal expansion coefficient of sea-water

!SAL_FLAG=1 for salinity, SAL_FLAG=0 for no salinity
logical, parameter :: SAL_FLAG=.true.
logical, parameter :: passive_salinity = .false.
real(rprec), parameter :: beta_sal = -0.0008_rprec !haline contraction coefficient of sea-water
!real(rprec), parameter :: beta_sal = 0.002_rprec !haline contraction coefficient of sea-water

! Temperature parameter
real(rprec), parameter :: inv_strength=0.010_rprec
real(rprec) :: dTdz_top
real(rprec), parameter :: T_scale=300._rprec
real(rprec), parameter :: T_mean = 290._rprec
real(rprec), parameter :: z_inv = 33._rprec
real(rprec) :: dSdz_top
real(rprec), parameter :: S_scale=30._rprec
real(rprec), parameter :: S_mean = 290._rprec

! Stability function in M-O similarity expressions
real(rprec), parameter :: nu_m = 1.84e-6_rprec
real(rprec), parameter :: kappa_tm = 1.38e-7_rprec
real(rprec), parameter :: kappa_sm = 9.0e-10_rprec
real(rprec), parameter ::pow = 2._rprec/3._rprec

! Location of ITP profile
real(rprec), parameter :: longitude = -150.4544_rprec
real(rprec), parameter :: latitude = 74.3316_rprec

contains

!************************************************************************
subroutine temperature_surface_bc ()
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, dz, vonk, u_b, &
                 lbz, zo, jt, jt_total,   &
                 coord, MPI_RPREC, comm, ierr, &
                 domain_nstart, domain_nend, domain_nskip, &
                 u_star, z_i
use sim_param, only: u, v, theta, tstar, t_flux, &
                     tstar_inst_avg, t_flux_inst_avg, &
                     sal_w
use test_filtermodule

implicit none

integer :: i, j

real(rprec),dimension(nx,ny)::ustar, u_avg
real(rprec),dimension(ld,ny)::u1,v1,theta1
real(rprec),dimension(ld,ny)::u_r1,v_r1
real(rprec) :: denom
real(rprec), dimension(ld,ny) :: phi_t
real(rprec), dimension(ld,ny) :: theta_f
real(rprec) :: m

!!Wall model for dsdz
!Calculate friction velocity ustar
u1=u(:,:,1)
v1=v(:,:,1)
u_r1 = u1 - u_b
v_r1 = v1
theta1=theta(:,:,1)
call test_filter ( u_r1 )
call test_filter ( v_r1 )
call test_filter ( theta1 ) 
denom=log(0.5_rprec*dz/zo)
u_avg=sqrt(u_r1(1:nx,1:ny)**2+v_r1(1:nx,1:ny)**2)
ustar=u_avg*vonk/denom

phi_t = 0.0_rprec
m = 0.054_rprec * S_scale/T_scale

do j=1,Ny
   do i=1,Nx
    theta_f(i,j) = (-(m*sal_w(i,j))*T_scale + 273.15_rprec)/T_scale
    phi_t(i,j) = (1._rprec/vonk)*log(0.5_rprec*dz/zo) + 1.57_rprec*sqrt(zo*z_i*ustar(i,j)*u_star/nu_m)*((nu_m/kappa_tm)**pow)
    t_flux(i,j) = ( (theta_f(i,j) - theta1(i,j))/phi_t(i,j) ) * ustar(i,j)
    tstar(i,j) = (theta_f(i,j) - theta1(i,j))/phi_t(i,j) 
   end do
end do
!if (coord==0) write(*,*) 'jt_total, sal_w(25,25), theta_f(25,25), t_flux(i,j)=', jt_total, sal_w(25,25), theta_f(25,25)*T_scale, t_flux(i,j)

!Calculate instantaneous average
t_flux_inst_avg = 0.0_rprec
tstar_inst_avg = 0.0_rprec

if (coord==0) then
   do j=1,Ny
      do i=1,Nx
         t_flux_inst_avg = t_flux_inst_avg + t_flux(i,j)
         tstar_inst_avg = tstar_inst_avg + tstar(i,j)
      end do
   end do
   t_flux_inst_avg = t_flux_inst_avg/(real(Nx*Ny))
   tstar_inst_avg = tstar_inst_avg/(real(Nx*Ny))
end if

call mpi_bcast(t_flux_inst_avg, 1, MPI_RPREC, 0 , comm, ierr)
call mpi_bcast(tstar_inst_avg, 1, MPI_RPREC, 0 , comm, ierr)

!if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and. ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
if (mod(jt,100)==0) then
   if (coord==0) write(*,*) 'Instantaneous average temperature flux=', jt_total, t_flux_inst_avg
   if (coord==0) write(*,*) 'Instantaneous average friction temperature=', jt_total, tstar_inst_avg
end if

end subroutine temperature_surface_bc



!************************************************************************
subroutine salinity_surface_bc ()
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, dz, vonk, u_b, u_star, z_i, &
                 lbz, zo, jt, jt_total,   &
                 coord, MPI_RPREC, comm, ierr, &
                 domain_nstart, domain_nend, domain_nskip
use sim_param, only: u, v, theta, sal, sstar, sal_flux, &
                     sal_flux_inst_avg, sstar_inst_avg, &
                     sal_w
use test_filtermodule

implicit none

integer :: i, j

real(rprec),dimension(nx,ny)::ustar, u_avg
real(rprec),dimension(ld,ny)::u1,v1,theta1,sal1
real(rprec),dimension(ld,ny) :: u_r1,v_r1
real(rprec) :: denom
real(rprec),dimension(ld,ny) :: phi_s, phi_t
real(rprec) :: sal_ice, L, m
real(rprec), dimension(nx,ny) :: A, B, C, T1
!real(rprec), dimension(nx,ny) :: sal_w1, sal_w2, sal_w
real(rprec), dimension(nx,ny) :: sal_w1, sal_w2
real(rprec), dimension(nx,ny) :: w_ice 

phi_s = 0.0_rprec
phi_t = 0.0_rprec
w_ice = 0.0_rprec

sal_ice = 3._rprec/S_scale
L = 67.9_rprec/T_scale
m = 0.054_rprec * S_scale/T_scale

!!Wall model for dsdz
!Calculate friction velocity ustar
u1=u(:,:,1)
v1=v(:,:,1)
u_r1 = u1 - u_b
v_r1 = v1
theta1 = theta(:,:,1)
sal1=sal(:,:,1) 
call test_filter ( u_r1)
call test_filter ( v_r1 )
call test_filter ( theta1 ) 
call test_filter ( sal1 ) 
denom=log(0.5_rprec*dz/zo)
u_avg=sqrt(u_r1(1:nx,1:ny)**2+v_r1(1:nx,1:ny)**2)
ustar=u_avg*vonk/denom

do j=1,Ny
   do i=1,Nx
    phi_t(i,j) = (1._rprec/vonk)*log(0.5_rprec*dz/zo) + 1.57_rprec*sqrt(zo*z_i*ustar(i,j)*u_star/nu_m)*((nu_m/kappa_tm)**pow)
    phi_s(i,j) = (1._rprec/vonk)*log(0.5_rprec*dz/zo) + 1.57_rprec*sqrt(zo*z_i*ustar(i,j)*u_star/nu_m)*((nu_m/kappa_sm)**pow)
   end do
end do

do j=1,Ny
   do i=1,Nx
      T1(i,j) = ((theta1(i,j)*T_scale)-273.15_rprec)/T_scale 
      A(i,j) = m
      B(i,j) = T1(i,j) - m*sal_ice + L*phi_t(i,j)/phi_s(i,j)
      C(i,j) = -( T1(i,j)*sal_ice + L*phi_t(i,j)*sal1(i,j)/phi_s(i,j) )
      sal_w1(i,j) = ( -B(i,j) + sqrt(B(i,j)**2._rprec - 4._rprec*A(i,j)*C(i,j)) )/(2._rprec*A(i,j))
      sal_w2(i,j) = ( -B(i,j) - sqrt(B(i,j)**2._rprec - 4._rprec*A(i,j)*C(i,j)) )/(2._rprec*A(i,j))
      if (sal_w1(i,j).gt.0._rprec) then
         sal_w(i,j) = sal_w1(i,j)
      elseif (sal_w2(i,j).gt.0._rprec) then
         sal_w(i,j) = sal_w2(i,j)
      end if
   end do
end do
!if (coord==0) write(*,*) 'sal_w1(25,25), sal_w2(25,25), sal_w(25,25)=', jt_total, sal_w1(25,25),sal_w2(25,25),sal_w(25,25)
!if (coord==0) write(*,*) 'A. jt_total, sal_w(25,25)=', jt_total, sal_w(25,25)

do j=1,Ny
   do i=1,Nx
      w_ice(i,j) = ustar(i,j)*(sal1(i,j)-sal_w(i,j))/( phi_s(i,j)*(sal_w(i,j)-sal1(i,j)) )
      sal_flux(i,j) = -w_ice(i,j) * (sal_w(i,j) - sal1(i,j))
   end do
end do

!if (coord==0) write(*,*) 'jt_total, w_ice(25,25), sal_flux(25,25)=', jt_total, w_ice(25,25), sal_flux(25,25)

!Calculate instantaneous average
sal_flux_inst_avg = 0.0_rprec

if (coord==0) then
   do j=1,Ny
      do i=1,Nx
         sal_flux_inst_avg = sal_flux_inst_avg + sal_flux(i,j)
      end do
   end do
   sal_flux_inst_avg = sal_flux_inst_avg/(real(Nx*Ny))
end if

call mpi_bcast(sal_flux_inst_avg, 1, MPI_RPREC, 0 , comm, ierr)

!if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and. ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
if (mod(jt,100)==0) then
   if (coord==0) write(*,*) 'Instantaneous average salinity flux=', jt_total, sal_flux_inst_avg
end if


end subroutine salinity_surface_bc

!************************************************************************
subroutine salinity_similarity_bc ()
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, dz, vonk, u_b, u_star, z_i, &
                 lbz, zo, jt, jt_total,   &
                 coord, MPI_RPREC, comm, ierr, &
                 domain_nstart, domain_nend, domain_nskip
use sim_param, only: u, v, sal, sstar, sal_flux, &
                     sal_flux_inst_avg, sstar_inst_avg
use test_filtermodule

implicit none

integer :: i, j

real(rprec),dimension(nx,ny)::ustar, u_avg
real(rprec),dimension(ld,ny)::u1,v1,sal1
real(rprec),dimension(ld,ny) :: u_r1,v_r1
real(rprec) :: denom
real(rprec),dimension(ld,ny) :: phi_s


!!Wall model for dsdz
!Calculate friction velocity ustar
u1=u(:,:,1)
v1=v(:,:,1)
u_r1 = u1 - u_b
v_r1 = v1
sal1=sal(:,:,1)
call test_filter ( u_r1 )
call test_filter ( v_r1 )
call test_filter ( sal1 ) 
denom=log(0.5_rprec*dz/zo)
u_avg=sqrt(u_r1(1:nx,1:ny)**2+v_r1(1:nx,1:ny)**2)
ustar=u_avg*vonk/denom

phi_s = 0.0_rprec

do j=1,Ny
   do i=1,Nx
    phi_s(i,j) = (1._rprec/vonk)*log(0.5_rprec*dz/zo) + 1.57_rprec*sqrt(zo*z_i*ustar(i,j)*u_star/nu_m)*((nu_m/kappa_sm)**pow)
    sal_flux(i,j) =  (5.0_rprec/30.0_rprec - sal1(i,j))/phi_s(i,j) * ustar(i,j)
    sstar(i,j) =  (5.0_rprec/30.0_rprec - sal1(i,j))/phi_s(i,j)
   end do
end do

!Calculate instantaneous average
sal_flux_inst_avg = 0.0_rprec
sstar_inst_avg = 0.0_rprec

if (coord==0) then
   do j=1,Ny
      do i=1,Nx
         sal_flux_inst_avg = sal_flux_inst_avg + sal_flux(i,j)
         sstar_inst_avg = sstar_inst_avg + sstar(i,j)
      end do
   end do
   sal_flux_inst_avg = sal_flux_inst_avg/(real(Nx*Ny))
   sstar_inst_avg = sstar_inst_avg/(real(Nx*Ny))
end if

call mpi_bcast(sal_flux_inst_avg, 1, MPI_RPREC, 0 , comm, ierr)
call mpi_bcast(sstar_inst_avg, 1, MPI_RPREC, 0 , comm, ierr)

!if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and. ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
if (mod(jt,100)==0) then
   if (coord==0) write(*,*) 'Instantaneous average salinity flux=', jt_total, sal_flux_inst_avg
   if (coord==0) write(*,*) 'Instantaneous average friction salinity=', jt_total, sstar_inst_avg
end if



end subroutine salinity_similarity_bc


end module scalars_base
