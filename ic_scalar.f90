module ic_scalar

use types, only: rprec
use param, only: ld, nx, ny, nz, z_i, dz, lbz, &
                 nproc, coord, nz_tot, dz
use sim_param, only: u,v,w
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use ocean_base, only: ml_base_on_grid

implicit none

contains

!********************************************************
subroutine read_init_scalar(fname_scalar, scalar, pressure, z1, scalar1)
!********************************************************
use functions, only: linear_interp, linreg

character(64), intent(in) :: fname_scalar
real(rprec), allocatable, dimension(:) :: temp_scalar
real(rprec), allocatable, dimension(:) :: temp_z
real(rprec), allocatable, dimension(:) :: temp_pressure
integer :: pair(2)
integer :: io, unit
integer :: i, n
real(rprec), dimension(nz) :: z 
integer :: jz, jz_abs
real(rprec), dimension(lbz:nz), intent(out) :: scalar, pressure
real(rprec), intent(out) :: z1, scalar1
real(rprec), dimension(1:2) :: reg_coeffs
real(rprec) :: m, b

!Read from file
unit = 12
open(unit, file=fname_scalar)
n=0

do 
  read(unit,*,iostat=io) pair
  if (io/=0) exit
  n=n+1
end do

rewind(unit)

allocate(temp_scalar(n))
allocate(temp_z(n))
allocate(temp_pressure(n))

do i=1,n
   read(unit,*) temp_scalar(i), temp_z(i), temp_pressure(i)
!   if (coord==0) print*, i, temp_scalar(i), temp_z(i), temp_pressure(i)  
end do
close(unit)

!print *, 'nz_tot',nz_tot
!print*, 'dz=', dz
!if (coord==0) then 
!   do i=1,n
!      print *, i, temp_scalar(i)
!      print *, i, temp_z(i)
!   end do
!end if

z1 = temp_z(1)
scalar1 = temp_scalar(1)

!Interpolate to vertical grid
do jz = 1,nz
   $if($MPI)
      z(jz) = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
   $else
      z(jz) = (jz - 0.5_rprec) * dz * z_i
   $endif
   !print*, 'coord, jz, z(jz)', coord, jz, z(jz)
end do

!Initialize scalar and pressure
scalar = 0.0_rprec
pressure = 0.0_rprec

!Linear regression to extrapolate pressure values
reg_coeffs=linreg(temp_z,temp_pressure,n)
m = reg_coeffs(1) !Slope of regression
b = reg_coeffs(2) !y-intercept of regression

do jz=1, nz
   if (z(jz).le.temp_z(1)) then 
      scalar(jz) = temp_scalar(1)
      pressure(jz) = m*z(jz) + b 
      if (pressure(jz).lt.0._rprec) pressure(jz) = 0._rprec
   else
      do i=1,n-1
         if ( (temp_z(i).lt.z(jz)) .and. (temp_z(i+1).ge.z(jz)) ) then
              scalar(jz) = linear_interp(temp_scalar(i), temp_scalar(i+1), temp_z(i+1)-temp_z(i), z(jz)-temp_z(i))
              pressure(jz) = linear_interp(temp_pressure(i), temp_pressure(i+1), temp_z(i+1)-temp_z(i), z(jz)-temp_z(i))
         end if
      end do
   end if
   !write(*,*) 'coord, z(jz), pressure(jz)', coord, z(jz), pressure(jz)
end do

end subroutine read_init_scalar


!********************************************************
subroutine initialize_temperature()
!********************************************************
use param, only: initt
use param, only: MPI_RPREC, comm, ierr
use sim_param, only: theta, pressure_z, sal
use scalars_base, only: dTdz_top, T_scale, S_scale, &
                        longitude, latitude

character(64) :: fname_init_theta
real(rprec) :: z1, theta_i1
real(rprec), dimension(lbz:nz) :: theta_i, pressure_i
real(rprec) :: z
integer :: jx, jy, jz, jz_abs
real(rprec), dimension(lbz:nz) :: thetai_avg, sali_avg
real(rprec) :: sa, psal, p_ref
real(rprec) :: gsw_sa_from_sp, gsw_pt_from_t

if (initt) then 

   fname_init_theta = "theta_itp77_1447"
   call read_init_scalar(fname_init_theta, theta_i, pressure_i, z1, theta_i1)

   p_ref = 0.0_rprec !Reference sea pressure 

   !Convert in-situ temperature to conservative temperature
   do jz=lbz,nz
      $if($MPI)
      z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
      $else
      z = (jz - 0.5_rprec) * dz * z_i
      $endif
      do jy=1,ny
         do jx=1,nx
            if (z.le.z1) then
                pressure_z(jx,jy,jz) = pressure_i(jz)
                psal = sal(jx,jy,jz)*S_scale
                sa = gsw_sa_from_sp(psal, pressure_i(jz), longitude, latitude)
                theta(jx,jy,jz) = (gsw_pt_from_t(sa,theta_i1,pressure_i(jz),p_ref) + 273.15_rprec)/T_scale
            else
                pressure_z(jx,jy,jz) = pressure_i(jz) 
                psal = sal(jx,jy,jz)*S_scale
                sa = gsw_sa_from_sp(psal, pressure_i(jz), longitude, latitude)
                theta(jx,jy,jz) = (gsw_pt_from_t(sa,theta_i(jz),pressure_i(jz),p_ref) + 273.15_rprec)/T_scale
            end if
         end do
       end do  
  end do

  !Top boundary condition
  if (coord==nproc-1) then
     dTdz_top = (theta_i(nz)/T_scale - theta_i(nz-1)/T_scale)/(dz*z_i)
  end if
  call mpi_bcast(dTdz_top, 1, MPI_RPREC, nproc-1, comm, ierr)

else
   call init_temperature
end if

$if ($MPI)
  call mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP )
  call mpi_sync_real_array(pressure_z, 0, MPI_SYNC_DOWNUP )
$endif

do jz=1,nz
   $if($MPI)
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
   $else
   jz_abs = jz
   z = (jz - 0.5_rprec) * dz * z_i
   $endif
   sali_avg(jz) = 0.0_rprec
   thetai_avg(jz) = 0.0_rprec
   do jy=1,ny
      do jx=1,nx
         sali_avg(jz) = sali_avg(jz) + sal(jx,jy,jz)
         thetai_avg(jz) = thetai_avg(jz) + theta(jx,jy,jz)
      end do
   end do
sali_avg(jz) = sali_avg(jz)/real(Nx*Ny)
thetai_avg(jz) = thetai_avg(jz)/real(Nx*Ny)
!write(*,*) jz_abs, z, sali_avg(jz)
end do

!if (coord==8) write(*,*) 'B. coord, pressure_z(25,25,nz)', coord, pressure_z(25,25,nz)
!if (coord==8) write(*,*) 'B. coord, pressure_z(25,25,nz-1)', coord, pressure_z(25,25,nz-1)
!if (coord==9) write(*,*) 'B. coord, pressure_z(25,25,0)', coord, pressure_z(25,25,0)
!if (coord==9) write(*,*) 'B. coord, pressure_z(25,25,1)', coord, pressure_z(25,25,1)

end subroutine initialize_temperature


!********************************************************
subroutine init_temperature()
!********************************************************
use sim_param, only :theta
use scalars_base, only: T_mean, T_scale, dTdz_top, &
                        inv_strength

integer :: jx, jy, jz, jz_abs
integer :: seed
real(rprec) :: z, z_ml
real(rprec) :: rms, noise, sigma_rv 

!Find the grid level closest to the base of the thermocline
call ml_base_on_grid(z_ml)

rms=0.1_rprec
sigma_rv = 0.289
do jz=1,nz
   $if ($MPI)
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
   $else
   jz_abs = jz
   z = (jz - 0.5_rprec) * dz * z_i
   $endif
   seed = -80 - jz_abs
   call random_seed(seed)
   do jy=1,ny
      do jx=1,nx
         if (z.le.z_ml) then
            !call random_number(noise)
            theta(jx,jy,jz)=T_mean/T_scale
         else
            !call random_number(noise)
            theta(jx,jy,jz)=(T_mean+ (z-z_ml)*inv_strength)/T_scale
         end if
      end do
   end do
   !write(*,*) jz_abs, z, theta(25,25,jz)
end do

!Upper boundary condition for theta
dTdz_top = inv_strength/T_scale
$if($MPI)
    if (coord == nproc-1) then
       theta(:,:,nz) = theta(:,:,nz-1) + dTdz_top*z_i*dz
    end if
$else
      theta(:,:,nz) = theta(:,:,nz-1) + dTdz_top*z_i*dz
$endif

return

end subroutine init_temperature


!********************************************************
subroutine initialize_salinity()
!********************************************************
use param, only: inits
use param, only: MPI_RPREC, comm, ierr
use sim_param, only: sal, pressure_z
use scalars_base, only: dSdz_top, S_scale, longitude, &
                        latitude 

character(64) :: fname_init_sal
real(rprec) :: z1, sal_i1
real(rprec), dimension(lbz:nz) :: sal_i, pressure_i
real(rprec) :: z
integer :: jx, jy, jz, jz_abs
real(rprec), dimension(lbz:nz) :: sali_avg

if (inits) then 

   fname_init_sal = "sal_itp77_1447"
   call read_init_scalar(fname_init_sal, sal_i, pressure_i, z1, sal_i1)
   
   !Convert practical salinity to absolute salinity
   do jz=lbz,nz
      $if($MPI)
      z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
      $else
      z = (jz - 0.5_rprec) * dz * z_i
      $endif
      do jy=1,ny
         do jx=1,nx
                if (z.le.z1) then
                   sal(jx,jy,jz) = sal_i1/S_scale
                else    
                   sal(jx,jy,jz) = sal_i(jz)/S_scale 
                end if
         end do
       end do  
  end do

  !Top boundary condition
  if (coord==nproc-1) then
     dSdz_top = (sal_i(nz)/S_scale - sal_i(nz-1)/S_scale)/(dz*z_i)
  end if
  call mpi_bcast(dSdz_top, 1, MPI_RPREC, nproc-1, comm, ierr)

else
   call init_salinity
end if

$if ($MPI)
  call mpi_sync_real_array(sal, 0, MPI_SYNC_DOWNUP )
  call mpi_sync_real_array(pressure_z, 0, MPI_SYNC_DOWNUP )
$endif

do jz=1,nz
   $if($MPI)
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
   $else
   jz_abs = jz
   z = (jz - 0.5_rprec) * dz * z_i
   $endif
   sali_avg(jz) = 0.0_rprec
   do jy=1,ny
      do jx=1,nx
         sali_avg(jz) = sali_avg(jz) + sal(jx,jy,jz)
      end do
   end do
sali_avg(jz) = sali_avg(jz)/real(Nx*Ny)
!write(*,*) jz_abs, z, sali_avg(jz)
end do

!if (coord==3) write(*,*) 'B. coord, pressure_z(25,25,nz)', coord, pressure_z(25,25,nz)
!if (coord==3) write(*,*) 'B. coord, pressure_z(25,25,nz-1)', coord, pressure_z(25,25,nz-1)
!if (coord==4) write(*,*) 'B. coord, pressure_z(25,25,0)', coord, pressure_z(25,25,0)
!if (coord==4) write(*,*) 'B. coord, pressure_z(25,25,1)', coord, pressure_z(25,25,1)

end subroutine initialize_salinity


!********************************************************
subroutine init_salinity()
!********************************************************
use sim_param, only: sal
use scalars_base, only: S_mean, S_scale, dSdz_top, &
                        inv_strength

integer :: jx, jy, jz, jz_abs
integer :: seed
real(rprec) :: z, z_ml
real(rprec) :: rms, noise, sigma_rv 

!Find the grid level closest to the base of the thermocline
call ml_base_on_grid(z_ml)

rms=0.1_rprec
sigma_rv = 0.289
do jz=1,nz
   $if ($MPI)
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
   $else
   jz_abs = jz
   z = (jz - 0.5_rprec) * dz * z_i
   $endif
   seed = -80 - jz_abs
   call random_seed(seed)
   do jy=1,ny
      do jx=1,nx
         if (z.le.z_ml) then
            !call random_number(noise)
            sal(jx,jy,jz)=S_mean/S_scale
         else
            !call random_number(noise)
            sal(jx,jy,jz)=(S_mean - (z-z_ml)*inv_strength)/S_scale
         end if
      end do
   end do
   !write(*,*) jz_abs, z, sal(25,25,jz)
end do
!if (coord==0) write(*,*) 'z_ml=', z_ml


!Upper boundary condition for theta
dSdz_top = inv_strength/S_scale
$if($MPI)
    if (coord == nproc-1) then
       sal(:,:,nz) = sal(:,:,nz-1) + dSdz_top*z_i*dz
    end if
$else
      sal(:,:,nz) = sal(:,:,nz-1) + dSdz_top*z_i*dz
$endif

return

end subroutine init_salinity


end module ic_scalar


