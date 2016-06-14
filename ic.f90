!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

subroutine ic()
  use types,only:rprec
  use param
  use sim_param, only : u, v, w
  use messages, only : error
#ifdef PPTURBINES
  use turbines, only: turbine_vel_init
#endif 
  
  implicit none
  
#ifdef PPVERBOSE
  character(*), parameter :: sub_name = 'ic'
#endif

#ifdef PPTURBINES
    real(rprec) :: zo_turbines
#endif    
  
  integer :: jx, jy, jz
  integer :: jz_abs

  real(rprec),dimension(nz) :: ubar
  real(rprec) :: rms, sigma_rv, noise, arg, arg2
  real(rprec) :: z

#ifdef PPTURBINES
  zo_turbines = 0._rprec
#endif

#ifdef PPCPS  
    call boundary_layer_ic()
#else
    if ( inflow ) then  !--no turbulence
        call uniform_ic()
    else
        call boundary_layer_ic()
    end if
#endif

  ! Display the mean vertical profiles of the initialized variables to std output
  do jz = lbz, nz
#ifdef PPMPI
     z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#else
     z = (jz - 0.5_rprec) * dz
#endif
     write(6,7780) jz,z,sum(u(1:nx,:,jz))/float(nx*ny),sum(v(1:nx,:,jz))/&
          float(nx*ny),sum(w(1:nx,:,jz))/float(nx*ny)
  end do
7780 format('jz,z,ubar,vbar,wbar:',(1x,I3,4(1x,F9.4)))

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine uniform_ic()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

u = inflow_velocity 
v = 0._rprec
w = 0._rprec

return
end subroutine uniform_ic

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine boundary_layer_ic()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! This subroutine produces an initial condition for the boundary layer.
! A log profile is used that flattens at z=z_i. Noise is added to 
! promote the generation of turbulence
!
#ifdef PPMPI
use mpi_defs
#endif
implicit none

if( coord == 0 ) write(*,*) '------> Creating modified log profile for IC'

do jz = 1, nz
#ifdef PPMPI
    z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#else
    z = (jz - 0.5_rprec) * dz
#endif

    ! IC in equilibrium with rough surface (rough dominates in effective zo)
    arg2 = z/zo
    arg = (1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

#ifdef PPLVLSET
    ! Kludge to adjust magnitude of velocity profile
    ! Not critical - may delete
    arg = 0.357*arg
#endif

#ifdef PPTURBINES
    call turbine_vel_init (zo_turbines)
    arg2 = z/zo_turbines
    arg = (1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z          
#endif        

    if (coriolis_forcing) then
        if (z > 0.5_rprec) then
            ubar(jz) = ug
        else
            ubar(jz) = arg/30._rprec
        end if
    else
        ubar(jz) = arg
    end if
end do

rms = 3._rprec
sigma_rv = 0.289_rprec

! Fill u, v, and w with uniformly distributed random numbers between 0 and 1
call init_random_seed
call random_number(u)
call random_number(v)
call random_number(w)

! Center random number about 0 and rescale
u = rms / sigma_rv * (u - 0.5_rprec)
v = rms / sigma_rv * (v - 0.5_rprec)
w = rms / sigma_rv * (w - 0.5_rprec)

! Rescale noise depending on distance from wall and mean log profile
do jz = 1, nz
#ifdef PPMPI
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i    !dimensions in meters
#else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i                        !dimensions in meters
#endif
    if (z <= z_i) then
        u(:,:,jz) = u(:,:,jz) * (1._rprec-z / z_i) + ubar(jz)
        v(:,:,jz) = v(:,:,jz) * (1._rprec-z / z_i)
        w(:,:,jz) = w(:,:,jz) * (1._rprec-z / z_i)
    else
        u(:,:,jz) = u(:,:,jz) * 0.01_rprec + ubar(jz)
        v(:,:,jz) = v(:,:,jz) * 0.01_rprec
        w(:,:,jz) = w(:,:,jz) * 0.01_rprec
    end if
end do

! Bottom boundary conditions
if (coord == 0) then
   w(:, :, 1) = 0._rprec
#ifdef PPMPI
   u(:, :, 0) = 0._rprec
   v(:, :, 0) = 0._rprec
   w(:, :, 0) = 0._rprec
#endif
end if

! Set upper boundary condition as zero gradient in u and v and no penetration in w
#ifdef PPMPI
if (coord == nproc-1) then
#endif    
   w(1:nx, 1:ny, nz) = 0._rprec
   u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
   v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
#ifdef PPMPI
end if
#endif

! Exchange ghost node information for u, v, and w
#ifdef PPMPI
call mpi_sync_real_array( u, lbz, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( v, lbz, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( w, lbz, MPI_SYNC_DOWNUP )
#endif

return
end subroutine boundary_layer_ic

end subroutine ic
