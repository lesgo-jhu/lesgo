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

!*******************************************************************************
subroutine initial()
!*******************************************************************************
use iwmles
use types,only:rprec
use param
use sim_param, only : u,v,w,RHSx,RHSy,RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
#ifdef PPDYN_TN
use sgs_param, only:F_ee2,F_deedt2,ee_past
#endif
#ifdef PPTURBINES
use sim_param,only:fxa,fya,fza
#endif
#ifdef PPLVLSET
use sim_param,only:fxa,fya,fza
use sim_param,only:fx,fy,fz
#endif
use string_util, only : string_concat
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#endif

implicit none

character (64) :: fname

#ifdef PPDYN_TN
logical :: exst
character (64) :: fname_dyn_tn
#endif

integer::jz

! Flag to identify if file exists
logical :: file_flag
logical :: iwm_file_flag !xiang: for iwm restart

#ifdef PPTURBINES
fxa=0._rprec
fya=0._rprec
fza=0._rprec
#endif

#ifdef PPLVLSET
fx=0._rprec;fy=0._rprec;fz=0._rprec
fxa=0._rprec; fya=0._rprec; fza=0._rprec
#endif

#ifdef PPDYN_TN
!Will be over-written if read from dyn_tn.out files
ee_past = 0.1_rprec; F_ee2 = 10.0_rprec; F_deedt2 = 10000.0_rprec
fname_dyn_tn = path // 'dyn_tn.out'
#ifdef PPMPI
  call string_concat( fname_dyn_tn, '.c', coord )
#endif
#endif

fname = checkpoint_file

#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif

!check iwm restart file
inquire (file='iwm_checkPoint.dat', exist=iwm_file_flag)
if (iwm_file_flag) then
    if (lbc_mom == 3) then
        if (coord==0) call iwm_read_checkPoint()
    end if
end if

! Check if file exists and read from it
! if not then create new IC
inquire( file=fname, exist=file_flag )
if (file_flag) then
    initu = .true.
    inilag = .false.
else
    initu = .false.
    inilag = .true.
end if

if (initu) then
    if (coord == 0) write(*,*) '--> Reading initial velocity field from file'
    call ic_file
#ifndef PPCPS
else if (inflow) then
        if (coord == 0) write(*,*) '--> Creating initial uniform velocity field'
        call ic_uniform
#endif
else if (lbc_mom==1) then
    if (coord == 0) write(*,*) '--> Creating initial boundary layer velocity ',&
        'field with DNS BCs'
    call ic_dns()
else
    if (coord == 0) write(*,*) '--> Creating initial boundary layer velocity ',&
    'field with LES BCs'
    call ic_les()
end if

#ifdef PPDYN_TN
! Read dynamic timescale running averages from file
if (cumulative_time) then
    inquire (file=fname_dyn_tn, exist=exst)
    if (exst) then
        open(13,file=fname_dyn_tn,form='unformatted', convert=read_endian)
        read(13) F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)
    else
        write(*,*) trim(fname_dyn_tn), ' not found - using default values'
    end if
    close(13)
end if
#endif

! Write averaged vertical profiles to standard output
do jz=1,nz
    write(6,7780) jz, sum (u(1:nx, :, jz)) / (nx * ny),                        &
                  sum (v(1:nx, :, jz)) / (nx * ny),                            &
                  sum (w(1:nx, :, jz)) / (nx * ny)
end do
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))

#ifdef PPMPI
! Exchange ghost node information for u, v, and w
call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP ) 
call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP ) 

!--set 0-level velocities to BOGUS
if (coord == 0) then
    u(:, :, lbz) = BOGUS
    v(:, :, lbz) = BOGUS
    w(:, :, lbz) = BOGUS
end if
#endif

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ic_uniform()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine creates a uniform initial condition without turbulence.
!
implicit none

u = inflow_velocity 
v = 0._rprec
w = 0._rprec

end subroutine ic_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ic_file()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine reads the initial conditions from the checkpoint file.
!
open(12, file=fname, form='unformatted', convert=read_endian)

read(12) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),                          &
         RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),                 &
         Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),                    &
         F_QN(:,:,1:nz), F_NN(:,:,1:nz)   
         
close(12)

end subroutine ic_file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ic_dns()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine produces an initial condition for the boundary layer when 
! using DNS boundary conditions.
!
use types,only:rprec
use param
use sim_param,only:u,v,w
implicit none

real(rprec),dimension(nz)::ubar
real(rprec)::rms,temp
integer::jx,jy,jz,z
real(rprec) :: dummy_rand

! Calculate the average streamwise velocity based on height of first uvp point 
! in wall units
do jz=1,nz
    z = int((real(jz)-.5_rprec)*dz) ! non-dimensional
    ubar(jz) = (u_star*z_i/nu_molec)*z*(1._rprec-.5_rprec*z) ! non-dimensional
end do

! Get random seeds to populate the initial condition with noise
call init_random_seed

! Add noise to the velocity field
! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms=0.2_rprec
do jz=1,nz
    do jy=1,ny
        do jx=1,nx
            call random_number(dummy_rand)
            u(jx,jy,jz)=ubar(jz)+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
            call random_number(dummy_rand)
            v(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
            call random_number(dummy_rand)
            w(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
        end do
    end do
end do

! make sure w-mean is 0
temp=0._rprec
do jz=1,nz
    do jy=1,ny
        do jx=1,nx
            temp=temp+w(jx,jy,jz)
        end do
    end do
end do
temp=temp/(nx*ny*nz)

do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         w(jx,jy,jz)=w(jx,jy,jz)-temp
      end do
   end do
end do

! Make sure field satisfies boundary conditions
w(:,:,1)=0._rprec
w(:,:,nz)=0._rprec
u(:,:,nz)=u(:,:,nz-1)
v(:,:,nz)=v(:,:,nz-1)

end subroutine ic_dns

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ic_les()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine produces an initial condition for the boundary layer.
! A log profile is used that flattens at z=z_i. Noise is added to 
! promote the generation of turbulence
!
use types,only:rprec
use param
use sim_param, only : u, v, w
use messages, only : error
#ifdef PPTURBINES
use turbines, only: turbine_vel_init
#endif 

implicit none
integer :: jz, jz_abs
real(rprec), dimension(nz) :: ubar
real(rprec) :: rms, sigma_rv, arg, arg2, z

character(*), parameter :: sub_name = 'ic'

#ifdef PPTURBINES
    real(rprec) :: zo_turbines = 0._rprec
#endif

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

end subroutine ic_les

end subroutine initial
