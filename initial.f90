!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
#ifdef PPTURBINES
use sim_param, only : fxa, fya, fza
#endif
#ifdef PPLVLSET
use sim_param, only : fxa, fya, fza
use sim_param, only : fx, fy, fz
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
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
#endif

#ifdef PPLVLSET
fx = 0._rprec; fy = 0._rprec; fz = 0._rprec
fxa = 0._rprec;  fya = 0._rprec; fza = 0._rprec
#endif

#ifdef PPDYN_TN
! Will be over-written if read from dyn_tn.out files
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

! check iwm restart file
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
    if (coord == 0) write(*,*) '--> Creating initial laminar profile ',&
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
do jz = 1, nz
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

!*******************************************************************************
subroutine ic_uniform()
!*******************************************************************************
! This subroutine creates a uniform initial condition without turbulence.
!
implicit none

u = inflow_velocity
v = 0._rprec
w = 0._rprec

end subroutine ic_uniform

!*******************************************************************************
subroutine ic_file()
!*******************************************************************************
! This subroutine reads the initial conditions from the checkpoint file.
!
open(12, file=fname, form='unformatted', convert=read_endian)

read(12) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),                          &
         RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),                 &
         Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),                    &
         F_QN(:,:,1:nz), F_NN(:,:,1:nz)

close(12)

end subroutine ic_file

!*******************************************************************************
subroutine ic_dns()
!*******************************************************************************
! This subroutine produces an initial condition for the boundary layer when
! using DNS boundary conditions.
!
use types,only:rprec
use param
use sim_param,only:u,v,w
use io, only:output_loop,energy
implicit none

logical :: visualize_initial

real(rprec), dimension(nz) :: ubar
real(rprec) :: temp, z
integer :: jx, jy, jz
real(rprec), dimension(:,:,:), allocatable :: upi,vpi,wpi
real(rprec) :: ke

integer :: pert_type
pert_type = 2
! pert_type = 1 --> random noise (fraction of ubar)
! pert_type = 2 --> localized perturbation (four vortical structures)

! Calculate the average streamwise velocity based on height of first uvp point
! in wall units
if (use_mean_p_force==.true.) then
if ( abs(ubot) > 0 .or. abs(utop) > 0 ) then
    !! linear laminar profile (couette) z and ubar are non-dimensional
    do jz = 1, nz
#ifdef PPMPI
        z = (coord*(nz-1) + real(jz,rprec) - 0.5_rprec) * dz
#else
        z = (real(jz,rprec) - 0.5_rprec) * dz ! non-dimensional
#endif
        ubar(jz)= (utop-ubot)/2*(2*z/L_z-1)**5 + (utop+ubot)/2
    end do
else
    !! parabolic laminar profile (channel) z and ubar are non-dimensional
    do jz = 1, nz
#ifdef PPMPI
        z = (coord*(nz-1) + real(jz,rprec) - 0.5_rprec) * dz
#else
        z = (real(jz,rprec) - 0.5_rprec) * dz
#endif
        ubar(jz) = (u_star*z_i/nu_molec) * z * (1._rprec - 0.5_rprec*z)        
    end do
endif
else
!! parabolic laminar profile based on ubulk=1 and constant mass flow assumed
do jz = 1, nz
    z = (coord*(nz-1) + real(jz,rprec) - 0.5_rprec) * dz
    ubar(jz) = 3._rprec*z - 1.5_rprec*z**2
end do
endif

! Gets perturbation to add to initial profile
call ic_pert(pert_type,upi,vpi,wpi,ubar)

! Add noise to the velocity field
do jz = 1, nz
    do jy = 1, ny
        do jx = 1, nx
            u(jx,jy,jz) = ubar(jz)+upi(jx,jy,jz)
            v(jx,jy,jz) = 0._rprec+vpi(jx,jy,jz)
            w(jx,jy,jz) = 0._rprec+wpi(jx,jy,jz)
        end do
    end do
end do

deallocate(upi,vpi,wpi)

!!!!!! to visualize initial velocity field. remove later.
visualize_initial = .true.
if (visualize_initial .eqv. .true.) then
    call output_loop()
    call energy(ke)
end if


! make sure w-mean is 0
temp=0._rprec
do jz = 1, nz
    do jy = 1, ny
        do jx = 1, nx
            temp = temp+w(jx,jy,jz)
        end do
    end do
end do
temp = temp/(nx*ny*nz)

do jz = 1, nz
   do jy = 1, ny
      do jx = 1, nx
         w(jx,jy,jz) = w(jx,jy,jz)-temp
      end do
   end do
end do

! Make sure field satisfies boundary conditions
w(:,:,1) = 0._rprec
w(:,:,nz) = 0._rprec
if (ubc_mom == 0) then
   u(:,:,nz) = u(:,:,nz-1)
   v(:,:,nz) = v(:,:,nz-1)
endif

end subroutine ic_dns


!*******************************************************************************
subroutine ic_pert(pert_type,upi,vpi,wpi,ubar)
!*******************************************************************************
! Inserts perturbation for DNS simulation.
! pert_type = 1 --> random noise (fraction of ubar)
! pert_type = 2 --> localized perturbation (four vortical structures)


use types,only:rprec
use param
use functions, only: interp_to_w_grid
implicit none

integer, intent(in) :: pert_type
real(rprec), dimension(:,:,:), allocatable, intent(inout) :: upi,vpi,wpi
real(rprec), dimension(:,:,:), allocatable :: wpi_uv
real(rprec), dimension(nz), intent(in) :: ubar
real(rprec) :: ubarfrac, dummy_rand
integer :: jx, jy, jz

real(rprec) :: theta, eps, lx, ly, p, q, xc, yc, zc, ubulk, ubulk_global
real(rprec), dimension(nx) :: x
real(rprec), dimension(ny) :: y
real(rprec), dimension(nz) :: z, f, dfdz, z_w
real(rprec), dimension(nx,ny) :: xp, yp, g, dgdy

allocate(upi(nx,ny,lbz:nz), vpi(nx,ny,lbz:nz), wpi(nx,ny,lbz:nz))
allocate(wpi_uv(nx,ny,lbz:nz))
upi(nx,ny,lbz:nz) = 0._rprec
vpi(nx,ny,lbz:nz) = 0._rprec
wpi(nx,ny,lbz:nz) = 0._rprec
wpi_uv(nx,ny,lbz:nz) = 0._rprec

! random noise
if (pert_type==1) then
    ubarfrac = 0.3_rprec
    call init_random_seed()
    do jz = 1, nz
        do jy = 1, ny
            do jx = 1, nx
                call random_number(dummy_rand)
                upi(jx,jy,jz) = ubarfrac*2._rprec*(dummy_rand-.5_rprec)*ubar(jz)/u_star
                call random_number(dummy_rand)
                vpi(jx,jy,jz) = ubarfrac*2._rprec*(dummy_rand-.5_rprec)*ubar(jz)/u_star
                call random_number(dummy_rand)
                wpi(jx,jy,jz) = ubarfrac*2._rprec*(dummy_rand-.5_rprec)*ubar(jz)/u_star
            end do
        end do
    end do

! localized perturbation
elseif (pert_type==2) then
    theta = 0._rprec ! rotates localized perturbation if not zero
    lx = 2._rprec ! new parameter, not L_x
    ly = 2._rprec ! new parameter, not L_y
    p = 2._rprec ! how quickly perturbation decays near the wall
    q = 2._rprec ! how quickly perturbation decays near the wall
    xc = L_x/2._rprec
    yc = L_y/2._rprec
    zc = L_z/2._rprec

    ! set up grid
    ! x(1) and y(1) are dx and dy, respectively
    do jx = 1, nx
        x(jx) = real(jx,rprec) * dx
    end do
    do jy = 1, ny
        y(jy) = real(jy,rprec) * dy
    end do
    ! z on uv-grid
    do jz = 1, nz
        z(jz) = (coord*(nz-1) + real(jz,rprec) - 0.5_rprec) * dz
    end do
    
    ! z_w on w-grid
    do jz = 1, nz
        z_w(jz) = (coord*(nz-1) + real(jz,rprec) - 1._rprec) * dz
    end do

    ! calculate ubulk
    ubulk = 0.0_rprec
    do jz=1, nz-1
        ubulk = ubulk + ubar(jz)*dz
    end do
    ubulk = ubulk / L_z
    call mpi_allreduce(ubulk,ubulk_global,1,MPI_RPREC,MPI_SUM, &
        MPI_COMM_WORLD,ierr)
    ubulk = ubulk_global
    
    eps = 0.20970_rprec*ubulk ! amplitude of perturbation    
    
    ! change orientation of localized perturbation
    do jx = 1, nx
        do jy = 1, ny
            xp(jx,jy) = (x(jx)-xc)*cos(theta) - (y(jy)-yc)*sin(theta)
            yp(jx,jy) = (x(jx)-xc)*sin(theta) + (y(jy)-yc)*cos(theta) 
        end do
    end do

    ! main loop for perturbation calculations
    do jx = 1, nx
        do jy = 1, ny
            do jz = 1, nz
                f(jz) = (1._rprec+z(jz)-zc)**p*(1._rprec-z(jz)+zc)**q
                g(jx,jy) = xp(jx,jy)/lx*yp(jx,jy)* &
                   exp(-(xp(jx,jy)/lx)**2._rprec-(yp(jx,jy)/ly)**2._rprec)
                dfdz(jz) = p*(1._rprec+z(jz)-zc)**(p-1._rprec)*(1._rprec-z(jz)+zc)**q- &
                   q*(1._rprec+z(jz)-zc)**p*(1._rprec-z(jz)+zc)**(q-1._rprec)
                dgdy(jx,jy) = (xp(jx,jy)/lx)* &
                   (1._rprec-2._rprec*(yp(jx,jy)/ly)**2._rprec)* &
                   exp(-(xp(jx,jy)/lx)**2._rprec-(yp(jx,jy)/ly)**2._rprec)

                upi(jx,jy,jz) = -eps*g(jx,jy)*dfdz(jz)*sin(theta)
                vpi(jx,jy,jz) = -eps*g(jx,jy)*dfdz(jz)*cos(theta)
               ! wpi_uv(jx,jy,jz) = eps*f(jz)*dgdy(jx,jy)
            end do
        end do
    end do 
    
    do jx = 1, nx
        do jy = 1, ny
            do jz = 1, nz
                f(jz) = (1._rprec+z_w(jz)-zc)**p*(1._rprec-z_w(jz)+zc)**q
                g(jx,jy) = xp(jx,jy)/lx*yp(jx,jy)* &
                   exp(-(xp(jx,jy)/lx)**2._rprec-(yp(jx,jy)/ly)**2._rprec)
                dfdz(jz) = p*(1._rprec+z_w(jz)-zc)**(p-1._rprec)* &
                   (1._rprec-z_w(jz)+zc)**q- &
                   q*(1._rprec+z_w(jz)-zc)**p*(1._rprec-z_w(jz)+zc)**(q-1._rprec)
                dgdy(jx,jy) = (xp(jx,jy)/lx)* &
                   (1._rprec-2._rprec*(yp(jx,jy)/ly)**2._rprec)* &
                   exp(-(xp(jx,jy)/lx)**2._rprec-(yp(jx,jy)/ly)**2._rprec)

                wpi(jx,jy,jz) = eps*f(jz)*dgdy(jx,jy)
            end do
        end do
    end do 

   ! wpi = interp_to_w_grid(wpi_uv(1:nx,1:ny,lbz:nz),lbz)
   ! write(*,*) 'proc:',coord,'wpi(1):',wpi(60,60,1),'wpi(nz):',wpi(60,60,nz)
    deallocate(wpi_uv)
endif


end subroutine ic_pert

!*******************************************************************************
subroutine ic_les()
!*******************************************************************************
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
real(rprec) :: angle

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

    ! For channel flow, choose closest wall
    if(lbc_mom  > 0 .and. ubc_mom > 0) z = min(z, dz*nproc*(nz-1) - z)
    ! For upside-down half-channel, choose upper wall
    if(lbc_mom == 0 .and. ubc_mom > 0) z = dz*nproc*(nz-1) - z

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

! Modify angle of bulk flow based on pressure gradient
if (mean_p_force_x == 0._rprec) then
    angle = 3.14159265358979323846_rprec/2._rprec
else
    angle = atan(mean_p_force_y/mean_p_force_x)
end if

! Rescale noise depending on distance from wall and mean log profile
! z is in meters
do jz = 1, nz
#ifdef PPMPI
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
#else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
#endif

    ! For channel flow, choose closest wall
    if(lbc_mom  > 0 .and. ubc_mom > 0) z = min(z, dz*nproc*(nz-1)*z_i - z)
    ! For upside-down half-channel, choose upper wall
    if(lbc_mom == 0 .and. ubc_mom > 0) z = dz*nproc*(nz-1)*z_i - z

    if (z <= z_i) then
        u(:,:,jz) = u(:,:,jz) * (1._rprec-z / z_i) + ubar(jz)*cos(angle)
        v(:,:,jz) = v(:,:,jz) * (1._rprec-z / z_i) + ubar(jz)*sin(angle)
        w(:,:,jz) = w(:,:,jz) * (1._rprec-z / z_i)
    else
        u(:,:,jz) = u(:,:,jz) * 0.01_rprec + ubar(jz)*cos(angle)
        v(:,:,jz) = v(:,:,jz) * 0.01_rprec + ubar(jz)*sin(angle)
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

! Set upper boundary condition as zero for u, v, and w
#ifdef PPMPI
if (coord == nproc-1) then
#endif
    w(1:nx, 1:ny, nz) = 0._rprec
    u(1:nx, 1:ny, nz) = 0._rprec
    v(1:nx, 1:ny, nz) = 0._rprec
#ifdef PPMPI
end if
#endif

end subroutine ic_les

end subroutine initial
