!!
!!  Copyright (C) 2010-2016  Johns Hopkins University
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
module turbines
!*******************************************************************************
! This module contains all of the subroutines associated with drag-disk turbines

use types, only : rprec
use param
use grid_m
use messages
use string_util
use stat_defs, only : wind_farm
use bi_pchip
use wake_model_estimator
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array
#endif

implicit none

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize,&
          generate_splines

character (*), parameter :: mod_name = 'turbines'

! The following values are read from the input file
! number of turbines in the x-direction
integer, public :: num_x
! number of turbines in the y-direction
integer, public :: num_y
! baseline diameter in meters
real(rprec), public :: dia_all
! baseline height in meters
real(rprec), public :: height_all
! baseline thickness in meters
real(rprec), public :: thk_all
! orientation of turbines
integer, public :: orientation
! stagger percentage from baseline
real(rprec), public :: stag_perc
! angle from upstream (CCW from above, -x dir is zero)
real(rprec), public :: theta1_all
! angle above horizontal
real(rprec), public :: theta2_all
! thrust coefficient (default 1.33)
real(rprec), public :: Ct_prime
! Read parameters from input_turbines/param.dat
logical, public :: read_param
! Dynamically change theta1 from input_turbines/theta1.dat
logical, public :: dyn_theta1
! Dynamically change theta2 from input_turbines/theta2.dat
logical, public :: dyn_theta2
! disk-avg time scale in seconds (default 600)
real(rprec), public :: T_avg_dim
! filter size as multiple of grid spacing
real(rprec), public :: alpha
! indicator function only includes values above this threshold
real(rprec), public :: filter_cutoff
! Number of timesteps between the output
integer, public :: tbase
! Cp_prime corrections
real(rprec), public :: phi_a, phi_b, phi_c, phi_d, phi_x0
! Air density
real(rprec), public :: rho
! Inertia (kg*m^2)
real(rprec), public :: inertia_all
! Torque gain (kg*m^2)
real(rprec), public :: torque_gain
! time constant for estimating freestream velocity [seconds]
real(rprec), public :: tau_U_infty = 300
! std. deviation of noise of velocity deficit
real(rprec), public :: sigma_du = 0.5
! std. deviation of noise of wake expansion coefficient
real(rprec), public :: sigma_k = 0.001
! std. deviation of noise of power measurements
real(rprec), public :: sigma_uhat = 1.0
! std. deviation of noise of rotational speed
real(rprec), public :: sigma_omega = 0.01
! Number of members in ensemble
integer, public :: num_ensemble = 50

! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: sx           ! spacing in the x-direction, multiple of diameter
real(rprec) :: sy           ! spacing in the y-direction

! Arrays for interpolating dynamic controls
real(rprec), dimension(:,:), allocatable :: theta1_arr
real(rprec), dimension(:), allocatable :: theta1_time
real(rprec), dimension(:,:), allocatable :: theta2_arr
real(rprec), dimension(:), allocatable :: theta2_time

! Arrays for interpolating power and thrust coefficients for LES
type(bi_pchip_t), public :: Cp_prime_spline, Ct_prime_spline

! Input files
character(*), parameter :: input_folder = 'input_turbines/'
character(*), parameter :: param_dat = path // input_folder // 'param.dat'
character(*), parameter :: theta1_dat = path // input_folder // 'theta1.dat'
character(*), parameter :: theta2_dat = path // input_folder // 'theta2.dat'
character(*), parameter :: Ct_dat = path // input_folder // 'Ct.dat'
character(*), parameter :: Cp_dat = path // input_folder // 'Cp.dat'
character(*), parameter :: lambda_dat = path // input_folder // 'lambda.dat'
character(*), parameter :: beta_dat = path // input_folder // 'beta.dat'

! Output files
character(*), parameter :: output_folder = 'turbine/'
character(*), parameter :: vel_top_dat = path // output_folder // 'vel_top.dat'
character(*), parameter :: u_d_T_dat = path // output_folder // 'u_d_T.dat'
integer, dimension(:), allocatable :: forcing_fid

! epsilon used for disk velocity time-averaging
real(rprec) :: eps

! Commonly used indices
integer :: i, j, k, i2, j2, k2, l, s
integer :: k_start, k_end

! Variables to keep track of which processors have turbines
integer, dimension(:), allocatable :: turbine_in_proc_array
logical :: turbine_in_proc = .false.
#ifdef PPMPI
integer :: turbine_in_proc_cnt = 0
logical :: buffer_logical
#endif

! Wake model
type(wake_model_estimator_t) :: wm
character(*), parameter :: wm_path = path // 'wake_model'
integer, dimension(:), allocatable :: wm_fid
type(bi_pchip_t), public :: wm_Cp_prime_spline, wm_Ct_prime_spline

contains

!*******************************************************************************
subroutine turbines_init()
!*******************************************************************************
!
! This subroutine creates the 'turbine' folder and starts the turbine forcing
! output files. It also creates the indicator function (Gaussian-filtered from
! binary locations - in or out) and sets values for turbine type
! (node locations, etc)
!
use open_file_fid_mod
implicit none

real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.turbines_init'
integer :: fid
real(rprec) :: T_avg_dim_file, delta2
logical :: test_logical, exst
character (100) :: string1

! Set pointers
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

! Allocate and initialize
nloc = num_x*num_y
nullify(wind_farm%turbine)
allocate(wind_farm%turbine(nloc))
allocate(turbine_in_proc_array(nproc-1))
allocate(forcing_fid(nloc))
turbine_in_proc_array = 0

! Create turbine directory
call system("mkdir -vp turbine")
call system("mkdir -vp wake_model")

! Non-dimensionalize length values by z_i
height_all = height_all / z_i
dia_all = dia_all / z_i
thk_all = thk_all / z_i

! Spacing between turbines (as multiple of mean diameter)
sx = L_x / (num_x * dia_all )
sy = L_y / (num_y * dia_all )

! Place the turbines and specify some parameters
call place_turbines

! Resize thickness to capture at least on plane of gridpoints
! and set baseline values for size
do k = 1, nloc
    wind_farm%turbine(k)%thk = max(wind_farm%turbine(k)%thk, dx * 1.01)
    wind_farm%turbine(k)%vol_c = dx*dy*dz/(pi/4.*(wind_farm%turbine(k)%dia)**2 &
        * wind_farm%turbine(k)%thk)
end do

! Specify starting and ending indices for the processor
#ifdef PPMPI
k_start = 1+coord*(nz-1)
k_end = nz-1+coord*(nz-1)
#else
k_start = 1
k_end = nz
#endif

! Find the center of each turbine
do k = 1,nloc
    wind_farm%turbine(k)%icp = nint(wind_farm%turbine(k)%xloc/dx)
    wind_farm%turbine(k)%jcp = nint(wind_farm%turbine(k)%yloc/dy)
    wind_farm%turbine(k)%kcp = nint(wind_farm%turbine(k)%height/dz + 0.5)

    ! Check if turbine is the current processor
    test_logical = wind_farm%turbine(k)%kcp >= k_start .and.                   &
           wind_farm%turbine(k)%kcp<=k_end
    if (test_logical) then
        wind_farm%turbine(k)%center_in_proc = .true.
    else
        wind_farm%turbine(k)%center_in_proc = .false.
    end if

    ! Make kcp the local index
    wind_farm%turbine(k)%kcp = wind_farm%turbine(k)%kcp - k_start + 1

end do

! Read dynamic control input files
call read_control_files

! Read power and thrust coefficient curves
call generate_splines

!Compute a lookup table object for the indicator function
delta2 = alpha**2 * (dx**2 + dy**2 + dz**2)
do s = 1, nloc
    call  wind_farm%turbine(s)%turb_ind_func%init(delta2,                      &
            wind_farm%turbine(s)%thk, wind_farm%turbine(s)%dia,                &
            max( max(nx, ny), nz) )
end do

! Find turbine nodes - including filtered ind, n_hat, num_nodes, and nodes for
! each turbine. Each processor finds turbines in its domain
call turbines_nodes

! Read the time-averaged disk velocities from file if available
if (coord == 0) then
    inquire (file=u_d_T_dat, exist=exst)
    if (exst) then
        write(*,*) 'Reading from file ', trim(u_d_T_dat)
        fid = open_file_fid( u_d_T_dat, 'rewind', 'formatted' )
        do i=1,nloc
            read(fid,*) wind_farm%turbine(i)%u_d_T, wind_farm%turbine(i)%omega
        end do
        read(fid,*) T_avg_dim_file
        if (T_avg_dim_file /= T_avg_dim) then
            write(*,*) 'Time-averaging window does not match value in ',   &
                       trim(u_d_T_dat)
        end if
        close (fid)
    else
        write (*, *) 'File ', trim(u_d_T_dat), ' not found'
        write (*, *) 'Assuming u_d_T = -8, omega = 1 for all turbines'
        do k=1,nloc
            wind_farm%turbine(k)%u_d_T = -1._rprec
            wind_farm%turbine(k)%omega = 1._rprec
        end do
    end if
end if

! Generate top of domain file
if (coord .eq. nproc-1) then
    fid = open_file_fid( vel_top_dat, 'rewind', 'formatted' )
    close(fid)
end if

! Generate the files for the turbine forcing output
if(coord==0) then
    do s=1,nloc
        call string_splice( string1, path // 'turbine/turbine_', s, '.dat' )
        forcing_fid(s) = open_file_fid( string1, 'append', 'formatted' )
    end do
end if

if (coord==0) call wake_model_init

nullify(x,y,z)

end subroutine turbines_init

!*******************************************************************************
subroutine turbines_nodes
!*******************************************************************************
!
! This subroutine locates nodes for each turbine and builds the arrays: ind,
! n_hat, num_nodes, and nodes
!
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_nodes'

real(rprec) :: rx,ry,rz,r,r_norm,r_disk

real(rprec), pointer :: p_xloc => null(), p_yloc => null(), p_height => null()
real(rprec), pointer :: p_dia => null(), p_thk => null()
real(rprec), pointer :: p_theta1 => null(), p_theta2 => null()
real(rprec), pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3 => null()
integer :: icp, jcp, kcp
integer :: imax, jmax, kmax
integer :: min_i, max_i, min_j, max_j, min_k, max_k
integer :: count_i, count_n
real(rprec), dimension(:), allocatable :: z_tot

#ifdef PPMPI
real(rprec), dimension(:), allocatable :: buffer_array
#endif
real(rprec), pointer, dimension(:) :: x, y, z

real(rprec) :: filt
real(rprec), dimension(:), allocatable :: sumA, turbine_vol

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

turbine_in_proc = .false.

allocate(sumA(nloc))
allocate(turbine_vol(nloc))
sumA = 0

! z_tot for total domain (since z is local to the processor)
allocate(z_tot(nz_tot))
do k = 1,nz_tot
    z_tot(k) = (k - 0.5_rprec) * dz
end do

do s=1,nloc

    count_n = 0    !used for counting nodes for each turbine
    count_i = 1    !index count - used for writing to array "nodes"

    !set pointers
    p_xloc => wind_farm%turbine(s)%xloc
    p_yloc => wind_farm%turbine(s)%yloc
    p_height => wind_farm%turbine(s)%height
    p_dia => wind_farm%turbine(s)%dia
    p_thk => wind_farm%turbine(s)%thk
    p_theta1 => wind_farm%turbine(s)%theta1
    p_theta2 => wind_farm%turbine(s)%theta2
    p_nhat1 => wind_farm%turbine(s)%nhat(1)
    p_nhat2 => wind_farm%turbine(s)%nhat(2)
    p_nhat3 => wind_farm%turbine(s)%nhat(3)

    !identify "search area"
    imax = int(p_dia/dx + 2)
    jmax = int(p_dia/dy + 2)
    kmax = int(p_dia/dz + 2)

    !determine unit normal vector for each turbine
    p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat3 = sin(pi*p_theta2/180.)

    !determine nearest (i,j,k) to turbine center
    icp = nint(p_xloc/dx)
    jcp = nint(p_yloc/dy)
    kcp = nint(p_height/dz + 0.5)

    !determine limits for checking i,j,k
    !due to spectral BCs, i and j may be < 1 or > nx,ny
    !the mod function accounts for this when these values are used
    min_i = icp-imax
    max_i = icp+imax
    min_j = jcp-jmax
    max_j = jcp+jmax
    min_k = max((kcp-kmax),1)
    max_k = min((kcp+kmax),nz_tot)
    wind_farm%turbine(s)%nodes_max(1) = min_i
    wind_farm%turbine(s)%nodes_max(2) = max_i
    wind_farm%turbine(s)%nodes_max(3) = min_j
    wind_farm%turbine(s)%nodes_max(4) = max_j
    wind_farm%turbine(s)%nodes_max(5) = min_k
    wind_farm%turbine(s)%nodes_max(6) = max_k

    ! check neighboring grid points
    ! update num_nodes, nodes, and ind for this turbine
    ! split domain between processors
    ! z(nz) and z(1) of neighboring coords match so each coord gets
    ! (local) 1 to nz-1
    wind_farm%turbine(s)%ind = 0._rprec
    wind_farm%turbine(s)%nodes = 0
    wind_farm%turbine(s)%num_nodes = 0
    count_n = 0
    count_i = 1

    do k=k_start,k_end  !global k
        do j=min_j,max_j
            do i=min_i,max_i
                ! vector from center point to this node is (rx,ry,rz)
                ! with length r
                if (i<1) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (x(i2)-L_x) - p_xloc
                elseif (i>nx) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (L_x+x(i2)) - p_xloc
                else
                    i2 = i
                    rx = x(i) - p_xloc
                end if
                if (j<1) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (y(j2)-L_y) - p_yloc
                elseif (j>ny) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (L_y+y(j2)) - p_yloc
                else
                    j2 = j
                    ry = y(j) - p_yloc
                end if
                rz = z_tot(k) - p_height
                r = sqrt(rx*rx + ry*ry + rz*rz)
                !length projected onto unit normal for this turbine
                r_norm = abs(rx*p_nhat1 + ry*p_nhat2 + rz*p_nhat3)
                !(remaining) length projected onto turbine disk
                r_disk = sqrt(r*r - r_norm*r_norm)
                ! get the filter value
                filt = wind_farm%turbine(s)%turb_ind_func%val(r_disk, r_norm)
                if ( filt > filter_cutoff ) then
                    wind_farm%turbine(s)%ind(count_i) = filt
                    wind_farm%turbine(s)%nodes(count_i,1) = i2
                    wind_farm%turbine(s)%nodes(count_i,2) = j2
                    wind_farm%turbine(s)%nodes(count_i,3) = k-coord*(nz-1)!local
                    count_n = count_n + 1
                    count_i = count_i + 1
                    turbine_in_proc = .true.
                    sumA(s) = sumA(s) + filt * dx * dy * dz
                end if
           end do
       end do
    end do
    wind_farm%turbine(s)%num_nodes = count_n

    ! Calculate turbine volume
    turbine_vol(s) = pi/4. * p_dia**2 * p_thk

end do

! Sum the indicator function across all processors if using MPI
#ifdef PPMPI
allocate(buffer_array(nloc))
buffer_array = sumA
call MPI_Allreduce(buffer_array, sumA, nloc, MPI_rprec, MPI_SUM, comm, ierr)
deallocate(buffer_array)
#endif

! Normalize the indicator function
do s = 1, nloc
    wind_farm%turbine(s)%ind=wind_farm%turbine(s)%ind(:)*turbine_vol(s)/sumA(s)
end do

!each processor sends its value of turbine_in_proc
!if false, disk-avg velocity will not be sent (since it will always be 0.)
#ifdef PPMPI
turbine_in_proc_cnt = 0
if (coord == 0) then
    do i=1,nproc-1
        call MPI_recv(buffer_logical, 1, MPI_logical, i, 2, comm, status, ierr )
        if (buffer_logical) then
            turbine_in_proc_cnt = turbine_in_proc_cnt + 1
            turbine_in_proc_array(turbine_in_proc_cnt) = i
        end if
    end do
else
    call MPI_send(turbine_in_proc, 1, MPI_logical, 0, 2, comm, ierr )
end if
#endif

! Cleanup
deallocate(sumA)
deallocate(turbine_vol)
nullify(x,y,z)
deallocate(z_tot)

end subroutine turbines_nodes

!*******************************************************************************
subroutine turbines_forcing()
!*******************************************************************************
!
! This subroutine applies the drag-disk forcing
!
use sim_param, only : u,v,w, fxa,fya,fza
use functions, only : linear_interp, interp_to_uv_grid, bilinear_interp
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_forcing'

real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec), pointer :: p_Ct_prime => null(), p_Cp_prime => null()
real(rprec), pointer :: p_omega => null()
integer, pointer :: p_icp => null(), p_jcp => null(), p_kcp => null()

real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vel, disk_force
real(rprec), dimension(nloc) :: u_vel_center, v_vel_center, w_vel_center
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), pointer, dimension(:) :: y,z
real(rprec) :: const
real(rprec), dimension(:), allocatable :: beta

real(rprec), dimension(4*nloc) :: send_array
#ifdef PPMPI
real(rprec), dimension(4*nloc) :: recv_array
#endif

nullify(y,z)
y => grid % y
z => grid % z

allocate(w_uv(ld,ny,lbz:nz))

write(*,*) "ENTER TURBINES_FORCING", coord

#ifdef PPMPI
!syncing intermediate w-velocities
call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)
#endif

w_uv = interp_to_uv_grid(w, lbz)

! Do interpolation for dynamically changing parameters
do s = 1, nloc
    if (dyn_theta1) wind_farm%turbine(s)%theta1 =                              &
        linear_interp(theta1_time, theta1_arr(s,:), total_time_dim)
    if (dyn_theta2) wind_farm%turbine(s)%theta2 =                              &
        linear_interp(theta2_time, theta2_arr(s,:), total_time_dim)
end do

! Recompute the turbine position if theta1 or theta2 can change
if (dyn_theta1 .or. dyn_theta2) call turbines_nodes

write(*,*) "turbines_forcing loc 1:", coord

! Each processor calculates the weighted disk-averaged velocity
send_array = 0._rprec
disk_avg_vel = 0._rprec
u_vel_center = 0._rprec
v_vel_center = 0._rprec
w_vel_center = 0._rprec
if (turbine_in_proc) then
    do s=1,nloc
        ! Calculate total disk-averaged velocity for each turbine
        ! (current, instantaneous) in the normal direction. The weighted average
        ! is calculated using "ind"
        do l=1,wind_farm%turbine(s)%num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            disk_avg_vel(s) = disk_avg_vel(s) + wind_farm%turbine(s)%ind(l)    &
                            * ( wind_farm%turbine(s)%nhat(1)*u(i2,j2,k2)       &
                              + wind_farm%turbine(s)%nhat(2)*v(i2,j2,k2)       &
                              + wind_farm%turbine(s)%nhat(3)*w_uv(i2,j2,k2) )
        end do

        ! Set pointers
        p_icp => wind_farm%turbine(s)%icp
        p_jcp => wind_farm%turbine(s)%jcp
        p_kcp => wind_farm%turbine(s)%kcp

        ! Calculate disk center velocity
        if (wind_farm%turbine(s)%center_in_proc) then
            u_vel_center(s) = u(p_icp, p_jcp, p_kcp)
            v_vel_center(s) = v(p_icp, p_jcp, p_kcp)
            w_vel_center(s) = w_uv(p_icp, p_jcp, p_kcp)
        end if

        ! write this value to the array (which will be sent to coord 0)
        send_array(s)        = disk_avg_vel(s)
        send_array(nloc+s)   = u_vel_center(s)
        send_array(2*nloc+s) = v_vel_center(s)
        send_array(3*nloc+s) = w_vel_center(s)
    end do
end if

write(*,*) "turbines_forcing loc 2:", coord

! send the disk-avg values to coord==0
#ifdef PPMPI
call mpi_barrier (comm,ierr)

if (coord == 0) then
    ! Add all received values from other processors
    do i=1,turbine_in_proc_cnt
        j = turbine_in_proc_array(i)
        recv_array = 0._rprec
        call MPI_recv( recv_array, 4*nloc, MPI_rprec, j, 3, comm, status, ierr )
        send_array = send_array + recv_array
    end do

    ! Place summed values into arrays
    do s = 1, nloc
        disk_avg_vel(s) = send_array(s)
        u_vel_center(s) = send_array(nloc+s)
        v_vel_center(s) = send_array(2*nloc+s)
        w_vel_center(s) = send_array(3*nloc+s)
    end do
elseif (turbine_in_proc) then
    call MPI_send( send_array, 4*nloc, MPI_rprec, 0, 3, comm, ierr )
end if
#endif

write(*,*) "turbines_forcing loc 3:", coord

!Coord==0 takes that info and calculates total disk force, then sends it back
if (coord == 0) then
    !update epsilon for the new timestep (for cfl_dt)
    if (T_avg_dim > 0.) then
        eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
    else
        eps = 1.
    end if

    do s=1,nloc
        !set pointers
        p_u_d => wind_farm%turbine(s)%u_d
        p_u_d_T => wind_farm%turbine(s)%u_d_T
        p_f_n => wind_farm%turbine(s)%f_n
        p_Ct_prime => wind_farm%turbine(s)%Ct_prime
        p_Cp_prime => wind_farm%turbine(s)%Cp_prime
        p_omega => wind_farm%turbine(s)%omega

        ! Calculate rotational speed. Power needs to be dimensional.
        ! Use the previous step's values.
        const = -p_Cp_prime*0.5*rho*pi*0.25*(wind_farm%turbine(s)%dia*z_i)**2
        p_omega = p_omega + dt_dim / inertia_all *                             &
                (const*(p_u_d_T*u_star)**3/p_omega - 0.5*torque_gain*p_omega**2)

        !volume correction:
        !since sum of ind is turbine volume/(dx*dy*dz) (not exactly 1.)
        p_u_d = disk_avg_vel(s) * wind_farm%turbine(s)%vol_c

        !add this current value to the "running average" (first order filter)
        p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

        ! Calculate Ct_prime and Cp_prime
        call Ct_prime_spline%interp(0._rprec, -p_omega * 0.5                   &
            * wind_farm%turbine(s)%dia * z_i / p_u_d_T / u_star, p_Ct_prime)
        call Cp_prime_spline%interp(0._rprec, -p_omega * 0.5                   &
            * wind_farm%turbine(s)%dia * z_i / p_u_d_T / u_star, p_Cp_prime)

        !calculate total thrust force for each turbine  (per unit mass)
        !force is normal to the surface (calc from u_d_T, normal to surface)
        !write force to array that will be transferred via MPI
        p_f_n = -0.5*p_Ct_prime*abs(p_u_d_T)*p_u_d_T/wind_farm%turbine(s)%thk
        disk_force(s) = p_f_n

        !write current step's values to file
        if (modulo (jt_total, tbase) == 0) then
            write( forcing_fid(s), *) total_time_dim, u_vel_center(s),         &
                v_vel_center(s), w_vel_center(s), -p_u_d, -p_u_d_T,            &
                wind_farm%turbine(s)%theta1, wind_farm%turbine(s)%theta2,      &
                p_Ct_prime, p_Cp_prime, p_omega
        end if

    end do
end if

write(*,*) "turbines_forcing loc 4:", coord

!send total disk force to the necessary procs (with turbine_in_proc==.true.)
#ifdef PPMPI
if (coord == 0) then
    do i=1,turbine_in_proc_cnt
        j = turbine_in_proc_array(i)
        call MPI_send( disk_force, nloc, MPI_rprec, j, 5, comm, ierr )
    end do
elseif (turbine_in_proc) then
    call MPI_recv( disk_force, nloc, MPI_rprec, 0, 5, comm, status, ierr )
end if
#endif

write(*,*) "turbines_forcing loc 5:", coord

!apply forcing to each node
if (turbine_in_proc) then
    do s=1,nloc
        do l=1,wind_farm%turbine(s)%num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            ind2 = wind_farm%turbine(s)%ind(l)
            fxa(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(1)*ind2
            fya(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(2)*ind2
            fza(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(3)*ind2
        end do
    end do
end if

write(*,*) "turbines_forcing loc 6:", coord

!spatially average velocity at the top of the domain and write to file
if (coord .eq. nproc-1) then
    open(unit=1,file=vel_top_dat,status='unknown',form='formatted',            &
    action='write',position='append')
    write(1,*) total_time, sum(u(:,:,nz-1))/(nx*ny)
    close(1)
end if

write(*,*) "turbines_forcing loc 7:", coord

! Update wake model
if (coord == 0) then
    allocate( beta(nloc) )
    beta = 0._rprec
    call wm%advance(-wind_farm%turbine%u_d_T*u_star,                           &
        wind_farm%turbine(:)%omega, beta,                                      &
        torque_gain*wind_farm%turbine(:)%omega**2, dt_dim)

    ! write values to file
    if (modulo (jt_total, tbase) == 0) then
        do s = 1, nloc
            write( wm_fid(s), *) total_time_dim, wm%wm%Ctp(s), wm%wm%Cpp(s),   &
                wm%wm%uhat(s), wm%wm%omega(s), wm%wm%Phat(s)
        end do
    end if
    deallocate(beta)
end if

! Cleanup
deallocate(w_uv)
nullify(y,z)
nullify(p_icp, p_jcp, p_kcp)

write(*,*) "EXIT TURBINES_FORCING", coord

end subroutine turbines_forcing

!*******************************************************************************
subroutine turbines_finalize ()
!*******************************************************************************
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_finalize'

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
call turbines_checkpoint

!deallocate
deallocate(wind_farm%turbine)

end subroutine turbines_finalize

!*******************************************************************************
subroutine turbines_checkpoint ()
!*******************************************************************************
!
!
!
use open_file_fid_mod
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_checkpoint'
integer :: fid

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
if (coord == 0) then
    fid = open_file_fid( u_d_T_dat, 'rewind', 'formatted' )
    do i=1,nloc
        write(fid,*) wind_farm%turbine(i)%u_d_T, wind_farm%turbine(i)%omega
    end do
    write(fid,*) T_avg_dim
    close (fid)
end if

call wm%write_to_file(wm_path)

end subroutine turbines_checkpoint

!*******************************************************************************
subroutine turbine_vel_init(zo_high)
!*******************************************************************************
!
! called from ic.f90 if initu, lbc_mom==1, S_FLAG are all false.
! this accounts for the turbines when creating the initial velocity profile.
!
use param, only: zo
implicit none
character (*), parameter :: sub_name = mod_name // '.turbine_vel_init'

real(rprec), intent(inout) :: zo_high
real(rprec) :: cft, nu_w, exp_KE, induction_factor, Ct_noprime

! Convert Ct' to Ct
! a = Ct'/(4+Ct'), Ct = 4a(1-a)
induction_factor = Ct_prime / (4._rprec + Ct_prime)
Ct_noprime = 4*(induction_factor) * (1 - induction_factor)

! friction coefficient, cft
cft = pi*Ct_noprime/(4.*sx*sy)

!wake viscosity
nu_w = 28.*sqrt(0.5*cft)

!turbine friction height, Calaf, Phys. Fluids 22, 2010
zo_high = height_all*(1.+0.5*dia_all/height_all)**(nu_w/(1.+nu_w))* &
  exp(-1.*(0.5*cft/(vonk**2) + (log(height_all/zo* &
  (1.-0.5*dia_all/height_all)**(nu_w/(1.+nu_w))) )**(-2) )**(-0.5) )

exp_KE =  0.5*(log(0.45/zo_high)/0.4)**2

if(.false.) then
    write(*,*) 'sx,sy,cft: ',sx,sy,cft
    write(*,*) 'nu_w: ',nu_w
    write(*,*) 'zo_high: ',zo_high
    write(*,*) 'approx expected KE: ', exp_KE
end if
end subroutine turbine_vel_init

!*******************************************************************************
subroutine place_turbines
!*******************************************************************************
!
! This subroutine places the turbines on the domain. It also sets the values for
! each individual turbine. After the subroutine is called, the following values
! are set for each turbine in wind_farm: xloc, yloc, height, dia, thk, theta1,
! theta2, and Ct_prime.
!
use param, only: pi, z_i
use open_file_fid_mod
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

real(rprec) :: sxx, syy, shift_base, const
real(rprec) :: dummy, dummy2
logical :: exst
integer :: fid

! Read parameters from file if needed
if (read_param) then
    ! Check if file exists and open
    inquire (file = param_dat, exist = exst)
    if (.not. exst) then
        call error (sub_name, 'file ' // param_dat // 'does not exist')
    end if

    ! Check that there are enough lines from which to read data
    nloc = count_lines(param_dat)
    if (nloc < num_x*num_y) then
        nloc = num_x*num_y
        call error(sub_name, param_dat // 'must have num_x*num_y lines')
    else if (nloc > num_x*num_y) then
        call warn(sub_name, param_dat // ' has more than num_x*num_y lines. '  &
                  // 'Only reading first num_x*num_y lines')
    end if

    ! Read from parameters file, which should be in this format:
    ! xloc [meters], yloc [meters], height [meters], dia [meters], thk [meters],
    ! theta1 [degrees], theta2 [degrees], Ct_prime [-]
    write(*,*) "Reading from", param_dat
    fid = open_file_fid(param_dat, 'rewind', 'formatted')
    do k = 1, nloc
        read(fid,*) wind_farm%turbine(k)%xloc, wind_farm%turbine(k)%yloc,      &
            wind_farm%turbine(k)%height, wind_farm%turbine(k)%dia,             &
            wind_farm%turbine(k)%thk, wind_farm%turbine(k)%theta1,             &
            wind_farm%turbine(k)%theta2, wind_farm%turbine(k)%Ct_prime
    end do
    close(fid)

    ! Make lengths dimensionless
    do k = 1, nloc
        wind_farm%turbine(k)%xloc = wind_farm%turbine(k)%xloc / z_i
        wind_farm%turbine(k)%yloc = wind_farm%turbine(k)%yloc / z_i
        wind_farm%turbine(k)%height = wind_farm%turbine(k)%height / z_i
        wind_farm%turbine(k)%dia = wind_farm%turbine(k)%dia / z_i
        wind_farm%turbine(k)%thk = wind_farm%turbine(k)%thk / z_i
    end do
else
    ! Set values for each turbine based on values in input file
    wind_farm%turbine(:)%height = height_all
    wind_farm%turbine(:)%dia = dia_all
    wind_farm%turbine(:)%thk = thk_all
    wind_farm%turbine(:)%theta1 = theta1_all
    wind_farm%turbine(:)%theta2 = theta2_all
    wind_farm%turbine(:)%Ct_prime = Ct_prime

    ! Set baseline locations (evenly spaced, not staggered aka aligned)
    k = 1
    sxx = sx * dia_all  ! x-spacing with units to match those of L_x
    syy = sy * dia_all  ! y-spacing
    do i = 1,num_x
        do j = 1,num_y
            wind_farm%turbine(k)%xloc = sxx*real(2*i-1)/2
            wind_farm%turbine(k)%yloc = syy*real(2*j-1)/2
            k = k + 1
        end do
    end do

    ! Place turbines based on orientation flag
    ! This will shift the placement relative to the baseline locations abive
    select case (orientation)
        ! Evenly-spaced, not staggered
        case (1)

        ! Evenly-spaced, horizontally staggered only
        ! Shift each row according to stag_perc
        case (2)
            do i = 2, num_x
                do k = 1+num_y*(i-1), num_y*i
                    shift_base = syy * stag_perc/100.
                    wind_farm%turbine(k)%yloc = mod( wind_farm%turbine(k)%yloc &
                                                    + (i-1)*shift_base , L_y )
                end do
            end do

        ! Evenly-spaced, only vertically staggered (by rows)
        case (3)
            ! Make even rows taller
            do i = 2, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.+stag_perc/100.)
                end do
            end do
            ! Make odd rows shorter
            do i = 1, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.-stag_perc/100.)
                end do
            end do

        ! Evenly-spaced, only vertically staggered, checkerboard pattern
        case (4)
            k = 1
            do i = 1, num_x
                do j = 1, num_y
                    ! this should alternate between 1, -1
                    const = 2.*mod(real(i+j),2.)-1.
                    wind_farm%turbine(k)%height = height_all                   &
                                                  *(1.+const*stag_perc/100.)
                    k = k + 1
                end do
            end do

        ! Aligned, but shifted forward for efficient use of simulation space
        ! during CPS runs
        case (5)
        ! Shift in spanwise direction: Note that stag_perc is now used
            k=1
            dummy=stag_perc                                                    &
                  *(wind_farm%turbine(2)%yloc - wind_farm%turbine(1)%yloc)
            do i = 1, num_x
                do j = 1, num_y
                    dummy2=dummy*(i-1)
                    wind_farm%turbine(k)%yloc=mod( wind_farm%turbine(k)%yloc   &
                                                  + dummy2,L_y)
                    k=k+1
                end do
            end do

        case default
            call error (sub_name, 'invalid orientation')

    end select
end if

end subroutine place_turbines

!*******************************************************************************
subroutine read_control_files
!*******************************************************************************
!
! This subroutine reads the input files for dynamic controls with theta1,
! theta2, and Ct_prime. This is calles from turbines_init.
!
use open_file_fid_mod
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

integer :: fid, i, num_t

! Read the theta1 input data
if (dyn_theta1) then
    ! Count number of entries and allocate
    num_t = count_lines(theta1_dat)
    allocate( theta1_time(num_t) )
    allocate( theta1_arr(nloc, num_t) )

    ! Read values from file
    fid = open_file_fid(theta1_dat, 'rewind', 'formatted')
    do i = 1, num_t
        read(fid,*) theta1_time(i), theta1_arr(:,i)
    end do
end if

! Read the theta2 input data
if (dyn_theta2) then
    ! Count number of entries and allocate
    num_t = count_lines(theta2_dat)
    allocate( theta2_time(num_t) )
    allocate( theta2_arr(nloc, num_t) )

    ! Read values from file
    fid = open_file_fid(theta2_dat, 'rewind', 'formatted')
    do i = 1, num_t
        read(fid,*) theta2_time(i), theta2_arr(:,i)
    end do
end if

end subroutine read_control_files

!*******************************************************************************
subroutine generate_splines
!*******************************************************************************
use open_file_fid_mod
use functions, only : linear_interp
use pchip
implicit none
integer :: N, fid, Nlp
real(rprec), dimension(:), allocatable :: lambda
real(rprec), dimension(:,:), allocatable :: Ct, Cp, a
real(rprec), dimension(:,:), allocatable :: iCtp, iCpp, ilp
real(rprec) :: dlp, phi
real(rprec), dimension(:,:), allocatable :: Cp_prime_arr
real(rprec), dimension(:,:), allocatable :: Ct_prime_arr
real(rprec), dimension(:), allocatable :: lambda_prime
real(rprec), dimension(:), allocatable :: beta
type(pchip_t) :: cspl

! Read lambda
N = count_lines(lambda_dat)
allocate( lambda(N) )
fid = open_file_fid(lambda_dat, 'rewind', 'formatted')
do i = 1, N
    read(fid,*) lambda(i)
end do

! Read beta
N = count_lines(beta_dat)
allocate( beta(N) )
fid = open_file_fid(beta_dat, 'rewind', 'formatted')
do i = 1, N
    read(fid,*) beta(i)
end do

! Read Ct
allocate( Ct(size(beta), size(lambda)) )
fid = open_file_fid(Ct_dat, 'rewind', 'formatted')
do i = 1, size(beta)
    read(fid,*) Ct(i,:)
end do

! Read Cp
allocate( Cp(size(beta), size(lambda)) )
fid = open_file_fid(Cp_dat, 'rewind', 'formatted')
do i = 1, size(beta)
    read(fid,*) Cp(i,:)
end do

! Ct_prime and Cp_prime are only really defined if 0<=Ct<=1
! Prevent negative power coefficients
do i = 1, size(beta)
    do j = 1, size(lambda)
        if (Ct(i,j) < 0._rprec) Ct(i,j) = 0._rprec
        if (Ct(i,j) > 1._rprec) Ct(i,j) = 1._rprec
        if (Cp(i,j) < 0._rprec) Cp(i,j) = 0._rprec
    end do
end do

! Calculate induction factor
allocate( a(size(beta), size(lambda)) )
a = 0.5*(1._rprec-sqrt(1._rprec-Ct))

! Calculate local Ct, Cp, and lambda
allocate( iCtp(size(beta), size(lambda)) )
allocate( iCpp(size(beta), size(lambda)) )
allocate( ilp(size(beta), size(lambda)) )
iCtp = Ct/((1._rprec-a)**2)
iCpp = Cp/((1._rprec-a)**3)
do i = 1, size(beta)
    do j = 1, size(lambda)
        ilp(i,j) = lambda(j)/(1._rprec - a(i,j))
    end do
end do

! Allocate arrays
Nlp = size(lambda)*3
allocate( lambda_prime(Nlp) )
allocate( Ct_prime_arr(size(beta), size(lambda_prime)) )
allocate( Cp_prime_arr(size(beta), size(lambda_prime)) )

! First save the uncorrected splines for use with the wake model
! Set the lambda_prime's onto which these curves will be interpolated
lambda_prime(1) = minval(lambda)
lambda_prime(Nlp) = maxval(2._rprec*lambda)
dlp = (lambda_prime(Nlp) - lambda_prime(1))
dlp = dlp / (Nlp - 1)
do i = 2, Nlp - 1
    lambda_prime(i) = lambda_prime(i-1) + dlp
end do

! Interpolate onto Ct_prime and Cp_prime arrays
do i = 1, size(beta)
    cspl = pchip_t(ilp(i,:), iCtp(i,:))
    call cspl%interp(lambda_prime, Ct_prime_arr(i,:))
    cspl = pchip_t(ilp(i,:), iCpp(i,:))
    call cspl%interp(lambda_prime, Cp_prime_arr(i,:))
end do

! Make sure the edges of Cp_prime are zero with zero gradient
Cp_prime_arr(1,:) = 0._rprec
Cp_prime_arr(2,:) = 0._rprec
Cp_prime_arr(:,1) = 0._rprec
Cp_prime_arr(:,2) = 0._rprec
Cp_prime_arr(size(beta),:) = 0._rprec
Cp_prime_arr(size(beta)-1,:) = 0._rprec
Cp_prime_arr(:,Nlp) = 0._rprec
Cp_prime_arr(:,Nlp-1) = 0._rprec

! For Ct_prime, low beta and lambda are zero. All edges have zero gradient
Ct_prime_arr(1,:) = 0._rprec
Ct_prime_arr(2,:) = 0._rprec
Ct_prime_arr(:,1) = 0._rprec
Ct_prime_arr(:,2) = 0._rprec
Ct_prime_arr(size(beta),:) = Ct_prime_arr(size(beta)-1,:)
Ct_prime_arr(:,Nlp) = Ct_prime_arr(:,Nlp-1)

! Now generate splines
wm_Ct_prime_spline = bi_pchip_t(beta, lambda_prime, Ct_prime_arr)
wm_Cp_prime_spline = bi_pchip_t(beta, lambda_prime, Cp_prime_arr)

! Now save the adjusted splines for LES
! Adjust the lambda_prime and Cp_prime to use the LES velocity
do i = 1, size(beta)
    do j = 1, size(lambda)
        if (iCtp(i,j) > phi_x0) then
            phi = phi_a*phi_x0**3 + phi_b*phi_x0**2 + phi_c*phi_x0 + phi_d
        else
            phi = phi_a*iCtp(i,j)**3 + phi_b*iCtp(i,j)**2                      &
                + phi_c*iCtp(i,j) + phi_d
        end if
        ilp(i,j) = ilp(i,j) * phi
        iCpp(i,j) = iCpp(i,j) * phi**3
    end do
end do

! Set the lambda_prime's onto which these curves will be interpolated
lambda_prime(1) = maxval(ilp(:,1))
lambda_prime(Nlp) = minval(ilp(:,size(lambda)))
dlp = (lambda_prime(Nlp) - lambda_prime(1))
dlp = dlp / (Nlp - 1)
do i = 2, Nlp - 1
    lambda_prime(i) = lambda_prime(i-1) + dlp
end do

! Interpolate onto Ct_prime and Cp_prime arrays
do i = 1, size(beta)
    cspl = pchip_t(ilp(i,:), iCtp(i,:))
    call cspl%interp(lambda_prime, Ct_prime_arr(i,:))
    cspl = pchip_t(ilp(i,:), iCpp(i,:))
    call cspl%interp(lambda_prime, Cp_prime_arr(i,:))
end do

! Now generate splines
Ct_prime_spline = bi_pchip_t(beta, lambda_prime, Ct_prime_arr)
Cp_prime_spline = bi_pchip_t(beta, lambda_prime, Cp_prime_arr)

! Cleanup
deallocate (lambda)
deallocate (Ct)
deallocate (Cp)
deallocate (a)
deallocate (iCtp)
deallocate (iCpp)
deallocate (ilp)
deallocate (Cp_prime_arr)
deallocate (Ct_prime_arr)
deallocate (lambda_prime)
deallocate (beta)

end subroutine generate_splines

!*******************************************************************************
function count_lines(fname) result(N)
!*******************************************************************************
!
! This function counts the number of lines in a file
!
use open_file_fid_mod
use messages
use param, only : CHAR_BUFF_LENGTH
implicit none
character(*), intent(in) :: fname
logical :: exst
integer :: fid, ios
integer :: N

character(*), parameter :: sub_name = mod_name // '.count_lines'

! Check if file exists and open
inquire (file = trim(fname), exist = exst)
if (.not. exst) then
    call error (sub_name, 'file ' // trim(fname) // 'does not exist')
end if
fid = open_file_fid(trim(fname), 'rewind', 'formatted')

! count number of lines and close
ios = 0
N = 0
do
    read(fid, *, IOstat = ios)
    if (ios /= 0) exit
    N = N + 1
end do

! Close file
close(fid)

end function count_lines

!*******************************************************************************
subroutine wake_model_init
!*******************************************************************************
use param, only : CHAR_BUFF_LENGTH
use open_file_fid_mod
implicit none
real(rprec) :: U_infty
real(rprec), dimension(:), allocatable :: wm_k, wm_sx, wm_sy
integer :: i, j, fid
logical :: exst
character (CHAR_BUFF_LENGTH) :: fstring

fstring = path // 'wake_model/wm_est.dat'
inquire (file=fstring, exist=exst)

write(*,*) "ENTER"

if (exst) then
    write(*,*) 'Reading wake model estimator data from wake_model/'
    wm = wake_model_estimator_t(wm_path, wm_Ct_prime_spline,                   &
        wm_Cp_prime_spline, torque_gain, sigma_du, sigma_k, sigma_omega,       &
        sigma_uhat, tau_U_infty)
else
    ! Set initial velocity
    U_infty = 8._rprec

    ! Specify spacing and wake expansion coefficients
    allocate( wm_k(nloc) )
    allocate( wm_sx(nloc) )
    allocate( wm_sy(nloc) )
    wm_k = 0.05_rprec
    do i = 1, nloc
        wm_sx(i) = wind_farm%turbine(i)%xloc * z_i
        wm_sy(i) = wind_farm%turbine(i)%yloc * z_i
    end do

    ! Create wake model
    wm = wake_model_estimator_t(num_ensemble, wm_sx, wm_sy, U_infty,           &
        0.25*dia_all*z_i, wm_k, dia_all*z_i, rho, inertia_all, 2*nx, 2*ny,     &
        wm_Ct_prime_spline, wm_Cp_prime_spline,  torque_gain, sigma_du,        &
        sigma_k, sigma_omega, sigma_uhat, tau_U_infty)
    call wm%generate_initial_ensemble()

    ! Create output files
    allocate( wm_fid(nloc) )
    do i = 1, nloc
        call string_splice( fstring, path // 'turbine/wm_turbine_', i, '.dat' )
        wm_fid(i) = open_file_fid( fstring, 'append', 'formatted' )
    end do

    ! Cleanup
    deallocate(wm_k)
    deallocate(wm_sx)
    deallocate(wm_sy)
end if

write(*,*) "EXIT"

end subroutine wake_model_init

end module turbines
