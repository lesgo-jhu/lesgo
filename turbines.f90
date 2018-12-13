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
use turbine_indicator
use stat_defs, only : wind_farm
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array
#endif

implicit none

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize

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
! Dynamically change Ct_prime from input_turbines/Ct_prime.dat
logical, public :: dyn_Ct_prime
! disk-avg time scale in seconds (default 600)
real(rprec), public :: T_avg_dim
! filter size as multiple of grid spacing
real(rprec), public :: alpha1
real(rprec), public :: alpha2
! indicator function only includes values above this threshold
real(rprec), public :: filter_cutoff
! Correct ADM for filtered indicator function
logical, public :: adm_correction
! Number of timesteps between the output
integer, public :: tbase

! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: sx           ! spacing in the x-direction, multiple of diameter
real(rprec) :: sy           ! spacing in the y-direction

! Arrays for interpolating dynamic controls
real(rprec), dimension(:,:), allocatable :: theta1_arr
real(rprec), dimension(:), allocatable :: theta1_time
real(rprec), dimension(:,:), allocatable :: theta2_arr
real(rprec), dimension(:), allocatable :: theta2_time
real(rprec), dimension(:,:), allocatable :: Ct_prime_arr
real(rprec), dimension(:), allocatable :: Ct_prime_time

! Input files
character(*), parameter :: input_folder = 'input_turbines/'
character(*), parameter :: param_dat = path // input_folder // 'param.dat'
character(*), parameter :: theta1_dat = path // input_folder // 'theta1.dat'
character(*), parameter :: theta2_dat = path // input_folder // 'theta2.dat'
character(*), parameter :: Ct_prime_dat = path // input_folder // 'Ct_prime.dat'

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
implicit none

real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.turbines_init'
integer :: fid
real(rprec) :: T_avg_dim_file, delta1, delta2
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
allocate(forcing_fid(nloc))

! Create turbine directory
call system("mkdir -vp turbine")

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

!Compute a lookup table object for the indicator function
delta1 = alpha1 * sqrt(dx**2 + dy**2 + dz**2)
delta2 = alpha2 * sqrt(dx**2 + dy**2 + dz**2)
do s = 1, nloc
    call  wind_farm%turbine(s)%turb_ind_func%init(delta1, delta2,              &
            wind_farm%turbine(s)%thk, wind_farm%turbine(s)%dia)
end do

! Find turbine nodes - including filtered ind, n_hat, num_nodes, and nodes for
! each turbine. Each processor finds turbines in its domain
call turbines_nodes

! Read the time-averaged disk velocities from file if available
if (coord == 0) then
    inquire (file=u_d_T_dat, exist=exst)
    if (exst) then
        write(*,*) 'Reading from file ', trim(u_d_T_dat)
        open(newunit=fid, file=u_d_T_dat, status='unknown', form='formatted',  &
            position='rewind')
        do i=1,nloc
            read(fid,*) wind_farm%turbine(i)%u_d_T
        end do
        read(fid,*) T_avg_dim_file
        if (T_avg_dim_file /= T_avg_dim) then
            write(*,*) 'Time-averaging window does not match value in ',       &
                       trim(u_d_T_dat)
        end if
        close (fid)
    else
        write (*, *) 'File ', trim(u_d_T_dat), ' not found'
        write (*, *) 'Assuming u_d_T = -1. for all turbines'
        do k=1,nloc
            wind_farm%turbine(k)%u_d_T = -1.
        end do
    end if
end if

! Generate top of domain file
if (coord .eq. nproc-1) then
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        position='rewind')
    close(fid)
end if

! Generate the files for the turbine forcing output
if(coord==0) then
    do s=1,nloc
        call string_splice( string1, path // 'turbine/turbine_', s, '.dat' )
        open(newunit=forcing_fid(s), file=string1, status='unknown',           &
            form='formatted', position='append')
    end do
end if

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
real(rprec), dimension(nz_tot) :: z_tot

#ifdef PPMPI
real(rprec), dimension(nloc) :: buffer_array
#endif
real(rprec), pointer, dimension(:) :: x, y, z

real(rprec) :: filt
real(rprec), dimension(nloc) :: sumA, turbine_vol

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

sumA = 0._rprec

! z_tot for total domain (since z is local to the processor)
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
                    sumA(s) = sumA(s) + filt * dx * dy * dz
                end if
           end do
       end do
    end do
    wind_farm%turbine(s)%num_nodes = count_n

    ! Calculate turbine volume
    turbine_vol(s) = 0.25 * pi* p_dia**2 * p_thk

end do

! Sum the indicator function across all processors if using MPI
#ifdef PPMPI
buffer_array = sumA
call MPI_Allreduce(buffer_array, sumA, nloc, MPI_rprec, MPI_SUM, comm, ierr)
#endif

! Normalize the indicator function integrate to 1
! All integration will simply use a sum, so include the volume of each cell
! in ind
do s = 1, nloc
    wind_farm%turbine(s)%ind=wind_farm%turbine(s)%ind(:)/sumA(s)
end do

! Cleanup
nullify(x,y,z)

end subroutine turbines_nodes

!*******************************************************************************
subroutine turbines_forcing()
!*******************************************************************************
!
! This subroutine applies the drag-disk forcing
!
use param, only : pi
use sim_param, only : u, v, w, fxa, fya, fza
use functions, only : linear_interp, interp_to_uv_grid
implicit none

character(*), parameter :: sub_name = mod_name // '.turbines_forcing'

real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec), pointer :: p_Ct_prime => null()
integer, pointer :: p_icp => null(), p_jcp => null(), p_kcp => null()

integer :: fid

real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vel
real(rprec), dimension(nloc) :: u_vel_center, v_vel_center, w_vel_center
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nloc) :: buffer_array

nullify(y,z)
y => grid % y
z => grid % z

allocate(w_uv(ld,ny,lbz:nz))

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
    if (dyn_Ct_prime) wind_farm%turbine(s)%Ct_prime =                          &
        linear_interp(Ct_prime_time, Ct_prime_arr(s,:), total_time_dim)
end do

! Recompute the turbine position if theta1 or theta2 can change
if (dyn_theta1 .or. dyn_theta2) call turbines_nodes

!Each processor calculates the weighted disk-averaged velocity
disk_avg_vel = 0._rprec
u_vel_center = 0._rprec
v_vel_center = 0._rprec
w_vel_center = 0._rprec
do s=1,nloc
    ! Calculate total disk-averaged velocity for each turbine
    ! (current, instantaneous) in the normal direction by integrating the
    ! velocity times the indicator function
    do l=1,wind_farm%turbine(s)%num_nodes
        i2 = wind_farm%turbine(s)%nodes(l,1)
        j2 = wind_farm%turbine(s)%nodes(l,2)
        k2 = wind_farm%turbine(s)%nodes(l,3)
        disk_avg_vel(s) = disk_avg_vel(s) + dx*dy*dz*wind_farm%turbine(s)%ind(l)&
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
end do

! Calculate disk velocities by summing all processors and multiplying by disk volume
#ifdef PPMPI
call MPI_Allreduce(disk_avg_vel, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
disk_avg_vel = buffer_array
call MPI_Allreduce(u_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
u_vel_center = buffer_array
call MPI_Allreduce(v_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
v_vel_center = buffer_array
call MPI_Allreduce(w_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
w_vel_center = buffer_array
#endif

! Update epsilon for the new timestep (for cfl_dt)
if (T_avg_dim > 0.) then
    eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
else
    eps = 1.
end if

! Calculate and apply disk force
do s=1,nloc
    !set pointers
    p_u_d => wind_farm%turbine(s)%u_d
    p_u_d_T => wind_farm%turbine(s)%u_d_T
    p_f_n => wind_farm%turbine(s)%f_n
    p_Ct_prime => wind_farm%turbine(s)%Ct_prime

    !add this current value to the "running average" (first order filter)
    p_u_d = disk_avg_vel(s)
    if (adm_correction) then
        p_u_d = p_u_d /(1 + 0.25_rprec                                         &
            * (1-wind_farm%turbine(s)%turb_ind_func%M)*p_Ct_prime)
    end if
    p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

    !calculate total thrust force for each turbine  (per unit mass)
    !force is normal to the surface (calc from u_d_T, normal to surface)
    !write force to array that will be transferred via MPI
    p_f_n = -0.5*p_Ct_prime*abs(p_u_d_T)*p_u_d_T*0.25*pi*wind_farm%turbine(s)%dia**2

    !write values to file
    if (modulo (jt_total, tbase) == 0 .and. coord == 0) then
        write( forcing_fid(s), *) total_time_dim, u_vel_center(s),         &
            v_vel_center(s), w_vel_center(s), -p_u_d, -p_u_d_T,            &
            wind_farm%turbine(s)%theta1, wind_farm%turbine(s)%theta2,      &
            p_Ct_prime
    end if

    do l=1,wind_farm%turbine(s)%num_nodes
        i2 = wind_farm%turbine(s)%nodes(l,1)
        j2 = wind_farm%turbine(s)%nodes(l,2)
        k2 = wind_farm%turbine(s)%nodes(l,3)
        ind2 = wind_farm%turbine(s)%ind(l)
        fxa(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(1)*ind2
        fya(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(2)*ind2
        fza(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(3)*ind2
    end do
end do

!spatially average velocity at the top of the domain and write to file
if (coord .eq. nproc-1) then
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        action='write', position='append')
    write(fid,*) total_time, sum(u(:,:,nz-1))/(nx*ny)
    close(fid)
end if

! Cleanup
deallocate(w_uv)
nullify(y,z)
nullify(p_icp, p_jcp, p_kcp)

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
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_checkpoint'
integer :: fid

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
if (coord == 0) then
    open(newunit=fid, file=u_d_T_dat, status='unknown', form='formatted',      &
        position='rewind')
    do i=1,nloc
        write(fid,*) wind_farm%turbine(i)%u_d_T
    end do
    write(fid,*) T_avg_dim
    close (fid)
end if

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
    open(newunit=fid, file=param_dat, status='unknown', form='formatted',      &
        position='rewind')
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
use param, only: pi
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
    open(newunit=fid, file=theta1_dat, status='unknown', form='formatted',     &
        position='rewind')
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
    open(newunit=fid, file=theta2_dat, status='unknown', form='formatted',     &
        position='rewind')
    do i = 1, num_t
        read(fid,*) theta2_time(i), theta2_arr(:,i)
    end do
end if

! Read the Ct_prime input data
if (dyn_Ct_prime) then
    ! Count number of entries and allocate
    num_t = count_lines(Ct_prime_dat)
    allocate( Ct_prime_time(num_t) )
    allocate( Ct_prime_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=Ct_prime_dat, status='unknown', form='formatted',   &
        position='rewind')
    do i = 1, num_t
        read(fid,*) Ct_prime_time(i), Ct_prime_arr(:,i)
    end do
end if

end subroutine read_control_files

!*******************************************************************************
function count_lines(fname) result(N)
!*******************************************************************************
!
! This function counts the number of lines in a file
!
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
open(newunit=fid, file=trim(fname), status='unknown', form='formatted',        &
    position='rewind')

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

end module turbines
