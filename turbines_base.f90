!!
!!  Copyright (C) 2012-2016  Johns Hopkins University
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

module turbines_base
use types, only : rprec
use stat_defs, only : wind_farm
use param, only : path, coord
#ifdef PPMPI
  use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array 
#endif

implicit none
character (*), parameter :: mod_name = 'turbines'

! The following values are read from the input file
integer :: num_x            ! number of turbines in the x-direction
integer :: num_y            ! number of turbines in the y-direction

real(rprec) :: dia_all      ! baseline diameter in meters
real(rprec) :: height_all   ! baseline height in meters
real(rprec) :: thk_all      ! baseline thickness in meters

integer :: orientation      ! orientation 1=aligned, 2=horiz stagger,
                            !  3=vert stagger by row, 4=vert stagger checkerboard,
                            !  5=2 pushed forward for cps
real(rprec) :: stag_perc    ! stagger percentage from baseline

real(rprec) :: theta1_all   ! angle from upstream (CCW from above, -x dir is zero)
real(rprec) :: theta2_all   ! angle above horizontal

real(rprec) :: Ct_prime     ! thrust coefficient (default 1.33)

logical :: read_param       ! Read parameters from input_turbines/param.dat

logical :: dyn_theta1       ! Dynamically change theta1 from input_turbines/theta1.dat
logical :: dyn_theta2       ! Dynamically change theta2 from input_turbines/theta2.dat
logical :: dyn_Ct_prime     ! Dynamically change Ct_prime from input_turbines/Ct_prime.dat

real(rprec) :: T_avg_dim    ! disk-avg time scale in seconds (default 600)

real(rprec) :: alpha        ! filter size as multiple of grid spacing
integer :: trunc            ! Gaussian filter truncated after this many gridpoints
real(rprec) :: filter_cutoff  ! indicator function only includes values above this threshold

logical :: turbine_cumulative_time ! Used to read in the disk averaged velocities of the turbines

integer :: tbase            ! Number of timesteps between the output

! Arrays for interpolating dynamic controls
real(rprec), dimension(:,:), allocatable :: theta1_arr
real(rprec), dimension(:), allocatable :: theta1_time
real(rprec), dimension(:,:), allocatable :: theta2_arr
real(rprec), dimension(:), allocatable :: theta2_time
real(rprec), dimension(:,:), allocatable :: Ct_prime_arr
real(rprec), dimension(:), allocatable :: Ct_prime_time

! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: sx           ! spacing in the x-direction, multiple of (mean) diameter
real(rprec) :: sy           ! spacing in the y-direction
real(rprec) :: dummy,dummy2 ! used to shift the turbine positions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets the values for wind_farm based on values
!   read from the input file.

subroutine turbines_base_init()
use param, only: L_x, L_y, dx, dy, dz, pi, z_i
use open_file_fid_mod
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.turbines_base_init'
character(*), parameter :: param_dat = path // 'input_turbines/param.dat'
character(*), parameter :: theta1_dat = path // 'input_turbines/theta1.dat'
character(*), parameter :: theta2_dat = path // 'input_turbines/theta2.dat'
character(*), parameter :: Ct_prime_dat = path // 'input_turbines/Ct_prime.dat'

integer :: i, j, k, num_t
real(rprec) :: sxx, syy, shift_base, const
logical :: exst
integer :: fid, ios

! set turbine parameters
! turbines are numbered as follows:
!   #1 = turbine nearest (x,y)=(0,0)
!   #2 = next turbine in the y-direction, etc. (go along rows)

! Allocate wind turbine array derived type
nloc = num_x*num_y
nullify(wind_farm%turbine)
allocate(wind_farm%turbine(nloc))

! Read parameters from file if needed
if (read_param) then
    ! Check if file exists and open
    inquire (file = param_dat, exist = exst)
    if (exst) then
        fid = open_file_fid(param_dat, 'rewind', 'formatted')
    else
        call error (sub_name, 'file ' // param_dat // 'does not exist')
    end if

    ! count number of lines and close
    ios = 0
    nloc = 0
    do 
        read(fid, *, IOstat = ios)
        if (ios /= 0) exit
        nloc = nloc + 1
    enddo
    close(fid)
    
    ! Check that there are enough lines from which to read data
    if (nloc < num_x*num_y) then
        nloc = num_x*num_y
        call error(sub_name, param_dat // 'must have num_x*num_y lines')
    else if (nloc > num_x*num_y) then
        call warn(sub_name, param_dat // ' has more than num_x*num_y lines. ' &
                  // 'Only reading first num_x*num_y lines')
    end if
    
    ! Read from parameters file, which should be in this format:
    !   xloc [meters], yloc [meters], height [meters], dia [meters], thk [meters], 
    !   theta1 [degrees], theta2 [degrees], Ct_prime [-]
    write(*,*) "Reading from", param_dat
    fid = open_file_fid(param_dat, 'rewind', 'formatted')
    do k = 1, nloc
        read(fid,*) wind_farm%turbine(k)%xloc, wind_farm%turbine(k)%yloc, wind_farm%turbine(k)%height, &
                    wind_farm%turbine(k)%dia, wind_farm%turbine(k)%thk, wind_farm%turbine(k)%theta1,   &
                    wind_farm%turbine(k)%theta2, wind_farm%turbine(k)%Ct_prime
    enddo
    close(fid)
    
    ! Make turbine locations dimensionless
    do k = 1, nloc
        wind_farm%turbine(k)%xloc = wind_farm%turbine(k)%xloc / z_i
        wind_farm%turbine(k)%yloc = wind_farm%turbine(k)%yloc / z_i
    end do
    
else
    ! This will not yet set the locations
    wind_farm%turbine(:)%height = height_all
    wind_farm%turbine(:)%dia = dia_all
    wind_farm%turbine(:)%thk = thk_all 
    wind_farm%turbine(:)%theta1 = theta1_all
    wind_farm%turbine(:)%theta2 = theta2_all
    wind_farm%turbine(:)%Ct_prime = Ct_prime
endif

! Non-dimensionalize length values by z_i
do k = 1, nloc
    height_all = height_all / z_i
    dia_all = dia_all / z_i
    thk_all = thk_all / z_i
    wind_farm%turbine(k)%height = wind_farm%turbine(k)%height / z_i
    wind_farm%turbine(k)%dia = wind_farm%turbine(k)%dia / z_i
    ! Resize thickness capture at least on plane of gridpoints
    wind_farm%turbine(k)%thk = max(wind_farm%turbine(k)%thk / z_i, dx * 1.01)
    ! Set baseline values for size 
    wind_farm%turbine(k)%vol_c = dx*dy*dz/(pi/4.*(wind_farm%turbine(k)%dia)**2 * wind_farm%turbine(k)%thk)
end do

! Place turbines based on orientation flag
if (.not. read_param) then
    ! Spacing between turbines (as multiple of mean diameter)
    sx = L_x / (num_x * dia_all )
    sy = L_y / (num_y * dia_all )

    ! Baseline locations (evenly spaced, not staggered aka aligned)
    !  x,y-locations
    k = 1
    sxx = sx * dia_all  ! x-spacing with units to match those of L_x
    syy = sy * dia_all  ! y-spacing
    do i = 1,num_x
        do j = 1,num_y
            wind_farm%turbine(k)%xloc = sxx*real(2*i-1)/2
            wind_farm%turbine(k)%yloc = syy*real(2*j-1)/2
            k = k + 1
        enddo
    enddo
    
    select case (orientation)
        ! Evenly-spaced, not staggered
        !  Use baseline as set above      
        case (1)
    
        ! Evenly-spaced, horizontally staggered only
        ! Shift each row according to stag_perc    
        case (2)
            do i = 2, num_x
                do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
                    shift_base = syy * stag_perc/100.
                    wind_farm%turbine(k)%yloc = mod( wind_farm%turbine(k)%yloc + (i-1)*shift_base , L_y )
                enddo
            enddo
 
        ! Evenly-spaced, only vertically staggered (by rows)
        case (3)
            ! Make even rows taller
            do i = 2, num_x, 2
                do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
                    wind_farm%turbine(k)%height = height_all*(1.+stag_perc/100.)
                enddo
            enddo
            ! Make odd rows shorter
            do i = 1, num_x, 2
                do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
                    wind_farm%turbine(k)%height = height_all*(1.-stag_perc/100.)
                enddo
            enddo
       
        ! Evenly-spaced, only vertically staggered, checkerboard pattern
        case (4)
            k = 1
            do i = 1, num_x 
                do j = 1, num_y
                    const = 2.*mod(real(i+j),2.)-1.  ! this should alternate between 1, -1
                    wind_farm%turbine(k)%height = height_all*(1.+const*stag_perc/100.)
                    k = k + 1
                enddo
            enddo

        ! Aligned, but shifted forward for efficient use of simulation space during CPS runs
        ! Usual placement is baseline as set above
        case (5)
        ! Shift in spanwise direction: Note that stag_perc is now used
            k=1
            dummy=stag_perc*(wind_farm%turbine(2)%yloc - wind_farm%turbine(1)%yloc)
            do i = 1, num_x
                do j = 1, num_y
                    dummy2=dummy*(i-1)         
                    wind_farm%turbine(k)%yloc=mod(wind_farm%turbine(k)%yloc +dummy2,L_y)
                    k=k+1
                enddo
            enddo      

        case default
            call error (sub_name, 'invalid orientation')
        
    end select
end if

! Read the theta1 input data
if (dyn_theta1) then
    ! Check if file exists and open
    inquire (file = theta1_dat, exist = exst)
    if (exst) then
        fid = open_file_fid(theta1_dat, 'rewind', 'formatted')
    else
        call error (sub_name, 'file ' // theta1_dat // 'does not exist')
    end if

    ! count number of lines and close
    ios = 0
    num_t = 0
    do 
        read(fid, *, IOstat = ios)
        if (ios /= 0) exit
        num_t = num_t + 1
    enddo
    close(fid)

    allocate( theta1_time(num_t) )
    allocate( theta1_arr(nloc, num_t) )

    fid = open_file_fid(theta1_dat, 'rewind', 'formatted')
    do i = 1, num_t
        read(fid,*) theta1_time(i), theta1_arr(:,i)
    end do
end if

! Read the theta2 input data
if (dyn_theta2) then
    ! Check if file exists and open
    inquire (file = theta2_dat, exist = exst)
    if (exst) then
        fid = open_file_fid(theta2_dat, 'rewind', 'formatted')
    else
        call error (sub_name, 'file ' // theta2_dat // 'does not exist')
    end if

    ! count number of lines and close
    ios = 0
    num_t = 0
    do 
        read(fid, *, IOstat = ios)
        if (ios /= 0) exit
        num_t = num_t + 1
    enddo
    close(fid)

    allocate( theta2_time(num_t) )
    allocate( theta2_arr(nloc, num_t) )

    fid = open_file_fid(theta2_dat, 'rewind', 'formatted')
    do i = 1, num_t
        read(fid,*) theta2_time(i), theta2_arr(:,i)
    end do
end if

! Read the Ct_prime input data
if (dyn_Ct_prime) then
    ! Check if file exists and open
    inquire (file = Ct_prime_dat, exist = exst)
    if (exst) then
        fid = open_file_fid(Ct_prime_dat, 'rewind', 'formatted')
    else
        call error (sub_name, 'file ' // Ct_prime_dat // 'does not exist')
    end if

    ! count number of lines and close
    ios = 0
    num_t = 0
    do 
        read(fid, *, IOstat = ios)
        if (ios /= 0) exit
        num_t = num_t + 1
    enddo
    close(fid)

    allocate( Ct_prime_time(num_t) )
    allocate( Ct_prime_arr(nloc, num_t) )

    fid = open_file_fid(Ct_prime_dat, 'rewind', 'formatted')
    do i = 1, num_t
        read(fid,*) Ct_prime_time(i), Ct_prime_arr(:,i)
    end do
end if

end subroutine turbines_base_init

end module turbines_base
