!!
!!  Copyright (C) 2012-2013  Johns Hopkins University
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
use stat_defs, only : wind_farm, upstream_sensors
use param, only : path
$if ($MPI)
  use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array 
$endif

implicit none
character (*), parameter :: mod_name = 'turbines'

! The following values are read from the input file
integer :: num_x            ! number of turbines in the x-direction
integer :: num_y            ! number of turbines in the y-direction

real(rprec) :: dia_all      ! baseline diameter in meters
real(rprec) :: height_all   ! baseline height in meters
real(rprec) :: thk_all      ! baseline thickness in meters

integer :: orientation      ! orientation 1=aligned, 2=horiz stagger,
                            !  3=vert stagger by row, 4=vert stagger checkerboard
real(rprec) :: stag_perc    ! stagger percentage from baseline

real(rprec) :: theta1_all   ! angle from upstream (CCW from above, -x dir is zero)
real(rprec) :: theta2_all   ! angle above horizontal

integer :: control          ! Control method (1 = Ct_prime for all turbines)
real(rprec), dimension(:,:), allocatable :: Pref_list
real(rprec), dimension(:,:), allocatable :: Pref_t_list
real(rprec), dimension(:,:), allocatable :: Ct_prime_list
real(rprec), dimension(:,:), allocatable :: Ct_prime_t_list

logical :: up_vel_calc      ! Calculate upstream velocity
real(rprec) :: up_vel_nd    ! number of diameters upstream to calculate

real(rprec) :: Ct_prime_all ! thrust coefficient (default 1.33)
real(rprec) :: Ct_noprime   ! thrust coefficient (default 0.75)

real(rprec) :: T_avg_dim    ! disk-avg time scale in seconds (default 600)

real(rprec) :: alpha        ! filter size as multiple of grid spacing
integer :: trunc            ! Gaussian filter truncated after this many gridpoints
real(rprec) :: filter_cutoff  ! indicator function only includes values above this threshold

logical :: turbine_cumulative_time ! Used to read in the disk averaged velocities of the turbines

integer :: tbase     ! Number of timesteps between the output

! Turbine spacings, multiple of (mean) diameter
! This is read in for orientation = 7, and calculated for orientation < 7
real(rprec) :: sx                   ! spacing in the x-direction
real(rprec) :: sy                   ! spacing in the y-direction

! Only used for orientaiton = 5
real(rprec) :: x1, y1               ! Location of first turbine
 
! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: dummy,dummy2 ! used to shift the turbine positions

! PI control for back row
real(rprec) :: Kp_pref      ! proportional term
real(rprec) :: Ki_pref      ! integral term
real(rprec) :: e_pref       ! error in Pref

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
character(*), parameter :: turbine_locations_dat = path // 'turbine_locations.dat'

integer :: i, j, k
real(rprec) :: sxx, syy, shift_base, const
logical :: exst
integer :: fid, ios
real(rprec) :: xl, yl, zl

character(100) :: Pref_dat
character(100) :: Ct_prime_dat
character(5) :: chari
integer :: num_list

! set turbine parameters
! turbines are numbered as follows:
!   #1 = turbine nearest (x,y)=(0,0)
!   #2 = next turbine in the y-direction, etc. (go along rows)

! Allocate wind turbine array derived type
if (orientation == 6) then
! Count number of turbines listed in turbine_locations.dat
    ! Check if file exists and open
    inquire (file = turbine_locations_dat, exist = exst)
    if (exst) then
        fid = open_file_fid(turbine_locations_dat, 'rewind', 'formatted')
    else
        call error (sub_name, 'file ' // turbine_locations_dat // 'does not exist')
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
else
    nloc = num_x*num_y      !number of turbines (locations)
endif 

$if ($VERBOSE)
write(*,*) "Number of turbines: ", nloc
$endif
nullify(wind_farm%turbine)
allocate(wind_farm%turbine(nloc))
if (up_vel_calc) then
    nullify(upstream_sensors%turbine)
    allocate(upstream_sensors%turbine(nloc))
endif

! Non-dimensionalize length values by z_i
dia_all = dia_all / z_i
height_all = height_all / z_i
thk_all = thk_all / z_i
! Resize thickness capture at least on plane of gridpoints
thk_all = max ( thk_all, dx*1.01 )

! Set baseline values for size
wind_farm%turbine(:)%height = height_all
wind_farm%turbine(:)%dia = dia_all
wind_farm%turbine(:)%thk = thk_all                      
wind_farm%turbine(:)%vol_c =  dx*dy*dz/(pi/4.*(dia_all)**2 * thk_all) 
if (up_vel_calc) then
    upstream_sensors%turbine(:)%dia = dia_all
    upstream_sensors%turbine(:)%thk = thk_all                      
    upstream_sensors%turbine(:)%vol_c =  dx*dy*dz/(pi/4.*(dia_all)**2 * thk_all)   
endif

! Spacing between turbines (as multiple of mean diameter)
if (orientation < 7) then
    sx = L_x / (num_x * dia_all )
    sy = L_y / (num_y * dia_all )
else
    sx = sx / ( dia_all * z_i )
    sy = sy / ( dia_all * z_i )
endif

! Location of first turbine
if (orientation == 7) then
    x1 = x1 / z_i
    y1 = y1 / z_i
endif

if (orientation /= 6) then
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
endif

! HERE PLACE TURBINES (x,y-positions) BASED ON 'ORIENTATION' FLAG
if (orientation == 1) then
! Evenly-spaced, not staggered
!  Use baseline as set above       
 
elseif (orientation == 2) then
! Evenly-spaced, horizontally staggered only
! Shift each row according to stag_perc
    do i = 2, num_x
        do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
            shift_base = syy * stag_perc/100.
            wind_farm%turbine(k)%yloc = mod( wind_farm%turbine(k)%yloc + (i-1)*shift_base , L_y )
        enddo
    enddo
 
elseif (orientation == 3) then 
! Evenly-spaced, only vertically staggered (by rows)
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
 
elseif (orientation == 4) then        
! Evenly-spaced, only vertically staggered, checkerboard pattern
    k = 1
    do i = 1, num_x 
        do j = 1, num_y
            const = 2.*mod(real(i+j),2.)-1.  ! this should alternate between 1, -1
            wind_farm%turbine(k)%height = height_all*(1.+const*stag_perc/100.)
            k = k + 1
        enddo
    enddo

elseif (orientation == 5) then        
! Aligned, but shifted forward for efficient use of simulation space during CPS runs
! Usual placement is baseline as set above

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
elseif (orientation == 6) then
! Read x, y, and z locations in meters directly from turbine_locations.dat

    write(*,*) "Reading turbine locations from file."
    inquire (file = turbine_locations_dat, exist = exst)
    if (exst) then
        fid = open_file_fid(turbine_locations_dat, 'rewind', 'formatted')
    else
        call error (sub_name, 'file ' // turbine_locations_dat // 'does not exist')
    end if

    do k = 1, nloc
        read(fid,*) xl, yl, zl
        wind_farm%turbine(k)%xloc = xl / z_i
        wind_farm%turbine(k)%yloc = yl / z_i
        wind_farm%turbine(k)%height = zl / z_i
    enddo
    close(fid)
elseif (orientation == 7) then
! Aligned, but not evenly spaced
    k = 1
    do i = 1, num_x
        do j = 1, num_y
            wind_farm%turbine(k)%xloc = x1 + (i-1)*sxx
            wind_farm%turbine(k)%yloc = y1 + (j-1)*syy
            k = k + 1
        enddo
    enddo
endif

if (up_vel_calc) then
    do k = 1,nloc
        upstream_sensors%turbine(k)%xloc = wind_farm%turbine(k)%xloc - up_vel_nd*dia_all
        upstream_sensors%turbine(k)%yloc = wind_farm%turbine(k)%yloc
        upstream_sensors%turbine(:)%height = wind_farm%turbine(k)%height
    enddo
endif

$if ($VERBOSE)
do k = 1, nloc
    write(*,*) "Turbine ", k, " located at: ", wind_farm%turbine(k)%xloc, wind_farm%turbine(k)%yloc, wind_farm%turbine(k)%height
    if (up_vel_calc) then
        write(*,*) "Upstream sensor ", k, " located at: ", upstream_sensors%turbine(k)%xloc, upstream_sensors%turbine(k)%yloc, upstream_sensors%turbine(k)%height
    endif
enddo
$endif

if (control == 2) then
    ! Count number of lines
    Pref_dat = path // 'input/row_1_Pref.dat'
    num_list = get_number_of_lines(Pref_dat)

    ! Allocate
    allocate(Pref_list(num_x,num_list))
    allocate(Pref_t_list(num_x,num_list))

    do k = 1, num_x
        ! Open file
        write(chari, '(I5)') k
        Pref_dat = path // 'input/row_' // trim(adjustl(chari)) // '_Pref.dat'

        ! Count number of lines
        write(chari, '(I5)') num_list
        if (num_list /= get_number_of_lines(Pref_dat)) then
            call error(sub_name, 'file ' // Pref_dat // 'must have ' // trim(adjustl(chari)) // ' entries')
        endif

        ! Read from file
        call read_values_from_file(Pref_dat, Pref_t_list(k,:), Pref_list(k,:))

        $if ($VERBOSE)
        write(*,*) "Pref_t_list for row ", k, " is: ", Pref_t_list(k,:)
        write(*,*) "Pref_list for row ", k, " is: ", Pref_list(k,:)
        $endif
    enddo
elseif (control == 3 .OR. control == 6 .OR. control == 7) then
    ! Count number of lines
    if (control == 3) then
        Pref_dat = path // 'input/Pref.dat'
    else
        Pref_dat = path // 'input/Pref_tcm.dat'
    endif 
    num_list = get_number_of_lines(Pref_dat)

    ! Allocate
    allocate(Pref_list(1,num_list))
    allocate(Pref_t_list(1,num_list))

    ! Read from file
    call read_values_from_file(Pref_dat, Pref_t_list(1,:), Pref_list(1,:))

    $if ($VERBOSE)
    write(*,*) "Pref_t_list  is: ", Pref_t_list(1,:)
    write(*,*) "Pref_list is: ", Pref_list(1,:)
    $endif
elseif (control == 4) then
    ! Count number of lines
    Ct_prime_dat = path // 'input/row_1_Ct_prime.dat'
    num_list = get_number_of_lines(Ct_prime_dat)

    ! Allocate
    allocate(Ct_prime_list(num_x,num_list))
    allocate(Ct_prime_t_list(num_x,num_list))

    do k = 1, num_x
        ! Open file
        write(chari, '(I5)') k
        Ct_prime_dat = path // 'input/row_' // trim(adjustl(chari)) // '_Ct_prime.dat'

        ! Count number of lines
        write(chari, '(I5)') num_list
        if (num_list /= get_number_of_lines(Ct_prime_dat)) then
            call error(sub_name, 'file ' // Ct_prime_dat // 'must have ' // trim(adjustl(chari)) // ' entries')
        endif

        ! Read from file
        call read_values_from_file(Ct_prime_dat, Ct_prime_t_list(k,:), Ct_prime_list(k,:))

        $if ($VERBOSE)
        write(*,*) "Ct_prime_t_list for row ", k, " is: ", Ct_prime_t_list(k,:)
        write(*,*) "Ct_prime_list for row ", k, " is: ", Ct_prime_list(k,:)
        $endif
    enddo
elseif (control == 5) then
    ! Count number of lines
    Ct_prime_dat = path // 'input/Ct_prime.dat'
    num_list = get_number_of_lines(Ct_prime_dat)

    ! Allocate
    allocate(Ct_prime_list(1,num_list))
    allocate(Ct_prime_t_list(1,num_list))

    ! Read from file
    call read_values_from_file(Ct_prime_dat, Ct_prime_t_list(1,:), Ct_prime_list(1,:))

    $if ($VERBOSE)
    write(*,*) "Ct_prime_t_list  is: ", Ct_prime_t_list(1,:)
    write(*,*) "Ct_prime_list is: ", Ct_prime_list(1,:)
    $endif
endif


            
! orientation (angles)
wind_farm%turbine(:)%theta1 = theta1_all
wind_farm%turbine(:)%theta2 = theta2_all

if (up_vel_calc) then
    wind_farm%turbine(:)%theta1 = theta1_all
    wind_farm%turbine(:)%theta2 = theta2_all
endif

! local thrust coefficient
wind_farm%turbine(:)%Ct_prime = Ct_prime_all

end subroutine turbines_base_init

!**********************************************************************
function interpolate(x, y, xi) result(yi)
!**********************************************************************
!  This function interpolates from given x and y values (length n) to yi 
!  values corresponding to provided xi (length ni)
use messages
implicit none

real(rprec), intent(in), dimension(:) :: x, y
real(rprec), intent(in) :: xi
real(rprec) :: yi, dx, t
integer :: i, n

n = size(x);
if (size(y) /= n) then
    call error('interpolate','x and y are not the same size')
end if

if ( xi <= x(1) ) then
    yi = y(1)
else if ( xi >= x(n) ) then
    yi = y(n)
else
    i = 1
    do while ( xi > x(i + 1) )
        i = i + 1
    end do
    dx = x(i + 1) - x(i)
    t = ( xi - x(i) ) / dx
    yi = (1 - t) * y(i) + t * y(i + 1)
end if

end function interpolate

!**********************************************************************
function interpolate_vec(x, y, xi) result(yi)
!**********************************************************************
!  This function interpolates from given x and y values (length n) to yi 
!  values corresponding to provided xi (length ni)
use messages
implicit none

real(rprec), intent(in), dimension(:) :: x, y, xi
real(rprec), dimension(size(xi)) :: yi
real(rprec) :: dx, t
integer :: i, j, n, k

n = size(x);
k = size(xi);
if (size(y) /= n) then
    call error('interpolate','x and y are not the same size')
end if

j = 1
do i = 1, k
    if ( xi(i) <= x(1) ) then
        yi(i) = y(1)
    else if ( xi(i) >= x(n) ) then
        yi(i) = y(n)
    else
        do while ( xi(i) > x(j + 1) )
            j = j + 1
        end do
        dx = x(j + 1) - x(j)
        t = ( xi(i) - x(j) ) / dx
        yi(i) = (1 - t) * y(j) + t * y(j + 1)
    end if
enddo

end function interpolate_vec

!**********************************************************************
function get_number_of_lines(filename) result(n)
!**********************************************************************
! This function counts the number of lines in a file
use open_file_fid_mod
use messages
implicit none

character(*), intent(in) :: filename
integer :: n, fid, ios
logical :: exst

! Open file
inquire (file = filename, exist = exst)
if (exst) then
	fid = open_file_fid(filename, 'rewind', 'formatted')
else
	call error ('get_number_of_lines', 'file ' // filename // 'does not exist')
end if

! count number of lines
ios = 0
n = 0
do
	read(fid, *, IOstat = ios)
	if (ios /= 0) exit
	n = n + 1
enddo

! Close
close(fid)

end function get_number_of_lines

!**********************************************************************
subroutine read_values_from_file(filename, v1, v2)
!**********************************************************************
! This function counts the number of lines in a file
use open_file_fid_mod
use messages
implicit none

character(*), intent(in) :: filename
real(rprec), dimension(:), intent(inout) :: v1, v2
integer :: fid, ios, i
logical :: exst

! Open file
inquire (file = filename, exist = exst)
if (exst) then
	fid = open_file_fid(filename, 'rewind', 'formatted')
else
	call error ('read_values_from_file', 'file ' // filename // 'does not exist')
end if

! Check that sizes are the same
if ( size(v1) /= size(v2) ) then
	call error ('read_values_from_file', 'v1 and v2 are not the same size')
end if

! read values
do i  = 1, size(v1)
	read(fid, *) v1(i), v2(i)
enddo

! Close
close(fid)

end subroutine read_values_from_file

end module turbines_base
