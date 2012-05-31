module turbines_base
use types, only:rprec
use stat_defs, only:wind_farm_t
$if ($MPI)
  use mpi_defs
$endif

implicit none

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

real(rprec) :: Ct_prime     ! thrust coefficient (default 1.33)
real(rprec) :: Ct_noprime   ! thrust coefficient (default 0.75)

real(rprec) :: T_avg_dim    ! disk-avg time scale in seconds (default 600)

real(rprec) :: alpha        ! filter size as multiple of grid spacing
integer :: trunc            ! Gaussian filter truncated after this many gridpoints
real(rprec) :: filter_cutoff  ! indicator function only includes values above this threshold

logical :: turbine_cumulative_time ! Used to read in the disk averaged velocities of the turbines

integer(rprec) :: tbase     ! Number of timesteps between the output
 
! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: sx           ! spacing in the x-direction, multiple of (mean) diameter
real(rprec) :: sy           ! spacing in the y-direction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets the values for wind_farm_t based on values
!   read from the input file.

subroutine turbines_base_init()
use param, only: L_x, L_y, dx, dy, dz, pi, z_i
implicit none

integer :: i, j, k
real(rprec) :: sxx, syy, shift_base, const

! set turbine parameters
! turbines are numbered as follows:
!   #1 = turbine nearest (x,y)=(0,0)
!   #2 = next turbine in the y-direction, etc. (go along rows)

    ! Allocate wind turbine array derived type
    nloc = num_x*num_y      !number of turbines (locations) 
    nullify(wind_farm_t%turbine_t)
    allocate(wind_farm_t%turbine_t(nloc))

    ! Non-dimensionalize length values by z_i
    dia_all = dia_all / z_i
    height_all = height_all / z_i
    thk_all = thk_all / z_i
   
    ! Resize thickness capture at least on plane of gridpoints
    thk_all = max ( thk_all, dx*1.01 )

    ! Set baseline values for size
    wind_farm_t%turbine_t(:)%height = height_all
    wind_farm_t%turbine_t(:)%dia = dia_all
    wind_farm_t%turbine_t(:)%thk = thk_all                      
    wind_farm_t%turbine_t(:)%vol_c =  dx*dy*dz/(pi/4.*(dia_all)**2 * thk_all)        

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
        wind_farm_t%turbine_t(k)%xloc = sxx*real(2*i-1)/2
        wind_farm_t%turbine_t(k)%yloc = syy*real(2*j-1)/2
        k = k + 1
      enddo
    enddo

    ! HERE PLACE TURBINES (x,y-positions) BASED ON 'ORIENTATION' FLAG
    if (orientation.eq.1) then
    ! Evenly-spaced, not staggered
    !  Use baseline as set above       
 
    elseif (orientation.eq.2) then
    ! Evenly-spaced, horizontally staggered only
      ! Shift each row according to stag_perc
      do i = 2, num_x
        do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
          shift_base = syy * stag_perc/100.
          wind_farm_t%turbine_t(k)%yloc = mod( wind_farm_t%turbine_t(k)%yloc + (i-1)*shift_base , L_y )
        enddo
      enddo
 
    elseif (orientation.eq.3) then 
    ! Evenly-spaced, only vertically staggered (by rows)
      ! Make even rows taller
      do i = 2, num_x, 2
        do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
          wind_farm_t%turbine_t(k)%height = height_all*(1.+stag_perc/100.)
        enddo
      enddo
      ! Make odd rows shorter
      do i = 1, num_x, 2
        do k = 1+num_y*(i-1), num_y*i         ! these are the numbers for turbines in row i
          wind_farm_t%turbine_t(k)%height = height_all*(1.-stag_perc/100.)
        enddo
      enddo
 
    elseif (orientation.eq.4) then        
    !Evenly-spaced, only vertically staggered, checkerboard pattern
      k = 1
      do i = 1, num_x 
        do j = 1, num_y
          const = 2.*mod(real(i+j),2.)-1.  ! this should alternate between 1, -1
          wind_farm_t%turbine_t(k)%height = height_all*(1.+const*stag_perc/100.)
          k = k + 1
        enddo
      enddo

    elseif (orientation.eq.5) then        
    !Aligned, but shifted forward for efficient use of simulation space during CPS runs

      ! Usual placement is baseline as set above

      ! Shift the turbines forward
      k=1
      do i = 1, num_x
        do j = 1, num_y
          wind_farm_t%turbine_t(k)%xloc=wind_farm_t%turbine_t(k)%xloc -wind_farm_t%turbine_t(1)%xloc/2
          k=k+1
        enddo
      enddo
    
    endif
        
    !orientation (angles)
        wind_farm_t%turbine_t(:)%theta1 = theta1_all
        wind_farm_t%turbine_t(:)%theta2 = theta2_all

 
end subroutine turbines_base_init

end module turbines_base
