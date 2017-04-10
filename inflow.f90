!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
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
module inflow
!*******************************************************************************
use types, only : rprec
use param, only : inflow_cond, use_cps
use fringe_region
use hit_inflow
use concurrent_precursor
implicit none

save
private

public :: inflow_init, inflow_apply

contains

!*******************************************************************************
subroutine inflow_init()
!*******************************************************************************
use param, only : inflow_velocity, fringe_region_end, fringe_region_len
use param, only : recycl_region_end
use messages
implicit none

! For all non-periodic inflow conditions, the fringe region is used
if (inflow_cond > 1) then
    fringe = fringe_region_t(fringe_region_end, fringe_region_len)
end if

if (inflow_cond == 1) then
    ! Periodic boundary conditions (default)
    write(*,*) "Using periodic boundary conditions"
else if (inflow_cond == 2) then
    ! Uniform inflow
    ! Specify velocities in the fringe region
    write(*,*) "Initializing uniform inflow"
    fringe%u = inflow_velocity
    fringe%v = 0._rprec
    fringe%w = 0._rprec
else if (inflow_cond == 3) then
    ! Homogenous isotropic turbulence (HIT)
    write(*,*) "Initializing homogenous isotropic turbulence (HIT) inflow"
    call HIT_init()
else if (inflow_cond == 4) then 
    ! Shifted periodic boundary conditions    
    write(*,*) "Initializing shifted periodic boundary conditions"
    recycl = fringe_region_t(recycl_region_end, fringe%n)
else if (inflow_cond == 5) then 
    ! Using concurrent precursor inflow
    write(*,*) "Using periodic boundary conditions"
else
    call error("inflow.inflow_init", "Invalid inflow condition.")
end if

! initialize CPS
if (use_cps) then
#ifdef PPMPI
    call cps_init()
#else
    call error("inflow.inflow_init", "CPS only allowed in MPI mode.")
#endif
end if

end subroutine inflow_init

!*******************************************************************************
subroutine inflow_apply()
!*******************************************************************************
use param, only : coord, nx
use sim_param, only : u
implicit none
integer :: i

if (use_cps) then
    call synchronize_cps()
end if

if (inflow_cond == 2) then
    ! Uniform inflow
    call fringe%apply_vel()
else if (inflow_cond == 3) then
    ! Homogenous isotropic turbulence (HIT)
    call inflow_HIT()
else if (inflow_cond == 4) then
    ! Shifted periodic boundary conditions
    call fringe%sample_vel()
    recycl%u = fringe%u
    recycl%v = fringe%v
    recycl%w = fringe%w
    call recycl%apply_vel()
else if (inflow_cond == 4) then 
    ! Concurrent precursor
    call fringe%apply_vel()
end if

end subroutine inflow_apply

end module inflow