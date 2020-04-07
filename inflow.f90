!!
!!  Copyright (C) 2016-2020  Johns Hopkins University
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
use param, only : use_inflow
use fringe
#ifdef PPCPS
use concurrent_precursor, only : synchronize_cps, inflow_cond_cps
#endif
#ifdef PPHIT
use hit_inflow
#endif
implicit none

private
public :: inflow_init, apply_inflow

type (fringe_t) :: uniform_fringe

contains

!*******************************************************************************
subroutine inflow_init
!*******************************************************************************
use param, only : fringe_region_end, fringe_region_len

uniform_fringe = fringe_t(fringe_region_end, fringe_region_len)

end subroutine inflow_init

!*******************************************************************************
subroutine apply_inflow
!*******************************************************************************

! Cases for CPS, Isotropic Turbulence and Uniform inflow
#ifdef PPCPS
call synchronize_cps()
if (use_inflow) call inflow_cond_cps()
#elif defined(PPHIT)
if (use_inflow) call inflow_HIT()
#else
if (use_inflow) call uniform_inflow_cond()
#endif

end subroutine apply_inflow

!*******************************************************************************
subroutine uniform_inflow_cond ()
!*******************************************************************************
!  Enforces prescribed inflow condition based on an uniform inflow
!  velocity.
use param, only : nx, ny, nz, inflow_velocity
use sim_param, only : u, v, w
integer :: i, i_w

!--skip istart since we know vel at istart, iend already
do i = 1, uniform_fringe%nx
    i_w = uniform_fringe%iwrap(i)
    u(i_w,1:ny,1:nz) = uniform_fringe%alpha(i) * u(i_w,1:ny,1:nz)              \
        + uniform_fringe%beta(i) * inflow_velocity
    v(i_w,1:ny,1:nz) = uniform_fringe%alpha(i) * v(i_w,1:ny,1:nz)
    w(i_w,1:ny,1:nz) = uniform_fringe%alpha(i) * w(i_w,1:ny,1:nz)
end do

end subroutine uniform_inflow_cond

end module inflow
