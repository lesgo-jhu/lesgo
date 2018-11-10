!!
!!  Copyright (C) 2009-2018  Johns Hopkins University
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
subroutine divstress()
!*******************************************************************************
use param, only : BOGUS
use sim_param
use sim_var_3d_m
implicit none

! Calculate derivatives
call txx_var%ddx(dtxxdx)
call txy_var%ddxy(dtxydx, dtxydy)
call tyy_var%ddy(dtyydy)
call txz_var%ddx(dtxzdx)
call txz_var%ddz(dtxzdz)
call tyz_var%ddy(dtyzdy)
call tyz_var%ddz(dtyzdz)
call tzz_var%ddz(dtzzdz)

! Calculate divergence of stress on uv grid
divtx_var%real = dtxxdx%real + dtxydy%real + dtxzdz%real
divty_var%real = dtxydx%real + dtyydy%real + dtyzdz%real

! On w-grid, we need to account for boundary conditions
! Assume dtzzdz at the walls is 0
if (grid%coord == 0) dtzzdz%real(:,:,1) = 0._rprec
if (grid%coord == grid%nproc-1) dtzzdz%real(:,:,grid%nz) = 0._rprec
divtz_var%real = dtxzdx%real + dtyzdy%real + dtzzdz%real

#ifdef PPSAFETYMODE
! Set BOGUS values
divtx_var%real(:,:,0) = BOGUS
divtx_var%real(:,:,grid%nz) = BOGUS
divty_var%real(:,:,0) = BOGUS
divty_var%real(:,:,grid%nz) = BOGUS
divtz_var%real(:,:,0) = BOGUS
#endif

end subroutine
