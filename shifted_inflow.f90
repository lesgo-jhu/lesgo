!!
!!  Copyright (C) 2020  Johns Hopkins University
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
module shifted_inflow
!*******************************************************************************
use types, only : rprec
use fringe

implicit none

private
public shifted_inflow_init, inflow_shifted

type(fringe_t) :: sample_fringe
type(fringe_t) :: apply_fringe

real(rprec), allocatable, dimension(:,:,:) :: u_s, v_s, w_s
#ifdef PPSCALARS
real(rprec), allocatable, dimension(:,:,:) :: theta_s
#endif

contains
!*******************************************************************************
subroutine shifted_inflow_init
!*******************************************************************************
use param, only : fringe_region_end, fringe_region_len, sampling_region_end
use param, only : shift_n, ny, nz

sample_fringe = fringe_t(sampling_region_end, fringe_region_len)
apply_fringe = fringe_t(fringe_region_end, fringe_region_len)

! Allocate the sample block
allocate(u_s(sample_fringe%nx, ny, nz ))
allocate(v_s(sample_fringe%nx, ny, nz ))
allocate(w_s(sample_fringe%nx, ny, nz ))
#ifdef PPSCALARS
allocate(theta_s(sample_fringe%nx, ny, nz))
#endif

! Only allow positive shifts than are less than ny/2
shift_n = modulo(abs(shift_n), ny)
if (shift_n > ny/2) shift_n = ny-shift_n
if (shift_n == 0) shift_n = 1

end subroutine shifted_inflow_init

!*******************************************************************************
subroutine inflow_shifted
!*******************************************************************************
use param, only : shift_n, ny, nz, coord
use sim_param, only : u, v, w
#ifdef PPSCALARS
use scalars, only : theta
#endif
integer :: i, i_w

! Sample and shift velocity
u_s(:,shift_n+1:ny,:) = u(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
u_s(:,1:shift_n,:) = u(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
v_s(:,shift_n+1:ny,:) = v(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
v_s(:,1:shift_n,:) = v(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
w_s(:,shift_n+1:ny,:) = w(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
w_s(:,1:shift_n,:) = w(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
#ifdef PPSCALARS
theta_s(:,shift_n+1:ny,:) = theta(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
theta_s(:,1:shift_n,:) = theta(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
#endif

! Apply inflow conditions
do i = 1, apply_fringe%nx
    i_w = apply_fringe%iwrap(i)
    u(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * u(i_w,1:ny,1:nz)                  &
        + apply_fringe%beta(i) * u_s(i,1:ny,1:nz)
    v(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * v(i_w,1:ny,1:nz)                  &
        + apply_fringe%beta(i) * v_s(i,1:ny,1:nz)
    w(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * w(i_w,1:ny,1:nz)                  &
        + apply_fringe%beta(i) * w_s(i,1:ny,1:nz)
#ifdef PPSCALARS
    theta(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * theta(i_w,1:ny,1:nz)          &
        + apply_fringe%beta(i) * theta_s(i,1:ny,1:nz)
#endif
end do

end subroutine inflow_shifted

end module shifted_inflow
