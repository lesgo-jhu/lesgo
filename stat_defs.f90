!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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
module stat_defs
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lh
#ifdef PPTURBINES
use turbine_indicator
#endif

save
public

type point_t
    integer :: istart, jstart, kstart, coord
    real(rprec) :: xdiff, ydiff, zdiff
    integer :: fid
end type point_t

type plane_t
    integer :: istart
    real(rprec) :: ldiff
end type plane_t

type zplane_t
    integer :: istart, coord
    real(rprec) :: ldiff
end type zplane_t

type tavg_t
    real(rprec) :: u, v, w, u_w, v_w, w_uv
    real(rprec) :: u2, v2, w2, uv, uw, vw
    real(rprec) :: txx, tyy, tzz, txy, txz, tyz
    real(rprec) :: p, fx, fy, fz
    real(rprec) :: cs_opt2
end type tavg_t

type rs_t
    real(rprec) :: up2, vp2, wp2, upvp, upwp, vpwp
end type rs_t

real(rprec) :: tavg_total_time
! Time between calls of tavg_compute, built by summing dt
real(rprec) :: tavg_dt
! Switch for determining if time averaging has been initialized
logical :: tavg_initialized = .false.

! Types for including wind turbines as drag disks
#ifdef PPTURBINES
! Single turbines
type turbine_t
    real(rprec) :: xloc, yloc, height, dia, thk
    ! term used for volume correction
    real(rprec) :: vol_c
    ! angle CCW(from above) from -x direction [degrees]
    real(rprec) :: theta1
    ! angle above the horizontal, from -x dir [degrees]
    real(rprec) :: theta2
    ! number of nodes associated with each turbine
    integer :: num_nodes
    ! location of turbine center (local k)
    integer :: icp, jcp, kcp
    ! true if the center is in the processor
    logical :: center_in_proc
    ! thrust coefficient
    real(rprec) :: Ct_prime
    ! running time-average of mean disk velocity
    real(rprec) :: u_d, u_d_T
    ! normal force on turbine disk
    real(rprec) :: f_n
    ! (nx,ny,nz) of unit normal for each turbine
    real(rprec), dimension(3) :: nhat
    ! indicator function - weighting of each node
    real(rprec), dimension(10000) :: ind
    ! object to calculate indicator function
    type(turb_ind_func_t) :: turb_ind_func
    ! (i,j,k) of each included node
    integer, dimension(10000,3) :: nodes
    ! search area for nearby nodes
    integer, dimension(6) :: nodes_max
end type turbine_t

! A collection of wind turbines
type wind_farm_t
    type(turbine_t), pointer, dimension(:) :: turbine
end type wind_farm_t

! The wind farm
type(wind_farm_t) :: wind_farm
#endif

! Create types for outputting data (instantaneous or averaged)
type(point_t), allocatable, dimension(:) :: point
type(plane_t), allocatable, dimension(:) :: xplane, yplane
type(zplane_t), allocatable, dimension(:) :: zplane
type(tavg_t), allocatable, dimension(:,:,:) :: tavg
type(rs_t), allocatable, dimension(:,:,:) :: rs

contains

!*******************************************************************************
function rs_compute(a, lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(rs_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % up2 = a % u2 - a % u * a % u
c % vp2 = a % v2 - a % v * a % v
c % wp2 = a % w2 - a % w * a % w
c % upvp = a % uv - a % u * a % v
!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
c % upwp = a % uw - a % u_w * a % w
c % vpwp = a % vw - a % v_w * a % w

end function rs_compute

end module stat_defs
