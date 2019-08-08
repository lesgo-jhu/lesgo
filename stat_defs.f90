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

! Types for including wind turbines as drag disks
#ifdef PPTURBINES
! Single turbines
type turbine_t
    real(rprec) :: xloc, yloc, height, dia, thk
    ! term used for volume correction
    ! real(rprec) :: vol_c
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
    real(rprec), dimension(50000) :: ind
    ! tangential indicator function - weighting of each node
    real(rprec), dimension(50000) :: ind_t
    ! vector of tangential direction
    real(rprec), dimension(50000,3) :: e_theta
    ! object to calculate indicator function
    type(turb_ind_func_t) :: turb_ind_func
    ! (i,j,k) of each included node
    integer, dimension(50000,3) :: nodes
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

contains

end module stat_defs
