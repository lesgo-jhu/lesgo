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

type rs_t
    real(rprec) :: up2, vp2, wp2, upvp, upwp, vpwp
end type rs_t

type spectra_t
    real(rprec), dimension(:), allocatable :: power
    integer :: istart, coord
    real(rprec) :: ldiff
end type spectra_t

real(rprec) :: spectra_total_time
real(rprec) :: tavg_total_time
#ifdef PPOUTPUT_EXTRA
real(rprec) :: tavg_total_time_sgs
#endif
! Time between calls of tavg_compute, built by summing dt
real(rprec) :: tavg_dt
! Switch for determining if time averaging has been initialized
logical :: tavg_initialized = .false.

!  Sums performed over time
type tavg_t
    real(rprec) :: u, v, w, u_w, v_w, w_uv
    real(rprec) :: u2, v2, w2, uv, uw, vw
    ! real(rprec) :: dudz, dvdz
    real(rprec) :: txx, tyy, tzz, txy, txz, tyz
    real(rprec) :: p, fx, fy, fz
    real(rprec) :: cs_opt2
end type tavg_t

!  Sums performed over time (for subgrid variables)
#ifdef PPOUTPUT_EXTRA
type tavg_sgs_t
    real(rprec) :: Nu_t
end type tavg_sgs_t
#endif

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
type(tavg_t), allocatable, dimension(:) :: tavg_zplane

#ifdef PPOUTPUT_EXTRA
type(tavg_sgs_t), allocatable, dimension(:,:,:) :: tavg_sgs
#endif

type(rs_t), allocatable, dimension(:,:,:) :: rs
type(rs_t), allocatable, dimension(:) :: rs_zplane, cnpy_zplane

! Overloaded operators for tavg and rs types
interface operator (.ADD.)
    module procedure tavg_add, tavg_scalar_add, rs_add
end interface

interface operator (.SUB.)
    module procedure tavg_sub, rs_sub
end interface

interface operator (.DIV.)
#ifdef PPOUTPUT_EXTRA
    module procedure tavg_scalar_div, rs_scalar_div, tavg_sgs_scalar_div
#else
    module procedure tavg_scalar_div, rs_scalar_div
#endif
end interface

interface operator (.MUL.)
    module procedure tavg_mul, tavg_scalar_mul
end interface

interface type_set
#ifdef PPOUTPUT_EXTRA
    module procedure tavg_set, rs_set, tavg_sgs_set
#else
    module procedure tavg_set, rs_set
#endif
end interface

interface type_zero_bogus
    module procedure tavg_zero_bogus_2D, tavg_zero_bogus_3D
end interface

contains

!///////////////////////////////////////////////////////////////////////////////
!/// TAVG operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function tavg_add(a, b) result(c)
!*******************************************************************************
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u + b % u
c % v = a % v + b % v
c % w = a % w + b % w
c % u_w  = a % u_w  + b % u_w
c % v_w  = a % v_w  + b % v_w
c % w_uv = a % w_uv + b % w_uv
c % u2 = a % u2 + b % u2
c % v2 = a % v2 + b % v2
c % w2 = a % w2 + b % w2
c % uv = a % uv + b % uv
c % uw = a % uw + b % uw
c % vw = a % vw + b % vw
!c % dudz = a % dudz + b % dudz
!c % dvdz = a % dvdz + b % dvdz
c % txx = a % txx + b % txx
c % tyy = a % tyy + b % tyy
c % tzz = a % tzz + b % tzz
c % txy = a % txy + b % txy
c % txz = a % txz + b % txz
c % tyz = a % tyz + b % tyz
c % p = a % p + b % p
c % fx = a % fx + b % fx
c % fy = a % fy + b % fy
c % fz = a % fz + b % fz
c % cs_opt2 = a % cs_opt2 + b % cs_opt2

end function tavg_add

!*******************************************************************************
function tavg_sub( a, b) result(c)
!*******************************************************************************
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u - b % u
c % v = a % v - b % v
c % w = a % w - b % w
c % u_w  = a % u_w  - b % u_w
c % v_w  = a % v_w  - b % v_w
c % w_uv = a % w_uv - b % w_uv
c % u2 = a % u2 - b % u2
c % v2 = a % v2 - b % v2
c % w2 = a % w2 - b % w2
c % uv = a % uv - b % uv
c % uw = a % uw - b % uw
c % vw = a % vw - b % vw
!c % dudz = a % dudz - b % dudz
!c % dvdz = a % dvdz - b % dvdz
c % txx = a % txx - b % txx
c % tyy = a % tyy - b % tyy
c % tzz = a % tzz - b % tzz
c % txy = a % txy - b % txy
c % txz = a % txz - b % txz
c % tyz = a % tyz - b % tyz
c % p = a % p - b % p
c % fx = a % fx - b % fx
c % fy = a % fy - b % fy
c % fz = a % fz - b % fz
c % cs_opt2 = a % cs_opt2 - b % cs_opt2

end function tavg_sub

!*******************************************************************************
function tavg_scalar_add( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u + b
c % v = a % v + b
c % w = a % w + b
c % u_w  = a % u_w  + b
c % v_w  = a % v_w  + b
c % w_uv = a % w_uv + b
c % u2 = a % u2 + b
c % v2 = a % v2 + b
c % w2 = a % w2 + b
c % uv = a % uv + b
c % uw = a % uw + b
c % vw = a % vw + b
!c % dudz = a % dudz + b
!c % dvdz = a % dvdz + b
c % txx = a % txx + b
c % tzz = a % tzz + b
c % tyy = a % tyy + b
c % txy = a % txy + b
c % txz = a % txz + b
c % tyz = a % tyz + b
c % p = a % p + b
c % fx = a % fx + b
c % fy = a % fy + b
c % fz = a % fz + b
c % cs_opt2 = a % cs_opt2 + b

end function tavg_scalar_add

!*******************************************************************************
subroutine tavg_zero_bogus_2D( c )
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), dimension(:,:), intent(inout) :: c

c % txx = 0._rprec
c % tyy = 0._rprec
c % tzz = 0._rprec
c % txy = 0._rprec
c % txz = 0._rprec
c % tyz = 0._rprec
c % fx = 0._rprec
c % fy = 0._rprec
c % fz = 0._rprec

end subroutine tavg_zero_bogus_2D

!*******************************************************************************
subroutine tavg_zero_bogus_3D( c )
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), dimension(:,:,:), intent(inout) :: c

c % txx = 0._rprec
c % tyy = 0._rprec
c % tzz = 0._rprec
c % txy = 0._rprec
c % txz = 0._rprec
c % tyz = 0._rprec
c % fx = 0._rprec
c % fy = 0._rprec
c % fz = 0._rprec

end subroutine tavg_zero_bogus_3D

!*******************************************************************************
function tavg_scalar_div( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u / b
c % v = a % v / b
c % w = a % w / b
c % u_w  = a % u_w  / b
c % v_w  = a % v_w  / b
c % w_uv = a % w_uv / b
c % u2 = a % u2 / b
c % v2 = a % v2 / b
c % w2 = a % w2 / b
c % uv = a % uv / b
c % uw = a % uw / b
c % vw = a % vw / b
!c % dudz = a % dudz / b
!c % dvdz = a % dvdz / b
c % txx = a % txx / b
c % tyy = a % tyy / b
c % tzz = a % tzz / b
c % txy = a % txy / b
c % txz = a % txz / b
c % tyz = a % tyz / b
c % p = a % p / b
c % fx = a % fx / b
c % fy = a % fy / b
c % fz = a % fz / b
c % cs_opt2 = a % cs_opt2 / b

end function tavg_scalar_div

!*******************************************************************************
function tavg_mul( a, b) result(c)
!*******************************************************************************
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u * b % u
c % v = a % v * b % v
c % w = a % w * b % w
c % u_w  = a % u_w  * b % u_w
c % v_w  = a % v_w  * b % v_w
c % w_uv = a % w_uv * b % w_uv
c % u2 = a % u2 * b % u2
c % v2 = a % v2 * b % v2
c % w2 = a % w2 * b % w2
c % uv = a % uv * b % uv
c % uw = a % uw * b % uw
c % vw = a % vw * b % vw
!c % dudz = a % dudz * b % dudz
!c % dvdz = a % dvdz * b % dvdz
c % txx = a % txx * b % txx
c % tyy = a % tyy * b % tyy
c % tzz = a % tzz * b % tzz
c % txy = a % txy * b % txy
c % txz = a % txz * b % txz
c % tyz = a % tyz * b % tyz
c % p = a % p * b % p
c % fx = a % fx * b % fx
c % fy = a % fy * b % fy
c % fz = a % fz * b % fz
c % cs_opt2 = a % cs_opt2 * b % cs_opt2

end function tavg_mul

!*******************************************************************************
function tavg_scalar_mul( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u * b
c % v = a % v * b
c % w = a % w * b
c % u_w  = a % u_w  * b
c % v_w  = a % v_w  * b
c % w_uv = a % w_uv * b
c % u2 = a % u2 * b
c % v2 = a % v2 * b
c % w2 = a % w2 * b
c % uv = a % uv * b
c % uw = a % uw * b
c % vw = a % vw * b
!c % dudz = a % dudz * b
!c % dvdz = a % dvdz * b
c % txx = a % txx * b
c % tyy = a % tyy * b
c % tzz = a % tzz * b
c % txy = a % txy * b
c % txz = a % txz * b
c % tyz = a % tyz * b
c % p = a % p * b
c % fx = a % fx * b
c % fy = a % fy * b
c % fz = a % fz * b
c % cs_opt2 = a % cs_opt2 * b

end function tavg_scalar_mul

#ifdef PPOUTPUT_EXTRA
!*******************************************************************************
function tavg_sgs_scalar_div( a, b ) result(c)
!*******************************************************************************
use types, only : rprec
implicit none

type(tavg_sgs_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_sgs_t) :: c

!c % Tn = a % Tn / b
c % Nu_t = a % Nu_t / b
!c % F_LM = a % F_LM / b
!c % F_MM = a % F_MM / b
!c % F_QN = a % F_QN / b
!c % F_NN = a % F_NN / b
!c % ee_now = a % ee_now / b
!#ifdef PPDYN_TN
!c % F_ee2 = a % F_ee2 / b
!c % F_deedt2 = a % F_deedt2 / b
!#endif

end function tavg_sgs_scalar_div
#endif

!*******************************************************************************
function tavg_interp_to_uv_grid( a ) result(c)
!*******************************************************************************
use param, only: lbz
use functions, only : interp_to_uv_grid
implicit none

type(tavg_t), dimension(:,:,lbz:), intent(in) :: a
type(tavg_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz:ubz))

c = a

c % fz = interp_to_uv_grid(a % fz, lbz )
c % w  = interp_to_uv_grid(a % w,lbz)
c % w2 = interp_to_uv_grid(a % w2,lbz)

end function tavg_interp_to_uv_grid

!*******************************************************************************
function tavg_interp_to_w_grid( a ) result(c)
!*******************************************************************************
use param, only: lbz
use functions, only : interp_to_w_grid
implicit none

type(tavg_t), dimension(:,:,lbz:), intent(in) :: a
type(tavg_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz:ubz))

c = a

c % txx =  interp_to_w_grid( a % txx, lbz )
c % tyy =  interp_to_w_grid( a % tyy, lbz )
c % tzz =  interp_to_w_grid( a % tzz, lbz )
c % txy =  interp_to_w_grid( a % txy, lbz )

c % p = interp_to_w_grid( a % p, lbz )
c % fx = interp_to_w_grid( a % fx, lbz )
c % fy = interp_to_w_grid( a % fy, lbz )

end function tavg_interp_to_w_grid

!///////////////////////////////////////////////////////////////////////////////
!/// RS operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function rs_add( a, b) result(c)
!*******************************************************************************
implicit none

type(rs_t), intent(in) :: a, b
type(rs_t) :: c

c % up2 = a % up2 + b % up2
c % vp2 = a % vp2 + b % vp2
c % wp2 = a % wp2 + b % wp2
c % upvp = a % upvp + b % upvp
c % upwp = a % upwp + b % upwp
c % vpwp = a % vpwp + b % vpwp

end function rs_add

!*******************************************************************************
function rs_sub( a, b) result(c)
!*******************************************************************************
implicit none

type(rs_t), intent(in) :: a, b
type(rs_t) :: c

c % up2 = a % up2 - b % up2
c % vp2 = a % vp2 - b % vp2
c % wp2 = a % wp2 - b % wp2
c % upvp = a % upvp - b % upvp
c % upwp = a % upwp - b % upwp
c % vpwp = a % vpwp - b % vpwp

end function rs_sub

!*******************************************************************************
function rs_scalar_div( a, b) result(c)
!*******************************************************************************
implicit none

type(rs_t), intent(in) :: a
real(rprec), intent(in) :: b
type(rs_t) :: c

c % up2 = a % up2 / b
c % vp2 = a % vp2 / b
c % wp2 = a % wp2 / b
c % upvp = a % upvp / b
c % upwp = a % upwp / b
c % vpwp = a % vpwp / b

end function rs_scalar_div

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral RS operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function rs_compute( a , lbz2) result(c)
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
c % upwp = a % uw - a % u_w * a % w   !!pj
c % vpwp = a % vw - a % v_w * a % w   !!pj

end function rs_compute

!*******************************************************************************
function cnpy_tavg_mul( a ) result(c)
!*******************************************************************************
!
! This performs one set of multiplication for the canopy stresses
!
implicit none

type(tavg_t), intent(in) :: a
type(rs_t) :: c

c % up2 = a % u * a % u
c % vp2 = a % v * a % v
c % wp2 = a % w * a % w
c % upvp = a % u * a % v
c % upwp = a % u * a % w
c % vpwp = a % v * a % w

end function cnpy_tavg_mul

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral TAVG operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
subroutine tavg_set( c, a )
!*******************************************************************************
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_t), intent(out) :: c

c % u = a
c % v = a
c % w = a
c % u_w  = a
c % v_w  = a
c % w_uv = a
c % u2 = a
c % v2 = a
c % w2 = a
c % uv = a
c % uw = a
c % vw = a
!c % dudz = a
!c % dvdz = a
c % txx = a
c % tyy = a
c % tzz = a
c % txy = a
c % txz = a
c % tyz = a
c % fx = a
c % fy = a
c % fz = a
c % cs_opt2 = a

end subroutine tavg_set

#ifdef PPOUTPUT_EXTRA
!*******************************************************************************
subroutine tavg_sgs_set( c, a )
!*******************************************************************************
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_sgs_t), intent(out) :: c

!c % Tn =  a
c % Nu_t =  a
!c % F_LM =  a
!c % F_MM =  a
!c % F_QN =  a
!c % F_NN =  a
!c % ee_now = a
!#ifdef PPDYN_TN
!c % F_ee2 = a
!c % F_deedt2 = a
!#endif

end subroutine tavg_sgs_set
#endif

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral RS subroutines
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
subroutine rs_set( c, a )
!*******************************************************************************
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(rs_t), intent(out) :: c

c % up2 = a
c % vp2 = a
c % wp2 = a
c % upvp = a
c % upwp = a
c % vpwp = a

end subroutine rs_set

end module stat_defs
