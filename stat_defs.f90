!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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

!**********************************************************************
module stat_defs
!**********************************************************************
use types, only : rprec
use param, only : nx,ny,nz,lh

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

type theta_stat_t     !Eshwan:  - temperature statistics
  real(rprec) :: wpthetap, thetap2
  real(rprec) :: total_tflux
end type theta_stat_t

type sal_stat_t     !Eshwan:  - salinity statistics
  real(rprec) :: wpsalp, salp2
  real(rprec) :: total_sflux
end type sal_stat_t

type sigma_t 
     real(rprec) :: sigma_theta, sigma_sal
end type sigma_t

type spectra_t
  real(rprec), dimension(:), allocatable :: power
  integer :: istart, coord
  real(rprec) :: ldiff 
end type spectra_t

real(rprec) :: spectra_total_time
real(rprec) :: tavg_total_time
!$if($OUTPUT_EXTRA) !- Eshwan
real(rprec) :: tavg_total_time_sgs
!$endif !- Eshwan
! Time between calls of tavg_compute, built by summing dt
real(rprec) :: tavg_dt
! Switch for determining if time averaging has been initialized
logical :: tavg_initialized = .false.

! Time between calls of spectra_compute, built by summing dt
real(rprec) :: spectra_dt
! Switch for determining if time averaging has been initialized
logical :: spectra_initialized = .false.

!  Sums performed over time
type tavg_t
  real(rprec) :: u, v, w
  real(rprec) :: u_uv, v_uv, w_uv
  real(rprec) :: u2, v2, w2, uv, uw, vw
  real(rprec) :: dudz, dvdz, dudz2             !Eshwan
  real(rprec) :: txx, tyy, tzz, txy, txz, tyz
  real(rprec) :: fx, fy, fz
  real(rprec) :: cs_opt2, Nu_t  
  real(rprec) :: divtx, divty, divtz
  real(rprec) :: buoyancy
end type tavg_t
 
type theta_avg_t    !Eshwan
  real(rprec) :: theta
  real(rprec) :: theta_w
  real(rprec) :: wtheta    
  real(rprec) :: theta2
  real(rprec) :: dTdx, dTdy, dTdz
  real(rprec) :: sgs_tflux
end type theta_avg_t
 
type sal_avg_t    !Eshwan
  real(rprec) :: sal
  real(rprec) :: sal_w
  real(rprec) :: wsal    
  real(rprec) :: sal2
  real(rprec) :: dSdx, dSdy, dSdz    
  real(rprec) :: sgs_sflux
end type sal_avg_t

!  Sums performed over time (for subgrid variables)
!$if($OUTPUT_EXTRA) 
type tavg_sgs_t
!  real(rprec) :: Tn, Nu_t !Eshwan
  real(rprec) :: Tn, Beta
  real(rprec) :: F_LM, F_MM, F_QN, F_NN
  real(rprec) :: ee_now
  real(rprec) :: cs2_clips
  $if ($DYN_TN)
  real(rprec) :: F_ee2, F_deedt2
  $endif
end type tavg_sgs_t 

type tavg_scalar_sgs_t !Eshwan
  real(rprec) :: Ds_opt2_t, Ds_opt2_s 
  real(rprec) :: kappa_tt, kappa_ts 
  real(rprec) :: s_Tn_t, s_Tn_s
  real(rprec) :: I_LM_t, I_MM_t, I_QN_t, I_NN_t
  real(rprec) :: I_LM_s, I_MM_s, I_QN_s, I_NN_s
  real(rprec) :: s_Beta_t, s_Beta_s
  real(rprec) :: ds2_clips_t, ds2_clips_s
end type tavg_scalar_sgs_t
!$endif 

! Types for including wind-turbines as drag disks
$if ($TURBINES)
! Single turbines
type turbine_t
  real(rprec) :: xloc, yloc, height, dia, thk
  real(rprec) :: vol_c                        ! term used for volume correction  
  real(rprec) :: theta1                       ! angle CCW(from above) from -x direction [degrees]
  real(rprec) :: theta2                       ! angle above the horizontal, from -x dir [degrees]
  real(rprec), dimension(3) :: nhat           ! (nx,ny,nz) of unit normal for each turbine
  integer :: num_nodes                        ! number of nodes associated with each turbine
  integer, dimension(5000,3) :: nodes         ! (i,j,k) of each included node
  integer, dimension(6) :: nodes_max          ! search area for nearby nodes
  real(rprec) :: u_d, u_d_T                   ! running time-average of mean disk velocity
  real(rprec) :: f_n                          ! normal force on turbine disk
  real(rprec), dimension(5000) :: ind         ! indicator function - weighting of each node
end type turbine_t

! A collection of wind-turbines
type wind_farm_t
  type(turbine_t), pointer, dimension(:) :: turbine
end type wind_farm_t
    
type(wind_farm_t) :: wind_farm
$endif

! Histogram (single)
type hist_t
    real(rprec) :: bmin, bmax, db             ! bin min, max, and spacing
    integer :: nbins                          ! number of bins
    real(rprec), allocatable, dimension(:) :: bins  ! bin centers
    real(rprec), allocatable, dimension(:) :: vals  ! count for each bin (may be normalized)
end type hist_t

! Collection of histograms (one for each zplane) for a single variable
type hist_zplanes_t  
    integer, allocatable, dimension(:) :: coord         ! processor where this plane exists
    integer, allocatable, dimension(:) :: istart        ! nearest node below plane (for interpolation)
    real(rprec), allocatable, dimension(:) :: ldiff     ! distance from istart to plane (for interpolation)
    type(hist_t), allocatable, dimension(:) :: hist     ! the histograms for each plane
end type hist_zplanes_t

! Create histogram groups here 
type(hist_zplanes_t) :: HISTcs2   ! SGS coefficient, squared
type(hist_zplanes_t) :: HISTtn    ! Lagrangian time scale
type(hist_zplanes_t) :: HISTnu    ! Eddy viscosity
type(hist_zplanes_t) :: HISTee    ! Error in SGS model

! Create types for outputting data (instantaneous or averaged)
type(point_t), allocatable, dimension(:) :: point
type(plane_t), allocatable, dimension(:) :: xplane, yplane
type(zplane_t), allocatable, dimension(:) :: zplane

type(tavg_t), allocatable, dimension(:,:,:) :: tavg
type(tavg_t), allocatable, dimension(:) :: tavg_zplane
type(theta_avg_t), allocatable, dimension(:,:,:) :: theta_avg        !Eshwan
type(theta_avg_t), allocatable, dimension(:) :: theta_avg_zplane     !Eshwan
type(sal_avg_t), allocatable, dimension(:,:,:) :: sal_avg        !Eshwan
type(sal_avg_t), allocatable, dimension(:) :: sal_avg_zplane     !Eshwan

real(rprec), allocatable, dimension(:,:) :: ustar_avg 
real(rprec), allocatable, dimension(:,:) :: t_flux_avg
real(rprec), allocatable, dimension(:,:) :: tstar_avg
real(rprec), allocatable, dimension(:,:) :: sal_flux_avg
real(rprec), allocatable, dimension(:,:) :: sstar_avg
real(rprec) :: ustar_avg_zplane
real(rprec) :: t_flux_avg_zplane
real(rprec) :: tstar_avg_zplane
real(rprec) :: sal_flux_avg_zplane
real(rprec) :: sstar_avg_zplane

!$if ($OUTPUT_EXTRA) !- Eshwan
type(tavg_sgs_t), allocatable, dimension(:,:,:) :: tavg_sgs
type(tavg_sgs_t), allocatable, dimension(:) :: tavg_sgs_zplane       !Eshwan
type(tavg_scalar_sgs_t), allocatable, dimension(:,:,:) :: tavg_scalar_sgs
type(tavg_scalar_sgs_t), allocatable, dimension(:) :: tavg_scalar_sgs_zplane       !Eshwan
!$endif !- Eshwan

type(rs_t), allocatable, dimension(:,:,:) :: rs
type(rs_t), allocatable, dimension(:) :: rs_zplane, cnpy_zplane
type(theta_stat_t), allocatable, dimension(:,:,:) :: theta_stat !Eshwan
type(theta_stat_t), allocatable, dimension(:) :: theta_stat_zplane !Eshwan
type(sal_stat_t), allocatable, dimension(:,:,:) :: sal_stat !Eshwan
type(sal_stat_t), allocatable, dimension(:) :: sal_stat_zplane !Eshwan
type(sigma_t), allocatable, dimension(:) :: sigma_avg_zplane !Eshwan
type(spectra_t), allocatable, dimension(:) :: spectra

! Overloaded operators for tavg and rs types
INTERFACE OPERATOR (.ADD.)
  MODULE PROCEDURE tavg_add, tavg_scalar_add, rs_add, tavg_sgs_add, tavg_scalar_sgs_add, & !Eshwan
                   theta_avg_add, theta_stat_add, sal_avg_add, sal_stat_add                 !Eshwan
END INTERFACE

INTERFACE OPERATOR (.SUB.)
  MODULE PROCEDURE tavg_sub, rs_sub
END INTERFACE

INTERFACE OPERATOR (.DIV.)
    MODULE PROCEDURE tavg_scalar_div, rs_scalar_div, tavg_sgs_scalar_div, tavg_scalar_sgs_scalar_div, &   !Eshwan
                     theta_avg_scalar_div, theta_stat_scalar_div, sigma_avg_zplane_scalar_div, &   !Eshwan
                     sal_avg_scalar_div, sal_stat_scalar_div !Eshwan
END INTERFACE

INTERFACE OPERATOR (.MUL.)
  MODULE PROCEDURE tavg_mul, tavg_scalar_mul
END INTERFACE

INTERFACE type_set
    MODULE PROCEDURE tavg_set, rs_set, tavg_sgs_set, tavg_scalar_sgs_set, theta_avg_set, &
                     theta_stat_set, sal_avg_set, sal_stat_set, sigma_avg_zplane_set  !Eshwan
END INTERFACE

INTERFACE type_zero_bogus
  MODULE PROCEDURE tavg_zero_bogus_2D, tavg_zero_bogus_3D
END INTERFACE

INTERFACE hist_binit
  MODULE PROCEDURE hist_binit_1D, hist_binit_2D, hist_binit_3D
END INTERFACE

contains

!//////////////////////////////////////////////////////////////////////
!/////////////////// TAVG OPERATORS ///////////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u + b % u
c % v = a % v + b % v
c % w = a % w + b % w
c % u_uv = a % u_uv + b % u_uv
c % v_uv = a % v_uv + b % v_uv
c % w_uv = a % w_uv + b % w_uv
c % u2 = a % u2 + b % u2
c % v2 = a % v2 + b % v2
c % w2 = a % w2 + b % w2
c % uv = a % uv + b % uv
c % uw = a % uw + b % uw
c % vw = a % vw + b % vw
c % dudz = a % dudz + b % dudz
c % dvdz = a % dvdz + b % dvdz
c % dudz2 = a % dudz2 + b % dudz2
c % txx = a % txx + b % txx
c % tyy = a % tyy + b % tyy
c % tzz = a % tzz + b % tzz
c % txy = a % txy + b % txy
c % txz = a % txz + b % txz
c % tyz = a % tyz + b % tyz
c % fx = a % fx + b % fx
c % fy = a % fy + b % fy
c % fz = a % fz + b % fz
c % cs_opt2 = a % cs_opt2 + b % cs_opt2
c % Nu_t = a % Nu_t + b % Nu_t !Eshwan
c % divtx = a % divtx + b % divtx
c % divty = a % divty + b % divty
c % divtz = a % divtz + b % divtz
c % buoyancy = a % buoyancy + b % buoyancy 

return
end function tavg_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_sub( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u - b % u
c % v = a % v - b % v
c % w = a % w - b % w
c % u_uv = a % u_uv - b % u_uv
c % v_uv = a % v_uv - b % v_uv
c % w_uv = a % w_uv - b % w_uv
c % u2 = a % u2 - b % u2
c % v2 = a % v2 - b % v2
c % w2 = a % w2 - b % w2
c % uv = a % uv - b % uv
c % uw = a % uw - b % uw
c % vw = a % vw - b % vw 
c % dudz = a % dudz - b % dudz
c % dvdz = a % dvdz - b % dvdz
c % dudz2 = a % dudz2 - b % dudz2
c % txx = a % txx - b % txx
c % tyy = a % tyy - b % tyy
c % tzz = a % tzz - b % tzz
c % txy = a % txy - b % txy
c % txz = a % txz - b % txz
c % tyz = a % tyz - b % tyz
c % fx = a % fx - b % fx
c % fy = a % fy - b % fy
c % fz = a % fz - b % fz
c % cs_opt2 = a % cs_opt2 - b % cs_opt2
c % Nu_t = a % Nu_t - b % Nu_t !Eshwan
c % divtx = a % divtx - b % divtx
c % divty = a % divty - b % divty
c % divtz = a % divtz - b % divtz
c % buoyancy = a % buoyancy - b % buoyancy 

return
end function tavg_sub

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_add( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u + b
c % v = a % v + b
c % w = a % w + b
c % u_uv = a % u_uv + b
c % v_uv = a % v_uv + b
c % w_uv = a % w_uv + b
c % u2 = a % u2 + b
c % v2 = a % v2 + b
c % w2 = a % w2 + b
c % uv = a % uv + b
c % uw = a % uw + b
c % vw = a % vw + b
c % dudz = a % dudz + b
c % dvdz = a % dvdz + b
c % dudz2 = a % dudz2 + b
c % txx = a % txx + b
c % tzz = a % tzz + b
c % tyy = a % tyy + b
c % txy = a % txy + b
c % txz = a % txz + b
c % tyz = a % tyz + b
c % fx = a % fx + b
c % fy = a % fy + b
c % fz = a % fz + b
c % cs_opt2 = a % cs_opt2 + b
c % Nu_t = a % Nu_t + b !Eshwan
c % divtx = a % divtx + b 
c % divty = a % divty + b
c % divtz = a % divtz + b 
c % buoyancy = a % buoyancy + b 

return
end function tavg_scalar_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_zero_bogus_2D( c )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

return
end subroutine tavg_zero_bogus_2D

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_zero_bogus_3D( c )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

return
end subroutine tavg_zero_bogus_3D


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_div( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u / b
c % v = a % v / b
c % w = a % w / b
c % u2 = a % u2 / b
c % u_uv = a % u_uv / b
c % v_uv = a % v_uv / b
c % w_uv = a % w_uv / b
c % v2 = a % v2 / b
c % w2 = a % w2 / b
c % uv = a % uv / b
c % uw = a % uw / b
c % vw = a % vw / b
c % dudz = a % dudz / b
c % dvdz = a % dvdz / b
c % dudz2 = a % dudz2 / b
c % txx = a % txx / b
c % tyy = a % tyy / b
c % tzz = a % tzz / b
c % txy = a % txy / b
c % txz = a % txz / b
c % tyz = a % tyz / b
c % fx = a % fx / b
c % fy = a % fy / b
c % fz = a % fz / b
c % cs_opt2 = a % cs_opt2 / b
c % Nu_t = a % Nu_t / b !Eshwan
c % divtx = a % divtx / b 
c % divty = a % divty / b
c % divtz = a % divtz / b 
c % buoyancy = a % buoyancy / b 

return
end function tavg_scalar_div

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_mul( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(tavg_t), intent(in) :: a, b
type(tavg_t) :: c

c % u = a % u * b % u
c % v = a % v * b % v
c % w = a % w * b % w
c % u2 = a % u2 * b % u2
c % u_uv = a % u_uv * b % u_uv
c % v_uv = a % v_uv * b % v_uv
c % w_uv = a % w_uv * b % w_uv
c % v2 = a % v2 * b % v2
c % w2 = a % w2 * b % w2
c % uv = a % uv * b % uv
c % uw = a % uw * b % uw
c % vw = a % vw * b % vw
c % dudz = a % dudz * b % dudz
c % dvdz = a % dvdz * b % dvdz
c % dudz2 = a % dudz2 * b % dudz2
c % txx = a % txx * b % txx
c % tyy = a % tyy * b % tyy
c % tzz = a % tzz * b % tzz
c % txy = a % txy * b % txy
c % txz = a % txz * b % txz
c % tyz = a % tyz * b % tyz
c % fx = a % fx * b % fx
c % fy = a % fy * b % fy
c % fz = a % fz * b % fz
c % cs_opt2 = a % cs_opt2 * b % cs_opt2
c % Nu_t = a % Nu_t * b % Nu_t !Eshwan
c % divtx = a % divtx * b % divtx
c % divty = a % divty * b % divty
c % divtz = a % divtz * b % divtz
c % buoyancy = a % buoyancy * b % buoyancy

return
end function tavg_mul

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_mul( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_t) :: c

c % u = a % u * b
c % v = a % v * b
c % w = a % w * b
c % u_uv = a % u_uv * b
c % v_uv = a % v_uv * b
c % w_uv = a % w_uv * b
c % u2 = a % u2 * b
c % v2 = a % v2 * b
c % w2 = a % w2 * b
c % uv = a % uv * b
c % uw = a % uw * b
c % vw = a % vw * b
c % dudz = a % dudz * b
c % dvdz = a % dvdz * b
c % dudz2 = a % dudz2 * b
c % txx = a % txx * b
c % tyy = a % tyy * b
c % tzz = a % tzz * b
c % txy = a % txy * b
c % txz = a % txz * b
c % tyz = a % tyz * b
c % fx = a % fx * b
c % fy = a % fy * b
c % fz = a % fz * b
c % cs_opt2 = a % cs_opt2 * b
c % Nu_t = a % Nu_t * b !Eshwan
c % divtx = a % divtx * b 
c % divty = a % divty * b
c % divtz = a % divtz * b 
c % buoyancy = a % buoyancy * b 

return
end function tavg_scalar_mul

!$if($OUTPUT_EXTRA) !- Eshwan
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_sgs_add( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_sgs_t), intent(in) :: a, b
type(tavg_sgs_t) :: c

c % Tn = a % Tn + b % Tn !Eshwan
c % Beta = a % Beta + b % Beta !Eshwan
!c % Nu_t = a % Nu_t + b % Nu_t !Eshwan
c % F_LM = a % F_LM + b % F_LM
c % F_MM = a % F_MM + b % F_LM
c % F_QN = a % F_QN + b % F_QN
c % F_NN = a % F_NN + b % F_NN
c % ee_now = a % ee_now + b % ee_now
c % cs2_clips = a % cs2_clips + b % cs2_clips
$if($DYN_TN)
c % F_ee2 = a % F_ee2 + b % F_ee2
c % F_deedt2 = a % F_deedt2 + b % F_ee2
$endif

return
end function tavg_sgs_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_sgs_add( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_scalar_sgs_t), intent(in) :: a, b
type(tavg_scalar_sgs_t) :: c

c % kappa_tt = a % kappa_tt + b % kappa_tt
c % Ds_opt2_t = a % Ds_opt2_t + b % Ds_opt2_t
c % s_Tn_t = a % s_Tn_t + b % s_Tn_t
c % I_LM_t = a % I_LM_t + b % I_LM_t
c % I_MM_t = a % I_MM_t + b % I_MM_t
c % I_QN_t = a % I_QN_t + b % I_QN_t
c % I_NN_t = a % I_NN_t + b % I_NN_t
c % ds2_clips_t = a % ds2_clips_t + b % ds2_clips_t
c % s_Beta_t = a % s_Beta_t + b % s_Beta_t

c % kappa_ts = a % kappa_ts + b % kappa_ts
c % Ds_opt2_s = a % Ds_opt2_s + b % Ds_opt2_s
c % s_Tn_s = a % s_Tn_s + b % s_Tn_s
c % I_LM_s = a % I_LM_s + b % I_LM_s
c % I_MM_s = a % I_MM_s + b % I_MM_s
c % I_QN_s = a % I_QN_s + b % I_QN_s
c % I_NN_s = a % I_NN_s + b % I_NN_s
c % ds2_clips_s = a % ds2_clips_s + b % ds2_clips_s
c % s_Beta_s = a % s_Beta_s + b % s_Beta_s

return
end function tavg_scalar_sgs_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_sgs_scalar_div( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_sgs_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_sgs_t) :: c

c % Tn = a % Tn / b !Eshwan
c % Beta = a % Beta / b  !Eshwan
!c % Nu_t = a % Nu_t / b !Eshwan
c % F_LM = a % F_LM / b
c % F_MM = a % F_MM / b
c % F_QN = a % F_QN / b
c % F_NN = a % F_NN / b
c % ee_now = a % ee_now / b
c % cs2_clips = a % cs2_clips
$if($DYN_TN)
c % F_ee2 = a % F_ee2 / b
c % F_deedt2 = a % F_deedt2 / b
$endif

return
end function tavg_sgs_scalar_div
!$endif !- Eshwan

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_sgs_scalar_div( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_scalar_sgs_t), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_scalar_sgs_t) :: c

c % kappa_tt = a % kappa_tt / b
c % Ds_opt2_t = a % Ds_opt2_t / b 
c % s_Tn_t = a % s_Tn_t / b
c % I_LM_t = a % I_LM_t / b
c % I_MM_t = a % I_MM_t / b
c % I_QN_t = a % I_QN_t / b
c % I_NN_t = a % I_NN_t / b
c % ds2_clips_t = a % ds2_clips_t 
c % s_Beta_t = a % s_Beta_t

c % kappa_ts = a % kappa_ts / b 
c % Ds_opt2_s = a % Ds_opt2_s / b 
c % s_Tn_s = a % s_Tn_s / b
c % I_LM_s = a % I_LM_s / b
c % I_MM_s = a % I_MM_s / b
c % I_QN_s = a % I_QN_s / b
c % I_NN_s = a % I_NN_s / b
c % ds2_clips_s = a % ds2_clips_s 
c % s_Beta_s = a % s_Beta_s

return
end function tavg_scalar_sgs_scalar_div

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_interp_to_uv_grid( a ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

!c % fz = interp_to_uv_grid( a % fz, lbz )

return

end function tavg_interp_to_uv_grid

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_interp_to_w_grid( a ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

c % fx = interp_to_w_grid( a % fx, lbz )
c % fy = interp_to_w_grid( a % fy, lbz )

return

end function tavg_interp_to_w_grid


!//////////////////////////////////////////////////////////////////////
!/////////////////// THETA_AVG OPERATORS //////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function theta_avg_scalar_div( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(theta_avg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(theta_avg_t) :: c

c % theta  = a % theta / b
c % theta_w  = a % theta_w / b
c % wtheta = a % wtheta / b
c % theta2 = a % theta2 / b
c % dTdx = a % dTdx / b
c % dTdy = a % dTdy / b
c % dTdz = a % dTdz / b
c % sgs_tflux = a % sgs_tflux / b

return
end function theta_avg_scalar_div


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function theta_avg_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(theta_avg_t), intent(in) :: a, b
type(theta_avg_t) :: c

c % theta = a % theta + b % theta
c % theta_w = a % theta_w + b % theta_w
c % wtheta = a % wtheta + b % wtheta
c % theta2 = a % theta2 + b % theta2
c % dTdx = a % dTdx + b % dTdx
c % dTdy = a % dTdy + b % dTdy
c % dTdz = a % dTdz + b % dTdz
c % sgs_tflux = a % sgs_tflux + b % sgs_tflux

return
end function theta_avg_add


!//////////////////////////////////////////////////////////////////////
!//////////////////// SAL_AVG OPERATORS ///////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function sal_avg_scalar_div( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(sal_avg_t), intent(in) :: a
real(rprec), intent(in) :: b
type(sal_avg_t) :: c

c % sal  = a % sal / b
c % sal_w  = a % sal_w / b
c % wsal = a % wsal / b
c % sal2 = a % sal2 / b
c % dSdx = a % dSdx / b 
c % dSdy = a % dSdy / b 
c % dSdz = a % dSdz / b 
c % sgs_sflux = a % sgs_sflux / b

return
end function sal_avg_scalar_div


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function sal_avg_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(sal_avg_t), intent(in) :: a, b
type(sal_avg_t) :: c

c % sal = a % sal + b % sal
c % sal_w = a % sal_w + b % sal_w
c % wsal = a % wsal + b % wsal
c % sal2 = a % sal2 + b % sal2
c % dSdx = a % dSdx + b % dSdx
c % dSdy = a % dSdy + b % dSdy
c % dSdz = a % dSdz + b % dSdz
c % sgs_sflux = a % sgs_sflux + b % sgs_sflux

return
end function sal_avg_add


!//////////////////////////////////////////////////////////////////////
!///////////////////// RS OPERATORS ///////////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(rs_t), intent(in) :: a, b
type(rs_t) :: c

c % up2 = a % up2 + b % up2
c % vp2 = a % vp2 + b % vp2
c % wp2 = a % wp2 + b % wp2
c % upvp = a % upvp + b % upvp
c % upwp = a % upwp + b % upwp
c % vpwp = a % vpwp + b % vpwp

return
end function rs_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_sub( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(rs_t), intent(in) :: a, b
type(rs_t) :: c

c % up2 = a % up2 - b % up2
c % vp2 = a % vp2 - b % vp2
c % wp2 = a % wp2 - b % wp2
c % upvp = a % upvp - b % upvp
c % upwp = a % upwp - b % upwp
c % vpwp = a % vpwp - b % vpwp

return
end function rs_sub

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_scalar_div( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

return
end function rs_scalar_div


!//////////////////////////////////////////////////////////////////////
!///////////////// THETA STAT OPERATORS ///////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function theta_stat_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(theta_stat_t), intent(in) :: a, b
type(theta_stat_t) :: c

c % wpthetap = a % wpthetap + b % wpthetap
c % thetap2  = a % thetap2 + b % thetap2
c % total_tflux = a % total_tflux + b % total_tflux

return
end function theta_stat_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function theta_stat_scalar_div( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(theta_stat_t), intent(in) :: a
real(rprec), intent(in) :: b
type(theta_stat_t) :: c

c % wpthetap = a % wpthetap / b
c % thetap2 = a % thetap2 / b 
c % total_tflux = a % total_tflux / b

return
end function theta_stat_scalar_div


!//////////////////////////////////////////////////////////////////////
!////////////////// SAL STAT OPERATORS ////////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function sal_stat_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(sal_stat_t), intent(in) :: a, b
type(sal_stat_t) :: c

c % wpsalp = a % wpsalp + b % wpsalp
c % salp2  = a % salp2 + b % salp2
c % total_sflux = a % total_sflux + b % total_sflux

return
end function sal_stat_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function sal_stat_scalar_div( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(sal_stat_t), intent(in) :: a
real(rprec), intent(in) :: b
type(sal_stat_t) :: c

c % wpsalp = a % wpsalp / b
c % salp2 = a % salp2 / b 
c % total_sflux = a % total_sflux / b

return
end function sal_stat_scalar_div


!//////////////////////////////////////////////////////////////////////
!/////////////// SIGMA AVG ZPLANE FUNCTIONS ///////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function sigma_avg_zplane_scalar_div( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(sigma_t), intent(in) :: a
real(rprec), intent(in) :: b
type(sigma_t) :: c

c % sigma_theta = a % sigma_theta / b
c % sigma_sal = a % sigma_sal / b

return
end function sigma_avg_zplane_scalar_div


!//////////////////////////////////////////////////////////////////////
!/////////////////// SPECIAL RS FUNCTIONS /////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_compute( a , lbz2) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
c % upwp = a % uw - a % u * a % w
c % vpwp = a % vw - a % v * a % w

return
end function rs_compute


!//////////////////////////////////////////////////////////////////////
!///////////////// SPECIAL THETA STAT FUNCTIONS ///////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function theta_stat_compute( a , s , lbz2) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(theta_avg_t), dimension(:,:,lbz2:), intent(in) :: s
type(theta_stat_t), allocatable, dimension(:,:,:) :: c
integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % wpthetap = s % wtheta - a % w * s % theta_w
c % thetap2 = s % theta2 - s % theta_w * s % theta_w
c % total_tflux = (s % wtheta - a %w * s % theta_w) - s % sgs_tflux

return
end function theta_stat_compute


!//////////////////////////////////////////////////////////////////////
!////////////////// SPECIAL SAL STAT FUNCTIONS ////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function sal_stat_compute( a , s , lbz2) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(sal_avg_t), dimension(:,:,lbz2:), intent(in) :: s
type(sal_stat_t), allocatable, dimension(:,:,:) :: c
integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % wpsalp = s % wsal - a % w * s % sal_w
c % salp2 = s % sal2 - s % sal_w * s % sal_w
c % total_sflux = (s % wsal - a % w * s % sal_w) - s % sgs_sflux

return
end function sal_stat_compute


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function cnpy_tavg_mul( a ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

return
end function cnpy_tavg_mul

!//////////////////////////////////////////////////////////////////////
!///////////////// SPECIAL TAVG SUBROUTINES ///////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_t), intent(out) :: c

c % u = a
c % v = a
c % w = a
c % u_uv = a
c % v_uv = a
c % w_uv = a
c % u2 = a
c % v2 = a
c % w2 = a
c % uv = a
c % uw = a
c % vw = a
c % dudz = a
c % dvdz = a
c % dudz2 = a
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
c % Nu_t = a !Eshwan
c % divtx = a
c % divty = a
c % divtz = a
c % buoyancy = a

return
end subroutine tavg_set

!$if($OUTPUT_EXTRA) !- Eshwan
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_sgs_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_sgs_t), intent(out) :: c

c % Tn =  a !Eshwan
c % Beta =  a !Eshwan
!c % Nu_t =  a !Eshwan
c % F_LM =  a
c % F_MM =  a
c % F_QN =  a
c % F_NN =  a
c % ee_now = a
c % cs2_clips = a
$if($DYN_TN)
c % F_ee2 = a
c % F_deedt2 = a
$endif

return
end subroutine tavg_sgs_set
!$endif !- Eshwan

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_scalar_sgs_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_scalar_sgs_t), intent(out) :: c

c % kappa_tt = a 
c % Ds_opt2_t = a 
c % s_Tn_t = a
c % I_LM_t = a
c % I_MM_t = a
c % I_QN_t = a
c % I_NN_t = a
c % ds2_clips_t = a
c % s_Beta_t = a

c % kappa_ts = a 
c % Ds_opt2_s = a 
c % s_Tn_s = a
c % I_LM_s = a
c % I_MM_s = a
c % I_QN_s = a
c % I_NN_s = a
c % ds2_clips_s = a
c % s_Beta_s = a

return
end subroutine tavg_scalar_sgs_set
!$endif !- Eshwan


!//////////////////////////////////////////////////////////////////////
!///////////////// SPECIAL THETA_AVG SUBROUTINES //////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine theta_avg_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(theta_avg_t), intent(out) :: c

c % theta = a
c % theta_w = a
c % wtheta = a
c % theta2 = a
c % dTdx = a
c % dTdy = a
c % dTdz = a
c % sgs_tflux = a

return
end subroutine theta_avg_set


!//////////////////////////////////////////////////////////////////////
!///////////////// SPECIAL SAL_AVG SUBROUTINES //////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sal_avg_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(sal_avg_t), intent(out) :: c

c % sal = a
c % sal_w = a
c % wsal = a
c % sal2 = a
c % dSdx = a
c % dSdy = a
c % dSdz = a
c % sgs_sflux = a

return
end subroutine sal_avg_set


!//////////////////////////////////////////////////////////////////////
!/////////////////// SPECIAL RS SUBROUTINES ///////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rs_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

return
end subroutine rs_set


!//////////////////////////////////////////////////////////////////////
!///////////// SPECIAL THETA STAT SUBROUTINES /////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine theta_stat_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(theta_stat_t), intent(out) :: c

c % wpthetap = a
c % thetap2 = a
c % total_tflux = a

return
end subroutine theta_stat_set



!//////////////////////////////////////////////////////////////////////
!////////////// SPECIAL SAL STAT SUBROUTINES //////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sal_stat_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(sal_stat_t), intent(out) :: c

c % wpsalp = a
c % salp2 = a
c % total_sflux = a

return
end subroutine sal_stat_set



!//////////////////////////////////////////////////////////////////////
!/////////////// SPECIAL SIGMA SUBROUTINES ////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sigma_avg_zplane_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(sigma_t), intent(out) :: c

c % sigma_theta = a 
c % sigma_sal =  a

return
end subroutine sigma_avg_zplane_set

!//////////////////////////////////////////////////////////////////////
!/////////////////// SPECIAL HIST SUBROUTINES /////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine hist_binit_1D( a, var, phi_ls )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine takes the values in var and bins them in the histogram
!   a only if the location is outside a body (level set function phi>0 )
!
! If phi_ls is not included as an input then all values are binned
!
! Inputs:
!   a   :: the histogram to be updated
!   var :: the values that will be used to update the histogram
!   phi_ls :: the level set function (phi from level_set_base module)

use types, only : rprec
implicit none

type(hist_t), intent(inout) :: a                
real(rprec), intent(in), dimension(:) :: var
real(rprec), intent(in), dimension(:), optional :: phi_ls ! phi_ls<0 is inside a body

integer :: dim1, i, ib
real(rprec) :: countme

! Determine length of input arrays
    dim1 = size(var,1)

    !! Check that phi_ls is the same length as var (if present)
    !if (present (phi_ls)) then
    !    if ( size(phi_ls,1) .ne. dim1 ) then
    !        write(*,*) 'In hist_binit_1D: size of phi_ls should match size of var'
    !        stop
    !    endif
    !endif

! Prepare temp array and counting variable
    countme = 1.0_rprec
    
$if ($LVLSET) 
    if (present (phi_ls)) then
        do i=1,dim1
            ! if phi<0 (inside body) don't count it!  (1=outside body, 0=inside)
            countme = 0.5_rprec * ( 1.0_rprec + sign(1.0_rprec,phi_ls(i)) )  

            ! Determine which bin and add 1.0 to that val
            ib = min( ceiling( max(var(i)-a%bmin,-0.5_rprec) /a%db ), a%nbins+1 )
            a%vals(ib) = a%vals(ib) + countme
        enddo
    else
$endif
        do i=1,dim1
            ! Determine which bin and add 1.0 to that val
            ib = min( ceiling( max(var(i)-a%bmin,-0.5_rprec) /a%db ), a%nbins+1 )
            a%vals(ib) = a%vals(ib) + countme
        enddo
$if ($LVLSET) 
    endif
$endif

return
end subroutine hist_binit_1D

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine hist_binit_2D( a, var, phi_ls )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine takes the values in var and bins them in the histogram
!   a only if the location is outside a body (level set function phi>0 )
!
! If phi_ls is not included as an input then all values are binned
!
! Inputs:
!   a   :: the histogram to be updated
!   var :: the values that will be used to update the histogram
!   phi_ls :: the level set function (phi from level_set_base module)

use types, only : rprec
implicit none

type(hist_t), intent(inout) :: a                
real(rprec), intent(in), dimension(:,:) :: var
real(rprec), intent(in), dimension(:,:), optional :: phi_ls 

integer :: dim1, dim2, i, j, ib
real(rprec) :: countme

! Determine length of input arrays
    dim1 = size(var,1)
    dim2 = size(var,2)

    !! Check that phi_ls is the same length as var (if present)
    !if (present (phi_ls)) then
    !    if (( size(phi_ls,1) .ne. dim1 ).or.( size(phi_ls,2) .ne. dim2 )) then
    !        write(*,*) 'In hist_binit_2D: size of phi_ls should match size of var'
    !        stop
    !    endif
    !endif

! Prepare temp array and counting variable
    countme = 1.0_rprec
    
$if ($LVLSET) 
    if (present (phi_ls)) then
        do j=1,dim2
        do i=1,dim1
            ! if phi<0 (inside body) don't count it!  (1=outside body, 0=inside)
            countme = 0.5_rprec * ( 1.0_rprec + sign(1.0_rprec,phi_ls(i,j)) )  

            ! Determine which bin and add 1.0 to that val
            ib = min( ceiling( max(var(i,j)-a%bmin,-0.5_rprec) /a%db ), a%nbins+1 )
            a%vals(ib) = a%vals(ib) + countme
        enddo
        enddo
    else
$endif
        do j=1,dim2
        do i=1,dim1
            ! Determine which bin and add 1.0 to that val
            ib = min( ceiling( max(var(i,j)-a%bmin,-0.5_rprec) /a%db ), a%nbins+1 )
            a%vals(ib) = a%vals(ib) + countme
        enddo
        enddo
$if ($LVLSET) 
    endif
$endif

return
end subroutine hist_binit_2D

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine hist_binit_3D( a, var, phi_ls )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine takes the values in var and bins them in the histogram
!   a only if the location is outside a body (level set function phi>0 )
!
! If phi_ls is not included as an input then all values are binned
!
! Inputs:
!   a   :: the histogram to be updated
!   var :: the values that will be used to update the histogram
!   phi_ls :: the level set function (phi from level_set_base module)

use types, only : rprec
implicit none

type(hist_t), intent(inout) :: a                
real(rprec), intent(in), dimension(:,:,:) :: var
real(rprec), intent(in), dimension(:,:,:), optional :: phi_ls 

integer :: dim1, dim2, dim3, i, j, k, ib
real(rprec) :: countme

! Determine length of input arrays
    dim1 = size(var,1)
    dim2 = size(var,2)
    dim3 = size(var,3)

    !! Check that phi_ls is the same length as var (if present)
    !if (present (phi_ls)) then
    !    if ( size(phi_ls,1) .ne. dim1 ) then
    !        write(*,*) 'In hist_binit_3D: size of phi_ls should match size of var (1)'
    !        stop
    !    endif
    !    if ( size(phi_ls,2) .ne. dim2 ) then
    !        write(*,*) 'In hist_binit_3D: size of phi_ls should match size of var (2)'
    !        stop
    !    endif
    !    if ( size(phi_ls,3) .ne. dim3 ) then
    !        write(*,*) 'In hist_binit_3D: size of phi_ls should match size of var (3)'
    !        stop
    !   endif
    !endif

! Prepare temp array and counting variable
    countme = 1.0_rprec
    
$if ($LVLSET) 
    if (present (phi_ls)) then
        do k=1,dim3
        do j=1,dim2
        do i=1,dim1
            ! if phi<0 (inside body) don't count it!  (1=outside body, 0=inside)
            countme = 0.5_rprec * ( 1.0_rprec + sign(1.0_rprec,phi_ls(i,j,k)) )  

            ! Determine which bin and add 1.0 to that val
            ib = min( ceiling( max(var(i,j,k)-a%bmin,-0.5_rprec) /a%db ), a%nbins+1 )
            a%vals(ib) = a%vals(ib) + countme
        enddo
        enddo
        enddo
    else
$endif
        do k=1,dim3
        do j=1,dim2
        do i=1,dim1
            ! Determine which bin and add 1.0 to that val
            ib = min( ceiling( max(var(i,j,k)-a%bmin,-0.5_rprec) /a%db ), a%nbins+1 )
            a%vals(ib) = a%vals(ib) + countme
        enddo
        enddo
        enddo
$if ($LVLSET) 
    endif
$endif

return
end subroutine hist_binit_3D

end module stat_defs

