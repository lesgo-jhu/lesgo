!**********************************************************************
module stat_defs
!**********************************************************************
use types, only : rprec
use param, only : nx,ny,nz,lh

save
public

type point
  integer :: istart, jstart, kstart, coord
  real(rprec) :: xdiff, ydiff, zdiff
  character(64) :: fname
end type point

type plane
  integer :: istart
  real(rprec) :: ldiff
end type plane

type zplane
  integer :: istart, coord
  real(rprec) :: ldiff
end type zplane  

type rs
  real(rprec) :: up2, vp2, wp2, upwp, vpwp, upvp
end type rs

type spectra
  real(rprec), dimension(lh) :: power
  integer :: istart, coord
  real(rprec) :: ldiff 
end type spectra

real(rprec) :: spectra_total_time
real(rprec) :: tavg_total_time

!  Sums performed over time
type tavg
  real(rprec) :: u, v, w, u2, v2, w2, uw, vw, uv
  real(rprec) :: dudz, dvdz
  real(rprec) :: txx, txy, tyy, txz, tyz, tzz
  real(rprec) :: fx, fy, fz
  real(rprec) :: cs_opt2
end type tavg

$if ($TURBINES)
type turbine 
  real(rprec) :: xloc, yloc, height, dia, thk
  real(rprec) :: vol_c                        !term used for volume correction  
  real(rprec) :: theta1                       !angle CCW(from above) from -x direction [degrees]
  real(rprec) :: theta2                    !angle above the horizontal, from -x dir [degrees]
  real(rprec), dimension(3) :: nhat           !(nx,ny,nz) of unit normal for each turbine
  integer :: num_nodes                        !number of nodes associated with each turbine
  integer, dimension(1500,3) :: nodes         !(i,j,k) of each included node
  integer, dimension(6) :: nodes_max          !search area for nearby nodes
  real(rprec) :: u_d, u_d_T                   !running time-average of mean disk velocity
  real(rprec) :: f_n                          !normal force on turbine disk
  real(rprec), dimension(1500) :: ind         !indicator function - weighting of each node
        
  !real(rprec), dimension (nx,ny,nz) :: u_cond_avg_lo, v_cond_avg_lo, w_cond_avg_lo
  !real(rprec), dimension (nx,ny,nz) :: u_cond_avg_hi, v_cond_avg_hi, w_cond_avg_hi
  !logical :: cond_avg_calc_lo,cond_avg_calc_hi
  !real(rprec) :: cond_avg_ud_lo,cond_avg_ud_hi,cond_avg_time_lo,cond_avg_time_hi
end type turbine

type wind_farm
  integer ::ifilter                           !Filter type: 2->Gaussian
  real(rprec) :: alpha                        !filter size is alpha*(grid spacing)
  integer :: trunc                            !truncated - # grid points to include in each dir.
  real(rprec) :: filter_cutoff                !ind only includes values above this cutoff
  type(turbine), pointer, dimension(:) :: turbine_t
        
  !logical, pointer, dimension(:) :: cond_avg_flag_lo,cond_avg_flag_hi
end type wind_farm
    
type(wind_farm)        :: wind_farm_t
$endif

type(point), allocatable, dimension(:) :: point_t
type(plane), allocatable, dimension(:) :: xplane_t, yplane_t
type(zplane), allocatable, dimension(:) :: zplane_t

type(tavg), allocatable, dimension(:,:,:) :: tavg_t
type(tavg), allocatable, dimension(:) :: tavg_zplane_t

type(rs), allocatable, dimension(:,:,:) :: rs_t
type(rs), allocatable, dimension(:) :: rs_zplane_t, cnpy_zplane_t
type(spectra), allocatable, dimension(:) :: spectra_t

! Overloaded operators for tavg and rs types
INTERFACE OPERATOR (.ADD.)
  MODULE PROCEDURE tavg_add, tavg_scalar_add, rs_add
END INTERFACE

INTERFACE OPERATOR (.SUB.)
  MODULE PROCEDURE tavg_sub, rs_sub
END INTERFACE

INTERFACE OPERATOR (.DIV.)
  MODULE PROCEDURE tavg_scalar_div, rs_scalar_div
END INTERFACE

INTERFACE OPERATOR (.MUL.)
  MODULE PROCEDURE tavg_mul, tavg_scalar_mul
END INTERFACE

INTERFACE type_set
  MODULE PROCEDURE tavg_set, rs_set
END INTERFACE

INTERFACE type_zero_bogus
  MODULE PROCEDURE tavg_zero_bogus_2D, tavg_zero_bogus_3D
END INTERFACE

contains

!//////////////////////////////////////////////////////////////////////
!/////////////////// TAVG OPERATORS ///////////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(tavg), intent(in) :: a, b
type(tavg) :: c

c % u = a % u + b % u
c % v = a % v + b % v
c % w = a % w + b % w
c % u2 = a % u2 + b % u2
c % v2 = a % v2 + b % v2
c % w2 = a % w2 + b % w2
c % uw = a % uw + b % uw
c % vw = a % vw + b % vw
c % uv = a % uv + b % uv
c % dudz = a % dudz + b % dudz
c % dvdz = a % dvdz + b % dvdz
c % txx = a % txx + b % txx
c % txy = a % txy + b % txy
c % tyy = a % tyy + b % tyy
c % txz = a % txz + b % txz
c % tyz = a % tyz + b % tyz
c % tzz = a % tzz + b % tzz
c % fx = a % fx + b % fx
c % fy = a % fy + b % fy
c % fz = a % fz + b % fz
c % cs_opt2 = a % cs_opt2 + b % cs_opt2

return
end function tavg_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_sub( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(tavg), intent(in) :: a, b
type(tavg) :: c

c % u = a % u - b % u
c % v = a % v - b % v
c % w = a % w - b % w
c % u2 = a % u2 - b % u2
c % v2 = a % v2 - b % v2
c % w2 = a % w2 - b % w2
c % uw = a % uw - b % uw
c % vw = a % vw - b % vw
c % uv = a % uv - b % uv
c % dudz = a % dudz - b % dudz
c % dvdz = a % dvdz - b % dvdz
c % txx = a % txx - b % txx
c % txy = a % txy - b % txy
c % tyy = a % tyy - b % tyy
c % txz = a % txz - b % txz
c % tyz = a % tyz - b % tyz
c % tzz = a % tzz - b % tzz
c % fx = a % fx - b % fx
c % fy = a % fy - b % fy
c % fz = a % fz - b % fz
c % cs_opt2 = a % cs_opt2 - b % cs_opt2

return
end function tavg_sub

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_add( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg) :: c

c % u = a % u + b
c % v = a % v + b
c % w = a % w + b
c % u2 = a % u2 + b
c % v2 = a % v2 + b
c % w2 = a % w2 + b
c % uw = a % uw + b
c % vw = a % vw + b
c % uv = a % uv + b
c % dudz = a % dudz + b
c % dvdz = a % dvdz + b
c % txx = a % txx + b
c % txy = a % txy + b
c % tyy = a % tyy + b
c % txz = a % txz + b
c % tyz = a % tyz + b
c % tzz = a % tzz + b
c % fx = a % fx + b
c % fy = a % fy + b
c % fz = a % fz + b
c % cs_opt2 = a % cs_opt2 + b

return
end function tavg_scalar_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_zero_bogus_2D( c )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg), dimension(:,:), intent(inout) :: c

c % txx = 0._rprec
c % txy = 0._rprec
c % tyy = 0._rprec
c % txz = 0._rprec
c % tyz = 0._rprec
c % tzz = 0._rprec
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

type(tavg), dimension(:,:,:), intent(inout) :: c

c % txx = 0._rprec
c % txy = 0._rprec
c % tyy = 0._rprec
c % txz = 0._rprec
c % tyz = 0._rprec
c % tzz = 0._rprec
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

type(tavg), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg) :: c

c % u = a % u / b
c % v = a % v / b
c % w = a % w / b
c % u2 = a % u2 / b
c % v2 = a % v2 / b
c % w2 = a % w2 / b
c % uw = a % uw / b
c % vw = a % vw / b
c % uv = a % uv / b
c % dudz = a % dudz / b
c % dvdz = a % dvdz / b
c % txx = a % txx / b
c % txy = a % txy / b
c % tyy = a % tyy / b
c % txz = a % txz / b
c % tyz = a % tyz / b
c % tzz = a % tzz / b
c % fx = a % fx / b
c % fy = a % fy / b
c % fz = a % fz / b
c % cs_opt2 = a % cs_opt2 / b

return
end function tavg_scalar_div

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_mul( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
type(tavg), intent(in) :: a, b
type(tavg) :: c

c % u = a % u * b % u
c % v = a % v * b % v
c % w = a % w * b % w
c % u2 = a % u2 * b % u2
c % v2 = a % v2 * b % v2
c % w2 = a % w2 * b % w2
c % uw = a % uw * b % uw
c % vw = a % vw * b % vw
c % uv = a % uv * b % uv
c % dudz = a % dudz * b % dudz
c % dvdz = a % dvdz * b % dvdz
c % txx = a % txx * b % txx
c % txy = a % txy * b % txy
c % tyy = a % tyy * b % tyy
c % txz = a % txz * b % txz
c % tyz = a % tyz * b % tyz
c % tzz = a % tzz * b % tzz
c % fx = a % fx * b % fx
c % fy = a % fy * b % fy
c % fz = a % fz * b % fz
c % cs_opt2 = a % cs_opt2 * b % cs_opt2

return
end function tavg_mul

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_scalar_mul( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg) :: c

c % u = a % u * b
c % v = a % v * b
c % w = a % w * b
c % u2 = a % u2 * b
c % v2 = a % v2 * b
c % w2 = a % w2 * b
c % uw = a % uw * b
c % vw = a % vw * b
c % uv = a % uv * b
c % dudz = a % dudz * b
c % dvdz = a % dvdz * b
c % txx = a % txx * b
c % txy = a % txy * b
c % tyy = a % tyy * b
c % txz = a % txz * b
c % tyz = a % tyz * b
c % tzz = a % tzz * b
c % fx = a % fx * b
c % fy = a % fy * b
c % fz = a % fz * b
c % cs_opt2 = a % cs_opt2 * b

return
end function tavg_scalar_mul

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_interp_to_uv_grid( a ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only: lbz
use functions, only : interp_to_uv_grid
implicit none

type(tavg), dimension(:,:,lbz:), intent(in) :: a
type(tavg), allocatable, dimension(:,:,:) :: c

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

type(tavg), dimension(:,:,lbz:), intent(in) :: a
type(tavg), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx = ubound(a,1)
uby = ubound(a,2)
ubz = ubound(a,3)

allocate(c(ubx,uby,lbz:ubz))

c = a

c % txx =  interp_to_w_grid( a % txx, lbz )
c % txy =  interp_to_w_grid( a % txy, lbz )
c % tyy =  interp_to_w_grid( a % tyy, lbz )
c % tzz =  interp_to_w_grid( a % tzz, lbz )

c % fx = interp_to_w_grid( a % fx, lbz )
c % fy = interp_to_w_grid( a % fy, lbz )

return

end function tavg_interp_to_w_grid


!//////////////////////////////////////////////////////////////////////
!///////////////////// RS OPERATORS ///////////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_add( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(rs), intent(in) :: a, b
type(rs) :: c

c % up2 = a % up2 + b % up2
c % vp2 = a % vp2 + b % vp2
c % wp2 = a % wp2 + b % wp2
c % upwp = a % upwp + b % upwp
c % vpwp = a % vpwp + b % vpwp
c % upvp = a % upvp + b % upvp

end function rs_add

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_sub( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(rs), intent(in) :: a, b
type(rs) :: c

c % up2 = a % up2 - b % up2
c % vp2 = a % vp2 - b % vp2
c % wp2 = a % wp2 - b % wp2
c % upwp = a % upwp - b % upwp
c % vpwp = a % vpwp - b % vpwp
c % upvp = a % upvp - b % upvp

end function rs_sub

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_scalar_div( a, b) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(rs), intent(in) :: a
real(rprec), intent(in) :: b
type(rs) :: c

c % up2 = a % up2 / b
c % vp2 = a % vp2 / b 
c % wp2 = a % wp2 / b
c % upwp = a % upwp / b 
c % vpwp = a % vpwp / b 
c % upvp = a % upvp / b

end function rs_scalar_div

!//////////////////////////////////////////////////////////////////////
!/////////////////// SPECIAL RS FUNCTIONS /////////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function rs_compute( a , lbz2) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer, intent(in) :: lbz2
type(tavg), dimension(:,:,lbz2:), intent(in) :: a
type(rs), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % up2 = a % u2 - a % u * a % u
c % vp2 = a % v2 - a % v * a % v
c % wp2 = a % w2 - a % w * a % w
c % upwp = a % uw - a % u * a % w
c % vpwp = a % vw - a % v * a % w
c % upvp = a % uv - a % u * a % v

return

end function rs_compute

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function cnpy_tavg_mul( a ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! This performs one set of multiplication for the canopy stresses
!
implicit none

type(tavg), intent(in) :: a
type(rs) :: c

c % up2 = a % u * a % u
c % vp2 = a % v * a % v
c % wp2 = a % w * a % w
c % upwp = a % u * a % w
c % vpwp = a % v * a % w
c % upvp = a % u * a % v

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
type(tavg), intent(out) :: c

c % u = a
c % v = a
c % w = a
c % u2 = a
c % v2 = a
c % w2 = a
c % uw = a
c % vw = a
c % uv = a
c % dudz = a
c % dvdz = a
c % txx = a
c % txy = a
c % tyy = a
c % txz = a
c % tyz = a
c % tzz = a
c % fx = a
c % fy = a
c % fz = a
c % cs_opt2 = a

return
end subroutine tavg_set

!//////////////////////////////////////////////////////////////////////
!/////////////////// SPECIAL RS SUBROUTINES ///////////////////////////
!//////////////////////////////////////////////////////////////////////

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rs_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(rs), intent(out) :: c

c % up2 = a
c % vp2 = a
c % wp2 = a
c % upwp = a
c % vpwp = a
c % upvp = a

return
end subroutine rs_set

end module stat_defs

