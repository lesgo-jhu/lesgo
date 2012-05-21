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
  integer :: fid
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
  real(rprec), dimension(:), allocatable :: power
  integer :: istart, coord
  real(rprec) :: ldiff 
end type spectra

real(rprec) :: spectra_total_time
real(rprec) :: tavg_total_time
$if($OUTPUT_EXTRA)
real(rprec) :: tavg_total_time_sgs
$endif
! Time step used for time averaging computations
real(rprec) :: dt_tavg
! Time stamp of last time averaging computation
real(rprec) :: tavg_time_stamp  
  
!  Sums performed over time
type tavg
  real(rprec) :: u, v, w, u2, v2, w2, uw, vw, uv
  real(rprec) :: dudz, dvdz
  real(rprec) :: txx, txy, tyy, txz, tyz, tzz
  real(rprec) :: fx, fy, fz
  real(rprec) :: cs_opt2  
end type tavg

!  Sums performed over time (for subgrid variables)
$if($OUTPUT_EXTRA)
type tavg_sgs
  real(rprec) :: Tn, Nu_t
  real(rprec) :: F_LM, F_MM, F_QN, F_NN
  real(rprec) :: ee_now
  $if ($DYN_TN)
  real(rprec) :: F_ee2, F_deedt2
  $endif
end type tavg_sgs
$endif

! Types for including wind-turbines as drag disks
$if ($TURBINES)
! Single turbines
type turbine 
  real(rprec) :: xloc, yloc, height, dia, thk
  real(rprec) :: vol_c                        ! term used for volume correction  
  real(rprec) :: theta1                       ! angle CCW(from above) from -x direction [degrees]
  real(rprec) :: theta2                       ! angle above the horizontal, from -x dir [degrees]
  real(rprec), dimension(3) :: nhat           ! (nx,ny,nz) of unit normal for each turbine
  integer :: num_nodes                        ! number of nodes associated with each turbine
  integer, dimension(1500,3) :: nodes         ! (i,j,k) of each included node
  integer, dimension(6) :: nodes_max          ! search area for nearby nodes
  real(rprec) :: u_d, u_d_T                   ! running time-average of mean disk velocity
  real(rprec) :: f_n                          ! normal force on turbine disk
  real(rprec), dimension(1500) :: ind         ! indicator function - weighting of each node
end type turbine

! A collection of wind-turbines
type wind_farm
  type(turbine), pointer, dimension(:) :: turbine_t
end type wind_farm
    
type(wind_farm) :: wind_farm_t
$endif

! Histogram (single)
type hist
    real(rprec) :: bmin, bmax, db             ! bin min, max, and spacing
    integer :: nbins                          ! number of bins
    real(rprec), allocatable, dimension(:) :: bins  ! bin centers
    real(rprec), allocatable, dimension(:) :: vals  ! count for each bin (may be normalized)
end type hist

! Collection of histograms (one for each zplane) for a single variable
type hist_zplanes  
    integer, allocatable, dimension(:) :: coord         ! processor where this plane exists
    integer, allocatable, dimension(:) :: istart        ! nearest node below plane (for interpolation)
    real(rprec), allocatable, dimension(:) :: ldiff     ! distance from istart to plane (for interpolation)
    type(hist), allocatable, dimension(:) :: hist_t     ! the histograms for each plane
end type hist_zplanes

! Create histogram groups here 
type(hist_zplanes) :: HISTcs2_t   ! SGS coefficient, squared
type(hist_zplanes) :: HISTtn_t    ! Lagrangian time scale
type(hist_zplanes) :: HISTnu_t    ! Eddy viscosity
type(hist_zplanes) :: HISTee_t    ! Error in SGS model

! Create types for outputting data (instantaneous or averaged)
type(point), allocatable, dimension(:) :: point_t
type(plane), allocatable, dimension(:) :: xplane_t, yplane_t
type(zplane), allocatable, dimension(:) :: zplane_t

type(tavg), allocatable, dimension(:,:,:) :: tavg_t
type(tavg), allocatable, dimension(:) :: tavg_zplane_t

$if ($OUTPUT_EXTRA)
type(tavg_sgs), allocatable, dimension(:,:,:) :: tavg_sgs_t
$endif

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
  $if($OUTPUT_EXTRA)
    MODULE PROCEDURE tavg_scalar_div, rs_scalar_div, tavg_sgs_scalar_div
  $else
    MODULE PROCEDURE tavg_scalar_div, rs_scalar_div
  $endif  
END INTERFACE

INTERFACE OPERATOR (.MUL.)
  MODULE PROCEDURE tavg_mul, tavg_scalar_mul
END INTERFACE

INTERFACE type_set
  $if($OUTPUT_EXTRA)
    MODULE PROCEDURE tavg_set, rs_set, tavg_sgs_set
  $else
    MODULE PROCEDURE tavg_set, rs_set
  $endif  
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

$if($OUTPUT_EXTRA)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tavg_sgs_scalar_div( a, b ) result(c)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none

type(tavg_sgs), intent(in) :: a
real(rprec), intent(in) :: b
type(tavg_sgs) :: c

c % Tn = a % Tn / b
c % Nu_t = a % Nu_t / b
c % F_LM = a % F_LM / b
c % F_MM = a % F_MM / b
c % F_QN = a % F_QN / b
c % F_NN = a % F_NN / b
c % ee_now = a % ee_now / b
$if($DYN_TN)
c % F_ee2 = a % F_ee2 / b
c % F_deedt2 = a % F_deedt2 / b
$endif

return
end function tavg_sgs_scalar_div
$endif

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

$if($OUTPUT_EXTRA)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_sgs_set( c, a )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
implicit none
real(rprec), intent(in) :: a
type(tavg_sgs), intent(out) :: c

c % Tn =  a
c % Nu_t =  a
c % F_LM =  a
c % F_MM =  a
c % F_QN =  a
c % F_NN =  a
c % ee_now = a
$if($DYN_TN)
c % F_ee2 = a
c % F_deedt2 = a
$endif

return
end subroutine tavg_sgs_set
$endif

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

type(hist), intent(inout) :: a                
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

type(hist), intent(inout) :: a                
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

type(hist), intent(inout) :: a                
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

