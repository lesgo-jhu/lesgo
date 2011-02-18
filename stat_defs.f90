!**********************************************************************
module stat_defs
!**********************************************************************
use types, only : rprec
use param, only : nx,ny,nz,lh

save
public

!integer, dimension(point_nloc) :: point_coord=-1, point_istart=-1, point_jstart=-1, point_kstart=-1
!real(rprec), dimension(point_nloc) :: point_xdiff, point_ydiff, point_zdiff
!character(64), dimension(point_nloc) :: point_fname

!integer, dimension(xplane_nloc) :: xplane_istart=-1
!real(rprec), dimension(xplane_nloc) :: xplane_ldiff

!integer, dimension(yplane_nloc) :: yplane_istart=-1
!real(rprec), dimension(yplane_nloc) :: yplane_ldiff

!integer, dimension(zplane_nloc) :: zplane_istart=-1, zplane_coord=-1
!real(rprec), dimension(zplane_nloc) :: zplane_ldiff

!integer, dimension(spectra_nloc) :: spectra_istart=-1, spectra_coord=-1
!real(rprec), dimension(spectra_nloc) :: spectra_ldiff

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
  real(rprec), dimension(nx) :: uhat
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

end module stat_defs

