!**********************************************************************
module stat_defs
!**********************************************************************
use types, only : rprec
use param, only : nx,ny,nz

save
public

type rs
  real(rprec) :: up2, vp2, wp2, upwp, vpwp, upvp
end type rs

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
		real(rprec) :: theta1                       !angle CCW(from above) from -x direction [degrees]
        real(rprec) :: theta2	                    !angle above the horizontal, from -x dir [degrees]
		real(rprec), dimension(3) :: nhat           !(nx,ny,nz) of unit normal for each turbine
		integer :: num_nodes                        !number of nodes associated with each turbine
		integer, dimension(1500,3) :: nodes         !(i,j,k) of each included node
        integer, dimension(6) :: nodes_max          !search area for nearby nodes
		real(rprec) :: u_d, u_d_T                   !running time-average of mean disk velocity
        real(rprec) :: f_n                          !normal force on turbine disk
		real(rprec), dimension(1500) :: ind         !indicator function - weighting of each node
        
        real(rprec), dimension (nx,ny,nz) :: u_cond_avg_lo, v_cond_avg_lo, w_cond_avg_lo
        real(rprec), dimension (nx,ny,nz) :: u_cond_avg_hi, v_cond_avg_hi, w_cond_avg_hi
        logical :: cond_avg_calc_lo,cond_avg_calc_hi
        real(rprec) :: cond_avg_ud_lo,cond_avg_ud_hi,cond_avg_time_lo,cond_avg_time_hi
	end type turbine

	type wind_farm
		integer ::ifilter                           !Filter type: 2->Gaussian
		real(rprec) :: alpha                        !filter size is alpha*(grid spacing)
        integer :: trunc                            !truncated - # grid points to include in each dir.
        real(rprec) :: filter_cutoff                !ind only includes values above this cutoff
		type(turbine), pointer, dimension(:) :: turbine_t
        
        logical, pointer, dimension(:) :: cond_avg_flag_lo,cond_avg_flag_hi
	end type wind_farm	
    
  type(wind_farm)	        :: wind_farm_t	
$endif

type(tavg), allocatable, dimension(:,:,:) :: tavg_t
type(tavg), allocatable, dimension(:) :: tavg_zplane_t

type(rs), allocatable, dimension(:,:,:) :: rs_t
type(rs), allocatable, dimension(:) :: rs_zplane_t, cnpy_zplane_t

end module stat_defs

