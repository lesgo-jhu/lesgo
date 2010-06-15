!**********************************************************************
module stat_defs
!**********************************************************************
use types, only : rprec

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
		real(rprec), dimension(1500) :: ind                !indicator function - weighting of each node
	end type turbine

	type wind_farm
		integer ::ifilter                           !Filter type: 2->Gaussian
		real(rprec) :: alpha                        !filter size is alpha*(grid spacing)
    integer :: trunc                            !truncated - # grid points to include in each dir.
    real(rprec) :: filter_cutoff                !ind only includes values above this cutoff
		type(turbine), pointer, dimension(:) :: turbine_t
	end type wind_farm	
    
  type(wind_farm)	        :: wind_farm_t	
$endif

type(tavg), pointer, dimension(:,:,:) :: tavg_t
type(tavg), pointer, dimension(:) :: tavg_zplane_t

type(rs), pointer, dimension(:,:,:) :: rs_t
type(rs), pointer, dimension(:) :: rs_zplane_t

$if($MPI)
type(rs), pointer, dimension(:) :: rs_zplane_tot_t
type(rs), pointer, dimension(:) :: rs_zplane_buf_t

type(tavg), pointer, dimension(:) :: tavg_zplane_tot_t
type(tavg), pointer, dimension(:) :: tavg_zplane_buf_t
$endif

end module stat_defs

