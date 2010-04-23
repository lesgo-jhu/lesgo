!**********************************************************************
module stat_defs
!**********************************************************************
use types, only : rprec
save
public

!  Reynolds stresses
type rs
  logical :: calc=.false.
  real(rprec), pointer, dimension(:,:,:) :: up2, vp2, wp2, & 
                                            upwp, vpwp, upvp
end type rs

!  Sums performed over time
type tstats
  logical :: calc=.false.
  logical :: started=.false.
  integer :: nstart, nend, nskip
  real(rprec) :: total_time
  real(rprec), pointer, dimension(:,:,:) :: u, v, w, &
    u2, v2, w2, uw, vw, uv, dudz
    
  $if($LVLSET)
  $if($RNS_LS)
  real(rprec), pointer, dimension(:,:,:) :: fx, fy, fz
  $endif
  $endif
  
  real(rprec), pointer, dimension(:) :: u_avg, v_avg, w_avg, u2_avg, v2_avg, w2_avg
  real(rprec), pointer, dimension(:) :: txx_avg, txy_avg, tyy_avg, txz_avg, tyz_avg, tzz_avg
  real(rprec), pointer, dimension(:) :: dudz_avg, dvdz_avg, uv_avg, uw_avg, vw_avg
end type tstats	
  
!  Instantaneous Variables Storage (Parameters for storing velocity 
!  component values each time step)
type point
  logical :: calc=.false.
  logical :: started=.false.
  integer :: nstart, nend, nloc, nskip
  integer, dimension(10) :: coord=-1, istart=-1, jstart=-1, kstart=-1
  real(rprec), dimension(10) :: xdiff, ydiff, zdiff
  real(rprec), dimension(3,10) :: xyz=-1._rprec
  character(64), dimension(10) :: fname
end type

!  Instantaneous velocity global declarations
type domain
  logical :: calc=.false.
  logical :: started=.false.
  integer :: nstart, nend, nskip
end type domain  
  
!  Planar stats/data
type plane
  logical :: calc=.false.
  logical :: started=.false.
  integer :: nloc, nstart, nend, nskip
  integer, dimension(10) :: istart=-1, coord=-1
!   real(rprec) :: fa
  real(rprec), dimension (10) :: loc, ldiff
!   real(rprec), pointer, dimension(:,:,:) :: ua, va, wa
end type plane

$if ($TURBINES)
	type turbine 
		real :: xloc, yloc, height, dia, thk		  
		real :: theta1                              !angle CCW(from above) from -x direction [degrees]
        real :: theta2	                            !angle above the horizontal, from -x dir [degrees]
		real, dimension(3) :: nhat                  !(nx,ny,nz) of unit normal for each turbine
		integer :: num_nodes                        !number of nodes associated with each turbine
		integer, dimension(1500,3) :: nodes         !(i,j,k) of each included node
        integer, dimension(6) :: nodes_max          !search area for nearby nodes
		real :: u_d, u_d_T                          !running time-average of mean disk velocity
        real :: f_n                                 !normal force on turbine disk
        integer :: u_d_flag
		real, dimension(1500) :: ind                !indicator function - weighting of each node
	end type turbine

	type wind_farm
		integer ::ifilter                           !Filter type: 2->Gaussian
		real :: alpha                               !filter size is alpha*(grid spacing)
        integer :: trunc                            !truncated - # grid points to include in each dir.
        real :: filter_cutoff                       !ind only includes values above this cutoff
		type(turbine), pointer, dimension(:) :: turbine_t
	end type wind_farm	
    
    type(wind_farm)	        :: wind_farm_t	
$endif

  
type(rs)            		:: rs_t
type(tstats)        		:: tavg_t, tstats_t
type(tstats)        	 	:: zplane_avg_t
type(point), target 	  :: point_t
type(domain)        		:: domain_t
type(plane)         		:: yplane_t, zplane_t

end module stat_defs

