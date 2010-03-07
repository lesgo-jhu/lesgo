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
  real(rprec), pointer, dimension(:,:,:) :: u, v, w, &
    u2, v2, w2, uw, vw, uv, dudz
end type tstats	
  
!  Instantaneous Variables Storage (Parameters for storing velocity 
!  component values each time step)
type point
  logical :: calc=.false.
  logical :: started=.false.
  integer :: nstart, nend, nloc, nskip
  integer, dimension(10) :: coord, istart, jstart, kstart
  real(rprec), dimension(10) :: xdiff, ydiff, zdiff
  real(rprec), dimension(3,10) :: xyz
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
  integer, dimension(10) :: istart, coord
!   real(rprec) :: fa
  real(rprec), dimension (10) :: loc, ldiff
!   real(rprec), pointer, dimension(:,:,:) :: ua, va, wa
end type plane
  
type(rs)            		:: rs_t
type(tstats)        		:: tavg_t
type(tstats)        	 	:: tsum_t
type(point), target 	:: point_t
type(domain)        		:: domain_t
type(plane)         		:: yplane_t, zplane_t

end module stat_defs

