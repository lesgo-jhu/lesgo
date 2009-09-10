!**********************************************************************
module stat_defs
!**********************************************************************

save
public

!  Reynolds stresses
type rs
  logical :: calc
  double precision, pointer, dimension(:,:,:) :: up2, vp2, wp2, & 
                                                     upwp, vpwp, upvp
end type rs

!  Sums performed over time
type tstats
  logical :: calc, started
  integer :: nstart, nend
  double precision, pointer, dimension(:,:,:) :: u, v, w, &
    u2, v2, w2, uw, vw, uv, dudz
end type tstats	
  
!  Instantaneous Variables Storage (Parameters for storing velocity 
!  component values each time step)
type point
  logical :: calc,started
  integer :: nstart, nend, nloc, nskip
  integer, dimension(10) :: coord, istart, jstart, kstart
  double precision, dimension(10) :: xdiff, ydiff, zdiff
  double precision, dimension(3,10) :: xyz
end type

!  Instantaneous velocity global declarations
type domain
  logical :: calc,started
  integer :: nstart, nend, nskip
end type domain  
  
!  Planar stats/data
type plane
  logical :: calc
  integer :: nloc, nstart, nend
  integer, dimension(10) :: istart, coord
  double precision :: fa
  double precision, dimension (10) :: loc, ldiff
  double precision, pointer, dimension(:,:,:) :: ua, va, wa
end type plane
  
type(rs)             :: rs_t
type(tstats)	     :: tavg_t
type(tstats)          :: tsum_t
type(point), target :: point_t
type(domain)         :: domain_t
type(plane)     :: yplane_t, zplane_t

end module stat_defs

