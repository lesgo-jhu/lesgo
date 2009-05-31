module stat_defs
!  Define parameters for writing statistics
type stats
  logical :: calc
  integer :: nstart, nend !  Time step when to start and stop averaging
end type stats

!  Reynolds stresses
type rs
  logical :: calc
  double precision, allocatable, dimension(:,:,:) :: up2, vp2, wp2, & 
                                                     upwp, vpwp, upvp
end type rs

!  Sums performed over time
type tavg
  logical :: calc,started
  integer :: nstart, nend
  double precision, allocatable, dimension(:,:,:) :: u, v, w, &
    u2, v2, w2, uw, vw, uv, dudz
end type tavg	
  
!  Instantaneous Variables Storage (Parameters for storing velocity 
!  componentsvalues each time step)
type point
  logical :: calc,started,global
  integer :: nstart, nend, nloc, nskip
  integer :: ijk(3,10) !  Can specify up to 10 points to record 
end type

!  Instantaneous velocity global declarations
type global
  logical :: calc,started
  integer :: nstart, nend, nskip
end type   
  
!  Planar stats/data
type plane
  logical :: avg
  integer :: na, nstart, nend
  integer, dimension(10) :: istart
  double precision :: fa
  double precision, dimension (10) :: la, ldiff
  double precision, allocatable, dimension(:,:,:) :: ua, va, wa
end type	
  
type(stats)          :: stats_t
type(rs)             :: rs_t
type(tavg)          :: tavg_t
type(point), target :: upoint_t
type(global)         :: uglobal_t
type(plane)		     :: yplane_t, zplane_t
  
contains
!***************************************************************
double precision function interp_to_uv_grid(var,i,j,k)
!***************************************************************
!  This function computes any values the read in value u1(k) and
!  u2(k+1) to the w grid location k
use param,only : nz
use sim_param, only : w, dudz

character(*), intent(IN) :: var
integer,intent(IN) :: i,j,k

if(trim(adjustl(var)) == 'w') then
  if(k==nz) then
    interp_to_uv_grid = 3./2.*w(i,j,k) - 0.5*w(i,j,k-1)
  else
    interp_to_uv_grid = 0.5*(w(i,j,k)+w(i,j,k+1))
  endif
elseif(trim(adjustl(var)) == 'dudz') then 
  if(k==nz) then
    interp_to_uv_grid = 3./2.*dudz(i,j,k) - 0.5*dudz(i,j,k-1)
  else
    interp_to_uv_grid = 0.5*(dudz(i,j,k)+dudz(i,j,k+1))
  endif
else
  write(*,*) 'Error: variable specification not specified properly!'
  stop
endif
return
end function

end module stat_defs

module Sij_defs
!  This module was created as a work around of a memory issue
!  previously seen for large grids
  use param, only : ld, nx, ny, nz, USE_MPI, coord, lbc_mom
  double precision, dimension (ld, ny, nz) :: S11, S12, S22, S33, S13, S23
end module Sij_defs

