!**********************************************************************
module stat_defs
!**********************************************************************
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
type taver
  logical :: calc,started
  integer :: nstart, nend
  double precision, allocatable, dimension(:,:,:) :: u, v, w, &
                                                     u2, v2, w2, &
      											     uw, vw, uv, dudz
end type taver	
  
!  Instantaneous Variables Storage (Parameters for storing velocity 
!  componentsvalues each time step)
type ui_pnt
  logical :: calc,started,global
  integer :: nstart, nend, nloc, nskip
  integer :: ijk(3,10) !  Can specify up to 10 points to record 
end type

!  Instantaneous velocity global declarations
type ui_gbl
  logical :: calc,started,global
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
  
logical :: aver_calc
 
type(stats)          :: stats_t
type(rs)             :: rs_t
type(taver)          :: taver_t
type(ui_pnt), target :: ui_pnt_t
type(ui_gbl)         :: ui_gbl_t
type(plane)		     :: yplane_t, zplane_t
  
contains
!***************************************************************
double precision function interp_to_uv_grid(var,i,j,k)
!***************************************************************
!  This function computes any values the read in value u1(k) and
!  u2(k+1) to the w grid location k
use param2,only : nz
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

!**********************************************************************  
module Sij_defs
!**********************************************************************
!  This module contains the allocatable arrays for the sgs_stag routine
  use param, only : USE_MPI, coord
  use param2, only : ld, nx, ny, nz, lbc_mom
  implicit none
  double precision, allocatable, dimension (:, :, :) :: S11, S12, S22, S33, S13, S23
end module Sij_defs

!**********************************************************************
module convec_defs
!**********************************************************************
!  This module contains the allocatable arrays for the subroutine convec 
!  in convec.f90
use types,only:rprec
implicit none
!  This module contains the allocatable arrays for the subroutine convec 
!  in convec.f90

!--save forces heap storage
real(kind=rprec), save, allocatable, dimension(:,:,:)::cc_big
!--save forces heap storage
real (rprec), save, allocatable, dimension (:,:,:) :: u1_big, u2_big, u3_big
!--MPI: only u1_big(0:nz-1), u2_big(0:nz-1), u3_big(1:nz) are used
!--save forces heap storage 
real (rprec), save, allocatable,dimension (:, :, :) :: vort1_big, vort2_big, vort3_big

end module convec_defs

!**********************************************************************
module test_filtermodule_defs
!**********************************************************************
!  This module contains the allocatable arrays from the test_filermodule
!  module
use types,only:rprec
implicit none

real(kind=rprec),allocatable, dimension(:,:) :: G_test,G_test_test

end module test_filtermodule_defs

!**********************************************************************
module immersedbc_defs
!**********************************************************************
!  This module contains the allocatable arrays from the immersedbc module
use types, only : rprec
implicit none

end module immersedbc_defs
