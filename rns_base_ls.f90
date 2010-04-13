!**********************************************************************
module rns_base_ls
!**********************************************************************
use types, only : rprec
use param, only : ld, ny, nz
implicit none

save
public

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

logical, parameter :: clforce_calc = .true.
integer, parameter :: clforce_nskip = 10
logical, parameter :: clforce_vel_write = .true.


logical, parameter :: brforce_calc = .false.



!  Flag for writing info (forces, velocity, etc.) on tree 1 (main) only
logical, parameter :: use_main_tree_only = .true.

!type rns
!  integer :: ntrees, nplanes
!  logical :: plane_u_calc
!end type rns

type ref_plane
  integer :: nzeta, neta ! discretization
  real(rprec), dimension(3) :: p1, p2, p3 !  3 ordered points
  real(rprec) :: u ! reference values
  real(rprec) :: area
end type ref_plane

type force
  real(rprec) :: fD
  real(rprec) :: CD
end type force

type(ref_plane), pointer, dimension(:) :: cl_ref_plane_t
type(force), pointer, dimension(:) :: brforce_t, clforce_t

integer :: brindx(ld, ny, $lbz:nz)
logical :: brindx_initialized = .false.

real (rprep) :: chi(ld, ny, $lbz:nz)
logical :: chi_initialized = .false.

end module rns_base_ls
