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

real(rprec), parameter :: chi_cutoff = 1.e-9_rprec



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
  integer :: parent !  parent CD; for resolved branches 
  real(rprec) :: fD
  real(rprec) :: CD
end type force

type indx_array
  integer :: npoint
  integer, pointer, dimension(:,:) :: iarray
end type

type(indx_array), pointer, dimension(:) :: cl_indx_array

type(ref_plane), pointer, dimension(:) :: cl_ref_plane_t
type(force), pointer, dimension(:) :: brforce_t, clforce_t

integer :: ncluster_unreslv ! total number of unresolved clusters
integer :: ncluster_reslv ! total number of resolved clusters
integer :: ncluster_reslv_ref
integer :: ncluster_ref ! number of clusters used for computing reference quantities (size of cl_ref_plane_t)
integer :: ncluster_tot ! total number of clusters for all trees
integer :: ntree_ref !  number of trees used in reference calculations 

integer, target :: brindx(ld, ny, $lbz:nz)
logical :: brindx_initialized = .false.

real (rprec) :: chi(ld, ny, $lbz:nz)
logical :: chi_initialized = .false.

end module rns_base_ls
