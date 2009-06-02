!**********************************************************************
module trees_base_ls
!**********************************************************************
use types, only : rp => rprec
!use precision
use messages
implicit none

save
public

private :: mod_name

character (*), parameter :: mod_name = 'trees_base_ls'

character (*), parameter :: branch_cross_section = 'circular'
                            !--'circular', 'square'
                            !--only for trees_pre: 'square+plate'

integer, parameter :: nd = 3  ! number of dimensions
integer, parameter :: nzone = 1  ! # of averaging zones for dyn. Cd
integer, parameter :: tree_node = 1  ! trees based on u, v, or w nodes

!--specify force model to be used
$define $fmodel nba
                !--d, d_germano, dls, nba
$if ($fmodel eq "d")
  $define $nfcoeff 1
$elsif ($fmodel eq "d_germano")
  $define $nfcoeff 1
$elsif ($fmodel eq "dls")
  $define $nfcoeff 3
$elsif ($fmodel eq "nba")
  $define $nfcoeff 3
$else
  $error (Invalid force model specification)
$endif

character (*), parameter :: fmodel = $str($fmodel)
integer, parameter :: nfcoeff = $nfcoeff

logical, parameter :: DEBUG = .false.
logical, parameter :: VERBOSE = .false.

!--beware, some of these are no longer in use
logical, parameter :: use_tecplot = .true.
!logical, parameter :: use_Cd_dynamic = .true.
!logical, parameter :: use_Cd_with_resolved = .false.
!logical, parameter :: use_measure_cyl = .false.
!logical, parameter :: use_Cd_apriori = .false.
!logical, parameter :: use_res_f_smoothed = .false.
!logical, parameter :: use_res_f_test = .true.
!logical, parameter :: use_unres_f_test = .false.  !--F for DNS, T for LES
!logical, parameter :: chop_base = .true.
logical, parameter :: add_cap = .true.
logical, parameter :: add_base = .true.
!logical, parameter :: use_term_check = .false.
!logical, parameter :: add_k_layer = .true.
logical, parameter :: sub_branches_outside = .true.

character (*), parameter :: cap_shape = 'rectangular'
                            !--'rectangular', 'hemispherical'
character (*), parameter :: base_shape = 'rectangular'
                            !--'rectangular', 'hemispherical'

real (rp), parameter :: flow_dir(nd) = (/ 1._rp, 0._rp, 0._rp /)
                        !--for now a parameter, could change later

integer :: n_tree = -1  !--to be read from trees.conf file
integer :: n_br = -1

!--use a lightweight definition of branches
type branch_type
  integer :: ident = -1  ! branch identification number
                         ! note: could encode useful info in here
  integer :: itree = -1  !--tell which tree this branch belongs to
  integer :: gen = -1  !--number of "general" parents this br. has
                       !--i.e. the branch generation
  integer :: n_sub_branch = -1  ! number of children this br. has
  integer :: zone = -1 ! zone of this branch, if using multi-zone avging

  integer :: nrespt = -1  !--number grid pts inside this (resolved) branch

  integer :: nbboxpt = -1  !--number of pts in bbox
  integer, pointer :: bboxpt(:, :) => NULL ()

  logical :: resolved  ! whether or not this br. is resolved

  !--defaults should be assigned now, any non-defaults will be read
  !  in read_trees_conf
  real (rp) :: l = 0._rp  !--length of branch
  real (rp) :: d = 0._rp        !--diameter of branch
  real (rp) :: taper = 0._rp    !--d_tip = (1-taper) * d
  real (rp) :: twist = 0._rp    !--extra rotation about default coord. sys.
  real (rp) :: width_bbox = 0._rp
  real (rp) :: height_bbox = 0._rp
  real (rp) :: root_height ! this br. root, fraction of parent hght

  real (rp) :: fcoeff(nfcoeff) = -1._rp  !--local value of force coefficients
  real (rp) :: fdist_coeff(nfcoeff) = -1.0_rp
  real (rp) :: A(nfcoeff) = -1._rp
  real (rp) :: fcoeff_dir(nd, nfcoeff) = 0._rp
  real (rp) :: fnorm(nd, nfcoeff) = -1.0_rp

  real (rp) :: Mdyn (nd, nzone) = 0._rp

  real (rp) :: x0(nd)      ! coords of root of br.

  real (rp) :: abs_dir(nd) ! unit vector in dir of br. (abs coords)
  real (rp) :: rel_dir(nd) ! unit vctr in dir of br. (parent coords)

  real (rp) :: resf(nd) = 0._rp    !--res. force this br.
  real (rp) :: unresf(nd) = 0._rp  !--unres. force this br.

  real (rp) :: resftot(nd) = 0._rp    !--total forces on this branch
  real (rp) :: unresftot(nd) = 0._rp  !  cumulative over sub-branches
  real (rp) :: ftot(nd) = 0._rp       !--ftot = resftot + unresftot

  real (rp) :: velscale(nd) = 0._rp ! avg. vel. in bounding box

  ! unit vectors holding this branch's coordinate directions
  ! in the absolute coordinate frame
  real (rp) :: x_hat(nd), y_hat(nd), z_hat(nd)

  type (branch_type), pointer :: sub_branch(:) => NULL ()

  type (branch_type), pointer :: parent_branch => NULL ()

end type branch_type

!--idea is to define defaults here, which may be optionally changed in
!  trees.conf file, these will then be applied to tree % trunk
type tree_type

  !--these defaults may be changed upon reading trees.conf
  integer :: n_gen = 0 !--number of branch levels/generations, trunk is 0
  integer :: n_sub_branch = 0 !--number of (direct) sub branches this tree has
  integer :: max_res_gen = 0 !--maximum resolved generation
  integer :: trunk_dir(nd) = (/ 0, 0, 1 /)
                      !--direction trunks grow in

  !--these are optionally specified in trees.conf file
  real (rp) :: l = 1._rp
  real (rp) :: d = 0.1_rp
  real (rp) :: ratio = 1._rp  !--this is ratio next-gen/this-gen, i.e. r
  real (rp), pointer :: rel_dir(:, :) => NULL ()  !--nd X n_sub_branch
  real (rp), pointer :: root_height(:) => NULL ()  !--size n_sub_branch
  real (rp) :: taper = 0._rp
  real (rp) :: trunk_twist = 0._rp
  real (rp), pointer :: twist(:) => NULL ()  !--size n_sub_branch
  real (rp) :: x0(nd) = 0._rp

  !type (branch_type) :: trunk
  type (branch_type), pointer :: trunk => NULL ()

end type tree_type

! tree array for simulations
!type (tree_type), target :: tree_array(n_tree)
$if ($XLF)
  type (tree_type), allocatable :: tree_array(:)  !--experimental

  !type (tree_type), save :: tree_array(n_tree)
                    !--xlf want save here, even though its at top of module
$else
  type (tree_type), allocatable :: tree_array(:)  !--experimental
  !type (tree_type) :: tree_array(n_tree)
$endif

!--branch array (more convenient to access than linked list)

type (branch_type), allocatable, target :: branch_array(:)

logical :: tree_array_initialized = .false.

! simple descriptor for the grid
! in general: x_i = x_min + (i-1) * dx, for values of i between 1 and nx
type grid_type

  logical :: initialized = .false.

  logical :: staggered(nd)  ! indicates which nodes are staggered

  integer :: nx(nd)

  real (rp) :: x_min(nd, nd)  ! x_min(:, 1) are x_mins of u nodes, etc.
  real (rp) :: dx(nd)

end type grid_type

! global copy of grid info
$if ($XLF)
  type (grid_type), save :: grid  !--xlf want save here
$else
  type (grid_type) :: grid
$endif


contains

!**********************************************************************
function cross_product (a, b)
!**********************************************************************
!
!  This function computes the cross product of vectors a and b (a x b) 
!  with dimensions nd.
!

implicit none

real (rp) :: cross_product(nd)

real (rp), intent (in) :: a(nd), b(nd)

!----------------------------------------------------------------------

cross_product(1) = a(2) * b(3) - a(3) * b(2)
cross_product(2) = a(3) * b(1) - a(1) * b(3)
cross_product(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product

!**********************************************************************
function delta (h, x)
!**********************************************************************
implicit none

real (rp) :: delta

real (rp), intent (in) :: h, x

real (rp), parameter :: pi = 3.14159265359_rp

!----------------------------------------------------------------------

if (abs (x) <= 2._rp * h) then

  delta = 0.25_rp * (1._rp + cos (0.5_rp * pi * x / h)) / h

else

  delta = 0._rp

end if

end function delta

!**********************************************************************
subroutine grid_initialize ()
!**********************************************************************
use param, only : nx, ny, nz, dx, dy, dz
implicit none

character (*), parameter :: sub_name = mod_name // '.grid_initialize'

integer :: i, j
integer :: tmp(nd), not_i(nd-1)

real (rp), parameter :: thresh = 10._rp * epsilon (1._rp)

!----------------------------------------------------------------------
write(*,*) 'From trees_base_ls.grid_initialize, dx, dy, dz =', dx,dy,dz

if (DEBUG) call enter_sub (sub_name)

!  Set grid dimensions to those of global values
grid % nx = (/ nx, ny, nz /)

! currently, we make assumption that dx=dy=dz, so we had better
! enforce it until this changes
if ((abs (dx-dy) > thresh ) .or. (abs (dx-dz) > thresh) .or.  &
    (abs (dy-dz) > thresh)) then

  call error (sub_name, 'tree module requires dx=dy=dz for now')

end if

!  Set grid spacing to that of global values
grid % dx = (/ dx, dy, dz /)

! u-nodes
grid % x_min (:, 1) = (/ 0._rp, 0._rp, dz / 2._rp /)
! v-nodes
grid % x_min (:, 2) = (/ 0._rp, 0._rp, dz / 2._rp /)
! w-nodes
grid % x_min (:, 3) = (/ 0._rp, 0._rp, 0._rp /)

!  Set the grid staggered flags; u, v staggered in
!  the z direction
grid % staggered(1:2) = .true.
grid % staggered(3) = .false.

! ! determine which nodes are staggered
! !  Loop over number of dimensions
! do i = 1, nd
!   write(*,*) 'j = ', j
!   tmp = (/ ( j, j=1, nd ) /)
!   write(*,*) 'tmp = ', tmp
!   not_i = pack (tmp, tmp /= i)
!   write(*,*) 'not_i = ', not_i
! 
!   ! this is the definition of staggered (i.e. our convention)
!   ! note (-) in front of dx/2 part--could arguably be a (+)
!   ! this would change (i, i+1) pairs for interp to (i-1, i), I think
!   if ((abs (grid % x_min(i, not_i(1)) -                                  &
!             grid % x_min(i, not_i(2))) < thresh) .and.                   &
!       (abs (grid % x_min(i, not_i(1)) -                                  &
!             (grid % x_min(i, i) - (grid % dx(i)) / 2._rp)) < thresh)) then
! 
!     grid % staggered(i) = .true.
! 
!   else
! 
!     grid % staggered(i) = .false.
! 
!   end if
! 
!   pause
! 
! end do


write(*,*) 'nd = ', nd
write(*,*) 'grid%staggered = ', grid%staggered
write(*,*) 'epsilon = ', epsilon(0.)

grid % initialized = .true.

if (DEBUG) call exit_sub (sub_name)

end subroutine grid_initialize


!**********************************************************************
function grid_nearest_of_pt (x, d, node)
!**********************************************************************
! 
! This function finds nearest grid point corresponding to real position 
! x must be a better way to combine the following 2 routines.
!

implicit none

integer :: grid_nearest_of_pt

real (rp), intent (in) :: x

integer, intent (in) :: d, node

character (*), parameter :: sub_name = mod_name // '.grid_nearest_of_pt'

!----------------------------------------------------------------------

if (.not. grid % initialized) then
  call error (sub_name, 'grid not initialized')
end if

if ((d < 1) .or. (d > nd)) then
  call error (sub_name, 'd out of bounds')
end if

if ((node < 1) .or. (node > nd)) then
  call error (sub_name, 'node out of bounds')
end if

grid_nearest_of_pt = nint ( (x - (grid % x_min(d, node))) /  &
                            (grid % dx(d)) ) + 1

end function grid_nearest_of_pt

!**********************************************************************
function grid_of_pt (x, d, node)
!**********************************************************************

implicit none

integer :: grid_of_pt

real (rp), intent (in) :: x

integer, intent (in) :: d, node

character (*), parameter :: sub_name = mod_name // '.grid_of_pt'

!----------------------------------------------------------------------

if (.not. grid % initialized) then
  call error (sub_name, 'grid not initialized')
end if

if ((d < 1) .or. (d > nd)) then
  call error (sub_name, 'd out of bounds')
end if

if ((node < 1) .or. (node > nd)) then
  call error (sub_name, 'node out of bounds')
end if

grid_of_pt = floor ( (x - (grid % x_min(d, node))) / (grid % dx(d)) ) + 1

end function grid_of_pt

!***************************************************************
function mag (a)
!***************************************************************
implicit none

real (rp) :: mag

real (rp), intent (in) :: a(nd)

!----------------------------------------------------------------------

mag = sqrt(dot_product(a, a))

end function mag

!***************************************************************
function pt_of_grid (i, d, node)
!***************************************************************
!  This subroutine computes the x,y,z location for the z-partitioned grid
!
use param, only : coord
implicit none

real (rp) :: pt_of_grid

integer, intent (in) :: i, d, node

character (*), parameter :: sub_name = mod_name // '.pt_of_grid'

!----------------------------------------------------------------------

if (.not. grid % initialized) then
  call error (sub_name, 'grid not initialized')
end if

if ((d < 1) .or. (d > nd)) then
  call error (sub_name, 'grid not intialized')
end if

if ((node < 1) .or. (node > nd)) then
  call error (sub_name, 'node out of bounds')
end if

! we do want to allow this, in periodic cases, at least
!if ( (i < 1) .or. ( i > (grid % nx(d)) ) ) then
!  write (*, *) sub_name // ': i out of bounds'
!  stop
!end if

$if ($MPI)

  if (d == nd) then
    pt_of_grid = (grid % x_min(d, node)) +                           &
                 (coord * (grid % nx(d) - 1) + i - 1) * (grid % dx(d))
  else
    pt_of_grid = (grid % x_min(d, node)) + (i - 1) * (grid % dx(d))
  end if

$else

  pt_of_grid = (grid % x_min(d, node)) + (i - 1) * (grid % dx(d))

$endif


end function pt_of_grid

end module trees_base_ls
