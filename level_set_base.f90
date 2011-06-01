module level_set_base
use types, only : rp => rprec
use types, only : rprec
use param, only : ld, ny, nz, dx
implicit none

!--this definition of lbz must be consistent with that in level_set module
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif


save

public
private :: rp, ld, ny, nz, dx

!private
!public :: phi

logical, parameter :: global_CD_calc = .true. ! Compute global CD based on inflow velocity
integer, parameter :: Ldir = 2
                      !--lift direction:
                      !  2 when lift direction is y
                      !  3 when lift direction is z


logical, parameter :: vel_BC = .false. !--means we are forcing velocity for
                                       !  level set BC
logical, parameter :: use_log_profile = .false.
logical, parameter :: use_enforce_un = .false.
logical, parameter :: physBC = .true.
logical, parameter :: use_smooth_tau = .true.
logical, parameter :: use_extrap_tau_log = .false.
logical, parameter :: use_extrap_tau_simple = .true.
logical, parameter :: use_modify_dutdn = .false.  !--only works w/interp_tau; not MPI compliant
                                                  !--wont work w/extra_tau_log

real (rp), parameter :: z0 = 0.0001_rp
                        !--nondimensional roughness length of surface

logical :: phi_cutoff_is_set = .false.
logical :: phi_0_is_set = .false.


real (rp) :: phi(ld, ny, $lbz:nz)

end module level_set_base
