module level_set_base
use types, only : rp => rprec
use types, only : rprec
use param, only : ld, ny, nz, dx, lbz
implicit none

save

public
private :: rp, ld, ny, nz, dx

!private
!public :: phi

logical, parameter :: global_CD_calc = .false. ! Compute global CD based on inflow velocity
integer, parameter :: Ldir = 2
                      !--lift direction:
                      !  2 when lift direction is y
                      !  3 when lift direction is z


logical, parameter :: vel_BC = .false. !--means we are forcing velocity for
                                       !  level set BC 
                                       !  (default = .false.)
logical, parameter :: use_log_profile = .false.       !  (default = .false.)
logical, parameter :: use_enforce_un = .false.        !  (default = .false.)
logical, parameter :: physBC = .true.                 !  (default = .true.)        
logical, parameter :: use_smooth_tau = .true.         !  (default = .true.)
logical, parameter :: use_extrap_tau_log = .false.    !  (default = .false.)
logical, parameter :: use_extrap_tau_simple = .true.  !  (default = .true.)
logical, parameter :: use_modify_dutdn = .false.  !--only works w/interp_tau; not MPI compliant
                                                  !--wont work w/extra_tau_log
                                                  !  (default = .false.)

! Enables scale dependent Cs evaluations (not dynamic evaluation)
! Used when model=4 in param module
logical, parameter :: lag_dyn_modify_beta = .true.

! Configures the mode in which SOR smoothing is applied in the IB
! 'xy' may be safely used in most cases (must be used for MPI cases)
! '3d' not MPI compliant
character (*), parameter :: smooth_mode = 'xy'  !--'xy', '3d'

real (rp), parameter :: z0 = 0.0001_rp !--nondimensional roughness length of surface
!real (rp), parameter :: z0 = 1.6e-3_rp ! rough

logical :: phi_cutoff_is_set = .false.
logical :: phi_0_is_set = .false.


real (rp) :: phi(ld, ny, lbz:nz)

end module level_set_base
