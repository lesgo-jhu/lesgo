!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

module level_set_base
use param, only : rp => rprec
use param, only : rprec
use param, only : ld, ny, nz, dx, lbz
implicit none

save

public
private :: rp, ld, ny, nz, dx, lbz

!private
!public :: phi

!logical, parameter :: global_CD_calc = .true. ! Compute global CD based on inflow velocity
logical :: global_CA_calc = .false. ! Compute global CA based on inflow velocity
integer :: global_CA_nskip = 10     ! Number of time steps to skip between global CA writes

!logical, parameter :: vel_BC = .false.
logical :: vel_BC = .false. !--means we are forcing velocity for
                            !  level set BC
                            !  (default = .false.)
! logical, parameter :: use_log_profile = .false.       !  (default = .false.)
! logical, parameter :: use_enforce_un = .false.        !  (default = .false.)
! logical, parameter :: physBC = .true.                 !  (default = .true.)
! logical, parameter :: use_smooth_tau = .true.         !  (default = .true.)
! logical, parameter :: use_extrap_tau_log = .false.    !  (default = .false.)
! logical, parameter :: use_extrap_tau_simple = .true.  !  (default = .true.)
! logical, parameter :: use_modify_dutdn = .false.  !--only works w/interp_tau; not MPI compliant
!                                                   !--wont work w/extra_tau_log
!                                                   !  (default = .false.)
logical :: use_log_profile = .false.       !  (default = .false.)
logical :: use_enforce_un = .false.        !  (default = .false.)
logical :: physBC = .true.                 !  (default = .true.)
logical :: use_smooth_tau = .true.         !  (default = .true.)
logical :: use_extrap_tau_log = .false.    !  (default = .false.)
logical :: use_extrap_tau_simple = .true.  !  (default = .true.)
logical :: use_modify_dutdn = .false.  !--only works w/interp_tau; not MPI compliant
                                                  !--wont work w/extra_tau_log
                                                  !  (default = .false.)

! ! Enables scale dependent Cs evaluations (not dynamic evaluation)
! ! Used when sgs_model=4 in param module
! logical, parameter :: lag_dyn_modify_beta = .true.

! Enables scale dependent Cs evaluations (not dynamic evaluation)
! Used when sgs_model=4 in param module
logical :: lag_dyn_modify_beta = .true.

! ! Configures the mode in which SOR smoothing is applied in the IB
! ! 'xy' may be safely used in most cases (must be used for MPI cases)
! ! '3d' not MPI compliant
! character (*), parameter :: smooth_mode = 'xy'  !--'xy', '3d'

! Configures the mode in which SOR smoothing is applied in the IB
! 'xy' may be safely used in most cases (must be used for MPI cases)
! '3d' not MPI compliant
character(25) :: smooth_mode = 'xy'  !--'xy', '3d'

!real (rp), parameter :: zo_level_set = 0.0001_rp !--nondimensional roughness length of surface
real (rp) :: zo_level_set = 0.0001_rp !--nondimensional roughness length of surface

logical :: phi_cutoff_is_set = .false.
logical :: phi_0_is_set = .false.


!real (rp) :: phi(ld, ny, lbz:nz)
real(rp), allocatable, dimension(:,:,:) :: phi

logical :: use_trees

#ifdef PPMPI
  ! Make sure all values (top and bottom) are less than Nz
  integer :: nphitop = 3
  integer :: nphibot = 2
  integer :: nveltop = 1
  integer :: nvelbot = 1
  integer :: ntautop = 3
  integer :: ntaubot = 2
  integer :: nFMMtop = 1
  integer :: nFMMbot = 1
#else
  integer, parameter :: nphitop = 0
  integer, parameter :: nphibot = 0
  integer, parameter :: nveltop = 0
  integer, parameter :: nvelbot = 0
  integer, parameter :: ntautop = 0
  integer, parameter :: ntaubot = 0
  integer, parameter :: nFMMtop = 0
  integer, parameter :: nFMMbot = 0
#endif


contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_base_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine initializes all arrays defined in level_set_base
!
implicit none

allocate( phi( ld, ny, lbz:nz ) )

return
end subroutine level_set_base_init


end module level_set_base
