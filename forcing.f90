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

!**********************************************************************
module forcing
!**********************************************************************
!
! Provides subroutines and functions for computing forcing terms on the
! velocity field. Provides driver routine for IBM forces
! (forcing_induced), driver routine for RNS, turbine, etc. forcing
! (forcing_applied), and for the projection step. Also included are
! routines for enforcing a uniform inflow and the fringe region
! treatment.
!

#ifdef PPHIT
use hit_inflow, only : inflow_HIT
#endif

implicit none

save

private

public :: forcing_applied, &
          forcing_induced, &
          project

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine forcing_applied()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This subroutine acts as a driver for applying pointwise body forces
!  into the domain. Subroutines contained here should modify f{x,y,z}a
!  which are explicitly applied forces. These forces are applied to RHS 
!  in the evaluation of u* so that mass conservation is preserved.
!
use types, only : rprec

#ifdef PPTURBINES
use sim_param, only : fxa, fya, fza
use turbines, only:turbines_forcing
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tony ATM
#ifdef PPATM
use sim_param, only : fxa, fya, fza ! The body force components
use atm_lesgo_interface, only : atm_lesgo_forcing
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tony ATM

implicit none

#ifdef PPTURBINES
! Reset applied force arrays
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
call turbines_forcing ()
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tony ATM
#ifdef PPATM
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
call atm_lesgo_forcing ()
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tony ATM
   
end subroutine forcing_applied

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine forcing_induced()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  
!  These forces are designated as induced forces such that they are 
!  chosen to obtain a desired velocity at time
!  step m+1. If this is not the case, care should be taken so that the forces
!  here are divergence free in order to preserve mass conservation. For 
!  non-induced forces such as explicitly applied forces they should be 
!  placed in forcing_applied.
!  
use types, only : rprec
#ifdef PPLVLSET
use level_set, only : level_set_forcing
use sim_param, only : fx, fy, fz
#endif
implicit none

#ifdef PPLVLSET
! Initialize
fx = 0._rprec
fy = 0._rprec
fz = 0._rprec
!  Compute the level set IBM forces
call level_set_forcing ()

#endif

return
end subroutine forcing_induced

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine project ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! provides u, v, w at 1:nz 
!
use param
use sim_param
use messages
use inflow
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#endif
implicit none

integer :: jx, jy, jz
integer :: jz_min

real (rprec) :: RHS, tconst

#ifdef PPVERBOSE
character(*), parameter :: sub_name='project'
#endif

! Caching
tconst = tadv1 * dt

do jz = 1, nz - 1
  do jy = 1, ny
    do jx = 1, nx
 
#ifdef PPLVLSET 
      RHS = -tadv1 * dpdx(jx, jy, jz)
      u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS + fx(jx, jy, jz)))
      RHS = -tadv1 * dpdy(jx, jy, jz)
      v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS + fy(jx, jy, jz))) 
#else
      RHS = -tadv1 * dpdx(jx, jy, jz)
      u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS                 ))
      RHS = -tadv1 * dpdy(jx, jy, jz)
      v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS                 )) 
#endif

    end do
  end do
end do

if (coord == 0) then
  jz_min = 2
else
  jz_min = 1
end if

do jz = jz_min, nz - 1
  do jy = 1, ny
    do jx = 1, nx

#ifdef PPLVLSET 
      RHS = -tadv1 * dpdz(jx, jy, jz)
      w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS + fz(jx, jy, jz)))
#else
      RHS = -tadv1 * dpdz(jx, jy, jz)
      w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS                 ))
#endif

    end do
  end do
end do

! Apply inflow conditions in necessary
call inflow_apply()

! This stuff must be last, so BCs are still enforced, no matter what
! inflow_apply() does
#ifdef PPMPI
! Exchange ghost node information (since coords overlap)                     
call mpi_sync_real_array(u, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(v, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)  
#endif

!--enfore bc at top
#ifdef PPMPI
if (coord == nproc-1) then
#endif
  ! no-stress top
  u(:,:,nz)=u(:,:,nz-1)
  ! no-stress top
  v(:,:,nz)=v(:,:,nz-1)
  ! no permeability
  w(:, :, nz)=0._rprec
#ifdef PPMPI
endif
#endif

if (coord == 0) then

  ! No modulation of u and v since if a stress free condition (lbc_mom=0) is
  ! applied, it is applied through the momentum equation. 

  ! no permeability
  w(:, :, 1)=0._rprec

end if

end subroutine project

end module forcing
