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

public :: forcing_random, &
          forcing_applied, &
          forcing_induced, &
          inflow_cond, &
          project

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine forcing_random()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine generates a random body force that is helpful to
! trigger transition at low Re DNS. The forces are applied to RHS in
! evaluation of u* (not at wall) so that mass conservation is preserved.
!
use types, only : rprec
use param, only : nx,ny,nz,coord
use sim_param, only : RHSy, RHSz

real (rprec) :: rms_forcing, dummy_rand ! TODO make this a user input
integer :: jx,jy,jz ! TODO subroutine this?

! Note: the "default" rms of a unif variable is 0.289
rms_forcing = 0.4_rprec ! TODO do not hard code this

call init_random_seed
do jz=2,nz-1 ! don't force too close to the wall
do jy=1,ny
do jx=1,nx
    call random_number(dummy_rand)
    RHSy(jx,jy,jz) = RHSy(jx,jy,jz) + (rms_forcing/.289_rprec)*(dummy_rand-.5_rprec)
    call random_number(dummy_rand)
    RHSz(jx,jy,jz) = RHSz(jx,jy,jz) + (rms_forcing/.289_rprec)*(dummy_rand-.5_rprec)
end do
end do
end do
if(coord == 0) write(*,*) 'Random forcing added.' ! TODO remove this

end subroutine forcing_random


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
subroutine inflow_cond ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Enforces prescribed inflow condition based on an uniform inflow
!  velocity.
!
use types, only : rprec
use param, only : inflow_velocity, nx, ny, nz
use sim_param, only : u, v, w
use messages, only : error
use fringe_util
implicit none

#ifdef PPVERBOSE
character (*), parameter :: sub_name = 'inflow_cond'
#endif

integer :: i, i_w
integer :: istart, istart_w
integer :: iplateau
integer :: iend, iend_w

real (rprec) :: alpha, beta

!--these may be out of 1, ..., nx
call fringe_init( istart, iplateau, iend )

!--wrapped versions
iend_w = modulo (iend - 1, nx) + 1
istart_w = modulo (istart - 1, nx) + 1

! Set end of domain
u(iend_w, :, :) = inflow_velocity
v(iend_w, :, :) = 0._rprec
w(iend_w, :, :) = 0._rprec

!--skip istart since we know vel at istart, iend already
do i = istart + 1, iend - 1

  i_w = modulo (i - 1, nx) + 1

  beta = fringe_weighting( i, istart, iplateau )
  alpha = 1.0_rprec - beta

  u(i_w, 1:ny, 1:nz) = alpha * u(i_w, 1:ny, 1:nz) + beta * inflow_velocity
  v(i_w, 1:ny, 1:nz) = alpha * v(i_w, 1:ny, 1:nz)
  w(i_w, 1:ny, 1:nz) = alpha * w(i_w, 1:ny, 1:nz)

end do

return
end subroutine inflow_cond


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine project ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! provides u, v, w at 1:nz 
!
use param
use sim_param
use messages
#ifdef PPMPI
  use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#ifdef PPCPS
  use concurrent_precursor, only : synchronize_cps, inflow_cond_cps
#endif
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

! Cases for CPS, Isotropic Turbulence and Uniform inflow
#ifdef PPCPS
call synchronize_cps()
if( inflow ) call inflow_cond_cps()
#elif defined(PPHIT)
if( inflow ) call inflow_HIT()
#else
if ( inflow ) call inflow_cond ()
#endif

!--left this stuff last, so BCs are still enforced, no matter what
!  inflow_cond does

#ifdef PPMPI

! Exchange ghost node information (since coords overlap)                     
call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP )  
  
#endif

!--enfore bc at top
#ifdef PPMPI
if (coord == nproc-1) then
#endif

  ! Note: for ubc_mom > 0, u and v and nz will be written to output as BOGUS
  if (ubc_mom == 0) then    ! no-stress top
     u(:,:,nz) = u(:,:,nz-1)
     v(:,:,nz) = v(:,:,nz-1) 
  endif
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
