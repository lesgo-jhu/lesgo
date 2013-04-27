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
implicit none

save

private

public :: forcing_applied, &
          forcing_induced, &
          inflow_cond, &
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

$if ($LVLSET)
$if ($RNS_LS)
use sim_param, only : fxa, fya, fza
use rns_ls, only : rns_forcing_ls
$endif
$endif

$if ($TURBINES)
use sim_param, only : fxa
use turbines, only:turbines_forcing
$endif

implicit none

$if ($LVLSET)
$if ($RNS_LS)
! Reset applied force arrays
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
call rns_forcing_ls()
$endif
$endif

$if ($TURBINES)
! Reset applied force arrays
fxa = 0._rprec
call turbines_forcing ()
$endif
   
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
$if ($LVLSET)
  use level_set, only : level_set_forcing
  use sim_param, only : fx, fy, fz
  $if($RNS_LS)
  use rns_ls, only : rns_elem_force_ls
  $endif
$endif
implicit none

$if($LVLSET)
! Initialize
fx = 0._rprec
fy = 0._rprec
fz = 0._rprec
!  Compute the level set IBM forces
call level_set_forcing ()

$endif

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
use param, only : inflow_velocity, nx, ny, nz, &
                  fringe_region_end, fringe_region_len
use param, only : coord
use sim_param, only : u, v, w, theta
use messages, only : error
use fringe_util
implicit none

character (*), parameter :: sub_name = 'inflow_cond'

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
$if($MPI)
  use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
  $if($CPS)
  use concurrent_precursor, only : synchronize_cps, inflow_cond_cps
  $endif
$endif
implicit none

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: jx, jy, jz
integer :: jz_min

real (rprec) :: RHS, tconst

character(*), parameter :: sub_name='project'

! Caching
tconst = tadv1 * dt

do jz = 1, nz - 1
  do jy = 1, ny
    do jx = 1, nx
 
$if( $LVLSET) 
      RHS = -tadv1 * dpdx(jx, jy, jz)
      u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS + fx(jx, jy, jz)))
      RHS = -tadv1 * dpdy(jx, jy, jz)
      v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS + fy(jx, jy, jz))) 
$else
      RHS = -tadv1 * dpdx(jx, jy, jz)
      u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS                 ))
      RHS = -tadv1 * dpdy(jx, jy, jz)
      v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS                 )) 
$endif

      !if (DEBUG) then
      !  if ( isnan (u(jx, jy, jz)) ) then
      !    write (*, *) $str($context_doc)
      !    write (*, *) 'nan in u at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !  if ( isnan (v(jx, jy, jz)) ) then
      !    write (*, *) $str($context_doc)
      !    write (*, *) 'nan in v at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !end if

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

$if( $LVLSET) 
      RHS = -tadv1 * dpdz(jx, jy, jz)
      w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS + fz(jx, jy, jz)))
$else
      RHS = -tadv1 * dpdz(jx, jy, jz)
      w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS                 ))
$endif
      !if (DEBUG) then
      !  if ( isnan (w(jx, jy, jz)) ) then
      !    write (*, *) 'nan in w at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !end if

    end do
  end do
end do

$if($CPS)
call synchronize_cps()
if( inflow ) call inflow_cond_cps()
$else
if ( inflow ) call inflow_cond ()
$endif

!--left this stuff last, so BCs are still enforced, no matter what
!  inflow_cond does

$if ($MPI)

! Exchange ghost node information (since coords overlap)                     
call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP )  
  
$endif

!--enfore bc at top
$if ($MPI)
if (coord == nproc-1) then
$endif

  if (force_top_bot .and. inflow) then
    u(:, :, nz) = inflow_velocity
    v(:, :, nz) = 0._rprec
  else
    ! no-stress top
    u(:,:,nz)=u(:,:,nz-1)
    ! no-stress top
    v(:,:,nz)=v(:,:,nz-1)
  end if

  w(:, :, nz)=0._rprec

$if ($MPI)
endif
$endif

if (coord == 0) then
  ! just a test
  !if (lbc_mom == 0) then
  !  if (force_top_bot) then
  !    u(:, :, 1) = inflow_velocity
  !    v(:, :, 1) = 0._rprec
  !  else
  !    u(:, :, 1) = u(:, :, 2)
  !    v(:, :, 1) = v(:, :, 2)
  !  end if
  !end if

  w(:, :, 1)=0._rprec

end if

end subroutine project

end module forcing
