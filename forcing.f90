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
use types, only : rprec
implicit none

save

private

public :: forcing_applied, &
          forcing_induced, &
          inflow_cond, &
          project, &
          mode_limit    !!jb , RNL mode limiting

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

$if ($USE_RNL)
use sim_param, only: u, v, w, fxml_rnl, fyml_rnl, fzml_rnl
use param, only: use_ml,jt_total,ml_start
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

$if ($USE_RNL)
if (use_ml .and. jt_total >= ml_start) then
   ! Reset applied force arrays
   fxml_rnl = 0._rprec
   fyml_rnl = 0._rprec
   fzml_rnl = 0._rprec
   call mode_limit(u,v,w)
endif
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
use param, only : inflow_velocity, nx, ny, nz
use sim_param, only : u, v, w
use messages, only : error
use fringe_util
implicit none

$if($VERBOSE)
character (*), parameter :: sub_name = 'inflow_cond'
$endif

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

$if ($VERBOSE)
character(*), parameter :: sub_name='project'
$endif

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

  ! no-stress top
  u(:,:,nz)=u(:,:,nz-1)
  ! no-stress top
  v(:,:,nz)=v(:,:,nz-1)
  ! no permeability
  w(:, :, nz)=0._rprec

$if ($MPI)
endif
$endif

if (coord == 0) then

  ! No modulation of u and v since if a stress free condition (lbc_mom=0) is
  ! applied, it is applied through the momentum equation. 

  ! no permeability
  w(:, :, 1)=0._rprec

end if

end subroutine project

!!$subroutine mode_limit()
!!$use types, only: rprec
!!$use param, only: nx,ny,nz,lbz,jt_total,wbase,L_z,coord,u_star,nproc,kx_allow
!!$use sim_param, only: u,v,w
!!$use grid_defs, only: grid
!!$use fft
!!$
!!$implicit none
!!$
!!$integer :: jx,jy,jz,cutHere
!!$real(rprec) :: const
!!$real(rprec), pointer, dimension(:) :: zw
!!$
!!$nullify(zw)
!!$zw => grid % zw
!!$
!!$cutHere = 2 * kx_allow
!!$
!!$const = 1._rprec / nx
!!$do jy=1,ny
!!$do jz=1,nz-1
!!$
!!$     u(:,jy,jz) = const * u(:,jy,jz)
!!$     v(:,jy,jz) = const * v(:,jy,jz)
!!$     w(:,jy,jz) = const * w(:,jy,jz)
!!$
!!$     call dfftw_execute_dft_r2c(forw_1d, u(:,jy,jz), u(:,jy,jz))
!!$     call dfftw_execute_dft_r2c(forw_1d, v(:,jy,jz), v(:,jy,jz))
!!$     call dfftw_execute_dft_r2c(forw_1d, w(:,jy,jz), w(:,jy,jz))
!!$
!!$     u(3:cutHere,jy,jz) = 0._rprec    !! zero out kx modes 1 to (kx_allow-1)
!!$     v(3:cutHere,jy,jz) = 0._rprec    
!!$     w(3:cutHere,jy,jz) = 0._rprec
!!$
!!$     u(cutHere+3: ,jy,jz) = 0._rprec  !! zero out kx modes of kx_allow+1 to end
!!$     v(cutHere+3: ,jy,jz) = 0._rprec    
!!$     w(cutHere+3: ,jy,jz) = 0._rprec
!!$
!!$     call dfftw_execute_dft_c2r(back_1d, u(:,jy,jz), u(:,jy,jz))
!!$     call dfftw_execute_dft_c2r(back_1d, v(:,jy,jz), v(:,jy,jz))
!!$     call dfftw_execute_dft_c2r(back_1d, w(:,jy,jz), w(:,jy,jz))
!!$
!!$enddo
!!$enddo
!!$
!!$end subroutine mode_limit

subroutine mode_limit(u1,u2,u3)
use types, only: rprec
use param, only: ld, nx, ny, nz, lbz, kx_allow
use param, only: wbase, jt_total, coord
use sim_param, only: fxml_rnl, fyml_rnl, fzml_rnl
use grid_defs, only: grid
use fft

implicit none

real(rprec),dimension(ld,ny,lbz:nz),intent(in) :: u1,u2,u3

integer :: jx,jy,jz,cutHere_p1
real(rprec) :: const

!!$real(rprec), pointer, dimension(:) :: zw
!!$
!!$nullify(zw)
!!$zw => grid % zw

cutHere_p1 = 2 * kx_allow + 1

const = 1._rprec / nx
do jy=1,ny
do jz=1,nz-1

     fxml_rnl(:,jy,jz) = const * u1(:,jy,jz)
     fyml_rnl(:,jy,jz) = const * u2(:,jy,jz)
     fzml_rnl(:,jy,jz) = const * u3(:,jy,jz)

     call dfftw_execute_dft_r2c(forw_1d, fxml_rnl(:,jy,jz), fxml_rnl(:,jy,jz))
     call dfftw_execute_dft_r2c(forw_1d, fyml_rnl(:,jy,jz), fyml_rnl(:,jy,jz))
     call dfftw_execute_dft_r2c(forw_1d, fzml_rnl(:,jy,jz), fzml_rnl(:,jy,jz))

     !! now only zero-out the mode you want to keep
     !! since the forcing will damp/cancel out the others
     fxml_rnl( cutHere_p1 : cutHere_p1 + 1 ,jy,jz) = 0._rprec
     fyml_rnl( cutHere_p1 : cutHere_p1 + 1 ,jy,jz) = 0._rprec    
     fzml_rnl( cutHere_p1 : cutHere_p1 + 1 ,jy,jz) = 0._rprec

     call dfftw_execute_dft_c2r(back_1d,fxml_rnl(:,jy,jz),fxml_rnl(:,jy,jz))
     call dfftw_execute_dft_c2r(back_1d,fyml_rnl(:,jy,jz),fyml_rnl(:,jy,jz))
     call dfftw_execute_dft_c2r(back_1d,fzml_rnl(:,jy,jz),fzml_rnl(:,jy,jz))

enddo
enddo

!! make negative (we want to damp/cancel out these modes)
fxml_rnl = -1._rprec * fxml_rnl
fyml_rnl = -1._rprec * fyml_rnl
fzml_rnl = -1._rprec * fzml_rnl

if (modulo (jt_total,wbase) == 0 .and. coord==0) print*, "Called mode_limit"

end subroutine mode_limit

end module forcing
