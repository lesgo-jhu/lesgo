!!
!!  Copyright (C) 2011-2013  Johns Hopkins University
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

!***********************************************************************
subroutine initialize()
!***********************************************************************
! 
! This subroutine is a driver that calls all top-level initialization
! subroutines.
!
use types, only : rprec
use param, only : coord,dt,jt_total,nsteps,path
use param, only : use_cfl_dt, cfl, cfl_f, dt_dim, z_i, u_star
use param, only : sgs_hist_calc
$if($MPI)
use param, only : MPI_COMM_WORLD, ierr
$endif
$if($USE_RNL)
use param, only: L_x, L_y, kx_num, kx_vec, pi
$endif

use cfl_util
use io, only : output_init
use sgs_param, only : sgs_param_init
use input_util, only : read_input_conf
use test_filtermodule, only : test_filter_init
use sim_param, only : sim_param_init
use grid_defs, only : grid_build
use fft, only : init_fft
use io, only : openfiles
use sgs_hist

$if ($MPI)
use mpi_defs, only : initialize_mpi
  $if($CPS)
  use concurrent_precursor, only : initialize_cps
  $endif
$endif

$if ($LVLSET)
use level_set_base, only : level_set_base_init 
use level_set, only : level_set_init
  $if ($RNS_LS) 
  use rns_base_ls, only : rns_base_init_ls
  $endif
  $if ($RNS_LS and $CYL_SKEW_LS)
  use rns_cyl_skew_ls, only : rns_init_ls
  $endif
$endif

$if ($TURBINES)
use turbines_base, only: turbines_base_init
use turbines, only : turbines_init, turbines_forcing
$endif

$if ($DEBUG)
use debug_mod
$endif

implicit none

character(*), parameter :: make_output_dir = 'mkdir -p ' // path // 'output'
$if ($TURBINES)
character(*), parameter :: make_turbine_dir = 'mkdir -p ' // path // 'turbine'
$endif

! Create output directory
if( coord == 0 ) call system( make_output_dir )
$if ($TURBINES)
if( coord == 0 ) call system( make_turbine_dir )
$endif

! Read input file
! This obtains all major data defined in param
call read_input_conf()

! Initialize MPI
$if ($MPI)
  call initialize_mpi()
$else
  if (nproc /= 1) then
    write (*, *) 'nproc /=1 for non-MPI run is an error'
    stop
  end if
  if (USE_MPI) then
    write (*, *) 'inconsistent use of USE_MPI and $MPI'
    stop
  end if
  chcoord = ''
$endif

! Open output files (total_time.dat and check_ke.out)  
call openfiles()

if( jt_total >= nsteps ) then

   if(coord == 0 ) write(*,'(a)') 'Full number of time steps reached'
   ! MPI:
   $if ($MPI)
   ! First make sure everyone in has finished
   call mpi_barrier( MPI_COMM_WORLD, ierr )
   call mpi_finalize (ierr)
   $endif 
   stop
endif

! Write simulation data to file 
! Commented out since we now have an input file and case information
! can be preserved via it; may still be useful for double checking that
! the input was read correctly and is sane.
if(coord == 0) call param_output()

! Define simulation parameters
call sim_param_init ()
! Initialize sgs variables
call sgs_param_init()

! Initialize uv grid (calculate x,y,z vectors)
call grid_build()

!  Initialize variables used for output statistics and instantaneous data
call output_init()

! Initialize turbines
$if ($TURBINES)
call turbines_base_init()
call turbines_init()    !must occur before initial is called
$endif

! If using level set method
$if ($LVLSET)
call level_set_base_init()
call level_set_init ()

  $if ($RNS_LS)
  call rns_base_init_ls()
  call rns_init_ls ()
  $endif
 
$endif

! Formulate the fft plans--may want to use FFTW_USE_WISDOM
! Initialize the kx,ky arrays
call init_fft()

! Initialize test filter(s)
! this is used for lower BC, even if no dynamic model
call test_filter_init( )

! Initialize velocity field
call initial()
        

! Initialize concurrent precursor stuff
$if($MPI and $CPS)
call initialize_cps()
$endif

! Initialize sgs variable histogram calc
if (sgs_hist_calc) then
  call sgs_hist_init()
endif

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (u(:, :, 1:nz), 'main.start.u')
  call DEBUG_write (v(:, :, 1:nz), 'main.start.v')
  call DEBUG_write (w(:, :, 1:nz), 'main.start.w')
end if
$endif

! Initialize dt if needed to force 1st order Euler
if( use_cfl_dt ) then
   if( jt_total == 0 .or. abs((cfl_f - cfl)/cfl) > 1.e-2_rprec ) then
      if( coord == 0) write(*,*) '--> Using 1st order Euler for first time step.' 
      dt = get_cfl_dt() 
      dt = dt * huge(1._rprec) ! Force Euler advection (1st order)
      dt_dim = dt * z_i / u_star
   endif
endif

return
end subroutine initialize
