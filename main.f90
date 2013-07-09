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
program main
!**********************************************************************
!
! Main file for lesgo solver 
! Contains main time-loop
! 

use types, only : rprec
use clocks 
use param
use sim_param
use grid_defs, only : grid_build
use io, only : energy, output_loop, output_final, jt_total
use fft
use derivatives, only : filt_da, ddz_uv, ddz_w
use test_filtermodule
use cfl_util
use sgs_hist
use sgs_stag_util, only : sgs_stag
use forcing

$if ($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_UP
$endif

$if ($LVLSET)
use level_set, only : level_set_global_CA, level_set_vel_err
use level_set_base, only : global_CA_calc
  
$if ($RNS_LS)
use rns_ls, only : rns_elem_force_ls
$endif
$endif

$if ($TURBINES)
use turbines, only : turbines_forcing, turbine_vel_init
$endif

$if ($DEBUG)
use debug_mod
$endif

use messages

implicit none

character (*), parameter :: prog_name = 'main'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: nstart
real(kind=rprec) rmsdivvel,ke, maxcfl
real (rprec):: tt

type(clock_t) :: clock, clock_total

$if($MPI)
! Buffers used for MPI communication
real(rprec) :: rbuffer
$endif

! Start the clocks, both local and total
call clock_start( clock )

! Initialize time variable
tt = 0
jt = 0
jt_total = 0

! Initialize all data
call initialize()

if(coord == 0) then
   call clock_stop( clock )
   write(*,'(1a,E15.7)') 'Initialization time: ', clock % time
endif

call clock_start( clock_total )

! Initialize starting loop index 
! If new simulation jt_total=0 by definition, if restarting jt_total
! provided by total_time.dat
nstart = jt_total+1

! BEGIN TIME LOOP
time_loop: do jt_total=nstart,nsteps   
  
   ! Get the starting time for the iteration
   call clock_start( clock )

   if( use_cfl_dt ) then
      
      dt_f = dt
      dt = get_cfl_dt()
      dt_dim = dt * z_i / u_star
    
      tadv1 = 1._rprec + 0.5_rprec * dt / dt_f
      tadv2 = 1._rprec - tadv1

   endif

   ! Advance time
   jt = jt + 1
   total_time = total_time + dt
   total_time_dim = total_time_dim + dt_dim
   tt=tt+dt
  
    $if ($DEBUG)
        if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
    $endif

    ! Save previous time's right-hand-sides for Adams-Bashforth Integration
    ! NOTE: RHS does not contain the pressure gradient
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.p.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.p.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.p.w')
    end if
    $endif
  
    ! Calculate velocity derivatives
    ! Calculate dudx, dudy, dvdx, dvdy, dwdx, dwdy (in Fourier space)
    call filt_da (u, dudx, dudy, lbz)
    call filt_da (v, dvdx, dvdy, lbz)
    call filt_da (w, dwdx, dwdy, lbz)
         
    ! Calculate dudz, dvdz using finite differences (for 1:nz on uv-nodes)
    !  except bottom coord, only 2:nz
    call ddz_uv(u, dudz, lbz)
    call ddz_uv(v, dvdz, lbz)
       
    ! Calculate dwdz using finite differences (for 0:nz-1 on w-nodes)
    !  except bottom coord, only 1:nz-1
    call ddz_w(w, dwdz, lbz)

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.q.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.q.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.q.w')
        call DEBUG_write (dudx(:, :, 1:nz), 'main.q.dudx')
        call DEBUG_write (dudy(:, :, 1:nz), 'main.q.dudy')
        call DEBUG_write (dudz(:, :, 1:nz), 'main.q.dudz')
        call DEBUG_write (dvdx(:, :, 1:nz), 'main.q.dvdx')
        call DEBUG_write (dvdy(:, :, 1:nz), 'main.q.dvdy')
        call DEBUG_write (dvdz(:, :, 1:nz), 'main.q.dvdz')
        call DEBUG_write (dwdx(:, :, 1:nz), 'main.q.dwdx')
        call DEBUG_write (dwdy(:, :, 1:nz), 'main.q.dwdy')
        call DEBUG_write (dwdz(:, :, 1:nz), 'main.q.dwdz')
    end if
    $endif

    ! Calculate wall stress and derivatives at the wall (txz, tyz, dudz, dvdz at jz=1)
    !   using the velocity log-law
    !   MPI: bottom process only
    if (dns_bc) then
        if (coord == 0) then
            call wallstress_dns ()
        end if
    else    ! "impose" wall stress 
        if (coord == 0) then
            call wallstress ()                            
        end if
    end if    

    ! Calculate turbulent (subgrid) stress for entire domain
    !   using the model specified in param.f90 (Smag, LASD, etc)
    !   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
    if (dns_bc .and. molec) then
        call dns_stress(txx,txy,txz,tyy,tyz,tzz)
    else        
        call sgs_stag()
    end if

    ! Exchange ghost node information (since coords overlap) for tau_zz
    !   send info up (from nz-1 below to 0 above)
    $if ($MPI)
        call mpi_sendrecv (tzz(:, :, nz-1), ld*ny, MPI_RPREC, up, 6,   &
                           tzz(:, :, 0), ld*ny, MPI_RPREC, down, 6,  &
                           comm, status, ierr)
    $endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (divtx(:, :, 1:nz), 'main.r.divtx')
        call DEBUG_write (divty(:, :, 1:nz), 'main.r.divty')
        call DEBUG_write (divtz(:, :, 1:nz), 'main.r.divtz')
        call DEBUG_write (txx(:, :, 1:nz), 'main.r.txx')
        call DEBUG_write (txy(:, :, 1:nz), 'main.r.txy')
        call DEBUG_write (txz(:, :, 1:nz), 'main.r.txz')
        call DEBUG_write (tyy(:, :, 1:nz), 'main.r.tyy')
        call DEBUG_write (tyz(:, :, 1:nz), 'main.r.tyz')
        call DEBUG_write (tzz(:, :, 1:nz), 'main.r.tzz')
    end if
    $endif    
    
    ! Compute divergence of SGS shear stresses     
    !   the divt's and the diagonal elements of t are not equivalenced in this version
    !   provides divtz 1:nz-1, except 1:nz at top process
    call divstress_uv (divtx, divty, txx, txy, txz, tyy, tyz) ! saves one FFT with previous version
    call divstress_w(divtz, txz, tyz, tzz)

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (divtx(:, :, 1:nz), 'main.s.divtx')
        call DEBUG_write (divty(:, :, 1:nz), 'main.s.divty')
        call DEBUG_write (divtz(:, :, 1:nz), 'main.s.divtz')
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.preconvec.RHSx')
    end if
    $endif

    ! Calculates u x (omega) term in physical space. Uses 3/2 rule for
    ! dealiasing. Stores this term in RHS (right hand side) variable
    call convec()

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.postconvec.RHSx')
    end if
    $endif

    ! Add div-tau term to RHS variable 
    !   this will be used for pressure calculation
    RHSx(:, :, 1:nz-1) = -RHSx(:, :, 1:nz-1) - divtx(:, :, 1:nz-1)
    RHSy(:, :, 1:nz-1) = -RHSy(:, :, 1:nz-1) - divty(:, :, 1:nz-1)
    RHSz(:, :, 1:nz-1) = -RHSz(:, :, 1:nz-1) - divtz(:, :, 1:nz-1)

    ! Coriolis: add forcing to RHS
    if (coriolis_forcing) then
        ! This is to put in the coriolis forcing using coriol,ug and vg as
        ! precribed in param.f90. (ug,vg) specfies the geostrophic wind vector
        ! Note that ug and vg are non-dimensional (using u_star in param.f90)
        RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) +                 &
                             coriol * v(:, :, 1:nz-1) - coriol * vg
        RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) -                 &
                             coriol * u(:, :, 1:nz-1) + coriol * ug
    end if

    !--calculate u^(*) (intermediate vel field)
    !  at this stage, p, dpdx_i are from previous time step
    !  (assumes old dpdx has NOT been added to RHSx_f, etc)
    !  we add force (mean press forcing) here so that u^(*) is as close
    !  to the final velocity as possible
    if (use_mean_p_force) then
        RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) + mean_p_force
    end if

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.a.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.a.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.a.w')
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.a.RHSx')
        call DEBUG_write (RHSy(:, :, 1:nz), 'main.a.RHSy')
        call DEBUG_write (RHSz(:, :, 1:nz), 'main.a.RHSz')
        call DEBUG_write (RHSz_f(:, :, 1:nz), 'main.a.RHSx_f')
        call DEBUG_write (RHSz_f(:, :, 1:nz), 'main.a.RHSy_f')
        call DEBUG_write (RHSz_f(:, :, 1:nz), 'main.a.RHSz_f')
        call DEBUG_write (dpdx(:, :, 1:nz), 'main.a.dpdx')
        call DEBUG_write (dpdy(:, :, 1:nz), 'main.a.dpdy')
        call DEBUG_write (dpdz(:, :, 1:nz), 'main.a.dpdz')
        call DEBUG_write (fxa(:, :, 1:nz), 'main.a.fxa')
        call DEBUG_write (force, 'main.a.force')
    end if
    $endif

    !//////////////////////////////////////////////////////
    !/// APPLIED FORCES                                 ///
    !//////////////////////////////////////////////////////
    !  In order to save memory the arrays fxa, fya, and fza are now only defined when needed.
    !  For Levelset RNS all three arrays are assigned. 
    !  For turbines at the moment only fxa is assigned.
    !  Look in forcing_applied for calculation of forces.
    !  Look in sim_param.f90 for the assignment of the arrays.
        
    !  Applied forcing (forces are added to RHS{x,y,z})
    call forcing_applied()

    !  Update RHS with applied forcing
    $if ($LVLSET)
    $if ($RNS_LS)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)    
    $endif
    $endif

    $if ($TURBINES)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    $endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tony ATM
    $if ($ATM)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)
    $endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tony ATM

    !//////////////////////////////////////////////////////
    !/// EULER INTEGRATION CHECK                        ///
    !////////////////////////////////////////////////////// 
    ! Set RHS*_f if necessary (first timestep) 
    if ((jt_total == 1) .and. (.not. initu)) then
      ! if initu, then this is read from the initialization file
      ! else for the first step put RHS_f=RHS
      !--i.e. at first step, take an Euler step
      RHSx_f=RHSx
      RHSy_f=RHSy
      RHSz_f=RHSz
    end if    

    !//////////////////////////////////////////////////////
    !/// INTERMEDIATE VELOCITY                          ///
    !//////////////////////////////////////////////////////     
    ! Calculate intermediate velocity field
    !   only 1:nz-1 are valid
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSx(:, :, 1:nz-1) +  &
                            tadv2 * RHSx_f(:, :, 1:nz-1) )
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                            tadv2 * RHSy_f(:, :, 1:nz-1) )
    w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                            tadv2 * RHSz_f(:, :, 1:nz-1) )

    ! Set unused values to BOGUS so unintended uses will be noticable
    $if ($SAFETYMODE)
    $if ($MPI)
        u(:, :, 0) = BOGUS
        v(:, :, 0) = BOGUS
        w(:, :, 0) = BOGUS
    $endif
    !--this is an experiment
    !--u, v, w at jz = nz are not useful either, except possibly w(nz), but that
    !  is supposed to zero anyway?
    !--this has to do with what bc are imposed on intermediate velocity    
    u(:, :, nz) = BOGUS
    v(:, :, nz) = BOGUS
    w(:, :, nz) = BOGUS
    $endif
    
    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.b.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.b.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.b.w')
    end if
    $endif

    !//////////////////////////////////////////////////////
    !/// PRESSURE SOLUTION                              ///
    !//////////////////////////////////////////////////////
    ! Solve Poisson equation for pressure
    !   div of momentum eqn + continuity (div-vel=0) yields Poisson eqn
    !   do not need to store p --> only need gradient
    !   provides p, dpdx, dpdy at 0:nz-1 and dpdz at 1:nz-1
    call press_stag_array()

    ! Add pressure gradients to RHS variables (for next time step)
    !   could avoid storing pressure gradients - add directly to RHS
    RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) - dpdx(:, :, 1:nz-1)
    RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) - dpdy(:, :, 1:nz-1)
    RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) - dpdz(:, :, 1:nz-1)

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (dpdx(:, :, 1:nz), 'main.dpdx')
        call DEBUG_write (dpdy(:, :, 1:nz), 'main.dpdy')
        call DEBUG_write (dpdz(:, :, 1:nz), 'main.dpdz')
    end if
    $endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.d.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.d.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.d.w')
    end if
    $endif

    !//////////////////////////////////////////////////////
    !/// INDUCED FORCES                                 ///
    !//////////////////////////////////////////////////////    
    ! Calculate external forces induced forces. These are
    ! stored in fx,fy,fz arrays. We are calling induced 
    ! forces before applied forces as some of the applied
    ! forces (RNS) depend on the induced forces and the 
    ! two are assumed independent
    call forcing_induced()

    !//////////////////////////////////////////////////////
    !/// PROJECTION STEP                                ///
    !//////////////////////////////////////////////////////   
    ! Projection method provides u,v,w for jz=1:nz
    !   uses fx,fy,fz calculated above
    !   for MPI: syncs 1 -> Nz and Nz-1 -> 0 nodes info for u,v,w    
    call project ()

    $if($LVLSET and $RNS_LS)
    !  Compute the relavent force information ( include reference quantities, CD, etc.)
    !  of the RNS elements using the IBM force; No modification to f{x,y,z} is
    !  made here.
    call rns_elem_force_ls()
    $endif
   
    ! Write ke to file
    if (modulo (jt_total, nenergy) == 0) call energy (ke)

    $if ($LVLSET)
      if( global_CA_calc ) call level_set_global_CA ()
    $endif

    ! Write output files
    call output_loop()  

    ! Check the total time of the simulation up to this point on the master node and send this to all
   
    if (modulo (jt_total, wbase) == 0) then
       
       ! Get the ending time for the iteration
       call clock_stop( clock )
       call clock_stop( clock_total )

       ! Calculate rms divergence of velocity
       ! only written to screen, not used otherwise
       call rmsdiv (rmsdivvel)
       maxcfl = get_max_cfl()

       if (coord == 0) then
          write(*,*)
          write(*,'(a)') '========================================'
          write(*,'(a)') 'Time step information:'
          write(*,'(a,i9)') '  Iteration: ', jt_total
          write(*,'(a,E15.7)') '  Time step: ', dt
          write(*,'(a,E15.7)') '  CFL: ', maxcfl
          write(*,'(a,2E15.7)') '  AB2 TADV1, TADV2: ', tadv1, tadv2
          write(*,*) 
          write(*,'(a)') 'Flow field information:'          
          write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
          write(*,'(a,E15.7)') '  Kinetic energy: ', ke
          write(*,*)
          write(*,'(1a)') 'Simulation wall times (s): '
          write(*,'(1a,E15.7)') '  Iteration: ', clock % time
          write(*,'(1a,E15.7)') '  Cumulative: ', clock_total % time
          write(*,'(a)') '========================================'
       end if

       ! Check if we are to check the allowable runtime
       if( runtime > 0 ) then

          ! Determine the processor that has used most time and communicate this.
          ! Needed to make sure that all processors stop at the same time and not just some of them
          $if($MPI)
          call mpi_allreduce(clock_total % time, rbuffer, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
          clock_total % time = rbuffer
          $endif
       
          ! If maximum time is surpassed go to the end of the program
          if ( clock_total % time >= real(runtime,rprec) ) then
             call mesg( prog_name, 'Specified runtime exceeded. Exiting simulation.')
             exit time_loop
          endif

       endif

    end if

end do time_loop
! END TIME LOOP

! Finalize
close(2)
    
! Write total_time.dat and tavg files
call output_final()

! Stop wall clock
call clock_stop( clock_total )
if( coord == 0 )  write(*,"(a,e15.7)") 'Simulation wall time (s) : ', clock_total % time

call finalize()

if(coord == 0 ) write(*,'(a)') 'Simulation complete'

end program main
