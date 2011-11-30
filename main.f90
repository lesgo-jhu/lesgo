!**********************************************************************
program main
!**********************************************************************
use types, only : rprec
use clocks 
use param
use sim_param
use grid_defs, only : grid_build
use io, only : output_loop, output_final, jt_total
use fft
use derivatives, only : filt_da, ddz_uv, ddz_w
use test_filtermodule
use cfl_mod

use sgs_stag_util, only : sgs_stag

$if ($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_UP
$endif

$if ($LVLSET)
use level_set, only : level_set_global_CD, level_set_smooth_vel, level_set_vel_err
use level_set_base, only : global_CD_calc
  
  $if ($RNS_LS)
  use rns_ls, only : rns_finalize_ls, rns_elem_force_ls
  $endif

$endif

$if ($TURBINES)
use turbines, only : turbines_forcing, turbine_vel_init, turbines_finalize, turbines_cond_avg
$endif

$if ($DEBUG)
use debug_mod
$endif

use messages
implicit none

character (*), parameter :: sub_name = 'main'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

real(kind=rprec) rmsdivvel,ke, maxcfl
real (rprec):: tt
real (rprec) :: force

type(clock_type) :: clock_t, clock_total_t


! Start the clocks, both local and total
call clock_start( clock_t )
clock_total_t = clock_t

! Initialize time variable
tt = 0

! Initialize all data
call initialize()


$if($MPI)
  if(coord == 0) then
     call clock_stop( clock_t )
     write(*,'(1a,E15.7)') 'Initialization time: ', clock_time( clock_t ) 
  endif
$else
  call clock_stop( clock_t )
  write(*,'(1a,E15.7)') 'Initialization time: ', clock_time( clock_t ) 
$endif

! BEGIN TIME LOOP
do jt=1,nsteps   
   
   ! Get the starting time for the iteration
   call clock_start( clock_t )

    $if($CFL_DT)
      dt_f = dt

      dt = get_cfl_dt()

      dt_dim = dt * z_i / u_star
    
      tadv1 = 1._rprec + 0.5_rprec * dt / dt_f
      tadv2 = 1._rprec - tadv1

    $endif

    ! Advance time
    jt_total = jt_total + 1 
    total_time = total_time + dt
    total_time_dim = total_time_dim + dt_dim
    tt=tt+dt
  
    ! Debug
    $if ($DEBUG)
        if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
    $endif

    ! Save previous time's right-hand-sides for Adams-Bashforth Integration
    ! NOTE: RHS does not contain the pressure gradient
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz

    ! Level set: smooth velocity
    $if ($LVLSET_SMOOTH_VEL)
      call level_set_smooth_vel (u, v, w)
    $endif

    ! Debug
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

    ! Debug
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
        if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
            call wallstress_dns ()
        end if
    else    ! "impose" wall stress 
        if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
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

    ! Debug
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
    call divstress_uv(divtx, txx, txy, txz)
    call divstress_uv(divty, txy, tyy, tyz)    
    call divstress_w(divtz, txz, tyz, tzz)

    ! Debug
    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (divtx(:, :, 1:nz), 'main.s.divtx')
        call DEBUG_write (divty(:, :, 1:nz), 'main.s.divty')
        call DEBUG_write (divtz(:, :, 1:nz), 'main.s.divtz')
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.preconvec.RHSx')
    end if
    $endif

    ! Calculates u x (omega) term in physical space
    !   uses 3/2 rule for dealiasing
    !   stores this term in RHS (right hand side) variable
    !call convec(RHSx,RHSy,RHSz)
    call convec()

    ! Debug
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
        force = mean_p_force
    else
        force = 0._rprec
    end if
  
    RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) + force

    ! Debug
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
    !  Applied forcing (forces are added to RHS{x,y,z})
    call forcing_applied()

    !  Update RHS with applied forcing
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)    

    !//////////////////////////////////////////////////////
    !/// EULER INTEGRATION CHECK                        ///
    !////////////////////////////////////////////////////// 
    ! Set RHS*_f if necessary (first timestep) 
    if ((jt == 1) .and. (.not. initu)) then
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

    ! Debug
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

    ! Debug
    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (dpdx(:, :, 1:nz), 'main.dpdx')
        call DEBUG_write (dpdy(:, :, 1:nz), 'main.dpdy')
        call DEBUG_write (dpdz(:, :, 1:nz), 'main.dpdz')
    end if
    $endif

    ! Debug
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
   
    ! Perform conditional averaging - for turbines
    $if ($TURBINES)
        call turbines_cond_avg()
    $endif       

    ! Write ke to file
    if (modulo (jt, nenergy) == 0) call energy (ke)

    $if ($LVLSET)
      if( global_CD_calc ) call level_set_global_CD ()
    $endif

    ! Write output files
    call output_loop (jt)  
    !RNS: Determine if instantaneous plane velocities are to be recorded
        
    ! Write "jt,dt,rmsdivvel,ke" (and) Coriolis/Scalar info to screen
    if (modulo (jt, wbase) == 0) then
       
       ! Get the ending time for the iteration
       call clock_stop( clock_t )
       call clock_stop( clock_total_t )

       
        ! Calculate rms divergence of velocity
       !   only written to screen, not used otherwise
       call rmsdiv (rmsdivvel)
       maxcfl = get_max_cfl()

       
       if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
          write(*,*)
          write(*,'(a)') '========================================'
          write(*,'(a)') 'Time step information:'
          write(*,'(a,i9)') '  Iteration: ', jt
          write(*,'(a,E15.7)') '  Time step: ', dt
          $if($CFL_DT)
          write(*,'(a,E15.7)') '  CFL: ', maxcfl
          $endif
          write(*,'(a,2E15.7)') '  AB2 TADV1, TADV2: ', tadv1, tadv2
          write(*,*) 
          write(*,'(a)') 'Flow field information:'          
          write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
          write(*,'(a,E15.7)') '  Kinetic energy: ', ke
          ! $if($CFL_DT)
          ! write (6, 7777) jt, dt, rmsdivvel, ke, maxcfl, tadv1, tadv2
          ! $else
          ! write (6, 7777) jt, dt, rmsdivvel, ke, maxcfl
          ! $endif  
          write(*,*)
          write(*,'(1a)') 'Simulation wall times (s): '
          write(*,'(1a,E15.7)') '  Iteration: ', clock_time( clock_t )
          write(*,'(1a,E15.7)') '  Cumulative: ', clock_time( clock_total_t )
          write(*,'(a)') '========================================'
       end if

    end if

    ! $if($CFL_DT)
    ! 7777 format ('jt,dt,rmsdivvel,ke,cfl,tadv1,tadv2:',1x,i6.6,4(1x,e9.4),2(1x,f9.4))
    ! $else
    ! 7777 format ('jt,dt,rmsdivvel,ke,cfl:',1x,i6.6,4(1x,e9.4))
    ! $endif
    ! 7778 format ('wt_s(K-m/s),Scalars,patch_flag,remote_flag,&
    !      &coriolis,Ug(m/s):',(f7.3,1x,L2,1x,i2,1x,i2,1x,L2,1x,f7.3))

end do
! END TIME LOOP

! Finalize
close(2)
    
! Write total_time.dat and tavg files
call output_final (jt)

! Level set:
$if ($LVLSET)

  $if ($RNS_LS)
  call rns_finalize_ls ()
  $endif
$endif

! Turbines:
$if ($TURBINES)
call turbines_finalize ()   ! must come before MPI finalize
$endif    

! Stop wall clock
call clock_stop( clock_total_t )
$if($MPI)
  if( coord == 0 )  write(*,"(a,e15.7)") 'Simulation wall time (s) : ', clock_time( clock_total_t )
$else
  write(*,"(a,e15.7)") 'Simulation wall time (s) : ', clock_time( clock_total_t )
$endif

! MPI:
$if ($MPI)
call mpi_finalize (ierr)
$endif

end program main
