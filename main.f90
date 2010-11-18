program main
use types,only:rprec
use param
use sim_param
use grid_defs, only : grid_build
use io, only : openfiles, output_loop, output_final, jt_total, inflow_write, stats_init
use fft
use derivatives, only : filt_da, ddz_uv, ddz_w
use immersedbc
use test_filtermodule
use topbc,only:setsponge,sponge
use bottombc,only:num_patch,avgpatch
use scalars_module,only:beta_scal,obukhov,theta_all_in_one,RHS_T,RHS_Tf
use scalars_module2,only:patch_or_remote

$if ($MPI)
  use mpi_defs, only : initialize_mpi, mpi_sync_real_array, MPI_SYNC_UP
$endif

$if ($LVLSET)
use level_set, only : level_set_init, level_set_cylinder_CD, level_set_smooth_vel
  
  $if ($CYL_SKEW_LS)
  !use cyl_skew_ls, only : cyl_skew_init_ls, cyl_skew_CD_ls
  $endif
  
  $if ($RNS_LS)
  use rns_ls, only : rns_finalize_ls
  
    $if ($CYL_SKEW_LS)
    use rns_cyl_skew_ls, only : rns_init_ls
    $endif
  
  $endif
  
$endif

$if ($TREES_LS)
use trees_ls, only : trees_ls_finalize, trees_ls_init
$endif

$if ($TURBINES)
use turbines, only : turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize, turbines_cond_avg
$endif

$if ($DEBUG)
use debug_mod
$endif

use messages
implicit none

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

character (*), parameter :: sub_name = 'main'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

real(kind=rprec) rmsdivvel,ke, maxcfl
real (rprec):: tt
real (rprec) :: force
real clock_start, clock_end

! Check if read inflow from file is being specified; currently does not work
if(inflow) then
  write(*,*) 'Error: read inflow conditions from file has been specified!'
  write(*,*) 'This capability does not currently work. Please set to false.'
  stop
endif

! Create output directory
call system("mkdir -vp output")

! INITIALIZATION
! Define simulation parameters
call sim_param_init ()

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

$if($MPI)
  if(coord == 0) then
    call cpu_time(clock_start)
    call param_output()
  endif
$else
  call cpu_time(clock_start)
  call param_output()
$endif

    ! Initialize uv grid (calculate x,y,z vectors)
        call grid_build()

    ! Initialize variables used for tavg and other data output
        call stats_init()

    ! Initialize time variable
        tt=0

    ! Determine bottom BC roughness and temperature
    !  only coord 0 needs to do this since at the bottom
    !  roughness information then needs to be broadcast
        call patch_or_remote ()

    ! Initialize turbines
        $if ($TURBINES)
          call turbines_init()    !must occur before initial is called
        $endif

    ! Initialize velocity field
        call initial()

    ! If using level set method
        $if ($LVLSET)
          call level_set_init ()
          
          $if ($CYL_SKEW_LS)
            !call cyl_skew_init_ls ()
          $endif
          
          $if ($RNS_LS)
            call rns_init_ls ()
          $endif
          
        $endif

    ! Initialize fractal trees
        $if ($TREES_LS)
          !--this must come after initial, since fx, fy, fz are set 0 there
          !  and this call may read fx, fy, fz from a file
          call trees_ls_init ()
        $endif

    ! Formulate the fft plans--may want to use FFTW_USE_WISDOM
    ! Initialize the kx,ky arrays
        call init_fft()
    
    ! Open output files (total_time.dat and check_ke.out)  
        call openfiles()

    ! Initialize test filter
    ! this is used for lower BC, even if no dynamic model
        call test_filter_init (2._rprec * filter_size, G_test)

        if (model == 3 .or. model == 5) then  !--scale dependent dynamic
          call test_filter_init (4._rprec * filter_size, G_test_test)
        end if

    ! Initialize sponge variable for top BC, if applicable
        if (ubc == 1) call setsponge()
    

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (u(:, :, 1:nz), 'main.start.u')
  call DEBUG_write (v(:, :, 1:nz), 'main.start.v')
  call DEBUG_write (w(:, :, 1:nz), 'main.start.w')
end if
$endif

$if($CFL_DT)
if( jt_total == 0 .or. abs((cfl_f - cfl)/cfl) > 1.e-2_rprec ) then
  if(.not. USE_MPI .or. ( USE_MPI .and. coord == 0)) write(*,*) '--> Using 1st order Euler for first time step.' 
  call cfl_set_dt(dt) 
  dt = dt * huge(1._rprec) ! Force Euler advection (1st order)
endif
$endif

! BEGIN TIME LOOP
do jt=1,nsteps   
    
    $if($CFL_DT)
      dt_f = dt

      call cfl_set_dt(dt)

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
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz

    ! Compute scalars and __? 
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) call obukhov (jt)

    ! Level set: smooth velocity
    $if ($LVLSET)
        !call level_set_smooth_vel (u, v, w)
    $endif

    ! Buildings: smooth velocity
    if (use_bldg) then      !--no MPI here yet
        call building_interp (u, v, w, .04_rprec, 3)
        call building_interp (dudx, dvdx, dwdx, .04_rprec, 3)
        call building_interp (dudy, dvdy, dwdy, .04_rprec, 3)
    end if

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
        call filt_da (u, dudx, dudy, $lbz)
        call filt_da (v, dvdx, dvdy, $lbz)
        call filt_da (w, dwdx, dwdy, $lbz)
         
        ! Calculate dudz, dvdz using finite differences (for 1:nz on uv-nodes)
        !  except bottom coord, only 2:nz
        call ddz_uv(u, dudz, $lbz)
        call ddz_uv(v, dvdz, $lbz)
       
        ! Calculate dwdz using finite differences (for 0:nz-1 on w-nodes)
        !  except bottom coord, only 1:nz-1
        call ddz_w(w, dwdz, $lbz)

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
        if(use_bldg) call walldudx_building
    end if    

    ! Calculate turbulent (subgrid) stress for entire domain
    !   using the model specified in param.f90 (Smag, LASD, etc)
    !   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
    if (dns_bc .and. molec) then
        call dns_stress(txx,txy,txz,tyy,tyz,tzz)
    else        
        call sgs_stag()
    end if
    if(use_bldg)then
        call wallstress_building(txy,txz,tyz)
        call building_mask(u,v,w)
    endif

    ! Update scalars
    if(S_FLAG.and.(jt.GE.SCAL_INIT))  then
        call theta_all_in_one(jt)
    else
        beta_scal=0._rprec  !buoyancy term set to zero
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
    call convec(RHSx,RHSy,RHSz)

    ! Debug
    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.postconvec.RHSx')
    end if
    $endif

    ! Buildings: set vel to 0. inside buildings
    if (use_bldg) call building_mask (u, v, w)

    ! Add div-tau term to RHS variable
    !   this will be used for pressure calculation
    RHSx(:, :, 1:nz-1) = -RHSx(:, :, 1:nz-1) - divtx(:, :, 1:nz-1)
    RHSy(:, :, 1:nz-1) = -RHSy(:, :, 1:nz-1) - divty(:, :, 1:nz-1)
    RHSz(:, :, 1:nz-1) = -RHSz(:, :, 1:nz-1) - divtz(:, :, 1:nz-1)

    ! Scalars: add buoyancy term to RHS...only valid for theta
    if (S_FLAG .and. (.not.sflux_flag)) then
        RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) + beta_scal(:, :, 1:nz-1)
        if (mod (jt, 1000) == 0) print *, 'Adding buoyancy_term'
        print *,'buoyancy_term=',jt,maxval(beta_scal),minval(beta_scal)
    end if

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

    ! Set RHS*_f if necessary (first timestep)
    if ((jt == 1) .and. (.not. initu)) then
        ! if initu, then this is read from the initialization file
        ! else for the first step put RHS_f=RHS
        !--i.e. at first step, take an Euler step
        RHSx_f=RHSx
        RHSy_f=RHSy
        RHSz_f=RHSz
    end if

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
        call DEBUG_write (fx(:, :, 1:nz), 'main.a.fx')
        call DEBUG_write (force, 'main.a.force')
    end if
    $endif

    ! Calculate intermediate velocity field
    !   only 1:nz-1 are valid
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSx(:, :, 1:nz-1) +  &
                            tadv2 * RHSx_f(:, :, 1:nz-1) + force )
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

    ! Solve Poisson equation for pressure
    !   div of momentum eqn + continuity (div-vel=0) yields Poisson eqn
    !   do not need to store p --> only need gradient
    !   provides p, dpdx, dpdy at 0:nz-1
    !
    ! COMMENTING FOR NOW NEED TO FIX
    ! call press_stag_array (p, dpdx, dpdy)

    ! Calculate dpdz
    !   note: p has additional level at z=-dz/2 for this derivative
    dpdz(1:nx, 1:ny, 1:nz-1) = (p(1:nx, 1:ny, 1:nz-1) -   &
                              p(1:nx, 1:ny, 0:nz-2)) / dz
    dpdz(:, :, nz) = BOGUS

    ! Add pressure gradients to RHS variables
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

    ! Calculate external forces (buildings, trees, turbines, etc)
    !   store in fx,fy,fz arrays
    call forcing ()

    ! Debug
    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.d.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.d.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.d.w')
    end if
    $endif

    ! Projection method provides u,v,w for jz=1:nz
    !   uses fx,fy,fz calculated above
    !   for MPI: syncs nz and 1 node info for u,v,w    
    call project ()
    
    ! Exchange ghost node information (since coords overlap)
    !   send info up (from nz-1 below to 0 above)    
    
    $if ($MPI)
    call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                       u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                       comm, status, ierr)

    call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)

    call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                       w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                       comm, status, ierr)    
    !    call mpi_sync_real_array( u, MPI_SYNC_UP )
    !    call mpi_sync_real_array( v, MPI_SYNC_UP )
    !    call mpi_sync_real_array( w, MPI_SYNC_UP )
    $endif
    
    ! Perform conditional averaging - for turbines
    $if ($TURBINES)
        call turbines_cond_avg()
    $endif       

    ! Write ke to file
    if (modulo (jt, nenergy) == 0) call energy (ke)

    ! Level set: ?
    $if ($LVLSET)
        !call level_set_cylinder_CD ()   
        $if ($CYL_SKEW_LS)
            !call cyl_skew_CD_ls()
        $endif    
        $if ($RNS_LS)
            !call rns_CD_ls()
        $endif
    $endif

    ! Write "jt,dt,rmsdivvel,ke" (and) Coriolis/Scalar info to screen
    if (modulo (jt, wbase) == 0) then
        ! Calculate rms divergence of velocity
        !   only written to screen, not used otherwise
        call rmsdiv (rmsdivvel)
        call cfl_max ( maxcfl )

        if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
          $if($CFL_DT)
          write (6, 7777) jt, dt, rmsdivvel, ke, maxcfl, tadv1, tadv2
          $else
          write (6, 7777) jt, dt, rmsdivvel, ke, maxcfl
          $endif  
          if ((S_FLAG) .or. (coriolis_forcing)) then
            write (6, 7778) wt_s, S_FLAG, patch_flag, remote_flag, &
              coriolis_forcing, ug*u_star
          end if
        end if
    end if
    $if($CFL_DT)
    7777 format ('jt,dt,rmsdivvel,ke,cfl,tadv1,tadv2:',1x,i6.6,4(1x,e9.4),2(1x,f9.4))
    $else
    7777 format ('jt,dt,rmsdivvel,ke,cfl:',1x,i6.6,4(1x,e9.4))
    $endif
    7778 format ('wt_s(K-m/s),Scalars,patch_flag,remote_flag,&
             &coriolis,Ug(m/s):',(f7.3,1x,L2,1x,i2,1x,i2,1x,L2,1x,f7.3))
          
    ! Write output files
        call output_loop (jt)  
        !RNS: Determine if instantaneous plane velocities are to be recorded
        $if ($RNS_LS)
         
        $endif
        
        ! Write inflow_BC file for future use
          if (write_inflow_file) call inflow_write () 

    ! Debug
    $if ($DEBUG)
    if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
    $endif

end do
! END TIME LOOP

! Finalize
    close(2)
    
    ! Write total_time.dat and tavg files
    call output_final (jt)

    ! Level set:
    $if ($LVLSET)
        $if ($TREES_LS)
            call trees_ls_finalize ()
        $endif
        $if ($RNS_LS)
            call rns_finalize_ls ()
        $endif
    $endif

    ! Turbines:
    $if ($TURBINES)
        call turbines_finalize ()   ! must come before MPI finalize
    $endif    

! Stop wall clock
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    call cpu_time (clock_end)
    write(*,"(a,e15.6)") 'Simulation wall time (s) : ', clock_end - clock_start
endif

    
    ! MPI:
    $if ($MPI)
        call mpi_finalize (ierr)
    $endif

end program main
