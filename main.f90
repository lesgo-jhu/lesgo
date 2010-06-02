program main
use types,only:rprec
use param
use sim_param
use grid_defs, only : grid_build
!!$use io, only : output_loop, output_final, inflow_write,  &
!!$            avg_stats
use io, only : openfiles, output_loop, output_final, jt_total, inflow_write, stats_init
use fft
use immersedbc
use test_filtermodule
use topbc,only:setsponge,sponge
use bottombc,only:num_patch,avgpatch
use scalars_module,only:beta_scal,obukhov,theta_all_in_one,RHS_T,RHS_Tf
use scalars_module2,only:patch_or_remote

$if ($MPI)
  use mpi_defs, only : initialize_mpi
$endif

$if ($LVLSET)
use level_set, only : level_set_init, level_set_cylinder_CD, level_set_smooth_vel
  
  $if ($CYL_SKEW_LS)
  use cyl_skew_ls, only : cyl_skew_init_ls, cyl_skew_CD_ls
  $endif
  
  $if ($RNS_LS)
  use rns_ls, only : rns_init_ls, rns_CD_ls, rns_finalize_ls
  $endif
  
$endif

$if ($TREES_LS)
use trees_ls, only : trees_ls_finalize, trees_ls_init
$endif

$if ($TURBINES)
use turbines, only : turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize
$endif

$if ($DEBUG)
use debug_mod  !--just for debugging
$endif

use messages
implicit none

integer,parameter::wbase=100  !--controls the frequency of screen diagnostics
integer, parameter :: nenergy = 1  !--frequency of writes to check_ke.dat

character (*), parameter :: sub_name = 'main'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

!--for trees_CV, we need the SGS stress, so it was moved to sim_params
!  also, the equivalence was removed, since this overwrites what we need
!  also for same reason cx,cy,cz put in sim_params
!real(kind=rprec),dimension(ld,ny,nz)::divtx,divty,divtz
!real(kind=rprec), dimension(ld,ny,nz)::cx,cy,cz,&
!     divtx,divty,divtz
!real(kind=rprec),dimension(ld,ny,nz)::cx,cy,cz,&
!     txx,txy,txz,tyy,tyz,tzz,divtx,divty,divtz
!equivalence (txx,divtx),(tyy,divty),(tzz,divtz)
real(kind=rprec) rmsdivvel,ke
real (rprec):: tt
real (rprec) :: force
real clock_start, clock_end

!  Start wall clock
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  call cpu_time (clock_start)
endif

!---------------------------------------------------------------------
!  Check if read inflow from file is being specified; it currently does
!  not work
if(inflow) then
  write(*,*) 'Error: read inflow conditions from file has been specified!'
  write(*,*) 'This capability does not currently work. Please set to false.'
  stop
endif

!  Create output directory
call system("mkdir -vp output")

call sim_param_init ()

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

  !--leave this blank or put in coord
  !write (chcoord, '(a,i0,a)') '(', coord, ')'  !--() make easier to use
  chcoord = ''

$endif

!  Initialize uv grid
call grid_build()
!  Initialized statics arrays
call stats_init()

if(coord == 0) then
  write(*,*) 'dz = ', dz
  write(*,*) 'nz = ', nz
  write(*,*) 'nz_tot = ', nz_tot
endif

tt=0

!--only coord 0 needs to do this since at the bottom
!--roughness information then needs to be broadcast
!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  call patch_or_remote ()
!end if


$if ($TURBINES)
  call turbines_init()    !must occur before initial is called
$endif

call initial()
!--could move this into something like initial ()
$if ($LVLSET)

  call level_set_init ()
  
  $if ($CYL_SKEW_LS)
  call cyl_skew_init_ls ()
  $endif
  
  $if ($RNS_LS)
  call rns_init_ls ()
  $endif
  
$endif


$if ($TREES_LS)
  !--this must come after initial, since fx, fy, fz are set 0 there
  !  and this call may read fx, fy, fz from a file
  call trees_ls_init ()
$endif

! formulate the fft plans--may want to use FFTW_USE_WISDOM
! initialize the kx,ky arrays
call init_fft()
! Open output files      
call openfiles()

!--initialize test filter
!--this is used for lower BC, even if no dynamic model
call test_filter_init (2._rprec * filter_size, G_test)

if (model == 3 .or. model == 5) then  !--scale dependent dynamic

  call test_filter_init (4._rprec * filter_size, G_test_test)
end if

if (ubc == 1) call setsponge()

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  print *, 'Number of timesteps', nsteps
  print *, 'dt = ', dt
  if (model == 1) print *, 'Co = ', Co
  print *, 'Nx, Ny, Nz = ', nx, ny, nz
  print *, 'Lx, Ly, Lz = ', L_x, L_y, L_z
  if (USE_MPI) print *, 'Number of processes = ', nproc
  print *, 'Number of patches = ', num_patch
  !print *, 'sampling stats every ', c_count, ' timesteps'
  !print *, 'writing stats every ', p_count, ' timesteps'
  if (molec) print*, 'molecular viscosity (dimensional) ', nu_molec
end if

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (u(:, :, 1:nz), 'main.start.u')
  call DEBUG_write (v(:, :, 1:nz), 'main.start.v')
  call DEBUG_write (w(:, :, 1:nz), 'main.start.w')
end if
$endif

!--MPI: u,v,w should be set at jz = 0:nz before getting here, except
!  bottom process which is BOGUS (starts at 1)
! time Loop

do jt=1,nsteps   

  jt_total = jt_total + 1  !--moved from io.f90
  total_time = total_time + dt
  total_time_dim = total_time_dim + dt_dim
  
  $if ($DEBUG)
  if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
  $endif

  ! advance total time
  tt=tt+dt

  ! save previous time's right-hand-sides for Adams-Bashforth Integration
  !  (In subroutine "STEP" use first order time advancement on first time 
  !  step).
  !if (jt > 1) then    
     RHSx_f = RHSx
     RHSy_f = RHSy
     RHSz_f = RHSz
  !end if

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) call obukhov (jt)

  $if ($LVLSET)
    !call level_set_smooth_vel (u, v, w)
  $endif

  !--no MPI here yet
  if (use_bldg) then
    call building_interp (u, v, w, .04_rprec, 3)
    call building_interp (dudx, dvdx, dwdx, .04_rprec, 3)
    call building_interp (dudy, dvdy, dwdy, .04_rprec, 3)
  end if

  $if ($DEBUG)
  if (DEBUG) then
    call DEBUG_write (u(:, :, 1:nz), 'main.p.u')
    call DEBUG_write (v(:, :, 1:nz), 'main.p.v')
    call DEBUG_write (w(:, :, 1:nz), 'main.p.w')
  end if
  $endif

  ! kill oddballs and calculate horizontal derivatives
  !--MPI: all these arrays are 0:nz
  !--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
  !  except on bottom process (0 level set to BOGUS, starts at 1)
  call filt_da (u, dudx, dudy)
  call filt_da (v, dvdx, dvdy)
  call filt_da (w, dwdx, dwdy)
   ! finite differences
   !--MPI: on exit of ddz_uv, have dudz, dvdz at 1:nz, except
   !  bottom process has 2:nz
   call ddz_uv(dudz,u)
   call ddz_uv(dvdz,v)
   !--MPI: on exit of ddz_w, have dwdz at 0:nz-1, except top process
   !  has 0:nz, and bottom process has 1:nz-1
   call ddz_w(dwdz,w)

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

!TS calculate wall stress and calculate derivatives at wall
   if (dns_bc) then
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       call wallstress_dns ()
     end if
   else
!TS "impose" wall stress and calculate derivatives at wall
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       call wallstress ()  !--provides txz, tyz, dudz, dvdz at jz=1
                           !--MPI: bottom process only
     end if
     if(use_bldg) call walldudx_building
   end if

! compute turbulent viscosity (const.)
  if (dns_bc .and. molec) then
    call dns_stress(txx,txy,txz,tyy,tyz,tzz)
  else
    !--MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz
    call sgs_stag()	    
  end if
  if(use_bldg)then
     call wallstress_building(txy,txz,tyz)
     call building_mask(u,v,w)
  endif

!xx----VK -ADDED FOR SCALARS !! --xxxxxx
!TS ADD sflux_flag
  if(S_FLAG.and.(jt.GE.SCAL_INIT))  then
  call theta_all_in_one(jt)
  else
  beta_scal=0._rprec
!  print *,'buoyancy term set to zero !!'
  end if
!xx----VK -ADDED FOR SCALARS !! --xxxxxx

  $if ($MPI)
     !--exchange ghost-node information for tij
     !--send stuff up to ghost nodes
     !--move this into sgs_stag?
     call mpi_sendrecv (tzz(:, :, nz-1), ld*ny, MPI_RPREC, up, 6,   &
                        tzz(:, :, 0), ld*ny, MPI_RPREC, down, 6,  &
                        comm, status, ierr)
  $endif

! compute divergence of SGS shear stresses     
! note: the divt's and the diagonal elements of t are equivalenced!
!--actually, they are not equivalenced in this version
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

  !--provides divtz 1:nz-1
  call divstress_uv(divtx, txx, txy, txz)
  call divstress_uv(divty, txy, tyy, tyz)
  !--provides divtz 1:nz-1, except 1:nz at top process
  call divstress_w(divtz, txz, tyz, tzz)

 ! --------- NOTE ----------
 !  When using the ifort compiler for a 
 !  grid resolution > 32^3, during execution
 !  a segmentation fault is recieved somewhere 
 !  just be for this section (JSG - 1/23/09)
 !
  $if ($DEBUG)
  if (DEBUG) then
    call DEBUG_write (divtx(:, :, 1:nz), 'main.s.divtx')
    call DEBUG_write (divty(:, :, 1:nz), 'main.s.divty')
    call DEBUG_write (divtz(:, :, 1:nz), 'main.s.divtz')
    call DEBUG_write (RHSx(:, :, 1:nz), 'main.preconvec.RHSx')
  end if
  $endif

  !--provides RHS{x,y,z} 1:nz-1
  call convec(RHSx,RHSy,RHSz)

  $if ($DEBUG)
  if (DEBUG) then
    call DEBUG_write (RHSx(:, :, 1:nz), 'main.postconvec.RHSx')
  end if
  $endif

  if (use_bldg) call building_mask (u, v, w)

! Compute preliminary RHS matrices for pressure calculation
  RHSx(:, :, 1:nz-1) = -RHSx(:, :, 1:nz-1) - divtx(:, :, 1:nz-1)
  RHSy(:, :, 1:nz-1) = -RHSy(:, :, 1:nz-1) - divty(:, :, 1:nz-1)
  RHSz(:, :, 1:nz-1) = -RHSz(:, :, 1:nz-1) - divtz(:, :, 1:nz-1)

  if (S_FLAG .and. (.not.sflux_flag)) then
    !--add buoyancy term...only valid for theta
    RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) + beta_scal(:, :, 1:nz-1)
    if (mod (jt, 1000) == 0) print *, 'Adding buoyancy_term'
    !print *,'buoyancy_term=',jt,maxval(beta_scal),minval(beta_scal)
  end if

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

  if ((jt == 1) .and. (.not. initu)) then
    ! if initu, then this is read from the initialization file
    ! else for the first step put RHS_f=RHS
    !--i.e. at first step, take an Euler step
    RHSx_f=RHSx
    RHSy_f=RHSy
    RHSz_f=RHSz
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
    call DEBUG_write (fx(:, :, 1:nz), 'main.a.fx')
    call DEBUG_write (force, 'main.a.force')
  end if
  $endif

   !--only 1:nz-1 are valid
   u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                           &
                     dt * ( tadv1 * RHSx(:, :, 1:nz-1) +         &
                            tadv2 * RHSx_f(:, :, 1:nz-1) + force )
   v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                    &
                     dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                            tadv2 * RHSy_f(:, :, 1:nz-1) )
   w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                    &
                     dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                            tadv2 * RHSz_f(:, :, 1:nz-1) )

  $if ($MPI)
    !--after this point, u,v,w at jz = 0 are not useful, until updated
    u(:, :, 0) = BOGUS
    v(:, :, 0) = BOGUS
    w(:, :, 0) = BOGUS
  $endif

  !--this is an experiment
  u(:, :, nz) = BOGUS
  v(:, :, nz) = BOGUS
  w(:, :, nz) = BOGUS

  !--this is an experiment
  !u(:, :, nz) = u(:, :, nz-1)  !BOGUS
  !v(:, :, nz) = v(:, :, nz-1)  !BOGUS
  !w(:, :, nz) = 0._rprec  !BOGUS

  !--u, v, w at jz = nz are not useful either, except possibly w(nz), but that
  !  is supposed to zero anyway?
  !--this has to do with what bc are imposed on intermediate velocity

  $if ($DEBUG)
  if (DEBUG) then
    call DEBUG_write (u(:, :, 1:nz), 'main.b.u')
    call DEBUG_write (v(:, :, 1:nz), 'main.b.v')
    call DEBUG_write (w(:, :, 1:nz), 'main.b.w')
  end if
  $endif

  !--solve Poisson equation for pressure
  !--do we ever need p itself, or only its gradient? -> probably
  !  do not need to store p
  !call press_stag (p, dpdx, dpdy)
  !--provides p, dpdx, dpdy at 0:nz-1
  call press_stag_array (p, dpdx, dpdy)

  !--calculate dpdz here
  !--careful, p is not dimensioned the same as the others
  dpdz(1:nx, 1:ny, 1:nz-1) = (p(1:nx, 1:ny, 1:nz-1) -   &
                              p(1:nx, 1:ny, 0:nz-2)) / dz
  dpdz(:, :, nz) = BOGUS

  !--if really wanted to, could avoid storing pressure gradients
  !  just add them directly to RHS in press_stag
  RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) - dpdx(:, :, 1:nz-1)
  RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) - dpdy(:, :, 1:nz-1)
  RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) - dpdz(:, :, 1:nz-1)

  $if ($DEBUG)
  if (DEBUG) then
    !--note: only 1:nz-1 valid here
    call DEBUG_write (dpdx(:, :, 1:nz), 'main.dpdx')
    call DEBUG_write (dpdy(:, :, 1:nz), 'main.dpdy')
    call DEBUG_write (dpdz(:, :, 1:nz), 'main.dpdz')
  end if
  $endif

  call forcing ()

  $if ($DEBUG)
  if (DEBUG) then
    call DEBUG_write (u(:, :, 1:nz), 'main.d.u')
    call DEBUG_write (v(:, :, 1:nz), 'main.d.v')
    call DEBUG_write (w(:, :, 1:nz), 'main.d.w')
  end if
  $endif

  !--provides u, v, w at 1:nz 
  call project ()

  $if ($MPI)
    !--exchange ghost-node information
    !--send stuff up to ghost layer (0) (for z-derivs)
    !--nz should already be in sync with 1 level: done in project()
    call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                       u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                       comm, status, ierr)

    call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)

    call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                       w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                       comm, status, ierr)
  $endif

  !--MPI: at this point, have u, v, w at 0:nz

  if (modulo (jt, nenergy) == 0) call energy (ke)

!!$  call avg_stats ()  !--only does something once every n_avg_stats steps

  $if ($LVLSET)
    ! call level_set_cylinder_CD ()
    
    $if ($CYL_SKEW_LS)
    !  call cyl_skew_CD_ls()
    $endif
    
    $if ($RNS_LS)
    call rns_CD_ls()
    $endif

  $endif

  if (modulo (jt, wbase) == 0) then

    call rmsdiv (rmsdivvel)

    if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then

      write (6, 7777) jt, dt, rmsdivvel, ke

      if ((S_FLAG) .or. (coriolis_forcing)) then
        write (6, 7778) wt_s, S_FLAG, patch_flag, remote_flag, &
                        coriolis_forcing, ug*u_star
      end if
    end if

  end if
   
7777 format ('jt,dt,rmsdivvel,ke:',1x,i6.6,3(1x,e9.4))
7778 format ('wt_s(K-m/s),Scalars,patch_flag,remote_flag,&
             &coriolis,Ug(m/s):',(f7.3,1x,L2,1x,i2,1x,i2,1x,L2,1x,f7.3))
          
  call output_loop (jt)
  
!  $if ($RNS_LS)
!!  Determine if instantaneous plane velocities are to be recorded
!  if(rns_t%plane_u_calc) call rns_u_write_ls()
!  $endif

  if (write_inflow_file) call inflow_write () !--for creating inflow_BC file

  $if ($DEBUG)
  if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
  $endif

end do  !--end time loop
close(2)
call output_final (jt)

$if ($LVLSET)

$if ($TREES_LS)
  call trees_ls_finalize ()
$endif

$if($RNS_LS)
  call rns_finalize_ls ()
$endif 
 
$endif

$if ($MPI)
  call mpi_finalize (ierr)
$endif

$if ($LVLSET)
$if ($RNS_LS)
  call rns_finalize_ls ()
$endif
$endif

$if ($TURBINES)
  call turbines_finalize ()
$endif

!  Stop wall clock
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  call cpu_time (clock_end)
  write(*,"(a,e)") 'Simulation wall time (s) : ', clock_end - clock_start
endif

stop

end program main
