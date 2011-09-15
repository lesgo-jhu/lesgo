module param
  use types,only:rprec,point3D
  $if ($MPI)
  use mpi
  $endif
  implicit none

  save

  private rprec  !--this is dumb.
  public
  
!---------------------------------------------------
! MPI PARAMETERS
!---------------------------------------------------

  $if ($MPI)
  $define $MPI_LOGICAL .true.
  $define $NPROC 4
  $else
  $define $MPI_LOGICAL .false.
  $define $NPROC 1
  $endif

  logical, parameter :: USE_MPI = $MPI_LOGICAL

  $undefine $MPI_LOGICAL

  $if ($MPI)
  integer :: status(MPI_STATUS_SIZE)
  $endif

  character (*), parameter :: path = './'

  !--this stuff must be defined, even if not using MPI
  character (8) :: chcoord  !--holds character representation of coord
  integer, parameter :: nproc = $NPROC  !--this must be 1 if no MPI
  integer :: ierr
  integer :: comm
  integer :: up, down
  integer :: global_rank
  integer :: MPI_RPREC, MPI_CPREC
  integer :: rank = -1   !--init to bogus (so its defined, even if no MPI)
  integer :: coord = -1  !--same here
  integer :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)
  !--end mpi stuff
  
!---------------------------------------------------
! COMPUTATIONAL DOMAIN PARAMETERS
!---------------------------------------------------  
! characteristic length is H=z_i and characteristic velocity is u_star  
!   L_x, L_y, L_z, dx, dy, dz are non-dim. using H

  integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
  real (rprec), parameter :: BOGUS = -1234567890._rprec
  real(rprec),parameter::pi=3.1415926535897932384626433_rprec

  integer,parameter:: nx=64,ny=64,nz=(64)/nproc + 1 
  
  integer, parameter :: nz_tot = (nz - 1) * nproc + 1
  integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
  integer,parameter:: lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

  ! this value is dimensional [m]:
  real(rprec),parameter::z_i=1000._rprec   !dimensions in meters, height of BL
    
  ! these values should be non-dimensionalized by z_i: 
  ! set as multiple of BL height (z_i) then non-dimensionalized by z_i
  real(rprec),parameter::L_x= 2*pi
  real(rprec),parameter::L_y= L_x
  real(rprec),parameter::L_z= 1.0_rprec
  !real(rprec),parameter::L_y=(1.*ny)/(1.*nx)*L_x               ! ensure dy=dx
  !real(rprec),parameter::L_z=(nz_tot - 1.)/nx*L_x  ! ensure dz = dx

  ! these values are also non-dimensionalized by z_i:
  real(rprec),parameter::dz=L_z/(nz_tot-1.) ! or (L_z/nproc)/(nz - 1)
  real(rprec),parameter::dx=L_x/nx,dy=L_y/ny
  
!---------------------------------------------------
! MODEL PARAMETERS
!---------------------------------------------------   

  ! Model type: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent
  !             4->Lagrangian scale-sim   5-> Lagragian scale-dep
  ! Models type: 1->static prandtl, 2->Dynamic
  integer,parameter::model=5,models=1,nnn=2

  ! timesteps between dynamic Cs updates           
  integer, parameter :: cs_count = 5
  ! When to start dynamic Cs calculations
  integer,parameter::DYN_init=100
  
  ! Cs is the Smagorinsky Constant
  ! Co and nnn are used in the mason model for smagorisky coeff
  real(kind=rprec),parameter::Co=0.16_rprec
  
  ! test filter type: 1->cut off 2->Gaussian 3->Top-hat
  integer,parameter::ifilter=1

  ! u_star=0.45 m/s if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
  real(rprec),parameter::u_star=0.45_rprec,Pr=.4_rprec

  ! von Karman constant     
  real(rprec),parameter::vonk=0.4_rprec   
  
  ! Coriolis stuff
  ! coriol=non-dim coriolis parameter,
  ! ug=horiz geostrophic vel, vg=transverse geostrophic vel
  logical,parameter::coriolis_forcing=.false.
  real(rprec),parameter::coriol=9.125E-05*z_i/u_star,      &
       ug=u_star/u_star,vg=0._rprec/u_star

  ! nu_molec is dimensional m^2/s
  real(rprec),parameter::nu_molec=1.14e-5_rprec   
    
  logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.  
  
!---------------------------------------------------
! TIMESTEP PARAMETERS
!---------------------------------------------------   

  integer, parameter :: nsteps = 150000
 
  $if($CFL_DT)
  
  real(rprec), parameter :: cfl = 0.05
  real(rprec) :: dt, dt_f, dt_dim, cfl_f
  
  ! time advance parameters (Adams-Bashforth, 2nd order accurate)
  real (rprec) :: tadv1, tadv2
  
  $else
  
  real (rprec), parameter :: dt = 2.0e-4_rprec              ! dt=2.e-4 usually works for 64^3
  real (rprec), parameter :: dt_dim = dt*z_i/u_star     ! dimensional time step in seconds
  
  ! time advance parameters (Adams-Bashforth, 2nd order accurate)
  real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1
  
  $endif
  
  logical, parameter :: cumulative_time = .true.        ! to use total_time.dat
  character (*), parameter :: fcumulative_time = path // 'total_time.dat'
  
  integer :: jt                 ! global time-step counter
  integer :: jt_total           ! used for cumulative time (see io module)
  real(rprec) :: total_time, total_time_dim
  
!---------------------------------------------------
! BOUNDARY/INITIAL CONDITION PARAMETERS
!---------------------------------------------------  

  ! initu = true to read from a file; false to create with random noise
  logical, parameter :: initu = .false.
  ! initlag = true to initialize cs, FLM & FMM; false to read from vel.out
  logical, parameter :: inilag = .true.

  ! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
  integer,parameter::ubc=0
  ! lbc: lower boundary condition:  'wall', 'stress free'
  character (*), parameter :: lbc_mom = 'wall'
  
  ! lower boundary condition, roughness length
  ! if use_default_patch is false, zo will be read from 'patch.dat'
  logical, parameter :: use_default_patch = .true.
  real (rprec), parameter :: zo_default = 0.0001_rprec  ! nondimensional  

  ! prescribed inflow:   
  logical,parameter::inflow=.false.
  ! if inflow is true the following should be set:
    logical, parameter :: use_fringe_forcing = .false.  
    ! position of right end of buffer region, as a fraction of L_x
    real (rprec), parameter :: buff_end = 1._rprec
    ! length of buffer region as a fraction of L_x
    real (rprec), parameter :: buff_len = 0.125_rprec  
    real (rprec), parameter :: face_avg = 1.0_rprec
    ! true to read from file; false to set as constant
    ! read from file is not working properly
    logical, parameter :: read_inflow_file = .false.
    logical, parameter :: write_inflow_file = .false.
    ! records at position jx_s
    integer, parameter :: jt_start_write = 6
    ! forcing along top and bottom bdrys
    ! if inflow is true and force_top_bot is true, then the top & bottom
    ! velocities are forced to the inflow velocity
    logical, parameter :: force_top_bot = .false.

    ! If true the inflow will be forced to the velocity at a sampled location
    ! Use instead of face_avg and read_inflow_file
    logical, parameter :: inflow_sample_velocity=.true.
    ! Sample location as a fraction of L_x
    real(rprec), parameter :: inflow_sample_location=0.5_rprec

  ! if true, imposes a pressure gradient in the x-direction to force the flow
  logical, parameter :: use_mean_p_force = .true.
  real (rprec), parameter :: mean_p_force = 1._rprec / L_z
  
!---------------------------------------------------
! DATA OUTPUT PARAMETERS
!---------------------------------------------------

  ! how often to display "jt,dt,rmsdivvel,ke,cfl" output
  integer,parameter::wbase=100
  
  ! how often to write ke to check_ke.out
  integer, parameter :: nenergy = 100

  ! how often to display Lagrangian CFL condition of 
  ! dynamic SGS models
  integer,parameter::cfl_count=1000

  ! records time-averaged data to files ./output/*_avg.dat
  logical, parameter :: tavg_calc = .true.
  integer, parameter :: tavg_nstart = 50000, tavg_nend = nsteps

  ! turns instantaneous velocity recording on or off
  logical, parameter :: point_calc = .true.
  integer, parameter :: point_nstart = 50000, point_nend = nsteps, point_nskip = 100
  integer, parameter :: point_nloc = 1
  type(point3D), dimension(point_nloc) :: point_loc = (/ &
        point3D( (/ L_x / 2, L_y / 2, L_z / 2 /) ) &
        /)

  ! domain instantaneous output
  logical, parameter :: domain_calc = .true.
  integer, parameter :: domain_nstart = 50000, domain_nend = nsteps, domain_nskip = 50000
  
  ! x-plane instantaneous output
  logical, parameter :: xplane_calc   = .true.
  integer, parameter :: xplane_nstart = 50000, xplane_nend = nsteps, xplane_nskip  = 10000
  integer, parameter :: xplane_nloc   = 2
  real(rprec), dimension(xplane_nloc) :: xplane_loc = (/ L_x/4._rprec, L_x/2._rprec /)

  ! y-plane instantaneous output
  logical, parameter :: yplane_calc   = .false.
  integer, parameter :: yplane_nstart = 50000, yplane_nend = nsteps, yplane_nskip  = 50000
  integer, parameter :: yplane_nloc   = 2
  real(rprec), dimension(yplane_nloc) :: yplane_loc = (/ L_y/4._rprec, L_y/2._rprec  /)  

  ! z-plane instantaneous output
  logical, parameter :: zplane_calc   = .true.
  integer, parameter :: zplane_nstart = 50000, zplane_nend = nsteps, zplane_nskip  = 10000
  integer, parameter :: zplane_nloc   = 3
  real(rprec), dimension(zplane_nloc) :: zplane_loc = (/ 0.1_rprec, 0.25_rprec, 0.5_rprec /)

  logical, parameter :: spectra_calc = .true.
  integer, parameter :: spectra_nstart = 50000, spectra_nend = nsteps
  integer, parameter :: spectra_nloc = 3
  real(rprec), dimension(spectra_nloc) :: spectra_loc = (/ 0.1_rprec, 0.25_rprec, 0.5_rprec /)

!--------------------------------------------------- 
! SCALAR PARAMETERS
!---------------------------------------------------
 
  ! S_FLAG=1 for Theta and q, =0 for no scalars
  ! logical,parameter::S_FLAG=.TRUE.,coupling_flag=.FALSE.,mo_flag=.TRUE.
  logical,parameter::S_FLAG=.false.
  ! integer,parameter::DYN_init=2, SCAL_init=5, no_days=1
  integer,parameter :: SCAL_init=5, no_days=1
  ! integer,parameter::DYN_init=1, SCAL_init=1, no_days=1
  integer,parameter::patch_flag=1, remote_flag=0, time_start=0
  ! initu=.TRUE. & initsc=.FALSE read velocity fields from a binary file
  ! initu=.TRUE. & initsc=.TRUE. read velocity & scalar fields from a binary file
  ! initu=.FALSE. & S_FLAG=.TRUE. initialize velocity & scalar fields using ran
  ! initu=.FALSE. & S_FLAG=.FALSE. initialize velocity fields using ran
  logical,parameter::initsc=.false.
  ! lbc=0: prescribed surface temperature, lbc=1 prescribed surface flux
  ! (wT=0.06 Km/s)
  integer,parameter :: lbc=0
  ! Added a new parameter - sflux_flag for passive scalars with bldngs
  logical,parameter :: sflux_flag=.false.
  logical,parameter :: wt_evolution_flag=.FALSE.
  logical,parameter :: test_phase=.FALSE., vec_map=.FALSE., smag_sc=.FALSE.
  logical,parameter :: check_dt=.TRUE.
  integer,parameter :: stencil_pts=4
  logical,parameter :: coarse_grain_flag=.FALSE.
  ! inversion strength (K/m)
  real(kind=rprec),parameter::g=9.81_rprec, inv_strength=0._rprec
  real(kind=rprec),parameter::theta_top=300._rprec,T_scale=300._rprec&
       ,wt_s=20._rprec, T_init=300._rprec
  real(kind=rprec),parameter::cap_thick=80._rprec, z_decay=1._rprec
  integer,parameter::c_count=10000,p_count=10000  

end module param
