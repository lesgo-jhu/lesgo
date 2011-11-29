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

  $if ($MPI)
  integer, parameter :: lbz = 0  ! overlap level for MPI transfer
  $else
  integer, parameter :: lbz = 1  ! no overlap level necessary
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
  integer :: jzmin, jzmax  ! levels that "belong" to this processor, set w/ grid
  !--end mpi stuff
  
!---------------------------------------------------
! COMPUTATIONAL DOMAIN PARAMETERS
!---------------------------------------------------  
! characteristic length is H=z_i and characteristic velocity is u_star  
!   L_x, L_y, L_z, dx, dy, dz are non-dim. using H

  integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
  real (rprec), parameter :: BOGUS = -1234567890._rprec
  real(rprec),parameter::pi=3.1415926535897932384626433_rprec

  !integer,parameter:: nx=64,ny=64,nz=(64)/nproc + 1 
  integer :: Nx, Ny, Nz
  
  !integer, parameter :: nz_tot = (nz - 1) * nproc + 1
  !integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
  !integer,parameter:: lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big
  integer :: nz_tot
  integer :: nx2, ny2
  integer :: lh, ld, lh_big, ld_big

  ! this value is dimensional [m]:
  !real(rprec),parameter::z_i=1000._rprec   !dimensions in meters, height of BL
  real(rprec) :: z_i
    
  ! these values should be non-dimensionalized by z_i: 
  ! set as multiple of BL height (z_i) then non-dimensionalized by z_i
  !real(rprec),parameter::L_x= 2*pi
  !real(rprec),parameter::L_y= L_x
  !real(rprec),parameter::L_z= 1.0_rprec
  !real(rprec),parameter::L_y=ny*L_x/nx               ! ensure dy=dx
  !real(rprec),parameter::L_z=(nz_tot - 1)*L_x/nx  ! ensure dz = dx
  real(rprec) :: L_x, L_y, L_z

  ! these values are also non-dimensionalized by z_i:
  !real(rprec),parameter::dz=L_z/(nz_tot-1.) ! or (L_z/nproc)/(nz - 1)
  !real(rprec),parameter::dx=L_x/nx,dy=L_y/ny
  real(rprec) :: dx, dy, dz
  
!---------------------------------------------------
! MODEL PARAMETERS
!---------------------------------------------------   

  ! Model type: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent
  !             4->Lagrangian scale-sim   5-> Lagragian scale-dep
  !integer,parameter::model=5, nnn=2
  integer :: sgs_model, wall_damp_exp

  ! timesteps between dynamic Cs updates           
  !integer, parameter :: cs_count = 5
  integer :: cs_count

  ! When to start dynamic Cs calculations
  !integer,parameter::DYN_init=100
  integer :: DYN_init
  
  ! Cs is the Smagorinsky Constant
  ! Co and nnn are used in the mason model for smagorisky coeff
  !real(kind=rprec),parameter::Co=0.16_rprec
  real(rprec) :: Co
  
  ! test filter type: 1->cut off 2->Gaussian 3->Top-hat
  !integer,parameter::ifilter=1
  integer :: ifilter

  ! u_star=0.45 m/s if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
  real(rprec),parameter::u_star=0.45_rprec

  ! von Karman constant     
  real(rprec),parameter::vonk=0.4_rprec   
  
  ! Coriolis stuff
  ! coriol=non-dim coriolis parameter,
  ! ug=horiz geostrophic vel, vg=transverse geostrophic vel
  !logical,parameter::coriolis_forcing=.false.
  !real(rprec),parameter::coriol=9.125E-05*z_i/u_star,      &
  !     ug=u_star/u_star,vg=0._rprec/u_star
  logical :: coriolis_forcing
  real(rprec) :: coriol, ug, vg

  ! nu_molec is dimensional m^2/s
  !real(rprec),parameter::nu_molec=1.14e-5_rprec   
  real(rprec) :: nu_molec
    
  !logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.  
  logical :: molec, sgs, dns_bc
  
!---------------------------------------------------
! TIMESTEP PARAMETERS
!---------------------------------------------------   

  !integer, parameter :: nsteps = 1500
  integer :: nsteps
 
  $if($CFL_DT)
  
  !real(rprec), parameter :: cfl = 0.05
  real(rprec) :: cfl
  real(rprec) :: dt, dt_f, dt_dim, cfl_f
  
  ! time advance parameters (Adams-Bashforth, 2nd order accurate)
  real (rprec) :: tadv1, tadv2
  
  $else
  
  !real (rprec), parameter :: dt = 2.0e-4_rprec              ! dt=2.e-4 usually works for 64^3
  !real (rprec), parameter :: dt_dim = dt*z_i/u_star     ! dimensional time step in seconds
  real(rprec) :: dt
  real(rprec) :: dt_dim
  
  ! time advance parameters (Adams-Bashforth, 2nd order accurate)
  real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1
  
  $endif
  
  logical, parameter :: cumulative_time = .false.        ! to use total_time.dat
  character (*), parameter :: fcumulative_time = path // 'total_time.dat'
  
  integer :: jt                 ! global time-step counter
  integer :: jt_total           ! used for cumulative time (see io module)
  real(rprec) :: total_time, total_time_dim
  
!---------------------------------------------------
! BOUNDARY/INITIAL CONDITION PARAMETERS
!---------------------------------------------------  

  ! initu = true to read from a file; false to create with random noise
  !logical, parameter :: initu = .false.
  logical :: initu
  ! initlag = true to initialize cs, FLM & FMM; false to read from vel.out
  !logical, parameter :: inilag = .true.
  logical :: inilag

  ! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
  !integer,parameter::ubc=0
  integer :: ubc
  ! lbc: lower boundary condition:  'wall', 'stress free'
  !character (*), parameter :: lbc_mom = 'wall'
  character(50) :: lbc_mom
  
  ! lower boundary condition, roughness length
  !real (rprec), parameter :: zo = 0.0001_rprec  ! nondimensional  
  real(rprec) :: zo

  ! prescribed inflow:   
  !logical,parameter::inflow=.false.
  logical :: inflow
  ! if inflow is true the following should be set:
    ! position of right end of fringe region, as a fraction of L_x
    !real (rprec), parameter :: fringe_region_end = 1._rprec
    real(rprec) :: fringe_region_end
    ! length of fringe region as a fraction of L_x
    !real (rprec), parameter :: fringe_region_len = 0.125_rprec  
    real(rprec) :: fringe_region_len

    ! Use uniform inflow instead of concurrent precursor inflow
    !logical, parameter :: uniform_inflow = .false.
    logical :: uniform_inflow
      !real (rprec), parameter :: inflow_velocity = 1.0_rprec
      real(rprec) :: inflow_velocity
      ! velocities are forced to the inflow velocity
      !logical, parameter :: force_top_bot = .false.
      logical :: force_top_bot

    ! Use concurrent precursor setup instead of uniform inflow
    ! Make sure USE_CPS=true in Makefile.in
    logical, parameter :: concurrent_precursor_inflow = .true. 

  ! if true, imposes a pressure gradient in the x-direction to force the flow
  !logical, parameter :: use_mean_p_force = .true.
  logical :: use_mean_p_force
  !real (rprec), parameter :: mean_p_force = 1._rprec / L_z
  real(rprec) :: mean_p_force
  
!---------------------------------------------------
! DATA OUTPUT PARAMETERS
!---------------------------------------------------

  ! how often to display "jt,dt,rmsdivvel,ke,cfl" output
  !integer,parameter::wbase=100
  integer :: wbase
  
  ! how often to write ke to check_ke.out
  !integer, parameter :: nenergy = 100
  integer :: nenergy

  ! how often to display Lagrangian CFL condition of 
  ! dynamic SGS models
  !integer,parameter::cfl_count=1000
  integer :: cfl_count

  ! records time-averaged data to files ./output/*_avg.dat
  !logical, parameter :: tavg_calc = .true.
  !integer, parameter :: tavg_nstart = 50000, tavg_nend = nsteps
  logical :: tavg_calc
  integer :: tavg_nstart, tavg_nend

  ! turns instantaneous velocity recording on or off
  !logical, parameter :: point_calc = .true.
  !integer, parameter :: point_nstart = 50000, point_nend = nsteps, point_nskip = 100
  !integer, parameter :: point_nloc = 1
  !type(point3D), parameter, dimension(point_nloc) :: point_loc = (/ &
  !      point3D( (/ L_x / 2, L_y / 2, L_z / 2 /) ) &
  !      /)
  logical :: point_calc
  integer :: point_nstart, point_nend, point_nskip
  integer :: point_nloc
  type(point3D), allocatable, dimension(:) :: point_loc

  ! domain instantaneous output
  !logical, parameter :: domain_calc = .true.
  !integer, parameter :: domain_nstart = 50000, domain_nend = nsteps, domain_nskip = 50000
  logical :: domain_calc
  integer :: domain_nstart, domain_nend, domain_nskip
  
  ! x-plane instantaneous output
  !logical, parameter :: xplane_calc   = .true.
  !integer, parameter :: xplane_nstart = 50000, xplane_nend = nsteps, xplane_nskip  = 10000
  !integer, parameter :: xplane_nloc   = 2
  !real(rprec), parameter, dimension(xplane_nloc) :: xplane_loc = (/ L_x/4._rprec, L_x/2._rprec /)
  logical :: xplane_calc
  integer :: xplane_nstart, xplane_nend, xplane_nskip
  integer :: xplane_nloc
  real(rprec), allocatable, dimension(:) :: xplane_loc

  ! y-plane instantaneous output
  ! logical, parameter :: yplane_calc   = .false.
  ! integer, parameter :: yplane_nstart = 50000, yplane_nend = nsteps, yplane_nskip  = 50000
  ! integer, parameter :: yplane_nloc   = 2
  ! real(rprec), parameter, dimension(yplane_nloc) :: yplane_loc = (/ L_y/4._rprec, L_y/2._rprec  /)  
  logical :: yplane_calc
  integer :: yplane_nstart, yplane_nend, yplane_nskip
  integer :: yplane_nloc
  real(rprec), allocatable, dimension(:) :: yplane_loc

  ! z-plane instantaneous output
  ! logical, parameter :: zplane_calc   = .true.
  ! integer, parameter :: zplane_nstart = 50000, zplane_nend = nsteps, zplane_nskip  = 10000
  ! integer, parameter :: zplane_nloc   = 3
  ! real(rprec), parameter, dimension(zplane_nloc) :: zplane_loc = (/ 0.1_rprec, 0.25_rprec, 0.5_rprec /)
  logical :: zplane_calc
  integer :: zplane_nstart, zplane_nend, zplane_nskip
  integer :: zplane_nloc
  real(rprec), allocatable, dimension(:) :: zplane_loc

  ! logical, parameter :: spectra_calc = .true.
  ! integer, parameter :: spectra_nstart = 50000, spectra_nend = nsteps
  ! integer, parameter :: spectra_nloc = 3
  ! real(rprec), parameter, dimension(spectra_nloc) :: spectra_loc = (/ 0.1_rprec, 0.25_rprec, 0.5_rprec /) 
  logical :: spectra_calc
  integer :: spectra_nstart, spectra_nend
  integer :: spectra_nloc
  real(rprec), allocatable, dimension(:) :: spectra_loc


end module param
