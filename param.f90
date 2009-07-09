module param
  use types,only:rprec
  $if ($MPI)
  use mpi
  $endif
  implicit none

  save

  private rprec  !--this is dumb.
  public

  !--mpi stuff
  $if ($MPI)
  $define $MPI_LOGICAL .true.
  $define $NPROC 2
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

  logical, parameter :: VERBOSE = .false.  !--prints small stuff to screen
  !--use DEBUG to write lots of data/files

  integer,parameter:: nx=64,ny=64,nz=(64+(nproc-1)-1)/nproc + 1
  integer, parameter :: nz_tot = (nz - 1) * nproc + 1
  integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
  integer,parameter:: lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

  integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
  real (rprec), parameter :: BOGUS = -1234567890._rprec
  
  integer, parameter :: nsteps = 100

  real(rprec),parameter::pi=3.1415926535897932384626433_rprec
    !real(rprec),parameter::z_i=1._rprec, L_z=(1._rprec * z_i)/nproc
  real(rprec),parameter::z_i=1._rprec
  real(rprec),parameter::L_x=4.*z_i, L_y=4.*z_i, L_z=4.*z_i/nproc
  !--L_z is not nondimensionalized by z_i yet
  ! set the aspect ratio of the box, already nondimensional
  real(rprec),parameter::dz=nproc*L_z/z_i/(nz_tot-1./2.)
  real(rprec),parameter::dx=L_x/(nx-1),dy=L_y/(ny-1)

  ! u_star=0.45 if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
  real(rprec),parameter::u_star=0.45_rprec,Pr=.4_rprec
  
  !  U intialization for non-log profile IC
  logical, parameter :: ic_const = .false.
  real (rprec), parameter :: u_ic = 10.0/u_star, v_ic=0., w_ic=0.

  real (rprec), parameter :: dt = 2.e-4
  real (rprec), parameter :: dt_dim = dt*z_i/u_star
  
!  real(rprec),parameter::dt_dim=0.1 !dimensional time step in seconds
!  real(rprec),parameter::dt=dt_dim*u_star/z_i
                                  !--dt=2.e-4 usually works for 64^3
  
  !--Coriolis stuff
  ! coriol=non-dim coriolis parameter,
  ! ug=horiz geostrophic vel, vg=transverse geostrophic vel
  logical,parameter::coriolis_forcing=.false.
  real(rprec),parameter::coriol=9.125E-05*z_i/u_star,      &
       ug=u_star/u_star,vg=0._rprec/u_star

  real(rprec),parameter::vonk=0.4_rprec 
  integer,parameter::c_count=10000,p_count=10000
 !integer, parameter :: cs_count = 1  !--tsteps between dynamic Cs updates
  integer, parameter :: cs_count = 5  !--tsteps between dynamic Cs updates
  logical,parameter:: output=.true.
  logical, parameter :: use_avgslice = .false.
  !  Set minimum time step to write averaged slices
  integer, parameter :: avgslice_start = 0
  
  
  !--initu = true to read from a file; false to create with random noise
  logical, parameter :: initu = .false.
  !--initlag = true to initialize cs, FLM & FMM; false to read from vel.out
  logical, parameter :: inilag = .true.

  ! nu_molec is dimensional m^2/s
  real(rprec),parameter::nu_molec=1.14e-5_rprec

  logical,parameter::use_bldg=.false.
  logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.

  !Model type: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent
  !            4->Lagrangian scale-sim   5-> Lagragian scale-dep
  !Models type: 1->static prandtl, 2->Dynamic
  !Cs is the Smagorinsky Constant
  !Co and nnn are used in the mason model for smagorisky coeff
  integer,parameter::model=1,models=1,nnn=2
  real(kind=rprec),parameter::Co=0.2_rprec
  !  This was not originally here
  real(kind=rprec),parameter::cs=0.16_rprec

  !Test filter type: 1->cut off 2->Gaussian 3->Top-hat
  integer,parameter::ifilter=2

  ! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
  integer,parameter::ubc=0

  character (*), parameter :: lbc_mom = 'wall'
  !--'wall', 'stress free'

  !--prescribed inflow: constant or read from file
  !  read from file is not working properly
  logical,parameter::inflow=.false.
  logical, parameter :: use_fringe_forcing = .false.

  real (rprec), parameter :: buff_end = 1._rprec
  !--position of right end of buffer region,
  !  as a fraction of L_x
  real (rprec), parameter :: buff_len = 0.25_rprec
  !--length of buffer region as a fraction of L_x
  !real (rprec), parameter :: face_avg = 0.0_rprec
  real (rprec), parameter :: face_avg = 1.0_rprec

  logical, parameter :: read_inflow_file = .false.
  logical, parameter :: write_inflow_file = .false.

  !--records at position jx_s
  integer, parameter :: jt_start_write = 6

  ! forcing along top and bottom bdrys
  ! if inflow is true and force_top_bot is true, then the top & bottom
  ! velocities are forced to the inflow velocity
  logical, parameter :: force_top_bot = .false.

  logical, parameter :: use_mean_p_force = .true.
  real (rprec), parameter :: mean_p_force = 1._rprec * z_i/(nproc*L_z)
  !--usually just z_i/L_z

  integer :: jt  ! global time-step counter
  integer :: jt_total  !--used for cumulative time (see io module)

  ! time advance parameters (AB2)
  real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1

  !------xxxxxxxxx--SCALARS_PARAMETERS--xxxxxxxxx---------------
  ! S_FLAG=1 for Theta and q, =0 for no scalars
  !logical,parameter::S_FLAG=.TRUE.,coupling_flag=.FALSE.,mo_flag=.TRUE.
  logical,parameter::S_FLAG=.false.
  !integer,parameter::DYN_init=2, SCAL_init=5, no_days=1
  integer,parameter::DYN_init=100, SCAL_init=5, no_days=1
  !integer,parameter::DYN_init=1, SCAL_init=1, no_days=1
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
  !inversion strength (K/m)
  real(kind=rprec),parameter::g=9.81_rprec, inv_strength=0._rprec
  real(kind=rprec),parameter::theta_top=300._rprec,T_scale=300._rprec&
       ,wt_s=20._rprec,T_init=300._rprec
  real(kind=rprec),parameter::cap_thick=80._rprec, z_decay=1._rprec


end module param
