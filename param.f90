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


logical, parameter :: tloop = .false.
logical, parameter :: VERBOSE = .false.  !--prints small stuff to screen
                      !--use DEBUG to write lots of data/files

integer, parameter :: nz_tot = (nz - 1) * nproc + 1
integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
integer,parameter:: lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
real (rprec), parameter :: BOGUS = -1234567890._rprec

real(rprec),parameter::pi=3.1415926535897932384626433_rprec
!real(rprec),parameter::z_i=1._rprec, L_z=(1._rprec * z_i)/nproc

                            !--L_z is not nondimensionalized by z_i yet
! set the aspect ratio of the box, already nondimensional
real(rprec),parameter::dz=L_z/z_i/(nz-1)
real(rprec),parameter::dx=L_x/nx,dy=L_y/ny

integer, parameter :: nsteps = 2000
!  Commented out for now; see below u_star declaration
!  for details
!real (rprec), parameter :: dt = 2.e-6_rprec / 1._rprec
                           !--dt=2.e-4 usually works for 64^3

! u_star=0.45 if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
real(rprec),parameter::u_star=0.45_rprec,Pr=.4_rprec

!  dt_dim is not present in this version but is present int
!  MARCELO_CODE; adding it here to see what we get
real(rprec),parameter::dt_dim=0.1_rprec !dimensional time step in seconds
real(rprec),parameter::dt=dt_dim*u_star/z_i

!--Coriolis stuff
! coriol=non-dim coriolis parameter,
! ug=horiz geostrophic vel, vg=transverse geostrophic vel
! u_star=0.45 if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.                         
real(rprec),parameter::coriol=9.125E-05*z_i/u_star,      &
                      ug=u_star/u_star,vg=0._rprec/u_star

real(rprec),parameter::vonk=.4_rprec

integer,parameter::c_count=200,p_count=200
!integer, parameter :: cs_count = 1  !--tsteps between dynamic Cs updates
integer, parameter :: cs_count = 5  !--tsteps between dynamic Cs updates
logical,parameter::output=.true.
logical, parameter :: use_avgslice = .true.

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
integer,parameter::model=5,models=1,nnn=2
real(kind=rprec),parameter::Co=0.16_rprec
!  This was not originally here
real(kind=rprec),parameter::cs=0.2_rprec

!Test filter type: 1->cut off 2->Gaussian 3->Top-hat
integer,parameter::ifilter=2

! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
integer,parameter::ubc=1

character (*), parameter :: lbc_mom = 'wall'
                            !--'wall', 'stress free'

!--prescribed inflow: constant or read from file
!  read from file is not working properly
logical,parameter::inflow=.true.
logical, parameter :: use_fringe_forcing = .false.

real (rprec), parameter :: buff_end = 1._rprec
                           !--position of right end of buffer region,
                           !  as a fraction of L_x
real (rprec), parameter :: buff_len = 0.25_rprec
                           !--length of buffer region as a fraction of L_x
!real (rprec), parameter :: face_avg = 0.0_rprec
real (rprec), parameter :: face_avg = 1.0_rprec

! forcing along top and bottom bdrys
! if inflow is true and force_top_bot is true, then the top & bottom
! velocities are forced to the inflow velocity
logical, parameter :: force_top_bot = .false.

logical, parameter :: use_mean_p_force = .false.
real (rprec), parameter :: mean_p_force = 1._rprec * z_i/L_z/nproc
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
integer,parameter::DYN_init=50, SCAL_init=5, no_days=1
!integer,parameter::DYN_init=1, SCAL_init=1, no_days=1
integer,parameter::patch_flag=1, remote_flag=0, time_start=0
! initu=.TRUE. & initsc=.FALSE read velocity fields from a binary file
! initu=.TRUE. & initsc=.TRUE. read velocity & scalar fields from a binary file
! initu=.FALSE. & S_FLAG=.TRUE. initialize velocity & scalar fields using ran
! initu=.FALSE. & S_FLAG=.FALSE. initialize velocity fields using ran
logical,parameter::initsc=.false.
! lbc=0: prescribed surface temperature, lbc=1 prescribed surface flux
! (wT=0.06 Km/s)
integer,parameter :: lbc=1
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
