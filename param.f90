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

  !----- Variable Declaration for Input File -------

  !----- Variable Declaration for Input File -------


logical, parameter :: VERBOSE = .false.  !--prints small stuff to screen
                      !--use DEBUG to write lots of data/files

integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
real (rprec), parameter :: BOGUS = -1234567890._rprec

real(rprec),parameter::pi=3.1415926535897932384626433_rprec
!real(rprec),parameter::z_i=1._rprec, L_z=(1._rprec * z_i)/nproc

!  U intialization for non-log profile IC
logical, parameter :: ic_const = .false.
real (rprec), parameter :: u_ic = 20.0/u_star, v_ic=0., w_ic=0.
                  


real(rprec),parameter::vonk=.4_rprec

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
