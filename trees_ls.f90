!--version of tree module to be used with level-set
!--supposed to be a lightweight version of tree module
!--this version has full MPI-support
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module trees_ls
use trees_base_ls
use trees_setup_ls, only : fill_tree_array
use trees_io_ls, only : write_ta_data
use trees_fmodel_ls, only : clip_fcoeff, def_d_dir, def_dls_dir, def_nba_dir,  &
                            do_tavg_fcoeff, do_tavg_fdist, fcoeff_ta,          &
                            fcoeff_d_germano, fcoeff, fcoeff_clipped,          &
                            fcoeff_divide, fill_unresf_br,                     &
                            vel_d, vel_dls, vel_nba
                             
use trees_global_fmask_ls, only : read_global_fmask, global_fmask
use param, only : ld, nx, ny, nz, dx, dy, dz, USE_MPI, coord, nproc,   &
                  ierr, MPI_RPREC, comm, rank_of_coord, jt, jt_total,  &
                  nsteps, BOGUS, chcoord
use level_set_base, only : phi
use messages

$if ($DEBUG)
use debug_mod
$endif

$if ($MPI)
  use mpi
$endif
implicit none

save
$if ($TEST)
  public
$else
  private
$endif

public :: trees_ls_calc, trees_ls_finalize, trees_ls_init
public :: trees_ls_calc_init, apriori
public :: brindex, brindex_initialized
public :: phi

character (*), parameter :: mod_name = 'trees_ls'

!--select what kind of test to perform: a-priori or a-posteriori
!character (*), parameter :: test = 'a-priori'
character (*), parameter :: test = 'a-posteriori'
!--force models availability is controlled by fmodel module
character (*), parameter :: file_fcoeff_last = 'fcoeff.last.out'
character (*), parameter :: file_ftot_last = 'ftot.last.out'
character (*), parameter :: file_fcoeff_dat = 'fcoeff.dat'
character (*), parameter :: file_fdist_last = 'fdist.last.out'
character (*), parameter :: file_num_den_fcoeff_last = 'num_den_fcoeff.last.out'
character (*), parameter :: file_num_den_fdist_last = 'num_den_fdist.last.out'

!character (*), parameter :: flocal = 'unresf_dir'
character (*), parameter :: flocal = 'unresf_dir'
    !--'unresf_dir', 'fmodel'
    !--only active when use_local_vel is true
    !--'unresf_dir': rule for calculating local force distribution has same
    !  form as that for the "bulk" force model that is used to calculate
    !  force coefficients, e.g., drag coefficients

!--counters for a-posteriori test
integer, parameter :: nCDupdate = 5
integer, parameter :: nfdistupdate = 5
integer, parameter :: noutput = 500
integer, parameter :: nglobalCD = 50

!--counters for a-priori test
integer, parameter :: napriori = 5

integer, parameter :: nfxyz_write = 2000

logical, parameter :: init_from_file = .false.
                      !--this is for both fcoeff and fdist

logical, parameter :: unresf_clip = .true. 
logical, parameter :: unresf_dir_proj = .false.
    !--whether to project local velocity when using unresf_dir

logical, parameter :: time_avg = .false.  !--time avg fcoeff
logical, parameter :: use_unif_fdist_br = .false.
                      !--measures fdist in RESOLVED br's as constant avg. val.
logical, parameter :: use_local_vel = .true.
logical, parameter :: use_global_fmask = .false.

logical :: brindex_initialized = .false.
integer :: brindex(ld, ny, nz)
integer :: wksp_size(nd)  !--only the first 3 dimensions of wksp
integer :: max_res_gen

real (rp), allocatable :: fdist_wksp(:, :, :, :)
real (rp), allocatable :: fdist_mask(:, :, :)

!--used to restart the time averaging when separate num, den are used
!--these are module-variables to allow use with trees_finalize_ls
real (rp) :: num_fcoeff(nfcoeff, nzone) = BOGUS
real (rp) :: den_fcoeff(nfcoeff, nzone) = BOGUS
real (rp), allocatable :: num_fdist(:, :, :, :)  !--wksp_size(:) X nzone
real (rp) :: den_fdist(nzone) = BOGUS

$if ($MPI)
  integer :: nbr  !--number of branches

  !--refer to these as "reduce arrays": store some branch info in 
  !  one place, only for arrays that are to be "reduced" with MPI
  !--these arrays should make MPI communication more efficient
  !  may want to extend this idea to non-MPI version, that way
  !  no packing/unpacking is needed for MPI version  (even better)
  integer, allocatable :: nvelscale(:), nvelscale_sum(:)  
  integer, allocatable :: nrespt(:), nrespt_sum(:)

  real (rp), allocatable :: fnorm(:, :, :), fnorm_sum(:, :, :)
  real (rp), allocatable :: resf(:, :), resf_sum(:, :)
  real (rp), allocatable :: velscale(:, :), velscale_sum(:, :)
$endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trees_ls_calc_init ()
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_ls_calc_init'

character (64) :: sn

logical, save :: first_call = .true.

!---------------------------------------------------------------------

sn = trim ( sub_name ) // trim ( chcoord )

$if ($VERBOSE)
call enter_sub ( sn )
$endif

if (.not. tree_array_initialized) then

  !--initialize tree array, etc.
  call grid_initialize ()
  call fill_tree_array ()
  if ( .not. brindex_initialized ) call brindex_init ()
      !--this allows for external initialization
      !  e.g. for post-processing, when need to merge MPI brindex files

  $if ($MPI)
      call reduce_array_init ()
  $endif
 
  if ( use_unif_fdist_br ) call fill_nrespt_ta ()
      !--MUST be after reduce_array_init, since array
      !  nrespt, nrespt_sum needs to be allocated

  tree_array_initialized = .true.
  
end if

!--other initialization not related to tree_array init
if (first_call) then
  call set_max_res_gen ()
  call mesg (sub_name, 'max_res_gen =', max_res_gen)
  call mesg (sub_name, 'fmodel =' // fmodel)
  call mesg (sub_name, 'init_from_file =', init_from_file)
  first_call = .false.
end if

$if ($VERBOSE)
call exit_sub ( sn )
$endif

end subroutine trees_ls_calc_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this will be called from main code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trees_ls_calc ()
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_ls_calc'

logical, save :: first_call = .true.

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

if ( first_call ) then
    call trees_ls_calc_init ()
    first_call = .false.
end if

select case ( test )
case ( 'a-posteriori' )
    call aposteriori ()
case ( 'a-priori' )
    if ( modulo (jt, napriori) == 0 ) call apriori ()
case default
    call error ( sub_name, 'invalid test =' // test )
end select

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine trees_ls_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine aposteriori ()
implicit none

character (*), parameter :: sub_name = mod_name // '.aposteriori'

character (128) :: foutput

!---------------------------------------------------------------------

!--this routine assumes there are SGS branches

!--update velscale at every timestep
!--fill velscale: currently fills all generations, but only
!  generations >= max_res_gen are used for dynamic CD
!  a priori test will require all scales
!--at every timestep we need velscale @ max_res_gen + 1 = n_gen
call fill_velscale_ta ()

if ((modulo (jt, nCDupdate) == 0) .or. (jt == 1)) then

  !--fill resf (all branch generations)
  !--may only want to calculate at max_res_gen @ nCDupdate,
  !  and calculate for all generations less frequently
  !  (could do this with optional argument)
  call fill_resf_ta ()

  !--calculate new dynamic CD
  !--the old unresf (using old CD and old vel)  will be used in calculation
  !  of non-germano CD unless fill_unresf_ta and fill_ftot_ta are called here
  !--on first time step there is no old CD to use, so these calls
  !  actually need to happen inside update_fcoeff, after the old CD has
  !  been read from a file
  !--WARNING: is nCDupdate /= nfdistupdate, make sure that fill_ftot_ta
  !  is called when using non-Germano
  call update_fcoeff ()  !--this includes special code to handle jt = 1

end if

!--calculate unresolved branch force to apply to RHS of momentum
!--total unresolved branch force, not the spatial distribution
!--if choose not to update velscale every time step, then this
!  should only be updated when fcoeff or velscale changes
call fill_unresf_ta ()

if ( (modulo ( jt, nfdistupdate ) == 0 ) .or. ( jt == 1 ) ) then

  !--update_fdist expects that ftot is updated before its called
  call fill_ftot_ta ()  !--right now, this is calculated for all gens

  if (use_unif_fdist_br) then  !--experimental
    !--some routines within update_fdist assume resf has been
    !  updated, so we need special call when nfdistupdate /= nCDupdate
    if (nfdistupdate /= nCDupdate) call fill_resf_ta ()
  end if
  
  call update_fdist ()  !--includes special code for jt = 1
  
end if

!--apply unresf using fdist
!--this must be done every time step
call apply_unresf_ta ()
!call mesg (sub_name, 'apply_unresf_ta has been disabled')  !--for debug

!--all procs must call this routine
!--output global CD, and split into resolved, unresolved parts
!--uses the info from the latest call to fill_ftot_ta above
if (modulo (jt, nglobalCD) == 0) call calc_globalCD ()

!--general output: always called on last time step
if ((modulo (jt, noutput) == 0) .or. (jt == nsteps)) then

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    
    if (.not. USE_MPI) then
      write (foutput, '(a,i6.6,a)') 'tree_data.', jt_total, '.dat' 
    else
      write (foutput, '(a,i6.6,a,i0)') 'tree_data.', jt_total, '.dat.c', coord 
    end if

    call write_ta_data (foutput)

  end if

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

    write (foutput, '(a,i6.6,a)') 'output/fcoeff.', jt_total, '.dat'
    call write_fcoeff (foutput)

    write (foutput, '(a,i6.6,a)') 'output/fdist.', jt_total, '.dat'
    call output_fdist_wksp (foutput)  !--this uses formatted i/o

    !--all procs need to do this
    !write (foutput, '(a,i6.6,a)') 'output/fxyz.', jt_total, '.dat'
    !call write_fxyz (foutput)  !--this uses unformatted i/o

  end if

  if ( modulo (jt, nfxyz_write) == 0 ) then
    !--all procs need to do this
    write (foutput, '(a,i6.6,a)') 'output/fxyz.', jt_total, '.dat'
    call write_fxyz (foutput)  !--this uses unformatted i/o
  end if

  !if (DEBUG .and. USE_MPI) then
  !  !--write fdist for each process
  !  write (foutput, '(a,i6.6,a,i0)') 'output/fdist.', jt_total,  &
  !                                   '.dat.c', coord
  !  call output_fdist_wksp (foutput)
  !end if

end if

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine aposteriori

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--tree_data is written out so that simulation results may be
!  post-processed with another zone setup without re-running simulation
!--this routine does not need to be called at every time step
!--the output section writes out whenever this routine is called, it is
!  NOT controlled by any other counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apriori ()
use param, only : dt  !--plus module-wide
implicit none

character (*), parameter :: sub_name = mod_name // '.apriori'
character (*), parameter :: fprefix = 'apriori-fcoeff-gen'
character (*), parameter :: fsuffix = '.dat'

integer, parameter :: lun = 1

character (128) :: foutput
character (64) :: fmt

integer :: i
integer :: gen
integer :: ngen

logical, save :: apri_init = .false.
logical, save :: file_init = .false.

real (rp), allocatable, save :: num(:, :, :), den(:, :, :),  &
                                fcoeff_gen(:, :, :)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name // trim (chcoord))
$endif

ngen = tree_array(n_tree) % n_gen  !--assumes all trees are same here

if (.not. apri_init) then

  !--perform sanity check here: make sure all branches are resolved
  do i = 1, n_tree
    if (tree_array(i) % n_gen /= tree_array(i) % max_res_gen) then
      call error (sub_name, 'must have n_gen = max_res_gen')
    end if
  end do

  allocate (num(nfcoeff, nzone, 0:ngen))  !--0 is the trunk
  allocate (den(nfcoeff, nzone, 0:ngen))  !--0 is the trunk
  allocate (fcoeff_gen(nfcoeff, nzone, 0:ngen))  !--0 is the trunk

  num = 0._rp
  den = 0._rp
  fcoeff_gen = 0._rp
  
  apri_init = .true.

end if

!--these routines must calculate for all branch gen
call fill_velscale_ta ()

call fill_resf_ta ()

call fill_ftot_ta ()

!--output global CD, and split into resolved, unresolved parts
!--uses the info from the latest call to fill_ftot_ta above
!--ALL PROCS MUST CALL THIS (calculation of velocity scale)
call calc_globalCD ()

!--MPI: only one proc needs to calculate fcoeff
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  !--may want to rewrite to do all generations at once
  !--careful, some models may fail with an error at certain generations
  !  e.g. if there are no branches in a zone at a generation
  !  probably should modify this behavior--since it can cause a BRS
  !  to fail simply because the a zone is not populated
  do gen = 0, ngen  !--0 is the trunk
 
    !--need to offer model choice here
    !--num, den not used here (for now), so no time averaging
    select case (fmodel)
    case ('d', 'dls', 'nba')
      call fcoeff_ta (gen, num(:, :, gen), den(:, :, gen))
    case ('d_germano')
      call error (sub_name, 'no a-priori available for germano')
    case default
      call error (sub_name, 'invalid fmodel = ' // fmodel)
    end select
    
    !--after fcoeff_d, fcoeff(:, :) contains force coefficient for
    !  generation gen, so we need to copy these before overwritten
    fcoeff_gen(:, :, gen) = fcoeff(:, :)
    
  end do

  !--output fcoeff_gen: separate file for each generation
  !--due to the number of zones and generations, its probably best
  !  to postprocess this output into separate files for plotting
  do gen = 0, ngen

    write (foutput, '(a,i0,a)') fprefix, gen, fsuffix

    if (file_init) then
    
      open (lun, file=foutput, action='write', position='append')

    else

      open (lun, file=foutput, action='write', position='rewind')

      !--write header
      write (lun, '(a,i0,a,i0)') '# nzone =', nzone, ', ngen =', ngen

      if (gen == ngen) file_init = .true.
      
    end if
    
    write (fmt, '(a,i0,a)') '(', nfcoeff * nzone + 2, '(es13.6,1x))'
    write (lun, fmt) (jt_total - 1) * dt, fcoeff_gen(:, :, gen)
  
    close (lun)

  end do

  if (.not. USE_MPI) then
    write (foutput, '(a,i6.6,a)') 'tree_data.', jt_total, '.dat' 
  else
    write (foutput, '(a,i6.6,a,i0)') 'tree_data.', jt_total, '.dat.c', coord 
  end if

  call write_ta_data (foutput)

end if

!--all procs need to do this
!--experimental: when using in non-simulation a-priori test,
!  this will overwrite the input to the a-priori test, which
!  is at best inefficient.  However, this should not introduce
!  any errors (check to make sure that fxyz is unchanged).
!--add optional args, etc later to fix this
if ( modulo (jt, nfxyz_write) == 0 ) then
    write (foutput, '(a,i6.6,a)') 'output/fxyz.', jt_total, '.dat'
    call write_fxyz (foutput)  !--this uses unformatted i/o
end if

$if ($VERBOSE)
 call exit_sub (sub_name // trim (chcoord))
$endif

end subroutine apriori

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_nrespt_ta ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_nrespt_ta'

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

do i = 1, n_tree
  call fill_nrespt_br (tree_array(i) % trunk)
end do

$if ($MPI)

  if (.not. allocated (nrespt)) call error (sub_name, 'nrespt not allocated')
  if (.not. allocated (nrespt_sum)) call error (sub_name,                 &
                                                'nrespt_sum not allocated')

  call mpi_allreduce (nrespt, nrespt_sum, nbr, MPI_INTEGER, MPI_SUM,  &
                      comm, ierr)
  
  do i = 1, n_tree
    call irepack_br (tree_array(i) % trunk, nrespt_sum, 'nrespt')
  end do
  
$endif

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine fill_nrespt_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_nrespt_br (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_nrespt_br'

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

if (br % resolved) then

  !--calculate nrespt for this branch
  call fill_nrespt (br)
  
  !--recur over sub-branches
  if (associated (br % sub_branch)) then

    do i = 1, br % n_sub_branch
      call fill_nrespt_br (br % sub_branch(i))
    end do
    
  end if

end if

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine fill_nrespt_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_nrespt (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_nrespt'

integer :: i, j, k
integer :: n

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

!--possible to loop over bounding box?
n = 0

do k = 1, nz-1
  do j = 1, ny
    do i = 1, nx

      $if ($DEBUG)
      if (DEBUG) call mesg (sub_name, 'testing (i,j,k)=', (/ i, j, k /))
      $endif
      
      if (brindex(i, j, k) == br % ident) n = n + 1
      
    end do
  end do
end do

$if ($MPI)

  if (allocated (nrespt)) then  !--more checks like this needed
  
    nrespt(br % ident) = n  !--just partial sum, see fill_nrespt_ta
    
  else
  
    call error (sub_name, 'nrespt not allocated')
    
  end if

$else

  br % nrespt = n

$endif

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine fill_nrespt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--only required for non-Germano CD determination (for continued runs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_ftot (fname)
implicit none

character (*), intent (in) :: fname

character (*), parameter :: sub_name = mod_name // '.output_ftot'

integer, parameter :: lun = 1

integer :: i

logical :: opn, exst

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

inquire (unit=lun, opened=opn, exist=exst)
if (.not. exst) then
  call error (sub_name, 'lun does not exist')
else if (opn) then
  call error (sub_name, 'lun is already open') 
end if

inquire (file=fname, opened=opn)
if (opn) call error (sub_name, 'file ' // trim (fname) // ' already open')

!open (lun, file=fname, action='write', position='rewind')
open (lun, file=fname, action='write', form='unformatted', position='rewind')

do i = 1, n_tree
  call output_ftot_br (tree_array(i) % trunk, lun)
end do

close (lun)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine output_ftot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is already connected to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine output_ftot_br (br, lun)
implicit none

type (branch_type), intent (in) :: br
integer, intent (in) :: lun

character (*), parameter :: sub_name = mod_name // '.output_ftot_br'
character (*), parameter :: fmt = '(i0,3(1x,es12.5))'

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

!write (lun, fmt) br % ident, br % ftot
write (lun) br % ident, br % ftot

if (associated (br % sub_branch)) then

  do i = 1, br % n_sub_branch
    call output_ftot_br (br % sub_branch(i), lun)
  end do
  
end if

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine output_ftot_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_ftot (fname)
implicit none

character (*), intent (in) :: fname

character (*), parameter :: sub_name = mod_name // '.read_ftot'

integer, parameter :: lun = 1

integer :: i

logical :: opn, exst

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

inquire (unit=lun, opened=opn, exist=exst)
if (.not. exst) then
  call error (sub_name, 'lun does not exist')
else if (opn) then
  call error (sub_name, 'lun is already open') 
end if

inquire (file=fname, opened=opn)
if (opn) call error (sub_name, 'file ' // trim (fname) // ' already open')

!open (lun, file=fname, action='read', position='rewind')
open (lun, file=fname, action='read', form='unformatted', position='rewind')

do i = 1, n_tree
  call read_ftot_br (tree_array(i) % trunk, lun)
end do

close (lun)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine read_ftot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine read_ftot_br (br, lun)
implicit none

type (branch_type), intent (inout) :: br
integer :: lun

character (*), parameter :: sub_name = mod_name // '.read_ftot_br'

$if ($DEBUG)
logical, parameter :: DEBUG = .true.
$endif

integer :: i
integer :: ident

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

!read (lun, *) ident, br % ftot
read (lun) ident, br % ftot

if (ident /= br %ident) call error (sub_name, 'ident /= br % ident')

$if ($DEBUG)
if (DEBUG) then
  call mesg (sub_name, 'br % ident =', br % ident)
  call mesg (sub_name, 'br % ftot =', br % ftot)
end if
$endif
if (associated (br % sub_branch)) then

  do i = 1, br % n_sub_branch
    call read_ftot_br (br % sub_branch(i), lun)
  end do
  
end if

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine read_ftot_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--initializes things that cannot be initialized properly in the
!  main call to the trees module:
!  an example is (fx, fy, fz), which will modified by level set stuff
!  before we get to the tree module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trees_ls_init ()
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_ls_init'

!---------------------------------------------------------------------

if (init_from_file) call read_fxyz ('fxyz.last.out')

end subroutine trees_ls_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--writes the final fcoeff.last.dat, fdist.last.dat, ftot.last.dat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trees_ls_finalize ()
implicit none

character (*), parameter :: sub_name = mod_name // '.trees_ls_finalize'

!---------------------------------------------------------------------

!--unformatted
open (1, file=file_fcoeff_last, action='write', form='unformatted',  &
      position='rewind')
write (1) fcoeff
close (1)

!--write out num, den for fcoeff and fdist
!  (regardless of whether we are using them)
open (1, file=file_num_den_fcoeff_last, action='write',  &
      form='unformatted', position='rewind')
write (1) num_fcoeff, den_fcoeff
close (1)

open (1, file=file_num_den_fdist_last, action='write',  &
      form='unformatted', position='rewind')
write (1) num_fdist, den_fdist
close (1)

call output_fdist_wksp (file_fdist_last)  !--now unformatted

call output_ftot (file_ftot_last)

call write_fxyz ('fxyz.last.out')  !--.out since unformatted

end subroutine trees_ls_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_fxyz (fname)
use immersedbc, only : fx, fy, fz
implicit none

character (*), intent (in) :: fname

character (*), parameter :: sub_name = mod_name // '.write_fxyz'

character (128) :: file
character (32) :: suff

!---------------------------------------------------------------------

$if ($MPI)
  write (suff, '(".c",i0)') coord
  file = trim (fname) // trim (suff)
$else
  file = trim (fname)
$endif

open (1, file=file, action='write', position='rewind',  &
      form='unformatted')
write (1) fx, fy, fz
close (1)

end subroutine write_fxyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_fxyz (fname)
use immersedbc, only : fx, fy, fz
implicit none

character (*), intent (in) :: fname

character (*), parameter :: sub_name = mod_name // '.read_fxyz'

character (128) :: file
character (32) :: suff

!---------------------------------------------------------------------

$if ($MPI)
  write (suff, '(".c",i0)') coord
  file = trim (fname) // trim (suff)
$else
  file = trim (fname)
$endif

open (1, file=file, action='read', position='rewind',  &
      form='unformatted')
read (1) fx, fy, fz
close (1)

end subroutine read_fxyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--only for MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
$if ($MPI)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--only needed when "reduce arrays" are in use
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine reduce_array_init ()
  use param, only : iBOGUS
  implicit none

  character (*), parameter :: sub_name = mod_name // '.reduce_array_init'

  integer :: i

  !-------------------------------------------------------------------

  !--count number of branches
  nbr = 0
  
  do i = 1, n_tree
    call calc_nbr (tree_array(i) % trunk, nbr)
  end do
  !--now we have nbr
  
  allocate (nvelscale(nbr), nvelscale_sum(nbr))
  allocate (nrespt(nbr), nrespt_sum(nbr))
  allocate (resf(nd, nbr), resf_sum(nd, nbr))
  allocate (velscale(nd, nbr), velscale_sum(nd, nbr))

  nvelscale = iBOGUS
  nvelscale_sum = iBOGUS
  nrespt = iBOGUS
  nrespt_sum = iBOGUS
  resf = BOGUS
  resf_sum = BOGUS
  velscale = BOGUS
  velscale_sum = BOGUS

  if ( use_local_vel ) then
      allocate ( fnorm(nd, nfcoeff, nbr), fnorm_sum(nd, nfcoeff, nbr) )
      fnorm = BOGUS
      fnorm_sum = BOGUS
  end if

  end subroutine reduce_array_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine calc_nbr (br, n)
  implicit none

  type (branch_type), intent (in) :: br
  integer, intent (inout) :: n

  integer :: i

  !-------------------------------------------------------------------

  n = n + 1

  if (associated (br % sub_branch)) then
    do i = 1, br % n_sub_branch
      call calc_nbr (br % sub_branch(i), n)
    end do
  end if

  end subroutine calc_nbr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--this is quite limited for now, maybe more general later if needed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine repack_br (br, array, tag)
  implicit none

  type (branch_type), intent (inout) :: br
  real (rp), intent (in) :: array(:, :)  !--assume this is nd X nbr
  character (*), intent (in) :: tag

  character (*), parameter :: sub_name = mod_name // '.repack_br'

  integer :: i

  !-------------------------------------------------------------------

  select case (tag)
    case ('velscale')
      br % velscale = array(:, br % ident)
    case ('resf')
      br % resf = array(:, br % ident)
    case default
      call error (sub_name, 'invalid tag')
  end select
  
  if (associated (br % sub_branch)) then
  
    do i = 1, br % n_sub_branch
      call repack_br (br % sub_branch(i), array, tag)
    end do
    
  end if
  
  end subroutine repack_br

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--this is quite limited for now, maybe more general later if needed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine repack_rank3_br (br, array, tag)
  implicit none

  type (branch_type), intent (inout) :: br
  real (rp), intent (in) :: array(:, :, :)  !--assume nd X nfcoeff X nbr
  character (*), intent (in) :: tag

  character (*), parameter :: sub_name = mod_name // '.repack_rank3_br'

  integer :: i

  !-------------------------------------------------------------------

  select case (tag)
    case ('fnorm')
      br % fnorm = array(:, :, br % ident)
    case default
      call error (sub_name, 'invalid tag')
  end select
  
  if (associated (br % sub_branch)) then
  
    do i = 1, br % n_sub_branch
      call repack_rank3_br (br % sub_branch(i), array, tag)
    end do
    
  end if
  
  end subroutine repack_rank3_br
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--integer version of repack_br above, except 1-d array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine irepack_br (br, array, tag)
  implicit none

  type (branch_type), intent (inout) :: br
  integer, intent (in) :: array(:)  !--assume this is nbr
  character (*), intent (in) :: tag

  character (*), parameter :: sub_name = mod_name // '.irepack_br'

  integer :: i

  !-------------------------------------------------------------------

  select case (tag)
    case ('nrespt')
      br % nrespt = array (br % ident)
    case default
      call error (sub_name, 'invalid tag')
  end select
  
  if (associated (br % sub_branch)) then
  
    do i = 1, br % n_sub_branch
      call irepack_br (br % sub_branch(i), array, tag)
    end do
    
  end if
  
  end subroutine irepack_br
$endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--MPI: all processes must call this for the mpi_reduce to work
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_Uinf ( Uinf )
use sim_param, only : u, v !, w
implicit none

real (rp), intent (out) :: Uinf(nd)

character (*), parameter :: sub_name = mod_name // '.calc_Uinf'

character (64) :: sn

$if ($MPI)
  real (rp) :: Uinf_sum(nd)
$endif

!---------------------------------------------------------------------

sn = trim ( sub_name ) // trim ( chcoord )

$if ($VERBOSE)
call enter_sub ( sn )
$endif

if ( nd /= 3 ) call error ( sn, 'expecting nd=3' )

Uinf(1) = sum ( u(1, 1:ny, 1:nz-1) ) / ( ny * (nz - 1) )
Uinf(2) = sum ( v(1:nx, 1, 1:nz-1) ) / ( nx * (nz - 1) )
Uinf(3) = 0.0_rp  !--with current BC, this should always be true
    !--if change bc so that this is no longer true, then the MPI part
    !  for Uinf(3) must be different from others, since the averaging plane
    !  is a z-plane, that does not cut across the processes.  Note also that
    !  it does not actually need to be reduced (avg can be calculated over
    !  any z-plane, at least when periodic BC are used)
    
$if ($MPI)
  call mpi_allreduce ( Uinf, Uinf_sum, nd, MPI_RPREC, MPI_SUM,  &
                       comm, ierr )
  Uinf = Uinf_sum / nproc
$endif
    
$if ($VERBOSE)
call exit_sub ( sn )
$endif
    
end subroutine calc_Uinf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--MPI: all processes must call this for the mpi_reduce to work
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_globalCD ()
use param, only : L_x, L_y, L_z, dt  !--see top of this module for more
use sim_param, only : u, txz  !--txz for wall stress @ bottom of domain
implicit none

character (*), parameter :: sub_name = mod_name // '.calc_globalCD'
character (*), parameter :: foutput = 'globalCD.dat'

integer, parameter :: lun = 1

integer :: i, j

logical :: exst, opn

real (rp) :: gCD, gCD_r, gCD_u, gCD_w
real (rp) :: A
real (rp) :: Uinf
$if ($MPI)
    real (rp) :: Uinf_sum
$endif
real (rp) :: f, f_r, f_u, f_w

!---------------------------------------------------------------------

Uinf = sum (u(1, :, 1:nz-1)) / (ny * (nz - 1))
    
$if ($MPI)
  call mpi_reduce (Uinf, Uinf_sum, 1, MPI_RPREC, MPI_SUM,  &
                   rank_of_coord(0), comm, ierr)
  Uinf = Uinf_sum / nproc
$endif

if (USE_MPI .and. coord /= 0) return
!--only coord 0 from here down

!--domain wall area (friction coefficient)
A = L_x * L_y

f = 0._rp
f_r = 0._rp
f_u = 0._rp

do i = 1, n_tree
  !--minus signs: force on tree
  f = f - tree_array(i) % trunk % ftot(1)
  f_r = f_r - tree_array(i) % trunk % resftot(1)
  f_u = f_u - tree_array(i) % trunk % unresftot(1)
end do

!--calculate force on bottom wall: we may need to add this for global CD
f_w = 0._rp

do j = 1, ny
  do i = 1, nx

    !--phi here is u-node, but txz at w-node (what else to do?)
    if (phi(i, j, 1) > 0._rp) f_w = f_w - txz(i, j, 1)

  end do
end do

f_w = f_w * dx * dy

!--calc global drag coefficients
gCD = (f + f_w) / (0.5_rp * Uinf**2 * A)  !--include wall here
gCD_r = f_r / (0.5_rp * Uinf**2 * A)
gCD_u = f_u / (0.5_rp * Uinf**2 * A)
gCD_w = f_w / (0.5_rp * Uinf**2 * A)  !--same A as others, not A-A_tree

!--append to existing datafile
inquire (unit=lun, exist=exst, opened=opn)

if (opn) call error (sub_name, 'unit is already open, unit=', lun)
if (.not. exst) call error (sub_name, 'unit does not exist, unit=', lun)

!--open and close every time to flush
open (lun, file=foutput, action='write', position='append')

write (lun, '(5(es13.6,1x))') dt * (jt_total-1), gCD, gCD_r, gCD_u, gCD_w

close (lun)

end subroutine calc_globalCD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--may want to change this to unformatted i/o, if so, make sure 
!  dependent routines are also changed (e.g. where file is read)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_fdist_wksp (foutput)
implicit none

character (*), intent (in) :: foutput

character (*), parameter :: sub_name = mod_name // '.output_fdist_wksp'

integer, parameter :: lun = 1

character (64) :: rfmt

logical :: opn, exst

!---------------------------------------------------------------------

!!--calculate the number format
!!--idea is to have precision () total digits:
!!  1 leading (ES) plus precision () - 1
!!--for ESw.d format, need w >= d + 8
!write (rfmt, "('es',i0,'.',i0)") precision (0._rp) + 7, precision (0._rp) - 1

inquire (unit=lun, exist=exst, opened=opn)
if (.not. exst) call error (sub_name, 'lun does not exist')
if (opn) call error (sub_name, 'lun already open')

open (lun, file=trim (foutput), action='write', form='unformatted',  &
      position='rewind')

write (lun) fdist_wksp

close (lun)

end subroutine output_fdist_wksp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fnorm_ta ()
implicit none

character (*), parameter :: sub_name = mod_name // '.calc_fnorm_ta'

integer :: i

logical, save :: do_init_global_fmask = .true.

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sub_name )
$endif

if ( use_global_fmask ) then
    if ( do_init_global_fmask ) then
        call read_global_fmask ()
        do_init_global_fmask = .false.
    end if
end if

do i = 1, n_tree
    call calc_fnorm_br ( tree_array( i ) % trunk )
end do

$if ($MPI)
    !--fnorm, fnorm_sum should be the same size
    call mpi_allreduce ( fnorm, fnorm_sum, size (fnorm), MPI_RPREC, MPI_SUM,  &
                         comm, ierr )

    do i = 1, n_tree
        call repack_rank3_br ( tree_array(i) % trunk, fnorm_sum, 'fnorm' )
    end do
$endif

$if ($VERBOSE)
call exit_sub ( sub_name )
$endif

end subroutine calc_fnorm_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine calc_fnorm_br ( br )
implicit none

type ( branch_type ), intent ( in out ) :: br

character (*), parameter :: sub_name = mod_name // '.calc_fnorm_br'

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sub_name )
$endif

if ( .not. br % resolved ) then

    call calc_fnorm ( br )

else

    if (.not. associated ( br % sub_branch ) ) then
        call error ( sub_name, 'expecting sub_branch to be associated' )
    end if

    do i = 1, br % n_sub_branch
        call calc_fnorm_br ( br % sub_branch( i ) )
    end do
    
end if

$if ($VERBOSE)
call exit_sub ( sub_name )
$endif
    
end subroutine calc_fnorm_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--expects fdist is valid here
!--this needs to be fully consistent with apply_unresf_br
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fnorm ( br )
use sim_param, only : u, v, w
implicit none

type (branch_type), intent (in out) :: br

character (*), parameter :: sn = mod_name // '.calc_fnorm'

real (rp), parameter :: eps = epsilon (eps)

integer :: ipt
integer :: i, j, k
integer :: m

real (rp) :: fmask
real (rp) :: sigma
real (rp) :: fn(nd, nfcoeff)
real (rp) :: vel(nd)
real (rp) :: dir(nd, nfcoeff)
real (rp) :: vv(nd, nfcoeff)
real (rp) :: vvmag(nfcoeff)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sn )
$endif

if ( .not. associated ( br % bboxpt ) ) call bboxpt_init ( br )
    !--after this all bbox will be associated
    !--some procs may have 0-sized arrays

!if ( DEBUG ) call mesg ( sn, 'br % nbboxpt=', br % nbboxpt )

select case ( flocal )
case ( 'fmodel' )
    fn = 0.0_rp
case ( 'unresf_dir' )
    fn = BOGUS
    fn(1, 1) = 0.0_rp
case default
    call error (sn, 'illegal flocal: ' // flocal)
end select

bboxpt_loop: do ipt = 1, br % nbboxpt

    i = br % bboxpt( 1, ipt )
    j = br % bboxpt( 2, ipt )
    k = br % bboxpt( 3, ipt )

    if (brindex(i, j, k) /= -1) cycle  !--(-1) is the "no branch" value
        !--check consistency with apply_unresf_br here

    !--(k+1) is in bounds since max(k) = nz - 1 (see bboxpt_init)
    !--MPI: w(k+1) may be BOGUS
    vel = (/ u( i, j, k ), v( i, j, k ),                 &
             0.5_rp * ( w( i, j, k ) + w( i, j, k+1 ) ) /)

    if ( use_global_fmask ) then
        fmask = global_fmask( i, j, k )
    else
        fmask = local_fmask( br, i, j, k ) 
        if ( fmask <= epsilon (fmask) ) cycle bboxpt_loop
            !--check consistency with apply_unresf_br
    end if

    !--check consistency with calc_flocal here, this is error-prone:
    !  better to unify this block with calc_flocal
    !  (actually tried a quickie, but failed)
    !--check for validity of flocal is just before loop, above
    select case ( flocal )
    case ( 'fmodel' )
        select case ( fmodel )
        case ( 'd', 'd_germano' )
            !call def_d_dir ( br )
            call vel_d ( vel, vvmag, vv )  !--not dir needed
        case ( 'dls' )
            call def_dls_dir ( br, vel, dir )
            call vel_dls ( vel, dir, vvmag, vv )
        case ( 'nba' )
            call def_nba_dir ( br, vel, dir )
            call vel_nba ( vel, dir, vvmag, vv )
        case default
            call error ( sn, 'invalid fmodel: ' // trim (fmodel) )
        end select

        !--this needs to be a sum similar to that in fill_unres_br
        do m = 1, nfcoeff
            fn(:, m) = fn(:, m) + vv(:, m) * vvmag(m) * fmask
        end do
    case ( 'unresf_dir' )
        if ( mag ( br % unresf ) > eps ) then
            !--could do > eps test outside of loop
            if ( unresf_dir_proj ) then 
                !--note minus sign
                fn(1, 1) = fn(1, 1) -                                        &
                           ( mag (vel)  * dot_product ( vel, br % unresf ) * &
                             fmask / mag ( br % unresf ) )
            else
                !--note plus sign
                sigma = 1.0_rp
                if ( unresf_clip ) then
                    if ( dot_product ( vel, br % unresf ) > 0.0_rp ) then
                        sigma = 0.0_rp
                    end if
                end if
                fn(1, 1) = fn(1, 1) + sigma * fmask * mag (vel)**2
            end if
        else
            fn(1, 1) = 0.0_rp  !--make sure do not use this value later
                !--fdist_coeff will be set to zero
                !-- check consistency with local_vel_fdist_coeff_unresf_dir
        end if
    end select

    !if (DEBUG) then
    !    call mesg (sn, '(i,j,k)=', (/ i, j, k /) )
    !    call mesg (sn, 'vel=', vel )
    !    call mesg (sn, 'vv=', pack (vv, .true.) )
    !    call mesg (sn, 'vvmag=', vvmag)
    !    call mesg (sn, 'fmask=', fmask)
    !    call mesg (sn, 'dir=', pack (dir, .true.) )
    !end if

end do bboxpt_loop

!--current convention does not include dV in fnorm at all
$if ($DEBUG)
if (DEBUG) then
    call mesg (sn, 'br % ident = ', br % ident)
    call mesg ( sn, 'fn = ', pack (fn, .true.) )
end if
$endif
$if ($MPI)

    if ( allocated ( fnorm ) ) then
        fnorm(:, :, br % ident ) = fn
    else
        call error ( sn, 'expecting fnorm allocated' )
    end if

$else

    br % fnorm = fn

    !--this applies only to 'd' for 'd_germano' case
    !--generally need to check if determinant of matrix is non-zero
    !  but have not implemented this yet, instead this will show up
    !  when the system is solved
    !do m = 1, nfcoeff
    !    if ( mag ( fn(:, m) ) < epsilon ( 1._rp ) )  &
    !        call error ( sn, ' at branch % ident=', br % ident,  &
    !                         ': small mag (fn) = ', mag (fn) )
    !    end if
    !end do

$endif

$if ($VERBOSE)
call exit_sub ( sn )
$endif

end subroutine calc_fnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_unresf_ta ()
use sim_param, only : u, v, w  !--for debug only
use immersedbc, only : fx, fy, fz  !--for debug only
implicit none

character (*), parameter :: sub_name = mod_name // '.apply_unresf_ta'

integer :: i
integer :: j, k

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sub_name )
$endif

!--this is an experiment
!--now we have fdist, clear the old unresolved force field
do k = 1, nz-1  !--not sure if want nz or nz-1 here
    do j = 1, ny
        do i = 1, nx
        
            if ( brindex( i, j, k ) == -1 ) then
                !--careful of u/w nodes
                fx( i, j, k ) = 0._rp
                fy( i, j, k ) = 0._rp
                fz( i, j, k ) = 0._rp
            end if
            
        end do
    end do
end do

$if ($DEBUG)
if ( DEBUG ) then
    call DEBUG_write ( fx, 'apply_unresf_ta.fx.a' )
    call DEBUG_write ( fy, 'apply_unresf_ta.fy.a' )
    call DEBUG_write ( fz, 'apply_unresf_ta.fz.a' )
    call DEBUG_write ( u(:, :, 1:nz), 'apply_unresf_ta.u.a')
    call DEBUG_write ( v(:, :, 1:nz), 'apply_unresf_ta.v.a')
    call DEBUG_write ( w(:, :, 1:nz), 'apply_unresf_ta.w.a')
    call DEBUG_write ( global_fmask(:, :, 1:nz),  &
                       'apply_unresf_ta.global_fmask.a')
end if
$endif

if ( use_local_vel ) call calc_fnorm_ta ()  !--inits also global_fmask

do i = 1, n_tree
    call apply_unresf ( tree_array( i ) % trunk )
end do

$if ($DEBUG)
if ( DEBUG ) then
    call DEBUG_write ( fx, 'apply_unresf_ta.fx.b' )
    call DEBUG_write ( fy, 'apply_unresf_ta.fy.b' )
    call DEBUG_write ( fz, 'apply_unresf_ta.fz.b' )
    call DEBUG_write ( u(:, :, 1:nz), 'apply_unresf_ta.u.b')
    call DEBUG_write ( v(:, :, 1:nz), 'apply_unresf_ta.v.b')
    call DEBUG_write ( w(:, :, 1:nz), 'apply_unresf_ta.w.b')
    call DEBUG_write ( global_fmask(:, :, 1:nz),  &
                       'apply_unresf_ta.global_fmask.b')
end if
$endif

$if ($VERBOSE)
call exit_sub ( sub_name )
$endif

end subroutine apply_unresf_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine apply_unresf (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.apply_unresf'

integer :: i

!---------------------------------------------------------------------

if (br % gen == tree_array(br % itree) % n_gen) then

  if (br % resolved) call error (sub_name, 'expecting unresolved branch')

  call apply_unresf_br (br)

else

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if
  
  do i = 1, br % n_sub_branch
    call apply_unresf (br % sub_branch(i))
  end do

end if
  
end subroutine apply_unresf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--sets fx, fy, fz within the bounding box
!--assumes no overlap of SGS forcing regions (i.e. that masks do not
!  overlap), since it overwrites
!--need to be careful not to overwrite fx, fy, fz if that point is
!  inside a resolved branch--this is likely to be an issue near the
!  base of the branch where it meets a resolved branch
!--needs to be fully consistent with calc_fnorm if 'use_local_vel' is true
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_unresf_br (br)
use immersedbc, only : fx, fy, fz
use sim_param, only : u, v, w
implicit none

type (branch_type), intent (inout) :: br  !--out since sets fdist_coeff

character (*), parameter :: sub_name = mod_name // '.apply_unresf_br'

logical, parameter :: allow_overlap = .true.
!logical, parameter :: DEBUG = .true.
!logical, parameter :: local_vel_drag_only = .true.

character (64) :: sn

integer :: i, j, k
integer :: m
integer :: ipt

real (rp) :: dV
real (rp) :: fmask
real (rp) :: f(nd)
real (rp) :: ftot(nd)
real (rp) :: vel(nd)
real (rp) :: x(nd)
real (rp) :: xp(nd)

!---------------------------------------------------------------------

sn = trim ( sub_name ) // trim ( chcoord )

$if ($DEBUG)
if ( DEBUG ) ftot = 0._rp
$endif
dV = dx * dy * dz

if ( use_local_vel ) then
    br % fdist_coeff = local_vel_fdist_coeff ( br )
        !--code in transition here, may make this function into a subroutine
        !  and have it modify its argument
end if

$if ($DEBUG)
if ( DEBUG ) then
    call mesg ( sn, 'br % ident=', br % ident )
    call mesg ( sn, 'br % nbboxpt=', br % nbboxpt )
end if
$endif
bboxpt_loop: do ipt = 1, br % nbboxpt

    i = br % bboxpt( 1, ipt )
    j = br % bboxpt( 2, ipt )
    k = br % bboxpt( 3, ipt )

    if ( brindex(i, j, k) /= -1 ) cycle  !--(-1) is the "no branch" value
        !--adding this means that we are not adding the full force
        !  instead we are skipping some points
        !--check consistency with calc_fnorm here
        
    if ( use_local_vel ) then
    
        !--(k+1) is in bounds since max(k) = nz - 1 (see bboxpt_init)
        !--MPI: w(k+1) may be BOGUS
        vel = (/ u(i, j, k), v(i, j, k),                 &
                 0.5_rp * (w(i, j, k) + w(i, j, k + 1)) /)

        $if ($MPI)
            if ((coord == nproc - 1) .and. (k == nz - 1)) then
                call error (sn, 'bboxpt trying to access w@k=nz' //          &
                            n_l // 'domain is probably too small to fit bbox')
            end if
        $endif

        if ( use_global_fmask ) then
            fmask = global_fmask( i, j, k )
        else
            fmask = local_fmask ( br, i, j, k )
            if ( fmask <= epsilon (fmask) ) cycle bboxpt_loop
                !--so we do not overwrite valid force from another
                !  branch with 0 force from this branch
                !--check consistency with calc_fnorm
        end if

        f = calc_flocal (br, vel, br % fdist_coeff, fmask)
        
    else  !--use_local_vel
    
        fmask = local_fmask ( br, i, j, k )
        if ( fmask <= epsilon (fmask) ) cycle bboxpt_loop
            !--so we do not overwrite valid force from another
            !  branch with 0 force from this branch

        !--not sure about u/w node stuff
        !--overwrites existing distribution in the bbox
        !--this assumes no overlap w/ other branches (overwrites f)
        !  if want to allow overlap, set allow_overlap = T
        !  but f needs to be cleared within the mask region once per time step

        f = fmask * ( br % unresf )
            !--additional normalization may be required if fdist is in the
            !  range [0,1]

    end if

    if (allow_overlap) then
        !--make sure this is consistent with apply_unresf_ta:
        !  forces w/ brindex = -1 should be zeroed there each t-step

        fx(i, j, k) = fx(i, j, k) + f(1)
        fy(i, j, k) = fy(i, j, k) + f(2)
        fz(i, j, k) = fz(i, j, k) + f(3)
        
    else
    
        fx(i, j, k) = f(1)
        fy(i, j, k) = f(2)
        fz(i, j, k) = f(3)

    end if

    if (DEBUG) ftot = ftot + (/ fx(i, j, k), fy(i, j, k), fz(i, j, k) /)

end do bboxpt_loop

$if ($DEBUG)
if (DEBUG) then
  !--measure how much force we are actually applying after interpolation
  !  and whatever other conditions above
  call mesg (sn, 'br % ident =', br % ident)
  call mesg (sn, 'br % ftot =', br % ftot)
  ftot = ftot * (dx * dy * dz)
  call mesg (sn, 'measured ftot =', ftot)
end if
$endif
end subroutine apply_unresf_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--coeff can be br % fdist_coeff, but does not have to be
!  e.g. if coeff = 1, then this routine can be used to calculate
!  the normalization fnorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function calc_flocal (br, vel, coeff, fmask)  result (f)
implicit none

real (rp) :: f(nd)

type (branch_type), intent (in) :: br
real (rp), intent (in) :: vel(nd)
real (rp), intent (in) :: coeff(nfcoeff)
real (rp), intent (in) :: fmask

character (*), parameter :: sn = mod_name // '.calc_flocal'

integer :: m

real (rp) :: sigma
real (rp) :: dir(nd, nfcoeff)
real (rp) :: vv(nd, nfcoeff)
real (rp) :: vvmag(nfcoeff)
real (rp ) :: uf_dir(nd)

!---------------------------------------------------------------------

select case ( flocal )
case ( 'fmodel' )
    select case ( fmodel )
    case ( 'd', 'd_germano' )
        !call def_d_dir ( br )
        call vel_d ( vel, vvmag, vv )  !--not dir needed
    case ( 'dls' )
        call def_dls_dir ( br, vel, dir )
        call vel_dls ( vel, dir, vvmag, vv )
    case ( 'nba' )
        call def_nba_dir ( br, vel, dir )
        call vel_nba ( vel, dir, vvmag, vv )
    case default
        call error ( sn, 'invalid fmodel: ' // trim (fmodel) )
    end select

    !--this needs to be a sum similar to that in fill_unres_br
    f = 0.0_rp
    do m = 1, nfcoeff
        f = f + coeff(m) * vv(:, m) * vvmag(m) * fmask
    end do
case ( 'unresf_dir' )
    if ( unresf_dir_proj ) then
        if ( mag ( br % unresf ) > epsilon ( 0.0_rp ) ) then
            !--note minus sign
            f = -coeff(1) * mag (vel) * dot_product ( vel, br % unresf ) *  &
                fmask * ( br % unresf ) / mag ( br % unresf )
        else
            f = 0.0_rp
        end if
    else
        !--note no minus sign   
        sigma = 1.0_rp
        if ( unresf_clip ) then
            if ( dot_product ( vel, br % unresf ) > 0.0_rp ) then
                sigma = 0.0_rp
            end if
        end if
        f = coeff(1) * sigma * ( mag (vel)**2 ) * fmask * ( br % unresf )
    end if
case default
    call error (sn, 'illegal flocal: ' // trim (flocal) )
end select

end function calc_flocal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--the sign of these coefficients is important because positive
!  values will cause exponential increases in local velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function local_vel_fdist_coeff ( br )  result ( coeff )
implicit none

real (rp) :: coeff(nfcoeff)

type ( branch_type ), intent (in) :: br

character (*), parameter :: sn = mod_name // '.local_vel_fdist_coeff'

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sn)
$endif

select case ( flocal )
case ( 'fmodel' )
    coeff = local_vel_fdist_coeff_fmodel (br)
case ( 'unresf_dir' )
    coeff = local_vel_fdist_coeff_unresf_dir (br)
case default
  call error ( sn, 'illegal flocal: ' // flocal )
end select

$if ($VERBOSE)
 call exit_sub (sn)
$endif

end function local_vel_fdist_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function local_vel_fdist_coeff_unresf_dir ( br )  result ( coeff )
implicit none

real (rp) :: coeff(nfcoeff)

type (branch_type), intent (in) :: br

character (*), parameter :: sn = mod_name //                       &
                                 '.local_vel_fdist_coeff_unresf_dir'

real (rp), parameter :: eps = epsilon (eps)

real (rp) :: dV
real (rp) :: Q

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sn)
$endif

dV = dx * dy * dz
Q = br % fnorm(1, 1)  !--all other components of fnorm are BOGUS

!--fnorm = sum ( |u|^2 * chi ), so should always be > 0
!  if unresf_dir_proj is true, then fnorm should still be > 0, unless
!  recirculation is dominating the sum (and that is BAD news)
!--with this convention, coeff(1) is >= 0
if ( Q < eps ) then
    !call error (sn, 'fnorm is too small: ', Q)
    coeff(1) = 0.0_rp
else
    coeff(1) = 1.0_rp / ( Q * dV )
end if
coeff(2:) = BOGUS

$if ($VERBOSE)
 call exit_sub (sn)
$endif

end function local_vel_fdist_coeff_unresf_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function local_vel_fdist_coeff_fmodel ( br )  result ( coeff )
use linear_simple, only : solve_linear
implicit none

real (rp) :: coeff(nfcoeff)

type ( branch_type ), intent (in) :: br

character (*), parameter :: sn = mod_name // '.local_vel_fdist_coeff_fmodel'

!logical, parameter :: DEBUG = .false.
logical, parameter :: clip_fdist_coeff = .true.
logical, parameter :: limit_ratio = .true.

real (rp), parameter :: ratio_max = 2.0_rp
real (rp), parameter :: eps = 0.05_rp

integer :: i, j

real (rp) :: A(nfcoeff, nfcoeff)
!real (rp) :: A2(2, 2)  !--experimental
real (rp) :: b(nfcoeff)
!real (rp) :: b2(2)
!real (rp) :: coeff2(2)
real (rp) :: det
real (rp) :: dV
real (rp) :: x(nfcoeff)
real (rp) :: Q(nd, nfcoeff)
real (rp) :: F(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sn)
$endif

dV = dx * dy * dz

!--shorthand to simplify formulas
Q = br % fnorm
F = br % unresf

$if ($DEBUG)
if (DEBUG) then
    call mesg (sn, 'Q = ', pack (Q, .true.) )
    call mesg (sn, 'F = ', F )
end if
$endif
select case ( fmodel )
case ( 'd', 'd_germano' )
    coeff(1) = dot_product ( Q(:, 1), F ) / ( mag (Q(:, 1))**2 * dV )
case ( 'dls' )

    call warn (sn, 'have not considered linear system in detail: ' // n_l //  &
                   'may be singular or may simplify' )

    !--can unify nba and dls if we use linear system approch
    do j = 1, nfcoeff
        do i = 1, nfcoeff
            A(i, j) = dot_product ( Q(:, i), Q(:, j) )
        end do
    end do

    !--could write as
    !  matmul ( transpose ( br % fnorm ), br % unresf )
    do i = 1, nfcoeff
        b(i) = dot_product ( F, Q(:, i) )
    end do

    call solve_linear ( A, b, coeff )

case ( 'nba' )
    !!--linear system: Ax = b
    !!--issue: this system is singular if include b-component
    !do j = 1, nfcoeff
    !    do i = 1, nfcoeff
    !        A(i, j) = dot_product ( Q(:, i), Q(:, j) )
    !    end do
    !end do

    !!--could write as
    !!  matmul ( transpose ( br % fnorm ), br % unresf )
    !do i = 1, nfcoeff
    !    b(i) = dot_product ( F, Q(:, i) ) / dV
    !end do

    !call mesg (sn, 'about to call solve_linear' )
    !call solve_linear ( A, b, coeff )
    !call mesg (sn, 'solve_linear gives coeff = ', coeff )

    !--since a is constant for a branch, and n is always normal to a,
    !  the terms decouple, i.e., Q(:, 1) .dot. Q(:, 3) = 0
    coeff(1) = dot_product ( Q(:, 1), F ) / ( mag (Q(:, 1))**2 * dV )
    coeff(2) = 0.0_rp
    coeff(3) = dot_product ( Q(:, 3), F ) / ( mag (Q(:, 3))**2 * dV )

    !!--neglect a component if it is too small (relative to the other)
    !if ( mag (Q(:, 1)) / mag (Q(:, 3)) < eps ) then
    !    !--neglect Q(:, 1)
    !    coeff(1) = 0.0_rp
    !    coeff(2) = 0.0_rp
    !    coeff(3) = dot_product ( Q(:, 3), F ) / ( mag (Q(:, 3))**2 * dV )
    !else if ( mag (Q(:, 3)) / mag (Q(:, 1)) < eps ) then
    !    !--neglect Q(:, 3)
    !    coeff(1) = dot_product ( Q(:, 1), F ) / ( mag (Q(:, 1))**2 * dV )
    !    coeff(2) = 0.0_rp
    !    coeff(3) = 0.0_rp
    !else
    !    !!--system without b-component
    !    !A2(1, 1) = dot_product ( Q(:, 1), Q(:, 1) )
    !    !A2(1, 2) = dot_product ( Q(:, 1), Q(:, 3) )
    !    !A2(2, 1) = dot_product ( Q(:, 3), Q(:, 1) )
    !    !A2(2, 2) = dot_product ( Q(:, 3), Q(:, 3) )

    !    !b2(1) = dot_product ( F, Q(:, 1) ) / dV
    !    !b2(2) = dot_product ( F, Q(:, 3) ) / dV

    !    !call mesg (sn, 'about to call solve_linear' )
    !    !call mesg (sn, 'A2 = ', pack( A2, .true. ) )
    !    !call mesg (sn, 'b2 = ', b2 )
    !    !call solve_linear ( A2, b2, coeff2 )
    !    !coeff(1) = coeff2(1)
    !    !coeff(2) = 0.0_rp
    !    !coeff(3) = coeff2(2)
    !    !call mesg (sn, 'solve_linear gives coeff = ', coeff )

    !    !--direct solution, using fnorm(:, 2) = 0 identically
    !    !  (causing an underdetermined linear system)
    !    det = ( mag (Q(:, 1)) * mag (Q(:, 3)) )**2 -  &
    !          ( dot_product (Q(:, 1), Q(:, 3)) )**2

    !    if ( abs (det) >= epsilon (det) ) then
    !        coeff(1) = ( mag (Q(:, 3))**2 * dot_product( Q(:, 1), F ) -  &
    !                     dot_product( Q(:, 1), Q(:, 3) ) *               &
    !                       dot_product( Q(:, 3), F )                     &
    !                   ) / ( det * dV )
    !        coeff(2) = 0.0_rp
    !        coeff(3) = ( mag (Q(:, 1))**2 * dot_product( Q(:, 3), F ) -  &
    !                     dot_product( Q(:, 1), Q(:, 3) ) *               &
    !                       dot_product( Q(:, 1), F )                     &
    !                   ) / ( det * dV )
    !        !call mesg (sn, 'direct solution gives ', coeff )
    !    else
    !        coeff = 0.0_rp  !--any better idea?
    !    end if

    !end if
case default
    call error (sn, 'invalid force model: ' // fmodel )
end select

if ( clip_fdist_coeff ) then
    coeff = min ( coeff, 0.0_rp )
        !--min since we want to keep negative values here
end if

if ( limit_ratio ) then
    if ( fmodel == 'nba' ) then
        !--coeff should be negative here
        !--increase the smaller coefficient if the ratio of coefficients is
        !  either too large, or too small
        if ( abs (coeff(1)) > ratio_max * abs (coeff(3)) ) then
            coeff(3) = coeff(1) / ratio_max
        else if ( abs (coeff(3)) > ratio_max * abs (coeff(1)) ) then
            coeff(1) = coeff(3) / ratio_max
        end if
    end if
end if

!--check signs
if ( any ( coeff > 0.0_rp ) ) then
    call mesg (sn, 'br % abs_dir =', br % abs_dir )
    call mesg (sn, 'Q =', pack (Q, .true.) )
    call mesg (sn, 'F =', F )
    !call error (sn, 'positive coeff: ', coeff )
        !--for now error, but could also prescribe 0 instead?
    call mesg (sn, 'positive coeff: ', coeff )
        !--see what happens if we leave it alone
end if

$if ($VERBOSE)
 call exit_sub (sn)
$endif

end function local_vel_fdist_coeff_fmodel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function local_fmask ( br, i, j, k )
implicit none

real (rp) :: local_fmask

type (branch_type) , intent (in) :: br
integer, intent (in) :: i, j, k

character (*), parameter :: sub_name = mod_name // '.local_fmask'

real (rp) :: fdist
real (rp) :: x(nd), xp(nd)

!---------------------------------------------------------------------

!--calculate local fmask
x(1) = pt_of_grid (i, 1, 1)
x(2) = pt_of_grid (j, 2, 1)
x(3) = pt_of_grid (k, 3, 1)

!--get branch-local coords
x = x - br % x0

xp(1) = dot_product (x, br % x_hat)
xp(2) = dot_product (x, br % y_hat)
xp(3) = dot_product (x, br % z_hat)

!--interpolate the force distribution to use from fdist_wksp
!--this may be a problem, since it will not strictly conserve force
call interp_fdist ( xp, br % zone, fdist )

local_fmask = fdist

end function local_fmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--performs bilinear interpolation to get fdist at location x
!--x is branch-local coordinate
!--not conservative? i.e. interpolation introduces small errors
!  and the force applied is less than the force that we meant to
!  apply.  is this significant?
!  if it is significant, then we probably need to make a temp copy of
!  the interpolated version and normalize that before application
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_fdist (x, iz, fdist)
implicit none

real (rp), intent (in) :: x(nd)
integer :: iz
real (rp), intent (out) :: fdist

character (*), parameter :: sub_name = mod_name // '.interp_fdist'

integer :: i, j, k
integer :: i1, j1, k1
integer :: io, jo, ko

real (rp) :: x1, x2, x3
real (rp) :: w1, w2, w3, w4, w5, w6, w7, w8

!---------------------------------------------------------------------

io = wksp_size(1) / 2 + 1
jo = wksp_size(2) / 2 + 1
ko = 1

i = floor (x(1) / dx) + io
j = floor (x(2) / dy) + jo
k = floor (x(3) / dz) + ko

!--the >= here means that we are checking i+1, j+1, k+1 in bounds
if ( (i < 1) .or. (i >= wksp_size(1)) ) then
  call error (sub_name, 'i out of bounds, i =', i)
else if ( (j < 1) .or. (j >= wksp_size(2)) ) then
  call error (sub_name, 'j out of bounds, j =', j)
else if ( (k < 1) .or. (k >= wksp_size(3)) ) then
  call error (sub_name, 'k out of bounds, k =', k)
end if

i1 = i + 1
j1 = j + 1
k1 = k + 1

!--interpolation weights: since io, jo, ko are integers
x1 = modulo (x(1), dx) / dx
x2 = modulo (x(2), dy) / dy 
x3 = modulo (x(3), dz) / dz  !--no need to worry about u/w here

w1 = (1._rp - x1) * (1._rp - x2) * (1._rp - x3)
w2 = (    x1    ) * (1._rp - x2) * (1._rp - x3)
w3 = (1._rp - x1) * (    x2    ) * (1._rp - x3)
w4 = (    x1    ) * (    x2    ) * (1._rp - x3)
w5 = (1._rp - x1) * (1._rp - x2) * (    x3    )
w6 = (    x1    ) * (1._rp - x2) * (    x3    )
w7 = (1._rp - x1) * (    x2    ) * (    x3    )
w8 = (    x1    ) * (    x2    ) * (    x3    )

!--experiment: direct copy
!fdist = fdist_wksp(i, j, k, iz)

fdist = w1 * fdist_wksp(i, j , k , iz) + w2 * fdist_wksp(i1, j , k , iz) +  &
        w3 * fdist_wksp(i, j1, k , iz) + w4 * fdist_wksp(i1, j1, k , iz) +  &
        w5 * fdist_wksp(i, j , k1, iz) + w6 * fdist_wksp(i1, j , k1, iz) +  &
        w7 * fdist_wksp(i, j1, k1, iz) + w8 * fdist_wksp(i1, j1, k1, iz)

if (.false.) then
!if (DEBUG) then
  if (fdist /= 0._rp) then
    write (101, *) 'i, j, k =', i, j, k
    write (101, *) 'x =', x
    write (101, *) 'x1 =', x1
    write (101, *) 'x2 =', x2
    write (101, *) 'x3 =', x3
    write (101, *) 'fdist =', fdist
    write (101, *) 'fdist_wksp(i , j , k , 1) =', fdist_wksp(i , j , k , 1)
    write (101, *) 'fdist_wksp(i1, j , k , 1) =', fdist_wksp(i1, j , k , 1)
    write (101, *) 'fdist_wksp(i , j1, k , 1) =', fdist_wksp(i , j1, k , 1)
    write (101, *) 'fdist_wksp(i1, j1, k , 1) =', fdist_wksp(i1, j1, k , 1)
    write (101, *) 'fdist_wksp(i , j , k1, 1) =', fdist_wksp(i , j , k1, 1)
    write (101, *) 'fdist_wksp(i1, j , k1, 1) =', fdist_wksp(i1, j , k1, 1)
    write (101, *) 'fdist_wksp(i , j1, k1, 1) =', fdist_wksp(i , j1, k1, 1)
    write (101, *) 'fdist_wksp(i1, j1, k1, 1) =', fdist_wksp(i1, j1, k1, 1)
    write (101, *) ' '
  end if
end if

end subroutine interp_fdist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--cumulative sum over this branch and its sub-branches
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_ftot_ta ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_ftot_ta'

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name // trim (chcoord))
$endif

do i = 1, n_tree
  call fill_ftot (tree_array(i) % trunk)
end do

$if ($VERBOSE)
 call exit_sub (sub_name // trim (chcoord))
$endif

end subroutine fill_ftot_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this calculates ftot at all generations
!--may want to add optional argument to pick a min. generation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_ftot (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_ftot'

integer :: i

!---------------------------------------------------------------------

!--must descend to highest branch generation first, this allows
!  easier recursive relationship if we want to operate on
!  arbitrary generation later
if (br % gen < tree_array(br % itree) % n_gen) then

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_ftot (br % sub_branch(i))
  end do
  
end if

call fill_ftot_br (br)

end subroutine fill_ftot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes ftot had already been calculated at sub-branches of br
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_ftot_br (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_ftot_br'

integer :: i

!---------------------------------------------------------------------

if (br % resolved) then

  !--contribution from this branch
  br % resftot = br % resf
  br % unresftot = 0._rp

  !--this is not an error for a-priori test
  !if (.not. associated (br % sub_branch)) then
  !  call error (sub_name, 'expecting associated sub-branches')
  !end if

  if (associated (br % sub_branch)) then
    do i = 1, br % n_sub_branch
      br % resftot = br % resftot + br % sub_branch(i) % resftot
      br % unresftot = br % unresftot + br % sub_branch(i) % unresftot
    end do
  end if

else

  br % resftot = 0._rp
  br % unresftot = br % unresf
  
end if

br % ftot = br % resftot + br % unresftot

end subroutine fill_ftot_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_unresf_ta ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_unresf_ta'

integer :: i
!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

do i = 1, n_tree
  call fill_unresf (tree_array(i) % trunk)
end do

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine fill_unresf_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--br % unresf for resolved branches is not touched (and since it is 
!  initialized to zero, it should still be zero)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_unresf (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_unresf'

integer :: i

!---------------------------------------------------------------------

if (.not. br % resolved) then

  call fill_unresf_br (br)  !--maybe put the body of this inline

else

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_unresf (br % sub_branch(i))
  end do

end if

end subroutine fill_unresf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fdist_init ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fdist_init'

!---------------------------------------------------------------------

call eval_wksp_size ()
  
if (allocated (fdist_wksp)) then
  call error (sub_name, 'fdist_wksp already allocated')
end if
  
allocate (fdist_wksp(wksp_size(1), wksp_size(2), wksp_size(3), nzone))
allocate ( fdist_mask(wksp_size(1), wksp_size(2), wksp_size(3)) )

call mask_init (fdist_mask)

$if ($VERBOSE)
 call mesg (sub_name, 'fdist_wksp, fdist_mask initialized')
$endif
  
end subroutine fdist_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes that ftot has been updated already (for fdist_num_den)
!--if wksp_size is large, may be helpful to write out some of the
!  array expressions as loops to keep stack requirements down
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_fdist ()
implicit none

character (*), parameter :: sub_name = mod_name // '.update_fdist'

logical, parameter :: static_fdist = .true.
logical, parameter :: filt_static_fdist = .true.
                      !--only filters if static_fdist is true as well

integer :: iz

logical, save :: do_fdist_init = .true.

real (rp) :: total

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sub_name )
$endif

if ( do_fdist_init ) then

    if ( use_local_vel .and. (.not. static_fdist) )  &
        call error ( sub_name, 'use_local_vel=T requires static_fdist=T' )

    call fdist_init ()  !--allocate arrays, mask, etc.
 
    if ( static_fdist ) then
  
        if ( filt_static_fdist ) then
            call mask_filter ( fdist_mask, fdist_wksp(:, :, :, 1) )
        else
            fdist_wksp(:, :, :, 1) = fdist_mask(:, :, :)
        end if

        do iz = 2, nzone
            fdist_wksp(:, :, :, iz) = fdist_wksp(:, :, :, 1)
        end do

        if ( .not. use_local_vel ) call normalize_fdist_wksp ()
            !--leave fdist as mask: [0,1] when use_local_vel
    
        !--at this point we do not need to update fdist_mask anymore
    
    end if

    do_fdist_init = .false.

end if

if (.not. static_fdist) then

  if ((init_from_file .and. (jt == 1)) .and. (nfdistupdate /= 1)) then
    !--when nfdistupdate = 1, this will give different results
    !  from a continuous run since fx, fy, fz will not be the same

    call read_fdist_file ()
  
  else  !--right now we are not prescribing any fdist, just let go by itself

    call dynamic_fdist ()  !--result in fdist_wksp

  end if
  
  call normalize_fdist_wksp ()  !--normalize, even if read from file
    
end if

$if ($DEBUG)
if (DEBUG) then
  !--check normalization
  do iz = 1, nzone
    total = sum (fdist_wksp(:, :, :, iz))
    call mesg (sub_name, 'total =', total)
  end do
  call mesg (sub_name, '1/dV =', 1._rp / (dx * dy * dz))
end if
$endif

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine update_fdist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--does not normalize the resulting fdist_wksp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dynamic_fdist ()
implicit none

character (*), parameter :: sub_name = mod_name // '.dynamic_fdist'

real (rp), parameter :: fdist_clip_min = 0._rp
real (rp), parameter :: fdist_clip_max = huge ( 1._rp )
                        !--experiment to get good value?
                        !  1000 * dx * dy * dz worked?

character (64) :: fmt

integer :: i, j, k
integer :: iz

real (rp) :: den(nzone)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

call fdist_num_den_ta (wksp_size, fdist_wksp, den)
     !--this does not calculate fdist, just num and den
     !--fdist_wksp will hold numerator (save storage)

if (do_tavg_fdist) then
  call tavg_fdist (init_from_file, wksp_size, fdist_wksp, den)
       !--this call: num is stored in fdist_wksp
end if
  
!--perform division: is this really needed in light of the normalization?
do iz = 1, nzone

  if (den(iz) > 0._rp) then

    do k = 1, wksp_size(3)
      do j = 1, wksp_size(2)
        do i = 1, wksp_size(1)
          
          fdist_wksp(i, j, k, iz) = fdist_wksp(i, j, k, iz) / den(iz)
            
        end do
      end do
    end do
      
  else
    
    write (fmt, '(a,i0,a)') '(a,i0,a,a,', nzone, '(es13.6,1x))'
    write (msg, fmt) 'encountered den <= 0. for zone = ',  &
                                 iz, n_l, 'den = ', den

    call warn (sub_name, msg)

    fdist_wksp(:, :, :, iz) = 1._rp  !--0 here will screw-up normalization

  end if

end do

!--clip fdist_wksp here, (recall we assume alpha > 0)
where (fdist_wksp < fdist_clip_min) fdist_wksp = fdist_clip_min

where (fdist_wksp > fdist_clip_max) fdist_wksp = fdist_clip_max

!--special treatment when CD is clipped (experimental)
!do iz = 1, nzone
!  if (fcoeff_clipped(iz)) fdist_wksp(:, :, :, iz) = 1._rp
!end do

!--multiply by mask
do iz = 1, nzone
  fdist_wksp(:, :, :, iz) = fdist_wksp(:, :, :, iz) * fdist_mask
end do

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine dynamic_fdist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_fdist_file ()
implicit none

character (*), parameter :: sub_name = mod_name // '.read_fdist_file'

integer, parameter :: lun = 1

logical :: exst

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

inquire (file=file_fdist_last, exist=exst)
if (.not. exst) call error (sub_name, 'file' //                       &
                            trim (file_fdist_last) // 'does not exist')

!--make sure format matches that in trees_ls_finalize
open (lun, file=file_fdist_last, action='read', form='unformatted',  &
      position='rewind')
read (lun) fdist_wksp
close (lun)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine read_fdist_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normalize_fdist_wksp ()
implicit none

character (*), parameter :: sub_name = mod_name // '.normalize_fdist_wksp'

integer :: iz

real (rp) :: total

!---------------------------------------------------------------------

!--perform normalization: integral(alpha, dV) = 1, so sum (alpha) = 1/dV
do iz = 1, nzone
  total = sum (fdist_wksp(:, :, :, iz))
  fdist_wksp(:, :, :, iz) = fdist_wksp(:, :, :, iz) /  &
                            (total * dx * dy * dz)
end do

end subroutine normalize_fdist_wksp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--not going to bother with branch cross section here, since that
!  should be taken care of by sampling the force
!--purpose is mainly to avoid cross-talk between overlapping bboxes
!  from different branches
!--set k = 1 part of mask to zero so our unresf cannot interfere with
!  the resolved force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mask_init (mask)
implicit none

real (rp), intent (out) :: mask(:, :, :)

character (*), parameter :: sub_name = mod_name // '.mask_init'
character (*), parameter :: fout = 'mask.dat'

integer :: addlgen
integer :: n_gen
integer :: i, j, k

logical :: opn

real (rp) :: l0, d0
real (rp) :: l, d
real (rp) :: r
real (rp) :: x0(nd)
real (rp) :: dir(nd)
real (rp) :: x_hat(nd), y_hat(nd), z_hat(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

l0 = tree_array(n_tree) % l
d0 = tree_array(n_tree) % d

!--could make this more detailed: add_base, etc.
if (add_cap) l0 = l0 + 0.5_rp * d0

r = tree_array(n_tree) % ratio
n_gen = tree_array(n_tree) % n_gen  !--mask is size of bbox at n_gen

!--branch dimensions at n_gen
l = l0 * (r**n_gen)
d = d0 * (r**n_gen)

addlgen = 0  !--indicates how many generations beyond n_gen we are at

!--these are relative to fdist_wksp coordinates
x0 = 0._rp
dir = (/ 0._rp, 0._rp, 1._rp /)
x_hat = (/ 1._rp, 0._rp, 0._rp /)
y_hat = (/ 0._rp, 1._rp, 0._rp /)
z_hat = (/ 0._rp, 0._rp, 1._rp /)

mask = 0._rp

!--this will flip mask to 1 in certain locations
call calc_mask (addlgen, l, d, x0, dir, x_hat, y_hat, z_hat, mask)

!--zero out k = 1 to disallow interference with resolved force
mask(:, :, 1) = 0._rp

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit already open')

open (1, file=fout, action='write', position='rewind')

write (1, '(a)') 'variables = "i", "j", "k", "mask"'
write (1, '(3(a,i0))') 'zone, f=point, i=', size (mask, 1),  &
                       ', j=', size (mask, 2),               &
                       ', k=', size (mask, 3)

do k = 1, size (mask, 3)
  do j = 1, size (mask, 2)
    do i = 1, size (mask, 1)

      write (1, '(3(i0,1x),es12.5)') i, j, k, mask(i, j, k)

    end do
  end do
end do

close (1)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine mask_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--dat holds the data and boundary padding, which can be zero, or wrapped
!  copies of the data, so that periodic or zero BC can be used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filt_gauss_1d ( npd, nhw, ndat, dat, fdat )
implicit none

integer, intent (in) :: npd
    !--# boundary pts for padding, on each side of data
integer, intent (in) :: nhw
    !--the half-width of the filter in terms of number of grid-cells
    !  i.e. Delta/2h
integer, intent (in) :: ndat  !--number of actual data points

real (rp), intent (in) :: dat(ndat+2*npd)
real (rp), intent (out) :: fdat(ndat)

integer :: i, im, ip
integer :: j

real (rp) :: pow
real (rp) :: w
real (rp) :: norm

!---------------------------------------------------------------------

pow = -6.0_rp / (4 * nhw**2) 
w = exp ( pow )

!--for Gaussian filter, npd sets the integration range, for npd small
!  this means that we calculate only a rough approx to a Gaussian filter
!  For this reason, we also keep track of the normalization factor, since
!  it may differ significantly from the infinite domain value.
norm = 1.0_rp + w**(npd**2)
do j = 1, npd-1
    norm = norm + 2.0_rp * w **(j**2)
end do

!--trapezoidal rule, for the Gaussian-weighted integration
do i = npd + 1, npd + ndat

    im = i - npd
    ip = i + npd
    
    fdat(i-npd) = 0.5_rp * ( dat(im) + dat(ip) ) * w**(npd**2)

    fdat(i-npd) = fdat(i-npd) + dat(i)

    do j = 1, npd-1
        fdat(i-npd) = fdat(i-npd) + ( dat(i+j) + dat(i-j) ) * w**(j**2)
    end do

end do

fdat = fdat / norm

end subroutine filt_gauss_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--dat holds the data and boundary padding, which can be zero, or wrapped
!  copies of the data, so that periodic or zero BC can be used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filt_box_1d ( npd, ndat, dat, fdat )
implicit none

integer, intent (in) :: npd
    !--# boundary pts for padding, on each side of data
    !--for the box filter, this also the half-width of the filter
    !  in terms of number of grid-cells
integer, intent (in) :: ndat  !--number of actual data points

real (rp), intent (in) :: dat(ndat+2*npd)
real (rp), intent (out) :: fdat(ndat)

integer :: i, im, ip

!---------------------------------------------------------------------

!--trapezoidal rule
do i = npd + 1, npd + ndat
    im = i - npd
    ip = i + npd
    fdat(i-npd) = 0.5_rp * ( dat(im) + dat(ip) ) + sum ( dat( im+1:ip-1 ) )
end do

fdat = fdat / (2 * npd)

end subroutine filt_box_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--apply 3-d physical-space box filter with length scale delta = 2*h
!  (h = grid spacing)
!--mask, filtmask must have same shape & size
!--enforce filtmask = 0 at b.box k = 1 to avoid interference with the
!  resolved forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mask_filter (mask, filtmask)
implicit none

real (rp), intent (in) :: mask(:, :, :)
real (rp), intent (out) :: filtmask(:, :, :)

character (*), parameter :: sub_name = mod_name // '.mask_filter'
character (*), parameter :: filter = 'gauss'
    !--'box', 'gauss'
    
integer, parameter :: ido2h = 2
    !--idelta_over_2h: Delta/(2h), h = dx, dy, or dz
    !--sets the filter width: 1 -> 2h, 2-> 4h, etc.
!integer, parameter :: npad = 2*ido2h
    !--for Gaussian, want this to be at least 2 * ido2h

logical, parameter :: zero_k1 = .true.

!real (rp), parameter :: wgt = 6._rp

integer :: i, j, k
integer :: m1, m2, m3
!integer :: im, ip, jm, jp, km, kp, imm, ipp, jmm, jpp, kmm, kpp
integer :: m
integer :: npad

real (rp) :: mask_m, mask_p
real (rp), allocatable :: tmp(:), ftmp(:)
    
!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

m1 = size (mask, 1)
m2 = size (mask, 2)
m3 = size (mask, 3)

if ( (m1 /= size (filtmask, 1)) .or.  &
     (m2 /= size (filtmask, 2)) .or.  &
     (m3 /= size (filtmask, 3)) ) then
  call error (sub_name, 'mask, filtmask size mis-match')
end if

m = max ( size (mask, 1), size (mask, 2), size (mask, 3) )

select case ( filter )
case ( 'box' )
    npad = ido2h
case ( 'gauss' )
    npad = 2*ido2h  !--very limited integration region
case default
    call error ( sub_name, 'invalid filter =' // filter )
end select

allocate ( tmp( m + 2 * npad ) )  !--rooom for data plus padding
allocate ( ftmp(m) )  !--just room for data

!--speed is not of primary concern here (only done once)
!--assume everthing outside the index range is zero

do k = 1, m3
  do j = 1, m2

    tmp = 0.0_rp
    tmp( npad+1 : npad+m1 ) = mask( 1:m1, j, k )

    select case ( filter )
    case ( 'box' )
        call filt_box_1d ( npad, m1, tmp, ftmp )
    case ( 'gauss' )
        call filt_gauss_1d ( npad, ido2h, m1, tmp, ftmp )
    end select
    
    filtmask( 1:m1, j, k ) = ftmp( 1:m1 )

    !do i = 1, m1
   
    !  im = max ( i-ido2h, 1 )
    !  imm = max ( im+1, 1 )
    !  
    !  ip = min ( i+ido2h, m1 )
    !  ipp = min ( ip-1, m1 )
    ! 
    !  if ( im < 1 ) then
    !    mask_m = 0.0_rp
    !  else
    !    mask_m = mask( im, j, k )
    !  end if
    !  
    !  if ( ip > m1 ) then
    !    mask_p = 0.0_rp
    !  else
    !    mask_p = mask( ip, j, k )
    !  end if
    !  
    !  filtmask( i, j, k ) =  ( 0.5_rp * mask_m +             &
    !                           0.5_rp * mask_p +             &
    !                           sum ( mask( imm:ipp, j, k ) ) ) 

    !end do

  end do
end do

!--note order of loops
do k = 1, m3
  do i = 1, m1

    tmp = 0.0_rp
    tmp( npad+1 : npad+m2 ) = filtmask( i, 1:m2, k )
    
    select case ( filter )
    case ( 'box' )
        call filt_box_1d ( npad, m2, tmp, ftmp )
    case ( 'gauss' )
        call filt_gauss_1d ( npad, ido2h, m2, tmp, ftmp )
    end select

    filtmask( i, 1:m2, k ) = ftmp( 1:m2 )
    
    !tmp(1:m2) = filtmask( i, 1:m2, k )
    !
    !do j = 1, m2
    !
    !  jm = max ( j-ido2h, 1 )
    !  jmm = max ( jm+1, 1 )
    !
    !  jp = min ( j+ido2h, m2 )
    !  jpp = min ( jp-1, m2 )

    !  if ( jm < 1 ) then
    !    mask_m = 0.0_rp
    !  else
    !    mask_m = tmp(jm)
    !  end if

    !  if ( jp > m2 ) then
    !    mask_p = 0.0_rp
    !  else
    !    mask_p = tmp(jp)
    !  end if

    !  filtmask( i, j, k ) = ( 0.5_rp * mask_m + 0.5_rp * mask_p +  &
    !                          sum ( tmp(jmm:jpp) ) )

    !end do

  end do
end do

!--note order of loops
do j = 1, m2
  do i = 1, m1
  
    tmp = 0.0_rp
    tmp( npad+1 : npad+m3 ) = filtmask( i, j, 1:m3 )
    
    select case ( filter )
    case ( 'box' )
        call filt_box_1d ( npad, m3, tmp, ftmp )
    case ( 'gauss' )
        call filt_gauss_1d ( npad, ido2h, m3, tmp, ftmp )
    end select

    filtmask( i, j, 1:m3 ) = ftmp( 1:m3 )

!    tmp(1:m3) = filtmask( i, j, 1:m3 )
!    
!    do k = 1, m3
!    
!      km = max ( k-ido2h, 1 )
!      kmm = max ( km+1, 1 )
!    
!      kp = min ( k+ido2h, m3 )
!      kpp = min ( kp-1, m3 )
!
!      if ( km < 1 ) then
!        mask_m = 0.0_rp
!      else
!        mask_m = tmp(km)
!      end if
!
!      if ( kp > m3 ) then
!        mask_p = 0.0_rp
!      else
!        mask_p = tmp(kp)
!      end if
!
!      filtmask( i, j, k ) = ( 0.5_rp * mask_m + 0.5_rp * mask_p +  &
!                              sum ( tmp(kmm:kpp) ) )
!
!    end do

  end do
end do

!filtmask = filtmask / ( 2 * ido2h )**3

!do k = 1, m3
!  do j = 1, m2
!    do i = 1, m1
!
!      if (i <= 1) then
!        mask_im1 = 0._rp
!      else
!        mask_im1 = mask(i - 1, j, k)
!      end if
!
!      if (i >= m1) then
!        mask_ip1 = 0._rp
!      else
!        mask_ip1 = mask(i + 1, j, k)
!      end if
!     
!      if (j <= 1) then
!        mask_jm1 = 0._rp
!      else
!        mask_jm1 = mask(i, j - 1, k)
!      end if
!
!      if (j >= m2) then
!        mask_jp1 = 0._rp
!      else
!        mask_jp1 = mask(i, j + 1, k)
!      end if
!
!      if (k <= 1) then
!        mask_km1 = 0._rp
!      else
!        mask_km1 = mask(i, j, k - 1)
!      end if
!
!      if (k >= m3) then
!        mask_kp1 = 0._rp
!      else
!        mask_kp1 = mask(i, j, k + 1)
!      end if
!      
!      filtmask(i, j, k) = ( mask_im1 + mask_ip1 +               &
!                            mask_jm1 + mask_jp1 +               &
!                            mask_km1 + mask_kp1 +               &
!                            wgt * mask(i, j, k) ) / (2._rp * wgt)
!
!    end do
!  end do
!end do

!--tmp, ftmp should be auto-deallocated
deallocate ( ftmp )
deallocate ( tmp )

if (zero_k1) filtmask(:, :, 1) = 0._rp  !--avoid interference w/ resolved f

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine mask_filter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--beware, the integration accuracy is low
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function box_kernel_3d ( x, delta )
implicit none

real (rp) :: box_kernel_3d

real (rp), intent (in) :: x(nd)
real (rp), intent (in) :: delta

!---------------------------------------------------------------------

if ( any ( abs (x) > delta / 2.0_rp ) ) then
    box_kernel_3d = 0.0_rp
else
    box_kernel_3d = 1.0_rp / delta**3
end if

end function box_kernel_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gaussian_kernel_3d ( x, delta )
implicit none

real (rp) :: gaussian_kernel_3d

real (rp), intent (in) :: x(nd)
real (rp), intent (in) :: delta

real (rp) :: pi

!---------------------------------------------------------------------

pi = acos( -1.0_rp )

gaussian_kernel_3d = ( 6.0_rp / pi / delta**2 )**(1.5_rp) *  &
                     exp ( -6.0_rp * mag (x)**2 / delta**2 )

end function gaussian_kernel_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes mask has been set to 0 before top-most call
!--only flips mask to 1 at certain locations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine calc_mask (addlgen, l, d, x0, dir, x_hat, y_hat,  &
                                z_hat, mask)
implicit none

integer, intent (in) :: addlgen
real (rp), intent (in) :: l, d
real (rp), intent (in) :: x0(nd)
real (rp), intent (in) :: dir(nd)
real (rp), intent (in) :: x_hat(nd), y_hat(nd), z_hat(nd)
real (rp), intent (inout) :: mask(:, :, :)

character (*), parameter :: sub_name = mod_name // '.calc_mask'

integer, parameter :: addlgen_max = 5  !--depth of recursion

logical, parameter :: use_loose_mask = .false.  !--expands mask slightly
    !--may want this false when filtering mask

real (rp), parameter :: thresh = 100._rp * epsilon (0._rp)

integer :: i, j, k
integer :: io, jo, ko
integer :: nsub

real (rp) :: d_para
real (rp) :: dir_sub(nd)
real (rp) :: eps
real (rp) :: r
real (rp) :: twist
real (rp) :: xtip(nd)
real (rp) :: x(nd)
real (rp) :: x_perp(nd)
real (rp) :: x0_sub(nd)
real (rp) :: x_hat_sub(nd), y_hat_sub(nd), z_hat_sub(nd)
real (rp) :: x_tmp(nd), y_tmp(nd)

!---------------------------------------------------------------------

eps = 10.0_rp * epsilon ( 1._rp )  !--this is the fudge factor in our mask

if ( use_loose_mask ) eps = eps + dx

r = tree_array(n_tree) % ratio  !--assume self similar

xtip = x0 + l * dir  !--mask this branch

!--these define the origin of the local coordinate (index space)
io = size (mask, dim = 1) / 2 + 1
jo = size (mask, dim = 2) / 2 + 1
ko = 1

do k = 1, size (mask, 3)
  do j = 1, size (mask, 2)
    do i = 1, size (mask, 1)

      x = (/ i - io, j - jo, k - ko /) * (grid % dx) - x0

      d_para = dot_product (x, dir)

      if ( (d_para  >= 0._rp) .and. (d_para <= l + eps) ) then

        x_perp = x - d_para * dir
        
        !--square cross-section test
        if (maxval (abs (x_perp)) <= 0.5_rp * d + eps) mask(i, j, k) = 1._rp

      end if

    end do
  end do
end do

!--recursion part
if (addlgen < addlgen_max) then

  nsub = tree_array(n_tree) % n_sub_branch

  do i = 1, nsub

    !--make sure this agrees with add_sub_branches
    dir_sub = tree_array(n_tree) % rel_dir(1, i) * x_hat +  &
              tree_array(n_tree) % rel_dir(2, i) * y_hat +  &
              tree_array(n_tree) % rel_dir(3, i) * z_hat

    z_hat_sub = dir_sub
    x_hat_sub = cross_product (z_hat_sub, dir)
    y_hat_sub = cross_product (z_hat_sub, x_hat_sub)

    !--if this coordinate system is degenerate (x/y_hat_sub = 0)
    !  then just use parents local coordinate system
    if (maxval (abs (y_hat_sub)) < thresh) then

      $if ($DEBUG)
      if (DEBUG) then
        call mesg (sub_name, "sub's coords degenerate, using parent's")
      end if
      $endif
      x_hat_sub = x_hat
      y_hat_sub = y_hat
      z_hat_sub = z_hat

    else  !--normalize

      x_hat_sub = x_hat_sub / mag (x_hat_sub)
      y_hat_sub = y_hat_sub / mag (y_hat_sub)
      z_hat_sub = z_hat_sub / mag (z_hat_sub)
      
    end if

    twist = tree_array(n_tree) % twist(i)

    !--apply twist
    x_tmp = cos (twist) * x_hat_sub + sin (twist) * y_hat_sub
    y_tmp = cos (twist) * y_hat_sub - sin (twist) * x_hat_sub
    x_hat_sub = x_tmp
    y_hat_sub = y_tmp

    !--may need to check to see if other options apply
    if ((add_cap) .and. (sub_branches_outside)) then

      x0_sub = x0 + (tree_array(n_tree) % root_height(i) *  &
                     (l - 0.5_rp * d) * dir)                &
                  + (0.5_rp * d) * dir_sub

    else

      call error (sub_name, 'current x0_sub only valid for' // n_l //  &
                            'add_cap .and. sub_branches_outside')

    end if

    call calc_mask (addlgen + 1, r * l, r * d, x0_sub, dir_sub,  &
                    x_hat_sub, y_hat_sub, z_hat_sub, mask)

  end do

end if

end subroutine calc_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this does not calculate an fdist value, only the numerator and
!  in the denominator in the expression for fdist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fdist_num_den_ta (nwksp, num, den)
$if ($MPI)
  use param, only : up, down, status
  use immersedbc, only : fz
$endif
implicit none

integer, intent (in) :: nwksp(nd)
real (rp), intent (out) :: num(nwksp(1), nwksp(2), nwksp(3), nzone)
real (rp), intent (out) :: den(nzone)

integer :: i

$if ($MPI)
  integer, parameter :: tag = 444

  !--if num_sum causes stack problems, make allocatable
  real (rp) :: num_sum(nwksp(1), nwksp(2), nwksp(3), nzone)
  real (rp) :: den_sum(nzone)
$endif

!---------------------------------------------------------------------

$if ($MPI)
  !--sync fz(nz) <-> fz(1'): needed in fdist_num_den
  !--may remove this if change fdist_num_den to not bother
  !  interpolating fz to u-nodes (since brindex stuff is also like that)
  call mpi_sendrecv (fz(1, 1, 1), ld*ny, MPI_RPREC, down, tag,  &
                     fz(1, 1, nz), ld*ny, MPI_RPREC, up, tag,   &
                     comm, status, ierr)
$endif

num = 0._rp
den = 0._rp

do i = 1, n_tree
  call fdist_num_den (tree_array(i) % trunk, nwksp, num, den)
end do

$if ($MPI)
  call mpi_allreduce (num, num_sum, product (nwksp) * nzone, MPI_RPREC,  &
                      MPI_SUM, comm, ierr)
  call mpi_allreduce (den, den_sum, nzone, MPI_RPREC, MPI_SUM, comm, ierr)

  num = num_sum
  den = den_sum
$endif

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (num(:, :, :, 1), 'fdist_num_den_ta.num')
  call DEBUG_write (den, 'fdist_num_den_ta.den')
end if
$endif

end subroutine fdist_num_den_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fdist_num_den (br, nwksp, num, den) 
implicit none

type (branch_type), intent (in) :: br
integer, intent (in) :: nwksp(nd)
real (rp), intent (inout) :: num(nwksp(1), nwksp(2), nwksp(3), nzone)
real (rp), intent (inout) :: den(nzone)

character (*), parameter :: sub_name = mod_name // '.fdist_num_den'

integer :: i

!---------------------------------------------------------------------

if (br % gen == tree_array(br % itree) % max_res_gen) then

  call fdist_num_den_br (br, nwksp, num, den)

else if (br % gen < tree_array(br % itree) % max_res_gen) then

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if
  
  do i = 1, br % n_sub_branch
    call fdist_num_den (br % sub_branch(i), nwksp, num, den)
  end do

end if

end subroutine fdist_num_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes br % ftot has already been calculated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fdist_num_den_br (br, nwksp, num, den)
use immersedbc, only : fx, fy, fz
implicit none

type (branch_type), intent (in) :: br
integer, intent (in) :: nwksp(nd)
real (rp), intent (inout) :: num(nwksp(1), nwksp(2), nwksp(3), nzone)
real (rp), intent (inout) :: den(nzone)

character (*), parameter :: sub_name = mod_name // '.fdist_num_den_br'


integer :: iz
integer :: ipt
integer :: i, j, k
integer :: ip, jp, kp

real (rp) :: r
real (rp) :: f(nd)
real (rp) :: x(nd)
real (rp) :: xp(nd)

!---------------------------------------------------------------------

iz = br % zone

r = tree_array(br % itree) % ratio

do ipt = 1, br % nbboxpt

  i = br % bboxpt(1, ipt)
  j = br % bboxpt(2, ipt)
  k = br % bboxpt(3, ipt)

  x(1) = pt_of_grid(i, 1, 1)  !--u-nodes
  x(2) = pt_of_grid(j, 2, 1)
  x(3) = pt_of_grid(k, 3, 1)

  x = x - br % x0

  !--convert to branch local coords
  xp(1) = dot_product (x, br % x_hat)
  xp(2) = dot_product (x, br % y_hat)
  xp(3) = dot_product (x, br % z_hat)
  
  !--rescale (so we are at n_gen scale)
  !--do this since this is the scale at which we are going to apply
  !  the forcing
  xp = r * xp

  !--calculate indices for wksp: origin of x,y coords are in middle of wksp
  ip = floor (xp(1) / dx) + (wksp_size(1) + 1) / 2
  jp = floor (xp(2) / dx) + (wksp_size(2) + 1) / 2
  kp = floor (xp(3) / dx) + 1

  call check_bounds ()

  !f = (/ fx(i, j, k), fy(i, j, k), fz(i, j, k) /)  !--debug

  if (use_unif_fdist_br) then
    !--special treatment of resolved branches to remove large
    !  spikes in force: smear out force over the branch
  
    !--make sure that br % resf is up to date
    if (brindex(i, j, k) == br % ident) then
      f = (br % resf) / ((br % nrespt) * dx * dy * dz)
    else
      !--MPI: k+1 here requires a resync of fz(k=nz)
      f = (/ fx(i, j, k), fy(i, j, k),                 &
             0.5_rp * (fz(i, j, k) + fz(i, j, k + 1)) /)
    end if
    
  else

    !--MPI: k+1 here requires a resync of fz(k=nz)
    f = (/ fx(i, j, k), fy(i, j, k),                 &
           0.5_rp * (fz(i, j, k) + fz(i, j, k + 1)) /)
           
  end if
  
  !--must have br % ftot already
  num(ip, jp, kp, iz) = num(ip, jp, kp, iz) +    &
                        dot_product (f, br % ftot)  

  den(iz) = den(iz) + dot_product (br % ftot, br % ftot)
 
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_bounds ()
  implicit none

  character (32) :: fmt

  !-------------------------------------------------------------------

  fmt = '(3a,3(i0,1x),2a,3(i0,1x))'
  if ( (ip < 1) .or. (ip > wksp_size(1)) ) then
    write (msg, fmt) 'ip is out of bounds', n_l,            &
                     '(ip,jp,kp)=', (/ ip, jp, kp /), n_l,  &
                     'wksp_size =', wksp_size
    call error (sub_name, msg)
  end if

  if ( (jp < 1) .or. (jp > wksp_size(2)) ) then
    write (msg, fmt) 'jp is out of bounds', n_l,            &
                     '(ip,jp,kp)=', (/ ip, jp, kp /), n_l,  &
                     'wksp_size =', wksp_size
    call error (sub_name, msg)
  end if

  if ( (kp < 1) .or. (kp > wksp_size(3)) ) then
    write (msg, fmt) 'kp is out of bounds', n_l,            &
                     '(ip,jp,kp)=', (/ ip, jp, kp /), n_l,  &
                     'wksp_size =', wksp_size
    call error (sub_name, msg)
  end if

  end subroutine check_bounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine fdist_num_den_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this is needed for unresf stuff
!--bounding box dimensions for bbox at n_gen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_wksp_size ()
implicit none

character (*), parameter :: sub_name = mod_name // '.eval_wksp_size'

type (branch_type), pointer :: brp => NULL()

real (rp) :: height, width
real (rp) :: margin

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

margin = 1.5_rp * dx  !--add extra space around the boundaries

!--assume dx = dy = dz here: since in local coords the bbox may be
!  aligned in any direction it is hard to do anything without this
!--this assumes all trees are identical
!--this assumes the bbox size of all the sub-branches are identical

brp => tree_array(n_tree) % trunk

do

  if (brp % gen == tree_array(n_tree) % n_gen) then
    width = brp % width_bbox + margin
    height = brp % height_bbox + margin
    exit
  else
    if (.not. associated (brp % sub_branch)) then
      call error (sub_name, 'expecting associated sub-branches')
    end if
    brp => brp % sub_branch(brp % n_sub_branch)
  end if
  
end do

brp => NULL ()  !--we are done with brp

wksp_size(1) = 2 * ceiling (width / 2._rp / dx) + 1  !--odd for nice centering
wksp_size(2) = 2 * ceiling (width / 2._rp / dx) + 1
wksp_size(3) = ceiling (height / dx)

call mesg (sub_name, 'wksp_size =', wksp_size)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine eval_wksp_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--no nfcoeff-dependence for tscale right now, may add later
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_tavg_timescale (n1, tscale)
implicit none

integer, intent (in) :: n1
real (rp), intent (out) :: tscale(n1, nzone)

character (*), parameter :: sub_name = mod_name // '.set_tavg_timescale'
!character (*), parameter :: method = 'L-V'
                            !--'L-V', 'N-D'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif
real (rp), parameter :: eps = 100._rp * epsilon (0._rp)
real (rp), parameter :: c1 = 5.0_rp
                        !--1/St = 1/0.13 = 7.69 square cyl.
                        !--1/St = 1/0.2 = 5.0 for circ. cyl.

integer :: i
integer :: gen

real (rp) :: vel_z(nzone)
real (rp) :: length
real (rp) :: r

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

gen = tree_array(n_tree) % max_res_gen
!gen = tree_array(n_tree) % n_gen

!--for now: height of bbox?
r = tree_array(n_tree) % ratio
!--use diam of branch at gen
length = (r**gen) * (tree_array(n_tree) % trunk % d)

!--use geometric mean of bbox dimensions at gen
!length = (r**gen) *  &
!         ( (tree_array(n_tree) % trunk % width_bbox)**2 *              &
!           (tree_array(n_tree) % trunk % height_bbox) )**(1._rp / 3._rp) 

!select case (method)
!case ('L-V')
  
  call mean_velscale_z (gen, vel_z)
  do i = 1, n1  !--not i dependence here yet
    tscale(i, :) = c1 * length / max (vel_z(:), eps)
  end do
    
!case ('N-D')
!
!  !--not sure whether to use c1 here
!  tscale = length**2 / max ( (num * den)**(1._rp/8._rp), eps )
!    
!  if (DEBUG) then
!    call mesg (sub_name, 'num =', num)
!    call mesg (sub_name, 'den =', den)
!    call mesg (sub_name, 'tscale =', tscale)
!  end if
!
!end select
$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine set_tavg_timescale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mean_velscale_z (gen, vel)
implicit none

integer, intent (in) :: gen

real (rp), intent (out) :: vel(nzone)

character (*), parameter :: sub_name = mod_name // '.mean_velscale'

integer :: i
integer :: nvel(nzone)

!---------------------------------------------------------------------

vel = 0._rp
nvel = 0

do i = 1, n_tree
  call mean_velscale_br (tree_array(i) % trunk, gen, nvel, vel)
end do

if (all (nvel /= 0)) then
  vel = vel / nvel
else
  call error (sub_name, 'cannot determine mean vel. scale: nvel is 0')
end if

end subroutine mean_velscale_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine mean_velscale_br (br, gen, nvel, vel)
implicit none

type (branch_type), intent (in) :: br
integer, intent (in) :: gen
integer, intent (inout) :: nvel(nzone)
real (rp), intent (inout) :: vel(nzone)

character (*), parameter :: sub_name = mod_name // '.mean_velscale_br'

integer :: i

!---------------------------------------------------------------------

if (br % gen == gen) then

  vel(br % zone) = vel(br % zone) + mag (br % velscale)
  nvel(br % zone) = nvel(br % zone) + 1
  
else if (br % gen < gen) then

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting sub_branch to be associated')
  end if

  do i = 1, br % n_sub_branch
    call mean_velscale_br (br % sub_branch(i), gen, nvel, vel)
  end do
  
end if

end subroutine mean_velscale_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_max_res_gen ()
implicit none

character (*), parameter :: sub_name = mod_name // '.set_max_res_gen'

integer :: i

!---------------------------------------------------------------------

max_res_gen = tree_array (n_tree) % max_res_gen  !--sets module variable

!--sanity check: max_res_gen for all trees must be the same
do i = 1, n_tree - 1
  if (tree_array(i) % max_res_gen /= max_res_gen) then
    call error (sub_name, 'expecting all trees to have same max_res_gen')
  end if
end do
    
end subroutine set_max_res_gen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_fcoeff (update)
implicit none

logical, intent (out) :: update

character (*), parameter :: sub_name = mod_name // '.init_fcoeff'

integer, parameter :: jt_min = 20
integer, parameter :: lun = 1

real (rp), parameter :: fcoeff_initial = 1._rp

logical :: exst

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

if ((.not. init_from_file) .and. (jt < jt_min)) then

  fcoeff = fcoeff_initial
  update = .false.

else if (init_from_file .and. (jt == 1)) then  !--read fcoeff from file

  call mesg (sub_name, 'reading fcoeff from file')
  
  inquire (file=file_fcoeff_last, exist=exst)
    
  if (.not. exst) call error (sub_name, 'file ' //                        &
                              trim (file_fcoeff_last) // ' does not exist')
                              
  open (lun, file=file_fcoeff_last, action='read', form='unformatted',  &
        position='rewind')
  read (lun) fcoeff
  close (lun)

  if (nCDupdate == 1) then

    if (fmodel /= 'd_germano') call read_ftot(file_ftot_last)

    !--fcoeff update still required, old fcoeff is just for ftot calc.

    update = .true.
    
  else
  
    update = .false.
    
  end if

else

  update = .true.
  
end if

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine init_fcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open
!--this includes time in first column
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_fcoeff_line (lun)
use param, only : dt  !--also see top of module
implicit none

integer, intent (in) :: lun

character (64) :: fmt

!---------------------------------------------------------------------

!--get precision decimal places out
write (fmt, '(3(a,i0),a)') '(', size (fcoeff)+1, '(es',   &
                           precision (1._rp) + 7, '.',    &
                           precision (1._rp) - 1, ',1x))'

write (lun, fmt) dt * (jt_total - 1), fcoeff

end subroutine write_fcoeff_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_fcoeff (file)
implicit none

character (*), intent (in) :: file

character (*), parameter :: sub_name = mod_name // '.write_fcoeff'

integer, parameter :: lun = 1

logical :: opn, exst

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

inquire (unit=lun, exist=exst, opened=opn)
if (opn) call error (sub_name, 'unit is already open, unit=', lun)
if (.not. exst) call error (sub_name, 'unit does not exist, unit=', lun)

!--open and close every time to flush
open (lun, file=file, action='write', position='append')

call write_fcoeff_line (lun)

close (lun)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine write_fcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_tscale (n, tscale)
use param, only : dt  !--also see top of module
implicit none

integer, intent (in) :: n
real (rp), intent (in) :: tscale(n)

character (*), parameter :: sub_name = mod_name // '.write_tscale'

integer, parameter :: lun = 1

character (64) :: fmt

logical :: opn

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

!--make sure this lun is not already open
inquire (lun, opened=opn)
if (opn) call error (sub_name, 'unit already open, lun=', lun)
  
open (lun, file='fcoeff_tscale.dat', action='write', position='append')

write (fmt, '(a,i0,a)') '(', size (tscale) + 1, '(es13.6,1x))'
write (lun, fmt) dt * (jt_total - 1), tscale

close (lun)

$if ($VERBOSE) 
 call exit_sub (sub_name)
$endif

end subroutine write_tscale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--resf at max_res_gen MUST have been updated before calling this
!--assumes max_res_gen is same for all trees
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_fcoeff ()
implicit none

character (*), parameter :: sub_name = mod_name // '.update_fcoeff'

logical :: update

real (rp) :: tscale(nfcoeff, nzone)
real (rp) :: num(nfcoeff, nzone), den(nfcoeff, nzone)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

!--determination of whether or not to update is somewhat complex
!  therefore we put that part, and reading from file in to one routine
call init_fcoeff (update)

if (update) then

  if (fmodel /= 'd_germano') then
  
    !--get "base" value of coefficients averaged over zone/generations
    !--do not actually apply coefficient yet until t-avg and clipping applied
    call fcoeff_ta (max_res_gen, num, den)

    !--apply time averaging to coefficients
    if (do_tavg_fcoeff) call tavg_fcoeff (init_from_file, num, den, tscale)
                        !--this resets fcoeff, num_fcoeff, den_fcoeff

    !--to avoid repeating division in above routines, can take out and
    !  put into a routine called from right here
    !  however, this may mean that above routines should be renamed
    !  since they will not then actually calculate the coefficients

  else  !--germano: does not fit same form as others
  
    call fcoeff_d_germano ()
    !--no time averaging available for Germano
    
  end if
  
  !--apply clipping to coefficients
  !  watch interface: 2d fcoeff -> 1d dummy
  call clip_fcoeff (size (fcoeff), fcoeff, fcoeff_clipped)

end if

!--output to file
if ((modulo (jt, nCDupdate) == 0) .or. update) then

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  
    call write_fcoeff (file_fcoeff_dat)

    if (do_tavg_fcoeff) call write_tscale (size (tscale), tscale)

  end if
  
end if

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine update_fcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tavg_fcoeff (file_init, num_in, den_in, tscale)
use param, only : dt  !--also see top of module
$if ($XLF)
  use ieee_arithmetic  !--for NAN checking
$endif
implicit none

logical, intent (in) :: file_init
real (rp), intent (in) :: num_in(nfcoeff, nzone), den_in(nfcoeff, nzone)
real (rp), intent (out) :: tscale(nfcoeff, nzone)

character (*), parameter :: sub_name = mod_name // '.tavg_fcoeff'

integer, parameter :: lun = 1

integer :: i, z

logical :: nan(nfcoeff, nzone)
logical, save :: init = .false.

!real (rp) :: tscale(nfcoeff, nzone)
!--using module variables num_fcoeff, den_fcoeff so can use trees_finalize_ls
!  perhaps there is a better way?
!real (rp), save :: num(nfcoeff, nzone), den(nfcoeff, nzone)
real (rp) :: wgt

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

if (.not. init) then

  !--init num, den here
  if (file_init) then
  
    open (lun, file='num_den_fcoeff.last.out', action='read',  &
          form='unformatted', position='rewind')
          
    read (lun) num_fcoeff, den_fcoeff
  
    close (lun)
    
  else

    num_fcoeff = 0._rp
    den_fcoeff = 0._rp
    
  end if
  
  init = .true.
  
end if

!--determine time scale for time averaging
call set_tavg_timescale (nfcoeff, tscale)  !--no fcoeff dependence right now

!--protect from NaN:
nan = .false.

!--may be better to write compiler-dependent isnan routine to contain mess
$if ($IFC || $IFORT)
  nan = (isnan (tscale) .or. isnan (num_fcoeff) .or. isnan (den_fcoeff))
$elsif ($XLF)
  if (ieee_support_datatype (1._rp)) then
    if (ieee_support_nan (1._rp)) then
        nan = (ieee_is_nan (tscale) .or. ieee_is_nan (num_fcoeff) .or.  &
               ieee_is_nan (den_fcoeff))
    end if
  end if
$endif
   
if (any (nan)) then
    
  call mesg (sub_name, 'encountered NAN:')
  call mesg (sub_name, 'tscale=', pack (tscale, mask=.true.))
  call mesg (sub_name, 'num=', pack (num_fcoeff, mask=.true.))
  call mesg (sub_name, 'den=', pack (den_fcoeff, mask=.true.))
      
end if

do z = 1, nzone
  do i = 1, nfcoeff
  
    if (.not. nan(i, z)) then
  
      wgt = (nCDupdate * dt) / tscale(i, z)
      !wgt = (nCDupdate * dt) / (tscale(i, z) + nCDupdate*dt)
    
      !--weighted time avg. as in Lagrangian dynamic model
      num_fcoeff(i, z) = wgt * num_in(i, z) + (1._rp - wgt) * num_fcoeff(i, z)
      den_fcoeff(i, z) = wgt * den_in(i, z) + (1._rp - wgt) * den_fcoeff(i, z)
      
    else  !--nan

      !--leave num, den, fcoeff set to their previous old values
      
      tscale(i, z) = huge (1._rp)  !--for purpose of writing to file only
      
    end if
    
  end do
end do

call fcoeff_divide (num_fcoeff, den_fcoeff)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine tavg_fcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--does not modify fdist_wksp (unless it is used as an argument, i.e.
!  storage for num_in)
!--do_file_init will probably always be init_from_file, but having it as
!  an argument means that we can move this into the fmodel module
!--for now there is only 1 fdist, num, den, not nfcoeff of them
!--perhaps could unify this with tavg_fcoeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tavg_fdist (do_file_init, nwksp, num_in, den_in)
use param, only : dt  !--also see top of module
implicit none

logical, intent (in) :: do_file_init
integer, intent (in) :: nwksp(3)
real (rp), intent (in) :: num_in(nwksp(1), nwksp(2), nwksp(3), nzone)
real (rp), intent (in) :: den_in(nzone)

character (*), parameter :: sub_name = mod_name // '.tavg_fdist'

integer, parameter :: lun = 1

integer :: i, j, k, z
integer, save :: nwksp_old(3)

logical, save :: do_init = .true.

real (rp) :: tscale(nzone)
!--now using module variable num_fdist, den_fdist so can use trees_finalize_ls
!  perhaps there is a better way to do it
!--num, is allocatable since it needs to be saved
!real (rp), allocatable, save :: num(:, :, :, :)
!real (rp), save :: den(nzone)
real (rp) :: wgt

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

if (do_init) then

  !--num_fdist is module variable
  allocate (num_fdist(nwksp(1), nwksp(2), nwksp(3), nzone))
  nwksp_old = nwksp

  !--init num, den here
  if (do_file_init) then

    open (lun, file=file_num_den_fdist_last, action='read',  &
          form='unformatted', position='rewind')
          
    read (lun) num_fdist, den_fdist
  
    close (lun)
    
  else

    num_fdist = 0._rp
    den_fdist = 0._rp
    
  end if
  
  do_init = .false.

else  !--check nwksp is consistent with nwksp_old

  if (any (nwksp /= nwksp_old)) call error (sub_name, 'nwksp size has changed')
  
end if

!--determine time scale for time averaging
!--the 1 here is first dimension of tscale: 1 X nzone as 2d-array
call set_tavg_timescale (1, tscale)

!--need to add nan-checking here
do z = 1, nzone

  wgt = (nCDupdate * dt) / tscale(z)
  !wgt = (nCDupdate * dt) / (tscale(z) + nCDupdate*dt)
    
  !--weighted time avg. as in Lagrangian dynamic model (sort of)
  do k = 1, nwksp(3)
    do j = 1, nwksp(2)
      do i = 1, nwksp(1)
        
        num_fdist(i, j, k, z) = wgt * num_in(i, j, k, z) +          &
                                (1._rp - wgt) * num_fdist(i, j, k, z)
      end do
    end do
  end do
    
  den_fdist(z) = wgt * den_in(z) + (1._rp - wgt) * den_fdist(z)

  !--moved the division into update_fdist: avoids repitition/inconsistencies
    
end do

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine tavg_fdist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--pay attention to whether ftot is updated before calling this
!!  (its not critical either way, but it uses ftot)
!!--placing time averaging in here, may want to put into coeff_d
!!  at a later time
!!--reports timescale used in the averaging: < 0 if no averaging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine fcoeff_d_iter (timescale)
!use param, only : dt  !--also see top of module
!$if ($XLF)
!  use ieee_arithmetic  !--for NAN checking
!$endif
!implicit none
!
!real (rp), intent (out), optional :: timescale(nzone)
!
!character (*), parameter :: sub_name = mod_name // '.fcoeff_d_iter'
!
!!logical, parameter :: nd_time_avg = .true.
!
!logical :: nan(nzone)
!logical, save :: init = .false.
!
!integer :: z
!
!!real (rp), save :: num(nzone) = 0._rp, den(nzone) = 0._rp
!real (rp) :: numold(nzone), denold(nzone)
!real (rp) :: CDold(nzone)
!real (rp) :: wgt(nzone)
!
!!---------------------------------------------------------------------
!
!if (VERBOSE) call enter_sub (sub_name)
!
!if (.not. init) then
!
!  !--init num, den here
!  if (init_from_file) then
!  
!    open (1, file='num_den_fcoeff.last.out', action='read',  &
!          form='unformatted', position='rewind')
!          
!    read (1) num, den
!  
!    close (1)
!    
!  else
!
!    num = 0._rp
!    den = 0._rp
!    
!  end if
!  
!  init = .true.
!  
!end if
!
!!--save old value
!!if (time_avg) then
!!  if (nd_time_avg) then
!!    numold = num
!!    denold = den
!!  else
!!    fcoeff = fcoeff
!!  end if
!!end if
!
!call fcoeff_d (max_res_gen, nzone, fcoeff, num, den)
!
!if (time_avg) then
!
!  call error (sub_name, 'time_avg disabled during code transition')
!  
!!  
!!  !--determine time scale for time averaging
!!  call set_timescale_fcoeff (timescale, num, den)
!!
!!  !--protect from NaN:
!!  nan = .false.
!!      
!!  $if ($IFC || $IFORT)
!!    if (nd_time_avg) then
!!      nan = (isnan (timescale) .or. isnan (num) .or. isnan (den))
!!    else
!!      nan = (isnan (timescale) .or. isnan (fcoeff))
!!    end if
!!  $endif
!!  $if ($XLF)
!!    if (ieee_support_datatype (1._rp)) then
!!      if (ieee_support_nan (1._rp)) then
!!        if (nd_time_avg) then
!!          nan = (ieee_is_nan (timescale) .or. ieee_is_nan (num) .or.  &
!!                 ieee_is_nan (den))
!!        else
!!          nan = (ieee_is_nan (timescale) .or. ieee_is_nan (fcoeff))
!!        end if
!!      end if
!!    end if
!!  $endif
!!   
!!  if (any (nan)) then
!!    
!!    call mesg (sub_name, 'encountered NAN:')
!!    call mesg (sub_name, 'timescale=', timescale)
!!      
!!    if (nd_time_avg) then
!!      call mesg (sub_name, 'num=', num)
!!      call mesg (sub_name, 'den=', den)
!!    else
!!      call mesg (sub_name, 'fcoeff=', fcoeff)
!!    end if
!!      
!!  end if
!!
!!  do z = 1, nzone
!!
!!    if (.not. nan(z)) then
!!  
!!      wgt(z) = (nCDupdate * dt) / timescale(z)
!!      !wgt = (nCDupdate * dt) / (timescale(z) + nCDupdate*dt)
!!    
!!      !--weighted time avg. as in Lagrangian dynamic model
!!      if (nd_time_avg) then
!!        num(z) = wgt(z) * num(z) + (1._rp - wgt(z)) * numold(z)
!!        den(z) = wgt(z) * den(z) + (1._rp - wgt(z)) * denold(z)
!!        fcoeff(:, z) = num(:, z) / den(:, z)
!!      else
!!        fcoeff(:, z) = wgt(z) * fcoeff(z) + (1._rp - wgt(z)) * CDold(z)
!!      end if
!!      
!!    else  !--nan
!!    
!!      if (nd_time_avg) then
!!        num(z) = numold(z)
!!        den(z) = denold(z)
!!      end if
!!      
!!      CDdyn(z) = CDold(z) 
!!      timescale(z) = huge (1._rp)  !--for purpose of writing to file only
!!      
!!    end if
!!    
!!  end do
!!
!else
!
!  timescale = -abs (BOGUS)  !--since not time averaing
!  
!end if
!
!if (VERBOSE) call exit_sub (sub_name)
!
!end subroutine fcoeff_d_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--sets velscale to value provided in argument
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine set_velscale_br ( br, velscale )
implicit none

type (branch_type), intent (inout) :: br
real (rp), intent (in) :: velscale(nd)

character (*), parameter :: sub_name = mod_name // '.set_velscale_br'

!logical, parameter :: DEBUG = .true.

character (64) :: sn

integer :: i

!---------------------------------------------------------------------

sn = trim ( sub_name ) // trim ( chcoord )

$if ($VERBOSE)
call enter_sub ( sn )
$endif

br % velscale = velscale

$if ($DEBUG)
if ( DEBUG ) then
    call mesg ( sn, 'br % ident=', br % ident )
    call mesg ( sn, 'br % velscale=', br % velscale )
end if
$endif

if ( associated ( br % sub_branch ) ) then
    do i = 1, br % n_sub_branch
        call set_velscale_br ( br % sub_branch(i), velscale )
    end do
end if

$if ($VERBOSE)
call exit_sub ( sn )
$endif
    
end subroutine set_velscale_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--possible major problem: right now, this is called at a stage when
!  the velocities are not really valid: they are partially updated
!  and still need to be projected.  MAY NEED TO SPLIT UP TREE CODE
!  TO ACCESS VELOCITYS FROM n-1 BEFORE UPDATE BEGINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_velscale_ta ()
$if ($MPI)
  use param, only : up, down, status  !--plus those at top of module
  !use sim_param, only : w
$endif
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_velscale_ta'

integer, parameter :: tag = 333

logical, parameter :: use_global_velscale = .true.
    !--uses one velscale for all branches: mean velocity in whole domain

!logical, parameter :: DEBUG = .true.

integer :: i

logical, save :: test_velscale_mask = .false.

real (rp) :: Uinf(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name // trim (chcoord))
$endif

$if ($MPI)
  !--resync w(nz) <-> w(1') here
  !--this may not always be required, depends on main code
  !--at top process, w(nz) will be BOGUS
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag,   &
                     comm, status, ierr)                     
$endif

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (u(:, :, 1:nz), 'fill_velscale_ta.u')
  call DEBUG_write (v(:, :, 1:nz), 'fill_velscale_ta.v')
  call DEBUG_write (w(:, :, 1:nz), 'fill_velscale_ta.w')
end if
$endif

if ( use_global_velscale ) then

    call calc_Uinf ( Uinf )
    
    do i = 1, n_tree
        call set_velscale_br ( tree_array(i) % trunk,  Uinf )
    end do

else  !--use a locally defined branch velocity scale

    do i = 1, n_tree

        call fill_velscale (tree_array(i) % trunk)
    
        if ( test_velscale_mask )  &
            call write_velscale_mask (tree_array(i) % trunk)
    
    end do

    if (test_velscale_mask) test_velscale_mask = .false.  !--only test once

    $if ($MPI)

        call mpi_allreduce (nvelscale, nvelscale_sum, nbr, MPI_INTEGER,  &
                            MPI_SUM, comm, ierr)
        call mpi_allreduce (velscale, velscale_sum, nd * nbr, MPI_RPREC,  &
                            MPI_SUM, comm, ierr)

        !--now perform division for the averaging
        do i = 1, nbr
        
            if (nvelscale_sum(i) > 0) then
                velscale_sum(:, i) = velscale_sum(:, i) / nvelscale_sum(i)
            else
                !--nvelscale_sum will be zero when the branch is not at 
                !  max_res_gen or ngen, so it is not always an error
                !--actually, this has changed: velscale should now be valid
                !  at all generations
                !call error (sub_name,  &
                !             'nvelscale_sum is zero for br % ident =', i)
                velscale_sum(:, i) = BOGUS  !--make sure it is not used
            end if

        end do

        do i = 1, n_tree
            call repack_br (tree_array(i) % trunk, velscale_sum, 'velscale')
        end do

    $endif

end if

$if ($VERBOSE)
 call exit_sub (sub_name // trim (chcoord))
$endif

end subroutine fill_velscale_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--only fills gen >= max_res_gen
!--for now, fills all gen (see commented parts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_velscale (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_velscale'

integer :: i

!---------------------------------------------------------------------

!--allow velscale to be calculated for all scales (for a priori)
!  can provide optional argument version if this is too slow
!if (br % gen >= tree_array(br % itree) % max_res_gen) then

  call velscale_br (br)

!end if

if (br % gen < tree_array(br % itree) % n_gen) then

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_velscale (br % sub_branch(i))
  end do

end if

end subroutine fill_velscale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--averages velocity over bounding box
!--MPI: only coord 0 has valid velscale after this
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine velscale_br (br)
use sim_param, only : u, v, w
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.velscale_br'
character (*), parameter :: f_vmask_pre = 'vel_mask.'
character (*), parameter :: f_vmask_suf = '.dat'

integer, parameter :: lun = 1  !--for velocity mask

logical, parameter :: include_solid_pts = .true.
logical, parameter :: use_shifted_difference = .false.
logical, parameter :: write_vmask = .false.

character (128) :: fname

integer :: ipt
integer :: i, j, k
integer :: iv, jv, kv
integer :: n_f  !--number of fluid pts in bbox, could save in needed
integer :: shft(nd)  !--for use with shifted difference

logical :: opn, exst

real (rp) :: vel(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sub_name )
$endif
!--option check
if ( use_local_vel .and. ( .not. include_solid_pts ) )  &
    call error ( sub_name, 'use_local_vel=T requires include_solid_pts=T' )

if ( .not. associated ( br % bboxpt ) ) call bboxpt_init ( br )
    !--after this all bbox will be associated
    !--some procs may have 0-sized arrays

if ( use_shifted_difference ) then
    shft = (/ -nint( 0.25_rp * br % width_bbox / dx ), 0, 0 /)
        !--assumes flow is in positive x-direction
        !--the shift here sort-of agrees with the upstream dimension
        !  in the udiam, and Y velocity scales
else
    shft = 0
end if

if ( write_vmask ) then
    inquire ( lun, opened=opn, exist=exst )
    if ( (.not. exst) .or. opn ) call error ( sub_name, 'problem with lun' )

    write ( fname, '(a,i0,a)' ) f_vmask_pre // 'ident',  &
                                br % ident, f_vmask_suf

    $if ($MPI)
        write ( fname, '(a,i0,a,i0)' ) f_vmask_pre // 'ident', br % ident,  &
                                       f_vmask_suf // '.c', coord
    $else
        write ( fname, '(a,i0,a)' ) f_vmask_pre // 'ident', br % ident,  &
                                    f_vmask_suf
    $endif

    inquire (file=fname, opened=opn)
    if ( opn ) call error ( sub_name, 'file already open' )

    open (lun, file=fname)
        !--no header for now, just scatter plot data
end if

vel = 0.0_rp
n_f = 0

!--MPI: some processes may have nbboxpt == 0, so they will skip this
do ipt = 1, br % nbboxpt

    i = br % bboxpt(1, ipt)
    j = br % bboxpt(2, ipt)
    k = br % bboxpt(3, ipt)

    $if ($MPI)
        if ((coord == nproc - 1) .and. (k == nz - 1)) then
            call error (sub_name, 'bboxpt trying to access w@k=nz' // n_l //  &
                                  'domain is probably too small to fit bbox')
        end if
    $endif

    if (.not. velscale_mask (i, j, k, br)) cycle  !--skip to next point

    !--point where velocity will be sampled
    iv = i + shft(1)
    jv = j + shft(2)
    kv = k + shft(3)

    !--now we must check that iv, jv, kv are still in the domain
    !--this does NOT apply periodic BC
    if (( iv < 1 .or. iv > nx ) .or. ( jv < 1 .or. jv > ny ) .or.         &
        ( kv < 1 .or. kv > nz-1 ) )                                      &
        call error ( sub_name, 'velocity point out of bounds' // n_l //  &
                               '(iv,jv,kv)=', (/ iv, jv, kv /) )

    if ( use_shifted_difference ) then
        if ( velscale_mask (iv, jv, kv, br) ) cycle
            !--excludes points if their shifted image is still inside
            !  the original velscale_mask
            !--hopefully shft is non-zero, or else we will not
            !  select any points
    end if

    !--could have made bboxpt consist of pts with phi > 0, but that
    !  is not convenient for other parts of code (unres. force distrib)
    if (.not. include_solid_pts) then
        if (phi(iv, jv, kv) <= 0._rp) cycle  !--skip to next point
    end if

    n_f = n_f + 1

    !--(k+1) is in bounds since max(k) = nz - 1 (see bboxpt_init)
    !--MPI: w(k+1) may be BOGUS
    vel = vel + (/ u(iv, jv, kv), v(iv, jv, kv),                 &
                   0.5_rp * (w(iv, jv, kv) + w(iv, jv, kv + 1)) /)

    if ( write_vmask )  &
        write (lun, '(3(i0,1x))' ) iv, jv, kv

end do

if ( write_vmask ) close (lun)

$if ($MPI)

  !--we do not divide here, this is done after reduction in 
  !  fill_velscale_ta
  velscale(:, br % ident) = vel  !--just a partial sum
  nvelscale(br % ident) = n_f    !--just a partial total

$else

  if (n_f > 0) then
    br % velscale = vel / n_f
  else  !--non-MPI case: all branches should have some bbox pts
    call error (sub_name, 'no fluid bboxpt for br % ident =', br % ident)
  end if
  
$endif

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine velscale_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine write_velscale_mask (br)
implicit none

type (branch_type), intent (inout) :: br  !--out since bbox may be init

character (*), parameter :: sub_name = mod_name // '.write_velscale_mask'

integer :: i

!---------------------------------------------------------------------

call write_velscale_mask_br (br)

if (br % gen < tree_array(br % itree) % n_gen) then

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call write_velscale_mask (br % sub_branch(i))
  end do

end if

end subroutine write_velscale_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this write this mask for this branch (1 = true, 0 = false)
!  for MPI, the whole mask for a branch may be spread across several files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_velscale_mask_br (br)
implicit none

type (branch_type), intent (inout) :: br  !--out, since may init bbox

character (*), parameter :: sub_name = mod_name // '.write_velscale_mask_br'
character (*), parameter :: fname_base = 'velscale_mask.ident'

integer, parameter :: lun = 1

character (128) :: fname

integer :: mask_value
integer :: i, j, k
integer :: ipt

logical :: ext, opn

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

!--need to check init in case velscale_br has not been called yet
if (.not. associated (br % bboxpt)) call bboxpt_init (br)
                                    !--after this all bbox will be associated
                                    !--some procs may have 0-sized arrays

if (.not. USE_MPI) then
  write (fname, '(a,i0,a)') fname_base, br % ident, '.dat'
else
  write (fname, '(a,i0,a,i0)') fname_base, br % ident, '.dat.c', coord
end if

inquire (unit=lun, opened=opn, exist=ext)
if (.not. ext .or. opn) call error (sub_name, 'lun open or nonexistant')

inquire (file=fname, opened=opn)
if (opn) call error (sub_name, 'file ' // trim (fname) // 'already open')

open (lun, file=fname, action='write', position='rewind')
write (lun, '(a)') 'variables = "i" "j" "k" "mask_value"'
write (lun, '(a,i0)') 'zone, f=point, i=', br % nbboxpt

!--MPI: some processes may have nbboxpt == 0, so they will skip this
do ipt = 1, br % nbboxpt

  i = br % bboxpt(1, ipt)
  j = br % bboxpt(2, ipt)
  k = br % bboxpt(3, ipt)

  $if ($MPI)
    if ((coord == nproc - 1) .and. (k == nz - 1)) then
      call error (sub_name, 'bboxpt trying to access w@k=nz' // n_l //  &
                            'domain is probably too small to fit bbox')
    end if
  $endif

  if (.not. velscale_mask (i, j, k, br)) then
    mask_value = 0
  else
    mask_value = 1
  end if

  write (lun, '(4(i0,1x))') i, j, k, mask_value

end do

close (lun)

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine write_velscale_mask_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--returns true when do want to include ipt, false when ipt to be exluded
!--may want to store mask so this calculation only needs to be done once
!  however, this is a first attempt so we do not (and the resulting code
!  is inefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function velscale_mask (i, j, k, br)
implicit none

logical :: velscale_mask

integer, intent (in) :: i, j, k
type (branch_type), intent (in) :: br

character (*), parameter :: sub_name = mod_name // '.velscale_mask'
character (*), parameter :: mask_type = 'Y_convex_hull_exact_upstream'
    !--'none', 'chop_top', 'upstream_bbox', 'upstream_diamond',
    !  'diamond', 'small_upstream_box', 'Y_box', 'Y_upstream_box'
    !  'Y_convex_hull_exact_upstream', 'Y_convex_hull_exact'

real (rp) :: r
real (rp) :: xmax, ymax, ymax_z
real (rp) :: zmin, zmax
real (rp) :: pivot(nd)  !--used as "center" to determine what is upstream
real (rp) :: x(nd)
real (rp) :: xi(nd), xi_p(nd)
real (rp) :: zeta(nd), zeta_p(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub ( sub_name )
$endif
    
!--do we need absolute k here?
x(1) = pt_of_grid ( i, 1, 1 )  !--u-nodes
x(2) = pt_of_grid ( j, 2, 1 )
x(3) = pt_of_grid ( k, 3, 1 )

xi = (x - br % x0)

xi_p(1) = dot_product ( xi, br % x_hat )
xi_p(2) = dot_product ( xi, br % y_hat )
xi_p(3) = dot_product ( xi, br % z_hat )

velscale_mask = .false.

select case ( mask_type )
case ( 'chop_top' )

    if ( xi_p(3) < (br % l) ) velscale_mask = .true.
  
case ( 'upstream_bbox' )

    !--the limits here are specific to 2d-cross geometry
    if ( xi_p(1) < 0.0_rp ) velscale_mask = .true.

case ( 'small_upstream_box' )

    !--the limits here are specific to 2d-cross geometry
    if ( ( xi_p(1) < 0._rp ) .and.                                &
         ( abs ( xi_p(2) ) < 0.25_rp * (br % width_bbox) ) .and.  &
         ( xi_p(3) < (br % l) )                                   &
       ) velscale_mask = .true.

case ( 'upstream_diamond' )

    zeta = xi - (br % l) * (br % z_hat)
      !--reset the "origin" to center of top of branch (no cap)
      
    zeta_p(1) = dot_product ( zeta, br % x_hat )
    zeta_p(2) = dot_product ( zeta, br % y_hat )
    zeta_p(3) = dot_product ( zeta, br % z_hat )
    
    if ( ( -0.25_rp * br % width_bbox <= zeta_p(1) ) .and.  &
         ( zeta_p(1) < 0.0_rp ) .and.                       &
         ( abs ( zeta_p(2) ) + abs ( zeta_p(3) ) <=         &
           0.5_rp * (br % width_bbox) )                     &
       ) velscale_mask = .true.

!    if ( ( zeta_p(1) < 0.0_rp ) .and.                 &
!         ( abs ( zeta_p(2) ) + abs ( zeta_p(3) ) <=   &
!           0.5_rp * (br % width_bbox) )               &
!       ) velscale_mask = .true.

case ( 'diamond' )

    zeta = xi - (br % l) * (br % z_hat)
        !--reset the "origin" to center of top of branch (no cap)
      
    zeta_p(1) = dot_product ( zeta, br % x_hat )
    zeta_p(2) = dot_product ( zeta, br % y_hat )
    zeta_p(3) = dot_product ( zeta, br % z_hat )
   
    !--includes any zeta_p(1)
    if ( abs ( zeta_p(2) ) + abs ( zeta_p(3) ) <= 0.5_rp * (br % width_bbox)  &
       ) velscale_mask = .true.

case ( 'Y_box' )

    r = tree_array ( br % itree ) % ratio

    !--coarse box that has branch tips on all four sides in 
    !  xi_(2), xi_p(3) plane
    !  l = br % l (****NO CAP****)
    !  zmax = max height above x0 = l + (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)
    !  zmin = min height above x0 \
    !       = l + (r*l+d/2) * (1-r**2) / sqrt(2) - r**3 * d/2 - r**4 * zmax
    !  ymax = half-width of box = (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)

    ymax = ( r * (br % l) + 0.5_rp * (br % d) ) *           &
           ( r + 1.0_rp / sqrt (2.0_rp) ) / ( 1.0_rp - r**2 )
    zmax = (br % l) + ymax
    zmin = (br % l) + ( r * (br % l) + 0.5_rp * (br % d) ) *  &
                      ( 1.0_rp - r**2 ) / sqrt (2.0_rp)       &
                    - 0.5_rp * (r**3 )* (br % d) - r**4 * zmax

    !--conforms to branch width in local x1 direction
    if ( ( abs ( xi_p(1) ) <= 0.5_rp * (br % d) ) .and.  &
         ( abs ( xi_p(2) ) <= ymax ) .and.               &
         ( zmin <= xi_p(3) .and. xi_p(3) <= zmax )       &
       ) velscale_mask = .true.
    
case ( 'Y_upstream_box' )

    r = tree_array ( br % itree ) % ratio

    !--coarse box that has branch tips on all four sides in 
    !  xi_(2), xi_p(3) plane
    !  l = br % l (no cap)
    !  zmax = max height above x0 = l + (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)
    !  zmin = min height above x0 \
    !       = l + (r*l+d/2) * (1-r**2) / sqrt(2) - r**3 * d/2 - r**4 * zmax
    !  ymax = half-width of box = (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)

    ymax = ( r * (br % l) + 0.5_rp * (br % d) ) *           &
           ( r + 1.0_rp / sqrt (2.0_rp) ) / ( 1.0_rp - r**2 )
    zmax = (br % l) + ymax
    zmin = (br % l) + ( r * (br % l) + 0.5_rp * (br % d) ) *  &
                      ( 1.0_rp - r**2 ) / sqrt (2.0_rp)       &
                    - 0.5_rp * (r**3 )* (br % d) - r**4 * zmax

    if ( ( -0.25_rp * br % width_bbox <= xi_p(1) ) .and.  &
         ( xi_p(1) < 0.0_rp ) .and.                       &
         ( abs ( xi_p(2) ) <= ymax ) .and.                &
         ( zmin <= xi_p(3) .and. xi_p(3) <= zmax )        &
       ) velscale_mask = .true.

case ( 'Y_convex_hull_exact' )
    !--this is only approximate convex hull at the moment
   
    !--takes upstream and downstream points
    r = tree_array ( br % itree ) % ratio

    !--coarse box that has branch tips on all four sides in 
    !  xi_(2), xi_p(3) plane
    !  l = br % l (no cap)
    !  zmax = max height above x0 = l + (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)
    !  zmin = min height above x0 \
    !       = l + (r*l+d/2) * (1-r**2) / sqrt(2) - r**3 * d/2 - r**4 * zmax
    !  ymax = half-width of box = (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)
   
    !xmax = 0.25_rp * br % width_bbox
    xmax = br % d  !--total width is then 2d
    ymax = ( r * (br % l) + 0.5_rp * (br % d) ) *           &
           ( r + 1.0_rp / sqrt (2.0_rp) ) / ( 1.0_rp - r**2 )
    zmax = (br % l) + ymax
    zmin = (br % l) + ( r * (br % l) + 0.5_rp * (br % d) ) *  &
                      ( 1.0_rp - r**2 ) / sqrt (2.0_rp)       &
                    - 0.5_rp * (r**3 )* (br % d) - r**4 * zmax
    
   if ( zmin <= xi_p(3) .and. xi_p(3) <= zmax ) then
    
       if ( ( abs ( xi_p(1) ) <= xmax ) .and.  &
            ( xi_p(1) < 0.0_rp ) .and.                       &
            ( abs ( xi_p(2) ) <= ymax )                      &
          ) velscale_mask = .true.
    
   else if ( 0.0_rp <= xi_p(3) .and. xi_p(3) < zmin ) then
    
       ymax_z = 0.5_rp * ( br % d ) +                           &
                ( ymax - 0.5_rp * (br % d) ) * ( xi_p(3) / zmin )
       if ( ( abs ( xi_p(1) ) <= 0.25_rp * br % width_bbox ) .and.  &
            ( xi_p(1) < 0.0_rp ) .and.                       &
            ( abs ( xi_p(2) ) <= ymax_z )                    &
          ) velscale_mask = .true.
    
   end if

case ( 'Y_convex_hull_exact_upstream' )
    !--this is only approximate convex hull at the moment
    
    pivot = (br % x0) + (br % l) * (br % z_hat)  !--in global coords
    if ( dot_product ( x - pivot, flow_dir ) <= 0.0_rp ) then

        r = tree_array ( br % itree ) % ratio

        !--coarse box that has branch tips on all four sides in 
        !  xi_(2), xi_p(3) plane
        !  l = br % l (no cap)
        !  zmax = max height above x0 = l + (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)
        !  zmin = min height above x0 \
        !       = l + (r*l+d/2) * (1-r**2) / sqrt(2) - r**3 * d/2 - r**4 * zmax
        !  ymax = half-width of box = (r*l+d/2) * (r+1/sqrt(2)) / (1-r**2)
   
        !xmax = 0.25_rp * br % width_bbox
        xmax = br % d  !--total width is then 2d
        ymax = ( r * (br % l) + 0.5_rp * (br % d) ) *           &
               ( r + 1.0_rp / sqrt (2.0_rp) ) / ( 1.0_rp - r**2 )
        zmax = (br % l) + ymax
        zmin = (br % l) + ( r * (br % l) + 0.5_rp * (br % d) ) *  &
                          ( 1.0_rp - r**2 ) / sqrt (2.0_rp)       &
                        - 0.5_rp * (r**3 )* (br % d) - r**4 * zmax
    
       if ( zmin <= xi_p(3) .and. xi_p(3) <= zmax ) then
    
           if ( ( abs ( xi_p(1) ) <= xmax ) .and.  &
                ( xi_p(1) < 0.0_rp ) .and.                       &
                ( abs ( xi_p(2) ) <= ymax )                      &
              ) velscale_mask = .true.
    
       else if ( 0.0_rp <= xi_p(3) .and. xi_p(3) < zmin ) then
    
           ymax_z = 0.5_rp * ( br % d ) +                           &
                    ( ymax - 0.5_rp * (br % d) ) * ( xi_p(3) / zmin )
           if ( ( abs ( xi_p(1) ) <= 0.25_rp * br % width_bbox ) .and.  &
                ( xi_p(1) < 0.0_rp ) .and.                       &
                ( abs ( xi_p(2) ) <= ymax_z )                    &
              ) velscale_mask = .true.
    
       end if

   end if

case ( 'none' )

    velscale_mask = .true.

case default

    call error ( sub_name, 'invalid mask_type=' // mask_type )
    
end select

$if ($VERBOSE)
call exit_sub ( sub_name )
$endif
    
end function velscale_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--Assumes the following branch variables are already set:
!  height_bbox, width_bbox, x0, x_hat, y_hat, z_hat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bboxpt_init (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.bboxpt_init'

integer :: pass
integer :: ipt
integer :: i, j, k

real (rp) :: x(nd), xp(nd)
real (rp) :: height, width

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

height = br % height_bbox
width = br % width_bbox

$if ($DEBUG)
if (DEBUG) then
  call mesg (sub_name, 'height =', height)
  call mesg (sub_name, 'width =', width)
end if
$endif

do pass = 1, 2

  ipt = 0

  do k = 1, nz - 1
    do j = 1, ny
      do i = 1, nx

        x(1) = pt_of_grid (i, 1, 1)  !--u-node
        x(2) = pt_of_grid (j, 2, 1) 
        x(3) = pt_of_grid (k, 3, 1)

        x = x - br % x0  !--relative to branch origin

        xp(1) = dot_product (x, br % x_hat)
        xp(2) = dot_product (x, br % y_hat)
        xp(3) = dot_product (x, br % z_hat)

        !--the bounding box is defined in branch-local coordinates
        !--inexact nature of floating point may result in ragged
        !  surfaces of our box--e.g. this occurs in the case where
        !  we have applied a twist to a branch and the trig errors
        !  result in x_hat, y_hat, z_hat being slightly out of alignment
        if ((abs (xp(1)) <= 0.5_rp * width) .and.      &
            (abs (xp(2)) <= 0.5_rp * width) .and.      &
            (0._rp <= xp(3) .and. xp(3) <= height)) then

          !--take all points, independent of value of phi
          ipt = ipt + 1
          if (pass == 2) br % bboxpt(:, ipt) = (/ i, j, k /)

        end if
        
      end do
    end do
  end do

  if (pass == 1) then
  
    br % nbboxpt = ipt

    $if ($VERBOSE)
    call mesg (sub_name, 'br % ident = ', br % ident)
    call mesg (sub_name, 'set br % nbboxpt = ', br % nbboxpt)
    $endif

    if (.not. associated (br % bboxpt)) then
      allocate (br % bboxpt(nd, br % nbboxpt))
    else
      call error (sub_name, 'bboxpt is already associated')
    end if
    
  end if

end do

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine bboxpt_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_resf_ta ()
use immersedbc, only : fx, fy, fz
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_resf_ta'

!logical, parameter :: DEBUG = .true.

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name // trim (chcoord))
$endif

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (fx(:, :, 1:nz), 'fill_resf_ta.fx')
  call DEBUG_write (fy(:, :, 1:nz), 'fill_resf_ta.fy')
  call DEBUG_write (fz(:, :, 1:nz), 'fill_resf_ta.fz')
end if
$endif

do i = 1, n_tree
  call fill_resf (tree_array(i) % trunk)
end do

$if ($MPI)

  call mpi_allreduce (resf, resf_sum, nd * nbr, MPI_RPREC, MPI_SUM,  &
                      comm, ierr)

  do i = 1, n_tree
    call repack_br (tree_array(i) % trunk, resf_sum, 'resf')
  end do
  
$endif

$if ($VERBOSE)
 call exit_sub (sub_name // trim (chcoord))
$endif

end subroutine fill_resf_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this fills resf for all generations
!--could add optional arg to specify a particular generation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_resf (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_resf'

integer :: i

!---------------------------------------------------------------------

call resf_br (br)  !--assumes br is resolved

!--only recur over sub-branches if they are resolved
if (br % gen < tree_array(br % itree) % max_res_gen) then

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_resf (br % sub_branch(i))
  end do
  
end if

end subroutine fill_resf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--calculate resolved force on each resolved branch
!--br must be resolved else issue error
!--multiplied by dV here
!--MPI: only coord 0 will have valid resf
!--probably better to do one loop over domain, and somehow directly
!  address each branch (array of derived type, whose only data are pointers?)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine resf_br (br)
use immersedbc, only : fx, fy, fz
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.resf_br'

integer :: i, j, k

real (rp) :: rf(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
 call enter_sub (sub_name)
$endif

if (.not. br % resolved) call error (sub_name, 'expecting resolved branch')

rf = 0._rp

!--this is not very efficient
!--probably better to replace with a loop over bboxpt, or
!  do all branches at once with one loop over nx, ny, nz
do k = 1, nz - 1
  do j = 1, ny
    do i = 1, nx

      if (brindex(i, j, k) == br % ident) then
        !--need to refine this for w-nodes: brindex is for u-nodes
        rf(1) = rf(1) + fx(i, j, k)
        rf(2) = rf(2) + fy(i, j, k)
        rf(3) = rf(3) + fz(i, j, k)
      end if
      
    end do
  end do
end do

rf = rf * (dx * dy * dz)  !--mult by dV

$if ($MPI)

  resf(:, br % ident) = rf  !--just partial sum, see fill_resf_ta

$else

  br % resf = rf

$endif

$if ($DEBUG)
if (DEBUG) then
  call mesg (sub_name, 'br % ident =', br % ident)
  call mesg (sub_name, 'br % resf =', br % resf)
end if
$endif

$if ($VERBOSE)
 call exit_sub (sub_name)
$endif

end subroutine resf_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine brindex_init ()
use param, only : iBOGUS
implicit none

character (*), parameter :: sub_name = mod_name // '.brindex_init'
character (*), parameter :: fbrindex_in = 'brindex.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fbrindex_in_MPI
$endif

integer :: ip

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fbrindex_in_MPI, '(a,a,i0)') fbrindex_in, MPI_suffix, coord
    
  inquire (file=fbrindex_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fbrindex_in_MPI)

  open (1, file=fbrindex_in_MPI, action='read', position='rewind',  &
         form='unformatted')
  read (1) brindex(:, :, 1:nz-1)
  close (1)

  brindex(:, :, nz) = iBOGUS

$else

  inquire (file=fbrindex_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fbrindex_in)

  open (1, file=fbrindex_in, action='read', position='rewind',  &
         form='unformatted')
  read (1) brindex
  close (1)

$endif

brindex_initialized = .true.

end subroutine brindex_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_ls
