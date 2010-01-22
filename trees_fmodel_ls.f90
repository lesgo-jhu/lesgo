!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this is supposed to be largely independent of param and sim_param
!  (and the main code), but this is not entirely possible since small
!  from param will be used via trees_base_ls
!--it is supposed to offer a lightweight way to do apriori test
!  (see trees_apri_ls)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module trees_fmodel_ls
use trees_base_ls
use messages
implicit none

save
$if ($TEST)
  public
$else
  private
$endif

public :: clip_fcoeff
public :: def_d_dir, def_dls_dir, def_nba_dir
public :: do_tavg_fcoeff, do_tavg_fdist
public :: fcoeff_ta, fcoeff_d_germano
public :: fcoeff, fcoeff_clipped
public :: fcoeff_divide
public :: fill_unresf_br
public :: vel_d, vel_dls, vel_nba
!public :: write_dls_line, write_nba_line

character (*), parameter :: mod_name = 'trees_fmodel_ls'

character (*), parameter :: vel_opt = 'component'
                            !--selects form of velocity in the model
                            !--values: 'absolute', 'component', 'mixed'

logical, parameter :: fcoeff_gen_avg = .true.
logical, parameter :: do_tavg_fcoeff = .false.
logical, parameter :: do_tavg_fdist = do_tavg_fcoeff

logical :: fcoeff_clipped(nfcoeff, nzone)  !--nfcoeff X nzone

real (rp) :: fcoeff(nfcoeff, nzone)  !--nfcoeff X nzone
             !--only used with averaging

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this assumes fcoeff has been updated
!--this assumes br is unresolved
!--this must be totally consistent with corresponding routine that
!  calculates fcoeff
!--fill_unresf_br parent routines in trees_ls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_unresf_br (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.fill_unresf_br'

!logical, parameter :: DEBUG = .true.

integer :: i

!real (rp) :: dir(nd, nfcoeff)  !--may get rid of br % fcoeff_dir
real (rp) :: cav2(nd)
real (rp) :: vmag(nfcoeff)
real (rp) :: v(nd, nfcoeff)

!---------------------------------------------------------------------

!--def_XXX_area only needed once per branch--if the area is flow-independent
select case (fmodel)
case ('d', 'd_germano')
  call def_d_dir (br, br % velscale, br % fcoeff_dir)
  call def_d_area (br)  !--define area
      !--this is only needed once per branch--if the area is flow-independent
  call vel_d (br % velscale, vmag, v)  !--no directions required here
case ('dls')
  call def_dls_dir (br, br % velscale, br % fcoeff_dir)
  call def_dls_area (br)  !--define area
  call vel_dls (br % velscale, br % fcoeff_dir, vmag, v)
case ('nba')
  call def_nba_dir (br, br % velscale, br % fcoeff_dir)
  call def_nba_area (br)  !--define area
  call vel_nba (br % velscale, br % fcoeff_dir, vmag, v)
case default
  call error (sub_name, 'invalid fmodel')
end select

cav2 = 0._rp
do i = 1, nfcoeff
  cav2 = cav2 + ( fcoeff(i, br % zone) * (br % A(i)) * vmag(i) * v(:, i) )
end do

br % unresf = -0.5_rp * cav2

end subroutine fill_unresf_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--resf at max_res_gen MUST have been updated before calling this
!--assumes max_res_gen is same for all trees
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fcoeff_d_germano ()
use linear_simple, only : solve_linear
implicit none

character (*), parameter :: sub_name = mod_name // '.fcoeff_germano'

integer :: i

real (rp) :: LHSmatrix(nzone, nzone)
real (rp) :: RHSvector(nzone)

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (nfcoeff /= 1) then
  call error (sub_name, 'expecting nfcoeff = 1')
end if

!--fill Mdyn
call fill_Mdyn_ta ()

!--set up linear system
LHSmatrix = 0._rp
RHSvector = 0._rp

do i = 1, n_tree
  call fill_LHSmatrix (tree_array(i) % trunk, LHSmatrix)
  call fill_RHSvector (tree_array(i) % trunk, RHSvector)
end do

!--solve the system for dynamic CD
call solve_linear (LHSmatrix, RHSvector, fcoeff(1, :))

$if ($DEBUG)
if (DEBUG) then
  call mesg (sub_name, 'LHSmatrix =', pack (LHSmatrix, mask=.true.))
  call mesg (sub_name, 'RHSvector =', RHSvector)
  call mesg (sub_name, '(unclipped) fcoeff =', pack (fcoeff, mask=.true.))
end if
$endif

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine fcoeff_d_germano

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--RHSvector may contain valuable info on entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_RHSvector (br, RHSvector)
implicit none

type (branch_type), intent (in) :: br
real (rp), intent (inout) :: RHSvector(nzone)

character (*), parameter :: sub_name = mod_name // '.fill_RHSvector'

integer :: iz
integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (br % gen == tree_array(br % itree) % max_res_gen) then

  do iz = 1, nzone

    !--dV is in resf already
    RHSvector(iz) = RHSvector(iz) - 2._rp * dot_product (br % resf,      &
                                                         br % Mdyn(:, iz))

  end do

else if (br % gen < tree_array(br % itree) % max_res_gen) then

  !--traverse sub-branches

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_RHSvector (br % sub_branch(i), RHSvector)
  end do

else

  call error (sub_name, 'invalid gen for this routine')

end if

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine fill_RHSvector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--LHSmatrix may contain valuable info on entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_LHSmatrix (br, LHSmatrix)

type (branch_type), intent (in) :: br
real (rp), intent (inout) :: LHSmatrix(nzone, nzone)

character (*), parameter :: sub_name = mod_name // '.fill_LHSmatrix'

integer :: i, j

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (br % gen == tree_array(br % itree) % max_res_gen) then

  !--could take advantage of symmetry
  do j = 1, nzone
    do i = 1, nzone

      LHSmatrix(i, j) = LHSmatrix(i, j) + dot_product (br % Mdyn(:, i),  &
                                                       br % Mdyn(:, j))
    end do
  end do
  
else if (br % gen < tree_array(br % itree) % max_res_gen) then

  !--traverse sub-branches

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub-branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_LHSmatrix (br % sub_branch(i), LHSmatrix)
  end do

else

  call error (sub_name, 'invalid gen for this routine')

end if
  
$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine fill_LHSmatrix 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_Mdyn_ta ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_Mdyn_ta'

integer :: i

!---------------------------------------------------------------------

do i = 1, n_tree
  call fill_Mdyn (tree_array(i) % trunk)
end do

end subroutine fill_Mdyn_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fill_Mdyn (br)
implicit none

type (branch_type) :: br

character (*), parameter :: sub_name = mod_name // '.fill_Mdyn'

integer :: i

!---------------------------------------------------------------------

if (br % gen == tree_array(br % itree) % max_res_gen) then

  call Mdyn_br (br)
  
else if (br % gen < tree_array(br % itree) % max_res_gen) then

  !--traverse sub-branches

  if (.not. associated (br % sub_branch)) then
    call error (sub_name, 'expecting associated sub_branches')
  end if

  do i = 1, br % n_sub_branch
    call fill_Mdyn (br % sub_branch(i))
  end do

else

   call error (sub_name, 'invalid gen for this routine')

end if

end subroutine fill_Mdyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--fills br % Mdyn(:, :)
!--br must be at max_res_gen
!--MPI: only proc w/ valid velscale should call this (i.e. coord 0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mdyn_br (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.Mdyn_br'

integer :: iz
integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (br % gen /= tree_array (br % itree) % max_res_gen) then
  call error (sub_name, 'expecting branch at max_res_gen')
end if

call def_d_area (br)

br % Mdyn(:, :) = 0._rp

!--the delta function part for this branch
br % Mdyn(:, br % zone) = mag (br % velscale) * (br % velscale) *  &
                          (br % A(1))

!do iz = 1, nzone
!
!  !--the delta function part for this branch
!  if (br % zone == iz) then
!    br % Mdyn(:, iz) = mag (br % velscale) * (br % velscale) * (br % Ap_bbox)
!  else
!    br % Mdyn(:, iz) = 0._rp
!  end if
!
!end do

!--sub-branch part
if (.not. associated (br % sub_branch)) then
  call error (sub_name, 'expecting associated sub_branches')
end if

do i = 1, br % n_sub_branch

  iz = br % sub_branch(i) % zone

  br % Mdyn(:, iz) = br % Mdyn(:, iz) -                     &
                     mag (br % sub_branch(i) % velscale) *  &
                     (br % sub_branch(i) % velscale) *      &
                     (br % sub_branch(i) % A(1))

end do

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine Mdyn_br
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--does not make use of globals, so that it can be used to clip
!  branch-local coefficients as well
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine clip_fcoeff (nc, c, active)
implicit none

integer, intent (in) :: nc  !--number of force coefficients
real (rp), intent (inout) :: c(nc)  !--force to be 1d array
logical, intent (out), optional :: active(nc)  !--indicates clipping

character (*), parameter :: sub_name = mod_name // '.clip_fcoeff'

real (rp), parameter :: fcoeff_clip_min = 0._rp
real (rp), parameter :: fcoeff_clip_max = 10000._rp

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

active = .false.

do i = 1, nc

  if (c(i) < fcoeff_clip_min) then
    c(i) = fcoeff_clip_min
    active(i) = .true.
  end if
  
  if (c(i) > fcoeff_clip_max) then
    c(i) = fcoeff_clip_max
    active(i) = .true.
  end if
  
end do

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine clip_fcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine write_nba_line (jt, gen, lun)
!implicit none
!
!integer, intent (in) :: jt, gen, lun
!
!character (*), parameter :: sub_name = mod_name // '.write_nba_line'
!
!integer :: i
!
!logical :: opn
!
!!---------------------------------------------------------------------
!
!if (VERBOSE) call enter_sub (sub_name)
!
!inquire (unit=lun, opened=opn)
!if (.not. opn) then
!  call error (sub_name, 'lun not open')
!end if
!
!do i = 1, n_tree
!
!  if (gen <= tree_array(i) % n_gen) then
!    call write_nba_line_br (jt, gen, lun, tree_array(i) % trunk)
!  else
!    call error (sub_name, 'gen > n_gen for tree i=', i)
!  end if
!  
!end do
!
!if (VERBOSE) call exit_sub (sub_name)
!
!end subroutine write_nba_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine write_dls_line (jt, gen, lun)
!implicit none
!
!integer, intent (in) :: jt, gen, lun
!
!character (*), parameter :: sub_name = mod_name // '.write_dls_line'
!
!integer :: i
!
!logical :: opn
!
!!---------------------------------------------------------------------
!
!if (VERBOSE) call enter_sub (sub_name)
!
!inquire (unit=lun, opened=opn)
!if (.not. opn) then
!  call error (sub_name, 'lun not open')
!end if
!
!do i = 1, n_tree
!
!  if (gen <= tree_array(i) % n_gen) then
!    call write_dls_line_br (jt, gen, lun, tree_array(i) % trunk)
!  else
!    call error (sub_name, 'gen > n_gen for tree i=', i)
!  end if
!  
!end do
!
!if (VERBOSE) call exit_sub (sub_name)
!
!end subroutine write_dls_line
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!recursive subroutine write_nba_line_br (jt, gen, lun, br)
!implicit none
!
!integer, intent (in) :: jt, gen, lun
!type (branch_type), intent (in) :: br
!
!character (*), parameter :: sub_name = mod_name // '.write_nba_line_br'
!character (*), parameter :: fmt = '(3(i0,1x),6(es13.6,1x))'
!
!integer :: i
!
!!---------------------------------------------------------------------
!
!if (br % gen == gen) then
!
!  write (lun, fmt) jt, br % gen, br % ident,                   &
!                   br % Cnormal, br % Cbinormal, br % Caxial,  &
!                   br % Fnormal, br % Fbinormal, br % Faxial
!
!else if (br % gen < gen) then
!
!  if (associated (br % sub_branch)) then
!
!    do i = 1, br % n_sub_branch
!      call write_nba_line_br (jt, gen, lun, br % sub_branch(i))
!    end do
!    
!  end if
!
!else
!
!  call error (sub_name, 'unexpected condition br % gen > gen')
!
!end if
!
!end subroutine write_nba_line_br
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!recursive subroutine write_dls_line_br (jt, gen, lun, br)
!implicit none
!
!integer, intent (in) :: jt, gen, lun
!type (branch_type), intent (in) :: br
!
!character (*), parameter :: sub_name = mod_name // '.write_dls_line_br'
!character (*), parameter :: fmt = '(3(i0,1x),6(es13.6,1x))'
!
!integer :: i
!
!!---------------------------------------------------------------------
!
!if (br % gen == gen) then
!
!  write (lun, fmt) jt, br % gen, br % ident,  &
!                   br % Cdrag, br % Clift, br % Cside,  &
!                   br % Fdrag, br % Flift, br % Fside
!
!else if (br % gen < gen) then
!
!  if (associated (br % sub_branch)) then
!
!    do i = 1, br % n_sub_branch
!      call write_dls_line_br (jt, gen, lun, br % sub_branch(i))
!    end do
!    
!  end if
!
!else
!
!  call error (sub_name, 'unexpected condition br % gen > gen')
!
!end if
!
!end subroutine write_dls_line_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--main reason for this routine is just to keep the method for avoiding
!  division by zero consistent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fcoeff_divide (num, den)
implicit none

real (rp), intent (in) :: num(nfcoeff, nzone), den(nfcoeff, nzone)

character (*), parameter :: sub_name = mod_name // '.fcoeff_divide'

real (rp), parameter :: eps = epsilon (eps)

!---------------------------------------------------------------------

where (abs (den) > eps)
  fcoeff = num / den
elsewhere
  fcoeff = 0._rp  !--clip when denominator is near zero
                  !  den will be near zero when vels are near zero
                  !  so c A v^2 will should also be near zero
                  !  so clipping this way is consistent
end where

end subroutine fcoeff_divide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fcoeff_ta (gen, num, den)
implicit none

integer, intent (in) :: gen
real (rp), intent (out) :: num(nfcoeff, nzone), den(nfcoeff, nzone)
                               
character (*), parameter :: sub_name = mod_name // '.fcoeff_ta'

real (rp), parameter :: eps = epsilon (eps)

character (32) :: fmt

integer :: i
integer :: nbr(nzone)  !--counter for number of branches in each zone

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

!--only used in the averaging process
nbr = 0
num = 0._rp
den = 0._rp

do i = 1, n_tree

  if (gen <= tree_array(i) % n_gen) then
 
    call fcoeff_br (tree_array(i) % trunk, gen, nbr, num, den)
    
  else
  
    call error (sub_name, 'gen > n_gen for tree i=', i)
    
  end if
  
end do

!--perform check to make sure there was a branch in each zone
if (any (nbr <= 0)) then  !--< just to be safe

  write (fmt, '(a,i0,a)') '(a,i0,a,', nzone, '(i0,1x))'
  write (msg, fmt) 'zero branches in at least one zone' // n_l //  &
                   'branch generation =', gen, n_l // 'nbr =', nbr
                   
  call error (sub_name, msg)

end if

!--calculate zonal averages
do i = 1, nfcoeff
  num(i, :) = num(i, :) / nbr(:)
  den(i, :) = den(i, :) / nbr(:)
end do

!--calculate average coefficients
call fcoeff_divide (num, den)

!--set coefficient at each branch to the average
!do i = 1, n_tree
!  call set_fcoeff_nba_br (gen, tree_array(i) % trunk, fcoeff)
!  !--set the unresolved branch level coefficients too!!!
!  call set_fcoeff_nba_br (gen+1, tree_array(i) % trunk, fcoeff)
!end do

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine fcoeff_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--this routine will set internal branch n,b,a coefficients
!!--the value of the coefficients can be set using the branch-local value
!!  or the average over a generation (depending on options below)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine fcoeff_nba (gen, num, den)
!implicit none
!
!integer, intent (in) :: gen
!real (rp), intent (out) :: num(nfcoeff, nzone), den(nfcoeff, nzone)
!                               !--nfcoeff must be 3 here
!                               
!character (*), parameter :: sub_name = mod_name // '.fcoeff_nba'
!
!real (rp), parameter :: eps = epsilon (eps)
!
!character (32) :: fmt
!
!integer :: i
!integer :: nbr(nzone)  !--counter for number of branches in each zone
!
!!---------------------------------------------------------------------
!
!if (VERBOSE) call enter_sub (sub_name)
!
!!--only used in the averaging process
!nbr = 0
!num = 0._rp
!den = 0._rp
!
!do i = 1, n_tree
!
!  if (gen <= tree_array(i) % n_gen) then
! 
!    call fcoeff_nba_br (tree_array(i) % trunk, gen, nbr, num, den)
!    
!  else
!  
!    call error (sub_name, 'gen > n_gen for tree i=', i)
!    
!  end if
!  
!end do
!
!!--perform check to make sure there was a branch in each zone
!if (any (nbr <= 0)) then  !--< just to be safe
!
!  write (fmt, '(a,i0,a)') '(a,i0,a,', nzone, '(i0,1x))'
!  write (msg, fmt) 'zero branches in at least one zone' // n_l //  &
!                   'branch generation =', gen, n_l // 'nbr =', nbr
!                   
!  call error (sub_name, msg)
!
!end if
!
!!--calculate zonal averages
!do i = 1, nfcoeff
!  num(i, :) = num(i, :) / nbr(:)
!  den(i, :) = den(i, :) / nbr(:)
!end do
!
!!--calculate average coefficients
!call fcoeff_divide (num, den)
!
!!--set coefficient at each branch to the average
!!do i = 1, n_tree
!!  call set_fcoeff_nba_br (gen, tree_array(i) % trunk, fcoeff)
!!  !--set the unresolved branch level coefficients too!!!
!!  call set_fcoeff_nba_br (gen+1, tree_array(i) % trunk, fcoeff)
!!end do
!  
!if (VERBOSE) call exit_sub (sub_name)
!
!end subroutine fcoeff_nba

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--this code is essentially the same as fcoeff_nba, only difference is
!!  the specific fcoeff_dls_br routine is called--perhaps unify the two
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine fcoeff_dls (gen, num, den)
!implicit none
!
!integer, intent (in) :: gen
!real (rp), intent (out) :: num(nfcoeff, nzone), den(nfcoeff, nzone)
!
!character (*), parameter :: sub_name = mod_name // '.fcoeff_dls'
!
!real (rp), parameter :: eps = epsilon (eps)
!
!character (32) :: fmt
!
!integer :: i
!integer :: nbr(nzone)  !--counter for number of branches in each zone
!
!!---------------------------------------------------------------------
!
!if (VERBOSE) call enter_sub (sub_name)
!
!!--only used in the averaging process
!nbr = 0
!num = 0._rp
!den = 0._rp
!
!do i = 1, n_tree
!
!  if (gen <= tree_array(i) % n_gen) then
! 
!    call fcoeff_dls_br (tree_array(i) % trunk, gen, nbr, num, den)
!    
!  else
!  
!    call error (sub_name, 'gen > n_gen for tree i=', i)
!    
!  end if
!  
!end do
!
!!--perform check to make sure there was a branch in each zone
!if (any (nbr <= 0)) then  !--< just to be safe
!
!  write (fmt, '(a,i0,a)') '(a,i0,a,', nzone, '(i0,1x))'
!  write (msg, fmt) 'zero branches in at least one zone' // n_l //  &
!                   'branch generation =', gen, n_l // 'nbr =', nbr
!                   
!  call error (sub_name, msg)
!
!end if
!
!!--calculate zonal averages
!do i = 1, nfcoeff
!  num(i, :) = num(i, :) / nbr(:)
!  den(i, :) = den(i, :) / nbr(:)
!end do
!
!!--calculate average coefficients
!call fcoeff_divide (num, den)
!  
!if (VERBOSE) call exit_sub (sub_name)
!
!end subroutine fcoeff_dls

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine set_fcoeff_nba_br (gen, br, fcoeff_in)
implicit none

integer, intent (in) :: gen
type (branch_type), intent (inout) :: br
real (rp), intent (in) :: fcoeff_in(nfcoeff, nzone)
                           !--storage for is n, b, a in first dimension

character (*), parameter :: sub_name = mod_name // '.set_fcoeff_nba_br'

integer :: i

!---------------------------------------------------------------------

if (br % gen == gen) then

  br % fcoeff(:) = fcoeff_in(:, br % zone)

else if (br % gen < gen) then

  if (associated (br % sub_branch)) then

    do i = 1, br % n_sub_branch
      call set_fcoeff_nba_br (gen, br % sub_branch(i), fcoeff_in)
    end do
  
  else
  
    !--this will also catch error when gen > ngen
    call error (sub_name, 'expecting br % sub_branch to associated' //  &
                          n_l // 'try checking gen < n_gen')

  end if

else

  call error (sub_name, 'unexpected br % gen > gen')

end if

end subroutine set_fcoeff_nba_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vel_d (vin, vmag, v)
implicit none

real (rp), intent (in) :: vin(nd)
real (rp), intent (out) :: vmag(nfcoeff)
real (rp), intent (out) :: v(nd, nfcoeff)

character (*), parameter :: sub_name = mod_name // '.vel_d'

!---------------------------------------------------------------------

if (nfcoeff /= 1) call error (sub_name, 'expecting nfcoeff = 1')

!--does not use vel_opt for now
!--no need to use fcoeff_dir for the original form of this model
vmag(1) = mag (vin)
v(:, 1) = vin

end subroutine vel_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine vel_d_br (br, vmag, v)
!implicit none
!
!type (branch_type), intent (in) :: br
!real (rp), intent (out) :: vmag(nfcoeff)
!real (rp), intent (out) :: v(nd, nfcoeff)
!
!character (*), parameter :: sub_name = mod_name // '.vel_d_br'
!
!!---------------------------------------------------------------------
!
!if (nfcoeff /= 1) call error (sub_name, 'expecting nfcoeff = 1')
!
!!--does not use vel_opt for now
!!--no need to use fcoeff_dir for the original form of this model
!vmag(1) = mag (br % velscale)
!v(:, 1) = br % velscale
!
!end subroutine vel_d_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--storage is n, b, a in nfcoeff dimensions
!!--expects nfcoeff = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine vel_nba_br (br, vmag, v)
!implicit none
!
!type (branch_type), intent (in) :: br
!real (rp), intent (out) :: vmag(nfcoeff)
!real (rp), intent (out) :: v(nd, nfcoeff)
!
!character (*), parameter :: sub_name = mod_name // '.vel_nba_br'
!
!!---------------------------------------------------------------------
!
!if (nfcoeff /= 3) call error (sub_name, 'expecting nfcoeff = 3')
!
!select case (vel_opt)
!case ('absolute')
!
!  vmag(1) = mag (br % velscale)
!  vmag(2) = mag (br % velscale)
!  vmag(3) = mag (br % velscale)
!    
!  v(:, 1) = vmag(1) * (br % fcoeff_dir(:, 1))
!  v(:, 2) = vmag(2) * (br % fcoeff_dir(:, 2))
!  v(:, 3) = vmag(3) * (br % fcoeff_dir(:, 3))
!
!case ('component')
!
!  v(:, 1) = dot_product (br % velscale, br % fcoeff_dir(:, 1)) *  &
!            (br % fcoeff_dir(:, 1))
!  v(:, 2) = dot_product (br % velscale, br % fcoeff_dir(:, 2)) *  &
!            (br % fcoeff_dir(:, 2))
!  v(:, 3) = dot_product (br % velscale, br % fcoeff_dir(:, 3)) *  &
!            (br % fcoeff_dir(:, 3))
! 
!  vmag(1) = mag (v(:, 1))
!  vmag(2) = mag (v(:, 2))
!  vmag(3) = mag (v(:, 3))
!    
!case ('mixed')
! 
!  v(:, 1) = dot_product (br % velscale, br % fcoeff_dir(:, 1)) *  &
!            (br % fcoeff_dir(:, 1))
!  v(:, 2) = dot_product (br % velscale, br % fcoeff_dir(:, 2)) *  &
!            (br % fcoeff_dir(:, 2))
!  v(:, 3) = dot_product (br % velscale, br % fcoeff_dir(:, 3)) *  &
!            (br % fcoeff_dir(:, 3))
!
!  vmag(1) = mag (br % velscale)
!  vmag(2) = mag (br % velscale)
!  vmag(3) = mag (br % velscale)
!
!case default
!  
!  call error (sub_name, 'invalid vel_opt=' // vel_opt)
!    
!end select
!
!end subroutine vel_nba_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--for now, the nba and dls models are so close that they call the
!  the same routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vel_nba (vin, dir, vmag, v)
implicit none

real (rp), intent (in) :: vin(nd)
real (rp), intent (in) :: dir(nd, nfcoeff)
real (rp), intent (out) :: vmag(nfcoeff)
real (rp), intent (out) :: v(nd, nfcoeff)

character (*), parameter :: sub_name = mod_name // '.vel_nba'

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

call vel_3fcoeff( vin, dir, vmag, v )

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine vel_nba

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--for now, the nba and dls models are so close that they call the
!  the same routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vel_dls (vin, dir, vmag, v)
implicit none

real (rp), intent (in) :: vin(nd)
real (rp), intent (in) :: dir(nd, nfcoeff)
real (rp), intent (out) :: vmag(nfcoeff)
real (rp), intent (out) :: v(nd, nfcoeff)

character (*), parameter :: sub_name = mod_name // '.vel_dls'

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif
call vel_3fcoeff( vin, dir, vmag, v )

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine vel_dls

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--expects nfcoeff = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vel_3fcoeff (vin, dir, vmag, v)
implicit none

real (rp), intent (in) :: vin(nd)
real (rp), intent (in) :: dir(nd, nfcoeff)
real (rp), intent (out) :: vmag(nfcoeff)
real (rp), intent (out) :: v(nd, nfcoeff)

character (*), parameter :: sub_name = mod_name // '.vel_3fcoeff'

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (nfcoeff /= 3) call error (sub_name, 'expecting nfcoeff = 3')

select case (vel_opt)
case ('absolute')

  do i = 1, nfcoeff
    vmag(i) = mag (vin)
    v(:, i) = vmag(i) * dir(:, i)
  end do

case ('component')

  do i = 1, nfcoeff
    v(:, i) = dot_product ( vin, dir(:, i) ) * dir(:, i)
    vmag(i) = mag (v(:, i))
  end do
    
case ('mixed')

  do i = 1, nfcoeff
    v(:, i) = dot_product ( vin, dir(:, i) ) * dir(:, i)
    vmag(i) = mag (vin)
  end do

case default
  
  call error (sub_name, 'invalid vel_opt=' // vel_opt)
    
end select

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine vel_3fcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--for now, this also evaluates areas to go with these coefficients
!--assumes gen <= n_gen
!--nbr, num, den must be zeroed before first call
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine fcoeff_br (br, gen, nbr, num, den)
implicit none

type (branch_type), intent (inout) :: br
integer, intent (in) :: gen
integer, intent (inout) :: nbr(nzone)
real (rp), intent (inout) :: num(nfcoeff, nzone), den(nfcoeff, nzone)
                             !--storage for is n, b, a in first
                             !  dimension

character (*), parameter :: sub_name = mod_name // '.fcoeff_br'

real (rp), parameter :: eps = epsilon (1._rp)

integer :: i

real (rp) :: v(nd, nfcoeff)
real (rp) :: vmag(nfcoeff)
real (rp) :: num_br(nfcoeff), den_br(nfcoeff)

!---------------------------------------------------------------------

if (br % gen == gen) then

  select case ( fmodel )
  case ( 'd' )
      call def_d_dir (br, br % velscale, br % fcoeff_dir)
      call def_d_area (br)  !--define area
      call vel_d (br % velscale, vmag, v)  !--no directions required here
  case ( 'd_germano')
      call error (sub_name, 'not expecting d_germano')
  case ( 'dls' )
      call def_dls_dir (br, br % velscale, br % fcoeff_dir)
      call def_dls_area (br)  !--define area
      call vel_dls (br % velscale, br % fcoeff_dir, vmag, v)
  case ( 'nba' )
      call def_nba_dir (br, br % velscale, br % fcoeff_dir)
      call def_nba_area (br)  !--define area
      call vel_nba (br % velscale, br % fcoeff_dir, vmag, v)
  case default
      call error (sub_name, 'invalid fmodel')
  end select

  do i = 1, nfcoeff
    !--numerator of coeff expression for this branch only
    !--the 2.0 is from the 1/2 in the force model (take it to be in numerator) 
    num_br(i) = -2.0_rp * dot_product (br % ftot, v(:, i)) *  &
                vmag(i) * (br % A(i))
    !--denominator
    den_br(i) = mag (v(:, i))**2 * vmag(i)**2 * (br % A(i))**2
  end do

  !--add this branches contribution to the global average
  !--watch the storage order: n, b, a along first dimension
  num(:, br % zone) = num(:, br % zone) + num_br(:)
  den(:, br % zone) = den(:, br % zone) + den_br(:)

  nbr(br % zone) = nbr(br % zone) + 1  !--branch counter

else if (br % gen < gen) then

  if (associated (br % sub_branch)) then

    do i = 1, br % n_sub_branch
      call fcoeff_br (br % sub_branch(i), gen, nbr, num, den)
    end do
  
  else
  
    !--this will also catch error when gen > ngen
    call error (sub_name, 'expecting br % sub_branch to associated' //  &
                          n_l // 'try checking gen < n_gen')

  end if

else

  call error (sub_name, 'unexpected br % gen > gen')

end if

end subroutine fcoeff_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--calculates the projected area of bbox as seen from flow_dir
!--br % height_bbox, br % width_bbox must already be set
!--br % x_hat, br % y_hat, br % z_hat must already be set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine def_d_area (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.def_d_area'

!---------------------------------------------------------------------

br % A = (br % width_bbox) * (br % height_bbox) *        &
         ( abs (dot_product (flow_dir, br % x_hat)) +    &
           abs (dot_product (flow_dir, br % y_hat)) ) +  &
         (br % width_bbox) * (br % width_bbox) *         &
         abs (dot_product (flow_dir, br % z_hat))

end subroutine def_d_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--define normal, axial, binormal areas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine def_nba_area (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.def_nba_area'

logical, parameter :: orient_indep = .true.  !--perhaps make module variable

real (rp) :: dir(nd)
real (rp) :: area

!---------------------------------------------------------------------

if (orient_indep) then  !--orientation independent: same all components

  !--for now, use volume^{2/3} for area scale
  br % A = ( (br % width_bbox)**2 * (br % height_bbox) )**(2._rp / 3._rp)

else  !--use projected bbox areas

  dir = br % fcoeff_dir(:, 1)
  call abbox_from_dir (dir, br, area)
  br % A(1) = area  !--normal

  dir = br % fcoeff_dir(:, 2)
  call abbox_from_dir (dir, br, area)
  br % A(2) = area  !--binormal
  
  dir = br % fcoeff_dir(:, 3)
  call abbox_from_dir (dir, br, area)
  br % A(3) = area  !--axial

end if

end subroutine def_nba_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--define drag, lift, side areas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine def_dls_area (br)
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.def_dls_area'

logical, parameter :: orient_indep = .true.  !--perhaps make module variable

real (rp) :: dir(nd)
real (rp) :: area

!---------------------------------------------------------------------

if (orient_indep) then  !--orientation independent: same all components

  !--for now, use volume^{2/3} for area scale
  br % A(:) = ( (br % width_bbox)**2 * (br % height_bbox) )**(2._rp / 3._rp)

else  !--use projected bbox areas

  dir = br % fcoeff_dir(:, 1)
  call abbox_from_dir (dir, br, area)
  br % A(1) = area  !--drag

  dir = br % fcoeff_dir(:, 2)
  call abbox_from_dir (dir, br, area)
  br % A(2) = area  !--lift

  dir = br % fcoeff_dir(:, 3)
  call abbox_from_dir (dir, br, area)
  br % A(3) = area  !--side

end if

end subroutine def_dls_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--perhaps move this into trees_setup module?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine abbox_from_dir (dir, br, area)
implicit none

real (rp), intent (in) :: dir(nd)
type (branch_type), intent (inout) :: br
real (rp), intent (out) :: area

!---------------------------------------------------------------------

area = (br % width_bbox) * (br % height_bbox) *   &
       ( abs (dot_product (dir, br % x_hat)) +    &
         abs (dot_product (dir, br % y_hat)) ) +  &
       (br % width_bbox) * (br % width_bbox) *    &
       abs (dot_product (dir, br % z_hat))
       
end subroutine abbox_from_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine def_d_dir (br, v, dir)
implicit none

type (branch_type), intent (in) :: br
real (rp), intent (in) :: v(nd)
real (rp), intent (out) :: dir(nd, nfcoeff)

character (*), parameter :: sn = mod_name // '.def_d_dir'

real (rp), parameter :: eps = epsilon (1._rp)

real (rp) :: d_hat(nd) 

!---------------------------------------------------------------------

if ( nfcoeff /= 1 ) call error ( sn, 'expecting nfcoeff = 1' )

if (mag (v) > eps) then

  d_hat = v / mag (v)

else
  !--just use branch-local x_hat, make sure that fcoeff is zero in this case

  d_hat = br % x_hat
  
  write (msg, '(a,i0)') 'velscale too small to form basis' // n_l //  &
                        'br %ident = ', br % ident
  call warn (sn, msg)
  
end if

dir(:, 1) = d_hat

end subroutine def_d_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this always fills fcoeff_dir with an orthogonal basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine def_nba_dir (br, v, dir)
implicit none

type (branch_type), intent (in) :: br
real (rp), intent (in) :: v(nd)
real (rp), intent (out) :: dir(nd, nfcoeff)

character (*), parameter :: sn = mod_name // '.def_nba_dir'

real (rp), parameter :: eps = epsilon (1._rp)

real (rp) :: n_hat(nd), b_hat(nd), a_hat(nd)

!---------------------------------------------------------------------

if ( nfcoeff /= 3 ) call error ( sn, 'expecting nfcoeff = 3' )

if (mag (v) > eps) then

  a_hat = br % abs_dir

  n_hat = v - ( dot_product (v, a_hat) * (a_hat) )

  !--possible problem when velscale = +/- a_hat
  if (mag (n_hat) > eps) then
    n_hat = (n_hat) / mag (n_hat)
  else  !--use branch-local x_hat
    n_hat = br % x_hat  !--this /is/ orthogonal to a_hat
  end if

  b_hat = cross_product (a_hat, n_hat)

else  !--just use branch-local (x,y,z)=(n,b,a) system

  n_hat = br % x_hat
  b_hat = br % y_hat
  a_hat = br % z_hat
  
  write (msg, '(a,i0)') 'velscale too small to form basis' // n_l //  &
                        'br %ident = ', br % ident
  call warn (sn, msg)
  
end if

dir(:, 1) = n_hat
dir(:, 2) = b_hat
dir(:, 3) = a_hat

end subroutine def_nba_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine def_dls_dir (br, v, dir)
implicit none

type (branch_type), intent (in) :: br
real (rp), intent (in) :: v(nd)
real (rp), intent (out) :: dir(nd, nfcoeff)

character (*), parameter :: sn = mod_name // '.def_dls_dir'

real (rp), parameter :: eps = epsilon (1._rp)

real (rp) :: d_hat(nd), l_hat(nd), s_hat(nd)

!---------------------------------------------------------------------

if ( nfcoeff /= 3 ) call error ( sn, 'expecting nfcoeff = 3' )

if (mag (v) > eps) then

  d_hat = v / mag (v)

else

  write (msg, '(a,i0,a)') 'br % ident = ', br % ident,     &
                ': velscale too close to zero to form basis'
  call error (sn, msg)
      !--figure out way to handle without error
  
end if

l_hat = br % abs_dir - dot_product (br % abs_dir, d_hat) * (d_hat)

if (mag (l_hat) > eps) then  !--normalize

  l_hat = (l_hat) / mag (l_hat)
  
else

  write (msg, '(a,i0,a)') 'br % ident = ', br % ident,     &
                ': l_hat too close to zero to form basis'
  call error (sn, msg)

end if

s_hat = cross_product (d_hat, l_hat)

dir(:, 1) = d_hat
dir(:, 2) = l_hat
dir(:, 3) = s_hat

end subroutine def_dls_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--would like to be able to call this from separate apriori program
!!--nfcoeff may vary between calls?
!!--this puts force coefficients in fcoeff (module-variable)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine fcoeff_d (gen, num, den)
!implicit none
!
!integer, intent (in) :: gen
!real (rp), intent (out) :: num(nfcoeff, nzone), den(nfcoeff, nzone)
!
!character (*), parameter :: sub_name = mod_name // '.fcoeff_d'
!
!!logical, parameter :: DEBUG = .true.
!
!real (rp), parameter :: eps = epsilon (1._rp)
!
!integer :: i
!integer :: navg(nzone)
!
!!type (branch_type), pointer :: b => NULL ()
!
!!---------------------------------------------------------------------
!
!if (VERBOSE) call enter_sub (sub_name)
!
!!--check that nfcoeff = 1
!if (nfcoeff /= 1) call error (sub_name, 'nfcoeff = 1 is required for ' //  &
!                              'fcoeff_d; nfcoeff =', nfcoeff)
!
!num = 0._rp
!den = 0._rp
!navg = 0
!
!do i = 1, n_tree
!
!  if ((gen < 0) .or. (gen > tree_array(i) % n_gen)) then
!    call error (sub_name, 'gen out of range')
!  end if
!  
!  call fcoeff_d_br (tree_array(i) % trunk, gen, navg, num, den)
!
!end do
!
!!--normalize the sums num, den to make them averages
!if (all (navg /= 0)) then
!
!  do i = 1, nzone
!    num(:, i) = num(:, i) / navg(i)
!    den(:, i) = den(:, i) / navg(i)
!  end do
!  
!else
!
!  call error (sub_name, 'navg has at least one zero element')
!
!end if
!
!call fcoeff_divide (num, den)
!
!if (DEBUG) then
!  call mesg (sub_name, 'num =', pack (num, mask=.true.))
!  call mesg (sub_name, 'den =', pack (den, mask=.true.))
!  call mesg (sub_name, 'fcoeff =', pack (fcoeff, mask=.true.))
!end if
!
!if (VERBOSE) call exit_sub (sub_name)
!
!end subroutine fcoeff_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--nfcoeff should be 1 here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!recursive subroutine fcoeff_d_br (br, gen, nbr, num, den)
!implicit none
!
!type (branch_type), intent (inout) :: br
!
!integer, intent (in) :: gen
!integer, intent (inout) :: nbr(nzone)
!real (rp), intent (inout) :: num(nfcoeff, nzone), den(nfcoeff, nzone)
!                             !--may be easier just to use assumed shape
!
!character (*), parameter :: sub_name = mod_name // '.fcoeff_d_br'
!
!!logical, parameter :: DEBUG = .true.
!
!integer :: i
!
!real (rp) :: num_br, den_br
!
!!---------------------------------------------------------------------
!
!if (nfcoeff /= 1) call error (sub_name, 'expecting nfcoeff = 1') 
!
!if (br % gen == gen) then
!
!  call def_d_dir (br)  !--does not really do anything, but here
!                       !  to make this parallel with 3 fcoeff models
!  call def_d_area (br)  !--really only need this on first call
!
!  call vel_d ( br % velscale, vmag, v )  !--no direction needed here
!
!  num_br = -2.0_rp * dot_product ( br % ftot, v(:, 1)) * vmag(1) * (br % A(1))
!  den_br = mag ( v(:, 1) )**2 * vmag(1)**2 * (br % A(1))**2
!
!  !--add this branches contribution to the global average
!  num(1, br % zone) = num(1, br % zone) + num_br
!  den(1, br % zone) = den(1, br % zone) + den_br
!
!  nbr(br % zone) = nbr(br % zone) + 1
!
!  if (DEBUG) then
!
!    call mesg (sub_name, 'br % ident =', br % ident)
!    call mesg (sub_name, 'br % ftot =', br % ftot)
!    call mesg (sub_name, 'br % resf =', br % resf)
!    call mesg (sub_name, 'br % unresf =', br % unresf)
!    call mesg (sub_name, 'br % velscale =', br % velscale)
!    call mesg (sub_name, 'br % A =', br % A)
!    
!  end if
!  
!else
!
!  if (associated (br % sub_branch)) then
!  
!    do i = 1, br % n_sub_branch
!      call fcoeff_d_br (br % sub_branch(i), gen, nbr, num, den)
!    end do
!
!  end if
!
!end if
!
!end subroutine fcoeff_d_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_fmodel_ls
