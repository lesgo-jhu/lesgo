module trees_post_mod_ls
use trees_base_ls
use messages
implicit none

save
private

public :: mean_ftot, write_drag_lift, write_CD  !--procedures
public :: itime, first_time  !--vars

character (*), parameter :: mod_name = 'trees_post_mod_ls'

real (rp), parameter :: BOGUS = -1234567890._rp

integer :: itime

logical :: first_time  !--write_drag_lift, keeps track of whether its
                       !  the first call to this routine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_CD (ofile, whichCD)
implicit none

character (*), intent (in) :: ofile
character (*), intent (in) :: whichCD

character (*), parameter :: sub_name = mod_name // '.write_CD'

integer, parameter :: lun = 1

integer :: i

logical :: opn, ext

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

inquire (unit=lun, opened=opn, exist=ext)
if (opn .or. (.not.ext)) then
  call error (sub_name, 'problem with lun=', lun)
end if

inquire (file=ofile, opened=opn)
if (opn) call error (sub_name, 'not expecting file'  // ofile // 'to be open')

open (lun, file=ofile)

do i = 1, n_tree
  call write_CD_br (tree_array(i)%trunk, lun, whichCD)
end do

close (lun)

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine write_CD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is already open
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine write_CD_br (br, lun, whichCD)
implicit none

type (branch_type), intent (in) :: br
integer, intent (in) :: lun
character (*), intent (in) :: whichCD  !--'resolved', 'total'

character (*), parameter :: sub_name = mod_name // '.write_CD_br'

integer :: i

logical :: opn

real (rp) :: CD
real (rp) :: A
real (rp) :: F(nd)
real (rp) :: num, den

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (br % resolved) then
  
  !A = (br % d) * ((br % l) + 0.5_rp * d)  !--frontal area
                              !--add cap contrib here
  A = br % A(1)  !--simulation A, not frontal
  
  select case (whichCD)
  case ('total')
    F = br % ftot
  case ('resolved')
    F = br % resf
  case default
    call error (sub_name, 'invalid whichCD = ' // trim(whichCD))
  end select

  num = -dot_product (F, br % velscale) / mag (br % velscale)
        !--(-) here since want force on branch, not the fluid

  den = 0.5_rp * (mag (br % velscale))**2 * A
  CD = num / den

  inquire (unit=lun, opened=opn)
  if (.not. opn) call error (sub_name, 'expecting open lun=', lun)
  
  write (lun, '(2(i0,1x),4(es13.6,1x))') br % ident, br % gen, CD,      &
                                         mag (br % velscale)**2, num, den

  !--recursion
  if (associated (br % sub_branch)) then

    do i = 1, br % n_sub_branch
      call write_CD_br (br % sub_branch(i), lun, whichCD)
    end do
  
  end if
  
end if

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine write_CD_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_drag_lift (gen)
implicit none

integer, intent (in) :: gen

character (*), parameter :: sub_name = mod_name // '.drag_lift'

logical :: exst, opn

integer :: i

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

!--check 0 <= gen <= ngen
if ((gen < 0) .or. (gen > tree_array (n_tree) % n_gen)) then
  call error (sub_name, 'invalid gen (out of range)')
end if

do i = 1, n_tree
  call write_drag_lift_br (gen, tree_array(i) % trunk)
end do

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes 0 <= gen <= ngen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine write_drag_lift_br (gen, br)
implicit none

integer, intent (in) :: gen
type (branch_type), intent (in) :: br

character (*), parameter :: sub_name = mod_name // '.write_drag_lift_br'

character (64) :: fname
character (64) :: fmt
character (8) :: pos

integer :: i
integer :: lun = 1

real (rp) :: drag, lift, side

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

if (br % gen == gen) then

  !--decompose ftot into drag, lift
  call decompose_drag_lift (br, drag, lift, side)  

  write (fname, '(2(a,i0),a)') 'drag_lift_gen', gen,      &
                               '_ident', br % ident, '.dat'

  !--first time, want to rewind, after want to append
  if (first_time) then
    pos = 'rewind'
  else
    pos = 'append'
  end if

  open (lun, file=fname, action='write', position=pos)
  
  fmt = '(1(i0,1x),3(es13.6,1x))'
  write (lun, fmt) itime, drag, lift, side

  close (lun)

  !--no recursion once we reach the desired generation
  
else if (br % gen < gen) then

  if (associated (br % sub_branch)) then

    do i = 1, br % n_sub_branch
      call write_drag_lift_br (gen, br % sub_branch(i))
    end do
    
  end if

else  !--br % gen > gen

  call error (sub_name, 'unexpected condition: br % gen > gen')
  
end if

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine write_drag_lift_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decompose_drag_lift (br, drag, lift, side)
implicit none

type (branch_type), intent (in) :: br
real (rp), intent (out) :: drag, lift, side

character (*), parameter :: sub_name = mod_name // '.decompose_drag_lift'

real (rp), parameter :: eps = epsilon (0._rp)

logical :: degen  !--will indicate if there is a problem with decomposition

real (rp) :: d_hat(nd), l_hat(nd), y_hat(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

degen = .false.

d_hat = (br % velscale)
if (mag (d_hat) >= eps) then
  d_hat = d_hat / mag (d_hat)
else
  d_hat = flow_dir
  degen = .true.
end if

l_hat = (br % abs_dir) - dot_product ( br % abs_dir, d_hat) * d_hat
if (mag (l_hat) >= eps) then
  l_hat = l_hat / mag (l_hat)
else
  l_hat = (/ 0._rp, 0._rp, 1._rp /)
  degen = .true.
end if

y_hat = cross_product (d_hat, l_hat)

if (.not. degen) then
  drag = dot_product (br % ftot, d_hat)
  lift = dot_product (br % ftot, l_hat)
  side = dot_product (br % ftot, y_hat)
else
  drag = BOGUS
  lift = BOGUS
  side = BOGUS
end if

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine decompose_drag_lift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--would like to be able to call this from separate program
!  this is why nzo is used instead of nzone (may want to vary)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mean_ftot (gen, nzo, ftot, navg)
implicit none

integer, intent (in) :: gen, nzo
real (rp), intent (out) :: ftot(nd, nzo)
integer, intent (out) :: navg(nzo)

character (*), parameter :: sub_name = mod_name // '.mean_ftot'

!logical, parameter :: DEBUG = .true.

!integer :: navg(nzo)
integer :: z

type (branch_type), pointer :: b => NULL ()

!---------------------------------------------------------------------
$if ($VERBOSE)
if (VERBOSE) call enter_sub (sub_name)
$endif

ftot = 0._rp
navg = 0

call sum_ftot_ta (gen, nzo, navg, ftot)

!--normalize the sums num, den to make them averages
do z = 1, nzo

  if (navg(z) /= 0) then
    ftot(1:nd, z) = ftot(1:nd, z) / navg(z)
  else
    ftot(1:nd, z) = BOGUS
  end if
  
end do

$if ($DEBUG)
if (DEBUG) then
  do z = 1, nzo
    call mesg (sub_name, 'ftot =', ftot(:, z))
  end do
end if
$endif

$if ($VERBOSE)
if (VERBOSE) call exit_sub (sub_name)
$endif

end subroutine mean_ftot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sum_ftot_ta (gen, nzo, navg, ftot)
implicit none

integer, intent (in) :: gen, nzo
integer, intent (inout) :: navg(nzo)
real (rp), intent (inout) :: ftot(nd, nzo)

character (*), parameter :: sub_name = mod_name // '.sum_ftot_taa'

integer :: i

!---------------------------------------------------------------------

do i = 1, n_tree

  if ((gen < 0) .or. (gen > tree_array(i) % n_gen)) then
    call error (sub_name, 'gen out of range')
  end if
  
  call sum_ftot (tree_array(i) % trunk, gen, nzo, navg, ftot)

end do

end subroutine sum_ftot_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine sum_ftot (br, gen, nzo, navg, ftot)
implicit none

type (branch_type), intent (in) :: br

integer, intent (in) :: gen, nzo
integer, intent (inout) :: navg(nzo)
real (rp), intent (inout) :: ftot(nd, nzo)

character (*), parameter :: sub_name = mod_name // '.sum_ftot'

!logical, parameter :: DEBUG = .true.

integer :: i

!---------------------------------------------------------------------

if (br % gen == gen) then

  ftot(:, br % zone) = ftot(:, br % zone) + br % ftot

  navg(br % zone) = navg(br % zone) + 1

  $if ($DEBUG)
  if (DEBUG) then

    call mesg (sub_name, 'br % ident =', br % ident)
    call mesg (sub_name, 'br % ftot =', br % ftot)
    call mesg (sub_name, 'br % resf =', br % resf)
    call mesg (sub_name, 'br % unresf =', br % unresf)
    call mesg (sub_name, 'br % velscale =', br % velscale)
    call mesg (sub_name, 'br % A(:) =', br % A(:))
    
  end if
  $endif
  
else

  if (associated (br % sub_branch)) then
  
    do i = 1, br % n_sub_branch
      call sum_ftot (br % sub_branch(i), gen, nzo, navg, ftot)
    end do

  end if

end if

end subroutine sum_ftot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_post_mod_ls
