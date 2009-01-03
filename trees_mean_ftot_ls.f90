module trees_mean_ftot_ls
use trees_base_ls
use messages
implicit none

save
private

public :: mean_ftot

character (*), parameter :: mod_name = 'trees_mean_ftot_ls'

real (rp), parameter :: BOGUS = -1234567890._rp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
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

if (VERBOSE) call enter_sub (sub_name)

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

if (DEBUG) then
  do z = 1, nzo
    call mesg (sub_name, 'ftot =', ftot(:, z))
  end do
end if

if (VERBOSE) call exit_sub (sub_name)

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

  if (DEBUG) then

    call mesg (sub_name, 'br % ident =', br % ident)
    call mesg (sub_name, 'br % ftot =', br % ftot)
    call mesg (sub_name, 'br % resf =', br % resf)
    call mesg (sub_name, 'br % unresf =', br % unresf)
    call mesg (sub_name, 'br % velscale =', br % velscale)
    call mesg (sub_name, 'br % Ap_bbox =', br % Ap_bbox)
    
  end if
  
else

  if (associated (br % sub_branch)) then
  
    do i = 1, br % n_sub_branch
      call sum_ftot (br % sub_branch(i), gen, nzo, navg, ftot)
    end do

  end if

end if

end subroutine sum_ftot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_mean_ftot_ls
