module level_set_base
use types, only : rp => rprec
use types, only : rprec
use param, only : ld, ny, nz, dx
implicit none

save

private
public :: phi


!--this definition of lbz must be consistent with that in level_set module
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real (rp) :: phi(ld, ny, $lbz:nz)

end module level_set_base
