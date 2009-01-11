module level_set_base
use types, rp => rprec
use param2, only : ld, ny, nz
implicit none

save
public


!--this definition of lbz must be consistent with that in level_set module
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real (rp) :: phi(ld, ny, $lbz:nz)

end module level_set_base
