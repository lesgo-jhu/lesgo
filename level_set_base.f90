!**********************************************************************
module level_set_base
!**********************************************************************
use types, rp => rprec
use param2, only : ld, ny, nz
implicit none

save
public

real (rp), allocatable, dimension(:,:,:) :: phi

contains
!**********************************************************************
subroutine alloc_level_set_base()
!**********************************************************************
implicit none

!--this definition of lbz must be consistent with that in level_set module
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

allocate(phi(ld,ny,$lbz:nz))

return
end subroutine alloc_level_set_base

end module level_set_base
