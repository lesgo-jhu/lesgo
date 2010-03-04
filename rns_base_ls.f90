!**********************************************************************
module rns_base_ls
!**********************************************************************
use types, only : rprec
implicit none

save
public

type rns
  integer :: ntrees, nplanes
end type rns

type rns_planes
  integer :: pindex
  real(rprec) :: bpoints(3,3)
end type rns_planes

type(rns) :: rns_t
type(rns_planes), allocatable, dimension(:) :: rns_planes_t

end module rns_base_ls
