!**********************************************************************
module rns_base_ls
!**********************************************************************
use types, only : rprec
implicit none

save
public

type rns
  integer :: ntrees, nplanes
  logical :: plane_u_calc
end type rns

type rns_planes
  integer :: indx
  real(rprec) :: bp(3,3)
  real(rprec) :: u, v, w
end type rns_planes

type(rns) :: rns_t
type(rns_planes), allocatable, dimension(:) :: rns_planes_t

end module rns_base_ls
