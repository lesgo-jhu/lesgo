!**********************************************************************
module rns_base_ls
!**********************************************************************
use types, only : rprec
implicit none

save
public

!type rns
!  integer :: ntrees, nplanes
!  logical :: plane_u_calc
!end type rns

type ref_plane
  integer :: nzeta, neta ! discretization
  real(rprec), dimension(3) :: p1, p2, p3 !  3 ordered points
  real(rprec) :: u, v, w ! reference values
end type ref_plane

type(ref_plane), pointer, dimension(:) :: ref_plane_t

end module rns_base_ls
