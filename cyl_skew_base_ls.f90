!**********************************************************************
module cylinder_skew_base_ls
!**********************************************************************
use types, only : rprec
use param, only : pi,nproc,nx,ny,nz,nz_tot,L_x,L_y,L_z,dx,dy,dz

implicit none

save
public
!private :: mod_name

character (*), parameter :: mod_name = 'cylinder_skew_base_ls'

real(rprec), parameter :: zrot_angle = -90._rprec*pi/180._rprec
real(rprec), parameter :: skew_angle = 45._rprec*pi/180._rprec

integer, parameter :: ntree = 7

integer, parameter :: ngen = 2
integer, parameter :: ngen_reslv = 2

integer, parameter :: nbranch = 3

real(rprec), parameter :: d = 28.8_rprec*4._rprec/185._rprec
real(rprec), parameter :: l = 50.4_rprec/cos(skew_angle)*4._rprec/185._rprec
real(rprec), parameter :: offset = 9._rprec*4._rprec/185._rprec
real(rprec), parameter :: scale_fact = 0.5_rprec

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
real(rprec), parameter :: z_bottom_surf = 0.125490193552785_rprec

real(rprec), parameter :: filt_width = 2.*dx  !  Filter width for filtered indicator function

integer, dimension(:), allocatable :: igen, kbottom, kbottom_inside, ktop, ktop_inside, lun
integer, dimension(:,:,:), allocatable :: itype
real(rprec), dimension(:), allocatable :: dz_bottom, dz_top

type cs0
     integer :: brindex, iset, itype
     real(rprec) :: phi, chi
     real(rprec), dimension(3) :: xyz
end type cs0

type cs1
    real(rprec), dimension(3) :: xyz
end type cs1

type cs2
  real(rprec), pointer, dimension(:,:) :: xyz
end type cs2

type rot
  real(rprec), pointer, dimension(:) :: angle
  real(rprec), pointer, dimension(:,:) :: axis
end type rot

type vector
  real(rprec) :: mag
  real(rprec), dimension(3) :: xyz
end type vector

type branch
    integer :: indx
    real(rprec) :: d, l, a, b, offset, skew_angle, angle
    real(rprec), dimension(3) :: skew_axis, bot, top
end type branch

type cluster
    integer :: nbranch, indx
    real(rprec), dimension(3) :: origin ! origin of center (at bottom)
    type(branch), pointer, dimension(:) :: br_t
end type cluster

type generation
    integer :: ncluster
    real(rprec) :: bplane, tplane !  assume all branches are the same height
    type(cluster), pointer, dimension(:) :: cl_t
end type generation    

type tree
    real(rprec), dimension(3) :: origin
    integer :: ngen, ngen_reslv
    type(generation), pointer, dimension(:) :: gen_t 
end type tree

integer, pointer, dimension(:,:,:) :: brindex
real(rprec), pointer, dimension(:,:,:) :: phi

type(tree), pointer, dimension(:) :: tr_t ! Tree type
integer, pointer, dimension(:,:) :: clindx_to_loc_id, brindx_to_loc_id

contains 

!**********************************************************************
subroutine set_tree_origin(nt,origin_out)
!**********************************************************************

implicit none

integer, intent(in) :: nt
real(rprec),  dimension(3), intent(out) :: origin_out
real(rprec), allocatable, dimension(:,:) :: origin

allocate(origin(3,ntree))

origin(:,1) = (/ L_x/2., L_y/2., z_bottom_surf /)
origin(:,2) = (/ 0._rprec, L_y, z_bottom_surf /)
origin(:,3) = (/ 0._rprec, 0._rprec, z_bottom_surf /)
origin(:,4) = (/ L_x, 0._rprec, z_bottom_surf /)
origin(:,5) = (/ L_x, L_y, z_bottom_surf /)
origin(:,6) = (/ L_x/2, 3./2.*L_y, z_bottom_surf /)
origin(:,7) = (/ L_x/2, -1./2.*L_y, z_bottom_surf /)

origin_out = origin(:,nt)

deallocate(origin)

return

end subroutine set_tree_origin

end module cylinder_skew_base_ls
