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

integer, parameter :: ntree = 1

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
    real(rprec) :: d
    real(rprec) :: l
    real(rprec) :: a, b
    real(rprec) :: offset
    real(rprec) :: skew_angle, angle
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

type(tree), allocatable, dimension(:) :: tr_t ! Tree type
integer, allocatable, dimension(:,:) :: clindx_to_gen_cl, brindx_to_gen_cl_br

!contains 

!!**********************************************************************
!subroutine set_tree_origin()
!!**********************************************************************

!implicit none



!return
!end subroutine set_tree_origin

end module cylinder_skew_base_ls
