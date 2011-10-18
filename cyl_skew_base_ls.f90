!**********************************************************************
module cyl_skew_base_ls
!**********************************************************************
use types, only : rprec, point3D
use param, only : pi,nproc,nx,ny,nz,nz_tot,L_x,L_y,L_z,dx,dy,dz,z_i

implicit none

save

public

private :: rprec
private :: pi,nproc,nx,ny,nz,nz_tot,L_x,L_y,L_z,dx,dy,dz
private :: mod_name

character (*), parameter :: mod_name = 'cyl_skew_base_ls'

!---------------------------------------------------
! CYL_SKEW TREE PARAMETERS
!--------------------------------------------------- 

real(rprec), parameter :: zrot_angle = -90._rprec*pi/180._rprec
real(rprec), parameter :: skew_angle = 45._rprec*pi/180._rprec

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
real(rprec), parameter :: z_bottom_surf = 0.6_rprec ! Already in non-dimensional units

integer, parameter :: ntree = 1

type(point3D), parameter, dimension(ntree) :: tree_location = (/ &
     point3D( (/ L_x / 2, L_y / 2, z_bottom_surf /) ) &
     /)

integer, parameter :: ngen = 1
integer, parameter :: ngen_reslv = 1

integer, parameter :: nbranch = 1

!  Make sure they are non-dimensional
!  dm = 28.8 mm
!  hm = 50.4 mm
!  offset = 9 mm
real(rprec), parameter :: d = 1._rprec
real(rprec), parameter :: l = 1._rprec
real(rprec), parameter :: offset = 0._rprec

real(rprec), parameter :: scale_fact = 0.5_rprec

logical, parameter :: filter_chi = .false.
real(rprec), parameter :: filt_width = 2._rprec*dx  !  Filter width for filtered indicator function

!---------------------------------------------------
!
!---------------------------------------------------

integer, dimension(:), allocatable :: igen, kbottom, kbottom_inside, ktop, ktop_inside, lun
integer, dimension(:,:,:), allocatable :: itype
real(rprec), dimension(:), allocatable :: dz_bottom, dz_top

type cs0
     integer :: clindx, brindx, iset, itype
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

type point_2d
    real(rprec), dimension(2) :: xy
end type point_2d

type point_3d
    real(rprec), dimension(3) :: xyz
end type point_3d

type branch
    integer :: indx
    real(rprec) :: d, l, a, b, offset, skew_angle, angle
    real(rprec), dimension(3) :: skew_axis, bot, top
end type branch

type cluster
    integer :: nbranch, indx, parent
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
    integer :: ngen, ngen_reslv, ncluster, nbranch
    type(generation), pointer, dimension(:) :: gen_t 
end type tree

type(tree), pointer, dimension(:) :: tr_t ! Tree type
integer, pointer, dimension(:,:) :: clindx_to_loc_id, brindx_to_loc_id
integer, pointer, dimension(:,:) :: reslv_clindx_to_loc_id, unreslv_clindx_to_loc_id

integer :: ncluster_reslv, ncluster_tot

end module cyl_skew_base_ls
