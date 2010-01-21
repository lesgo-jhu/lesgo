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

double precision, parameter :: zrot_angle = -90._rprec*pi/180._rprec
double precision, parameter :: skew_angle = 45._rprec*pi/180._rprec

integer, parameter :: ntree = 7
integer, parameter :: ntrunk = 3
integer, parameter :: ngen = 3
real(rprec), parameter :: d = 28.8_rprec*4._rprec/185._rprec
real(rprec), parameter :: l = 50.4_rprec/cos(skew_angle)*4._rprec/185._rprec
real(rprec), parameter :: offset = 9._rprec*4._rprec/185._rprec
real(rprec), parameter :: scale_fact = 0.5_rprec

logical, parameter :: use_bottom_surf = .true. !  True for making a bottom surface
real(rprec), parameter :: z_bottom_surf = 4.*dz

integer, dimension(:), allocatable :: igen, kbottom, kbottom_inside, ktop, ktop_inside, lun
integer, dimension(:,:,:), allocatable :: itype
real(rprec), dimension(:), allocatable :: dz_bottom, dz_top

type cs0
     integer :: brindex, iset, itype
     real(rprec) :: phi
     real(rprec), dimension(3) :: xyz
end type cs0

type cs1
    real(rprec), dimension(3) :: xyz
end type cs1

type cs2
    !double precision, allocatable, dimension(:,:) :: xyz
  real(rprec), pointer, dimension(:,:) :: xyz
end type cs2

type rot
  real(rprec), pointer, dimension(:) :: angle
  real(rprec), pointer, dimension(:,:) :: axis
end type rot

type vector
  real(rprec), dimension(3) :: xyz
end type vector

integer, pointer, dimension(:,:,:) :: brindex
real(rprec), pointer, dimension(:,:,:) :: phi

end module cylinder_skew_base_ls
