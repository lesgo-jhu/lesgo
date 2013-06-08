!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!**********************************************************************
module cyl_skew_base_ls
!**********************************************************************
use types, only : rprec, point3D_t
use param, only : pi

implicit none

save

public

private :: rprec, point3D_t
private :: pi
private :: mod_name

character (*), parameter :: mod_name = 'cyl_skew_base_ls'

!---------------------------------------------------
! CYL_SKEW TREE PARAMETERS
!--------------------------------------------------- 
real(rprec) :: zrot_angle = -90._rprec*pi/180._rprec
real(rprec) :: skew_angle = 45._rprec*pi/180._rprec

! Bottom and top surfaces
logical :: use_bottom_surf = .false. !  True for making a bottom surface
real(rprec) :: z_bottom_surf = 0.0_rprec ! Non-dimensional
logical :: use_top_surf = .false.
real(rprec) :: z_top_surf = 0.0_rprec ! Non-dimensional
logical :: use_left_surf = .false.
real(rprec) :: y_left_surf = 0.0_rprec ! Non-dimensional
logical :: use_right_surf = .false.
real(rprec) :: y_right_surf = 0.0_rprec ! Non-dimensional


! Tree settings
integer :: ntree = 1
type(point3D_t), allocatable, dimension(:) :: tree_location 

! Number of generations (total and resolved)
integer :: ngen = 5
integer :: ngen_reslv = 2
! Number of branches of generator
integer :: nbranch = 3

!  Make sure they are non-dimensional
! Generator branch diameter and length
real(rprec) :: d = 1._rprec
real(rprec) :: l = 1._rprec
real(rprec) :: offset = 0.1_rprec

! Fractal scale factor
real(rprec) :: scale_fact = 0.5_rprec

! Compute the filter indicator function (if true)
logical :: filter_chi = .false.
real(rprec) :: filt_width = 0.125  !  Filter width for filtered indicator function

!---------------------------------------------------
!
!---------------------------------------------------

integer, dimension(:), allocatable :: igen, kbottom, kbottom_inside, ktop, ktop_inside
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
