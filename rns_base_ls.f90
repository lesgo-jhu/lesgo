!!
!!  Copyright (C) 2010-2013  Johns Hopkins University
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
module rns_base_ls
!**********************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz
implicit none

save

public

private :: rprec
private :: ld, ny, nz, lbz

!public r_elem_t, beta_elem_t, b_elem_t

!---------------------------------------------------
! RNS PARAMETERS
!---------------------------------------------------  
! kc-3
!integer, parameter :: rns_ntree = 4 ! Number of unique trees 
! vtree-3
!integer, parameter :: rns_ntree = 3 ! Number of unique trees
integer :: rns_ntree = 1 

!  Options: 'default', 'kc-3'
! kc-3
!character(*), parameter :: rns_tree_layout = 'kc-3'
! vtree-3
!character(*), parameter :: rns_tree_layout = 'default'
character(25) :: rns_tree_layout = 'default'

!  Weighting/averaging: Off - 0, On - 1
!integer, parameter :: temporal_weight = 0
integer :: temporal_weight = 1
  
  ! Time weighting constant
  !real(rprec), parameter :: Tconst = 1.0_rprec
  real(rprec) :: Tconst = 1.0_rprec
  ! Time step to start weighting
  !integer, parameter :: weight_nstart = 20000
  integer :: weight_nstart = 10000

!  Explict - 1, Implicit - 2
!integer, parameter :: temporal_model = 1
integer :: temporal_model = 1

!  Local - 1, Global - 2
!integer, parameter :: spatial_model = 2
integer :: spatial_model = 1

!integer, parameter :: output_nskip = 10
!integer, parameter :: CD_ramp_nstep = 10000
integer :: output_nskip = 10
integer :: CD_ramp_nstep = 5000

!  Parameters for setting reference regions
!real(rprec), parameter :: alpha_width = 2.0_rprec
!real(rprec), parameter :: alpha_dist = 0.571428571_rprec
real(rprec) :: alpha_width = 1.0
real(rprec) :: alpha_dist = 0.5

!real(rprec), parameter :: chi_cutoff = 1.0e-9_rprec
real(rprec) :: chi_cutoff = 1.0e-9

! Number of spatial dimensions (MUST BE 2)
integer, parameter :: ndim = 2

!---------------------------------------------------
! 
!---------------------------------------------------  

integer :: nr_elem, nbeta_elem, nb_elem
!real(rprec) :: chi(ld, ny, lbz:nz)
real(rprec), allocatable, dimension(:,:,:) :: chi
logical :: chi_initialized = .false.

! ---- Secondary structures ----
type ref_region
  integer :: npoint
  real(rprec), pointer, dimension(:,:) :: points !  3 ordered points
  real(rprec) :: u, v ! reference values
  real(rprec) :: area
end type ref_region

type force_type_1
  !integer :: parent !  parent CD; for resolved branches 
  real(rprec) :: fx, fy
  real(rprec) :: CD
  real(rprec) :: kappa ! Used for unresolved branches
end type force_type_1

type force_type_2
  !integer :: parent !  parent CD; for resolved branches 
  real(rprec) :: fx, fy
  real(rprec) :: CD
  real(rprec) :: error, error_norm
  real(rprec) :: LAB, LBB ! Used for temporal weighting
end type force_type_2

type indx_array
  integer :: npoint
  integer, pointer, dimension(:,:) :: iarray
end type indx_array

type child_elem
  integer :: nelem
  integer, pointer, dimension(:) :: indx
end type child_elem

! ---- Secondary structures ----

! ---- Primary structures ----
type primary_struct_type_1
  type(ref_region)   :: ref_region_t
  type(force_type_1) :: force_t
  type(indx_array)   :: indx_array_t
end type primary_struct_type_1

type primary_struct_type_2
  type(ref_region)   :: ref_region_t
  type(force_type_2) :: force_t
  type(child_elem)   :: r_child_t
  type(child_elem)   :: beta_child_t
end type primary_struct_type_2

type(primary_struct_type_1), pointer, dimension(:) :: r_elem_t
type(primary_struct_type_1), pointer, dimension(:) :: beta_elem_t
type(primary_struct_type_2), pointer, dimension(:) :: b_elem_t
! ---- Primary structures ----

!integer, pointer, dimension(:) :: reslv_to_rbeta_map
!integer, pointer, dimension(:) :: beta_to_rbeta_map

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rns_base_init_ls
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : path
implicit none

character(*), parameter :: mkdir_cmd = 'mkdir -p ' // path // 'output/rns'

! Make output directory for RNS module
call system(mkdir_cmd)

! Allocate space for signed distance function
allocate( chi( ld, ny, lbz:nz ) )

return
end subroutine rns_base_init_ls

end module rns_base_ls
