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

module sim_param
use types, only : rprec
use param, only : ld, ny, nz, lbz
implicit none

save
public

logical :: sim_param_initialized = .false.

!--still testing allocatable array implementation:
!  ifc 7.1 segfaults
!  ifort 8.1 ok
!  xlf segfaults when in MPI mode 256^3/32 cpu (need to test other combos)
    
real (rprec), dimension (:, :, :), allocatable :: u, v, w
real (rprec), dimension (:, :, :), allocatable :: dudx, dudy, dudz,  &
                                                  dvdx, dvdy, dvdz,  &
                                                  dwdx, dwdy, dwdz,  &
                                                  RHSx, RHSy, RHSz,  &
                                                  RHSx_f, RHSy_f, RHSz_f

real (rprec), dimension (:, :, :), allocatable :: dpdx, dpdy, dpdz

real (rprec), dimension (:, :, :), allocatable :: txx, txy, tyy
real (rprec), dimension (:, :, :), allocatable :: txz, tyz, tzz

real (rprec), target, dimension (:, :, :), allocatable :: p

real (rprec), dimension (:, :, :), allocatable :: divtx, divty, divtz

real (rprec), dimension (:, :, :), allocatable :: fx, fy, fz, &
                                                  fxa, fya, fza

$if ($USE_RNL)
real (rprec), dimension (:, :, :), allocatable :: u_rnl, v_rnl, w_rnl
real (rprec), dimension (:, :, :), allocatable :: dudx_rnl, dudy_rnl, dudz_rnl,  &
                                                  dvdx_rnl, dvdy_rnl, dvdz_rnl,  &
                                                  dwdx_rnl, dwdy_rnl, dwdz_rnl,  &
                                                  RHSx_rnl, RHSy_rnl, RHSz_rnl
real (rprec), dimension (:, :, :), allocatable :: fxml_rnl, fyml_rnl, fzml_rnl
!!jb - clean up and remove these ^ when no longer needed for testing
$endif

contains

!
! This subroutine initilizes all global arrays defined in the sim_param
! module. Here they are allocated and initialized to zero.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sim_param_init ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

allocate ( u(ld, ny, lbz:nz) ); u = 0.0_rprec
allocate ( v(ld, ny, lbz:nz) ); v = 0.0_rprec
allocate ( w(ld, ny, lbz:nz) ); w = 0.0_rprec
allocate( dudx(ld, ny, lbz:nz) ); dudx = 0.0_rprec    !!jb - clean up and remove these when no longer needed for testing
allocate( dudy(ld, ny, lbz:nz) ); dudy = 0.0_rprec
allocate( dudz(ld, ny, lbz:nz) ); dudz = 0.0_rprec
allocate( dvdx(ld, ny, lbz:nz) ); dvdx = 0.0_rprec
allocate( dvdy(ld, ny, lbz:nz) ); dvdy = 0.0_rprec
allocate( dvdz(ld, ny, lbz:nz) ); dvdz = 0.0_rprec
allocate( dwdx(ld, ny, lbz:nz) ); dwdx = 0.0_rprec
allocate( dwdy(ld, ny, lbz:nz) ); dwdy = 0.0_rprec
allocate( dwdz(ld, ny, lbz:nz) ); dwdz = 0.0_rprec
allocate( RHSx(ld, ny, lbz:nz) ); RHSx = 0.0_rprec
allocate( RHSy(ld, ny, lbz:nz) ); RHSy = 0.0_rprec
allocate( RHSz(ld, ny, lbz:nz) ); RHSz = 0.0_rprec
allocate( RHSx_f(ld, ny, lbz:nz) ); RHSx_f = 0.0_rprec
allocate( RHSy_f(ld, ny, lbz:nz) ); RHSy_f = 0.0_rprec
allocate( RHSz_f(ld, ny, lbz:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
allocate ( txx(ld, ny, lbz:nz) ); txx = 0.0_rprec
allocate ( txy(ld, ny, lbz:nz) ); txy = 0.0_rprec
allocate ( tyy(ld, ny, lbz:nz) ); tyy = 0.0_rprec
allocate ( txz(ld, ny, lbz:nz) ); txz = 0.0_rprec
allocate ( tyz(ld, ny, lbz:nz) ); tyz = 0.0_rprec
allocate ( tzz(ld, ny, lbz:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
allocate ( divtx(ld, ny, lbz:nz) ); divtx = 0.0_rprec
allocate ( divty(ld, ny, lbz:nz) ); divty = 0.0_rprec
allocate ( divtz(ld, ny, lbz:nz) ); divtz = 0.0_rprec

$if ($USE_RNL)
allocate ( u_rnl(ld, ny, lbz:nz) ); u_rnl = 0.0_rprec
allocate ( v_rnl(ld, ny, lbz:nz) ); v_rnl = 0.0_rprec
allocate ( w_rnl(ld, ny, lbz:nz) ); w_rnl = 0.0_rprec
!! note dudx_rnl, dvdy_rnl, and dwdz_rnl are not here
!! since they are not needed for convec subroutine
allocate( dudx_rnl(ld, ny, lbz:nz) ); dudx_rnl = 0.0_rprec
allocate( dudy_rnl(ld, ny, lbz:nz) ); dudy_rnl = 0.0_rprec
allocate( dudz_rnl(ld, ny, lbz:nz) ); dudz_rnl = 0.0_rprec
allocate( dvdx_rnl(ld, ny, lbz:nz) ); dvdx_rnl = 0.0_rprec
allocate( dvdz_rnl(ld, ny, lbz:nz) ); dvdz_rnl = 0.0_rprec
allocate( dwdx_rnl(ld, ny, lbz:nz) ); dwdx_rnl = 0.0_rprec
allocate( dwdy_rnl(ld, ny, lbz:nz) ); dwdy_rnl = 0.0_rprec
allocate( RHSx_rnl(ld, ny, lbz:nz) ); RHSx_rnl = 0.0_rprec
allocate( RHSy_rnl(ld, ny, lbz:nz) ); RHSy_rnl = 0.0_rprec
allocate( RHSz_rnl(ld, ny, lbz:nz) ); RHSz_rnl = 0.0_rprec
allocate( fxml_rnl(ld, ny, nz) ); fxml_rnl = 0.0_rprec
allocate( fyml_rnl(ld, ny, nz) ); fyml_rnl = 0.0_rprec
allocate( fzml_rnl(ld, ny, nz) ); fzml_rnl = 0.0_rprec
$endif

$if($TURBINES)
allocate ( fxa(ld, ny, nz) ); fxa = 0.0_rprec
$endif

$if($LVLSET)
allocate ( fx(ld, ny, nz) ); fx = 0.0_rprec
allocate ( fy(ld, ny, nz) ); fy = 0.0_rprec
allocate ( fz(ld, ny, nz) ); fz = 0.0_rprec
! May already be allocated if using TURBINES
if( .not. allocated(fxa) ) allocate ( fxa(ld, ny, nz) ); fxa = 0.0_rprec
allocate ( fya(ld, ny, nz) ); fya = 0.0_rprec
allocate ( fza(ld, ny, nz) ); fza = 0.0_rprec
$endif 

sim_param_initialized = .true.

return
end subroutine sim_param_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module sim_param
