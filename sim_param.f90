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

real (rprec), dimension (:, :, :), allocatable :: theta, q
    !--Added for scalars

$if($BINARY)    
real (rprec), dimension (:, :, :), allocatable :: u_avg
$endif
    
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sim_param_init ( array_list_opt )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
!--this is needed to make some post-processing code compilable, since
!  is only allocates the space to be used
!--initialized all data to zero
!
implicit none

character (*), intent (in), optional :: array_list_opt
    !--comma separated list of arrays to init

integer, parameter :: narray_max = 128
integer, parameter :: narray_name_len = 8  !--can increase if needed

$if ($TURBINES)
character (*), parameter :: array_list_def = 'u, v, w,' //                 &
                                             'dudx, dudy, dudz,' //        &
                                             'dvdx, dvdy, dvdz,' //        &
                                             'dwdx, dwdy, dwdz,' //        &
                                             'RHSx, RHSy, RHSz,' //        &
                                             'RHSx_f, RHSy_f, RHSz_f,' //  &
                                             'dpdx, dpdy, dpdz,' //        &
                                             'txx, txy, tyy,' //           &
                                             'txz, tyz, tzz,' //           &
                                             'p,' //                       &
                                             'divtx, divty, divtz,' //     &
                                             'fxa,'           //           &
                                             'theta, q,' //                &
                                             'u_avg'                                             
$elseif ($LVLSET)
character (*), parameter :: array_list_def = 'u, v, w,' //                 &
                                             'dudx, dudy, dudz,' //        &
                                             'dvdx, dvdy, dvdz,' //        &
                                             'dwdx, dwdy, dwdz,' //        &
                                             'RHSx, RHSy, RHSz,' //        &
                                             'RHSx_f, RHSy_f, RHSz_f,' //  &
                                             'dpdx, dpdy, dpdz,' //        &
                                             'txx, txy, tyy,' //           &
                                             'txz, tyz, tzz,' //           &
                                             'p,' //                       &
                                             'divtx, divty, divtz,' //     &
                                             'fx, fy, fz,' //              &
                                             'fxa, fya, fza,' //           &
                                             'theta, q,' //                &
                                             'u_avg'                                             
$else
character (*), parameter :: array_list_def = 'u, v, w,' //                 &
                                             'dudx, dudy, dudz,' //        &
                                             'dvdx, dvdy, dvdz,' //        &
                                             'dwdx, dwdy, dwdz,' //        &
                                             'RHSx, RHSy, RHSz,' //        &
                                             'RHSx_f, RHSy_f, RHSz_f,' //  &
                                             'dpdx, dpdy, dpdz,' //        &
                                             'txx, txy, tyy,' //           &
                                             'txz, tyz, tzz,' //           &
                                             'p,' //                       &
                                             'divtx, divty, divtz,' //     &
                                             'theta, q,' //                &
                                             'u_avg'                                             
$endif

$if ($DEBUG)
logical, parameter :: DEBUG = .true.
$endif

character (narray_name_len * narray_max) :: array_list
character (narray_name_len) :: array(narray_max)
character (narray_name_len) :: alloced_array(narray_max)

integer :: i
integer :: ios

!---------------------------------------------------------------------

array = ''
alloced_array = ''

if ( present ( array_list_opt ) ) then
    array_list = array_list_opt
else  !--use default list: all arrays allocated
    array_list = array_list_def
end if

read ( array_list, *, iostat=ios ) array
!--comma or space separated

!--ios < 0 indicates array_list is too small  (or too large) to fit
!  into array, which for now is OK (at least when its too small)
if ( ios > 0 ) then
    write (*, *) 'sim_param_init: error reading array list'
    stop
end if

do i = 1, size ( array )
    
    if ( trim (array(i)) == '' ) exit
        !--assume rest of array(i:) is blank as well

    select case ( trim ( array(i) ) )
    case ( 'u' )
        allocate ( u(ld, ny, lbz:nz) )
        u = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'v' )
        allocate ( v(ld, ny, lbz:nz) )
        v = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'w' )
        allocate ( w(ld, ny, lbz:nz) )
        w = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dudx' )
        allocate( dudx(ld, ny, lbz:nz) )
        dudx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dudy' )
        allocate( dudy(ld, ny, lbz:nz) )
        dudy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dudz' )
        allocate( dudz(ld, ny, lbz:nz) )
        dudz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dvdx' )
        allocate( dvdx(ld, ny, lbz:nz) )
        dvdx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dvdy' )
        allocate( dvdy(ld, ny, lbz:nz) )
        dvdy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dvdz' )
        allocate( dvdz(ld, ny, lbz:nz) )
        dvdz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dwdx' )
        allocate( dwdx(ld, ny, lbz:nz) )
        dwdx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dwdy' )
        allocate( dwdy(ld, ny, lbz:nz) )
        dwdy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dwdz' )
        allocate( dwdz(ld, ny, lbz:nz) )
        dwdz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSx' )
        allocate( RHSx(ld, ny, lbz:nz) )
        RHSx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSy' )
        allocate( RHSy(ld, ny, lbz:nz) )
        RHSy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSz' )
        allocate( RHSz(ld, ny, lbz:nz) )
        RHSz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSx_f' )
        allocate( RHSx_f(ld, ny, lbz:nz) )
        RHSx_f = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSy_f' )
        allocate( RHSy_f(ld, ny, lbz:nz) )
        RHSy_f = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSz_f' )
        allocate( RHSz_f(ld, ny, lbz:nz) )
        RHSz_f = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dpdx' )
        allocate ( dpdx(ld, ny, nz) )
        dpdx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dpdy' )
        allocate ( dpdy(ld, ny, nz) )
        dpdy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dpdz' )
        allocate ( dpdz(ld, ny, nz) )
        dpdz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'txx' )
        allocate ( txx(ld, ny, lbz:nz) )
        txx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'txy' )
        allocate ( txy(ld, ny, lbz:nz) )
        txy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'tyy' )
        allocate ( tyy(ld, ny, lbz:nz) )
        tyy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'txz' )
        allocate ( txz(ld, ny, lbz:nz) )
        txz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'tyz' )
        allocate ( tyz(ld, ny, lbz:nz) )
        tyz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'tzz' )
        allocate ( tzz(ld, ny, lbz:nz) )
        tzz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'p' )
        allocate ( p(ld, ny, 0:nz) )
        p = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'divtx' )
        allocate ( divtx(ld, ny, lbz:nz) )
        divtx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'divty' )
        allocate ( divty(ld, ny, lbz:nz) )
        divty = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'divtz' )
        allocate ( divtz(ld, ny, lbz:nz) )
        divtz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fx' )
        allocate ( fx(ld, ny, nz) )        
        fx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fy' )
        allocate ( fy(ld, ny, nz) )        
        fy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )  
    case ( 'fz' )
        allocate ( fz(ld, ny, nz) )        
        fz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fxa' )
        allocate ( fxa(ld, ny, nz) )
        fxa = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fya' )
        allocate ( fya(ld, ny, nz) )
        fya = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )  
    case ( 'fza' )
        allocate ( fza(ld, ny, nz) )
        fza = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'theta' )
        allocate ( theta(ld, ny, nz) )
        theta = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'q' )
        allocate ( q(ld, ny, nz) )
        q = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
$if($BINARY)
    case ( 'u_avg' )
        allocate ( u_avg(ld, ny, nz) )
        u_avg = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
$endif
    case default
        write (*, *) 'sim_param_init: invalid array ' // trim ( array(i) )
        stop
    end select

end do

write(*,*) 'Setting sim_param_initialized true'
sim_param_initialized = .true.
    !--beware this does NOT guarentee that a particular set of arrays
    !  was initialized!

$if ($DEBUG)
if ( DEBUG ) then

    write (*, *) 'sim_param_init: the following arrays were initialized'

    do i = 1, size (alloced_array)
        if ( trim ( alloced_array(i) ) == '' ) cycle
        write ( *, '(3a)' ) 'array=<', trim (alloced_array(i)), '>'
    end do
    
end if
$endif

end subroutine sim_param_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module sim_param
