!!
!!  Copyright (C) 2009-2018  Johns Hopkins University
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

!*******************************************************************************
module data_writer
!*******************************************************************************
use param, only : rprec, write_endian
use messages
#ifdef PPMPI
use mpi
use param, only : nproc, coord, comm, ierr, MPI_RPREC
#endif
#ifdef PPCGNS
use cgns
#endif
implicit none

private
public :: data_writer_t

!  Sums performed over time
type data_writer_t
    integer :: fid
    logical :: opened = .false.
    integer :: nx, ny, nz
    integer :: counter
#ifdef PPMPI
    integer :: nz_tot, nz_below
#ifdef PPCGNS
    integer :: base = 1
    integer :: zone = 1
    integer :: sol = 1
    integer(cgsize_t) :: start_n(3)
    integer(cgsize_t) :: end_n(3)
#endif
#endif
    integer :: num_fields
contains
    procedure, public :: open_file
    procedure, public :: write_field
    procedure, public :: close_file
end type data_writer_t

contains

!*******************************************************************************
subroutine open_file(this, fname, i_nx, i_ny, i_nz, x, y, z, num_fields)
!*******************************************************************************
use mpi
use string_util
class(data_writer_t), intent(inout) :: this
character(*), intent(in) :: fname
! Size of the plane to write
integer, intent(in) :: i_nx, i_ny, i_nz
! Coordinate system
real(rprec), intent(in), dimension(:) :: x, y, z
! Number of fields
integer :: num_fields
! Full file name with extension
character(64) :: full_fname
! Extension
character(64) :: ext
#ifdef PPMPI
integer, dimension(nproc) :: cum_nz, buffer
#ifdef PPCGNS
! Sizes
integer(cgsize_t) :: sizes(3,3)
! Local mesh
real(rprec), dimension(:,:,:), allocatable :: xyz
! Loop through arrays
integer :: i
#endif
#endif

! Set the number of fields to write
this%num_fields = num_fields

! Concatenate extension onto filename
full_fname = fname
#ifdef PPCGNS
ext = '.cgns'
#else
ext = '.bin'
#endif
call string_concat(full_fname, ext)

! Set record counter
this%counter = 1

! Set size of record
this%nx = i_nx
this%ny = i_ny
this%nz = i_nz

#ifdef PPMPI
! Figure out nz_tot and nz_below
cum_nz = 0
cum_nz(coord+1:) = this%nz
call mpi_allreduce(cum_nz, buffer, nproc, MPI_INTEGER, MPI_SUM, comm, ierr)
this%nz_tot = buffer(nproc)
this%nz_below = buffer(coord+1)-this%nz
#endif

write(*,*) coord, this%nz, this%nz_tot, this%nz_below

! Open file if not already open
if (this%opened) call error('data_writer_t%open_file', 'File already opened')

#ifdef PPMPI
#ifdef PPCGNS
! Open CGNS file
call cgp_open_f(full_fname, CG_MODE_WRITE, this%fid, ierr)
if (ierr .ne. CG_OK) call cgp_error_exit_f
#else
! Open MPI file
call mpi_file_open(comm, full_fname, MPI_MODE_WRONLY + MPI_MODE_CREATE,        &
    MPI_INFO_NULL, this%fid, ierr)
#endif
#else
! Open simple fortran file
open(newunit=this%fid, file=full_fname, form='unformatted',                    &
    convert=write_endian, access='direct', recl=this%nx*this%ny*this%nz*rprec)
#endif

#ifdef PPCGNS
    ! Write this%base
    call cg_base_write_f(this%fid, 'Base', 3, 3, this%base, ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f

    ! Sizes, used to create this%zone
    sizes(:,1) = (/int(this%nx, cgsize_t), int(this%ny, cgsize_t),             &
        int(this%nz_tot, cgsize_t)/)
    sizes(:,2) = (/int(this%nx-1, cgsize_t), int(this%ny-1, cgsize_t),         &
        int(this%nz_tot-1, cgsize_t)/)
    sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

    ! Write this%zone
    call cg_zone_write_f(this%fid, this%base, 'Zone', sizes, Structured,       &
        this%zone, ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f

    ! Create data nodes for coordinates
    call cgp_coord_write_f(this%fid, this%base, this%zone, RealDouble,         &
        'CoordinateX', this%nx*this%ny*this%nz, ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f
    call cgp_coord_write_f(this%fid, this%base, this%zone, RealDouble,         &
        'CoordinateY', this%nx*this%ny*this%nz, ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f
    call cgp_coord_write_f(this%fid, this%base, this%zone, RealDouble,         &
        'CoordinateZ', this%nx*this%ny*this%nz,ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f

    ! Set start and end points
    if (this%nz == 0) then
        this%start_n = int(1, cgsize_t)
        this%end_n(1) = int(this%nx, cgsize_t)
        this%end_n(2) = int(this%ny, cgsize_t)
        this%end_n(3) = int(this%nz_tot, cgsize_t)
    else
        this%start_n(1) = int(1, cgsize_t)
        this%start_n(2) = int(1, cgsize_t)
        this%start_n(3) = int(this%nz_below + 1, cgsize_t)
        this%end_n(1) = int(this%nx, cgsize_t)
        this%end_n(2) = int(this%ny, cgsize_t)
        this%end_n(3) = int(this%nz_below + this%nz, cgsize_t)
    end if

    ! Write mesh
    if (this%nz > 0) then

        ! Write local x-mesh
        allocate(xyz(this%nx, this%ny, this%nz))
        do i = 1, this%nx
            xyz(i,:,:) = x(i)
        end do
        call cgp_coord_write_data_f(this%fid, this%base, this%zone, 1,         &
            this%start_n, this%end_n, xyz(:,:,:), ierr)
        if (ierr .ne. CG_OK) call cgp_error_exit_f

        ! Write local y-mesh
        do i = 1, this%ny
            xyz(:,i,:) = y(i)
        end do
        call cgp_coord_write_data_f(this%fid, this%base, this%zone, 2,         &
            this%start_n, this%end_n, xyz(:,:,:), ierr)
        if (ierr .ne. CG_OK) call cgp_error_exit_f

        ! Write local z-mesh
        do i = 1, this%nz
            xyz(:,:,i) = z(i)
        end do
        call cgp_coord_write_data_f(this%fid, this%base, this%zone, 3,         &
            this%start_n, this%end_n, xyz(:,:,:), ierr)
        if (ierr .ne. CG_OK) call cgp_error_exit_f

        ! Create a centered solution
        call cg_sol_write_f(this%fid, this%base, this%zone, 'Solution',        &
            Vertex, this%sol, ierr)
        if (ierr .ne. CG_OK) call cgp_error_exit_f
    end if
#endif

! Mark file as opened
this%opened = .true.

end subroutine open_file

!*******************************************************************************
subroutine write_field(this, field, field_name)
!*******************************************************************************
use mpi
class(data_writer_t), intent(inout) :: this
real(rprec), intent(inout), dimension(:,:,:) :: field
character(*), intent(in) :: field_name
#ifdef PPMPI
#ifdef PPCGNS
integer :: sec
#else
integer(MPI_OFFSET_KIND) :: offset
#endif
#endif

! Check record size
if (size(field, 1) /= this%nx .or. size(field, 2) /= this%ny .or.              &
    size(field, 3) /= this%nz) call error('data_writer_t%write_field',         &
    'Invalid record size')

! Check field counter
if (this%counter > this%num_fields) call error('data_write_t%write_field',     &
    'All records already recorded.')

#ifdef PPMPI
if (this%nz > 0) then
#ifdef PPCGNS
    ! Create the field
    call cgp_field_write_f(this%fid, this%base, this%zone, this%sol,           &
        RealDouble, field_name, sec, ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f

    ! Write data field
    call cgp_field_write_data_f(this%fid, this%base, this%zone, this%sol, sec, &
        this%start_n, this%end_n, field(:,:,:), ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f
#else
    ! Set the offset for each record
    offset = (this%counter-1)*this%nx*this%ny*this%nz_tot*rprec

    call mpi_file_write_at(this%fid,                                           &
        offset+this%nx*this%ny*this%nz_below*rprec, field(:,:,:),              &
        this%nx*this%ny*this%nz, MPI_RPREC, MPI_STATUS_IGNORE, ierr)
#endif
end if
#else
write(this%fid, rec=this%counter) field
#endif

! Increment counter
this%counter = this%counter+1

end subroutine write_field

!*******************************************************************************
subroutine close_file(this)
!*******************************************************************************
class(data_writer_t), intent(inout) :: this

! Close file if open
if (this%opened) then
#ifdef PPMPI
#ifdef PPCGNS
    ! Close the file
    call cgp_close_f(this%fid, ierr)
    if (ierr .ne. CG_OK) call cgp_error_exit_f
#else
    call mpi_file_close(this%fid, ierr)
#endif
#else
    close(this%fid)
#endif
    this%opened = .false.
end if

end subroutine close_file

end module data_writer
