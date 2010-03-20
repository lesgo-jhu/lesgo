!**********************************************************************
module mpi_defs
!**********************************************************************
use mpi
implicit none

save
private

public :: initialize_mpi

contains

!**********************************************************************
subroutine initialize_mpi()
!**********************************************************************
use types, only : rprec
use param
implicit none

integer :: ip, np, coords(1)

!--check for consistent preprocessor & param.f90 definitions of 
!  MPI and $MPI
if (.not. USE_MPI) then
  write (*, *) 'inconsistent use of USE_MPI and $MPI'
  stop
end if

call mpi_init (ierr)
call mpi_comm_size (MPI_COMM_WORLD, np, ierr)
call mpi_comm_rank (MPI_COMM_WORLD, global_rank, ierr)

  !--check if run-time number of processes agrees with nproc parameter
if (np /= nproc) then
  write (*, *) 'runtime number of procs = ', np,  &
               ' not equal to nproc = ', nproc
  stop
end if

  !--set up a 1d cartesian topology 
call mpi_cart_create (MPI_COMM_WORLD, 1, (/ nproc /), (/ .false. /),  &
  .true., comm, ierr)
  !--slight problem here for ghost layers:
  !  u-node info needs to be shifted up to proc w/ rank "up",
  !  w-node info needs to be shifted down to proc w/ rank "down"
call mpi_cart_shift (comm, 0, 1, down, up, ierr)
call mpi_comm_rank (comm, rank, ierr)
call mpi_cart_coords (comm, rank, 1, coords, ierr)
coord = coords(1)  !--use coord (NOT rank) to determine global position

write (chcoord, '(a,i0,a)') '(', coord, ')'  !--() make easier to use

  !--rank->coord and coord->rank conversions
do ip = 0, np-1
  call mpi_cart_rank (comm, (/ ip /), rank_of_coord(ip), ierr)
  call mpi_cart_coords (comm, ip, 1, coord_of_rank(ip), ierr)
end do

write (*, *) 'Hello! from process with coord = ', coord

  !--set the MPI_RPREC variable
if (rprec == kind (1.e0)) then
  MPI_RPREC = MPI_REAL
  MPI_CPREC = MPI_COMPLEX
else if (rprec == kind (1.d0)) then
  MPI_RPREC = MPI_DOUBLE_PRECISION
  MPI_CPREC = MPI_DOUBLE_COMPLEX
else
  write (*, *) 'error defining MPI_RPREC/MPI_CPREC'
  stop
end if

return
end subroutine initialize_mpi

end module mpi_defs

