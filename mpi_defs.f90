!**********************************************************************
module mpi_defs
!**********************************************************************
use mpi
implicit none

save
private

public :: initialize_mpi, mpi_sync_real_array

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

!**********************************************************************
subroutine mpi_sync_real_array(var)
!**********************************************************************
!  This subroutine syncs data for the:
!    k = nz-1 at coord to k=0 at coord+1
!    k = 1 at coord+1  to k=nz at coord
!  nodes; these are the ghost and interprocessor overlap nodes. 
use types, only : rprec
use mpi
use param, only : MPI_RPREC, down, up, comm, status, ierr, nproc, coord

implicit none

real(rprec), dimension(:,:,:), intent(INOUT) :: var

integer :: lbx,ubx,lby,uby,lbz,ubz
integer :: mpi_datasize

!  Get bounds of var array
lbx=lbound(var,1); ubx=ubound(var,1)
lby=lbound(var,2); uby=ubound(var,2)
lbz=lbound(var,3); ubz=ubound(var,3)

!  Set mpi data size
mpi_datasize = (ubx-lbx+1)*(uby-lby+1)

!  ----- Need to get all overlapping values -----
if(coord < nproc - 1) then
  call mpi_send (var(:,:,ubz-1), mpi_datasize, MPI_RPREC, up, 1, comm, ierr)
  call mpi_recv (var(:,:,ubz), mpi_datasize, MPI_RPREC, up, 2, comm, status, ierr)
endif

if(coord > 0) then
  call mpi_recv(var(:,:,lbz), mpi_datasize, MPI_RPREC, down, 1, comm, status, ierr)
  call mpi_send (var(:,:,lbz+1), mpi_datasize, MPI_RPREC, down, 2, comm, ierr)
endif

return

end subroutine mpi_sync_real_array

end module mpi_defs

