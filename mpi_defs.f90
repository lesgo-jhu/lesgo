!**********************************************************************
module mpi_defs
!**********************************************************************
use mpi
implicit none

save
private

public :: initialize_mpi, mpi_sync_real_array
public :: MPI_SYNC_DOWN, MPI_SYNC_UP, MPI_SYNC_DOWNUP

character (*), parameter :: mod_name = 'mpi_defs'

integer, parameter :: MPI_SYNC_DOWN=1
integer, parameter :: MPI_SYNC_UP=2
integer, parameter :: MPI_SYNC_DOWNUP=3

contains

!**********************************************************************
subroutine initialize_mpi()
!**********************************************************************
use types, only : rprec
use param
$if($CPS)
use concurrent_precursor
$endif
implicit none

integer :: ip, np, coords(1)

$if($CPS)
integer :: world_np, world_rank
integer :: remoteLeader
integer :: memberKey
integer :: localComm
$endif

!--check for consistent preprocessor & param.f90 definitions of 
!  MPI and $MPI
if (.not. USE_MPI) then
  write (*, *) 'inconsistent use of USE_MPI and $MPI'
  stop
end if

call mpi_init (ierr)

$if($CPS)

! Get number of processors in world comm
call mpi_comm_size (MPI_COMM_WORLD, world_np, ierr)
call mpi_comm_rank (MPI_COMM_WORLD, world_rank, ierr)

! Set color and remote leader for intercommunicator interComm
if( world_rank < world_np / 2 ) then
   color = RED
   remoteLeader = world_np / 2
else
   color = BLACK
   remoteLeader = 0
endif

memberKey=modulo(world_rank, world_np / 2)

! Split the world communicator into intracommunicators localComm
call MPI_Comm_split(MPI_COMM_WORLD, color, memberKey, localComm, ierr)
call mpi_comm_size (localComm, np, ierr)
call mpi_comm_rank (localComm, global_rank, ierr)

! Create intercommunicator interComm
call mpi_intercomm_create( localComm, 0, MPI_COMM_WORLD, remoteLeader, 1, interComm, ierr)

!--set up a 1d cartesian topology with localComm creates
call mpi_cart_create (localComm, 1, (/ nproc /), (/ .false. /),  &
  .true., comm, ierr)

!write(*,'(a,4i3)') 'program black, world_np, world_rank, color, remoteLeader : ', world_np, world_rank, color, remoteLeader

$else

call mpi_comm_size (MPI_COMM_WORLD, np, ierr)
call mpi_comm_rank (MPI_COMM_WORLD, global_rank, ierr)

  !--set up a 1d cartesian topology 
call mpi_cart_create (MPI_COMM_WORLD, 1, (/ nproc /), (/ .false. /),  &
  .true., comm, ierr)

$endif

!--check if run-time number of processes agrees with nproc parameter
if (np /= nproc) then
  write (*, *) 'runtime number of procs = ', np,  &
               ' not equal to nproc = ', nproc
  stop
end if

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
subroutine mpi_sync_real_array( var, lbz, isync )
!**********************************************************************
! 
! This subroutine provides a generic method for syncing arrays in
! lesgo. This method applies to arrays indexed in the direction starting
! from both 0 and 1. For arrays starting from index of 1, only the
! SYNC_DOWN procedure may be performed. No assumption is made about the
! dimensions of the other directions (x and y) and the bounds of these
! indices are obtained directly.
!
! The synchronization is provided according to the following rules:
!
! SYNC_DOWN : Send data from k = 1 at coord+1 to k=nz at coord
! SYNC_UP   : Send data from k = nz-1 at coord to k=0 at coord+1
!
! Arguments:
!
! var   : three dimensional array to be sync'd accross processors 
! lbz   : the lower bound of var for the z index; its specification resolves
!         descrepencies between arrays indexed starting at 0 from those at 1
! isync : flag for determining the type of synchronization and takes on values, 
!         MPI_SYNC_DOWN, MPI_SYNC_UP, or MPI_SYNC_DOWNUP from the MPI_DEFS 
!         module. 
!
use types, only : rprec
use mpi
use param, only : MPI_RPREC, down, up, comm, status, ierr, nproc, coord, nz
use messages

implicit none

character (*), parameter :: sub_name = mod_name // '.mpi_sync_real_array'

real(rprec), dimension(:,:,lbz:), intent(INOUT) :: var
integer, intent(in) :: lbz
integer, intent(in) :: isync

!integer :: lbx,ubx
!integer :: lby,uby
integer :: sx, sy
integer :: ubz
integer :: mpi_datasize

!lbx=lbound(var,1); ubx=ubound(var,1)
!lby=lbound(var,2); uby=ubound(var,2)
!lbz=lbound(var,3); ubz=ubound(var,3)

! Get size of var array in x and y directions
sx=size(var,1)
sy=size(var,2)
! Get upper bound of z index; the lower bound is specified
ubz=ubound(var,3)

! We are assuming that the array var has nz as the upper bound - checking this
if( ubz .ne. nz ) call error( sub_name, 'Input array must lbz:nz z dimensions.')

!  Set mpi data size
mpi_datasize = sx*sy

if(isync == MPI_SYNC_DOWN) then

   call sync_down()

elseif( isync == MPI_SYNC_UP) then

   if( lbz /= 0 ) call error( sub_name, 'Cannot SYNC_UP with variable with non-zero starting index')
   call sync_up()

elseif( isync == MPI_SYNC_DOWNUP) then

   if( lbz /= 0 ) call error( sub_name, 'Cannot SYNC_DOWNUP with variable with non-zero starting index')
   call sync_down()
   call sync_up()

else

   call error( sub_name, 'isync not specified properly')

endif

if(ierr .ne. 0) call error( sub_name, 'Error occured during mpi sync; recieved mpi error code :', ierr)

! Enforce globally synchronous MPI behavior. Most likely safe to comment
! out, but can be enabled to ensure absolute safety.
!call mpi_barrier( comm, ierr )

return

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sync_down()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

call mpi_sendrecv (var(:,:,1), mpi_datasize, MPI_RPREC, down, 1,  &
  var(:,:,ubz), mpi_datasize, MPI_RPREC, up, 1,   &
  comm, status, ierr)
    
return

end subroutine sync_down

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sync_up()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

call mpi_sendrecv (var(:,:,ubz-1), mpi_datasize, MPI_RPREC, up, 2,  &
  var(:,:,0), mpi_datasize, MPI_RPREC, down, 2,   &
  comm, status, ierr)

return

end subroutine sync_up 

end subroutine mpi_sync_real_array

end module mpi_defs

