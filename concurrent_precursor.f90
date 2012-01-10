!*********************************************************************
module concurrent_precursor
!*********************************************************************
use types, only : rprec
implicit none

save
private

public :: interComm, color, RED, BLUE
public :: vel_sample_t
public :: create_mpi_comms_cps, &
          initialize_cps, &
          synchronize_cps, &
          inflow_cond_cps 

character (*), parameter :: mod_name = 'concurrent_precursor'

integer, parameter :: RED=0 ! Upstream domain (producer)
integer, parameter :: BLUE=1 ! Downstream domain (consumer) 

integer :: interComm, color

type vel_sample_type
   integer :: nx
   integer :: istart
   integer :: iplateau
   integer :: iend
   real(rprec), allocatable, dimension(:,:,:) :: u, v, w
end type vel_sample_type

type(vel_sample_type), target :: vel_sample_t 

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine create_mpi_comms_cps( localComm )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine does two things. It first splits the MPI_COMM_WORLD
! communicator into two communicators (localComm). The two new
! communicators are then bridged to create an intercommunicator
! (interComm).
! 
use mpi
use param, only : ierr
implicit none

integer, intent(out) :: localComm

integer :: world_np, world_rank
integer :: remoteLeader
integer :: memberKey

! Get number of processors in world comm
call mpi_comm_size (MPI_COMM_WORLD, world_np, ierr)
call mpi_comm_rank (MPI_COMM_WORLD, world_rank, ierr)

! Set color and remote leader for intercommunicator interComm
if( world_rank < world_np / 2 ) then
   color = RED
   remoteLeader = world_np / 2
else
   color = BLUE
   remoteLeader = 0
endif

! Generate member key
memberKey=modulo(world_rank, world_np / 2)

! Split the world communicator into intracommunicators localComm
call MPI_Comm_split(MPI_COMM_WORLD, color, memberKey, localComm, ierr)

! Create intercommunicator interComm
call mpi_intercomm_create( localComm, 0, MPI_COMM_WORLD, remoteLeader, 1, interComm, ierr)

return
end subroutine create_mpi_comms_cps

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine initialize_cps()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : nx, ny, nz
use param, only : coord, rank_of_coord, status, ierr
use param, only : fringe_region_end, fringe_region_len
use messages
use mpi
use fringe_util, only : fringe_init
implicit none

character (*), parameter :: sub_name = mod_name // '.initialize_cps'

!integer :: rankTest, coordTest

real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
integer, pointer :: nx_p, istart_p, iplateau_p, iend_p

nullify( u_p, v_p, w_p )
nullify( nx_p, istart_p, iplateau_p, iend_p )

istart_p   => vel_sample_t % istart
iplateau_p => vel_sample_t % iplateau
iend_p     => vel_sample_t % iend
nx_p       => vel_sample_t % nx
u_p        => vel_sample_t % u
v_p        => vel_sample_t % v
w_p        => vel_sample_t % w

if( color == BLUE ) then

   call fringe_init( istart_p, iplateau_p, iend_p )

   ! Sample size same as buffer region (omitting istart from the block
   ! since velocity is already set there)
   nx_p = iend_p - istart_p
   
   ! Send size of the sample block to upstream domain (RED)
   call mpi_send( nx_p , 1, MPI_INT, &
        rank_of_coord(coord), 1, interComm, ierr )

elseif( color == RED ) then

   ! Receive from downstream domain (BLUE) the length of the sample block
   call mpi_recv( nx_p , 1, MPI_INT, &
        rank_of_coord(coord), 1, interComm, status, ierr)
   
   ! Should end up as nx + 1 (this eventually gets wrapped) 
   iend_p = floor ( 1.0_rprec * nx + 1.0_rprec) 
   ! Plateau location not used since no fringe treatment on the RED domain, but setting so it is at
   ! least initialized.
   iplateau_p = iend_p
   ! Set istart based on the size of the sample block
   istart_p = iend_p - nx_p

else

  call error(sub_name,'Erroneous color specification')

endif

! Allocate the sample block
allocate( u_p( nx_p, ny, nz ) )
allocate( v_p( nx_p, ny, nz ) )
allocate( w_p( nx_p, ny, nz ) )

nullify( u_p, v_p, w_p )
nullify( nx_p, istart_p, iplateau_p, iend_p )

return
end subroutine initialize_cps

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine synchronize_cps()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use messages
use param, only : nx, ny, nz
use param, only : coord, rank_of_coord, status, ierr, MPI_RPREC
use sim_param, only : u,v,w
implicit none

character (*), parameter :: sub_name = mod_name // '.synchronize_cps'

real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
integer, pointer :: nx_p, istart_p, iend_p

integer :: i, i_w, index
integer :: sendsize, recvsize

nullify( u_p, v_p, w_p )
nullify( nx_p, istart_p, iend_p )

istart_p => vel_sample_t % istart
iend_p   => vel_sample_t % iend
nx_p     => vel_sample_t % nx
u_p      => vel_sample_t % u
v_p      => vel_sample_t % v
w_p      => vel_sample_t % w

index=0

sendsize = nx_p * ny * nz
recvsize = sendsize

if( color == BLUE ) then

   ! Recieve sampled velocities from upstream (RED)
   call mpi_recv( u_p(1,1,1) , recvsize, MPI_RPREC, &
        rank_of_coord(coord), 1, interComm, status, ierr)
   call mpi_recv( v_p(1,1,1) , recvsize, MPI_RPREC, &
        rank_of_coord(coord), 2, interComm, status, ierr)
   call mpi_recv( w_p(1,1,1) , recvsize, MPI_RPREC, &
        rank_of_coord(coord), 3, interComm, status, ierr)
   
elseif( color == RED ) then

   ! Remember, omitting istart since velocity is already known there
   do i=istart_p+1, iend_p

      index=index+1
      i_w = modulo( i - 1, nx ) + 1

      ! Sample velocity and copy to buffers
      u_p(index,:,:) = u( i_w, 1:ny, 1:nz)
      v_p(index,:,:) = v( i_w, 1:ny, 1:nz)
      w_p(index,:,:) = w( i_w, 1:ny, 1:nz)

   enddo

   ! Send sampled velocities to downstream domain (BLUE)
   call mpi_send( u_p(1,1,1), sendsize, MPI_RPREC, &
        rank_of_coord(coord), 1, interComm, ierr )
   call mpi_send( v_p(1,1,1), sendsize, MPI_RPREC, &
        rank_of_coord(coord), 2, interComm, ierr )
   call mpi_send( w_p(1,1,1), sendsize, MPI_RPREC, &
        rank_of_coord(coord), 3, interComm, ierr )
   
else

   call error( sub_name, 'Erroneous color specification')
   
endif

nullify( u_p, v_p, w_p )
nullify( nx_p, istart_p, iend_p )

return
end subroutine synchronize_cps

!**********************************************************************
subroutine inflow_cond_cps ()
!**********************************************************************
!
!  Enforces prescribed inflow condition from an inlet velocity field
!  generated from a precursor simulation. The inflow condition is
!  enforced by direct modulation on the velocity in the fringe region.
!
use types, only : rprec
use param, only : nx, ny, nz
use sim_param, only : u, v, w
use messages, only : error
use fringe_util, only : fringe_weighting
implicit none

character (*), parameter :: sub_name = 'inflow_cond_cps'

integer :: i, i_w, index
real (rprec) :: alpha, beta

integer, pointer :: nx_p, istart_p, iplateau_p, iend_p
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p

nullify( u_p, v_p, w_p )
nullify( nx_p, istart_p, iplateau_p, iend_p )

u_p        => vel_sample_t % u
v_p        => vel_sample_t % v
w_p        => vel_sample_t % w
istart_p   => vel_sample_t % istart
iplateau_p => vel_sample_t % iplateau
iend_p     => vel_sample_t % iend

index=0

!--skip istart since we know vel at istart already (where beta=0, alpha=1)
do i = istart_p + 1, iend_p 

  i_w = modulo (i - 1, nx) + 1

  beta = fringe_weighting( i, istart_p, iplateau_p )
  alpha = 1.0_rprec - beta

  index = index + 1

  u(i_w, 1:ny, 1:nz) = alpha * u(i_w, 1:ny, 1:nz) + beta * u_p(index, 1:ny, 1:nz) 
  v(i_w, 1:ny, 1:nz) = alpha * v(i_w, 1:ny, 1:nz) + beta * v_p(index, 1:ny, 1:nz)
  w(i_w, 1:ny, 1:nz) = alpha * w(i_w, 1:ny, 1:nz) + beta * w_p(index, 1:ny, 1:nz)
  
end do

if( index .ne. nx_p ) call error( sub_name, 'Mismatch in expected sample size')

nullify( u_p, v_p, w_p )
nullify( nx_p, istart_p, iplateau_p, iend_p )

return
end subroutine inflow_cond_cps


end module concurrent_precursor
