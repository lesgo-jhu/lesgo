!**********************************************************************
module rns_ls
!**********************************************************************
use rns_base_ls

implicit none

save
private

public :: rns_forcing_ls ! Apply forcing
public :: rns_finalize_ls
public :: rns_force_init_ls

character (*), parameter :: mod_name = 'rns_ls'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine rns_forcing_ls()
!**********************************************************************
!  This subroutine computes the forces on the unresolved branches
!
use types, only : rprec
use sim_param, only : u
use immersedbc, only : fx

$if($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use mpi
use param, only : ld, ny, nz, MPI_RPREC, down, up, comm, status, ierr
$endif

use param, only : dx, dy, dz, coord, jt, USE_MPI

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_forcing_ls'

integer :: np, n

integer, pointer :: i, j, k

integer, pointer :: npoint_p
real(rprec), pointer :: kappa_p

nullify(i,j,k)
nullify(npoint_p)
nullify(kappa_p)

!  Compute the relavent force information ( include reference quantities, CD, etc.)
call rns_elem_force()

!  Apply the RNS forcing to appropriate nodes
do n = 1, nbeta_elem
 
    !  Loop over number of points used in beta region
    npoint_p => beta_elem_t(n) % indx_array_t % npoint    
    kappa_p  => beta_elem_t(n) % force_t % kappa
  
    do np = 1, npoint_p
  
      i => beta_elem_t( n ) % indx_array_t % iarray(1,np)
      j => beta_elem_t( n ) % indx_array_t % iarray(2,np)
      k => beta_elem_t( n ) % indx_array_t % iarray(3,np)
    
      fx(i,j,k) = - kappa_p * dabs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
      nullify(i,j,k)
      
    enddo
    
    nullify( npoint_p, kappa_p )
    
enddo

!$if($MPI)
!call mpi_sync_real_array( fx, MPI_SYNC_DOWNUP )
!$endif

$if($MPI)
!  Sync fx; can't use mpi_sync_real_array since its not allocated from 0 -> nz
call mpi_sendrecv (fx(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
  fx(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
  comm, status, ierr)
$endif

!endif

if(modulo (jt, output_nskip) == 0) call rns_elem_output()

return

end subroutine rns_forcing_ls

!**********************************************************************
subroutine rns_elem_output()
!**********************************************************************
!
use param, only : jt, USE_MPI, coord
use messages
!!$if($CYL_SKEW_LS)
!!use cyl_skew_base_ls, only : ngen, ngen_reslv
!!$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_elem_output'
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then

  call r_elem_data_write()
  call beta_elem_data_write()
  call b_elem_data_write()
      
endif
    
return
end subroutine rns_elem_output


!**********************************************************************
subroutine rns_elem_force()
!**********************************************************************
!  This subroutine computes the CD of the beta elements
!
use types, only : rprec
use messages
use sim_param, only : u
use immersedbc, only : fx
use functions, only : points_avg_3D
use param, only : nx, nz, dx, dy, dz, coord, jt
$if($MPI)
use mpi
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_elem_force'

integer ::  n, ns
integer, pointer :: nelem_p
integer, pointer, dimension(:) :: indx_p

real(rprec), allocatable, dimension(:) :: beta_gamma
real(rprec), allocatable, dimension(:) :: beta_gamma_sum
!real(rprec), allocatable, dimension(:) :: beta_gamma_CD_sum
real(rprec), allocatable, dimension(:) :: b_gamma
real(rprec), allocatable, dimension(:) :: b_r_force, b_force, b_m

real(rprec), pointer :: area_p, u_p
real(rprec), pointer, dimension(:,:) :: points_p

real(rprec) ::  CD_num, CD_denom

integer, pointer :: i,j,k
integer, pointer :: npoint_p

real(rprec), pointer :: kappa_p, CD_p

real(rprec) :: beta_int

$if($MPI)
real(rprec) :: beta_int_global
$endif

$if($MPI)
real(rprec) :: fD
$endif

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

nullify(nelem_p, indx_p)

nullify(area_p, u_p)
nullify(points_p)

allocate(beta_gamma(nbeta_elem))
allocate(b_gamma(nb_elem))
allocate(beta_gamma_sum(nb_elem))
!allocate(beta_gamma_CD_sum(nb_elem))

beta_gamma(:)=0._rprec
b_gamma(:)=0._rprec
beta_gamma_sum(:)=0._rprec
!beta_gamma_CD_sum(:)=0._rprec

$if($MPI)
!  Make sure intermediate velocity is sync'd
call mpi_sync_real_array( u, MPI_SYNC_DOWNUP)
$endif

!  Get the force for the resolved elements
call r_elem_force()

do n=1, nbeta_elem

  u_p      => beta_elem_t(n) % ref_region_t % u
  area_p   => beta_elem_t(n) % ref_region_t % area
  points_p => beta_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3D( u(1:nx,:,1:nz), points_p ) 
  
  beta_gamma(n) = dabs( u_p ) * u_p * area_p
  
  nullify(points_p, area_p, u_p)

enddo

do n=1, nb_elem

  u_p      => b_elem_t(n) % ref_region_t % u
  area_p   => b_elem_t(n) % ref_region_t % area
  points_p => b_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3D( u(1:nx,:,1:nz), points_p ) 
  
  b_gamma(n) = dabs( u_p ) * u_p * area_p
  
  nullify(points_p, area_p, u_p)
  
  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx
  
  do ns=1, nelem_p
    beta_gamma_sum(n) = beta_gamma_sum(n) + beta_gamma( indx_p(ns) )
    !beta_gamma_CD_sum(n) = beta_gamma_CD_sum(n) + beta_elem_t( indx_p(ns) ) % force_t % CD * beta_gamma( indx_p(ns) )
  enddo
  
  nullify(indx_p, nelem_p)
  
enddo

!  Compute the total resolved force of each b_elem
allocate(b_r_force( nb_elem )) 
b_r_force(:) = 0._rprec

do n=1, nb_elem

  nelem_p => b_elem_t(n) % r_child_t % nelem
  indx_p  => b_elem_t(n) % r_child_t % indx
  
  do ns=1, nelem_p
    b_r_force(n) = b_r_force(n) + r_elem_t( indx_p(ns) )% force_t % fD 
  enddo 
  
  !  Perform partial sum for b_elem force
  b_elem_t(n) % force_t % fD = b_r_force(n) 	
  
  nullify(indx_p, nelem_p)
	
enddo

if( time_model == 1 ) then

  allocate(b_force( nb_elem ))
  !b_force = 0._rprec
  
  do n=1, nb_elem
    	
    b_force(n) = b_r_force(n) - 0.5_rprec * b_elem_t(n) % force_t % CD * beta_gamma_sum(n)
    !b_force(n) = b_r_force(n) - 0.5_rprec * beta_gamma_CD_sum(n)
	
  enddo
  
  if( spatial_model == 1 ) then
  
    do n=1, nb_elem
	
      b_elem_t(n) % force_t % CD = -2._rprec * b_force(n) / b_gamma(n)
	  
    enddo
	
  elseif( spatial_model == 2) then
   
    CD_num=0._rprec
	CD_denom=0._rprec
    
    do n=1, nb_elem
	
      CD_num = CD_num + b_force(n) 
      CD_denom = CD_denom + b_gamma(n)
      
    enddo
    
    b_elem_t(:) % force_t % CD = -2._rprec * CD_num / CD_denom    
  
  elseif( spatial_model == 3) then ! compute global CD
  
    CD_num=0._rprec
	CD_denom=0._rprec
    
    do n=1, nb_elem
	
      CD_num = CD_num + b_force(n) * b_gamma(n) 
      CD_denom = CD_denom + b_gamma(n) * b_gamma(n)
      
    enddo
    
    b_elem_t(:) % force_t % CD = -2._rprec * CD_num / CD_denom
	
  else
  
    call error( sub_name, 'spatial_model not specified correctly.')
	
  endif
  
  deallocate(b_force)

elseif( time_model == 2) then ! use implicit formulation

  allocate( b_m( nb_elem ) ) 
  
  do n=1, nb_elem  
    !b_m(n) = beta_gamma_sum(n) - b_gamma(n)
    b_m(n) = b_gamma(n) - beta_gamma_sum(n)
  enddo

  if( spatial_model == 1 ) then
  
    do n=1, nb_elem
	
      b_elem_t(n) % force_t % CD = 2._rprec * b_r_force(n) / b_m(n)
	  
    enddo
	
  elseif( spatial_model == 2 ) then
  
    CD_num=0._rprec
    CD_denom=0._rprec
    
    do n=1, nb_elem
	
      CD_num = CD_num + b_r_force(n) 
      CD_denom = CD_denom + b_m(n)
	  
    enddo

    b_elem_t(:) % force_t % CD = 2._rprec * CD_num / CD_denom    
	
  elseif( spatial_model == 3 ) then
  
    CD_num=0._rprec
    CD_denom=0._rprec
    
    do n=1, nb_elem
	
      CD_num = CD_num + b_r_force(n) * b_m(n)
      CD_denom = CD_denom + b_m(n) * b_m(n)
	  
    enddo

    b_elem_t(:) % force_t % CD = 2._rprec * CD_num / CD_denom
	
  else
  
    call error( sub_name, 'spatial_model not specified correctly.')	
	
  endif
  
  deallocate( b_m )
  
else
  
  call error( sub_name, 'time_model not specified correctly.')  
 
endif

deallocate(b_r_force)

!  Check if CD is to be modulated
if( jt < CD_ramp_nstep ) b_elem_t(:) % force_t % CD = b_elem_t(:) % force_t % CD * jt / CD_ramp_nstep

!  Now update the CD of all the beta_elem 
do n=1, nb_elem

  !  Check if b_elem CD < 0
  if( b_elem_t(n) % force_t % CD < 0._rprec ) b_elem_t(n) % force_t % CD = 0._rprec

  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx
  
  do ns=1, nelem_p
    beta_elem_t( indx_p(ns) ) % force_t % CD = b_elem_t(n) % force_t % CD
  enddo 
	
  nullify(indx_p, nelem_p)
	
enddo

!  Compute the total force for beta_elem
do n=1, nbeta_elem
  beta_elem_t(n) % force_t % fD = - 0.5_rprec * beta_elem_t(n) % force_t % CD * beta_gamma(n)
enddo

!  Compute the total force for b_elem; r_elem has already been accounted for from above
do n=1, nb_elem
  
  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx

  !  Perform secondary sum over beta_elem for b_elem force
  do ns = 1, nelem_p
    b_elem_t(n) % force_t % fD = b_elem_t(n) % force_t % fD + beta_elem_t( indx_p(ns) ) % force_t % fD
  enddo 

  nullify( nelem_p, indx_p )
  
enddo

!  Now need to compute kappa; each beta region gets its own kappa value
do n = 1, nbeta_elem

  !  Compute beta_int over each region beta
  beta_int = 0._rprec
    
  $if($MPI)
  beta_int_global = 0._rprec
  $endif
  
  !  Loop over number of points used in beta region
  npoint_p => beta_elem_t(n) % indx_array_t % npoint
  
  do ns = 1, npoint_p
  
    i => beta_elem_t(n) % indx_array_t % iarray(1,ns)
    j => beta_elem_t(n) % indx_array_t % iarray(2,ns)
    k => beta_elem_t(n) % indx_array_t % iarray(3,ns)
    
    beta_int = beta_int + dabs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
    nullify(i,j,k)
      
  enddo
  
  beta_int = beta_int * dx * dy * dz
    
  nullify( npoint_p )
    
  $if($MPI)
  call mpi_allreduce (beta_int, beta_int_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  beta_int = beta_int_global
  $endif
    
  kappa_p => beta_elem_t(n) % force_t % kappa
  CD_p    => beta_elem_t(n) % force_t % CD
    
  kappa_p = CD_p * beta_gamma(n) / ( 2._rprec * beta_int )
    
  if(coord == 0 .and. (modulo (jt, output_nskip) == 0)) write(*,'(1a,i3,3f18.6)') 'beta_indx, kappa, CD, beta_int : ', n, kappa_p, CD_p, beta_int
    
  nullify(kappa_p, CD_p)

enddo
   
deallocate(beta_gamma, beta_gamma_sum)
deallocate(b_gamma)


return
end subroutine rns_elem_force

!**********************************************************************
subroutine r_elem_force()
!**********************************************************************
!  This subroutine computes the CD of the all the resolved elements.
!
use types, only : rprec
use messages
use param, only : nx, ny, nz, dx, dy, dz, coord, USE_MPI, jt
$if($MPI)
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif
use sim_param, only : u
use functions, only : points_avg_3D
use immersedbc, only : fx
implicit none

character (*), parameter :: sub_name = mod_name // '.r_elem_force'

integer :: n, np

integer, pointer :: i, j, k
integer, pointer :: npoint_p
integer, pointer, dimension(:,:) :: iarray_p

real(rprec), pointer :: fD_p

$if ($MPI)
real(rprec) :: fD
$endif

type(ref_region), pointer :: ref_region_t_p
type(indx_array), pointer :: indx_array_t_p

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

!  Comment starts here 
nullify(ref_region_t_p)
nullify(indx_array_t_p)
nullify(npoint_p, iarray_p)
nullify(i,j,k)
nullify(fD_p)

do n = 1, nr_elem

  !  Get the reference velocity
  ref_region_t_p => r_elem_t( n ) % ref_region_t
  ref_region_t_p % u = points_avg_3D( u(1:nx,:,1:nz), ref_region_t_p % points)
  
  indx_array_t_p => r_elem_t( n ) % indx_array_t
     
  npoint_p => indx_array_t_p % npoint
  iarray_p => indx_array_t_p % iarray
  
  $if($MPI)
  fD = 0._rprec
  $endif

  fD_p => r_elem_t( n ) % force_t % fD
  fD_p = 0._rprec
  
  do np=1, npoint_p
  
    i => iarray_p(1,np)
    j => iarray_p(2,np)
    k => iarray_p(3,np)
    
    if( k == nz ) call error( sub_name, 'Summing over bogus fx')
  
    $if($MPI)
    fD = fD + fx(i,j,k)
    $else
    fD_p = fD_p + fx(i,j,k)
    $endif
    
    nullify(i,j,k)
    
  enddo
  
  $if($MPI)
  call mpi_allreduce (fD, fD_p, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  $endif

  fD_p = fD_p * dx * dy * dz
   
  !  Compute CD
  r_elem_t(n) % force_t % CD = -fD_p / (0.5_rprec * ref_region_t_p % area * (ref_region_t_p % u)**2)
  
  if(coord == 0 .and. (modulo (jt, output_nskip) == 0)) write(*,'(1a,i3,3f18.6)') 'r_indx, fD, CD : ', n, -fD_p, r_elem_t(n) % force_t % CD

  nullify(fD_p)
  nullify(npoint_p, iarray_p)
  nullify(indx_array_t_p)
  nullify(ref_region_t_p)
 
enddo

return
end subroutine r_elem_force

!**********************************************************************
subroutine r_elem_data_write()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : total_time, path
use strmod

character(*), parameter :: sub_name = mod_name // '.r_elem_data_write'
character(*), parameter :: fname_CD = path // 'output/rns_r_elem_CD.dat'
character(*), parameter :: fname_fD = path // 'output/rns_r_elem_fD.dat'
character(*), parameter :: fname_vel = path // 'output/rns_r_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n

inquire (file=fname_CD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nr_elem
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_CD, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_CD, 'append', 'formatted', nr_elem + 1, (/ total_time, r_elem_t(:) % force_t % CD /))

inquire (file=fname_fD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nr_elem
    !  Create variable list name:
    call strcat(var_list, ',"fD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_fD, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nr_elem+1, (/ total_time, -r_elem_t(:) % force_t % fD /))

inquire (file=fname_vel, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nr_elem
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_vel, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_vel, 'append', 'formatted', nr_elem+1, (/ total_time, r_elem_t(:) % ref_region_t % u /))

return
end subroutine r_elem_data_write

!**********************************************************************
subroutine beta_elem_data_write()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : total_time, path
use strmod

character(*), parameter :: sub_name = mod_name // '.beta_elem_data_write'
character(*), parameter :: fname_CD = path // 'output/rns_beta_elem_CD.dat'
character(*), parameter :: fname_fD = path // 'output/rns_beta_elem_fD.dat'
character(*), parameter :: fname_kappa = path // 'output/rns_beta_elem_kappa.dat'
character(*), parameter :: fname_vel = path // 'output/rns_beta_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n

inquire (file=fname_CD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_CD, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_CD, 'append', 'formatted', nbeta_elem + 1, (/ total_time, beta_elem_t(:) % force_t % CD /))

inquire (file=fname_fD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"fD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_fD, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nbeta_elem + 1, (/ total_time, -beta_elem_t(:) % force_t % fD /))

inquire (file=fname_kappa, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"<greek>k</greek><sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_kappa, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_kappa, 'append', 'formatted', nbeta_elem + 1, (/ total_time, beta_elem_t(:) % force_t % kappa /))

inquire (file=fname_vel, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_vel, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_vel, 'append', 'formatted', nbeta_elem+1, (/ total_time, beta_elem_t(:) % ref_region_t % u /))

return
end subroutine beta_elem_data_write

!**********************************************************************
subroutine b_elem_data_write()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : total_time, path
use strmod

character(*), parameter :: sub_name = mod_name // '.b_elem_data_write'
character(*), parameter :: fname_CD = path // 'output/rns_b_elem_CD.dat'
character(*), parameter :: fname_fD = path // 'output/rns_b_elem_fD.dat'
character(*), parameter :: fname_vel = path // 'output/rns_b_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n
inquire (file=fname_CD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nb_elem
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_CD, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_CD, 'append', 'formatted', nb_elem + 1, (/ total_time, b_elem_t(:) % force_t % CD /))

inquire (file=fname_fD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nb_elem
    !  Create variable list name:
    call strcat(var_list, ',"fD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_fD, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nb_elem+1, (/ total_time, -b_elem_t(:) % force_t % fD /))

inquire (file=fname_vel, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nb_elem
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname_vel, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_vel, 'append', 'formatted', nb_elem+1, (/ total_time, b_elem_t(:) % ref_region_t % u /))

return
end subroutine b_elem_data_write

!**********************************************************************
subroutine rns_force_init_ls ()
!**********************************************************************
!  
!  This subroutine reads the last BETA force data from a previous simulation
!
use param, only : coord, USE_MPI
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_force_init'
character (*), parameter :: fname_in = 'rns_force.out'
character (128) :: fname
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

$endif

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') fname_in, MPI_suffix, coord
$else
fname = trim(adjustl(fname_in))
$endif

inquire (file=fname, exist=exst)

if (.not. exst) then
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,*) ' '
    write(*,*)'No previous RNS force data - starting from scratch.'
  endif
  return ! Do nothing if not present
endif 

$if ($READ_BIG_ENDIAN)
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted', convert='big_endian')
$elseif ($READ_LITTLE_ENDIAN)
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted', convert='little_endian')  
$else
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted')
$endif

read (1) r_elem_t(:) % force_t 
read (1) beta_elem_t(:) % force_t 
read (1) b_elem_t(:) % force_t 
close (1)

end subroutine rns_force_init_ls

!**********************************************************************
subroutine rns_finalize_ls()
!**********************************************************************
! 
!  This subroutine writes all restart data to file
!
use param, only : coord
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_finalize_ls'
character (*), parameter :: fname_out = 'rns_force.out'

character (128) :: fname
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

$endif

logical :: opn

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') fname_out, MPI_suffix, coord
$else
fname = trim(adjustl(fname_out))
$endif

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, action='write', position='rewind',  &
  form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, action='write', position='rewind',  &
  form='unformatted', convert='little_endian')  
$else
open (1, file=fname, action='write', position='rewind',  &
  form='unformatted')
$endif

write(1) r_elem_t(:) % force_t 
write(1) beta_elem_t(:) % force_t 
write(1) b_elem_t(:) % force_t
close (1)

deallocate(r_elem_t)
deallocate(beta_elem_t)
deallocate(b_elem_t)

return
end subroutine rns_finalize_ls

end module rns_ls


