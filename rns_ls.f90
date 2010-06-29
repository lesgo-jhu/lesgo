!**********************************************************************
module rns_ls
!**********************************************************************
use rns_base_ls
!!$if($CYL_SKEW_LS)
!!use cyl_skew_base_ls, only : tr_t, clindx_to_loc_id, brindx_to_loc_id, ntree
!!use cyl_skew_ls, only : cyl_skew_fill_tree_array_ls, ngen, ngen_reslv
!!!use cyl_skew_ls, only : cyl_skew_fill_cl_ref_plane_array_ls
!!!use cyl_skew_ls, only : cyl_skew_get_branch_id_ls
!!$endif

implicit none

save
private

public :: rns_CD_ls
public :: rns_forcing_ls ! Apply forcing
public :: rns_finalize_ls

character (*), parameter :: mod_name = 'rns_ls'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine rns_CD_ls()
!**********************************************************************
!  This subroutine handles all CD calculation within the RNS module; 
!  all CD and force calculations associated with
!
!  r_elem
!  beta_elem
!
!  are handled here
!
use param, only : jt, USE_MPI, coord
use messages
!!$if($CYL_SKEW_LS)
!!use cyl_skew_base_ls, only : ngen, ngen_reslv
!!$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_CD_ls'

call r_elem_force_ls()     !  Get CD, force etc for resolved regions
call beta_elem_force_ls()  !  Get force of 
    
  !if(ngen > ngen_reslv) call rns_cl_unreslv_CD_ls()
    
  if(modulo (jt, clforce_nskip) == 0) then
    
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then

      call write_cl_CD_ls()
      call write_cl_fD_ls()
      call write_cl_vel_ls()
      
      call write_beta_CD_ls()
      call write_beta_fD_ls()
      call write_beta_vel_ls()
      call write_beta_kappa_ls()
      
    endif
    
  endif
  
!endif

!if(coord == 0) call mesg(sub_name, 'Exiting ' // sub_name)

return
end subroutine rns_CD_ls


!**********************************************************************
subroutine r_elem_force_ls()
!**********************************************************************
!  This subroutine computes the CD of the all the resolved elements.
!
use types, only : rprec
use messages
use param, only : nx, ny, nz, dx, dy, dz, coord, USE_MPI
$if($MPI)
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif
use sim_param, only : u
use functions, only : points_avg_3D
use immersedbc, only : fx
implicit none

character (*), parameter :: sub_name = mod_name // '.r_elem_force_ls'

integer, pointer :: i, j, k
integer, pointer :: npoint_p
integer, pointer, dimension(:,:) :: iarray_p

real(rprec), pointer :: fD_p

integer :: ncluster_tot
integer :: nt, ng, nc, np

$if ($MPI)
real(rprec) :: fD
$endif

type(ref_region) :: ref_region_t_p
type(indx_array) :: indx_array_t_p

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

!  Comment starts here 
nullify(ref_region_t_p)
nullify(indx_array_t_p)
nullify(npoint_p, iarray_p)
nullify(i,j,k)
nullify(fD_p)

!!$if ($MPI)
!!allocate (cl_fD ( ncluster_reslv_ref ) )
!!cl_fD = 0._rprec
!!$endif

do n = 1, nr_elem

  !  Get the reference velocity
  ref_region_t_p => r_elem_t( n ) % ref_region_t
  ref_region_t_p % u = points_avg_3D( u(1:nx,:,1:nz), ref_region_t_p % points)
  
  indx_array_t_p => r_elem_t( n ) % indx_array_t
     
  npoint_p => indx_array_t_p % npoint
  iarray_p => cindx_array_t_p % iarray
  
  $if($MPI)
  fD = 0._rprec
  $endif

  fD_p => r_elem_t( n ) % force_t % fD
  fD_p = 0._rprec
  
  do np=1, npoint_p
  
    i => iarray_p(1,np)
    j => iarray_p(2,np)
    k => iarray_p(3,np)
  
    $if($MPI)
    fD = fD + fx(i,j,k) * dx * dy * dz
    $else
    fD_p = fD_p + fx(i,j,k) * dx * dy * dz
    $endif
    
    nullify(i,j,k)
    
  enddo
  
  $if($MPI)
  call mpi_allreduce (fD, fD_p, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  $endif
  
  !  Compute CD
  r_elem_t(n) % force_t % CD = -fD_p / (0.5_rprec * ref_region_t_p(n)%area * (ref_region_t_p(n)%u)**2)
  
  nullify(fD_p)
  nullify(npoint_p, iarray_p)
  nullify(indx_array_t_p)
  nullify(ref_region_t_p)
 
enddo

return
end subroutine r_elem_force_ls

!**********************************************************************
subroutine rns_elem_force_ls()
!**********************************************************************
!  This subroutine computes the CD of the beta elements
!
use types, only : rprec
use messages
use sim_param, only : u
use immersedbc, only : fx
use functions, only : points_avg_3D
$if($MPI)
use mpi
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_elem_force_ls'

integer, pointer :: i,j,k, n, ns
integer, pointer :: npoint_p
integer, pointer :: clindx_p
integer, pointer :: nelem_p, indx_p

real(rprec), allocatable, dimension(:) :: beta_gamma
real(rprec), allocatable, dimension(:) :: beta_gamma_sum
real(rprec), allocatable, dimension(:) :: b_gamma

real(rprec), pointer :: area_p, u_p
real(rprec), pointer :: kappa_p, CD_p
real(rprec), pointer, dimension(:,:) :: points_p

!real(rprec) :: sigma
!real(rprec), allocatable, dimension(:) :: fD_dir

real(rprec) :: CD, Lint

real(rprec), allocatable, dimension(:) ::  fD_tot
real(rprec), allocatable, dimension(:) ::  CD_num, CD_denom

$if($MPI)
real(rprec) :: Lint_global
$endif

$if($MPI)
real(rprec) :: fD
$endif

!if(coord == 0) call mesg(sub_name, 'Entered ' // sub_name)

nullify(i,j,k)
nullify(npoint_p)
nullify(clindx_p)
nullify(nelem_p, indx_p)

nullify(area_p, u_p)
nullify(kappa_p, CD_p)

allocate(beta_gamma(nbeta_elem))
allocate(b_gamma(nb_elem))
allocate(beta_gamma_sum(nb_elem))
beta_gamma_sum=0

do n=1, nbeta_elem

  u_p      => beta_elem_t(n) % ref_region_t % u
  area_p   => beta_elem_t(n) % ref_region_t % area
  points_p => beta_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3D( u(1:nx,:,1:nz), points_p ) 
  
  beta_gamma(n) = dabs( u_p ) * u_p * area_p
  
  nullify(u_p, area_p, points_p)
  
enddo

do n=1, nb_elem

  u_p      => b_elem_t(n) % ref_region_t % u
  area_p   => b_elem_t(n) % ref_region_t % area
  points_p => b_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3D( u(1:nx,:,1:nz), points_p ) 
  
  b_gamma(n) = dabs( u_p ) * u_p * area_p
  
  nullify(u_p, area_p, points_p)
  
  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx
  
  do ns=1, nelem_p
    beta_gamma_sum(n) = beta_gamma_sum(n) + beta_gamma( indx_p(ns) )
  enddo
  
  nullify(indx_p, nelem_p)
  
enddo

!  Compute the total resolved force of each b_elem
allocate(b_r_force( nb_elem )) 
b_r_force = 0._rprec
do n=1, nb_elem

  nelem_p => b_elem_t(n) % r_child_t % nelem
  indx_p  => b_elem_t(n) % r_child_t % indx
  
  do ns=1, nelem_p
    b_r_force(n) = b_r_force(n) + r_elem_t( indx_p(ns) )% force_t % fD 
  enddo 
	
  nullify(indx_p, nelem_p)
	
enddo

if( use_explicit_formulation ) then

  allocate(b_force( nb_elem ))
  b_force = 0._rprec
  
  do n=1, nb_elem
    	
	b_force(n) = b_r_force(n) - 0.5_rprec * b_elem_t(n) % force_t % CD * beta_gamma_sum(n)
	
  enddo
  
  if( use_local_CD ) then
  
    do n=1, nb_elem
	
	  b_elem_t(n) % force_t % CD = 2._rprec * b_force(n) / b_gamma(n)
	  
	enddo
  
  else ! compute global CD
  
    CD_num=0._rprec
	CD_denom=0._rprec
    do n=1, nb_elem
	
	  CD_num = CD_num + b_force(n) * b_gamma(n) 
	  CD_denom = b_gamma(n) * b_gamma(n)
	  
	enddo

	b_elem_t(:) % force_t % CD = 2._rprec * CD_num / CD_denom
	
  endif
  
  deallocate(b_force)

else ! use implicit formulation

  allocate( b_m( nb_elem )) 
  do n=1, nb_elem  
    b_m(n) = 0.5_rprec * (beta_gamma_sum(n) - b_gamma(n))
  enddo

  if( use_local_CD ) then
  
    do n=1, nb_elem
	
	  b_elem_t(n) % force_t % CD = b_r_force(n) / b_m(n)
	  
	enddo
	
  else
  
    CD_num=0._rprec
	CD_denom=0._rprec
    do n=1, nb_elem
	
	  CD_num = CD_num + b_r_force(n) * b_m(n)
	  CD_denom = b_m(n) * b_m(n)
	  
	enddo

	b_elem_t(:) % force_t % CD = CD_num / CD_denom
	
  endif
  
  deallocate( b_m )
 
endif

deallocate(b_r_force)

!  Now update the CD of all the beta_elem 
do n=1, nb_elem

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

!  Compute the total force for b_elem
do n=1, nb_elem
  b_elem_t(n) % force_t % fD = - 0.5_rprec * b_elem_t(n) % force_t % CD * b_gamma(n)
enddo

!  Now need to compute kappa
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
    
  if(coord == 0 .and. (modulo (jt, 10) == 0)) write(*,'(1a,i3,3f18.6)') 'beta_indx, kappa, CD, beta_int : ', n, kappa_p, CD_p, beta_int
    
  nullify(kappa_p, CD_p)
  nullify(u_p, area_p)
        
enddo
   
deallocate(beta_gamma, beta_gamma_sum)
deallocate(b_gamma)

return
end subroutine rns_elem_force_ls

!**********************************************************************
subroutine rns_forcing_ls()
!**********************************************************************
!  This subroutine computes the forces on the unresolved branches
!
use types, only : rprec
use sim_param, only : u
use immersedbc, only : fx
$if($MPI)
use mpi
use param, only : MPI_RPREC, up, down, comm, status, ierr, ld, ny, nz, nproc
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif
use param, only : dx, dy, dz, coord, jt, USE_MPI

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_forcing_ls'

integer :: ib, np

integer, pointer :: i, j, k

integer, pointer :: npoint_p

real(rprec), pointer :: kappa_p

nullify(i,j,k)
nullify(npoint_p)
nullify(kappa_p)


do n = 1, nbeta_elem
 
    !  Loop over number of points used in beta region
    npoint_p => beta_force_t( n ) % indx_array_t % npoint    
    kappa_p  => beta_force_t( n ) % force_t % kappa
  
    do np = 1, npoint_p
  
      i => beta_force_t( n ) % indx_array_t % iarray(1,np)
      j => beta_force_t( n ) % indx_array_t % iarray(2,np)
      k => beta_force_t( n ) % indx_array_t % iarray(3,np)
    
      fx(i,j,k) = - kappa_p * dabs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
      nullify(i,j,k)
      
    enddo
    
    nullify( npoint_p, kappa_p )
    
enddo

$if($MPI)
call mpi_sync_real_array( fx, MPI_SYNC_DOWNUP )
$endif

return

end subroutine rns_forcing_ls

!**********************************************************************
subroutine write_r_elem_data()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : total_time, path
use strmod

character(*), parameter :: sub_name = mod_name // '.write_r_elem_data'
character(*), parameter :: fname_CD = path // 'output/rns_r_elem_CD.dat'
character(*), parameter :: fname_fD = path // 'output/rns_r_elem_fD.dat'
character(*), parameter :: fname_vel = path // 'output/rns_r_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n, nvar

inquire (file=fname_CD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nr_elem
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
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
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nr_elem+1, (/ total_time, r_elem_t(:) % force_t % fD /))

inquire (file=fname_vel, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nr_elem
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', 'formatted', nr_elem+1, (/ total_time, r_elem_t(:) % ref_region_t % u /))

return
end subroutine write_r_elem_data

!**********************************************************************
subroutine write_beta_elem_data()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : total_time, path
use strmod

character(*), parameter :: sub_name = mod_name // '.write_beta_elem_data'
character(*), parameter :: fname_CD = path // 'output/rns_beta_elem_CD.dat'
character(*), parameter :: fname_fD = path // 'output/rns_beta_elem_fD.dat'
character(*), parameter :: fname_kappa = path // 'output/rns_beta_elem_kappa.dat'
character(*), parameter :: fname_vel = path // 'output/rns_beta_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n, nvar

inquire (file=fname_CD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
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
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nbeta_elem + 1, (/ total_time, beta_elem_t(:) % force_t % fD /))

inquire (file=fname_fD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"<greek>k</greek><sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nbeta_elem + 1, (/ total_time, beta_elem_t(:) % force_t % kappa /))

inquire (file=fname_vel, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nbeta_elem
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', 'formatted', nbeta_elem+1, (/ total_time, beta_elem_t(:) % ref_region_t % u /))

return
end subroutine write_beta_elem_data

!**********************************************************************
subroutine write_b_elem_data()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : total_time, path
use strmod

character(*), parameter :: sub_name = mod_name // '.write_b_elem_data'
character(*), parameter :: fname_CD = path // 'output/rns_b_elem_CD.dat'
character(*), parameter :: fname_fD = path // 'output/rns_b_elem_fD.dat'
character(*), parameter :: fname_vel = path // 'output/rns_b_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n, nvar

inquire (file=fname_CD, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nb_elem
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
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
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname_fD, 'append', 'formatted', nb_elem+1, (/ total_time, b_elem_t(:) % force_t % fD /))

inquire (file=fname_vel, exist=exst)
if (.not. exst) then
  var_list = '"t"'
  do n = 1, nb_elem
    !  Create variable list name:
    call strcat(var_list, ',"u<sub>')
    call strcat(var_list, n)
    call strcat(var_list, '</sub>"')
  enddo
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

call write_real_data(fname, 'append', 'formatted', nb_elem+1, (/ total_time, r_elem_t(:) % ref_region_t % u /))

return
end subroutine write_r_elem_data

!**********************************************************************
subroutine rns_force_init_ls ()
!**********************************************************************
!  
!  This subroutine reads the last BETA force data from a previous simulation
!
use param, only : coord, USE_MPI
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_force_init_ls'
character (*), parameter :: fname_in = 'rns_force_ls.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fname
$endif

integer :: ip

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

open (1, file=fname, action='read', position='rewind',  &
  form='unformatted')
read (1) beta_force_t
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
character (*), parameter :: fname_out = 'rns_force_ls.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fname
$endif

integer :: ip

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') fname_out, MPI_suffix, coord
$else
fname = trim(adjustl(fname_out))
$endif

open (1, file=fname, action='write', position='rewind',  &
  form='unformatted')
write (1) beta_force_t
close (1)

return
end subroutine rns_finalize_ls

end module rns_ls


