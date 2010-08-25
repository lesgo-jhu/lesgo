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

contains

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
    
      fx(i,j,k) = - kappa_p * abs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
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
use param, only : nx, nz, dx, dy, dz, coord, jt, jt_total
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

nullify(nelem_p, indx_p)
nullify(area_p, u_p)
nullify(points_p)

allocate(beta_gamma(nbeta_elem))
allocate(b_gamma(nb_elem))
allocate(beta_gamma_sum(nb_elem))

beta_gamma(:)=0._rprec
b_gamma(:)=0._rprec
beta_gamma_sum(:)=0._rprec

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
  
  u_p = points_avg_3D( u(1:nx,:,1:nz), beta_elem_t(n) % ref_region_t % npoint, points_p ) 
  
  beta_gamma(n) = abs( u_p ) * u_p * area_p
  
  nullify(points_p, area_p, u_p)

enddo

do n=1, nb_elem

  u_p      => b_elem_t(n) % ref_region_t % u
  area_p   => b_elem_t(n) % ref_region_t % area
  points_p => b_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3D( u(1:nx,:,1:nz), b_elem_t(n) % ref_region_t % npoint, points_p ) 
  
  b_gamma(n) = abs( u_p ) * u_p * area_p
  
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

if( temporal_weight == 0 ) then

  if( temporal_model == 1 ) then

    allocate(b_force( nb_elem ))
    b_force(:) = b_r_force(:) - 0.5_rprec * b_elem_t(:) % force_t % CD * beta_gamma_sum(:)
  
    if( spatial_model == 1 ) then
	
	  call b_elem_CD_LE()
  
    elseif( spatial_model == 2) then
	
	  call b_elem_CD_GED()
    
    elseif( spatial_model == 3) then
	
	  call b_elem_CD_GELS()
        
    else
  
      call error( sub_name, 'spatial_model not specified correctly.')

    endif
  
    deallocate(b_force)

  elseif( temporal_model == 2) then ! use implicit formulation

    allocate( b_m( nb_elem ) )  
    b_m(:) = beta_gamma_sum(:) - b_gamma(:)

    if( spatial_model == 1 ) then
	  
	  call b_elem_CD_LI()
 
    elseif( spatial_model == 2 ) then
	
	  call b_elem_CD_GID() ! Global, implicit, direct summation (GID)
  
    elseif( spatial_model == 3 ) then
	
	  call b_elem_CD_GILS()
  
    else
  
      call error( sub_name, 'spatial_model not specified correctly.')

    endif
  
    deallocate( b_m )
  
  else
  
    call error( sub_name, 'temporal_model not specified correctly.')  
 
  endif

  deallocate(b_r_force)

elseif( temporal_weight == 1 ) then

  if( temporal_model == 1 ) then
    
    call error( sub_name, 'temporal_method not specified correctly.')
	
  elseif( temporal_model == 2 ) then
  
    allocate( b_m( nb_elem ) ) 
    b_m(:) = beta_gamma_sum(:) - b_gamma(:)  
	
	if( spatial_model == 1 ) then
	
	  call b_elem_CD_LITW()
	
	elseif( spatial_model == 2 ) then
	
	  call b_elem_CD_GITW()
	
	else
	  
	  call error( sub_name, 'spatial_method not specified correctly.')
	
	endif
	
  else
   
    call error( sub_name, 'temporal_method not specified correctly.')
	
  endif
  
else  

  call error( sub_name, 'temporal_weight not specified correctly.')

endif

!  Check if CD is to be modulated (based on jt_total so can span across multiple runs)
if( jt_total < CD_ramp_nstep ) b_elem_t(:) % force_t % CD = b_elem_t(:) % force_t % CD * jt_total / CD_ramp_nstep

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
beta_elem_t(:) % force_t % fD = - 0.5_rprec * beta_elem_t(:) % force_t % CD * beta_gamma(:)

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
call beta_elem_kappa()
  
deallocate(beta_gamma, beta_gamma_sum)
deallocate(b_gamma)


return

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LE()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the explicit 
!  formulation
!
!  Used variable declarations from contained subroutine rns_elem_force
!
implicit none

b_elem_t(:) % force_t % CD = -2._rprec * b_force(:) / b_gamma(:)   

return
end subroutine b_elem_CD_LE
     

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GED()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the explicit 
!  formulation with direct summation
!
!  Used variable declarations from contained subroutine rns_elem_force
!
implicit none

CD_num=0._rprec
CD_denom=0._rprec
    
do n=1, nb_elem
      
  CD_num = CD_num + b_force(n) 
  CD_denom = CD_denom + b_gamma(n)
    
enddo
    
b_elem_t(:) % force_t % CD = -2._rprec * CD_num / CD_denom    

return
end subroutine b_elem_CD_GED

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GELS()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the explicit 
!  formulation with least squares summation
!
!  Used variable declarations from contained subroutine rns_elem_force
!
implicit none

CD_num=0._rprec
CD_denom=0._rprec
    
do n=1, nb_elem

  CD_num = CD_num + b_force(n) * b_gamma(n) 
  CD_denom = CD_denom + b_gamma(n) * b_gamma(n)
     
enddo
    
b_elem_t(:) % force_t % CD = -2._rprec * CD_num / CD_denom

return
end subroutine b_elem_CD_GELS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LI()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the explicit 
!  formulation
!
!  Used variable declarations from contained subroutine rns_elem_force
!
implicit none

b_elem_t(:) % force_t % CD = 2._rprec * b_r_force(:) / b_m(:)

return
end subroutine b_elem_CD_LI

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GID()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the implicit 
!  formulation with direct summation
!
!  Used variable declarations from contained subroutine rns_elem_force
!
implicit none

CD_num=0._rprec
CD_denom=0._rprec
    
do n=1, nb_elem

  CD_num = CD_num + b_r_force(n) 
  CD_denom = CD_denom + b_m(n)
  
enddo

b_elem_t(1) % force_t % CD = 2._rprec * CD_num / CD_denom 

b_elem_t(:) % force_t % CD = b_elem_t(1) % force_t % CD

return
end subroutine b_elem_CD_GID

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GILS()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the implicit 
!  formulation with direct summation
!
!  Used variable declarations from contained subroutine rns_elem_force
!
implicit none

CD_num=0._rprec
CD_denom=0._rprec
    
do n=1, nb_elem

  CD_num = CD_num + b_r_force(n) * b_m(n)
  CD_denom = CD_denom + b_m(n) * b_m(n)
    
enddo

b_elem_t(:) % force_t % CD = 2._rprec * CD_num / CD_denom

return
end subroutine b_elem_CD_GILS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LITW()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the implicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force
!
use param, only : wbase
implicit none

real(rprec) :: b_r_fsum, b_m_wsum, b_m_wsum2, lambda

real(rprec), pointer :: LRM_p, LMM_p, CD_p
real(rprec), parameter :: thresh = 1.e-4_rprec

integer :: mu_iter, iter
real(rprec) :: b_m_msum, mu_rms, rms
real(rprec), allocatable, dimension(:) :: CD_f, mu, mu_f
!real(rprec) :: b_m_psum, b_m_qsum, denom, 
!

nullify(LRM_p, LMM_p, CD_p)

!  Get LRM and LMM for all b_elem
do n=1,nb_elem
  
  LRM_p    => b_elem_t(n) % force_t % LRM
  LMM_p    => b_elem_t(n) % force_t % LMM
  
  call Lsim( b_r_force(n) * b_m(n), LRM_p )
  call Lsim( b_m(n) * b_m(n), LMM_p )
  
  nullify( LRM_p, LMM_p )
  
enddo

if( jt < weight_nstart ) then

  call b_elem_CD_GID()
  
else

  b_r_fsum  = sum( b_r_force(:) )
  b_m_wsum  = sum( b_elem_t(:) % force_t % LRM * b_m(:) / b_elem_t(:) % force_t % LMM )
  b_m_wsum2 = sum( b_m(:) * b_m(:) / b_elem_t(:) % force_t % LMM ) 

  !  Compute the Lagrange multiplier
  lambda = 2._rprec * ( b_r_fsum  - b_m_wsum ) / b_m_wsum2 
  
  !  Compute CD
  do n = 1, nb_elem
  
    LRM_p => b_elem_t(n) % force_t % LRM
    LMM_p => b_elem_t(n) % force_t % LMM
    CD_p  => b_elem_t(n) % force_t % CD
	
    CD_p = (2._rprec * LRM_p + lambda * b_m(n) ) / LMM_p

    nullify( LRM_p, LMM_p, CD_p )   

  enddo
 

  if( minval( b_elem_t(:) % force_t % CD ) < 0._rprec ) then

    allocate(CD_f(nb_elem))
	allocate(mu(nb_elem))
	allocate(mu_f(nb_elem))

    iter=0
    rms=1.

    do while ( rms > thresh )
	  
	  iter=iter+1
	  CD_f(:) = b_elem_t(:) % force_t % CD
	  
	  mu_rms = 1.
	  mu_iter = 0.
	  mu(:) = 0._rprec
	  do while ( mu_rms > thresh )
	  
	    mu_iter = mu_iter + 1
		mu_f = mu
		b_m_msum = 0._rprec
	    !  Initalize mu
	    do n=1, nb_elem
		  if( b_elem_t(n) % force_t % CD  < 0._rprec ) then
		    mu(n) = b_elem_t(n) % force_t % LRM + 0.5_rprec * lambda * b_m(n)
		  else
		    mu(n) = 0._rprec
		  endif
		  b_m_msum = b_m_msum + mu(n) * b_m(n) / b_elem_t(n) % force_t % LMM
		enddo
		
		lambda = 2._rprec * ( b_r_fsum  - b_m_wsum  - b_m_msum ) / b_m_wsum2 
		
		mu_rms = sqrt( sum( ( mu(:) - mu_f(:) )**2 ) )
	  enddo
	  if(coord == 0) write(*,*) 'mu_iter : ', mu_iter
	  
	  !  Compute new CD
	  b_m_msum = 0._rprec
	  do n=1, nb_elem
        b_m_msum = b_m_msum + mu(n) * b_m(n) / b_elem_t(n) % force_t % LMM
	  enddo		

      b_elem_t(:) % force_t % CD = ( 2._rprec * b_elem_t(:) % force_t % LRM + &
	    lambda * b_m(:) - 2._rprec * mu(:) ) / b_elem_t(:) % force_t % LMM
		
      rms = sqrt( sum( ( b_elem_t(:) % force_t % CD - CD_f(:) )**2 ) )
	  
      if(iter == 1 .and. coord == 0) then

         write(*,'(1a,i6)') 'iter : ', iter
         write(*,'(1a,6f9.4)') 'CD : ', b_elem_t(:) % force_t % CD
         write(*,'(1a,6f9.4)') 'LRM : ', b_elem_t(:) % force_t % LRM
         write(*,'(1a,6f9.4)') 'LMM : ', b_elem_t(:) % force_t % LMM
         write(*,'(1a,6f9.4)') 'b_r_force : ', b_r_force(:)
         write(*,'(1a,6f9.4)') 'b_m : ', b_m(:)

      endif	  
	  
	enddo
	if(coord == 0) write(*,*) 'iter : ', iter
	
	deallocate(CD_f, mu, mu_f )
	
  endif
  
  if(modulo(jt,wbase)==0 .and. coord == 0) then
    write(*,*) '--> Computing LITW CD'
    !write(*,*) '--> lambda : ', lambda
  endif  

endif



return
end subroutine b_elem_CD_LITW

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GITW()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the implicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force
!
use param, only : wbase
implicit none

real(rprec) :: b_r_fsum, b_m_sum
real(rprec) :: lambda
real(rprec), pointer :: LRM_p, LMM_p, CD_p

!integer :: iter
!real(rprec) ::  rms, CD_f, sigma
!real(rprec), parameter :: sigmult = 1.0_rprec
!real(rprec), parameter :: thresh = 1.0e-6_rprec

nullify(LRM_p, LMM_p, CD_p)

LRM_p => b_elem_t(1) % force_t % LRM
LMM_p => b_elem_t(1) % force_t % LMM
CD_p  => b_elem_t(1) % force_t % CD

b_r_fsum = sum( b_r_force(:) )
b_m_sum = sum( b_m(:) )

call Lsim( b_r_fsum * b_m_sum, LRM_p )
call Lsim( b_m_sum * b_m_sum, LMM_p )

!  Update all b elements
b_elem_t(:) % force_t % LRM = LRM_p
b_elem_t(:) % force_t % LMM = LMM_p

if( jt < weight_nstart ) then

  call b_elem_CD_GID()
    
else

  !!  Compute the Lagrange multiplier
  !lambda = 2._rprec * ( b_r_fsum  - LRM_p / LMM_p * b_m_sum ) / ( b_m_sum * b_m_sum / LMM_p ) 

  !  Compute CD
  !CD_p = ( 2._rprec *  LRM_p + lambda * b_m_sum ) / LMM_p
  CD_p = 2._rprec *  LRM_p / LMM_p
 
  !if( CD_p < 0._rprec ) then

  !  sigma = 0.25 * (sum( b_elem_t(:) % force_t % LMM ) / nb_elem) / sigmult
  !  iter=0
  !  rms=1.

  !  do while ( rms > thresh )

  !    CD_f = CD_p
  
  !    sigma = sigmult * sigma
  !    iter = iter + 1
      
      !LRM_p => b_elem_t(1) % force_t % LRM
      !LMM_p => b_elem_t(1) % force_t % LMM
      !CD_p  => b_elem_t(1) % force_t % CD
	
   !   CD_p = (2._rprec * LRM_p - 2._rprec * sigma * minval((/ 0._rprec,  2._rprec * CD_f /))) / LMM_p

   !   rms = abs( CD_p - CD_f )

   !   if(coord == 0) write(*,*) iter, CD_p
    
   ! enddo
 
  !endif 
  
  !  Update all b elements
  b_elem_t(:) % force_t % CD = CD_p

  if(modulo(jt,wbase)==0 .and. coord == 0) then
    write(*,*) '--> Computing GITW CD'
    !write(*,*) '--> lambda : ', lambda
  endif 

endif

nullify( LRM_p, LMM_p, CD_p )

return
end subroutine b_elem_CD_GITW

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine beta_elem_kappa()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the implicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force
!
use param, only : wbase
implicit none

real(rprec) :: beta_int

$if($MPI)
real(rprec) :: beta_int_global
$endif

real(rprec), pointer :: kappa_p, CD_p

nullify(i,j,k)
nullify(npoint_p)
nullify(kappa_p, CD_p)

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
    
    beta_int = beta_int + abs( u(i,j,k) ) * u(i,j,k) * chi(i,j,k) 
 
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
    
  if(coord == 0 .and. (modulo (jt, wbase) == 0)) write(*,'(1a,i3,3f18.6)') 'beta_indx, kappa, CD, beta_int : ', n, kappa_p, CD_p, beta_int
    
  nullify(kappa_p, CD_p)
  
enddo

return
end subroutine beta_elem_kappa

end subroutine rns_elem_force

!**********************************************************************
subroutine Lsim(F, L)
!**********************************************************************
!  This subroutine updates the time weighted function L when using
!  temporal weighted. This advances dL/dt where
!
!  L = \int_{-\infty}^t F * W(t-t') dt'
! 
!  and W(t - t') = 1/Tconst * exp( (t-t')/Tconst ).
!
use types, only : rprec
use param, only : dt

implicit none

real(rprec), intent(IN) :: F
real(rprec), intent(INOUT) :: L

real(rprec) :: Tratio

Tratio = dt / Tconst

L = ( 1._rprec - Tratio ) * L + Tratio * F 

return
end subroutine Lsim

!**********************************************************************
subroutine r_elem_force()
!**********************************************************************
!  This subroutine computes the CD of the all the resolved elements.
!
use types, only : rprec
use messages
use param, only : nx, ny, nz, dx, dy, dz, coord, USE_MPI, jt, wbase
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
  ref_region_t_p % u = points_avg_3D( u(1:nx,:,1:nz), ref_region_t_p % npoint, ref_region_t_p % points)
  
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
  
  if(coord == 0 .and. (modulo (jt, wbase) == 0)) write(*,'(1a,i3,3f18.6)') 'r_indx, fD, CD : ', n, -fD_p, r_elem_t(n) % force_t % CD

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


