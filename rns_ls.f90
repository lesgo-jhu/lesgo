!!
!!  Copyright (C) 2010-2013  Johns Hopkins University
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

!**********************************************************************
module rns_ls
!**********************************************************************
use rns_base_ls

implicit none

save
private

public :: rns_forcing_ls ! Apply forcing
public :: rns_elem_force_ls
public :: rns_finalize_ls
public :: rns_force_init_ls

character (*), parameter :: mod_name = 'rns_ls'

logical :: r_elem_data_init = .false.
logical :: beta_elem_data_init = .false.
logical :: b_elem_data_init = .false.

! Output file id's
integer :: r_elem_cd_fid, &
           r_elem_force_fid, &
           r_elem_vel_fid

integer :: beta_elem_cd_fid, &
           beta_elem_force_fid, &
           beta_elem_kappa_fid, &
           beta_elem_vel_fid

integer :: b_elem_cd_fid, &
           b_elem_error_fid, &
           b_elem_error_norm_fid, &
           b_elem_force_fid, &
           b_elem_vel_fid

contains

!**********************************************************************
subroutine rns_forcing_ls()
!**********************************************************************
!
!  This subroutine computes the forces on the unresolved branches. These
!  forces are based on the velocity field at u^m and not the intermediate
!  velocity u*. The values for kappa is from the m-1 time step
!  as the values from time step m are unknown until after the IBM is applied
!  and rns_elem_force_ls is called.
!
use types, only : rprec
use sim_param, only : u, v
use sim_param, only : fxa, fya
use messages

$if($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP, MPI_SYNC_DOWN
use mpi
use param, only : ld, ny, nz, MPI_RPREC, down, up, comm, status, ierr
$endif

use param, only : dx, dy, dz, coord

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_forcing_ls'

integer :: np, n

integer, pointer :: i, j, k

integer, pointer :: npoint_p
real(rprec) :: uc, vc, cache
real(rprec), pointer :: kappa_p

nullify(i,j,k)
nullify(npoint_p)
nullify(kappa_p)

$if($VERBOSE)
call enter_sub(sub_name)
$endif

!  Apply the RNS forcing to appropriate nodes
do n = 1, nbeta_elem
 
    !  Loop over number of points used in beta region
    npoint_p => beta_elem_t(n) % indx_array_t % npoint    
    kappa_p  => beta_elem_t(n) % force_t % kappa
  
    do np = 1, npoint_p
  
      i => beta_elem_t( n ) % indx_array_t % iarray(1,np)
      j => beta_elem_t( n ) % indx_array_t % iarray(2,np)
      k => beta_elem_t( n ) % indx_array_t % iarray(3,np)
    
      !  Cache values
      uc = u(i,j,k)
      vc = v(i,j,k)
      cache = -kappa_p * sqrt( uc**2 + vc**2 ) * chi(i,j,k)  
      
      !  Modify the RHS to include forces in the evaluation
      !  of the intermediate velocity u*
      fxa(i,j,k) = cache * uc
      fya(i,j,k) = cache * vc

      nullify(i,j,k)
      
    enddo
    
    nullify( npoint_p, kappa_p )
    
enddo

$if($MPI)
! Sync applied forces
call mpi_sync_real_array( fxa, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fya, 1, MPI_SYNC_DOWN )
$endif

$if($VERBOSE)
call exit_sub(sub_name)
$endif

return

end subroutine rns_forcing_ls

!**********************************************************************
subroutine rns_elem_output()
!**********************************************************************
!
use param, only : coord
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_elem_output'

$if($VERBOSE)
call enter_sub(sub_name)
$endif
  
if (coord == 0) then

  call r_elem_data_write()
  call beta_elem_data_write()
  call b_elem_data_write()
      
endif
    
$if($VERBOSE)
call exit_sub(sub_name)
$endif

return
end subroutine rns_elem_output


!**********************************************************************
subroutine rns_elem_force_ls()
!**********************************************************************
!  This subroutine computes the CD and kappa of the beta elements
!
use types, only : rprec
use messages
use sim_param, only : u, v
use sim_param, only : fx, fy
use functions, only : points_avg_3d
use param, only : nx, nz, dx, dy, dz, coord, jt_total
$if($MPI)
use mpi
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_elem_force_ls'

integer ::  n, ns
integer, pointer :: nelem_p
integer, pointer, dimension(:) :: indx_p

real(rprec) :: cache

real(rprec), allocatable, dimension(:,:) :: beta_gamma
real(rprec), allocatable, dimension(:,:) :: b_beta_gamma_sum
real(rprec), allocatable, dimension(:,:) :: b_gamma
real(rprec), allocatable, dimension(:,:) :: b_r_force
real(rprec), allocatable, dimension(:,:) :: b_force, b_m

real(rprec), pointer :: area_p, u_p, v_p
real(rprec), pointer, dimension(:,:) :: points_p

type(force_type_2), pointer :: force_t_p

integer, pointer :: i,j,k
integer, pointer :: npoint_p

$if($VERBOSE)
call enter_sub(sub_name)
$endif

nullify(nelem_p, indx_p)
nullify(area_p, u_p, v_p)
nullify(points_p)
nullify(force_t_p)

allocate(beta_gamma(ndim, nbeta_elem))
allocate(b_gamma(ndim, nb_elem))
allocate(b_beta_gamma_sum(ndim, nb_elem))

beta_gamma=0._rprec
b_gamma=0._rprec
b_beta_gamma_sum=0._rprec

!  Get the force for the resolved elements
call r_elem_force()

do n=1, nbeta_elem

  u_p      => beta_elem_t(n) % ref_region_t % u
  v_p      => beta_elem_t(n) % ref_region_t % v
  area_p   => beta_elem_t(n) % ref_region_t % area
  npoint_p => beta_elem_t(n) % ref_region_t % npoint
  points_p => beta_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3d( u(1:nx,:,1:nz), 1, npoint_p, points_p ) 
  v_p = points_avg_3d( v(1:nx,:,1:nz), 1, npoint_p, points_p )
  
  cache = 0.5_rprec * sqrt( u_p**2 + v_p**2 ) * area_p
  beta_gamma(:,n) = cache * (/ u_p, v_p /)

  nullify(npoint_p, points_p, area_p, u_p, v_p)

enddo

do n=1, nb_elem

  u_p      => b_elem_t(n) % ref_region_t % u
  v_p      => b_elem_t(n) % ref_region_t % v
  area_p   => b_elem_t(n) % ref_region_t % area
  npoint_p => b_elem_t(n) % ref_region_t % npoint
  points_p => b_elem_t(n) % ref_region_t % points
  
  u_p = points_avg_3d( u(1:nx,:,1:nz), 1, npoint_p, points_p ) 
  v_p = points_avg_3d( v(1:nx,:,1:nz), 1, npoint_p, points_p )
  
  cache = 0.5_rprec * sqrt( u_p**2 + v_p**2 ) * area_p
  b_gamma(:,n) = cache * (/ u_p, v_p /)
  
  nullify(npoint_p, points_p, area_p, u_p, v_p)
  
  !  Sum over 
  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx
  b_beta_gamma_sum(:,n) = 0._rprec
  do ns=1, nelem_p
    b_beta_gamma_sum(:,n) = b_beta_gamma_sum(:,n) + beta_gamma(:, indx_p(ns) )
  enddo
  
  nullify(indx_p, nelem_p)
  
enddo

!  Compute the total resolved force of each b_elem
allocate(b_r_force( ndim, nb_elem )) 
b_r_force(:,:) = 0._rprec

do n=1, nb_elem

  nelem_p => b_elem_t(n) % r_child_t % nelem
  indx_p  => b_elem_t(n) % r_child_t % indx
  
  do ns=1, nelem_p
    b_r_force(:,n) = b_r_force(:,n) + (/ r_elem_t( indx_p(ns) )% force_t % fx, &
                                         r_elem_t( indx_p(ns) )% force_t % fy /)
  enddo 
  
  ! Perform partial sum for b_elem force This only sums the resolved
  ! forces. The contribution due to the unresolved (subgrid) forces is
  ! performed after the Cd is calculated.
  b_elem_t(n) % force_t % fx = b_r_force(1,n)
  b_elem_t(n) % force_t % fy = b_r_force(2,n)
  
  nullify(indx_p, nelem_p)

enddo

if( temporal_model == 1 ) then

  !  Compute F_b^n (CD^{n-1})
  allocate(b_force( ndim, nb_elem ))
  do n=1, nb_elem
     b_force(:,n) = b_r_force(:,n) - b_elem_t(n) % force_t % CD * b_beta_gamma_sum(:,n)
  enddo

  if( temporal_weight == 0 ) then

    if( spatial_model == 1 ) then

      call b_elem_CD_LE()
  
    elseif( spatial_model == 2) then

      call b_elem_CD_GE()
        
    else
  
      call error( sub_name, 'spatial_model not specified correctly.')

    endif

  elseif( temporal_weight == 1 ) then 

    if( spatial_model == 1 ) then

      call b_elem_CD_LETW()

    elseif( spatial_model == 2 ) then
  
      call b_elem_CD_GETW()

    else

      call error( sub_name, 'spatial_model not specified correctly.')

    endif 

  else

    call error( sub_name, 'temporal_weight not specified correctly.')

  endif

  !  Check if CD is to be modulated (based on jt_total so can span across multiple runs)
  if( jt_total < CD_ramp_nstep ) b_elem_t(:) % force_t % CD = b_elem_t(:) % force_t % CD * jt_total / CD_ramp_nstep  

  !  Compute the RNS error for the b elem (e_b = F_b + CD_b * gamma_b)
  do n=1, nb_elem
    force_t_p => b_elem_t(n) % force_t
    force_t_p % error = sum( ( b_force(:,n) + force_t_p % CD * b_gamma(:,n) )**2 )
    !force_t_p % error_norm = sum( b_force(:,n)**2 )
    nullify(force_t_p)                            
  enddo

  deallocate(b_force)  

elseif( temporal_model == 2) then ! use implicit formulation

  allocate( b_m( ndim, nb_elem ) )  
  do n=1,nb_elem
    b_m(:,n) = b_gamma(:,n) - b_beta_gamma_sum(:,n)
  enddo  

  if( temporal_weight == 0 ) then

    if( spatial_model == 1 ) then
 
      call b_elem_CD_LI()
 
    elseif( spatial_model == 2 ) then

      call b_elem_CD_GI() ! Global implicit
  
    else
  
      call error( sub_name, 'spatial_model not specified correctly.')

    endif

  elseif( temporal_weight == 1 ) then

    if( spatial_model == 1 ) then

      call b_elem_CD_LITW()

    elseif( spatial_model == 2 ) then
  
      call b_elem_CD_GITW()

    else

      call error( sub_name, 'spatial_model not specified correctly.')

    endif

  else

    call error( sub_name, 'temporal_weight not specified correctly.')

  endif    

  !  Check if CD is to be modulated (based on jt_total so can span across multiple runs)
  if( jt_total < CD_ramp_nstep ) b_elem_t(:) % force_t % CD = b_elem_t(:) % force_t % CD * jt_total / CD_ramp_nstep  

  !  Compute the RNS error for the b elem (e_b = R_b + CD_b * M_b)
  do n=1, nb_elem
    force_t_p => b_elem_t(n) % force_t
    force_t_p % error = sum(( b_r_force(:,n) + force_t_p % CD * b_m(:,n))**2 )
    !force_t_p % error_norm = sum( b_r_force(:,n)**2 )
    nullify(force_t_p)
  enddo  

  deallocate( b_m )

else
  
  call error( sub_name, 'temporal_model not specified correctly.')  
 
endif  

!  Now update the CD of all the beta_elem 
do n=1, nb_elem

  !  Check if b_elem CD < 0
  !if( b_elem_t(n) % force_t % CD < 0._rprec ) b_elem_t(n) % force_t % CD = 0._rprec

  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx
  
  do ns=1, nelem_p
    beta_elem_t( indx_p(ns) ) % force_t % CD = b_elem_t(n) % force_t % CD
  enddo 

  nullify(indx_p, nelem_p)

enddo

!  Compute the total force for beta_elem
do n=1, nbeta_elem
  cache = - beta_elem_t(n) % force_t % CD
  beta_elem_t(n) % force_t % fx = cache * beta_gamma(1,n)
  beta_elem_t(n) % force_t % fy = cache * beta_gamma(2,n)
enddo

!  Compute the total force for b_elem; r_elem has already been accounted for from above
do n=1, nb_elem
  
  nelem_p => b_elem_t(n) % beta_child_t % nelem
  indx_p  => b_elem_t(n) % beta_child_t % indx

  force_t_p => b_elem_t(n) % force_t

  !  Perform secondary sum over beta_elem for b_elem force and the RNS error
  do ns = 1, nelem_p
    force_t_p % fx = force_t_p % fx + beta_elem_t( indx_p(ns) ) % force_t % fx
    force_t_p % fy = force_t_p % fy + beta_elem_t( indx_p(ns) ) % force_t % fy
  enddo 

  ! Compute the normalization factor for the RNS error
  force_t_p % error_norm = force_t_p % fx**2 + force_t_p % fy**2

  nullify( nelem_p, indx_p, force_t_p )
  
enddo

!  Now need to compute kappa; each beta region gets its own kappa value
call beta_elem_kappa()
  
deallocate(beta_gamma, b_beta_gamma_sum)
deallocate(b_gamma)

if(modulo (jt_total, output_nskip) == 0) call rns_elem_output()

$if($VERBOSE)
call exit_sub(sub_name)
$endif

return

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LE()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the explicit 
!  formulation
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

do n=1, nb_elem
  b_elem_t(n) % force_t % CD = - sum( b_force(:,n) * b_gamma(:,n) ) / &
    sum ( b_gamma(:,n) * b_gamma(:,n) )
enddo

if(modulo(jt_total,wbase)==0 .and. coord == 0) then
  write(*,*) '--> Computing LE CD'
endif

return
end subroutine b_elem_CD_LE
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GE()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the explicit 
!  formulation with least squares summation
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec) :: CD_num, CD_denom

CD_num=0._rprec
CD_denom=0._rprec
    
do n=1, nb_elem

  CD_num = CD_num + sum( b_force(:,n) * b_gamma(:,n) )
  CD_denom = CD_denom + sum( b_gamma(:,n) * b_gamma(:,n) )
     
enddo
    
b_elem_t(:) % force_t % CD = - CD_num / CD_denom

if(modulo(jt_total,wbase)==0 .and. coord == 0) then
  write(*,*) '--> Computing GE CD'
endif

return
end subroutine b_elem_CD_GE

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LI()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the implicit
!  formulation
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

do n=1, nb_elem
  b_elem_t(n) % force_t % CD = - sum( b_r_force(:,n) * b_m(:,n) ) / &
    sum( b_m(:,n) * b_m(:,n) )
enddo

if(modulo(jt_total,wbase)==0 .and. coord == 0) then
  write(*,*) '--> Computing LI CD'
endif

return
end subroutine b_elem_CD_LI

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GI()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the implicit 
!  formulation with direct summation
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec) :: CD_num, CD_denom

CD_num=0._rprec
CD_denom=0._rprec
    
do n=1, nb_elem

  CD_num = CD_num + sum( b_r_force(:,n) * b_m(:,n) )
  CD_denom = CD_denom + sum ( b_m(:,n) * b_m(:,n) )
    
enddo

b_elem_t(:) % force_t % CD = - CD_num / CD_denom

if(modulo(jt_total,wbase)==0 .and. coord == 0) then
  write(*,*) '--> Computing GI CD'
endif

return
end subroutine b_elem_CD_GI

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LETW()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the explicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec), dimension(ndim) :: b_gamma_gsum, b_force_gsum

real(rprec), pointer :: LAB_p, LBB_p, CD_p

nullify(LAB_p, LBB_p, CD_p)


!  Get LAB and LBB for all b_elem
do n=1,nb_elem
  
  LAB_p    => b_elem_t(n) % force_t % LAB
  LBB_p    => b_elem_t(n) % force_t % LBB
  
  call Lsim( sum(b_force(:,n) * b_gamma(:,n)), LAB_p )
  call Lsim( sum(b_gamma(:,n) * b_gamma(:,n)), LBB_p )
  
  nullify( LAB_p, LBB_p )
  
enddo

if( jt_total < weight_nstart ) then

  call b_elem_CD_GE()
  
else

  !  Compute CD
  do n = 1, nb_elem
  
    LAB_p => b_elem_t(n) % force_t % LAB
    LBB_p => b_elem_t(n) % force_t % LBB
    CD_p  => b_elem_t(n) % force_t % CD

    CD_p = - LAB_p / LBB_p

    nullify( LAB_p, LBB_p, CD_p )   

  enddo
  
  if(modulo(jt_total,wbase)==0 .and. coord == 0) then
    write(*,*) '--> Computing LETW CD'
  endif  

endif

return
end subroutine b_elem_CD_LETW

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_GETW()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the explicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec) :: b_force_gsum, b_gamma_gsum
real(rprec), pointer :: LAB_p, LBB_p, CD_p

nullify(LAB_p, LBB_p, CD_p)

LAB_p => b_elem_t(1) % force_t % LAB
LBB_p => b_elem_t(1) % force_t % LBB
CD_p  => b_elem_t(1) % force_t % CD

b_force_gsum = 0._rprec
b_gamma_gsum = 0._rprec

!  Least squares summation
do n=1, nb_elem
  b_force_gsum = b_force_gsum + sum( b_force(:,n) * b_gamma(:,n) )
  b_gamma_gsum = b_gamma_gsum + sum( b_gamma(:,n) * b_gamma(:,n) )
enddo

call Lsim( b_force_gsum, LAB_p )
call Lsim( b_gamma_gsum, LBB_p )

!  Update all b elements
b_elem_t(:) % force_t % LAB = LAB_p
b_elem_t(:) % force_t % LBB = LBB_p

if( jt_total < weight_nstart ) then

  call b_elem_CD_GE()
    
else

  !  Compute CD
  CD_p = - LAB_p / LBB_p
 
  !  Update all b elements
  b_elem_t(:) % force_t % CD = CD_p

  if(modulo(jt_total,wbase)==0 .and. coord == 0) then
    write(*,*) '--> Computing GETW CD'
  endif 

endif

nullify( LAB_p, LBB_p, CD_p )

return
end subroutine b_elem_CD_GETW

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine b_elem_CD_LITW()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the local b_elem CD using the implicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec), dimension(ndim) :: b_r_fsum, b_m_wsum, b_m_wsum2
real(rprec), pointer :: LAB_p, LBB_p, CD_p

nullify(LAB_p, LBB_p, CD_p)

!  Get LAB and LBB for all b_elem
do n=1,nb_elem
  
  LAB_p    => b_elem_t(n) % force_t % LAB
  LBB_p    => b_elem_t(n) % force_t % LBB
  
  call Lsim( sum(b_r_force(:,n) * b_m(:,n)), LAB_p )
  call Lsim( sum(b_m(:,n) * b_m(:,n)), LBB_p )
  
  nullify( LAB_p, LBB_p )
  
enddo

if( jt_total < weight_nstart ) then

  call b_elem_CD_GI()
  
else

  !  Compute CD
  do n = 1, nb_elem
  
    LAB_p => b_elem_t(n) % force_t % LAB
    LBB_p => b_elem_t(n) % force_t % LBB
    CD_p  => b_elem_t(n) % force_t % CD

    CD_p = - LAB_p / LBB_p

    nullify( LAB_p, LBB_p, CD_p )   

  enddo
  
  if(modulo(jt_total,wbase)==0 .and. coord == 0) then
    write(*,*) '--> Computing LITW CD'
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
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec) :: b_r_msum, b_m_msum
real(rprec), pointer :: LAB_p, LBB_p, CD_p

nullify(LAB_p, LBB_p, CD_p)

LAB_p => b_elem_t(1) % force_t % LAB
LBB_p => b_elem_t(1) % force_t % LBB
CD_p  => b_elem_t(1) % force_t % CD

b_r_msum = 0._rprec
b_m_msum = 0._rprec

!  Least squares summation
do n=1, nb_elem
  b_r_msum = b_r_msum + sum( b_r_force(:,n) * b_m(:,n) )
  b_m_msum = b_m_msum + sum( b_m(:,n) * b_m(:,n) )
enddo

call Lsim( b_r_msum, LAB_p )
call Lsim( b_m_msum, LBB_p )

!  Update all b elements
b_elem_t(:) % force_t % LAB = LAB_p
b_elem_t(:) % force_t % LBB = LBB_p

if( jt_total < weight_nstart ) then

  call b_elem_CD_GI()
    
else

  !  Compute CD
  CD_p = - LAB_p / LBB_p
 
  !  Update all b elements
  b_elem_t(:) % force_t % CD = CD_p

  if(modulo(jt_total,wbase)==0 .and. coord == 0) then
    write(*,*) '--> Computing GITW CD'
  endif 

endif

nullify( LAB_p, LBB_p, CD_p )

return
end subroutine b_elem_CD_GITW

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine beta_elem_kappa()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine computes the global b_elem CD using the implicit 
!  formulation with temporal weighting
!
!  Used variable declarations from contained subroutine rns_elem_force_ls
!
use param, only : wbase
implicit none

real(rprec) :: uc, vc
real(rprec), dimension(ndim) :: beta_int

$if($MPI)
real(rprec), dimension(ndim) :: beta_int_global
$endif

real(rprec), pointer :: kappa_p, CD_p

nullify(i,j,k)
nullify(npoint_p)
nullify(kappa_p, CD_p)

do n = 1, nbeta_elem

  !  Compute beta_int over each region beta
  beta_int(:) = 0._rprec
    
  $if($MPI)
  beta_int_global(:) = 0._rprec
  $endif
  
  !  Loop over number of points used in beta region
  npoint_p => beta_elem_t(n) % indx_array_t % npoint
  do ns = 1, npoint_p
  
    i => beta_elem_t(n) % indx_array_t % iarray(1,ns)
    j => beta_elem_t(n) % indx_array_t % iarray(2,ns)
    k => beta_elem_t(n) % indx_array_t % iarray(3,ns)
   
    uc = u(i,j,k)
    vc = v(i,j,k) 
    cache = sqrt( uc**2 + vc**2 ) * chi(i,j,k)
    beta_int(:) = beta_int(:) + cache * (/ uc, vc /)  
 
    nullify(i,j,k)
      
  enddo
  
  beta_int(:) = beta_int(:) * dx * dy * dz
    
  nullify( npoint_p )
    
  $if($MPI)
  call mpi_allreduce (beta_int(1), beta_int_global(1), 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce (beta_int(2), beta_int_global(2), 1, MPI_RPREC, MPI_SUM, comm, ierr)
  beta_int(:) = beta_int_global(:)
  $endif
    
  kappa_p => beta_elem_t(n) % force_t % kappa
  CD_p    => beta_elem_t(n) % force_t % CD
    
  kappa_p = CD_p * sum( beta_gamma(:,n) * beta_int(:) ) / ( sum( beta_int(:) * beta_int(:) ) )
    
  if(coord == 0 .and. (modulo (jt_total, wbase) == 0)) write(*,'(1a,(1x,i3),4(1x,f9.4))') 'beta_indx, kappa, CD, beta_int : ', n, kappa_p, CD_p, beta_int
    
  nullify(kappa_p, CD_p)
  
enddo

return
end subroutine beta_elem_kappa

end subroutine rns_elem_force_ls

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
use param, only : nx, ny, nz, dx, dy, dz, coord, jt_total, wbase
$if($MPI)
use param, only : MPI_RPREC, MPI_SUM, comm, ierr
$endif
use sim_param, only : u, v
use functions, only : points_avg_3d
use sim_param, only : fx, fy
implicit none

character (*), parameter :: sub_name = mod_name // '.r_elem_force'

integer :: n, np

integer, pointer :: i, j, k
integer, pointer :: npoint_p
integer, pointer, dimension(:,:) :: iarray_p

real(rprec) :: cache
real(rprec), pointer :: fx_p, fy_p

$if ($MPI)
real(rprec) :: fx_l, fy_l
$endif

type(ref_region), pointer :: ref_region_t_p
type(indx_array), pointer :: indx_array_t_p

$if($VERBOSE)
call enter_sub(sub_name)
$endif

!  Comment starts here 
nullify(ref_region_t_p)
nullify(indx_array_t_p)
nullify(npoint_p, iarray_p)
nullify(i,j,k)
nullify(fx_p, fy_p)

do n = 1, nr_elem

  !  Get the reference velocity
  ref_region_t_p => r_elem_t( n ) % ref_region_t
  ref_region_t_p % u = points_avg_3d( u(1:nx,:,1:nz), 1, ref_region_t_p % npoint, ref_region_t_p % points)
  ref_region_t_p % v = points_avg_3d( v(1:nx,:,1:nz), 1, ref_region_t_p % npoint, ref_region_t_p % points)
  
  indx_array_t_p => r_elem_t( n ) % indx_array_t
     
  npoint_p => indx_array_t_p % npoint
  iarray_p => indx_array_t_p % iarray
  
  $if($MPI)
  fx_l = 0._rprec
  fy_l = 0._rprec
  $endif

  fx_p => r_elem_t( n ) % force_t % fx
  fy_p => r_elem_t( n ) % force_t % fy
  fx_p = 0._rprec
  fy_p = 0._rprec
  
  do np=1, npoint_p
  
    i => iarray_p(1,np)
    j => iarray_p(2,np)
    k => iarray_p(3,np)
    
    if( k == nz ) call error( sub_name, 'Summing over bogus fx')
  
    $if($MPI)
    fx_l = fx_l + fx(i,j,k)
    fy_l = fy_l + fy(i,j,k)
    $else
    fx_p = fx_p + fx(i,j,k)
    fy_p = fy_p + fy(i,j,k)
    $endif
    
    nullify(i,j,k)
    
  enddo
  
  $if($MPI)
  call mpi_allreduce (fx_l, fx_p, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce (fy_l, fy_p, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  $endif

  cache = dx * dy * dz
  fx_p = fx_p * cache
  fy_p = fy_p * cache
   
  !  Compute CD
  r_elem_t(n) % force_t % CD = -(fx_p * ref_region_t_p % u + fy_p * ref_region_t_p % v)/ &
    (0.5_rprec * ref_region_t_p % area * sqrt((ref_region_t_p % u)**2 + (ref_region_t_p % v)**2)**3)
  
  if(coord == 0 .and. (modulo (jt_total, wbase) == 0)) write(*,'(1a,(1x,i3),4(1x,f9.4))') 'r_indx, fx, fy, CD : ', n, -fx_p, -fy_p, r_elem_t(n) % force_t % CD

  nullify(fx_p, fy_p)
  nullify(npoint_p, iarray_p)
  nullify(indx_array_t_p)
  nullify(ref_region_t_p)
 
enddo

$if($VERBOSE)
call exit_sub(sub_name)
$endif

return
end subroutine r_elem_force

!**********************************************************************
subroutine r_elem_data_write()
!**********************************************************************
use param, only : total_time, path
use string_util
use open_file_fid_mod

implicit none

include 'tecryte.h'

character(*), parameter :: sub_name = mod_name // '.r_elem_data_write'
character(*), parameter :: fname_CD = path // 'output/rns/r_elem_cd.dat'
character(*), parameter :: fname_force = path // 'output/rns/r_elem_force.dat'
character(*), parameter :: fname_vel = path // 'output/rns/r_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n

if( .not. r_elem_data_init ) then

  inquire (file=fname_CD, exist=exst)
  r_elem_cd_fid = open_file_fid( fname_cd, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nr_elem
      !  Create variable list name:
      call string_concat(var_list, ',"CD<sub>',n,'</sub>"')
    enddo
    call write_tecplot_header_xyline( r_elem_cd_fid, trim(adjustl(var_list)) )
  endif

  inquire (file=fname_force, exist=exst)
  r_elem_force_fid = open_file_fid( fname_force, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nr_elem
      !  Create variable list name:
      call string_concat(var_list, ',"fx<sub>',n,'</sub>"')
    enddo
    ! Write the total fx name
    call string_concat(var_list, ',"fx_tot"')

    do n = 1, nr_elem
      !  Create variable list name:
      call string_concat(var_list, ',"fy<sub>',n,'</sub>')
    enddo  
    ! Write the total fy name
    call string_concat(var_list, ',"fy_tot"')
    call write_tecplot_header_xyline( r_elem_force_fid, trim(adjustl(var_list)) )

  endif

  inquire (file=fname_vel, exist=exst)
  r_elem_vel_fid = open_file_fid( fname_vel, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nr_elem
      !  Create variable list name:
      call string_concat(var_list, ',"u<sub>',n,'</sub>"')
    enddo
    do n = 1, nr_elem
      !  Create variable list name:
      call string_concat(var_list, ',"v<sub>',n,'</sub>"')
    enddo

    call write_tecplot_header_xyline( r_elem_vel_fid, trim(adjustl(var_list)) )

  endif

  r_elem_data_init = .true.

endif  

call write_real_data(r_elem_cd_fid, 'formatted', nr_elem + 1, &
     (/ total_time, r_elem_t(:) % force_t % CD /))

call write_real_data(r_elem_force_fid, 'formatted', ndim*(nr_elem+1)+1, &
     (/ total_time, -r_elem_t(:) % force_t % fx, &
     -sum( r_elem_t(:) % force_t % fx ), &
     -r_elem_t(:) % force_t % fy, &
     -sum( r_elem_t(:) % force_t % fy ) /))

call write_real_data(r_elem_vel_fid, 'formatted', ndim*nr_elem+1, &
     (/ total_time, r_elem_t(:) % ref_region_t % u, &
     r_elem_t(:) % ref_region_t % v /))

return
end subroutine r_elem_data_write

!**********************************************************************
subroutine beta_elem_data_write()
!**********************************************************************
use param, only : total_time, path
use string_util
use open_file_fid_mod

implicit none

include 'tecryte.h'

character(*), parameter :: sub_name = mod_name // '.beta_elem_data_write'
character(*), parameter :: fname_CD = path // 'output/rns/beta_elem_cd.dat'
character(*), parameter :: fname_force = path // 'output/rns/beta_elem_force.dat'
character(*), parameter :: fname_kappa = path // 'output/rns/beta_elem_kappa.dat'
character(*), parameter :: fname_vel = path // 'output/rns/beta_elem_vel.dat'


logical :: exst
character(5000) :: var_list
integer :: n


if( .not. beta_elem_data_init ) then

   ! Initialize files and headers
   inquire (file=fname_CD, exist=exst)
   beta_elem_cd_fid = open_file_fid( fname_cd, 'append', 'formatted' )

   if (.not. exst) then
      var_list = '"t"'
      do n = 1, nbeta_elem
         !  Create variable list name:
         call string_concat(var_list, ',"CD<sub>',n,'</sub>"')
      enddo

      call write_tecplot_header_xyline( beta_elem_cd_fid, trim(adjustl(var_list)) )

   endif

   inquire (file=fname_force, exist=exst)
   beta_elem_force_fid = open_file_fid( fname_force, 'append', 'formatted' )

   if (.not. exst) then
      var_list = '"t"'
      do n = 1, nbeta_elem
         !  Create variable list name:
         call string_concat(var_list, ',"fx<sub>',n,'</sub>"')
      enddo
      ! Write the total fx name
      call string_concat(var_list, ',"fx_tot"')

      do n = 1, nbeta_elem
         !  Create variable list name:
         call string_concat(var_list, ',"fy<sub>',n,'</sub>"')
      enddo
      ! Write the total fy name
      call string_concat(var_list, ',"fy_tot"')

      call write_tecplot_header_xyline( beta_elem_force_fid, trim(adjustl(var_list)) )

   endif

   inquire (file=fname_kappa, exist=exst)
   beta_elem_kappa_fid = open_file_fid( fname_kappa, 'append', 'formatted' )

   if (.not. exst) then
      var_list = '"t"'
      do n = 1, nbeta_elem
         !  Create variable list name:
         call string_concat(var_list, ',"<greek>k</greek><sub>',n,'</sub>"')
      enddo
      call write_tecplot_header_xyline( beta_elem_kappa_fid, trim(adjustl(var_list)) )

   endif

   inquire (file=fname_vel, exist=exst)
   beta_elem_vel_fid = open_file_fid( fname_vel, 'append', 'formatted' )

   if (.not. exst) then
      var_list = '"t"'
      do n = 1, nbeta_elem
         !  Create variable list name:
         call string_concat(var_list, ',"u<sub>',n,'</sub>"')
      enddo
      do n = 1, nbeta_elem
         !  Create variable list name:
         call string_concat(var_list, ',"v<sub>',n,'</sub>"')
      enddo
      call write_tecplot_header_xyline(beta_elem_vel_fid, trim(adjustl(var_list)) )

   endif

   beta_elem_data_init = .true. 

endif  

! Write all the data to file
call write_real_data(beta_elem_cd_fid, 'formatted', nbeta_elem + 1, &
     (/ total_time, beta_elem_t(:) % force_t % CD /))

call write_real_data(beta_elem_force_fid, 'formatted', ndim*(nbeta_elem + 1) + 1, &
     (/ total_time, -beta_elem_t(:) % force_t % fx, &
     -sum(beta_elem_t(:) % force_t % fx), &
     -beta_elem_t(:) % force_t % fy, &
     -sum(beta_elem_t(:) % force_t % fy) /))

call write_real_data(beta_elem_kappa_fid, 'formatted', nbeta_elem + 1, &
     (/ total_time, beta_elem_t(:) % force_t % kappa /))

call write_real_data(beta_elem_vel_fid, 'formatted', ndim*nbeta_elem+1, &
     (/ total_time, beta_elem_t(:) % ref_region_t % u, &
     beta_elem_t(:) % ref_region_t % v /))

return
end subroutine beta_elem_data_write

!**********************************************************************
subroutine b_elem_data_write()
!**********************************************************************
use param, only : total_time, path
use string_util
use open_file_fid_mod
implicit none

include 'tecryte.h'

character(*), parameter :: sub_name = mod_name // '.b_elem_data_write'
character(*), parameter :: fname_CD = path // 'output/rns/b_elem_cd.dat'
character(*), parameter :: fname_force = path // 'output/rns/b_elem_force.dat'
character(*), parameter :: fname_vel = path // 'output/rns/b_elem_vel.dat'
character(*), parameter :: fname_error = path // 'output/rns/b_elem_error.dat'
character(*), parameter :: fname_error_norm = path // 'output/rns/b_elem_error_norm.dat'

logical :: exst
character(5000) :: var_list
integer :: n

if( .not. b_elem_data_init ) then

  inquire (file=fname_CD, exist=exst)
  b_elem_cd_fid = open_file_fid( fname_cd, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"CD<sub>',n,'</sub>"')
    enddo
    call write_tecplot_header_xyline(b_elem_cd_fid, trim(adjustl(var_list)))

  endif

  inquire (file=fname_force, exist=exst)
  b_elem_force_fid = open_file_fid( fname_force, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"fx<sub>',n,'</sub>"')
    enddo
    call string_concat(var_list, ',"fx_tot"')
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"fy<sub>',n,'</sub>"')
    enddo  
    call string_concat(var_list, ',"fy_tot"')
    call write_tecplot_header_xyline(b_elem_force_fid, trim(adjustl(var_list)))

  endif

  inquire (file=fname_error, exist=exst)
  b_elem_error_fid = open_file_fid( fname_error, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"error<sub>',n,'</sub>"')
    enddo
    call write_tecplot_header_xyline(b_elem_error_fid, trim(adjustl(var_list)))

  endif 

  inquire (file=fname_error_norm, exist=exst)
  b_elem_error_norm_fid = open_file_fid( fname_error_norm, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"error_norm<sub>',n,'</sub>"')
    enddo
    call write_tecplot_header_xyline(b_elem_error_norm_fid, trim(adjustl(var_list)))

  endif  

  inquire (file=fname_vel, exist=exst)
  b_elem_vel_fid = open_file_fid( fname_vel, 'append', 'formatted' )

  if (.not. exst) then
    var_list = '"t"'
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"u<sub>',n,'</sub>"')
    enddo
    do n = 1, nb_elem
      !  Create variable list name:
      call string_concat(var_list, ',"v<sub>',n,'</sub>"')
    enddo  
    call write_tecplot_header_xyline(b_elem_vel_fid, trim(adjustl(var_list)))
    
  endif

  b_elem_data_init = .true. 

endif       


call write_real_data(b_elem_cd_fid, 'formatted', nb_elem + 1, &
     (/ total_time, b_elem_t(:) % force_t % CD /))

call write_real_data(b_elem_force_fid, 'formatted', ndim*(nb_elem+1)+1, &
     (/ total_time, -b_elem_t(:) % force_t % fx, &
     -sum(b_elem_t(:) % force_t % fx), &
     -b_elem_t(:) % force_t % fy, &
     -sum(b_elem_t(:) % force_t % fy) /))

call write_real_data(b_elem_error_fid, 'formatted', nb_elem+1, &
     (/ total_time, b_elem_t(:) % force_t % error /))

call write_real_data(b_elem_error_norm_fid, 'formatted', nb_elem+1, &
                     (/ total_time, b_elem_t(:) % force_t % error_norm /))

call write_real_data(b_elem_vel_fid, 'formatted', ndim*nb_elem+1, &
     (/ total_time, b_elem_t(:) % ref_region_t % u, &
     b_elem_t(:) % ref_region_t % v /))

return
end subroutine b_elem_data_write

!**********************************************************************
subroutine rns_force_init_ls ()
!**********************************************************************
!  
!  This subroutine reads the last BETA force data from a previous simulation
!
use types, only : rprec
use param, only : coord, path
use string_util, only : string_splice
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_force_init'
character (*), parameter :: fname_in = path // 'rns.out'
character (128) :: fname
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

$endif

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
call string_splice( fname, fname_in // MPI_suffix, coord )
$else
fname = fname_in
$endif

inquire (file=fname, exist=exst)

if (.not. exst) then
  if (coord == 0) then
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

read(1) r_elem_t(:) % force_t, beta_elem_t(:) % force_t, b_elem_t(:) % force_t

close(1)

end subroutine rns_force_init_ls

!**********************************************************************
subroutine rns_finalize_ls()
!**********************************************************************
! 
!  This subroutine writes all restart data to file
!
use param, only : coord, path
$if($MPI)
use param, only : comm, ierr
$endif
use messages
use string_util, only : string_splice
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_finalize_ls'
character (*), parameter :: fname_out = path // 'rns.out'

character (128) :: fname
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

$endif

logical :: opn

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
call string_splice( fname, fname_out // MPI_suffix, coord )
$else
fname = fname_out
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

write(1) r_elem_t(:) % force_t, beta_elem_t(:) % force_t, b_elem_t(:) % force_t
close (1)

deallocate(r_elem_t)
deallocate(beta_elem_t)
deallocate(b_elem_t)

! Close all opened files
close( r_elem_cd_fid )
close( r_elem_force_fid )
close( r_elem_vel_fid )

close( beta_elem_cd_fid )
close( beta_elem_force_fid )
close( beta_elem_kappa_fid )
close( beta_elem_vel_fid )

close( b_elem_cd_fid )
close( b_elem_error_fid )
close( b_elem_error_norm_fid )
close( b_elem_force_fid )
close( b_elem_vel_fid )

$if($MPI)
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
$endif

return
end subroutine rns_finalize_ls

end module rns_ls


