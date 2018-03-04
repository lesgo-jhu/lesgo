!!!
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

! this is the w-node version
!--provides Ds_opt2 1:nz
!--MPI: requires u,v 0:nz, except bottom process only 1:nz

module s_lagrange_Sdep

implicit none

contains

subroutine scalars_lagrange_Sdep(scalar,dsdx,dsdy,dsdz,Ds_opt2,s_Tn_all,s_Beta,I_LM,I_MM,I_QN,I_NN,ds2_clips,I_LM_MM_init,I_QN_NN_init,sigma_s)
! Lagrangian scale-dependent dynamic model to calculate the 
! equivalent of the Smagorinsky coefficient for the eddy diffusivity 
! for scalarsthis is done layer-by-layer to save memory.
! note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary.  mem requirement ~ n^2, not n^3
use types,only:rprec
use param
use sim_param,only:u,v,w
!use sgs_param,only: I_LM,I_MM,I_QN,I_NN,s_beta,opftime,s_lagran_dt
use sgs_param,only:S11,S12,S13,S22,S23,S33,delta,S
use sgs_param,only:S_bar,S11_bar,S12_bar,S13_bar,S22_bar,S23_bar,S33_bar
use sgs_param,only:S_S11_bar,S_S12_bar,S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
use sgs_param,only:u_bar,v_bar,w_bar,scalar_bar
use sgs_param,only:u_hat,v_hat,w_hat,scalar_hat
use sgs_param,only:S_hat,S11_hat,S12_hat,S13_hat,S22_hat,S23_hat,S33_hat
use sgs_param,only:S_S11_hat,S_S12_hat,S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat
use sgs_param,only:L1,L2,L3,M1,M2,M3
use sgs_param,only:Q1,Q2,Q3,N1,N2,N3
use sgs_param,only:dsdx_bar, dsdy_bar, dsdz_bar, S_dsdx_bar, S_dsdy_bar, S_dsdz_bar
use sgs_param,only:dsdx_hat, dsdy_hat, dsdz_hat, S_dsdx_hat, S_dsdy_hat, S_dsdz_hat
use sgs_param,only:opftime,s_lagran_dt
!use sgs_param,only:s_Tn_all, ds2_clips
use test_filtermodule
$if ($MPI)
use mpi_defs, only:mpi_sync_real_array,MPI_SYNC_DOWNUP
$endif

implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: scalar
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdx
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdy
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdz
real(rprec), dimension(ld,ny,nz), intent (out) :: Ds_opt2
real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: I_LM, I_MM, I_QN, I_NN
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: s_Tn_all, s_Beta, ds2_clips
logical, intent(inout) :: I_LM_MM_init
logical, intent(inout) :: I_QN_NN_init
!real(rprec), intent(in) :: scalarstar
real(rprec), dimension(lbz:nz), intent(inout) :: sigma_s

integer :: jx,jy,jz

real(rprec):: tf1,tf2,tf1_2,tf2_2 ! Size of the second test filter
real(rprec) :: fractus
real(rprec) :: Betaclip  !--scalar to save mem., otherwise (ld,ny,nz)
real(rprec), dimension(ld,ny) :: Ds_opt2_2d,Ds_opt2_4d

real(rprec), dimension(ld,ny) :: LM,MM,QN,NN,s_Tn,epsi,dumfac

real(rprec) :: const
real(rprec) :: opftdelta,powcoeff

real(rprec), parameter :: zero=1.e-24_rprec ! zero = infimum(0)

!logical, save :: I_LM_MM_init = .false.
!logical, save :: I_QN_NN_init = .false.

integer :: i, j
integer :: count_ds2_all, count_ds2_clip

!---------------------------------------------------------------------
! Set coefficients
opftdelta = opftime*delta
powcoeff = -1._rprec/4._rprec
fractus= 1._rprec/real(ny*nx,kind=rprec)
const = (delta**2) !Note difference from const in eddy viscosity model
tf1=2._rprec
tf2=4._rprec
tf1_2=tf1**2
tf2_2=tf2**2

!Get standard deviation of scalar quantity for the calculation of Tn
call get_scalar_std(scalar,sigma_s)

!if (coord==6) write(*,*) 'A. jt, coord, sigma_s(nz-1)', jt, coord, sigma_s(nz-1)
!if (coord==6) write(*,*) 'A. jt, coord, sigma_s(nz)', jt, coord, sigma_s(nz)
!if (coord==7) write(*,*) 'A. jt, coord, sigma_s(0)', jt, coord, sigma_s(0)
!if (coord==7) write(*,*) 'A. jt, coord, sigma_s(1)', jt, coord, sigma_s(1)

! "Rearrange" I_* (running averages) so that their new positions (i,j,k)
!  correspond to the current (i,j,k) particle 
call scalars_interpolag_Sdep(I_LM,I_MM,I_QN,I_NN)

! For each horizontal level, calculate L_i(:,:), Q_i(:,:), M_i(:,:), and N_i(:,:).
! Then update the running averages, I_*(:,:,jz), which are used to 
! calculate Ds_opt2(:,:,jz).
do jz=1,nz
   count_ds2_all = 0
   count_ds2_clip = 0

   ! Calculate L_i
   !Interp u,v,w,scalar onto w-nodes and store result as u_bar,v_bar,w_bar
   !(except for very first level which should be on uvp-nodes
   if ( (coord == 0) .and. (jz == 1) ) then
      u_bar(:,:) = u(:,:,1)     
      v_bar(:,:) = v(:,:,1)     
      w_bar(:,:) = 0.25_rprec*w(:,:,2)     
      scalar_bar(:,:) = scalar(:,:,1)     
   else !w-nodes
      u_bar(:,:) = 0.5_rprec*( u(:,:,jz) + u(:,:,jz-1) ) 
      v_bar(:,:) = 0.5_rprec*( v(:,:,jz) + v(:,:,jz-1) ) 
      w_bar(:,:) = w(:,:,jz)
      scalar_bar(:,:) = 0.5_rprec*( scalar(:,:,jz) + scalar(:,:,jz-1) ) 
   end if
   u_hat = u_bar
   v_hat = v_bar
   w_hat = w_bar
   scalar_hat = scalar_bar  

   !First term before filtering
   L1 = u_bar*scalar_bar
   L2 = v_bar*scalar_bar
   L3 = w_bar*scalar_bar

   Q1 = L1
   Q2 = L2 
   Q3 = L3

   !Filter first term and add the second term to get the final value
   call test_filter ( u_bar )   ! in-place filtering
   call test_filter ( v_bar )
   call test_filter ( w_bar )
   call test_filter ( scalar_bar )
   call test_filter ( L1 )  
   L1 = L1 - u_bar*scalar_bar  
   call test_filter ( L2 )
   L2 = L2 - v_bar*scalar_bar
   call test_filter ( L3 )
   L3 = L3 - w_bar*scalar_bar
   
   call test_test_filter ( u_hat )
   call test_test_filter ( v_hat )
   call test_test_filter ( w_hat )
   call test_test_filter ( scalar_hat )
   call test_test_filter ( Q1 )
   Q1 = Q1 - u_hat*scalar_hat
   call test_test_filter ( Q2 )
   Q2 = Q2 - v_hat*scalar_hat
   call test_test_filter ( Q3 )
   Q3 = Q3 - w_hat*scalar_hat
   
   !calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +& 
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

   !already on w-nodes
   S11_bar(:,:) = S11(:,:,jz)  
   S12_bar(:,:) = S12(:,:,jz)  
   S13_bar(:,:) = S13(:,:,jz)  
   S22_bar(:,:) = S22(:,:,jz)  
   S23_bar(:,:) = S23(:,:,jz)  
   S33_bar(:,:) = S33(:,:,jz)  
   
   S11_hat = S11_bar
   S12_hat = S12_bar
   S13_hat = S13_bar
   S22_hat = S22_bar
   S23_hat = S23_bar
   S33_hat = S33_bar
   
   call test_filter ( S11_bar )
   call test_filter ( S12_bar )
   call test_filter ( S13_bar )
   call test_filter ( S22_bar )
   call test_filter ( S23_bar )
   call test_filter ( S33_bar )

   call test_test_filter ( S11_hat )
   call test_test_filter ( S12_hat )
   call test_test_filter ( S13_hat )
   call test_test_filter ( S22_hat )
   call test_test_filter ( S23_hat )
   call test_test_filter ( S33_hat )

   !Calculate |S_bar| (the test_filtered Sij)
   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   !Calculate |S_hat| (the test_test_filtered Sij)
   S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
        2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

   !Filter the temperature gradients
   if ( (coord==0).and.(jz==1) ) then !store on uvp-nodes 
      dsdx_bar(:,:) = dsdx(:,:,1)
      dsdy_bar(:,:) = dsdy(:,:,1)
      dsdz_bar(:,:) = dsdz(:,:,2) !Try this
   else
      dsdx_bar(:,:) = 0.5_rprec*( dsdx(:,:,jz) + dsdx(:,:,jz-1) )
      dsdy_bar(:,:) = 0.5_rprec*( dsdy(:,:,jz) + dsdy(:,:,jz-1) )
      dsdz_bar(:,:) = dsdz(:,:,jz) 
   end if 

   dsdx_hat = dsdx_bar
   dsdy_hat = dsdy_bar
   dsdz_hat = dsdz_bar

   !Calculate |S|ds_i
   S_dsdx_bar = S * dsdx_bar
   S_dsdy_bar = S * dsdy_bar
   S_dsdz_bar = S * dsdz_bar

   S_dsdx_hat = S_dsdx_bar
   S_dsdy_hat = S_dsdy_bar
   S_dsdz_hat = S_dsdz_bar

   !Test filter temperature gradients
   call test_filter ( dsdx_bar ) 
   call test_filter ( dsdy_bar ) 
   call test_filter ( dsdz_bar ) 

   call test_test_filter ( dsdx_hat ) 
   call test_test_filter ( dsdy_hat ) 
   call test_test_filter ( dsdz_hat ) 

   call test_filter ( S_dsdx_bar ) 
   call test_filter ( S_dsdy_bar ) 
   call test_filter ( S_dsdz_bar ) 

   call test_test_filter ( S_dsdx_hat ) 
   call test_test_filter ( S_dsdy_hat ) 
   call test_test_filter ( S_dsdz_hat ) 

   !Calculate M_i and N_i
   M1 = const*(S_dsdx_bar - tf1_2*S_bar*dsdx_bar)
   M2 = const*(S_dsdy_bar - tf1_2*S_bar*dsdy_bar)
   M3 = const*(S_dsdz_bar - tf1_2*S_bar*dsdz_bar)

   N1 = const*(S_dsdx_hat - tf2_2*S_hat*dsdx_hat)
   N2 = const*(S_dsdy_hat - tf2_2*S_hat*dsdy_hat)
   N3 = const*(S_dsdz_hat - tf2_2*S_hat*dsdz_hat)

   !Calculate L_i*M_i, M_i*M_i, Q_i*N_i, and N_i*N_i for each 
   !point in the plane 
   LM = L1*M1 + L2*M2 + L3*M3
   MM = M1*M1 + M2*M2 + M3*M3
   QN = Q1*N1 + Q2*N2 + Q3*N3
   NN = N1*N1 + N2*N2 + N3*N3

   !Initialize sgs quantities
   if (inilag_scalar) then 
      if ((.not.I_LM_MM_init).and.(jt==cs_count.or.jt==ssgs_init)) then
         !print*,'I_LM and I_MM initialized'
         I_MM (:,:,jz) = MM
         I_LM (:,:,jz) = 0.03_rprec * MM
         I_MM (ld-1:ld,:,jz) = 1._rprec
         I_LM (ld-1:ld,:,jz) = 1._rprec
   
         if (jz==nz) I_LM_MM_init=.true.
      end if
   end if

   !Update running averages (I_LM, I_MM)
   !Determine averaging timescale (for 2-delta filter)
   s_Tn = max( I_LM(:,:,jz)*I_MM(:,:,jz),zero )
   !s_Tn = opftdelta*(s_Tn**powcoeff)
   s_Tn = opftdelta*sigma_s(jz)*(s_Tn**powcoeff)
   !Clip if necessary
   s_Tn(:,:) = max( zero, s_Tn(:,:) )

   !Calculate new running average 
   dumfac = s_lagran_dt/s_Tn
   epsi = dumfac / (1._rprec+dumfac)
 
   I_LM(:,:,jz) = epsi*LM + (1._rprec-epsi)*I_LM(:,:,jz)
   I_MM(:,:,jz) = epsi*MM + (1._rprec-epsi)*I_MM(:,:,jz)
   !Clip if necessary
   I_LM(:,:,jz) = max( zero, I_LM(:,:,jz) )

   !Calculate Ds_opt2 (for 2-delta filter)
   ! Add +zero in denominator to avoid division by identically zero
   Ds_opt2_2d(:,:) = I_LM(:,:,jz)/(I_MM(:,:,jz)+zero)
   Ds_opt2_2d(ld,:) = zero 
   Ds_opt2_2d(ld-1,:) = zero 
   !Clip if necessary 
   Ds_opt2_2d(:,:) = max( zero, Ds_opt2_2d(:,:) ) 

   !Initialize sgs quantities
   if (inilag_scalar) then
      if ((.not.I_QN_NN_init).and.(jt==cs_count.or.jt==ssgs_init)) then
         !print*, 'I_NN and I_QN initialized'
         I_NN(:,:,jz) = NN
         I_QN(:,:,jz) = 0.03_rprec*NN 
         I_NN(ld-1:ld,:,jz) = 1._rprec
         I_QN(ld-1:ld,:,jz) = 1._rprec

         if (jz==nz) I_QN_NN_init=.true.
      end if
   end if

   !Update running averages 
   s_Tn = max( I_QN(:,:,jz)*I_NN(:,:,jz), zero )
   !s_Tn = opftdelta*(s_Tn**powcoeff)
   s_Tn = opftdelta*sigma_s(jz)*(s_Tn**powcoeff)
   !Clip, if necessary
   s_Tn(:,:) = max( zero, s_Tn(:,:) )

   !Calculate new running average
   dumfac = s_lagran_dt/s_Tn
   epsi = dumfac / (1._rprec+dumfac)
 
   I_QN(:,:,jz) = epsi*QN + (1._rprec-epsi)*I_QN(:,:,jz)
   I_NN(:,:,jz) = epsi*NN + (1._rprec-epsi)*I_NN(:,:,jz)
   !Clip if necessary
   I_QN(:,:,jz) = max( zero, I_QN(:,:,jz) )

   !Calculate Ds_opt2 (for 4-delta filter)
   !Add +zero in denominator to avoid division by identically zero
   Ds_opt2_4d(:,:) = I_QN(:,:,jz)/(I_NN(:,:,jz) + zero)
   Ds_opt2_4d(ld,:) = zero
   Ds_opt2_4d(ld-1,:) = zero
   !Clip if necessary
   Ds_opt2_4d(:,:) = max( zero, Ds_opt2_4d(:,:) )    

   s_Beta(:,:,jz) = &
       (Ds_opt2_4d(:,:)/Ds_opt2_2d(:,:))**(log(tf1)/(log(tf2)-log(tf1)))
    
   !--MPI
   $if ($MPI)
       if ((coord==nproc-1).and.(jz==nz)) then
          s_Beta(:,:,jz) = 1._rprec
       end if
   $else
       if (jz==nz) then
          s_Beta(:,:,jz) = 1._rprec
       end if
   $endif

   !Clip s_Beta and 
   do jy=1,ny
      do jx=1,ld 
         Betaclip = max(s_Beta(jx,jy,jz),1._rprec/8._rprec)
         Ds_opt2(jx,jy,jz) = Ds_opt2_2d(jx,jy)/Betaclip
      end do
   end do
   Ds_opt2(ld,:,jz) = zero
   Ds_opt2(ld-1,:,jz) = zero


   !Count how often Ds is clipped
   !if (coord == 0) then
   do i=1,nx
      do j=1,ny
         if (real(Ds_opt2(i,j,jz)).lt.real(zero)) then
            count_ds2_clip = count_ds2_clip + 1
            ds2_clips(i,j,jz) = 1._rprec 
         end if    
         count_ds2_all = count_ds2_all + 1
       end do
       !print*, 'jt, coord, count_ds2_clip, sum(ds2_clips(:,:,jz))', jt, coord, count_ds2_clip, sum(ds2_clips(:,:,jz))
   end do
   !end if
   
   !Clip if necessary
   Ds_opt2(:,:,jz) = max( zero, Ds_opt2(:,:,jz) )

   !Save s_Tn to 3D array for use with tavg_scalar_sgs
   s_Tn_all(:,:,jz) = s_Tn(:,:)

end do

!Share new data between overlapping nodes
$if($MPI) 
   call mpi_sync_real_array( I_LM, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( I_MM, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( I_QN, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( I_NN, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( s_Tn_all, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( ds2_clips, 0, MPI_SYNC_DOWNUP )
   call mpi_sync_real_array( s_Beta, 0, MPI_SYNC_DOWNUP )
$endif

!Reset variable for use during next set of cs_count timesteps
if (use_cfl_dt ) s_lagran_dt = 0._rprec

end subroutine scalars_lagrange_Sdep





subroutine get_scalar_std(scalar,sigma_s)

use types, only: rprec
use param, only: ld, nx, ny, nz, lbz, coord, jt
use functions, only: interp_to_w_grid
$if($MPI)
use mpi_defs, only: mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: scalar
real(rprec), dimension(lbz:nz), intent(out) :: sigma_s
real(rprec), dimension(lbz:nz) :: scalar_avg_zplane, temp
real(rprec), dimension(ld,ny,lbz:nz) :: scalar_w
real(rprec):: favg
integer :: jx, jy, jz

scalar_w = interp_to_w_grid(scalar,lbz)

favg = real(nx*ny,kind=rprec)


do jz=lbz,nz
   scalar_avg_zplane(jz) = 0._rprec
   do jy=1,ny
      do jx=1,nx
         scalar_avg_zplane(jz) = scalar_avg_zplane(jz) + scalar_w(jx,jy,jz)
      end do
   end do 
   scalar_avg_zplane(jz) = scalar_avg_zplane(jz)/favg
end do

do jz=lbz,nz
   temp(jz) = 0.0_rprec
   do jy=1,ny
      do jx=1,nx
         temp(jz) = temp(jz) + (scalar_w(jx,jy,jz) - scalar_avg_zplane(jz))**2._rprec  
      end do
   end do
   sigma_s(jz) = sqrt(temp(jz)/favg)
end do

end subroutine get_scalar_std



end module s_lagrange_Sdep
