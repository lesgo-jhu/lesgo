!!
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
!--provides Ds_1D 1:nz
!--MPI: requires u,v 0:nz, except bottom process only 1:nz
subroutine scalars_scaledep_dynamic(jt,scalar,dsdx,dsdy,dsdz,Ds_1D)
! standard dynamic model to calculate the equivalent of the 
! Smagorinsky coefficient for the eddy diffusivity for scalars
! this is done layer-by-layer to save memory
! note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary.  mem requirement ~ n^2, not n^3
use types,only:rprec
use param,only:path,ld,nx,ny,nz,dx,dy,dz,coord, lbz
use sim_param,only:u,v,w
use sgs_stag_util,only:rtnewt
use sgs_param,only:S11,S12,S13,S22,S23,S33,delta,S
use sgs_param,only:S_bar,S11_bar,S12_bar,S13_bar,S22_bar,S23_bar,S33_bar
use sgs_param,only:S_S11_bar,S_S12_bar,S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
use sgs_param,only:u_bar,v_bar,w_bar,scalar_bar
use sgs_param,only:u_hat,v_hat,w_hat,scalar_hat
use sgs_param,only:S_hat,S11_hat,S12_hat,S13_hat,S22_hat,S23_hat,S33_hat
use sgs_param,only:S_S11_hat,S_S12_hat,S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat
use sgs_param,only:L1,L2,L3
use sgs_param,only:dsdx_bar, dsdy_bar, dsdz_bar, S_dsdx_bar, S_dsdy_bar, S_dsdz_bar
use sgs_param,only:dsdx_hat, dsdy_hat, dsdz_hat, S_dsdx_hat, S_dsdy_hat, S_dsdz_hat
use test_filtermodule
use sgs_param, only:ee_now
implicit none

integer, intent(in) :: jt
integer :: jz
real(rprec), dimension(nz), intent (out) :: Ds_1D
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: scalar
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdx
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdy
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdz

real(rprec), save, allocatable, target, dimension(:,:) :: Q1,Q2,Q3
real(rprec), pointer, dimension(:,:) :: M1,M2,M3

real(rprec), save, dimension(:), allocatable :: beta

logical, save :: arrays_allocated = .false. 

real(rprec), parameter :: zero = 1e-24
real(rprec) :: const
real(rprec), dimension(0:5) :: A
real(rprec) :: a1, b1, c1, d1, e1, a2, b2, c2, d2, e2

!TSreal(kind=rprec) :: rtnewt
real(rprec), dimension(ld,ny) :: LM,MM

! Allocate arrays
if( .not. arrays_allocated ) then

   allocate ( Q1(ld,ny), Q2(ld,ny), Q3(ld,ny) )

   allocate ( beta(nz) )

   arrays_allocated = .true. 

endif

! Associate pointers
M1 => Q1
M2 => Q2
M3 => Q3


do jz=1,nz
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

   ! Filter first term and add the second term to get the final value
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
   
! calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +& 
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

! already on w-nodes
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

   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

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
      dsdz_bar(:,:) = dsdz(:,:,jz) !Try this
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

! note: check that the Nyquist guys are zero!
! the 1./(nx*ny) is not really neccessary, but in practice it does 
! affect the results.
   a1 = -(delta**2) * 4._rprec * sum( S_bar*(dsdx_bar*L1 + dsdy_bar*L2 +&
        dsdz_bar*L3) ) / (nx*ny)
   b1 = -(delta**2) * sum ( S_dsdx_bar*L1 + S_dsdy_bar*L2 + S_dsdz_bar*L3 ) / (nx*ny)
   c1 = (delta**2)**2 * sum( S_dsdx_bar**2 + S_dsdy_bar**2 + S_dsdz_bar**2 ) / (nx*ny)
   d1 = (delta**2)**2 * 16._rprec * sum ( S_bar**2 * (dsdx_bar**2 + &
        dsdy_bar**2 + dsdz_bar**2) ) / (nx*ny)
   e1 = 2._rprec * (delta**2)**2 * 4._rprec * sum( S_bar*(dsdx_bar*S_dsdx_bar + &
        dsdy_bar*S_dsdy_bar + dsdz_bar*S_dsdz_bar) )/ (nx*ny)                               
     
   a2 = -(delta**2) * 16._rprec * sum( S_hat*(dsdx_hat*Q1 + dsdy_hat*Q2 +&
        dsdz_hat*Q3) ) / (nx*ny)
   b2 = -(delta**2) * sum( S_dsdx_hat*Q1 + S_dsdy_hat*Q2 + S_dsdz_hat*Q3 ) / (nx*ny)
   c2 = (delta**2)**2 * sum (S_dsdx_hat**2 + S_dsdy_hat**2 + S_dsdz_hat*Q3 ) / (nx*ny)
   d2 = (delta**2)**2 * 256._rprec * sum (S_hat**2 * (dsdx_hat**2 + & 
        dsdy_hat**2 + dsdz_hat**2) ) / (nx*ny)
   e2 = 2._rprec * (delta**2)**2 * 16._rprec *sum( S_hat*(dsdx_hat*S_dsdx_hat + &
        dsdy_hat*S_dsdy_hat + dsdz_hat*S_dsdz_hat) ) / (nx*ny)

   A(0) = b2*c1 - b1*c2
   A(1) = a1*c2 - b2*e1
   A(2) = b2*d1 + b1*e2 - a2*c1
   A(3) = a2*e1 - a1*e2
   A(4) = -a2*d1 - b1*d2
   A(5) = a1*d2

! note: this is not the temperature-beta
! I think we might need to do something to ensure that we always get the 
!  largest root, as in the JFM paper
   beta(jz) = rtnewt(A,jz,jt)  ! subroutine at end of this file

! now put beta back into M_ij: using Q_ij as storage
   const = delta**2
! it might be faster just to go ahead and not store M_ij; just code directly
   M1 = const*(S_dsdx_bar - 4._rprec*beta(jz)*S_bar*dsdx_bar)
   M2 = const*(S_dsdy_bar - 4._rprec*beta(jz)*S_bar*dsdy_bar)
   M3 = const*(S_dsdz_bar - 4._rprec*beta(jz)*S_bar*dsdz_bar)
     
   if ( (coord.le.2).and.(jz.ge.1) ) then 
   !print*, 'jt, coord, jz, L1(25,25)', jt, coord, jz, L1(25,25)
   !print*, 'jt, coord, jz, L2(25,25)', jt, coord, jz, L2(25,25)
   !print*, 'jt, coord, jz, L3(25,25)', jt, coord, jz, L3(25,25)
   !print*, 'jt, coord, jz, M1(25,25)', jt, coord, jz, M1(25,25)
   !print*, 'jt, coord, jz, M2(25,25)', jt, coord, jz, M2(25,25)
   !print*, 'jt, coord, jz, M3(25,25)', jt, coord, jz, M3(25,25)
   !print*, 'jt, coord, jz, S_bar(25,25)', jt, coord, jz, S_bar(25,25)
   !print*, 'jt, coord, jz, dsdx_bar(25,25)', jt, coord, jz, dsdx_bar(25,25)
   !print*, 'jt, coord, jz, dsdy_bar(25,25)', jt, coord, jz, dsdy_bar(25,25)
   !print*, 'jt, coord, jz, dsdz_bar(25,25)', jt, coord, jz, dsdz_bar(25,25)
   !print*, 'jt, coord, jz, S_dsdx_bar(25,25)', jt, coord, jz, S_dsdx_bar(25,25)
   !print*, 'jt, coord, jz, S_dsdy_bar(25,25)', jt, coord, jz, S_dsdy_bar(25,25)
   !print*, 'jt, coord, jz, S_dsdz_bar(25,25)', jt, coord, jz, S_dsdz_bar(25,25)
   !print*, 'jt, coord, jz, dsdx(25,25,jz)', jt, coord, jz, dsdx(25,25,jz)
   !print*, 'jt, coord, jz, dsdy(25,25,jz)', jt, coord, jz, dsdy(25,25,jz)
   !print*, 'jt, coord, jz, dsdz(25,25,jz)', jt, coord, jz, dsdz(25,25,jz)
   
   !print*, 'jt, coord, jz, scalar(25,25,jz)', jt, coord, jz, scalar(25,25,jz)
   end if

   Ds_1D(jz) = sum(L1*M1 + L2*M2 + L3*M3)/sum(M1**2 + M2**2 + M3**2 +zero)
 
   Ds_1D(jz) = max(0._rprec, real(Ds_1D(jz),kind=rprec))

    ! Calculate ee_now (the current value of eij*eij) 
    !LM=L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
    !MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)
    !ee_now(:,:,jz) = L11**2+L22**2+L33**2+2._rprec*(L12**2+L13**2+L23**2) &
    !                -2._rprec*LM*Cs_1D(jz) + MM*Cs_1D(jz)**2  
end do   

! Nullify pointers
nullify( M1, M2, M3 )

end subroutine scalars_scaledep_dynamic
