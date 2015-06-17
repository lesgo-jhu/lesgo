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

!***********************************************************************
module sgs_stag_util
!***********************************************************************
implicit none

save 
private

public sgs_stag, rtnewt

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sgs_stag ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Calculates turbulent (subgrid) stress for entire domain
!   using the model specified in param.f90 (Smag, LASD, etc)
!   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)
!   Sij values are stored on w-nodes (1:nz)
!
!   module is used to share Sij values b/w subroutines
!     (avoids memory error when arrays are very large)
!
! put everything onto w-nodes, follow original version

use types,only:rprec
use param
use sim_param,only: u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,  &
                    txx, txy, txz, tyy, tyz, tzz
use sgs_param
use messages

$if ($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
$endif

$if ($DEBUG)
use debug_mod
$endif

$if ($LVLSET)
  use level_set, only : level_set_BC, level_set_Cs
$endif

!use sgs_hist, only: sgs_hist_update_vals

implicit none

$if ($VERBOSE)
character (*), parameter :: sub_name = 'sgs_stag'
$endif

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

real(kind=rprec),dimension(nz)::l,ziko,zz
!real(kind=rprec),dimension(ld,ny) :: txzp, tyzp
real(kind=rprec) :: const
real(kind=rprec) :: const2 !RICHARD: USED FOR OPTIMIZATION
real(kind=rprec) :: const3 !RICHARD: USED FOR OPTIMIZATION
real(kind=rprec) :: const4 !RICHARD: USED FOR OPTIMIZATION

integer::jx,jy,jz
integer :: jz_min

$if ($VERBOSE)
call enter_sub (sub_name)
$endif

! Cs is Smagorinsky's constant. l is a filter size (non-dim.)  

$if ($DEBUG)
if (DEBUG) then
    call DEBUG_write (dudx(:, :, 1:nz), 'sgs_stag.dudx.a')
    call DEBUG_write (dudy(:, :, 1:nz), 'sgs_stag.dudy.a')
    call DEBUG_write (dudz(:, :, 1:nz), 'sgs_stag.dudz.a')
    call DEBUG_write (dvdx(:, :, 1:nz), 'sgs_stag.dvdx.a')
    call DEBUG_write (dvdy(:, :, 1:nz), 'sgs_stag.dvdy.a')
    call DEBUG_write (dvdz(:, :, 1:nz), 'sgs_stag.dvdz.a')
    call DEBUG_write (dwdx(:, :, 1:nz), 'sgs_stag.dwdx.a')
    call DEBUG_write (dwdy(:, :, 1:nz), 'sgs_stag.dwdy.a')
    call DEBUG_write (dwdz(:, :, 1:nz), 'sgs_stag.dwdz.a')
end if
$endif

! Calculate S12, S13, S23, etc.
call calc_Sij ()

$if ($DEBUG)
if (DEBUG) then
    call DEBUG_write (S11(:, :, 1:nz), 'sgs_stag.S11.b')
    call DEBUG_write (S12(:, :, 1:nz), 'sgs_stag.S12.b')
    call DEBUG_write (S13(:, :, 1:nz), 'sgs_stag.S13.b')
    call DEBUG_write (S22(:, :, 1:nz), 'sgs_stag.S22.b')
    call DEBUG_write (S23(:, :, 1:nz), 'sgs_stag.S23.b')
    call DEBUG_write (S33(:, :, 1:nz), 'sgs_stag.S33.b')
end if
$endif

! This approximates the sum displacement during cs_count timesteps
! This is used with the lagrangian model only
$if ($CFL_DT)
    if (sgs_model == 4 .OR. sgs_model==5) then
      if ( ( jt .GE. DYN_init-cs_count + 1 ) .OR.  initu ) then
        lagran_dt = lagran_dt + dt
      endif
    endif
$else
    lagran_dt = cs_count*dt
$endif


if (sgs) then 
    if((sgs_model == 1))then  ! Traditional Smagorinsky model

        $if ($LVLSET)
            l = delta
            call level_set_Cs (delta)
        $else
            ! Parameters (Co and nn) for wallfunction defined in param.f90
            Cs_opt2 = Co**2  ! constant coefficient
            
            if (lbc_mom == 0) then
               ! Stress free
                l = delta        
            else       
                ! The variable "l" calculated below is l_sgs/Co where l_sgs is from JDA eqn(2.30)
                if (coord == 0) then

                    ! z's nondimensional, l here is on uv-nodes
                    zz(1) = 0.5_rprec * dz                    
                    l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp) +  &
                                (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                
                    do jz = 2, nz
                        ! z's nondimensional, l here is on w-nodes
                        zz(jz) = (jz - 1) * dz                        
                        l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp) +  &
                                (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                    end do
                else
                    do jz = 1, nz
                        ! z's nondimensional, l here is on w-nodes
                        zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz                        
                        l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp) +  &
                                (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                    end do 
                end if

            end if
        $endif

    else    ! Dynamic procedures: modify/set Sij and Cs_opt2 (specific to sgs_model)
   
        l = delta       ! recall: l is the filter size
        
        if ((jt == 1) .and. (inilag)) then 
        ! Use the Smagorinsky model until DYN_init timestep
            print *,'CS_opt2 initialiazed'
            Cs_opt2 = 0.03_rprec

        elseif ( ((jt.GE.DYN_init).OR.(initu)) .AND. (mod(jt_total,cs_count)==0) ) then
        ! Update Sij, Cs every cs_count timesteps (specified in param)
        
            if (jt == DYN_init) print *,'running dynamic sgs_model = ',sgs_model
            
            if (sgs_model == 2) then        ! Standard dynamic model
                call std_dynamic(ziko)
                forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)
            else if (sgs_model==3) then     ! Plane average dynamic model
                call scaledep_dynamic(ziko)
                do jz = 1, nz
                    Cs_opt2(:, :, jz) = ziko(jz)
                end do
            else if (sgs_model==4.) then    ! Lagrangian scale similarity model
                call lagrange_Ssim()
            elseif (sgs_model==5) then      ! Lagrangian scale dependent model
                call lagrange_Sdep()
            end if       
        end if
 
    end if 
 
end if 


$if ($DEBUG)
if (DEBUG) then
    call DEBUG_write (Cs_opt2, 'sgs_stag.Cs_opt2')
end if
$endif

! Define |S| and eddy viscosity (nu_t= c_s^2 l^2 |S|) for entire domain
!   stored on w-nodes (on uvp node for jz=1 and 'wall' BC only) 
do jz=1,nz
do jy=1,ny
do jx=1,nx
    S(jx,jy) = sqrt( 2._rprec*(S11(jx,jy,jz)**2 +           S22(jx,jy,jz)**2 +&
                               S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +&
                               S13(jx,jy,jz)**2 +           S23(jx,jy,jz)**2 )))
    Nu_t(jx,jy,jz)=S(jx,jy)*Cs_opt2(jx,jy,jz)*l(jz)**2
end do
end do
end do

! Update the values for the sgs-variable histograms
!  if (sgs_hist_calc) then
!  if ( (jt_total .ge. sgs_hist_nstart) .and. (mod(jt_total,sgs_hist_nskip).eq.0) ) then
!    call sgs_hist_update_vals( )
!  endif
!  endif
  
! Calculate txx, txy, tyy, tzz for bottom level: jz=1 node (coord==0 only)
if (coord == 0) then

    select case (lbc_mom)

    case (0) ! Stress free
    ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
    !   recall: for this case, Sij are stored on w-nodes      
        
        if (sgs) then
            do jy=1,ny
            do jx=1,nx        
               ! Total viscosity
                const=0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2)) + nu
!                txx(jx,jy,1)=-2._rprec*(const+nu)*0.5_rprec*(S11(jx,jy,1) + S11(jx,jy,2))
!                txy(jx,jy,1)=-2._rprec*(const+nu)*0.5_rprec*(S12(jx,jy,1) + S12(jx,jy,2))
!                tyy(jx,jy,1)=-2._rprec*(const+nu)*0.5_rprec*(S22(jx,jy,1) + S22(jx,jy,2))
!                tzz(jx,jy,1)=-2._rprec*(const+nu)*0.5_rprec*(S33(jx,jy,1) + S33(jx,jy,2))
                txx(jx,jy,1)=-const*(S11(jx,jy,1) + S11(jx,jy,2))
                txy(jx,jy,1)=-const*(S12(jx,jy,1) + S12(jx,jy,2))
                tyy(jx,jy,1)=-const*(S22(jx,jy,1) + S22(jx,jy,2))
                tzz(jx,jy,1)=-const*(S33(jx,jy,1) + S33(jx,jy,2))
            end do
            end do
        else    
            const = 0._rprec
            do jy=1,ny
            do jx=1,nx      
                ! const = 0. (therefore, removed from expressions below)                
!                txx(jx,jy,1)=-2._rprec*(nu)*0.5_rprec*(S11(jx,jy,1) + S11(jx,jy,2))
!                txy(jx,jy,1)=-2._rprec*(nu)*0.5_rprec*(S12(jx,jy,1) + S12(jx,jy,2))
!                tyy(jx,jy,1)=-2._rprec*(nu)*0.5_rprec*(S22(jx,jy,1) + S22(jx,jy,2))
!                tzz(jx,jy,1)=-2._rprec*(nu)*0.5_rprec*(S33(jx,jy,1) + S33(jx,jy,2))
                txx(jx,jy,1)=-(nu)*(S11(jx,jy,1) + S11(jx,jy,2))
                txy(jx,jy,1)=-(nu)*(S12(jx,jy,1) + S12(jx,jy,2))
                tyy(jx,jy,1)=-(nu)*(S22(jx,jy,1) + S22(jx,jy,2))
                tzz(jx,jy,1)=-(nu)*(S33(jx,jy,1) + S33(jx,jy,2))
            end do
            end do
        end if                     

    case (1) ! Wall
       
    ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
    !   recall: for this case, Sij are stored on uvp-nodes        
        
        if (sgs) then
            do jy=1,ny
            do jx=1,nx         
!                const = Nu_t(jx,jy,1)                   
!                txx(jx,jy,1) = -2._rprec*(const+nu)*S11(jx,jy,1)
!                txy(jx,jy,1) = -2._rprec*(const+nu)*S12(jx,jy,1)
!                tyy(jx,jy,1) = -2._rprec*(const+nu)*S22(jx,jy,1)
!                tzz(jx,jy,1) = -2._rprec*(const+nu)*S33(jx,jy,1)
                const = -2._rprec*(Nu_t(jx,jy,1)+nu)                   
                txx(jx,jy,1) = const*S11(jx,jy,1)
                txy(jx,jy,1) = const*S12(jx,jy,1)
                tyy(jx,jy,1) = const*S22(jx,jy,1)
                tzz(jx,jy,1) = const*S33(jx,jy,1)
            end do
            end do
        else    
            const = 0._rprec
            do jy=1,ny
            do jx=1,nx      
                ! const = 0. (therefore, removed from expressions below)                
                txx(jx,jy,1) = -2._rprec*(nu)*S11(jx,jy,1)
                txy(jx,jy,1) = -2._rprec*(nu)*S12(jx,jy,1)
                tyy(jx,jy,1) = -2._rprec*(nu)*S22(jx,jy,1)
                tzz(jx,jy,1) = -2._rprec*(nu)*S33(jx,jy,1)
            end do
            end do
        end if       
       
    end select
  
    jz_min = 2      ! since first level already calculated
    
else   
    jz_min = 1
    
end if

! Calculate all tau for the rest of the domain
!   txx, txy, tyy, tzz not needed at nz (so they aren't calculated)
!     txz, tyz at nz will be done later
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)

if (sgs) then 
    const3=-2._rprec*(nu)*0.5_rprec
    const4=-2._rprec*(nu)
    do jz=jz_min, nz-1
    do jy=1,ny
    do jx=1,nx
              
       const =-0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))  
       const2=-2._rprec*Nu_t(jx,jy,jz)                

       txx(jx,jy,jz)=(const+const3)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
       txy(jx,jy,jz)=(const+const3)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
       tyy(jx,jy,jz)=(const+const3)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
       tzz(jx,jy,jz)=(const+const3)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))        
       txz(jx,jy,jz)=(const2+const4)* S13(jx,jy,jz)
       tyz(jx,jy,jz)=(const2+const4)* S23(jx,jy,jz)        
                
    end do
    end do
    end do

else    
    
    const=0._rprec  ! removed from tij expressions below since it's zero
    
    do jz=jz_min, nz-1
    do jy=1,ny
    do jx=1,nx
       
!       txx(jx,jy,jz)=-2._rprec*(nu)*0.5_rprec*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
!       txy(jx,jy,jz)=-2._rprec*(nu)*0.5_rprec*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
!       tyy(jx,jy,jz)=-2._rprec*(nu)*0.5_rprec*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
!       tzz(jx,jy,jz)=-2._rprec*(nu)*0.5_rprec*(S33(jx,jy,jz) + S33(jx,jy,jz+1))            
       txx(jx,jy,jz)=-(nu)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
       txy(jx,jy,jz)=-(nu)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
       tyy(jx,jy,jz)=-(nu)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
       tzz(jx,jy,jz)=-(nu)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))            

        
       txz(jx,jy,jz)=-2._rprec*(nu) * S13(jx,jy,jz)
       tyz(jx,jy,jz)=-2._rprec*(nu) * S23(jx,jy,jz)        
                
    end do
    end do
    end do

end if    
      
$if ($DEBUG)
if (DEBUG) then
    call DEBUG_write (txx(:, :, 1:nz), 'sgs_stag.txx.a')
    call DEBUG_write (txy(:, :, 1:nz), 'sgs_stag.txy.a')
    call DEBUG_write (txz(:, :, 1:nz), 'sgs_stag.txz.a')
end if
$endif

$if ($LVLSET)
  !--at this point tij are only set for 1:nz-1
  !--at this point u, v, w are set for 0:nz, except bottom process is 1:nz
  !--some MPI synchronizing may be done in here, but this will be kept
  !  separate from the rest of the code (at the risk of some redundancy)
  call level_set_BC ()
$endif

$if ($DEBUG)
if (DEBUG) then
    call DEBUG_write (txx(:, :, 1:nz), 'sgs_stag.txx.b')
    call DEBUG_write (txy(:, :, 1:nz), 'sgs_stag.txy.b')
    call DEBUG_write (txz(:, :, 1:nz), 'sgs_stag.txz.b')
end if
$endif

$if ($MPI)
    ! txz,tyz calculated for 1:nz-1 (on w-nodes) except bottom process
    ! (only 2:nz-1) exchange information between processors to set
    ! values at nz from jz=1 above to jz=nz below
    call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
    call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
    ! Set bogus values (easier to catch if there's an error)
$if (SAFETYMODE)
    txx(:, :, 0) = BOGUS
    txy(:, :, 0) = BOGUS
    txz(:, :, 0) = BOGUS
    tyy(:, :, 0) = BOGUS
    tyz(:, :, 0) = BOGUS
    tzz(:, :, 0) = BOGUS 
$endif
$endif

! Set bogus values (easier to catch if there's an error)
$if ($SAFETYMODE)
txx(:, :, nz) = BOGUS
txy(:, :, nz) = BOGUS
tyy(:, :, nz) = BOGUS
tzz(:, :, nz) = BOGUS
$endif 

$if ($MPI) 
  if (coord == nproc-1) then  !assuming stress-free lid?
    txz(:,:,nz)=0._rprec
    tyz(:,:,nz)=0._rprec
  end if
$else
  txz(:,:,nz)=0._rprec
  tyz(:,:,nz)=0._rprec
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine sgs_stag

!**********************************************************************
subroutine calc_Sij()
!**********************************************************************
! Calculate the resolved strain rate tensor, Sij = 0.5(djui - diuj)
!   values are stored on w-nodes (1:nz)

use types,only:rprec
use param
use sim_param,only: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use sgs_param
$if ($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
$endif

implicit none

integer::jx,jy,jz
integer :: jz_min

real (rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz

! Calculate Sij for jz=1 (coord==0 only)
!   stored on uvp-nodes (this level only) for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == 0) then

    select case (lbc_mom)

    case (0) ! Stress free

        do jy=1,ny
        do jx=1,nx              ! Sij values are supposed to be on w-nodes for this case
                                !   does that mean they (Sij) should all be zero?
            ux=dudx(jx,jy,1)    ! check this, WAS 0.5_rprec*(dudx(jx,jy,1) + dudx(jx,jy,1))
            uy=dudy(jx,jy,1)    ! check this
            uz=dudz(jx,jy,1)  
            vx=dvdx(jx,jy,1)    ! check this
            vy=dvdy(jx,jy,1)    ! check this
            vz=dvdz(jx,jy,1) 
            wx=dwdx(jx,jy,1)  
            wy=dwdy(jx,jy,1)  
            wz=0.5_rprec*(dwdz(jx,jy,1) + 0._rprec)     ! check this
                              
            ! these values are stored on w-nodes
            S11(jx,jy,1)=ux          
            S12(jx,jy,1)=0.5_rprec*(uy+vx) 
            S13(jx,jy,1)=0.5_rprec*(uz+wx) 
            S22(jx,jy,1)=vy          
            S23(jx,jy,1)=0.5_rprec*(vz+wy) 
            S33(jx,jy,1)=wz          
        end do
        end do

    case (1) ! Wall
    ! recall dudz and dvdz are stored on uvp-nodes for first level only, 'wall' only
    ! recall dwdx and dwdy are stored on w-nodes (always)
        do jy=1,ny
        do jx=1,nx              
!            ux=dudx(jx,jy,1)  
!            uy=dudy(jx,jy,1)  
!            uz=dudz(jx,jy,1)    ! dudz on uvp-node for jz==1 (wallstress.f90)
!            vx=dvdx(jx,jy,1)  
!            vy=dvdy(jx,jy,1)  
!            vz=dvdz(jx,jy,1)    ! dvdz on uvp-node for jz==1 (wallstress.f90)
!            wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2)) 
!            wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  
!            wz=dwdz(jx,jy,1) 
         
            ! these values stored on uvp-nodes
!            S11(jx,jy,1)=ux         
!            S12(jx,jy,1)=0.5_rprec*(uy+vx) 
!            S13(jx,jy,1)=0.5_rprec*(uz+wx) 
!            S22(jx,jy,1)=vy          
!            S23(jx,jy,1)=0.5_rprec*(vz+wy) 
!            S33(jx,jy,1)=wz         
          
            ! these values stored on uvp-nodes
            S11(jx,jy,1)=dudx(jx,jy,1)         
            S12(jx,jy,1)=0.5_rprec*(dudy(jx,jy,1)+dvdx(jx,jy,1)) 
            wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2)) 
            S13(jx,jy,1)=0.5_rprec*(dudz(jx,jy,1)+wx) 
            S22(jx,jy,1)=dvdy(jx,jy,1)          
            wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))             
            S23(jx,jy,1)=0.5_rprec*(dvdz(jx,jy,1)+wy) 
            S33(jx,jy,1)=dwdz(jx,jy,1)         

        end do
        end do
  
    end select
  
    jz_min = 2      ! since first level already calculated

else
    jz_min = 1

end if

$if ($MPI)
    ! dudz calculated for 0:nz-1 (on w-nodes) except bottom process
    ! (only 1:nz-1) exchange information between processors to set
    ! values at nz from jz=1 above to jz=nz below
    call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
$endif

! Calculate Sij for the rest of the domain
!   values are stored on w-nodes
!   dudz, dvdz, dwdx, dwdy are already stored on w-nodes
do jz=jz_min, nz
do jy=1,ny
do jx=1,nx              
!    ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1)) 
!    uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  
!    uz=dudz(jx,jy,jz)  
!    vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  
!    vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  
!    vz=dvdz(jx,jy,jz)  
!    wx=dwdx(jx,jy,jz)  
!    wy=dwdy(jx,jy,jz)  
!    wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1)) 
!    S11(jx,jy,jz)=ux         
!    S12(jx,jy,jz)=0.5_rprec*(uy+vx) 
!    S13(jx,jy,jz)=0.5_rprec*(uz+wx) 
!    S22(jx,jy,jz)=vy          
!    S23(jx,jy,jz)=0.5_rprec*(vz+wy) 
!    S33(jx,jy,jz)=wz          
    S11(jx,jy,jz)=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))          
    uy=(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  
    vx=(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  
    S12(jx,jy,jz)=0.25_rprec*(uy+vx) 
    S13(jx,jy,jz)=0.5_rprec*(dudz(jx,jy,jz) + dwdx(jx,jy,jz)) 
    S22(jx,jy,jz)=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))          
    S23(jx,jy,jz)=0.5_rprec*(dvdz(jx,jy,jz) + dwdy(jx,jy,jz)) 
    S33(jx,jy,jz)=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))          
end do
end do
end do

end subroutine calc_Sij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------
real(kind=rprec) function rtnewt(A, jz)
!-----------------------------------------------------------------------
use types,only:rprec
integer,parameter :: jmax=100
real(kind=rprec) :: x1,x2,xacc
integer :: j, jz
real(kind=rprec) :: df,dx,f
real(kind=rprec), dimension(0:5) :: A
x1 = 0._rprec
x2 = 15._rprec  ! try to find the largest root first....hmm
xacc = 0.001_rprec ! doesn't need to be that accurate...
rtnewt = 0.5_rprec*(x1+x2)
do j=1,jmax
   f = A(0)+rtnewt*(A(1)+rtnewt*(A(2)+rtnewt*(A(3)+rtnewt*(A(4)+rtnewt*A(5)))))
   df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) +&
        rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))
   dx=f/df
   rtnewt = rtnewt - dx
!        if ((x1-rtnewt)*(rtnewt-x2) < 0.) STOP 'rtnewt out of bounds'
   if (abs(dx) < xacc) return
end do
rtnewt = 1._rprec  ! if dont converge fast enough
write(6,*) 'using beta=1 at jz= ', jz
end function rtnewt
!-----------------------------------------------------------------------

end module sgs_stag_util
