!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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

!*******************************************************************************
module sgs_stag_util
!*******************************************************************************
implicit none

save
private

public sgs_stag, rtnewt

contains

!*******************************************************************************
subroutine sgs_stag ()
!*******************************************************************************
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

use types, only : rprec
use param
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
use sgs_param
use messages

#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif

#ifdef PPLVLSET
use level_set, only : level_set_BC, level_set_Cs
#endif

implicit none

character (*), parameter :: sub_name = 'sgs_stag'

real(rprec), dimension(nz) :: l, ziko, zz
real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: jz_min, jz_max

! Cs is Smagorinsky's constant. l is a filter size (non-dim.)
call calc_Sij ()

! This approximates the sum displacement during cs_count timesteps
! This is used with the lagrangian model only
#ifdef PPCFL_DT
if (sgs_model == 4 .OR. sgs_model==5) then
    if ( ( jt .GE. DYN_init-cs_count + 1 ) .OR.  initu ) then
        lagran_dt = lagran_dt + dt
    endif
endif
#else
lagran_dt = cs_count*dt
#endif


if (sgs) then
    ! Traditional Smagorinsky model
    if (sgs_model == 1) then

#ifdef PPLVLSET
        l = delta
        call level_set_Cs (delta)
#else
        ! Parameters (Co and nn) for wallfunction defined in param.f90
        Cs_opt2 = Co**2  ! constant coefficient

        ! both Stress free
        if (lbc_mom == 0 .and. ubc_mom == 0) then
            l = delta

        ! top Stress free, bottom wall
        else if (lbc_mom > 0 .and. ubc_mom == 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == 0) then
                ! z's nondimensional, l here is on uv-nodes
                zz(1) = 0.5_rprec * dz
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)    &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_min = 2
            else
                jz_min = 1
            end if

            do jz = jz_min, nz
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)  &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! both top and bottom walls, zz is distance to nearest wall
        else if (lbc_mom > 0 .and. ubc_mom > 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == 0) then
                ! z's nondimensional, l here is on uv-nodes
                zz(1) = 0.5_rprec * dz
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)&
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)

                jz_min = 2
            else
                jz_min = 1
            end if

            if (coord == nproc-1) then
                ! z's nondimensional, l here is on uv-nodes
                zz(nz) = 0.5_rprec * dz
                l(nz) = (Co**(wall_damp_exp)*(vonk*zz(nz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_max = nz-1
            else
                jz_max = nz
            end if

            do jz = jz_min, jz_max
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                zz(jz) = min( zz(jz), (nz-1)*nproc*dz - zz(jz) )
                l(jz) = (Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! top wall, bottom stress free, zz is distance to top
        else if (lbc_mom == 0 .and. ubc_mom > 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == nproc-1) then
                ! z's nondimensional, l here is on uv-nodes
                zz(nz) = 0.5_rprec * dz
                l(nz) = (Co**(wall_damp_exp)*(vonk*zz(nz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_max = nz-1
            else
                jz_max = nz
            end if

            do jz = 1, jz_max
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((nproc - coord)*(nz - 1) - (jz - 1)) * dz
                l(jz) = (Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! Invalid combination
        else
            call error (sub_name, 'invalid b.c. combination')
        end if
#endif

    ! Dynamic procedures: modify/set Sij and Cs_opt2 (specific to sgs_model)
    else
        ! recall: l is the filter size
        l = delta

        ! Use the Smagorinsky model until DYN_init timestep
        if ((jt == 1) .and. (inilag)) then
            write(*,*) 'CS_opt2 initialiazed'
            Cs_opt2 = 0.03_rprec

        ! Update Sij, Cs every cs_count timesteps (specified in param)
        elseif ( ((jt.GE.DYN_init).OR.(initu)) .AND.                           &
            (mod(jt_total,cs_count)==0) ) then

            if (jt == DYN_init) then
                write(*,*) 'running dynamic sgs_model = ', sgs_model
            end if

            ! Standard dynamic model
            if (sgs_model == 2) then
                call std_dynamic(ziko)
                forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)

            ! Plane average dynamic model
            else if (sgs_model==3) then
                call scaledep_dynamic(ziko)
                do jz = 1, nz
                    Cs_opt2(:, :, jz) = ziko(jz)
                end do
            ! Lagrangian scale similarity model
            else if (sgs_model==4.) then
                call lagrange_Ssim()

            ! Lagrangian scale dependent model
            elseif (sgs_model==5) then
                call lagrange_Sdep()
            end if

        end if
    end if
end if


! Define |S| and eddy viscosity (nu_t= c_s^2 l^2 |S|) for entire domain
!   stored on w-nodes (on uvp node for jz=1 or nz for 'wall' BC only)
do jz = 1, nz
do jy = 1, ny
do jx = 1, nx
    S(jx,jy) = sqrt( 2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +           &
        S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +                        &
        S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2 )))
    Nu_t(jx,jy,jz) = S(jx,jy)*Cs_opt2(jx,jy,jz)*l(jz)**2
end do
end do
end do

! Calculate txx, txy, tyy, tzz for bottom level: jz=1 node (coord==0 only)
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on w-nodes
        case (0)
            if (sgs) then
                do jy = 1, ny
                do jx = 1, nx
                   ! Total viscosity
                    const = 0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2)) + nu
                    txx(jx,jy,1) = -const*(S11(jx,jy,1) + S11(jx,jy,2))
                    txy(jx,jy,1) = -const*(S12(jx,jy,1) + S12(jx,jy,2))
                    tyy(jx,jy,1) = -const*(S22(jx,jy,1) + S22(jx,jy,2))
                    tzz(jx,jy,1) = -const*(S33(jx,jy,1) + S33(jx,jy,2))
                end do
                end do
            else
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,1) = -(nu)*(S11(jx,jy,1) + S11(jx,jy,2))
                    txy(jx,jy,1) = -(nu)*(S12(jx,jy,1) + S12(jx,jy,2))
                    tyy(jx,jy,1) = -(nu)*(S22(jx,jy,1) + S22(jx,jy,2))
                    tzz(jx,jy,1) = -(nu)*(S33(jx,jy,1) + S33(jx,jy,2))
                end do
                end do
            end if

        ! Wall
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on uvp-nodes
        case (1:)
            if (sgs) then
                do jy = 1, ny
                do jx = 1, nx
                    const = -2._rprec*(Nu_t(jx,jy,1)+nu)
                    txx(jx,jy,1) = const*S11(jx,jy,1)
                    txy(jx,jy,1) = const*S12(jx,jy,1)
                    tyy(jx,jy,1) = const*S22(jx,jy,1)
                    tzz(jx,jy,1) = const*S33(jx,jy,1)
                end do
                end do
            else
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,1) = -2._rprec*(nu)*S11(jx,jy,1)
                    txy(jx,jy,1) = -2._rprec*(nu)*S12(jx,jy,1)
                    tyy(jx,jy,1) = -2._rprec*(nu)*S22(jx,jy,1)
                    tzz(jx,jy,1) = -2._rprec*(nu)*S33(jx,jy,1)
                end do
                end do
            end if

    end select

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if


! Calculate txx, txy, tyy, tzz for bottom level: jz=nz node (coord==nproc-1)
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on w-nodes
        case (0)

            if (sgs) then
                do jy = 1, ny
                do jx = 1, nx
                    ! Total viscosity
                    const = 0.5_rprec*(Nu_t(jx,jy,nz-1) + Nu_t(jx,jy,nz)) + nu
                    const2 = 2._rprec*(Nu_t(jx,jy,nz-1) + nu)

                    ! for top wall, it is nz-1 on the uv-grid
                    txx(jx,jy,nz-1) = -const*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                    txy(jx,jy,nz-1) = -const*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                    tyy(jx,jy,nz-1) = -const*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                    tzz(jx,jy,nz-1) = -const*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1) = -const2*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1) = -const2*S23(jx,jy,nz-1)
                end do
                end do
            else
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,nz-1) = -(nu)*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                    txy(jx,jy,nz-1) = -(nu)*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                    tyy(jx,jy,nz-1) = -(nu)*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                    tzz(jx,jy,nz-1) = -(nu)*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1) = -2._rprec*(nu)*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1) = -2._rprec*(nu)*S23(jx,jy,nz-1)
                end do
                end do
            end if

        ! Wall
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on uvp-nodes
        case (1:)
            if (sgs) then
                do jy = 1, ny
                do jx = 1, nx
                    const = -2._rprec*(Nu_t(jx,jy,nz)+nu)
                    const2 = -2._rprec*(Nu_t(jx,jy,nz-1) + nu)

                    ! Note: Sij(nz) is on uvp-node at nz-1
                    txx(jx,jy,nz-1) = const*S11(jx,jy,nz)
                    txy(jx,jy,nz-1) = const*S12(jx,jy,nz)
                    tyy(jx,jy,nz-1) = const*S22(jx,jy,nz)
                    tzz(jx,jy,nz-1) = const*S33(jx,jy,nz)
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1)= const2*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1)= const2*S23(jx,jy,nz-1)
                end do
                end do
            else
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,nz-1) = -2._rprec*(nu)*S11(jx,jy,nz-1)
                    txy(jx,jy,nz-1) = -2._rprec*(nu)*S12(jx,jy,nz-1)
                    tyy(jx,jy,nz-1) = -2._rprec*(nu)*S22(jx,jy,nz-1)
                    tzz(jx,jy,nz-1) = -2._rprec*(nu)*S33(jx,jy,nz-1)
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1)=-2._rprec*(nu)*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1)=-2._rprec*(nu)*S23(jx,jy,nz-1)
                end do
                end do
            end if

    end select

    ! since last level already calculated
    jz_max = nz-2
else
    jz_max = nz-1
end if

! Calculate all tau for the rest of the domain
!   txx, txy, tyy, tzz not needed at nz (so they aren't calculated)
!     txz, tyz at nz will be done later
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)

if (sgs) then
    const3=-2._rprec*(nu)*0.5_rprec
    const4=-2._rprec*(nu)
    do jz=jz_min, jz_max
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

    do jz = jz_min, jz_max
    do jy = 1, ny
    do jx = 1, nx
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

#ifdef PPLVLSET
!--at this point tij are only set for 1:nz-1
!--at this point u, v, w are set for 0:nz, except bottom process is 1:nz
!--some MPI synchronizing may be done in here, but this will be kept
!  separate from the rest of the code (at the risk of some redundancy)
call level_set_BC ()
#endif


#ifdef PPMPI
! txz,tyz calculated for 1:nz-1 (on w-nodes) except bottom process
! (only 2:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
#ifdef PPSAFETYMODE
! Set bogus values (easier to catch if there's an error)
txx(:, :, 0) = BOGUS
txy(:, :, 0) = BOGUS
txz(:, :, 0) = BOGUS
tyy(:, :, 0) = BOGUS
tyz(:, :, 0) = BOGUS
tzz(:, :, 0) = BOGUS
#endif
#endif

! Set bogus values (easier to catch if there's an error)
#ifdef PPSAFETYMODE
txx(:, :, nz) = BOGUS
txy(:, :, nz) = BOGUS
tyy(:, :, nz) = BOGUS
tzz(:, :, nz) = BOGUS
#endif

end subroutine sgs_stag

!*******************************************************************************
subroutine calc_Sij()
!*******************************************************************************
! Calculate the resolved strain rate tensor, Sij = 0.5(djui - diuj)
!   values are stored on w-nodes (1:nz)

use types, only : rprec
use param
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dwdz_var
use sgs_param
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
implicit none

real(rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
integer :: jx, jy, jz
integer :: jz_min, jz_max


! Calculate Sij for jz=1 (coord==0 only)
!   stored on uvp-nodes (this level only) for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free
        case (0)
            do jy = 1, ny
            do jx = 1, nx
                ! Sij values are supposed to be on w-nodes for this case
                !   does that mean they (Sij) should all be zero?
                ! Check ux, uy, vx, vy, and qz
                ux = dudx(jx,jy,1) ! was 0.5_rprec*(dudx(jx,jy,1) + dudx(jx,jy,1))
                uy = dudy(jx,jy,1)
                uz = dudz(jx,jy,1)
                vx = dvdx(jx,jy,1)
                vy = dvdy(jx,jy,1)
                vz = dvdz(jx,jy,1)
                wx = dwdx(jx,jy,1)
                wy = dwdy(jx,jy,1)
                wz = 0.5_rprec*(dwdz(jx,jy,1) + 0._rprec)

                ! these values are stored on w-nodes
                S11(jx,jy,1) = ux
                S12(jx,jy,1) = 0.5_rprec*(uy+vx)
                S13(jx,jy,1) = 0.5_rprec*(uz+wx)
                S22(jx,jy,1) = vy
                S23(jx,jy,1) = 0.5_rprec*(vz+wy)
                S33(jx,jy,1) = wz
            end do
            end do

        ! Wall
        ! recall dudz and dvdz are stored on uvp-nodes for first level only,
        !   'wall' only
        ! recall dwdx and dwdy are stored on w-nodes (always)
        case (1:)
            do jy=1,ny
            do jx=1,nx
                ! these values stored on uvp-nodes
                S11(jx,jy,1) = dudx(jx,jy,1)
                S12(jx,jy,1) = 0.5_rprec*(dudy(jx,jy,1)+dvdx(jx,jy,1))
                wx = 0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))
                S13(jx,jy,1) = 0.5_rprec*(dudz(jx,jy,1)+wx)
                S22(jx,jy,1) = dvdy(jx,jy,1)
                wy = 0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))
                S23(jx,jy,1) = 0.5_rprec*(dvdz(jx,jy,1)+wy)
                S33(jx,jy,1) = dwdz(jx,jy,1)
            end do
            end do

    end select

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if

! Calculate Sij for jz=nz (coord==nproc-1 only)
!   stored on uvp-nodes (this level only nz on w-grid --> nz-1 on uvp-grid)
!       for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        case (0)

            do jy=1,ny
            do jx=1,nx
                ! Sij values are supposed to be on w-nodes for this case
                !   does that mean they (Sij) should all be zero?
                ux = dudx(jx,jy,nz-1)
                uy = dudy(jx,jy,nz-1)
                uz = dudz(jx,jy,nz)   ! this comes from wallstress() i.e. zero
                vx = dvdx(jx,jy,nz-1)
                vy = dvdy(jx,jy,nz-1)
                vz = dvdz(jx,jy,nz)   ! this comes from wallstress() i.e. zero
                wx = dwdx(jx,jy,nz)
                wy = dwdy(jx,jy,nz)
                wz = 0.5_rprec*(dwdz(jx,jy,nz-1) + 0._rprec)

                ! these values are stored on w-nodes
                S11(jx,jy,nz) = ux
                S12(jx,jy,nz) = 0.5_rprec*(uy+vx)
                S13(jx,jy,nz) = 0.5_rprec*(uz+wx)
                S22(jx,jy,nz) = vy
                S23(jx,jy,nz) = 0.5_rprec*(vz+wy)
                S33(jx,jy,nz) = wz
            end do
            end do

        ! Wall
        ! recall dudz and dvdz are stored on uvp-nodes for first level only,
        !   'wall' only
        ! recall dwdx and dwdy are stored on w-nodes (always)
        case (1:)
            do jy = 1, ny
            do jx = 1, nx
                ! these values stored on uvp-nodes
                S11(jx,jy,nz) = dudx(jx,jy,nz-1)
                S12(jx,jy,nz) = 0.5_rprec*(dudy(jx,jy,nz-1)+dvdx(jx,jy,nz-1))
                wx = 0.5_rprec*(dwdx(jx,jy,nz-1)+dwdx(jx,jy,nz))
                ! dudz from wallstress()
                S13(jx,jy,nz) = 0.5_rprec*(dudz(jx,jy,nz)+wx)
                S22(jx,jy,nz) = dvdy(jx,jy,nz-1)
                wy = 0.5_rprec*(dwdy(jx,jy,nz-1)+dwdy(jx,jy,nz))
                ! dvdz from wallstress()
                S23(jx,jy,nz) = 0.5_rprec*(dvdz(jx,jy,nz)+wy)
                S33(jx,jy,nz) = dwdz(jx,jy,nz-1)
            end do
            end do

    end select

    ! since last level already calculated
    jz_max = nz-1
else
    jz_max = nz
end if

#ifdef PPMPI
! dudz calculated for 0:nz-1 (on w-nodes) except bottom process
! (only 1:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
call dwdz_var%sync_down
! call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#endif

! Calculate Sij for the rest of the domain
!   values are stored on w-nodes
!   dudz, dvdz, dwdx, dwdy are already stored on w-nodes
do jz = jz_min, jz_max
do jy = 1, ny
do jx = 1, nx
    S11(jx,jy,jz) = 0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))
    uy = (dudy(jx,jy,jz) + dudy(jx,jy,jz-1))
    vx = (dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))
    S12(jx,jy,jz) = 0.25_rprec*(uy+vx)
    S13(jx,jy,jz) = 0.5_rprec*(dudz(jx,jy,jz) + dwdx(jx,jy,jz))
    S22(jx,jy,jz) = 0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))
    S23(jx,jy,jz) = 0.5_rprec*(dvdz(jx,jy,jz) + dwdy(jx,jy,jz))
    S33(jx,jy,jz) = 0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))
end do
end do
end do

end subroutine calc_Sij

!*******************************************************************************
real(rprec) function rtnewt(A, jz)
!*******************************************************************************
use types, only : rprec
integer, parameter :: jmax=100
real(rprec) :: x1, x2, xacc
real(rprec) :: df,dx,f
integer :: j, jz
real(rprec), dimension(0:5) :: A
x1 = 0._rprec
x2 = 15._rprec                  ! try to find the largest root first
xacc = 0.001_rprec              ! doesn't need to be that accurate
rtnewt = 0.5_rprec*(x1+x2)
do j = 1, jmax
    f = A(0)+rtnewt*(A(1)+rtnewt*(A(2)+rtnewt*(A(3)+rtnewt*(A(4)+rtnewt*A(5)))))
    df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) +                &
        rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))
    dx = f/df
    rtnewt = rtnewt - dx
    if (abs(dx) < xacc) return
end do
rtnewt = 1._rprec  ! if don't converge fast enough
write(6,*) 'using beta=1 at jz= ', jz

end function rtnewt

end module sgs_stag_util
