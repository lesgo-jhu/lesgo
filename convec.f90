!!
!!  Copyright (C) 2009-2018  Johns Hopkins University
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
subroutine convec
!*******************************************************************************
!
! Computes the rotation convective term in physical space
!       c = - (u X vort)
! Stores the convective term in RHS
! Uses 3/2-rule for dealiasing for more info see Canuto 1991 Spectral Methods,
! Chapter 7
!
use types, only : rprec
use param
use sim_param, only : u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy
use sim_param, only : RHSx, RHSy, RHSz
use sim_param, only : RHSx_var, RHSy_var, RHSz_var
use sim_param, only : u_big, v_big, w_big, u_big_var, v_big_var, w_big_var
use sim_param, only : u_var, v_var, w_var
use sim_param, only : vortx_big_var, vorty_big_var, vortz_big_var
use sim_param, only : vortx_big, vorty_big, vortz_big
use sim_param, only : RHSx_var, RHSy_var, RHSz_var
use sim_param, only : cc_big, cc_big_var
use fft

implicit none

integer :: jz
integer :: jz_min
integer :: jzLo, jzHi, jz_max  ! added for full channel capabilities

real(rprec) :: const

if (sgs) then
    jzLo = 2        !! necessary for LES or else blows up ....?
    jzHi = nz-1     !! can remove after testing
else
    jzLo = 1        !! for DNS
    jzHi = nz-1     !! can remove after testing
endif

! Calculate vorticity, using RHSx, RHSy, and RHSz as temporary storage locations
RHSx(:,:,1:nz) = dwdy(:,:,1:nz) - dvdz(:,:,1:nz)
RHSy(:,:,1:nz) = dudz(:,:,1:nz) - dwdx(:,:,1:nz)
RHSz(:,:,1:nz) = dvdx(:,:,1:nz) - dudy(:,:,1:nz)

! Bottom boundary condition
if (coord == 0) then
    if (lbc_mom == 0) then
        ! Stress free
        RHSx(:,:,1) = 0._rprec
        RHSy(:,:,1) = 0._rprec
    else
        ! Wall
        ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
        ! definition of the vorticity in that case which also interpolates
        ! dwdx, dwdy to the u-node at jz=1
        ! dwdy(1) and dwdx(1) should be 0
        RHSx(:,:,1) = 0.5_rprec*(dwdy(:,:,1) + dwdy(:,:,2)) - dvdz(:,:,1)
        RHSy(:,:,1) = dudz(:,:,1) - 0.5_rprec*(dwdx(:,:,1) + dwdx(:,:,2))
    endif
end if

! Top boundary condition
if (coord == nproc-1) then
    if (ubc_mom == 0) then
        ! Stress free
        RHSx(:,:,nz) = 0._rprec
        RHSy(:,:,nz) = 0._rprec
    else
        ! Wall
        ! if dudz, dvdz are on u-nodes for jz=nz, then we need a special
        ! definition of the vorticity in that case which also interpolates
        ! dwdx, dwdy to the u-node at jz=nz
        ! dwdy(1) and dwdx(1) should be 0
        RHSx(:,:,nz) = 0.5_rprec*(dwdy(:,:,nz-1) + dwdy(:,:,nz)) - dvdz(:,:,nz-1)
        RHSy(:,:,nz) = dudz(:,:,nz-1) - 0.5_rprec*(dwdx(:,:,nz-1) + dwdx(:,:,nz))
    endif
end if

! Transfer onto larger grid
call u_var%padd(u_big_var)
call v_var%padd(v_big_var)
call w_var%padd(w_big_var)
call RHSx_var%padd(vortx_big_var)
call RHSy_var%padd(vorty_big_var)
call RHSz_var%padd(vortz_big_var)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute (u x omega)_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc_big(:,:,0:nz-1) = -v_big(:,:,0:nz-1) * vortz_big(:,:,0:nz-1)                &
    + 0.5_rprec * (w_big(:,:,1:nz)*(vorty_big(:,:,1:nz))                       &
    + w_big(:,:,0:nz-1) * (vorty_big(:,:,0:nz-1)))

! Boundary conditions
if (coord == 0) then
    ! vort(1) is located on uvp-node
    ! the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    ! above the wall (could arguably be 0.25 * w(:,:,2))
    cc_big(:,:,1) = -v_big(:,:,1) * vortz_big(:,:,1)                           &
       + 0.5_rprec * w_big(:,:,2) * vorty_big(:,:,jzLo)
end if
if (coord == nproc-1) then
    ! vort(nz-1) is located on uvp-node
    ! the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    cc_big(:,:,nz-1) = -v_big(:,:,nz-1) * vortz_big(:,:,nz-1)                 &
       + 0.5_rprec * w_big(:,:,nz-1) * vorty_big(:,:,jzHi)
end if

! Dealias
call cc_big_var%unpadd(RHSx_var)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute (u x omega)_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc_big(:,:,0:nz-1) = u_big(:,:,0:nz-1) * vortz_big(:,:,0:nz-1)                &
    - 0.5_rprec * (w_big(:,:,1:nz)*(vortx_big(:,:,1:nz))                       &
    - w_big(:,:,0:nz-1) * (vortx_big(:,:,0:nz-1)))

! Boundary conditions
if (coord == 0) then
    ! vort(1) is located on uvp-node
    ! the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    ! above the wall (could arguably be 0.25 * w(:,:,2))
    cc_big(:,:,1) = u_big(:,:,1) * vortz_big(:,:,1)                           &
       - 0.5_rprec * w_big(:,:,2) * vortx_big(:,:,jzLo)
end if
if (coord == nproc-1) then
    ! vort(nz-1) is located on uvp-node
    ! the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    cc_big(:,:,nz-1) = u_big(:,:,nz-1) * vortz_big(:,:,nz-1)                 &
       - 0.5_rprec * w_big(:,:,nz-1) * vortx_big(:,:,jzHi)
end if

! Dealias
call cc_big_var%unpadd(RHSy_var)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute (u x omega)_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc_big(:,:,0:nz-1) = 0.5_rprec * (                                             &
    - (u_big(:,:,0:nz-1) + u_big(:,:,1:nz)) * vorty_big(:,:,0:nz-1)            &
    + (v_big(:,:,0:nz-1) + v_big(:,:,1:nz)) * vortx_big(:,:,0:nz-1) )

! Boundary conditions
! Some previous musings on this:
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
    !--earlier comment is also wrong, it is true that RHS = 0 at both walls and
    ! slip BC
if (coord == 0) then
    cc_big(:,:,1) = 0._rprec
end if
if (coord == nproc-1) then
    cc_big(:,:,nz-1) = 0._rprec
end if

! Dealias
call cc_big_var%unpadd(RHSz_var)

!

#ifdef PPSAFETYMODE
#ifdef PPMPI
RHSx(:, :, 0) = BOGUS
RHSy(:, :, 0) = BOGUS
RHSz(: ,:, 0) = BOGUS
#endif

!--top level is not valid
RHSx(:, :, nz) = BOGUS
RHSy(:, :, nz) = BOGUS
if (coord < nproc-1) RHSz(:, :, nz) = BOGUS
#endif

end subroutine convec
