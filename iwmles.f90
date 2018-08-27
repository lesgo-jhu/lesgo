!!
!!  Copyright (C) 2016  Johns Hopkins University
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
module iwmles
!*******************************************************************************
!
! This module contains the procedures for using the integral wall model.
! See Yang et al. 2015 for more details
!

use types, only : rprec

implicit none

private
public iwm_wallstress, iwm_init, iwm_finalize,                                 &
    iwm_checkpoint, iwm_read_checkpoint

! the instantaneous field
real(rprec), dimension(:,:), allocatable :: u_inst, v_inst, w_inst, p_inst
! u_tau,x  u_tau,y
real(rprec), dimension(:,:), allocatable :: iwm_utx, iwm_uty
real(rprec), dimension(:,:), allocatable :: iwm_utx_m, iwm_uty_m
real(rprec), dimension(:,:), allocatable :: iwm_utx_m2, iwm_uty_m2
! tau_wall,x tau_wall,y
real(rprec), dimension(:,:), allocatable :: iwm_tauwx, iwm_tauwy
! tau_top,x tau_top,y
real(rprec), dimension(:,:), allocatable :: iwm_tautx, iwm_tauty
! filtered tangential velocity, current and previous
real(rprec), dimension(:,:,:), allocatable :: iwm_flt_tagvel, iwm_flt_tagvel_m
! filtered pressure
real(rprec), dimension(:,:), allocatable :: iwm_flt_p
! direction x
integer :: iwm_dirx = 1
! direction y
integer :: iwm_diry = 2
! dimension of a surface, wall model always deal with 2D surfaces
! (because the world is 3D)
integer :: iwm_DN  = 2
! integrated profiles, current and previous
real(rprec), dimension(:,:,:), allocatable :: iwm_inte, iwm_inte_m
integer :: iwm_Lu = 1   ! index for integral of u
integer :: iwm_Luu = 2  ! index for integral of uu
integer :: iwm_Lv = 3   ! etc.
integer :: iwm_Lvv = 4
integer :: iwm_Luv = 5
! the total number of integrals that need to be calculated
integer :: iwm_LN  = 5
! unsteady, convective, pressure gradient
real(rprec), dimension(:,:,:), allocatable :: iwm_unsdy, iwm_conv, iwm_PrsGrad
! turbulent diffusion, RHS
real(rprec), dimension(:,:,:), allocatable :: iwm_diff, rhs, rhs_m, rhs_tot
! dudz at z=dz/2, dudz at z=zo
real(rprec), dimension(:,:,:), allocatable :: iwm_dudzT, iwm_dudzB
! filtered friction velocity, filtering time scale
real(rprec), dimension(:,:), allocatable :: iwm_flt_us, iwm_tR
! HALF cell height, zo, linear correction in x, y directions
real(rprec), dimension(:,:), allocatable :: iwm_Dz, iwm_Ax, iwm_Ay
real(rprec), dimension(:,:), allocatable :: iwm_z0

! number of time steps to skip between wall stress calculations
integer :: iwm_ntime_skip = 5
! time step size seen by the wall model
real(rprec) :: iwm_dt, iwm_dtp
real(rprec), dimension(:,:), allocatable :: delta_i_array
real(rprec) :: delta_i

real(rprec), dimension(:,:), allocatable :: u_sum,v_sum
integer :: time_steps

! flags for extraneous solutions of newton solver
integer, dimension(:,:), allocatable :: flag1, flag2, flag3, flag4, flag5, flag6
integer :: flag4_max
integer :: flag4_count

! storage of Ax,Ay before modifying to equilibrium solution
real(rprec), dimension(:,:), allocatable :: Axp, Ayp

! deverivatives in convective term
real(rprec), dimension(:,:), allocatable :: Luux, Lvvy, Luvx, Luvy, Lux, Lvy

contains

!*******************************************************************************
subroutine iwm_wallstress
!*******************************************************************************
use param, only : jt, nx, ny, dt
use sim_param , only : dudz, dvdz, txz, tyz
implicit none

integer :: iwm_i, iwm_j

! Calculate the time step used in the integral wall model
!! DO NOT USE iwm_ntime_skip=1 !! !! this number is hard coded to prevent any
!! mis-use...
if (mod(jt,iwm_ntime_skip)==1) then
    iwm_dt = dt
else
    iwm_dt = iwm_dt+dt
end if

! Compute the wall stress
if(mod(jt,iwm_ntime_skip)==0) then
    if (jt == iwm_ntime_skip) then
        iwm_dtp = iwm_dt
    end if
    ! gather flow status, update the integrated unsteady term, convective term,
    ! turbulent diffusion term etc.
    call iwm_calc_rhs()
    ! the subroutine to calculate wall stress
    call iwm_calc_wallstress()
    ! this is to monitor any quantity from the iwm, useful debugging tool
    call iwm_monitor()
    ! store previous time step for AB scheme
    iwm_dtp = iwm_dt
end if


! Imposing txz, tyz, dudz, dvdz every time step even iwm_* are not computed
! every time step.
do iwm_i = 1, nx
do iwm_j = 1, ny
    ! wall stress, use the value calculated in iwm, note the negative sign
    ! equilibrium wall model used in spanwise direction
    txz(iwm_i,iwm_j,1) = -iwm_tauwx(iwm_i,iwm_j)
    tyz(iwm_i,iwm_j,1) = -iwm_tauwy(iwm_i,iwm_j)

    ! Siang: I think those quantities should be consistent with the iwm.
    !        Users could switch back to equilibrium values, which I have tested
    !        and is fine.
    ! Note: Using the value calculated in iwm, note the positive sign
    dudz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
    dvdz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_diry)
end do
end do

end subroutine iwm_wallstress

!xiang:
!*******************************************************************************
subroutine iwm_init
!*******************************************************************************
!
! This subroutine allocates memory and initializes everything with plug flow
! conditions
!
use types, only : rprec
use param, only : nx, ny, ld, dz, vonk, zo, cfl, L_x, nu_molec, smooth

implicit none

real(rprec) :: usinit, uinit, vinit, Dzp

allocate(u_inst(ld,ny))
allocate(v_inst(ld,ny))
allocate(w_inst(ld,ny))
allocate(p_inst(ld,ny))

! at the height of the first grid point.
Dzp=dz/2._rprec

! initial value for us (the friction velocity)
usinit= 1._rprec
! initial value for the x-velocity at first grid point
if (smooth) then
    if (Dzp < 11._rprec*nu_molec/usinit) then
        uinit = usinit**2._rprec*Dzp/nu_molec
    else
        !uinit = usinit*(5._rprec + 1._rprec/vonk*log(usinit*dz/2._rprec/nu_molec))
        uinit = 1._rprec
    end if
else
    uinit = usinit/vonk*log(dz/2._rprec/zo)
end if
! initial value for the y-velocity at first grid point
vinit = 0._rprec

! us in x, y directions
allocate(iwm_utx(nx,ny))
allocate(iwm_uty(nx,ny))
allocate(iwm_utx_m(nx,ny))
allocate(iwm_uty_m(nx,ny))
iwm_utx = usinit
iwm_uty = 0._rprec
iwm_utx_m = iwm_utx
iwm_uty_m = iwm_uty

! wall stress in x, y directions
allocate(iwm_tauwx(nx,ny))
allocate(iwm_tauwy(nx,ny))
iwm_tauwx = usinit**2._rprec
iwm_tauwy = 0._rprec

allocate(iwm_tautx(nx,ny))
allocate(iwm_tauty(nx,ny))
if (smooth) then
    if (Dzp < 11._rprec*nu_molec/usinit) then
        iwm_tautx = usinit**2._rprec/nu_molec
        iwm_tauty = 0._rprec
    else
        iwm_tautx = 0._rprec
        iwm_tauty = 0._rprec
    end if
else
    iwm_tautx = 0._rprec
    iwm_tauty = 0._rprec
end if

! filtered velocity at the first grid point in x, y directions
allocate(iwm_flt_tagvel  (nx,ny,iwm_DN))
allocate(iwm_flt_tagvel_m(nx,ny,iwm_DN))
iwm_flt_tagvel  (:,:,iwm_dirx) = uinit
iwm_flt_tagvel  (:,:,iwm_diry) = vinit
iwm_flt_tagvel_m(:,:,iwm_dirx) = uinit
iwm_flt_tagvel_m(:,:,iwm_diry) = vinit

! pressure at first grid point
allocate(iwm_flt_p(nx,ny))
iwm_flt_p = 0._rprec

! integrals of Lu, Lv, etc.
allocate(iwm_inte  (nx,ny,iwm_LN))
allocate(iwm_inte_m(nx,ny,iwm_LN))
iwm_inte(:,:,iwm_Lu) = uinit*Dzp
iwm_inte(:,:,iwm_Lv) = 0._rprec
iwm_inte(:,:,iwm_Luu) = uinit**2._rprec*Dzp
iwm_inte(:,:,iwm_Lvv) = 0._rprec
iwm_inte(:,:,iwm_Luv) = 0._rprec
iwm_inte_m(:,:,iwm_Lu) = uinit*Dzp
iwm_inte_m(:,:,iwm_Lv) = 0._rprec
iwm_inte_m(:,:,iwm_Luu) = uinit**2._rprec*Dzp
iwm_inte_m(:,:,iwm_Lvv) = 0._rprec
iwm_inte_m(:,:,iwm_Luv) = 0._rprec

! each term in the integral equation and top/bottom derivatives
allocate(iwm_unsdy  (nx,ny,iWM_DN))
allocate(iwm_conv   (nx,ny,iWM_DN))
allocate(iwm_PrsGrad(nx,ny,iWM_DN))
allocate(iwm_diff   (nx,ny,iwm_DN))
allocate(rhs        (nx,ny,iWM_DN))
allocate(rhs_m      (nx,ny,iWM_DN))
allocate(rhs_tot    (nx,ny,iWM_DN))
allocate(iwm_dudzT  (nx,ny,iwm_DN))
allocate(iwm_dudzB  (nx,ny,iwm_DN))
iWM_unsdy   = 0._rprec
iWM_conv    = 0._rprec
iWM_PrsGrad = 0._rprec
iwm_diff(:,:,iwm_dirx) = iwm_tautx - iwm_tauwx
iwm_diff(:,:,iwm_diry) = iwm_tauty - iwm_tauwy
rhs         = Dzp
rhs_m       = rhs
rhs_tot     = 0._rprec
if (smooth) then
    if (Dzp < 11._rprec*nu_molec/usinit) then
        iwm_dudzT(:,:,iwm_dirx) = usinit**2._rprec/nu_molec
        iwm_dudzB(:,:,iwm_dirx) = usinit**2._rprec/nu_molec
    else
        iwm_dudzT(:,:,iwm_dirx) = usinit/vonk/Dzp
        iwm_dudzB(:,:,iwm_dirx) = usinit**2._rprec/nu_molec
    end if
else
    iwm_dudzT(:,:,iwm_dirx) = usinit/vonk/Dzp
    iwm_dudzB(:,:,iwm_dirx) = usinit/vonk/zo
end if
iwm_dudzT(:,:,iwm_diry) = 0._rprec
iwm_dudzB(:,:,iwm_diry) = 0._rprec

! filtered friction velocity and the filtering time scale, tR<1
allocate(iwm_flt_us(nx,ny))
allocate(iwm_tR    (nx,ny))
iwm_flt_us = usinit
iwm_tR = (cfl*L_x/nx/uinit)/(dz/2._rprec/vonk/usinit)

! cell height and imposed roughness length
allocate(iwm_Dz(nx,ny))
allocate(iwm_z0(nx,ny))
iwm_Dz = dz/2._rprec
iwm_z0 = zo !we leave the possibility of zo as a function of x-y

! linear correction to the log profile
allocate(iwm_Ax(nx,ny))
allocate(iwm_Ay(nx,ny))
iwm_Ax = 0._rprec
iwm_Ay = 0._rprec

! time step seen by the iwm
iwm_dt=iwm_ntime_skip*cfl*L_x/nx/uinit

allocate(delta_i_array(nx,ny))
delta_i_array = dz/2._rprec

allocate(u_sum(nx,ny))
allocate(v_sum(nx,ny))
u_sum = 0._rprec
v_sum = 0._rprec

time_steps = 0

allocate(flag1(nx,ny))
allocate(flag2(nx,ny))
allocate(flag3(nx,ny))
allocate(flag4(nx,ny))
allocate(flag5(nx,ny))
allocate(flag6(nx,ny))
flag1 = 0
flag2 = 0
flag3 = 0
flag4 = 0
flag5 = 0
flag6 = 0

allocate(Axp(nx,ny))
allocate(Ayp(nx,ny))
Axp = 0._rprec
Ayp = 0._rprec

allocate(Luux(nx,ny))
allocate(Lvvy(nx,ny))
allocate(Luvx(nx,ny))
allocate(Luvy(nx,ny))
allocate(Lux(nx,ny))
allocate(Lvy(nx,ny))
Luux = 0._rprec
Lvvy = 0._rprec
Luvx = 0._rprec
Luvy = 0._rprec
Lux = 0._rprec
Lvy = 0._rprec

end subroutine iwm_init

!*******************************************************************************
subroutine iwm_finalize
!*******************************************************************************
!
! This subroutine deallocates memory used for iwm
!
implicit none

deallocate(u_inst)
deallocate(v_inst)
deallocate(w_inst)
deallocate(p_inst)

deallocate(iwm_utx)
deallocate(iwm_uty)
deallocate(iwm_utx_m)
deallocate(iwm_uty_m)

deallocate(iwm_tauwx)
deallocate(iwm_tauwy)

deallocate(iwm_tautx)
deallocate(iwm_tauty)

deallocate(iwm_flt_tagvel  )
deallocate(iwm_flt_tagvel_m)

deallocate(iwm_inte  )
deallocate(iwm_inte_m)

deallocate(iwm_unsdy  )
deallocate(iwm_conv   )
deallocate(iwm_PrsGrad)
deallocate(iwm_diff   )
deallocate(rhs        )
deallocate(rhs_m      )
deallocate(rhs_tot    )
deallocate(iwm_dudzT  )
deallocate(iwm_dudzB  )

deallocate(iwm_flt_us)
deallocate(iwm_tR    )

deallocate(iwm_Dz)
deallocate(iwm_z0)
deallocate(iwm_Ax)
deallocate(iwm_Ay)
deallocate(Axp)
deallocate(Ayp)

deallocate(flag1)
deallocate(flag2)
deallocate(flag3)
deallocate(flag4)
deallocate(flag5)
deallocate(flag6)

deallocate(Luux)
deallocate(Lvvy)
deallocate(Luvx)
deallocate(Luvy)
deallocate(Lux)
deallocate(Lvy)

end subroutine iwm_finalize


!*******************************************************************************
subroutine iwm_calc_rhs()
!*******************************************************************************
!
! Ths subroutine calculates the right hand side of the iwm system.
!
use grid_m, only : grid
use types,only : rprec
use param,only : nx,ny,dx,dy,ld,coord
use sim_param,only : u,v,w,p
use test_filtermodule
implicit none

! useful array for autowraped index
integer, pointer, dimension(:) :: autowrap_i, autowrap_j
integer :: iwm_i,iwm_j
! the mean pressure at first grid point
real(rprec) :: p_bar
! for temporary storage of derivativex of integrals like dLudx, dLvdx...
real(rprec) :: phip, phim

nullify(autowrap_i, autowrap_j)
autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

! update the u, v for previous time step
do iwm_i = 1, nx
do iwm_j = 1, ny
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx) =                                   &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry) =                                   &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
end do
end do

! get the instantaneous field
u_inst = u(:,:,1)
! let us do not worry about the half cell displacement in v
v_inst = v(:,:,1)
! w is quadrantic near the wall
w_inst = w(:,:,2)*0.25_rprec

! the real pressure is needed, this step is CODE SPECIFIC!
p_inst = p(:,:,1)
do iwm_i = 1, nx
do iwm_j = 1, ny
    p_inst(iwm_i,iwm_j)= p_inst(iwm_i,iwm_j)                                   &
        -0.5_rprec*( (u_inst(iwm_i,iwm_j))**2._rprec                           &
        + (v_inst(iwm_i,iwm_j))**2._rprec                                      &
        +(w_inst(iwm_i,iwm_j))**2._rprec   )
end do
end do

! obtain the pressure fluctuations
p_bar = 0._rprec
do iwm_i = 1, nx
do iwm_j = 1, ny
    p_bar = p_bar+p_inst(iwm_i,iwm_j)
end do
end do
p_bar=p_bar/nx/ny

do iwm_i = 1, nx
do iwm_j = 1, ny
    p_inst(iwm_i,iwm_j) = p_inst(iwm_i,iwm_j)-p_bar
end do
end do

! all the data enters must be filtered (see Anderson & Meneveau 2011 JFM)
call test_filter ( u_inst )
call test_filter ( v_inst )
call test_filter ( w_inst )
call test_filter ( p_inst )

!time_steps = time_steps + 1;
!temporal filtering
do iwm_i=1,nx
do iwm_j=1,ny
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx) =                                     &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)*(1._rprec-iwm_tR(iwm_i,iwm_j))    &
        + u_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) =                                     &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)*(1._rprec-iwm_tR(iwm_i,iwm_j))    &
        + v_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
    iwm_flt_p(iwm_i,iwm_j) = iwm_flt_p(iwm_i,iwm_j)                            &
        * (1._rprec-iwm_tR(iwm_i,iwm_j))                                       &
        + p_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
end do
end do

!write(*,*) iwm_flt_tagvel(int(nx/2._rprec),int(ny/2._rprec),iwm_dirx),&
!           iwm_flt_tagvel(int(nx/2._rprec),int(ny/2._rprec),iwm_diry)

! calculate RHS, calculation of the integrals is done from the last time step
! in the subroutine iwm_calc_wallstress, so is iwm_diff
do iwm_i = 1, nx
do iwm_j = 1, ny
    ! store previous LHS calculation for Adams-Bashforth solver
    rhs_m(iwm_i,iwm_j,iwm_dirx) = rhs(iwm_i,iwm_j,iwm_dirx)
    rhs_m(iwm_i,iwm_j,iwm_diry) = rhs(iwm_i,iwm_j,iwm_diry)

    ! the unsteady term
    iwm_unsdy(iwm_i,iwm_j,iwm_dirx) =                                          &
        (iwm_inte(iwm_i,iwm_j,iwm_Lu)-iwm_inte_m(iwm_i,iwm_j,iwm_Lu))/iwm_dt
    iwm_unsdy(iwm_i,iwm_j,iwm_diry) =                                          &
        (iwm_inte(iwm_i,iwm_j,iwm_Lv)-iwm_inte_m(iwm_i,iwm_j,iwm_Lv))/iwm_dt

    ! the convective term
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luu)
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luu)
    Luux(iwm_i,iwm_j) = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Luv)
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Luv)
    Luvy(iwm_i,iwm_j) = (phip-phim)/dy/2._rprec
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luv)
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luv)
    Luvx(iwm_i,iwm_j) = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lvv)
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lvv)
    Lvvy(iwm_i,iwm_j) = (phip-phim)/dy/2._rprec
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Lu )
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Lu )
    Lux(iwm_i,iwm_j) = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lv )
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lv )
    Lvy(iwm_i,iwm_j) = (phip-phim)/dy/2._rprec
    iwm_conv(iwm_i,iwm_j,iwm_dirx) = Luux(iwm_i,iwm_j) + Luvy(iwm_i,iwm_j)     &
        -iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_dirx)*(Lux(iwm_i,iwm_j)+Lvy(iwm_i,iwm_j))
    iwm_conv(iwm_i,iwm_j,iwm_diry) = Luvx(iwm_i,iwm_j) + Lvvy(iwm_i,iwm_j)     &
        -iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_diry)*(Lux(iwm_i,iwm_j)+Lvy(iwm_i,iwm_j))

    ! the pressure gradient term
    phip = iwm_flt_p(autowrap_i(iwm_i+1),iwm_j)
    phim = iwm_flt_p(autowrap_i(iwm_i-1),iwm_j)
    ! including the mean unit pressure gradient
    iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx) = (phip-phim)/dx/2._rprec                &
        * iwm_Dz(iwm_i,iwm_j) - 1._rprec*iwm_Dz(iwm_i,iwm_j)
    phip = iwm_flt_p(iwm_i,autowrap_j(iwm_j+1))
    phim = iwm_flt_p(iwm_i,autowrap_j(iwm_j-1))
    iwm_PrsGrad(iwm_i,iwm_j,iwm_diry) = (phip-phim)/dy/2._rprec                &
        * iwm_Dz(iwm_i,iwm_j)

    ! the right hand side
    ! this is the integrated momentum equation before AB
    rhs(iwm_i,iwm_j,iwm_dirx) = -iwm_conv(iwm_i,iwm_j,iwm_dirx)                &
        - iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)                                    &
        + iwm_diff(iwm_i,iwm_j,iwm_dirx)
    rhs(iwm_i,iwm_j,iwm_diry) = -iwm_conv(iwm_i,iwm_j,iwm_diry)                &
        - iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)                                    &
        + iwm_diff(iwm_i,iwm_j,iwm_diry)

    ! this is the integrated momentum equation after AB
    rhs_tot(iwm_i,iwm_j,iwm_dirx) =  iwm_inte(iwm_i,iwm_j,iwm_Lu)              &
        + iwm_dt*( (1._rprec+iwm_dt/iwm_dtp/2._rprec)*rhs(iwm_i,iwm_j,iwm_dirx)&
                 - iwm_dt/iwm_dtp/2._rprec*rhs_m(iwm_i,iwm_j,iwm_dirx) )
    rhs_tot(iwm_i,iwm_j,iwm_diry) =  iwm_inte(iwm_i,iwm_j,iwm_Lv)              &
        + iwm_dt*( (1._rprec+iwm_dt/iwm_dtp/2._rprec)*rhs(iwm_i,iwm_j,iwm_diry)&
                 - iwm_dt/iwm_dtp/2._rprec*rhs_m(iwm_i,iwm_j,iwm_diry) )
end do
end do

nullify(autowrap_i, autowrap_j)

end subroutine iwm_calc_rhs

!*******************************************************************************
subroutine iwm_slv(rhsx,rhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy,Ax,Ay)
!*******************************************************************************
use types, only : rprec
use param, only : vonk, nu_molec, smooth, nx, ny
implicit none

real(rprec), intent(in)  :: rhsx, rhsy, Ux, Uy, Dz, utx, uty, z0
real(rprec), intent(out) :: fx,fy
real(rprec), intent(out) :: Ax, Ay
real(rprec) :: Lu, Lv, Cx, Cy, delta_nu
real(rprec) :: utau, vel

if (smooth) then

    utau = (utx**4.0 + uty**4.0)**0.25
    vel = sqrt(Ux**2.0+Uy**2.0)

    ! Eq. C25 in Yang et al. 2015
    delta_i = 11._rprec*nu_molec/utau

    ! Eq. C26 in Yang et al. 2015
    delta_nu = nu_molec/((utx**4.0 + uty**4.0)/(utx**2.0 + uty**2.0))**(1.0/2.0)

    ! Eq. C27 in Yang et al. 2015 (modified)
    Ax = (vel/utau + 1._rprec/vonk*log(delta_i/Dz)                            &
        - utx/utau*delta_i/delta_nu*vel/Ux)                                   &
        /(1._rprec - delta_i/delta_nu)
    Ay = (vel/utau + 1._rprec/vonk*log(delta_i/Dz)                            &
        - uty/utau*delta_i/delta_nu*vel/Uy)                                   &
        /(1._rprec - delta_i/delta_nu)

    ! Eq. C28 in Yang et al. 2015 (modified)
    Cx = vel/utau - Ax
    Cy = vel/utau - Ay

    ! Eq. E27 in Yang et al. 2015 (modified)
    Lu = 1._rprec/2._rprec*utx*delta_i**2.0/delta_nu + utau*Dz                &
        *(1._rprec/2._rprec*Ux/vel*Ax*(1._rprec - (delta_i/Dz)**2.0)          &
        + Ux/vel*Cx*(1._rprec - delta_i/Dz)                                   &
        - Ux/vel/vonk*(1._rprec - delta_i/Dz + delta_i/Dz*log(delta_i/Dz)))
    Lv = 1._rprec/2._rprec*uty*delta_i**2.0/delta_nu + utau*Dz                &
        *(1._rprec/2._rprec*Uy/vel*Ay*(1._rprec - (delta_i/Dz)**2.0)          &
        + Uy/vel*Cy*(1._rprec - delta_i/Dz)                                   &
        - Uy/vel/vonk*(1._rprec - delta_i/Dz + delta_i/Dz*log(delta_i/Dz)))

else ! rough wall

    Ax = (Ux/utx-1.0/vonk*log(Dz/z0))/(1.0-z0/Dz)
    Ay = (Uy/uty-1.0/vonk*log(Dz/z0))/(1.0-z0/Dz)

    Lu = 1.0/2.0*utx*Dz*Ax*(1.0-z0/Dz)**2.0+1.0/vonk*utx*Dz*(z0/Dz-1.0+log(Dz/z0))
    Lv = 1.0/2.0*uty*Dz*Ay*(1.0-z0/Dz)**2.0+1.0/vonk*uty*Dz*(z0/Dz-1.0+log(Dz/z0))

end if

fx = Lu-rhsx
fy = Lv-rhsy

end subroutine iwm_slv

!*******************************************************************************
subroutine visc_slv(Ux,Uy,Dz,utx,uty,fx2,fy2)
!*******************************************************************************

use types, only : rprec
use param, only : nu_molec
implicit none

real(rprec), intent(in)  :: Ux, Uy, Dz, utx, uty
real(rprec), intent(out) :: fx2,fy2
real(rprec) :: delta_nu

delta_i = Dz
delta_nu = nu_molec/((utx**4.0 + uty**4.0)/(utx**2.0 + uty**2.0))**0.5

fx2 = Ux - utx*delta_i/delta_nu
fy2 = Uy - uty*delta_i/delta_nu


end subroutine 

!*******************************************************************************
subroutine iwm_calc_wallstress
!*******************************************************************************
use types, only : rprec
use param, only : vonk, nx, ny, nu_molec, smooth, coord, jt
use test_filtermodule

implicit none

integer :: iwm_i, iwm_j
real(rprec) :: fx, fy, fxp, fyp
real(rprec) :: fx2, fy2, fx2p, fy2p
real(rprec) :: iwm_tol, iwm_eps
real(rprec) :: a11, a12, a21, a22
real(rprec) :: iwmutxP,iwmutyP
integer :: iter, MaxIter
real(rprec) :: equilWMpara,equilutx,equiluty
real(rprec) :: Ax, Ay, Cx, Cy, z0, Dz
real(rprec) :: Ux, Uy, dVelzT, dVelzB, Vel
real(rprec) :: utx, uty, utx_m, uty_m, utxp, utyp, utau, delta_nu
real(rprec) :: rhsx, rhsy
real(rprec) :: Lu, Lv
integer :: equil_flag

MaxIter=1500

iwm_tol = 0.000001_rprec
iwm_eps = 0.000000001_rprec

flag4_max = 0
flag4_count = 0

do iwm_i=1,nx
do iwm_j=1,ny

    ! Rename for simplification
    Ux = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
    Uy = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
    utx = iwm_utx(iwm_i,iwm_j)
    uty = iwm_uty(iwm_i,iwm_j)
    utx_m = iwm_utx_m(iwm_i,iwm_j)
    uty_m = iwm_uty_m(iwm_i,iwm_j)
    rhsx = rhs_tot(iwm_i,iwm_j,iwm_dirx)
    rhsy = rhs_tot(iwm_i,iwm_j,iwm_diry)
    Dz = iwm_Dz(iwm_i,iwm_j)
    z0 = iwm_z0(iwm_i,iwm_j)

    !initial guess before solver
    utx = utx_m
    uty = uty_m
!    utx = 1._rprec*sign(1._rprec,Ux)
!    uty = 0.1_rprec*sign(1._rprec,Uy)


    iter = 0
    equil_flag = 0
    flag1(iwm_i,iwm_j) = 0
    flag2(iwm_i,iwm_j) = 0
    flag3(iwm_i,iwm_j) = 0
    flag4(iwm_i,iwm_j) = 0
    flag5(iwm_i,iwm_j) = 0
    flag6(iwm_i,iwm_j) = 0

    call iwm_slv(rhsx,rhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy,Ax,Ay)

    ! Check if in viscous sublayer and solve if it is
    if (Dz <= 11._rprec*nu_molec/utx) then
        call visc_slv(Ux,Uy,Dz,utx,uty,fx2,fy2)
        do while (max(abs(fx2),abs(fy2))>iwm_tol)
            utxp = utx+iwm_eps
            utyp = uty
            call visc_slv(Ux,Uy,Dz,utxp,utyp,fx2p,fy2p)
            a11 = (fx2p-fx2)/iwm_eps
            a21 = (fy2p-fy2)/iwm_eps
            utxp = utx
            utyp = uty+iwm_eps
            call visc_slv(Ux,Uy,Dz,utxp,utyp,fx2p,fy2p)
            a12 = (fx2p-fx2)/iwm_eps
            a22 = (fy2p-fy2)/iwm_eps
            utx = utx - 0.5_rprec*( a22*fx2-a12*fy2)/(a11*a22-a12*a21)
            uty = uty - 0.5_rprec*(-a21*fx2+a11*fy2)/(a11*a22-a12*a21)
            call visc_slv(Ux,Uy,Dz,utx,uty,fx2,fy2)
            iter = iter+1
            if (iter>MaxIter) then
                flag1(iwm_i,iwm_j) = 1
                equil_flag = 1
                exit
            end if
        end do
        fx = 0._rprec
        fy = 0._rprec
    end if

    ! use Newton method to solve the system
    do while (max(abs(fx),abs(fy))>iwm_tol)
        utxp=utx+iwm_eps
        utyp=uty
        call iwm_slv(rhsx,rhsy,Ux,Uy,Dz,z0,utxp,utyp,fxp,fyp,Ax,Ay)
        a11 = (fxp-fx)/iwm_eps
        a21 = (fyp-fy)/iwm_eps
        utxp = utx
        utyp = uty+iwm_eps
        call iwm_slv(rhsx,rhsy,Ux,Uy,Dz,z0,utxp,utyp,fxp,fyp,Ax,Ay)
        a12 = (fxp-fx)/iwm_eps
        a22 = (fyp-fy)/iwm_eps
        utx = utx - 0.50*( a22*fx-a12*fy)/(a11*a22-a12*a21)
        uty = uty - 0.50*(-a21*fx+a11*fy)/(a11*a22-a12*a21)
        call iwm_slv(rhsx,rhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy,Ax,Ay)
        iter = iter+1
        ! maximum iteration reached
        if (iter>MaxIter) then
            flag1(iwm_i,iwm_j) = 1
            equil_flag = 1
            exit
        end if
    end do

    if (iwm_i == int(nx/2._rprec) .and. iwm_j == int(ny/2._rprec)) then
        write(*,*) utx, uty, Ax, iter
    end if

    ! store Ax,Ay before adjusting to equilibrium solution
    Axp(iwm_i,iwm_j) = Ax
    Ayp(iwm_i,iwm_j) = Ay

    ! large utx or uty
    if (abs(utx) > 100._rprec .or. abs(uty) > 100._rprec) then
        flag2(iwm_i,iwm_j) = 1
        equil_flag = 1
    end if

    ! large Ax or Ay
    if (abs(Ax) > 10._rprec .or. abs(Ay) > 10._rprec) then
        flag3(iwm_i,iwm_j) = 1
        equil_flag = 1
    end if

    ! infinity check
    if (utx - 1.0 == utx .or. uty - 1.0 == uty) then
        flag5(iwm_i,iwm_j) = 1
        equil_flag = 1
    end if

    ! NaN check
    if (utx /= utx .or. uty /= uty) then
        flag6(iwm_i,iwm_j) = 1
        equil_flag = 1
    end if

    ! calculate equilibrium us
    if (equil_flag /=0) then
        if (smooth) then
             if (Dz < 11._rprec*nu_molec) then
                 utx = sqrt(abs(Ux)*nu_molec/Dz)*sign(1._rprec,Ux)
                 uty = Uy*utx/Ux 
             else
                 delta_i = 11._rprec*nu_molec/1._rprec
                 z0 = delta_i*exp(-11._rprec*vonk)
                 utx = vonk*Ux/log(Dz/z0)
                 uty = Uy*utx/Ux
             end if
        else
            utx = vonk*Ux/log(Dz/z0)
            uty = vonk*Uy/log(Dz/z0)
        end if
    end if

    ! Total quantities
    utau = (utx**4._rprec + uty**4._rprec)**0.25_rprec
    Vel = sqrt(Ux**2.0 + Uy**2.0)

    if (smooth) then
        ! Eq. C25 in Yang et al. 2015
        if(Dz < 11._rprec*nu_molec/utau) then
            delta_i = Dz
        else
            delta_i = 11._rprec*nu_molec/utau
        end if
        ! Eq. C26 in Yang et al. 2015
        delta_nu = nu_molec/((utx**4.0 + uty**4.0)                         &
            /(utx**2.0 + uty**2.0))**(1.0/2.0)
    end if

    !calculate Ax, Ay
    if (smooth) then
        if (delta_i == Dz) then
            Ax = 0._rprec
            Ay = 0._rprec
        else
            ! Eq. C27 in Yang et al. 2015 (modified)
            Ax = (vel/utau + 1._rprec/vonk*log(delta_i/Dz)                            &
                - utx/utau*delta_i/delta_nu*vel/Ux)                                   &
                /(1._rprec - delta_i/delta_nu)
            Ay = (vel/utau + 1._rprec/vonk*log(delta_i/Dz)                            &
                - uty/utau*delta_i/delta_nu*vel/Uy)                                   &
                /(1._rprec - delta_i/delta_nu)
        end if
    else
        ! eq. D2 in Yang et al. 2015
        Ax = (Ux/utx-1._rprec/vonk*log(Dz/z0))/(1._rprec-z0/Dz)
        Ay = (Uy/uty-1._rprec/vonk*log(Dz/z0))/(1._rprec-z0/Dz)
    end if
    
    ! check if in viscous sublayer
    if(Dz <= delta_i) then
        flag4(iwm_i,iwm_j) = 1
        flag4_count = flag4_count + 1
    end if
        
    if (flag4(iwm_i,iwm_j) > flag4_max) then
        flag4_max = flag4(iwm_i,iwm_j)
    end if

    ! store the friction velocities last time step
    iwm_utx_m(iwm_i,iwm_j) = iwm_utx(iwm_i,iwm_j)
    iwm_uty_m(iwm_i,iwm_j) = iwm_uty(iwm_i,iwm_j)

    ! store the friction velocities
    iwm_utx(iwm_i,iwm_j) = utx
    iwm_uty(iwm_i,iwm_j) = uty

    ! store the linear correction
    iwm_Ax(iwm_i,iwm_j) = Ax
    iwm_Ay(iwm_i,iwm_j) = Ay

    ! store delta_i_array
    delta_i_array(iwm_i,iwm_j) = delta_i

    ! update integral for last time step
    iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu)
    iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv)
    iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv)
    iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu)
    iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv)

    ! calculate the needed integrals
    if (smooth) then
        
        if (delta_i == Dz) then !smooth wall, grid point inside viscous sublayer

        iwm_inte(iwm_i,iwm_j,iwm_Lu) = 1._rprec/2._rprec*utx                   &
            *delta_i**2.0/delta_nu
        iwm_inte(iwm_i,iwm_j,iwm_Lv) = 1._rprec/2._rprec*uty                   &
            *delta_i**2.0/delta_nu
        iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1._rprec/3._rprec*utx*uty              &
            *delta_i**3.0/delta_nu**2.0
        iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1._rprec/3._rprec*utx**2.0             &
            *delta_i**3.0/delta_nu**2.0
        iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1._rprec/3._rprec*uty**2.0             &
            *delta_i**3.0/delta_nu**2.0
       
        ! Velocity derivatives at top of profile 
        iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = utx/delta_nu
        iwm_dudzT(iwm_i,iwm_j,iwm_diry) = uty/delta_nu

        ! Calculate top shear stress
        ! Mixing length is zero in viscous sublayer
        iwm_tautx(iwm_i,iwm_j) = nu_molec*iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
        iwm_tauty(iwm_i,iwm_j) = nu_molec*iwm_dudzT(iwm_i,iwm_j,iwm_diry)
        
        ! Calculate wall stress
        iwm_tauwx(iwm_i,iwm_j) = utx**2.0*sign(1._rprec,utx)
        iwm_tauwy(iwm_i,iwm_j) = uty**2.0*sign(1._rprec,uty)

        ! Calculate turbulent diffusion (difference in shear stress)
        ! Eq. C15 in Yang et al. 2015
        iwm_diff(iwm_i,iwm_j,iwm_dirx) = iwm_tautx(iwm_i,iwm_j) - iwm_tauwx(iwm_i,iwm_j)
        iwm_diff(iwm_i,iwm_j,iwm_diry) = iwm_tautx(iwm_i,iwm_j) - iwm_tauwy(iwm_i,iwm_j)
        
        else !smooth wall, grid point outside viscous sublayer

        ! Eq. C28 in Yang et al. 2015 (modified)
        Cx = vel/utau - Ax
        Cy = vel/utau - Ay

        ! Eq. E27 in Yang et al. 2015 (modified)
        iwm_inte(iwm_i,iwm_j,iwm_Lu) =                                            &
            1._rprec/2._rprec*utx*delta_i**2.0/delta_nu + utau*Dz                 &
            *(1._rprec/2._rprec*Ux/vel*Ax*(1._rprec - (delta_i/Dz)**2.0)          &
            + Ux/vel*Cx*(1._rprec - delta_i/Dz)                                   &
            - Ux/vel/vonk*(1._rprec - delta_i/Dz + delta_i/Dz*log(delta_i/Dz)))
        iwm_inte(iwm_i,iwm_j,iwm_Lv) =                                            &
            1._rprec/2._rprec*uty*delta_i**2.0/delta_nu + utau*Dz                 &
            *(1._rprec/2._rprec*Uy/vel*Ay*(1._rprec - (delta_i/Dz)**2.0)          &
            + Uy/vel*Cy*(1._rprec - delta_i/Dz)                                   &
            - Uy/vel/vonk*(1._rprec - delta_i/Dz + delta_i/Dz*log(delta_i/Dz)))

        ! Eq. C21 in Yang et al. 2015 (modified)
        iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1._rprec/3._rprec*utx*uty              &
            *delta_i**3.0/delta_nu**2.0 + utau**2.0*Dz*(                       &
            - Ux*Uy/vel**2.0/vonk*(Ax + Ay)                                    &
            *(1._rprec/4._rprec-1._rprec/4._rprec*(delta_i/Dz)**2.0            &
            + 1._rprec/2._rprec*(delta_i/Dz)**2.0*log(delta_i/Dz))             &
            - Ux*Uy/vel**2.0/vonk*(Cx + Cy)*(1._rprec - delta_i/Dz             &
            + delta_i/Dz*log(delta_i/Dz)) - Ux*Uy/Vel**2.0/vonk**2.0           &
            *(delta_i/Dz - 2._rprec + delta_i/Dz*(log(delta_i/Dz) - 1._rprec)**2.0)&
            + Ux*Uy/vel**2.0/3._rprec*Ax*Ay*(1._rprec - (delta_i/Dz)**3.0)     &
            + Ux*Uy/vel**2.0/2._rprec*(Ax*Cy + Ay*Cx)*(1._rprec - (delta_i/Dz)**2.0)&
            + Ux*Uy/vel**2.0*Cx*Cy*(1._rprec - delta_i/Dz))

        ! Eq. C20 in Yang et al. 2015 (modified)
        iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1._rprec/3._rprec*utx**2.0             &
            *delta_i**3.0/delta_nu**2.0 + utau**2.0*Dz                         &
            *(-Ax*Ux**2.0/vel**2.0/vonk*(delta_i/Dz)**2.0*log(delta_i/Dz)      &
            + Ux/vel*Ax*(Ux/vel*Cx - 1._rprec/2._rprec/vonk*Ux/vel)            &
            *(1._rprec - (delta_i/Dz)**2.0)                                    &
            + (Ux/vel*Ax)**2.0/3._rprec*(1._rprec - (delta_i/Dz)**3.0)         &
            + (Ux/vel*Cx - 1._rprec/vonk*Ux/vel)**2.0 - delta_i/Dz             &
            *(Ux/vel*Cx - Ux/vel/vonk + Ux/vel/vonk*log(delta_i/Dz))**2.0      &
            + (Ux/vel/vonk)**2.0*(1._rprec - delta_i/Dz) )
        iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1._rprec/3._rprec*uty**2.0             &
            *delta_i**3.0/delta_nu**2.0 + utau**2.0*Dz                         &
            *(-Ay*Uy**2.0/vel**2.0/vonk*(delta_i/Dz)**2.0*log(delta_i/Dz)      &
            + Uy/vel*Ay*(Uy/vel*Cy - 1._rprec/2._rprec/vonk*Uy/vel)            &
            *(1._rprec - (delta_i/Dz)**2.0)                                    &
            + (Uy/vel*Ay)**2.0/3._rprec*(1._rprec - (delta_i/Dz)**3.0)         &
            + (Uy/vel*Cy - 1._rprec/vonk*Uy/vel)**2.0 - delta_i/Dz             &
            *(Uy/vel*Cy - Uy/vel/vonk + Uy/vel/vonk*log(delta_i/Dz))**2.0      &
            + (Uy/vel/vonk)**2.0*(1._rprec - delta_i/Dz) )

        iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0._rprec
        iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 0._rprec
        iwm_inte(iwm_i,iwm_j,iwm_Luv) = 0._rprec

        ! Calculate top derivatives
        ! Eq. C14 in Yang et al 2015 (modified)
        iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = utau*Ux/vel/Dz*(1._rprec/vonk + Ax)
        iwm_dudzT(iwm_i,iwm_j,iwm_diry) = utau*Uy/vel/Dz*(1._rprec/vonk + Ay)
        dVelzT = abs(Ux/Vel*iwm_dudzT(iwm_i,iwm_j,iwm_dirx)                    &
            + Uy/Vel*iwm_dudzT(iwm_i,iwm_j,iwm_diry))

        ! Calculate the stress at the top
        ! Eq. C12 and C13 in Yang et al 2015
        iwm_tautx(iwm_i,iwm_j) = abs((nu_molec + (vonk*Dz)**2._rprec*dVelzT)   &
            *iwm_dudzT(iwm_i,iwm_j,iwm_dirx))*sign(1._rprec,utx)
        iwm_tauty(iwm_i,iwm_j) = abs((nu_molec + (vonk*Dz)**2._rprec*dVelzT)   &
            *iwm_dudzT(iwm_i,iwm_j,iwm_diry))*sign(1._rprec,uty)

        ! Calculate wall stress
        iwm_tauwx(iwm_i,iwm_j) = utx**2.0*sign(1._rprec,utx)
        iwm_tauwy(iwm_i,iwm_j) = uty**2.0*sign(1._rprec,uty)

        ! Calculate turbulent diffusion
        iwm_diff(iwm_i,iwm_j,iwm_dirx) = iwm_tautx(iwm_i,iwm_j) - iwm_tauwx(iwm_i,iwm_j)
        iwm_diff(iwm_i,iwm_j,iwm_diry) = iwm_tauty(iwm_i,iwm_j) - iwm_tauwy(iwm_i,iwm_j)
        
        end if

    else ! rough wall

        ! Eq. D7 in Yang et al. 2015
        iwm_inte(iwm_i,iwm_j,iwm_Lu) = 1._rprec/2._rprec*Dz*Ax                 &
            * (1._rprec - z0/Dz)**2.0 + 1._rprec/vonk*utx*Dz                   &
            * (z0/Dz - 1._rprec + log(Dz/z0))
        iwm_inte(iwm_i,iwm_j,iwm_Lv) = 1._rprec/2._rprec*Dz*Ay                 &
            * (1._rprec - z0/Dz)**2.0 + 1._rprec/vonk*uty*Dz                   &
            * (z0/Dz - 1._rprec + log(Dz/z0))

        ! Eq. D8 in Yang et al 2015
        iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1._rprec/vonk**2.0*utx*uty*Dz          &
            * (1._rprec - 2*z0/Dz + (1._rprec - log(Dz/z0))**2.0)              &
            + 1._rprec/3._rprec*Ax*Ay*Dz*(1._rprec - z0/Dz)**3.0               &
            - 1._rprec/4._rprec/vonk*(Ax*uty + Ay*utx)*Dz                      &
            * (1._rprec - 4._rprec*z0/Dz + 3._rprec*z0**2.0/Dz**2.0            &
            - 2._rprec*log(Dz/z0) + 4._rprec*z0/Dz*log(Dz/z0))

        ! Eq. D9 in Yang et al 2015
        iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1._rprec/vonk**2.0*utx**2.0*Dz         &
            * ((log(Dz/z0) - 1._rprec)**2.0 - 2._rprec*z0/Dz + 1._rprec)       &
            + 1._rprec/3._rprec*Ax**2.0*Dz*(1._rprec - z0/Dz)**3.0             &
            - 1._rprec/2._rprec/vonk*utx*Ax*Dz                                 &
            * (1._rprec - 4._rprec*z0/Dz+3._rprec*z0**2.0/Dz**2.0              &
            - 2._rprec*log(Dz/z0) + 4._rprec*z0/Dz*log(Dz/z0))
        iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1._rprec/vonk**2.0*uty**2.0*Dz         &
            * ((log(Dz/z0) - 1._rprec)**2.0 - 2._rprec*z0/Dz + 1._rprec)       &
            + 1._rprec/3._rprec*Ay**2.0*Dz*(1._rprec - z0/Dz)**3.0             &
            - 1._rprec/2._rprec/vonk*uty*Ay*Dz                                 &
            * (1._rprec - 4._rprec*z0/Dz+3._rprec*z0**2.0/Dz**2.0              &
            - 2._rprec*log(Dz/z0) + 4._rprec*z0/Dz*log(Dz/z0))

        ! calculate top and bottom derivatives
        ! Eq. D5 (a) in Yang et al 2015
        iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = 1._rprec/Dz*(Ax+utx/vonk)
        iwm_dudzT(iwm_i,iwm_j,iwm_diry) = 1._rprec/Dz*(Ay+uty/vonk)
        ! Eq. D5 (b) in Yang et al 2015
        iwm_dudzB(iwm_i,iwm_j,iwm_dirx) = 1._rprec/Dz*Ax+utx/vonk/z0
        iwm_dudzB(iwm_i,iwm_j,iwm_diry) = 1._rprec/Dz*Ay+uty/vonk/z0
        ! Eq. D6 in Yang et al 2015
        dVelzT = abs(Ux/Vel*iwm_dudzT(iwm_i,iwm_j,iwm_dirx)                    &
            + Uy/Vel*iwm_dudzT(iwm_i,iwm_j,iwm_diry))
        dvelzB = sqrt(iwm_dudzB(iwm_i,iwm_j,iwm_dirx)**2._rprec                &
            + iwm_dudzB(iwm_i,iwm_j,iwm_diry)**2._rprec)

        ! shear stress at the top
        ! Eq. D4(a) in Yang et al 2015
        iwm_tautx(iwm_i,iwm_j) = (vonk*Dz)**2._rprec*dVelzT                    &
            *iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
        iwm_tauty(iwm_i,iwm_j) = (vonk*Dz)**2._rprec*dVelzT                    &
            *iwm_dudzT(iwm_i,iwm_j,iwm_diry)

        ! calculate the wall stress
        ! Eq. D4(b) in Yang et al 2015
        if (flag1(iwm_i,iwm_j)/=0 .or. flag2(iwm_i,iwm_j)/=0 .or.              &
            flag2(iwm_i,iwm_j)/=0) then
            equilWMpara = Vel*vonk**2._rprec/(log(Dz/z0))**2._rprec
            iwm_tauwx(iwm_i,iwm_j) = equilWMpara*Ux
            iwm_tauwy(iwm_i,iwm_j) = equilWMpara*Uy
        else
            ! Eq. D4
            iwm_tauwx(iwm_i,iwm_j) = (vonk*z0)**2._rprec*dVelzB                &
                *iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
            iwm_tauwy(iwm_i,iwm_j) = (vonk*z0)**2._rprec*dVelzB                &
                *iwm_dudzB(iwm_i,iwm_j,iwm_diry)
        end if

        ! Calculate turbulent diffusion
        ! Difference in shear stresses
        iwm_diff(iwm_i,iwm_j,iwm_dirx) =                                       &
            iwm_tautx(iwm_i,iwm_j) - iwm_tauwx(iwm_i,iwm_j)
        iwm_diff(iwm_i,iwm_j,iwm_diry) =                                       &
            iwm_tauty(iwm_i,iwm_j) - iwm_tauwy(iwm_i,iwm_j)

    end if 

    ! the filtered friction velocity used for filtering time scale
    iwm_flt_us(iwm_i,iwm_j) = iwm_flt_us(iwm_i,iwm_j)                          &
        *(1._rprec-iwm_tR(iwm_i,iwm_j))+utau*iwm_tR(iwm_i,iwm_j)

    ! update the filtering time scale
    ! Eq. 26
    iwm_tR(iwm_i,iwm_j) = iwm_dt/(Dz/iwm_flt_us(iwm_i,iwm_j)/vonk)

    ! filtering time scale can only be larger than the time step,
    ! if not, then just use the instantaneous flow field to do the model
    if (iwm_tR(iwm_i,iwm_j)>1._rprec) then
        iwm_tR(iwm_i,iwm_j) = 1._rprec
    end if

end do
end do

end subroutine iwm_calc_wallstress


!*******************************************************************************
subroutine iwm_monitor
!*******************************************************************************
!
! This subroutine is to monitor the parameters at one point, do not call this
! subroutine if you are not interested in how the model works
!
use param, only : nx,ny,jt_total
use open_file_fid_mod
use sim_param, only : u,v,p
implicit none

integer :: iwm_i,iwm_j,dmpPrd,fid
character*50 :: fname

dmpPrd = iwm_ntime_skip
iwm_i = int(nx/2._rprec)
iwm_j = int(ny/2._rprec)

write(fname,'(A,i5.5,A)') 'iwm_track.dat'
fid = open_file_fid(fname, 'append', 'formatted' )
if (jt_total == dmpPrd) then
write(fid,*) 'u_inst    ','v_inst    ','U    ','V    ','rhsx    ','rhsy    ',  &
    'utx    ','uty    ','u*_filt    ','Lu    ','Lv    ','Luv    ','Luu    ',   &
    'Lvv    ','Ax    ','Ay    ','tR    ','delta_i    ','tauwx    ','tauwy    ',&
    'tautx    ','tauty    ','flag1    ','flag2    ','flag3    ','flag4    ',   &
    'flag4_count    ','flag5    ','Axp    ','Ayp    '
end if
if( mod(jt_total,dmpPrd)==0)then
write(fid,*) u(iwm_i,iwm_j,1), v(iwm_i,iwm_j,1),                         &
    iwm_flt_tagvel(iwm_i,iwm_j,:),rhs(iwm_i,iwm_j,iwm_dirx),                   &
    rhs(iwm_i,iwm_j,iwm_diry),iwm_utx(iwm_i,iwm_j),                            &
    iwm_uty(iwm_i,iwm_j), iwm_flt_us(iwm_i,iwm_j),                             &
    iwm_inte(iwm_i,iwm_j,iwm_Lu),iwm_inte(iwm_i,iwm_j,iwm_Lv),                 &
    iwm_inte(iwm_i,iwm_j,iwm_Luv),iwm_inte(iwm_i,iwm_j,iwm_Luu),               &
    iwm_inte(iwm_i,iwm_j,iwm_Lvv),iwm_Ax(iwm_i,iwm_j), iwm_Ay(iwm_i,iwm_j),    &
    iwm_tR(iwm_i,iwm_j),delta_i_array(iwm_i,iwm_j),iwm_tauwx(iwm_i,iwm_j),     &
    iwm_tauwy(iwm_i,iwm_j),iwm_tautx(iwm_i,iwm_j),iwm_tauty(iwm_i,iwm_j),      &
    flag1(iwm_i,iwm_j),flag2(iwm_i,iwm_j),flag3(iwm_i,iwm_j),flag4(iwm_i,iwm_j),&
    flag4_count,flag5(iwm_i,iwm_j),Axp(iwm_i,iwm_j),Ayp(iwm_i,iwm_j),          &
    rhs_m(iwm_i,iwm_j,:),rhs_tot(iwm_i,iwm_j,:),iwm_z0(iwm_i,iwm_j),           &
    p(iwm_i,iwm_j,1),iwm_flt_p(iwm_i,iwm_j),iwm_PrsGrad(iwm_i,iwm_j,:),        &
    flag6(iwm_i,iwm_j),Luux(iwm_i,iwm_j),Lvvy(iwm_i,iwm_j),Luvx(iwm_i,iwm_j),  &
    Luvy(iwm_i,iwm_j),Lux(iwm_i,iwm_j),Lvy(iwm_i,iwm_j),iwm_conv(iwm_i,iwm_j,:)

end if
close(fid)

end subroutine iwm_monitor


!*******************************************************************************
subroutine iwm_checkPoint()
!*******************************************************************************
!
! This subroutine checkpoints the integral wall model. It is called after making
! sure lbc_mom=3
!
use open_file_fid_mod
implicit none

integer :: fid

fid = open_file_fid('iwm_checkPoint.dat', 'rewind', 'unformatted' )
write(fid) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:),         &
    iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN),              &
    iwm_flt_p(:,:), iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),          &
    iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN),&
    iwm_diff(:,:,1:iwm_DN), rhs(:,:,1:iwm_DN), iwm_dudzT(:,:,1:iwm_DN),        &
    iwm_dudzB(:,:,1:iwm_DN), iwm_flt_us(:,:), iwm_tR(:,:), iwm_Dz(:,:),        &
    iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), iwm_dt,                             &
    rhs_m(:,:,1:iwm_DN), rhs_tot(:,:,1:iwm_DN), iwm_utx(:,:), iwm_uty(:,:)
close(fid)

end subroutine iwm_checkPoint

!*******************************************************************************
subroutine iwm_read_checkPoint()
!*******************************************************************************
!
! This subroutine reads the check point data for the integral wall model. It is
! called after making sure lbc_mom=3
!
use open_file_fid_mod
implicit none

integer :: fid

fid = open_file_fid('iwm_checkPoint.dat', 'rewind', 'unformatted' )
read(fid) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:),         &
    iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN),              &
    iwm_flt_p(:,:), iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),          &
    iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN),&
    iwm_diff(:,:,1:iwm_DN), rhs(:,:,1:iwm_DN), iwm_dudzT(:,:,1:iwm_DN),    &
    iwm_dudzB(:,:,1:iwm_DN), iwm_flt_us(:,:), iwm_tR(:,:), iwm_Dz(:,:),        &
    iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), iwm_dt,                             &
    rhs_m(:,:,1:iwm_DN), rhs_tot(:,:,1:iwm_DN), iwm_utx(:,:), iwm_uty(:,:)
close(fid)

end subroutine iwm_read_checkPoint

end module iwmles
