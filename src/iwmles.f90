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

! u_tau,x  u_tau,y
real(rprec), dimension(:,:), allocatable :: iwm_utx, iwm_uty
! tau_wall,x tau_wall,y
real(rprec), dimension(:,:), allocatable :: iwm_tauwx, iwm_tauwy
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
! turbulent diffusion, LHS
real(rprec), dimension(:,:,:), allocatable :: iwm_diff, iwm_LHS
! dudz at z=dz/2, dudz at z=zo
real(rprec), dimension(:,:,:), allocatable :: iwm_dudzT, iwm_dudzB
! filtered friction velocity, filtering time scale
real(rprec), dimension(:,:), allocatable :: iwm_flt_us, iwm_tR
! HALF cell height, zo, linear correction in x, y directions
real(rprec), dimension(:,:), allocatable :: iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay

! number of time steps to skip between wall stress calculations
integer :: iwm_ntime_skip = 5
! time step size seen by the wall model
real(rprec) :: iwm_dt

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
    ! gather flow status, update the integrated unsteady term, convective term,
    ! turbulent diffusion term etc.
    call iwm_calc_lhs()
    ! the subroutine to calculate wall stress
    call iwm_calc_wallstress()
    ! this is to monitor any quantity from the iwm, useful debugging tool
    call iwm_monitor()
end if

! Imposing txz, tyz, dudz, dvdz every time step even iwm_* are not computed
! every time step.
do iwm_i = 1, nx
do iwm_j = 1, ny
    ! wall stress, use the value calculated in iwm, note the negative sign
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
use param, only : nx, ny, dz, vonk, zo, cfl, L_x

implicit none

real(rprec) :: usinit, uinit, vinit, Dzp

! initial value for us (the friction velocity)
usinit= 1._rprec
! initial value for the x-velocity at first grid point
uinit = usinit/vonk*log(dz/2._rprec/zo)
! initial value for the y-velocity at first grid point
vinit = 0._rprec
! at the height of the first grid point.
Dzp=dz/2._rprec

! us in x, y directions
allocate(iwm_utx(nx,ny))
allocate(iwm_uty(nx,ny))
iwm_utx = usinit
iwm_uty = 0._rprec

! wall stress in x, y directions
allocate(iwm_tauwx(nx,ny))
allocate(iwm_tauwy(nx,ny))
iwm_tauwx = usinit**2._rprec
iwm_tauwy = 0._rprec

! filitered velocity at the first grid point in x, y directions
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
allocate(iwm_LHS    (nx,ny,iWM_DN))
allocate(iwm_dudzT  (nx,ny,iwm_DN))
allocate(iwm_dudzB  (nx,ny,iwm_DN))
iWM_unsdy   = 0._rprec
iWM_conv    = 0._rprec
iWM_PrsGrad = 0._rprec
iwm_diff    = 0._rprec
iWM_LHS     = -uinit*Dzp
iwm_dudzT(:,:,iwm_dirx) = usinit/vonk/Dzp
iwm_dudzT(:,:,iwm_diry) = 0._rprec
iwm_dudzB(:,:,iwm_dirx) = usinit/vonk/zo
iwm_dudzB(:,:,iwm_diry) = 0._rprec

! filtered friction velocity and the filtering time scale, tR<1
allocate(iwm_flt_us(nx,ny))
allocate(iwm_tR    (nx,ny))
iwm_flt_us = usinit
iwm_tR = (cfl*L_x/nx/uinit)/(dz/2._rprec/vonk/usinit)

! cell height and imposed roughness length
allocate(iwm_Dz(nx,ny))
allocate(iwm_z0(nx,ny))
iWM_Dz = dz/2._rprec
iWM_z0 = zo !we leave the possibility of zo as a function of x-y

! linear correction to the log profile
allocate(iwm_Ax(nx,ny))
allocate(iwm_Ay(nx,ny))
iwm_Ax = 0._rprec
iwm_Ay = 0._rprec

! time step seen by the iwm
iwm_dt=iwm_ntime_skip*cfl*L_x/nx/uinit

end subroutine iwm_init

!*******************************************************************************
subroutine iwm_finalize
!*******************************************************************************
!
! This subroutine deallocates memory used for iwm
!
implicit none

deallocate(iwm_utx)
deallocate(iwm_uty)

deallocate(iwm_tauwx)
deallocate(iwm_tauwy)

deallocate(iwm_flt_tagvel  )
deallocate(iwm_flt_tagvel_m)

deallocate(iwm_inte  )
deallocate(iwm_inte_m)

deallocate(iwm_unsdy  )
deallocate(iwm_conv   )
deallocate(iwm_PrsGrad)
deallocate(iwm_diff   )
deallocate(iwm_LHS    )
deallocate(iwm_dudzT  )
deallocate(iwm_dudzB  )

deallocate(iwm_flt_us)
deallocate(iwm_tR    )

deallocate(iwm_Dz)
deallocate(iwm_z0)
deallocate(iwm_Ax)
deallocate(iwm_Ay)

end subroutine iwm_finalize


!*******************************************************************************
subroutine iwm_calc_lhs()
!*******************************************************************************
!
! Ths subroutine calculates the left hand side of the iwm system.
!
use grid_m, only : grid
use types,only : rprec
use param,only : nx,ny,dx,dy,ld
use sim_param,only : u,v,w,p
use test_filtermodule
implicit none

! useful array for autowraped index
integer, pointer, dimension(:) :: autowrap_i, autowrap_j
integer :: iwm_i,iwm_j
! the instantaneous field
real(rprec), dimension(ld,ny) :: u_inst, v_inst, w_inst, p_inst
! the mean pressure at first grid point
real(rprec) :: p_bar
! for temporary storage of derivativex of integrals like dLudx, dLvdx...
real(rprec) :: Luux, Luvx, Luvy, Lvvy, Lux, Lvy
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

! calculate LHS, calculation of the integrals is done from the last time step
! in the subroutine iwm_calc_wallstress, so is iwm_diff
do iwm_i = 1, nx
do iwm_j = 1, ny
    ! the unsteady term
    iwm_unsdy(iwm_i,iwm_j,iwm_dirx) =                                          &
        (iwm_inte(iwm_i,iwm_j,iwm_Lu)-iwm_inte_m(iwm_i,iwm_j,iwm_Lu))/iwm_dt
    iwm_unsdy(iwm_i,iwm_j,iwm_diry) =                                          &
        (iwm_inte(iwm_i,iwm_j,iwm_Lv)-iwm_inte_m(iwm_i,iwm_j,iwm_Lv))/iwm_dt

    ! the convective term
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luu)
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luu)
    Luux = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Luv)
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Luv)
    Luvy = (phip-phim)/dy/2._rprec
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luv)
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luv)
    Luvx = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lvv)
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lvv)
    Lvvy = (phip-phim)/dy/2._rprec
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Lu )
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Lu )
    Lux = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lv )
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lv )
    Lvy = (phip-phim)/dy/2._rprec
    iwm_conv(iwm_i,iwm_j,iwm_dirx) = Luux + Luvy                               &
        - iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_dirx)*(Lux+Lvy)
    iwm_conv(iwm_i,iwm_j,iwm_diry) = Luvx + Lvvy                               &
        - iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_diry)*(Lux+Lvy)

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

    ! the left hand side
    ! this is the integrated momentum equation, except for the Lu term
    iwm_lhs(iwm_i,iwm_j,iwm_dirx) = -iwm_inte(iwm_i,iwm_j,iwm_Lu)              &
        + iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_dirx)                              &
        + iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)                                    &
        - iwm_diff(iwm_i,iwm_j,iwm_dirx) )
    ! this is the integrated momentum equation, except for the Lv term
    iwm_lhs(iwm_i,iwm_j,iwm_diry) = -iwm_inte(iwm_i,iwm_j,iwm_Lv)              &
        + iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_diry)                              &
        + iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)                                    &
        - iwm_diff(iwm_i,iwm_j,iwm_diry) )
end do
end do

nullify(autowrap_i, autowrap_j)

end subroutine iwm_calc_lhs

!*******************************************************************************
subroutine iwm_slv(lhsx,lhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy)
!*******************************************************************************
use types, only : rprec
use param, only : vonk
implicit none

real(rprec), intent(in)  :: lhsx, lhsy, Ux, Uy, Dz, z0, utx, uty
real(rprec), intent(out) :: fx,fy
real(rprec) :: Ax, Ay, Vel, inteLu, inteLv

Ax = (Ux-utx/vonk*log(Dz/z0))/((1.0-z0/Dz))
Ay = (Uy-uty/vonk*log(Dz/z0))/((1.0-z0/Dz))
Vel = sqrt(Ux**2.0+Uy**2.0)
inteLu = 1.0/2.0*Dz*Ax*(1.0-z0/Dz)**2.0+1.0/vonk*utx*Dz*(z0/Dz-1.0+log(Dz/z0))
inteLv = 1.0/2.0*Dz*Ay*(1.0-z0/Dz)**2.0+1.0/vonk*uty*Dz*(z0/Dz-1.0+log(Dz/z0))
fx = inteLu+lhsx
fy = inteLv+lhsy

end subroutine iwm_slv

!*******************************************************************************
subroutine iwm_calc_wallstress
!*******************************************************************************
use types, only : rprec
use param, only : vonk, nx, ny
use test_filtermodule

implicit none

integer :: iwm_i, iwm_j
real(rprec) :: fx, fy, fxp, fyp
real(rprec) :: iwm_tol, iwm_eps
real(rprec) :: a11, a12, a21, a22
real(rprec) :: iwmutxP,iwmutyP
integer :: iter, MaxIter, equil_flag, div_flag
real(rprec) :: equilWMpara,equilutx,equiluty
real(rprec) :: iwmpAx, iwmpAy, iwmputx,iwmputy,iwmpz0,iwmpDz
real(rprec) :: utaup
real(rprec) :: dVelzT, dVelzB, Vel

MaxIter=1500

iwm_tol = 0.000001_rprec
iwm_eps = 0.000000001_rprec

do iwm_i=1,nx
do iwm_j=1,ny

    ! use Newton method to solve the system
    iwm_utx(iwm_i,iwm_j) = 1._rprec                                            &
        * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx))
    iwm_uty(iwm_i,iwm_j) = 0.1_rprec                                           &
        * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))

    call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx), iwm_lhs(iwm_i,iwm_j,iwm_diry), &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),                                  &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),                                  &
        iwm_Dz(iwm_i,iwm_j), iwm_z0(iwm_i,iwm_j),                              &
        iwm_utx(iwm_i,iwm_j), iwm_uty(iwm_i,iwm_j),fx,fy )

    iter = 0
    equil_flag = 0
    div_flag = 0
    do while (max(abs(fx),abs(fy))>iwm_tol)
        iwmutxP=iwm_utx(iwm_i,iwm_j)+iWM_eps
        iwmutyP=iwm_uty(iwm_i,iwm_j)
        call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),                            &
            iwm_lhs(iwm_i,iwm_j,iwm_diry),iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),&
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),iwm_Dz(iwm_i,iwm_j),          &
            iwm_z0(iwm_i,iwm_j), iwmutxP, iwmutyP, fxp, fyp )
        a11 = (fxp-fx)/iwm_eps
        a21 = (fyp-fy)/iwm_eps
        iwmutxP = iwm_utx(iwm_i,iwm_j)
        iwmutyP = iwm_uty(iwm_i,iwm_j)+iwm_eps
        call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),                            &
            iwm_lhs(iwm_i,iwm_j,iwm_diry),                                     &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),                              &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),                              &
            iwm_Dz(iwm_i,iwm_j), iwm_z0(iwm_i,iwm_j),                          &
            iwmutxP, iwmutyP, fxp, fyp)
        a12 = (fxp-fx)/iwm_eps
        a22 = (fyp-fy)/iwm_eps
        iwm_utx(iwm_i,iwm_j) = iwm_utx(iwm_i,iwm_j)                            &
            - 0.50*( a22*fx-a12*fy)/(a11*a22-a12*a21)
        iwm_uty(iwm_i,iwm_j) = iwm_uty(iwm_i,iwm_j)                            &
            - 0.50*(-a21*fx+a11*fy)/(a11*a22-a12*a21)
        call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),                            &
            iwm_lhs(iwm_i,iwm_j,iwm_diry),                                     &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),                              &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),                              &
            iwm_Dz(iwm_i,iwm_j), iwm_z0(iwm_i,iwm_j),                          &
            iwm_utx(iwm_i,iwm_j), iwm_uty(iwm_i,iwm_j), fx, fy)
        iter = iter+1

        ! maximum iteration reached
        if (iter>MaxIter) then
            equil_flag = 1
            div_flag = 1;
            exit
        end if
    end do

    ! infinity check
    if (iwm_utx(iwm_i,iwm_j)-1.0==iwm_utx(iwm_i,iwm_j) .or.                    &
        iwm_uty(iwm_i,iwm_j)-1.0==iwm_uty(iwm_i,iwm_j)) then
        equil_flag=1
        div_flag  =1
    end if

    ! calculate equilibrium us for equil_flag=1 use
    equilutx = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                       &
    / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
    equiluty = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                       &
        / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
    if (equil_flag==1) then
        iwm_utx(iwm_i,iwm_j) = equilutx
        iwm_uty(iwm_i,iwm_j) = equiluty
    end if

    !calculate Ax, Ay
    if(equil_flag==1)then
        iwmpAx = 0._rprec
        iwmpAy=0._rprec
    else
        ! eq. D2 in Yang et al. 2015
        iwmpAx = ( iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                        &
            -iwm_utx(iwm_i,iwm_j)/vonk*log(iwm_Dz(iwm_i,iwm_j)                 &
            /iwm_z0(iwm_i,iwm_j)))                                             &
            / ((1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j)))
        iwmpAy = ( iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                        &
            -iwm_uty(iwm_i,iwm_j)/vonk*log(iwm_Dz(iwm_i,iwm_j)                 &
            /iwm_z0(iwm_i,iwm_j)))                                             &
            /((1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j)))
    end if

    ! check for excessive linear term correction
    ! after first 100 step this check is rarely invoked
    if (abs(iwmpAx)>1._rprec .or. abs(iwmpAy)>1._rprec) then
        equil_flag = 1
        iwm_utx(iwm_i,iwm_j) = equilutx
        iwm_uty(iwm_i,iwm_j) = equiluty
        iwmpAx = 0._rprec
        iwmpAy = 0._rprec
    end if

    ! store the linear correction
    iwm_Ax(iwm_i,iwm_j) = iwmpAx
    iwm_Ay(iwm_i,iwm_j) = iwmpAy

    ! update integral for last time step
    iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu )
    iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv )
    iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv)
    iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu)
    iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv)

    !those temporary variables are used for convenient reference
    iwmputx = iwm_utx(iwm_i,iwm_j)
    iwmputy = iwm_uty(iwm_i,iwm_j)
    iwmpDz  = iwm_Dz (iwm_i,iwm_j)
    iwmpz0  = iwm_z0 (iwm_i,iwm_j)

    ! calculate the needed integrals

    ! Eq. D7 in Yang et al. 2015
    iwm_inte(iwm_i,iwm_j,iwm_Lu) = 1.0/2.0*iwmpDz*iwmpAx                       &
        *(1.0-iwmpz0/iwmpDz)**2.0                                              &
        +1.0/vonk*iwmputx*iwmpDz*(iwmpz0/iwmpDz-1.0+log(iwmpDz/iwmpz0))
    iwm_inte(iwm_i,iwm_j,iwm_Lv) = 1.0/2.0*iwmpDz*iwmpAy                       &
        *(1.0-iwmpz0/iwmpDz)**2.0                                              &
        +1.0/vonk*iwmputy*iwmpDz*(iwmpz0/iwmpDz-1.0+log(iwmpDz/iwmpz0))

    ! Eq. D8 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1.0/vonk**2.0*iwmputx*iwmputy*iwmpDz       &
        *(1.0-2*iwmpz0/iwmpDz+(1.0-log(iwmpDz/iwmpz0))**2.0)                   &
        +1.0/3.0*iwmpAx*iwmpAy*iwmpDz*(1.0 - iwmpz0/iwmpDz)**3.0               &
        -1.0/4.0/vonk*(iwmpAx*iwmputy+iwmpAy*iwmputx)*iwmpDz                   &
        *(1.0-4.0*iwmpz0/iwmpDz+3.0*iwmpz0**2.0/iwmpDz**2.0                    &
        -2.0*log(iwmpDz/iwmpz0)+4.0*iwmpz0/iwmpDz*log(iwmpDz/iwmpz0))
    iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1.0/vonk**2.0*iwmputx**2.0*iwmpDz          &
        *((log(iwmpDz/iwmpz0)-1.0)**2.0-2.0*iwmpz0/iwmpDz+1.0)                 &
        +1.0/3.0*iwmpAx**2.0*iwmpDz*(1.0-iwmpz0/iwmpDz)**3.0                   &
        -1.0/2.0/vonk*iwmputx*iwmpAx*iwmpDz                                    &
        *(1.0-4.0*iwmpz0/iwmpDz+3.0*iwmpz0**2.0/iwmpDz**2.0                    &
        -2.0*log(iwmpDz/iwmpz0)+4.0*iwmpz0/iwmpDz*log(iwmpDz/iwmpz0))

    ! Eq. D9 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1.0/vonk**2.0*iwmputy**2.0*iwmpDz*         &
        ((log(iwmpDz/iwmpz0)-1.0)**2.0-2.0*iwmpz0/iwmpDz+1.0)                  &
        +1.0/3.0*iwmpAy**2.0*iwmpDz*(1.0-iwmpz0/iwmpDz)**3.0                   &
        -1.0/2.0/vonk*iwmputy*iwmpAy*iwmpDz                                    &
        *(1.0-4.0*iwmpz0/iwmpDz-3.0*iwmpz0**2.0/iwmpDz**2.0                    &
        -2.0*log(iwmpDz/iwmpz0)+4.0*iwmpz0/iwmpDz*log(iwmpDz/iwmpz0))

    ! calculate top and bottom derivatives
    ! Eq. D5 (a)
    iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*(iwmpAx+iwmputx/vonk)
    iwm_dudzT(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*(iwmpAy+iwmputy/vonk)
    ! Eq. D5 (b)
    iwm_dudzB(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*iwmpAx+iwmputx/vonk/iwmpz0
    iwm_dudzB(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*iwmpAy+iwmputy/vonk/iwmpz0

    ! calculte the turbulent diffusion term
    !total velocity
    Vel = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)**2._rprec                  &
        +iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)**2._rprec)
    ! Eq. D6
    dVelzT=abs(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel                        &
        *iwm_dudzT(iwm_i,iwm_j,iwm_dirx)+iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)  &
        /Vel*iwm_dudzT(iwm_i,iwm_j,iwm_diry))
    dvelzB = sqrt(iwm_dudzB(iwm_i,iwm_j,iwm_dirx)**2._rprec                    &
        +iwm_dudzB(iwm_i,iwm_j,iwm_diry)**2._rprec)

    ! Eq. D4, the eddy viscosity is nu_T=(vonk*y)^2*dudy, hence the formula
    iwm_diff(iwm_i,iwm_j,iwm_dirx) = (vonk*iwmpDz)**2._rprec*dVelzT            &
        *iwm_dudzT(iwm_i,iwm_j,iwm_dirx)-(vonk*iwmpz0)**2._rprec               &
        *dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
    iwm_diff(iwm_i,iwm_j,iwm_diry) = (vonk*iwmpDz)**2._rprec*dVelzT            &
        *iwm_dudzT(iwm_i,iwm_j,iwm_diry)-(vonk*iwmpz0)**2._rprec               &
        *dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_diry)

    ! calculate the wall stress
    if (equil_flag==1) then
        equilWMpara = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)**2._rprec      &
            + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)**2._rprec)                  &
            *vonk**2._rprec/(log(iwm_Dz(iwm_i,iwm_j)                           &
            /iwm_z0(iwm_i,iwm_j)))**2._rprec
        iwm_tauwx(iwm_i,iwm_j) = equilWMpara                                   &
            *iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
        iwm_tauwy(iwm_i,iwm_j) = equilWMpara                                   &
            *iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
    else
        ! Eq. D4
        iwm_tauwx(iwm_i,iwm_j) = (vonk*iwmpz0)**2._rprec*dVelzB                &
            *iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
        iwm_tauwy(iwm_i,iwm_j) = (vonk*iwmpz0)**2._rprec*dVelzB                &
            *iwm_dudzB(iwm_i,iwm_j,iwm_diry)
    end if

    ! calculate the friciton velocity
    ! definition of friction velocity
    utaup = (iwm_tauwx(iwm_i,iwm_j)**2._rprec                                  &
        +iwm_tauwy(iwm_i,iwm_j)**2._rprec )**0.25_rprec
    ! the filtered friction velocity used for filtering time scale
    iwm_flt_us(iwm_i,iwm_i) = iwm_flt_us(iwm_i,iwm_j)                          &
        *(1._rprec-iwm_tR(iwm_i,iwm_j))+utaup*iwm_tR(iwm_i,iwm_j)

    ! update the filtering time scale
    ! Eq. 26
    iwm_tR(iwm_i,iwm_j) = iwm_dt/(iwm_Dz(iwm_i,iwm_j)                          &
        /iwm_flt_us(iwm_i,iwm_j)/vonk)
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
implicit none

integer :: iwm_i,iwm_j,dmpPrd,fid
character*50 :: fname

dmpPrd = iwm_ntime_skip
iwm_i = int(nx/2._rprec)
iwm_j = int(ny/2._rprec)
write(fname,'(A,i5.5,A)') 'iwm_track.dat'
open(newunit=fid, file=fname, status='unknown', form='formatted',              &
    position='append')
if( mod(jt_total,dmpPrd)==0)then
write(fid,*) iwm_flt_tagvel(iwm_i,iwm_j,:), iwm_utx(iwm_i,iwm_j),              &
    iwm_uty(iwm_i,iwm_j),  iwm_Ax(iwm_i,iwm_j),                                &
    iwm_Ay(iwm_i,iwm_j), iwm_tR(iwm_i,iwm_j)
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
implicit none

integer :: fid

open(newunit=fid, file='iwm_checkPoint.dat', status='unknown',                 &
    form='unformatted', position='rewind')
write(fid) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:),         &
    iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN),              &
    iwm_flt_p(:,:), iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),          &
    iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN),&
    iwm_diff(:,:,1:iwm_DN), iwm_LHS(:,:,1:iwm_DN), iwm_dudzT(:,:,1:iwm_DN),    &
    iwm_dudzB(:,:,1:iwm_DN), iwm_flt_us(:,:), iwm_tR(:,:), iwm_Dz(:,:),        &
    iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), iwm_dt
close(fid)

end subroutine iwm_checkPoint

!*******************************************************************************
subroutine iwm_read_checkPoint()
!*******************************************************************************
!
! This subroutine reads the check point data for the integral wall model. It is
! called after making sure lbc_mom=3
!
implicit none

integer :: fid

open(newunit=fid, file='iwm_checkPoint.dat', status='unknown',                 &
    form='unformatted', position='rewind')
read(fid) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:),          &
    iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN),              &
    iwm_flt_p(:,:), iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),          &
    iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN),&
    iwm_diff(:,:,1:iwm_DN), iwm_LHS(:,:,1:iwm_DN), iwm_dudzT(:,:,1:iwm_DN),    &
    iwm_dudzB(:,:,1:iwm_DN), iwm_flt_us(:,:), iwm_tR(:,:), iwm_Dz(:,:),        &
    iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), iwm_dt
close(fid)

end subroutine iwm_read_checkPoint

end module iwmles