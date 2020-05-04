!!
!!  Copyright (C) 2019  Johns Hopkins University
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
module coriolis
!*******************************************************************************
! This module contains all of the subroutines associated with scalar transport
use types, only : rprec
use pid_m
implicit none

save
private

public :: coriolis_init, coriolis_finalize, coriolis_calc


! Coriolis forcing
!   0->Off; 1->Fixed geostrophic wind;
!   2-> PID control of angle at specific height;
!   3-> interpolate from input file coriolis.dat
integer, public :: coriolis_forcing = 0

! fc -> coriolis parameter fc (dimensional)
! G -> geostrophic velocity (dimensional)
! alpha -> angle of geostrophic velocity
!   (in radians measured counter-clockwise from x-direction)
real(rprec), public :: fc = 0.0001
real(rprec), public :: G = 8.0
real(rprec), public :: alpha = 0.78539816339

! phi_set -> angle of planar-averaged velocity
!   (in radians measured counter-clockwise from x-direction)
! height_set -> height of angle set point (dimensional)
! Kp, Ki, Kd -> PID controller gains (dimensional)
integer, public :: pid_time = 0
real(rprec), public :: phi_set = 0.0
real(rprec), public :: height_set = 100
real(rprec), public :: Kp = 1e-4
real(rprec), public :: Ki = 3.8e-8
real(rprec), public :: Kd = 0.0

! Repeat_interval
real(rprec), public :: repeat_interval = 86400.0

! Geostrphic wind component
real(rprec) :: ug, vg

! PID controller
type(pid_t) :: pid

! Rotation rate
real(rprec) :: pid_rot_rate = 0._rprec
real(rprec), public :: phi_actual = 0._rprec

! Determining location of planar-averaging
logical :: plane_in_coord = .false.
integer :: k1, k2       ! Weighting locations
real(rprec) :: w1, w2   ! Weight values

! Interpolation of geostrophic wind
real(rprec), dimension(:), allocatable :: t_interp, alpha_interp, G_interp

contains

!*******************************************************************************
subroutine coriolis_init
!*******************************************************************************
! This subroutine initializes the variables for the coriolis
use param, only : z_i, u_star, L_z, nz, coord, nproc, dz, read_endian
use grid_m
use functions, only : binary_search, count_lines
logical :: exst
real(rprec) :: e_int
integer :: num_t, fid, i

! Non-dimensionalize
G = G/u_star
fc = fc*z_i/u_star
Ki = Ki*z_i/u_star
Kd = Kd*u_star/z_i
height_set = height_set/z_i

if (coriolis_forcing == 2) then

    ! Create PID controller
    pid = pid_t(Kp, Ki, Kd, phi_set)

    inquire (file='coriolis_pid.out', exist=exst)
    if (exst) then
        open(12, file='coriolis_pid.out', form='unformatted', convert=read_endian)
        read(12) e_int, alpha
        close(12)
        pid%e_int = e_int
    end if

    ! Create
    if (height_set < 0.5*dz ) then
        if (coord == 0) then
            plane_in_coord = .true.
            k1 = 1
            k2 = 2
            w1 = 1._rprec
            w2 = 0._rprec
        end if
    else if (height_set > L_z) then
        if (coord == nproc-1) then
            plane_in_coord = .true.
            k1 = Nz-2
            k2 = Nz-1
            w1 = 0._rprec
            w2 = 1._rprec
        end if
    else
        if (height_set >= grid%z(1) .and. height_set <= grid%z(nz)) then
            plane_in_coord = .true.
            k1 = binary_search(grid%z(1:), height_set)
            k2 = k1 + 1
            w1 = (height_set - grid%z(k1)) / dz
            w2 = 1._rprec - w1
        end if
    end if

else if (coriolis_forcing == 3) then
    ! Count number of entries and allocate
    num_t = count_lines('coriolis.dat')
    allocate( t_interp(num_t) )
    allocate( alpha_interp(num_t) )
    allocate( G_interp(num_t) )

    ! Read values from file
    open(newunit=fid, file='coriolis.dat', status='unknown', form='formatted', &
        position='rewind')
    do i = 1, num_t
        read(fid,*) t_interp(i), G_interp(i), alpha_interp(i)
    end do

end if

! Set components
ug = G*cos(alpha)
vg = G*sin(alpha)

end subroutine coriolis_init

!*******************************************************************************
subroutine coriolis_finalize
!*******************************************************************************
use param, only : read_endian

if (coriolis_forcing == 2) then
    open(12, file='coriolis_pid.out', form='unformatted', convert=read_endian)
    write(12) pid%e_int, alpha
    close(12)
end if

end subroutine coriolis_finalize

!*******************************************************************************
subroutine coriolis_calc
!*******************************************************************************
use param, only : MPI_RPREC, comm, ierr, dt, total_time_dim, u_star, jt_total
use sim_param, only : u, v, RHSx, RHSy, nx, ny, nz
use functions, only : linear_interp
use mpi
real(rprec) :: ubar = 0, vbar = 0, temp

if (coriolis_forcing == 2) then
    if (plane_in_coord) then
        ubar = w1*sum(u(1:nx,:,k1))/(nx*ny) + w2*sum(u(1:nx,:,k2))/(nx*ny)
        vbar = w1*sum(v(1:nx,:,k1))/(nx*ny) + w2*sum(v(1:nx,:,k2))/(nx*ny)
    else
        ubar = 0._rprec
        vbar = 0._rprec
    end if  

#ifdef PPMPI
    call MPI_AllReduce(ubar, temp, 1, MPI_RPREC, MPI_SUM, comm, ierr)
    ubar = temp
    call MPI_AllReduce(vbar, temp, 1, MPI_RPREC, MPI_SUM, comm, ierr)
    vbar = temp
#endif

    phi_actual = atan2(vbar,ubar)
    if (jt_total < pid_time) then
        ! Use PID to get new angle
        pid_rot_rate = -pid%advance(phi_actual, dt)
        alpha = alpha - pid_rot_rate*dt

        ! Set components
        ug = G*cos(alpha)
        vg = G*sin(alpha)
    else 
        pid_rot_rate = 0._rprec 
    end if

else if( coriolis_forcing == 3) then
    ! interpolate
    alpha = linear_interp(t_interp, alpha_interp,                              &
        mod(total_time_dim, repeat_interval))
    G = linear_interp(t_interp, G_interp,                                      &
        mod(total_time_dim, repeat_interval))/u_star

    ! Set components
    ug = G*cos(alpha)
    vg = G*sin(alpha)

    write(*,*) total_time_dim, G, alpha, ug, vg
end if

! Coriolis: add forcing to RHS
if (coriolis_forcing > 0) then
    ! This is to put in the coriolis forcing using coriol,ug and vg as
    ! precribed in param.f90. (ug,vg) specfies the geostrophic wind vector
    ! Note that ug and vg are non-dimensional (using u_star in param.f90)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fc * v(:,:,1:nz-1) - fc * vg
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) - fc * u(:,:,1:nz-1) + fc * ug
    if (coriolis_forcing == 2) then
        RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + pid_rot_rate * v(:,:,1:nz-1)
        RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) - pid_rot_rate * u(:,:,1:nz-1)
    end if

end if

end subroutine coriolis_calc



end module coriolis
