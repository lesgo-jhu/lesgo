!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by: 
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   The Actuator Turbine Model is free software: you can redistribute it 
!!   and/or modify it under the terms of the GNU General Public License as 
!!   published by the Free Software Foundation, either version 3 of the 
!!   License, or (at your option) any later version.
!!
!!   The Actuator Turbine Model is distributed in the hope that it will be 
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************************
module atm_lesgo_interface
!*******************************************************************************
! This module interfaces actuator turbine module with lesgo
! It is a lesgo specific module, unlike the atm module

! Remember to always dimensinoalize the variables from LESGO
! Length is non-dimensionalized by z_i

! Lesgo data used regarding the grid (LESGO)
use param, only : dt ,nx,ny,nz,jzmin,jzmax,dx,dy,dz,coord,lbz,nproc, z_i, u_star
! nx, ny, nz - nodes in every direction
! z_i - non-dimensionalizing length
! dt - time-step 

! These are the forces on x,y, and z respectively
use sim_param, only : fxa, fya, fza

! Grid definition (LESGO)
use grid_defs, only : grid

! Interpolating function for interpolating the velocity field to each
! actuator point
use functions, only : trilinear_interp

! Actuator Turbine Model module
use atm_base
use actuator_turbine_model
use atm_input_util, only : rprec, turbineArray, turbineModel

implicit none

private
public atm_lesgo_initialize, atm_lesgo_forcing

! This is a dynamic list that stores all the points in the domain with a body 
! force due to the turbines. 
type bodyForce_t
!     i,j,k stores the index for the point in the domain
    integer i,j,k
    real(rprec), dimension(3) :: force ! Force vector
    real(rprec), dimension(3) :: location ! Position vector
end type bodyForce_t

! This is the body force vector which has the forces by the atm
!type(bodyForce_t), allocatable, dimension(:) :: bodyForce

! These are the pointers to the grid arrays
real(rprec), pointer, dimension(:) :: x,y,z

! Vector used to store x, y, z locations
real(rprec), dimension(3) :: vector_point

! Integers to loop through x, y, and z  
integer :: i,j,k

! Integer to loop through all turbines
integer :: m

! Counter to establish number of points which are influenced by body forces
integer :: c0

! Body force field
type(bodyForce_t), allocatable, dimension(:) :: forceField

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_initialize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model

! This is the list of body forces for all the cells 
!type(DynamicList_), pointer :: bodyForceField
real(rprec) :: zeroVector(3)     ! A vector of zeros used for initiallizing vec
integer :: i, j, k

! Define all elements of zero vector to zero
do i=1,3
    zeroVector(i)=0.
enddo


! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

call atm_initialize () ! Initialize the atm (ATM)

! This will find all the locations that are influenced by each turbine
! It depends on a sphere centered on the rotor that extends beyond the blades
c0=0  ! Initialize conuter
do k=1,nz ! Loop through grid points in z
    do j=1,ny ! Loop through grid points in y
        do i=1,nx ! Loop through grid points in x
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i
            vector_point(3)=z(k)*z_i
            do m=1,numberOfTurbines
                if (distance(vector_point,turbineArray(m) % baseLocation) &
                    .le. turbineArray(m) % sphereRadius ) then
                    c0=c0+1
!                    call appendVector(bodyForceField,zeroVector)
                end if
            end do 
        end do
    end do
end do

allocate(forceField(c0))
write(*,*) 'Numebr of cells being affected by ATM = ', c0

c0=0
! Run the same loop and save all variables
do k=1,nz ! Loop through grid points in z
    do j=1,ny ! Loop through grid points in y
        do i=1,nx ! Loop through grid points in x
            do m=1,numberOfTurbines
                if (distance(vector_point,turbineArray(m) % baseLocation) &
                    .le. turbineArray(m) % sphereRadius ) then
                    c0=c0+1
                    forceField(c0) % i = i
                    forceField(c0) % j = j
                    forceField(c0) % k = k
                    forceField(c0) % location = vector_point
                end if
            end do 
        end do
    end do
end do



end subroutine atm_lesgo_initialize


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_forcing ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutines calls the update function from the ATM Library
! and calculates the body forces needed in the domain

! Reset applied force arrays to zero
!fxa = 0._rprec
!fya = 0._rprec
!fza = 0._rprec

call atm_update(dt)

do i=1,numberOfTurbines
    call atm_lesgo_force(i)
enddo

end subroutine atm_lesgo_forcing

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will feed the velocity at all the actuator points into the atm
! This is done by using trilinear interpolation from lesgo
! Force will be calculated based on the velocities

! Use the velocity field (LESGO)
use sim_param,only:u,v,w ! Load the velocity components

integer, intent(in) :: i ! The turbine number
integer :: m,n,q,j
real(rprec), dimension(3) :: velocity
real(rprec), dimension(3) :: xyz    ! Point onto which to interpolate velocity
!real(rprec) :: u_i, v_i, w_i ! Interpolated velocities
j=turbineArray(i) % turbineTypeID ! The turbine type ID

do q=1, turbineArray(i) %  numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
                       
            ! Actuator point onto which to interpolate the velocity
            xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

            ! Non-dimensionalizes the point location 
            xyz=vector_divide(xyz,z_i)
            ! Interpolate velocities
!            u_i=trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)  ! Vel in x
!            v_i=trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)  ! Vel in y
!            w_i=trilinear_interp(w(1:nx,1:ny,lbz:nz),lbz,xyz)  ! Vel in z
            velocity(1)=trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)!*u_star
            velocity(2)=trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)!*u_star
            velocity(3)=trilinear_interp(w(1:nx,1:ny,lbz:nz),lbz,xyz)!*u_star
write(*,*) 'Error NOT Here'

            call atm_computeBladeForce(i,m,n,q,velocity)
        enddo
    enddo
enddo

end subroutine atm_lesgo_force


end module atm_lesgo_interface
