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
use param, only : dt ,nx,ny,nz,dx,dy,dz,coord,nproc, z_i, u_star, lbz, jt_total
! nx, ny, nz - nodes in every direction
! z_i - non-dimensionalizing length
! dt - time-step 

! These are the forces, and velocities on x,y, and z respectively
use sim_param, only : fxa, fya, fza, u, v, w

! Grid definition (LESGO)
use grid_defs, only : grid

! MPI implementation from LESGO
$if ($MPI)
  use mpi_defs
  use mpi
  use param, only : ierr, mpi_rprec, comm, coord
$endif

! Interpolating function for interpolating the velocity field to each
! actuator point
use functions, only : trilinear_interp, interp_to_uv_grid

! Actuator Turbine Model module
use atm_base
use actuator_turbine_model
use atm_input_util, only : rprec, turbineArray, turbineModel, eat_whitespace

implicit none

! Variable for interpolating the velocity in w onto the uv grid
real(rprec), allocatable, dimension(:,:,:) :: w_uv 

private
public atm_lesgo_initialize, atm_lesgo_forcing

! This is a dynamic list that stores all the points in the domain with a body 
! force due to the turbines. 
type bodyForce_t
!     i,j,k stores the index for the point in the domain
    integer i,j,k
    real(rprec), dimension(3) :: force ! Force vector on uv grid
    real(rprec), dimension(3) :: location ! Position vector on uv grid
end type bodyForce_t

! Body force field
type(bodyForce_t), allocatable, dimension(:) :: forceFieldUV, forceFieldW

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_initialize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

! Counter to establish number of points which are influenced by body forces
integer :: cUV, cW  ! Counters for number of points affected on UV and W grids
integer :: i, j, k, m
! Vector used to store x, y, z locations
real(rprec), dimension(3) :: vector_point
! These are the pointers to the grid arrays
real(rprec), pointer, dimension(:) :: x,y,z,zw

! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! Allocate space for the w_uv variable
allocate(w_uv(nx,ny,lbz:nz))

call atm_initialize () ! Initialize the atm (ATM)

! This will find all the locations that are influenced by each turbine
! It depends on a sphere centered on the rotor that extends beyond the blades
cUV=0  ! Initialize conuter
cW=0  ! Initialize conuter
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz-1 ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i

            ! Take into account the UV grid
            vector_point(3)=z(k)*z_i
            do m=1,numberOfTurbines
                if (distance(vector_point,turbineArray(m) % towerShaftIntersect) &
                    .le. turbineArray(m) % sphereRadius ) then
                    cUV=cUV+1
                end if
            enddo 

            ! Take into account the W grid
            vector_point(3)=zw(k)*z_i
            do m=1,numberOfTurbines
                if (distance(vector_point,turbineArray(m) % towerShaftIntersect) &
                    .le. turbineArray(m) % sphereRadius ) then
                    cW=cW+1
                end if
            enddo 

        enddo
    enddo
enddo

! Allocate space for the force fields in UV and W grids
allocate(forceFieldUV(cUV))
allocate(forceFieldW(cW))

write(*,*) 'Number of cells being affected by ATM cUV, cW = ', cUV, cW

cUV=0
cW=0
! Run the same loop and save all variables
! The forceField arrays include all the forces which affect the domain
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz-1 ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i
            vector_point(3)=z(k)*z_i
            do m=1,numberOfTurbines
                if (distance(vector_point,turbineArray(m) % towerShaftIntersect) &
                    .le. turbineArray(m) % sphereRadius ) then
                    cUV=cUV+1
                    forceFieldUV(cUV) % i = i
                    forceFieldUV(cUV) % j = j
                    forceFieldUV(cUV) % k = k
                    forceFieldUV(cUV) % location = vector_point
                endif
            enddo 
            vector_point(3)=zw(k)*z_i
            do m=1,numberOfTurbines
                if (distance(vector_point,turbineArray(m) % towerShaftIntersect) &
                    .le. turbineArray(m) % sphereRadius ) then
                    cW=cW+1
                    forceFieldW(cW) % i = i
                    forceFieldW(cW) % j = j
                    forceFieldW(cW) % k = k
                    forceFieldW(cW) % location = vector_point
                endif
            enddo 
        enddo
    enddo
enddo

! Write the points which are affected
!open(unit=1787,file='./output/forceUV')
!do k=1, size(forceFieldUV)
!            write(1787,*) forceFieldUV(k) % location       
!enddo
!close(1787)
!open(unit=1787,file='./output/forceW')
!do k=1, size(forceFieldW)
!            write(1787,*) forceFieldW(k) % location       
!enddo
!close(1787)

end subroutine atm_lesgo_initialize


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_forcing ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutines calls the update function from the ATM Library
! and calculates the body forces needed in the domain
integer :: i, c

! Establish all turbine forces as zero
do i=1,numberOfTurbines
    turbineArray(i) % bladeForces = 0._rprec
enddo

! Get the velocity from w onto the uv grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

! If statement is for running code only if grid points affected are in this 
! processor. If not, no code is executed at all.
if (size(forceFieldUV) .gt. 0 .or. size(forceFieldW) .gt. 0) then

!    write(*,*) 'dt = ', dt



    ! Update the blade positions based on the time-step
    ! Time needs to be dimensionalized
    call atm_update(dt*z_i/u_star)

    ! Initialize force field to zero
    do c=1,size(forceFieldUV)
        forceFieldUV(c) % force = 0._rprec     
    enddo
    do c=1,size(forceFieldW)
        forceFieldW(c) % force = 0._rprec
    enddo
endif

    ! Calculate forces for all turbines
    do i=1,numberOfTurbines
        call atm_lesgo_force(i)
    enddo

    ! This will gather all the blade forces from all processors
    $if ($MPI)
    do i=1,numberOfTurbines
        call mpi_allreduce(turbineArray(i) % bladeForcesDummy,         &
                           turbineArray(i) % bladeForces,             &
                           size(turbineArray(i) % bladeForces),       &
                           mpi_rprec, mpi_sum, comm, ierr) 
    enddo
    $endif

if (size(forceFieldUV) .gt. 0 .or. size(forceFieldW) .gt. 0) then

    ! Convolute force onto the domain
    do i=1,numberOfTurbines
        call atm_lesgo_convolute_force(i)
    enddo

! This will apply body forces onto the flow field if there are forces within
! this domain
    call atm_lesgo_apply_force()
endif

!write(*,*) 'Bounds of fxa = ',nz,ubound(fxa(:,:,lbz:),3)
!write(*,*) 'Bounds of fya = ',nz,ubound(fya(:,:,lbz:),3)
!write(*,*) 'Bounds of fza = ',nz,ubound(fza(:,:,lbz:),3)
! This will take care of sync for the forces on all processors
$if ($MPI)
    call mpi_sync_real_array( fxa, 1, MPI_SYNC_DOWN )
    call mpi_sync_real_array( fya, 1, MPI_SYNC_DOWN )
    call mpi_sync_real_array( fza, 1, MPI_SYNC_DOWN )
$endif
call atm_lesgo_write_output() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VTK

!write(*,*) 'Passed Time'

end subroutine atm_lesgo_forcing

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will feed the velocity at all the actuator points into the atm
! This is done by using trilinear interpolation from lesgo
! Force will be calculated based on the velocities and stored on forceField
implicit none

integer, intent(in) :: i ! The turbine number
integer :: m,n,q,j,c
real(rprec), dimension(3) :: velocity
real(rprec), dimension(3) :: xyz    ! Point onto which to interpolate velocity
real(rprec), pointer, dimension(:) :: x,y,z,zw
character(128) :: coord_char,filename, time_char

j=turbineArray(i) % turbineTypeID ! The turbine type ID

! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! The forces need to be initialized at zero for every time-step
turbineArray(i) % bladeForces(:,:,:,:) = 0._rprec
turbineArray(i) % bladeForcesDummy(:,:,:,:) = 0._rprec

! This loop goes through all the blade points and calculates the respective
! body forces then imposes it onto the force field
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
                       
            ! Actuator point onto which to interpolate the velocity
            xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

            ! Non-dimensionalizes the point location 
            xyz=vector_divide(xyz,z_i)
!write(*,*) 'size z = ', size(z)
            ! Interpolate velocities if inside the domain
            if (  z(lbz) <= xyz(3) .and. xyz(3) < z(nz-1) ) then
                velocity(1)=trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                velocity(2)=trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                velocity(3)=trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

                ! This will compute the blade force for the specific point
                call atm_computeBladeForce(i,m,n,q,velocity)

            endif

        enddo
    enddo
enddo
!write(*,*) 'Z is = ',z(1),z(nz-1),z(nz)

end subroutine atm_lesgo_force



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_convolute_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will convolute the forces for each turbine

implicit none

integer, intent(in) :: i
integer :: j, m, n, q, c 

j=turbineArray(i) % turbineTypeID ! The turbine type ID

do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
            ! This will convolute the blade force onto the grid points
            ! affected by the turbines on both grids
            ! Only if the distance is less than specified value
            do c=1,size(forceFieldUV)
                if (distance(forceFieldUV(c) % location,                       &
                   turbineArray(i) % bladePoints(m,n,q,:)) >=                  &
                   3.*turbineArray(i) % epsilon) cycle
                forceFieldUV(c) % force = forceFieldUV(c) % force +            &
                atm_convoluteForce(i, m, n, q, forceFieldUV(c) % location)     &
                *dx*dy*dz*z_i**3./((u_star**2.)) !*z_i
            enddo
            do c=1,size(forceFieldW)
                if (distance(forceFieldW(c) % location,                        &
                   turbineArray(i) % bladePoints(m,n,q,:)) >=                  &
                   3.*turbineArray(i) % epsilon) cycle
                forceFieldW(c) % force = forceFieldW(c) % force +              &
                atm_convoluteForce(i, m, n, q, forceFieldW(c) % location)      &
                *dx*dy*dz*z_i**3./((u_star**2.)) !*z_i
            enddo
        enddo
    enddo
enddo

end subroutine atm_lesgo_convolute_force


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_apply_force()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will apply the blade force onto the CFD grid by using the convolution
! function in the ATM library
implicit none

integer :: c
integer :: i,j,k



! Impose force field onto the flow field variables
! The forces are non-dimensionalized here as well
do c=1,size(forceFieldUV)
    i=forceFieldUV(c) % i
    j=forceFieldUV(c) % j
    k=forceFieldUV(c) % k

    fxa(i,j,k) = fxa(i,j,k) + forceFieldUV(c) % force(1) 
    fya(i,j,k) = fya(i,j,k) + forceFieldUV(c) % force(2) 
enddo

do c=1,size(forceFieldW)
    i=forceFieldW(c) % i
    j=forceFieldW(c) % j
    k=forceFieldW(c) % k

    fza(i,j,k) = fza(i,j,k) + forceFieldW(c) % force(3) 
!write(*,*) 'forceFieldW = ',i,j,k,fza(i,j,k)
enddo


end subroutine atm_lesgo_apply_force

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_write_output()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will write the output of the code from within. No need to access lesgo
! Format is VTK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!VTK
implicit none

real(rprec), pointer, dimension(:) :: x,y,z,zw
character(128) :: coord_char,filename, time_char

! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

!comm_char=char(comm)
write(coord_char,'(i5)') coord
write(time_char,'(i5)') jt_total
!filename='output/data'// trim(time_char) //'.vtk.c'//trim(coord_char)
!call write_3D_Field_VTK(x, y, z, u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w(1:nx,1:ny,1:nz), nx, ny, nz,filename)

! Write plane
call eat_whitespace(time_char,' ')
call eat_whitespace(coord_char,' ')
filename='output/plane_velocity'// trim(time_char) //'.vtk.c'//trim(coord_char)
call write_3D_Field_VTK(x(32:33), y(1:ny), z(1:nz), u(32,1:ny,1:nz), v(32,1:ny,1:nz), w(32,1:ny,1:nz), 1, ny, nz,filename)
!write(*,*) 'Sizes ares = ',size(x),size(y),size(z)
filename='output/plane_force'// trim(time_char) //'.vtk.c'//trim(coord_char)
call write_3D_Field_VTK(x(32:33), y(1:ny), z(1:nz), fxa(32,1:ny,1:nz), fya(32,1:ny,1:nz), fza(32,1:ny,1:nz), 1, ny, nz,filename)

end subroutine atm_lesgo_write_output


end module atm_lesgo_interface





