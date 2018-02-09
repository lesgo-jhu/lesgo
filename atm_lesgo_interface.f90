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
! The MPI management is done only in this section of the code
! This is very code dependent and will have to be modified according to
! the code being used. In this case LESGO has its own MPI details
! Look into mpi_defs.f90 for the details

! Remember to always dimensinoalize the variables from LESGO
! Length is non-dimensionalized by z_i

! Lesgo data used regarding the grid (LESGO)
use param, only : dt ,nx,ny,nz,nz_tot,dx,dy,dz,coord,nproc, z_i, u_star, lbz,  &
                  total_time, jt_total
! nx, ny, nz - nodes in every direction
! z_i - non-dimensionalizing length
! dt - time-step

! These are the forces, and velocities on x,y, and z respectively
use sim_param, only : fxa, fya, fza, u, v, w

! Grid definition (LESGO)
use grid_m, only : grid

! MPI implementation from LESGO
#ifdef PPMPI
  use mpi_defs
  use mpi
  use param, only : ierr, mpi_rprec, comm, coord
#endif

! Interpolating function for interpolating the velocity field to each
! actuator point
use functions, only : trilinear_interp, interp_to_uv_grid

! Actuator Turbine Model module
use atm_base
use actuator_turbine_model
use atm_input_util, only : rprec, turbineArray, turbineModel, eat_whitespace, &
                           atm_print_initialize, updateInterval

! Used for testing time
! use clock_m

implicit none

! Variable for interpolating the velocity in w onto the uv grid
real(rprec), allocatable, dimension(:,:,:) :: w_uv

private
public atm_lesgo_initialize, atm_lesgo_forcing, atm_lesgo_finalize

! This is a list that stores all the points in the domain with a body
! force due to the turbines.
type bodyForce_t
    integer :: c ! Number of cells
    ! i,j,k stores the index for the point in the domain
    integer, allocatable :: ijk(:,:)
    real(rprec), allocatable :: force(:,:) ! Force vector on uv grid
    real(rprec), allocatable :: location(:,:) ! Position vector on uv grid
end type bodyForce_t

! Body force field
type(bodyForce_t), allocatable, target, dimension(:) :: forceFieldUV, forceFieldW

! The very crucial parameter pi
real(rprec), parameter :: pi=acos(-1._rprec)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_initialize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

! Counter to establish number of points which are influenced by body forces
integer ::  m

! Allocate space for the w_uv variable
allocate(w_uv(nx,ny,lbz:nz))

call atm_initialize () ! Initialize the atm (ATM)

! Allocate the body force variables. It is an array with one per turbine.
allocate(forceFieldUV(numberOfTurbines))
allocate(forceFieldW(numberOfTurbines))

do m=1, numberOfTurbines
    call atm_lesgo_findCells(m)
enddo

#ifdef PPMPI
    call mpi_barrier( comm, ierr )

    ! This will create the output files and write initialization to the screen
    if (coord == 0) then
        call atm_initialize_output()
    endif
#else
    call atm_initialize_output()
#endif

end subroutine atm_lesgo_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_finalize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

! Counter for turbines
integer ::  i

! Write if on main node
if (coord == 0) then
    write(*,*) 'Finalizing ATM...'
endif

    ! Loop through all turbines and finalize
    do i = 1, numberOfTurbines
        if (coord == turbineArray(i) % master) then
            call atm_write_restart(i) ! Write the restart file
        endif
    end do

if (coord == 0) then
    write(*,*) 'Done finalizing ATM'
endif

end subroutine atm_lesgo_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_findCells (m)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine finds all the cells that surround the turbines

! The awkward if statements are to only consider points in front and behind
! the turbine without having to

implicit none

! The turbine number
integer, intent(in) :: m

! Counter to establish number of points which are influenced by body forces
integer :: cUV, cW  ! Counters for number of points affected on UV and W grids
integer :: i, j, k

! Vector used to store x, y, z locations
real(rprec), dimension(3) :: vector_point
! These are the pointers to the grid arrays
real(rprec), pointer, dimension(:) :: x,y,z,zw

! Variables for MPI implementation
#ifdef PPMPI
integer :: base_group ! The base group from comm --> MPI_COMM_WORLD (all processors)
integer :: local_group  ! The local group of processors
integer :: member !  (1 or 0) yes or no
integer :: num_of_members  ! total number of members

! List of all the cores that belong to this turbine
! This variable gets allocated for each turbine
integer, allocatable, dimension(:) :: ls_of_cores
#endif

nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! Initialize internal counter to zero
forceFieldUV(m) % c = 0

! This will find all the locations that are influenced by each turbine
! It depends on a sphere centered on the rotor that extends beyond the blades
cUV=0  ! Initialize conuter
cW=0  ! Initialize conuter
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i

            ! Take into account the UV grid
            vector_point(3)=z(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cUV=cUV+1
!~ endif

                end if
                ! Take into account the W grid
                vector_point(3)=zw(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cW=cW+1
!~ endif
                end if
        enddo
    enddo
enddo

! Allocate space for the force fields in UV and W grids
forceFieldUV(m) % c = cUV  ! Counter
allocate(forceFieldUV(m) % force(3,cUV))
allocate(forceFieldUV(m) % location(3,cUV))
allocate(forceFieldUV(m) % ijk(3,cUV))

forceFieldW(m) % c = cW  ! Counter
allocate(forceFieldW(m) % force(3,cW))
allocate(forceFieldW(m) % location(3,cW))
allocate(forceFieldW(m) % ijk(3,cW))

#ifdef PPMPI
call mpi_barrier( comm, ierr )
write(*,*) 'Number of cells being affected by ATM in turbine', m,              &
           ' cUV, cW = ', cUV, cW
call mpi_barrier( comm, ierr )
#endif

cUV=0
cW=0
! Run the same loop and save all variables
! The forceField arrays include all the forces which affect the domain
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i
            vector_point(3)=z(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cUV=cUV+1
                    forceFieldUV(m) % ijk(1,cUV) = i
                    forceFieldUV(m) % ijk(2,cUV) = j
                    forceFieldUV(m) % ijk(3,cUV) = k
                    forceFieldUV(m) % location(1:3,cUV) = vector_point(1:3)
                    forceFieldUV(m) % force(1:3,cUV) = 0_rprec
                endif
!~ endif
            vector_point(3)=zw(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cW=cW+1
                    forceFieldW(m) % ijk(1,cW) = i
                    forceFieldW(m) % ijk(2,cW) = j
                    forceFieldW(m) % ijk(3,cW) = k
                    forceFieldW(m) % location(1:3,cW) = vector_point(1:3)
                    forceFieldW(m) % force(:,cW) = 0_rprec
                endif
!~ endif
        enddo
    enddo
enddo


! MPI distribution
! This will create new communicator for each turbine
#ifdef PPMPI

! Store the base group from the global communicator mpi_comm_world
call MPI_COMM_GROUP(comm, base_group, ierr)

! Assign member
member = 0
! Flag to know if this turbine is operating or not
turbineArray(m) % operate = .FALSE.

! Assign proper values if turbine affects processors in this region
if (cUV > 0 .or. cW >0) then
member = 1
turbineArray(m) % operate = .TRUE.
endif

! Find the total number of processors for each turbine
call mpi_allreduce(member, num_of_members, 1, MPI_INTEGER , MPI_SUM, comm, ierr)

if (turbineArray(m) % operate) then
! Find the master processor for each turbine
    call mpi_allreduce(coord, turbineArray(m) % master, 1, MPI_INTEGER ,       &
                       MPI_MIN, comm, ierr)
else
    ! This is bogus since nz will always be less than number of processors
    ! This is done to ensure that the master is part of the processors
    ! that hold the turbine model
    call mpi_allreduce(nz_tot, turbineArray(m) % master, 1, MPI_INTEGER ,       &
                       MPI_MIN, comm, ierr)
endif

allocate(ls_of_cores(num_of_members))
ls_of_cores(1) = turbineArray(m) % master

! Notice this list is valid only for decomposition in 1 direction
do i = 2, num_of_members
    ls_of_cores(i) = ls_of_cores(i-1) + 1
enddo

! Write if this processor is the master
if (coord == turbineArray(m) % master) then
    write(*,*) 'Master for turbine',m, 'is processor', turbineArray(m) % master
endif

! Create the new communicator and group for this turbine
call MPI_GROUP_INCL(base_group, num_of_members, ls_of_cores, local_group, ierr)
call MPI_COMM_CREATE(comm, local_group, turbineArray(m) % TURBINE_COMM_WORLD, ierr)

if (turbineArray(m) % operate) then
    write(*,*) 'Processor', coord, 'has elements in turbine', m
else
    write(*,*) 'Processor', coord, 'does NOT have elements in turbine', m
endif

    call mpi_barrier( comm, ierr )

#endif

end subroutine atm_lesgo_findCells

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_forcing ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutines calls the update function from the ATM Library
! and calculates the body forces needed in the domain
implicit none

integer :: i

!~ real(rprec) :: integrateNacelleForce, totForce
!~ integer :: c

!~ type(clock_t) :: myClock

!~ call myClock % start()
! Get the velocity from w onto the uv grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)


! Update the blade positions based on the time-step
! Time needs to be dimensionalized
! All processors carry the blade points
!~ call myCock%start_time();
!~ call atm_update(dt*z_i/u_star)

! Loop through all turbines and rotate the blades
do i = 1, numberOfTurbines
    ! If statement is for running code only with the processors on that turbine
    if (turbineArray(i) % operate) then
        ! Time is dimensionalize using velocity and length scale
        call atm_update(i, dt*z_i/u_star, total_time*z_i/u_star)
    endif
enddo

!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Update ', myClock % time


! Only calculate new forces if interval is correct
if ( mod(jt_total-1, updateInterval) == 0) then

    ! Establish all turbine properties as zero
    ! This is essential for paralelization
    do i=1,numberOfTurbines
        turbineArray(i) % bladeForces = 0._rprec
        turbineArray(i) % integratedBladeForces = 0._rprec
        turbineArray(i) % torqueRotor = 0._rprec
        turbineArray(i) % thrust = 0._rprec
        turbineArray(i) % alpha = 0._rprec
        turbineArray(i) % Cd = 0._rprec
        turbineArray(i) % Cl = 0._rprec
        turbineArray(i) % lift = 0._rprec
        turbineArray(i) % drag = 0._rprec
        turbineArray(i) % Vmag = 0._rprec
        turbineArray(i) % windVectors = 0._rprec
        turbineArray(i) % nacelleForce = 0._rprec
        turbineArray(i) % induction_a = 0._rprec
        turbineArray(i) % u_infinity = 0._rprec
        turbineArray(i) % bladeAlignedVectors = 0._rprec
        turbineArray(i) % VelNacelle_sampled = 0._rprec
        turbineArray(i) % VelNacelle_corrected = 0._rprec

        ! If statement is for running code only if grid points affected are in
        ! this processor. If not, no code is executed at all.
!~         if (forceFieldUV(i) % c .gt. 0 .or. forceFieldW(i) % c .gt. 0) then
        if (turbineArray(i) % operate) then

            ! Set body forces to zero
            forceFieldUV(i) % force = 0._rprec
            forceFieldW(i) % force = 0._rprec
            ! Calculate forces for all turbines
            call atm_lesgo_force(i)


        endif

    enddo
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Forces ', myClock % time


!~  call myClock % start()
! This will gather all the blade forces from all processors
#ifdef PPMPI
    ! This will gather all values used in MPI
!~     call mpi_barrier( MPI_COMM_WORLD, ierr )
    call atm_lesgo_mpi_gather()
!~     call mpi_barrier( MPI_COMM_WORLD, ierr )

#endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  MPI Gather ', myClock % time


!~  call myClock % start()
    do i=1,numberOfTurbines
!~         if ( forceFieldUV(i) % c .gt. 0 .or. forceFieldW(i) % c .gt. 0) then
        if (turbineArray(i) % operate) then
            ! Convolute force onto the domain
            call atm_lesgo_convolute_force(i)
        endif

!~         ! Sync the nacelle force
!~         integrateNacelleForce=0.
!~
!~         do c=1,forceFieldUV(i) % c
!~             if (turbineArray(i) % nacelle) then
!~                 integrateNacelleForce = integrateNacelleForce +  &
!~                     forceFieldUV(i) % force(1,c) * dx *dy * dz * z_i**2*u_star**2
!~             endif
!~         enddo

        ! Compute the correction for the Cl coefficient
        call atm_compute_cl_correction(i)

    enddo

!~         totForce=0.
!~         call mpi_allreduce( integrateNacelleForce,  totForce, 1,   &
!~                              mpi_rprec, mpi_sum, comm, ierr)

       !write(*,*) 'Integrated Nacelle Force is: ', integrateNacelleForce
!~         if (coord == 0) then
!~             write(*,*) 'Integrated Total Force is: ', totForce
!~         endif
endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Convolute force ', myClock % time

    ! This will apply body forces onto the flow field if there are forces within
    ! this domain
!~  call myClock % start()
    call atm_lesgo_apply_force()
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Apply force ', myClock % time

!!! Sync the integrated forces (used for debugging)
!do i=1,numberOfTurbines
!    j=turbineArray(i) % turbineTypeID ! The turbine type ID
!    ! Sync all the integrated blade forces
!    turbineArray(i) % bladeVectorDummy=turbineArray(i) % integratedBladeForces
!    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                   &
!                       turbineArray(i) % integratedBladeForces,              &
!                       size(turbineArray(i) % bladeVectorDummy),             &
!                       mpi_rprec, mpi_sum, comm, ierr)


!    if (coord==0) then

!    do q=1, turbineArray(i) % numBladePoints
!        do n=1, turbineArray(i) % numAnnulusSections
!            do m=1, turbineModel(j) % numBl
!                write(*,*) 'blade ',m,'section ',q, 'force ratio', &
!                turbineArray(i) % integratedBladeForces(m,n,q,1) /  &
!                turbineArray(i) % bladeForces(m,n,q,1) , &
!                turbineArray(i) % integratedBladeForces(m,n,q,2) /  &
!                turbineArray(i) % bladeForces(m,n,q,2) , &
!                turbineArray(i) % integratedBladeForces(m,n,q,3) /  &
!                turbineArray(i) % bladeForces(m,n,q,3)
!            enddo
!        enddo
!    enddo
!    endif

!enddo

do i=1, numberOfTurbines
    if (coord == turbineArray(i) % master) then
    !~  call myClock % start()

        call atm_output(i, jt_total, total_time*z_i/u_star)
    !~     call myClock % stop()
    !~     write(*,*) 'coord ', coord, '  Output ', myClock % time
    endif
enddo

#ifdef PPMPDI
! Make sure all processors stop wait for the output to be completed
call mpi_barrier( comm, ierr )
#endif

end subroutine atm_lesgo_forcing


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Complie this subroutines only if MPI will be used
#ifdef PPMPI

subroutine atm_lesgo_mpi_gather()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will gather the necessary outputs from the turbine models
! so all processors have acces to it. This is done by means of all reduce SUM
implicit none
integer :: i
real(rprec) :: torqueRotor, thrust, VelNacelle_sampled, VelNacelle_corrected
real(rprec), dimension(3) :: nacelleForce

! Pointer for MPI communicator
integer, pointer :: TURBINE_COMMUNICATOR

do i=1,numberOfTurbines

    ! Only do MPI sums if processors are operating in this turbine
    if (turbineArray(i) % operate) then

        TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

        turbineArray(i) % bladeVectorDummy = turbineArray(i) % bladeForces
        ! Sync all the blade forces
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeForces,                      &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync bladeAlignedVectors
        turbineArray(i) % bladeVectorDummy =                                   &
                              turbineArray(i) % bladeAlignedVectors(:,:,:,1,:)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeAlignedVectors(:,:,:,1,:),   &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        turbineArray(i) % bladeVectorDummy =                                   &
                              turbineArray(i) % bladeAlignedVectors(:,:,:,2,:)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeAlignedVectors(:,:,:,2,:),   &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        turbineArray(i) % bladeVectorDummy =                                   &
                              turbineArray(i) % bladeAlignedVectors(:,:,:,3,:)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeAlignedVectors(:,:,:,3,:),   &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)


        ! Sync alpha
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % alpha
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % alpha,                            &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        ! Sync lift
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % lift
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % lift,                             &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        ! Sync drag
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % drag
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % drag,                             &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        ! Sync Cl
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cl
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % Cl,                               &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync Cd
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cd
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % Cd,                               &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync Vmag
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % Vmag
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % Vmag,                             &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync axialForce
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % axialForce
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % axialForce,                       &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync tangentialForce
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % tangentialForce
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % tangentialForce,                                 &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync wind Vectors (Vaxial, Vtangential, Vradial)
        turbineArray(i) % bladeVectorDummy = turbineArray(i) %                 &
                                             windVectors(:,:,:,1:3)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % windVectors(:,:,:,1:3),                                 &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync induction factor
        turbineArray(i) % bladeScalarDummy = turbineArray(i) %                 &
                                             induction_a(:,:,:)
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % induction_a(:,:,:),                                 &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync u infinity
        turbineArray(i) % bladeScalarDummy = turbineArray(i) %                 &
                                             u_infinity(:,:,:)
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % u_infinity(:,:,:),                                 &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Store the torqueRotor.
        ! Needs to be a different variable in order to do MPI Sum
        torqueRotor=turbineArray(i) % torqueRotor
        thrust=turbineArray(i) % thrust
        nacelleForce=turbineArray(i) % nacelleForce
        VelNacelle_sampled=turbineArray(i) % VelNacelle_sampled
        VelNacelle_corrected=turbineArray(i) % VelNacelle_corrected

        ! Sum all the individual torqueRotor from different blade points
        call mpi_allreduce( torqueRotor, turbineArray(i) % torqueRotor,        &
                           1, mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sum all the individual thrust from different blade points
        call mpi_allreduce( thrust, turbineArray(i) % thrust,                  &
                           1, mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync the nacelle force
        call mpi_allreduce( nacelleForce, turbineArray(i) % nacelleForce,      &
                           size(turbineArray(i) % nacelleForce), mpi_rprec,    &
                                mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync the nacelle sampled velocity
        call mpi_allreduce( VelNacelle_sampled,                                &
                               turbineArray(i) % VelNacelle_sampled,  1,       &
                                mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync the nacelle corrected velocity
        call mpi_allreduce( VelNacelle_corrected,                              &
                               turbineArray(i) % VelNacelle_corrected, 1,      &
                                mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    endif
enddo
end subroutine atm_lesgo_mpi_gather
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will feed the velocity at all the actuator points into the atm
! This is done by using trilinear interpolation from lesgo
! Force will be calculated based on the velocities and stored on forceField
implicit none

integer, intent(in) :: i ! The turbine number
integer :: m,n,q,j
! mpi_velocity only used for Spalart method
real(rprec), dimension(3) :: velocity, mpi_velocity
real(rprec), dimension(3) :: xyz    ! Point onto which to interpolate velocity
real(rprec), pointer, dimension(:) :: x,y,z,zw

! The MPI turbine communcator
integer, pointer :: TURBINE_COMM

TURBINE_COMM => turbineArray(i) % TURBINE_COMM_WORLD

j=turbineArray(i) % turbineTypeID ! The turbine type ID

! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

if (turbineArray(i) % sampling == 'Spalart') then
    ! This loop goes through all the blade points and calculates the respective
    ! body forces then imposes it onto the force field
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl

                ! Actuator point onto which to interpolate the velocity
                xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

                velocity = 0._rprec
                mpi_velocity = 0._rprec

                call atm_lesgo_compute_spalart_u(i, xyz, velocity)

                mpi_velocity = velocity

                ! Complie this subroutines only if MPI will be used
#ifdef PPMPI
!~                     call mpi_barrier( TURBINE_COMM, ierr )
                    ! Sync all the blade forces
                    call mpi_allreduce(mpi_velocity, velocity, size(velocity), &
                           mpi_rprec, mpi_sum, TURBINE_COMM , ierr)
#endif

                ! This will compute the blade force for the specific point
                if (  z(1) <= xyz(3)/z_i .and. xyz(3)/z_i < z(nz) ) then
                    call atm_computeBladeForce(i,m,n,q,velocity)
                else
                    velocity = 0._rprec
                endif

            enddo
        enddo
    enddo


else if (turbineArray(i) % sampling == 'atPoint') then
    ! This loop goes through all the blade points and calculates the respective
    ! body forces then imposes it onto the force field
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl

                ! Actuator point onto which to interpolate the velocity
                xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

                ! Non-dimensionalizes the point location
                xyz=xyz/z_i

                ! Interpolate velocities if inside the domain
                if (  z(1) <= xyz(3) .and. xyz(3) < z(nz) ) then
                    velocity(1)=                                               &
                    trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                    velocity(2)=                                               &
                    trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                    velocity(3)=                                               &
                    trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

                    ! This will compute the blade force for the specific point
                    call atm_computeBladeForce(i,m,n,q,velocity)

                endif

            enddo
        enddo
    enddo
endif

    ! Calculate Nacelle force
    if (turbineArray(i) % nacelle) then
        xyz=turbineArray(i) % nacelleLocation
        xyz=xyz/z_i
        if (  z(1) <= xyz(3) .and. xyz(3) < z(nz) ) then

            velocity(1)=                                                   &
            trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
            velocity(2)=                                                   &
            trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
            velocity(3)=                                                   &
            trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

            call atm_computeNacelleForce(i,velocity)

        endif
    endif

!~ ! Compute the correction for the Cl coefficient
!~ call atm_compute_cl_correction(i)


end subroutine atm_lesgo_force

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_compute_Spalart_u(i, xyz, velocity)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will calculate the sampling velocity using the proposed method
! from Spalart
! n turbine number
! xyz actuator point position vector
! velocity reference velocity for computing lift and drag

implicit none

integer, intent(in) :: i
real(rprec), intent(in) :: xyz(3)
real(rprec), intent(inout) :: velocity(3)

integer :: c, m, n, q

! Pointers for mesh
real(rprec), pointer, dimension(:) :: z,zw

! Test for time optimization
real(rprec) :: dist, a(3), projectradius, epsilon

nullify(z,zw)
z => grid % z
zw => grid % zw

! Value of epsilon
epsilon=turbineArray(i) % epsilon

! Projection radius
projectradius = turbineArray(i) % projectionRadius

! Set the velocity to zero
velocity = 0._rprec


do c=1,forceFieldUV(i) % c

    a = forceFieldUV(i) %  location(1:3, c)
    m = forceFieldUV(i) %  ijk(1, c)
    n = forceFieldUV(i) %  ijk(2, c)
    q = forceFieldUV(i) %  ijk(3, c)

    dist=((a(1)-xyz(1))**2+(a(2)-xyz(2))**2+(a(3)-xyz(3))**2)**0.5
    if (dist .le. projectradius * z_i) then
        if ( z(1) <= a(3)/z_i .and. a(3)/z_i < z(nz)) then

        ! The value of the kernel. This is the actual smoothing function
        velocity(1) = velocity(1) + u(m,n,q) * exp(-(dist/epsilon)**2)    &
                                 /  ((epsilon**3.)*(pi**1.5))
        velocity(2) = velocity(2) + v(m,n,q) * exp(-(dist/epsilon)**2)    &
                                 /  ((epsilon**3.)*(pi**1.5))
        endif
    endif
enddo

do c=1,forceFieldW(i) % c
    a = forceFieldW(i) %  location(1:3, c)
    m = forceFieldW(i) %  ijk(1, c)
    n = forceFieldW(i) %  ijk(2, c)
    q = forceFieldW(i) %  ijk(3, c)

    dist=((a(1)-xyz(1))**2+(a(2)-xyz(2))**2+(a(3)-xyz(3))**2)**0.5

    if (dist .le. projectradius) then
        if ( z(1) <= a(3)/z_i .and. a(3)/z_i < z(nz)) then

        ! The value of the kernel. This is the actual smoothing function
        velocity(3) = velocity(3) + w(m,n,q) * exp(-(dist/epsilon)**2)    &
                                 /  ((epsilon**3.)*(pi**1.5))
        endif
    endif
enddo

velocity = velocity * u_star * z_i * dx * z_i * dy *z_i * dz

end subroutine atm_lesgo_compute_Spalart_u

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_convolute_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will convolute the forces for each turbine

implicit none

!~ type(clock_t) :: myClock

integer, intent(in) :: i
integer :: j, m, n, q, c,mmend,nnend,qqend

integer :: ii, jj, kk  ! Indices for lesgo fields

! Test for time optimization
real(rprec) :: dist,a(3),b(3),projectradius,epsilon,const1,const2,const3
real(rprec) :: nacelleEpsilon

! Variables for convolution force
real(rprec) :: kernel, force(3)

! Pointers for the turbineArray quantities
real(rprec), pointer, dimension(:,:,:,:) :: bladeForces, bladePoints

real(rprec), pointer, dimension(:,:) :: bodyForceUV, bodyForceW

nullify(bladeForces)
nullify(bladePoints)
nullify(bodyForceUV)
nullify(bodyForceW)

bladeForces => turbineArray(i) % bladeForces
bladePoints => turbineArray(i) % bladePoints

bodyForceUV => forceFieldUV(i) % force
bodyForceW =>  forceFieldW(i) % force

!real(rprec) :: dummyForce(3)  ! Debugging

j=turbineArray(i) % turbineTypeID ! The turbine type ID

! This will convolute the blade force onto the grid points
! affected by the turbines on both grids
! Only if the distance is less than specified value
mmend=turbineModel(j) % numBl
nnend=turbineArray(i) % numAnnulusSections
qqend=turbineArray(i) % numBladePoints
projectradius=turbineArray(i) % projectionRadius
epsilon=turbineArray(i) % epsilon
nacelleEpsilon = turbineArray(i) % nacelleEpsilon
const1=1./ ((epsilon**3.)*(pi**1.5))
const2= z_i/(u_star**2.)
const3=const1*const2

! Body Force implementation using velocity sampling at the actuator point
if (turbineArray(i) % sampling == 'atPoint') then

    !~  call myClock % start()
    do c=1,forceFieldUV(i) % c
        a= forceFieldUV(i) %  location(1:3,c)
        force=0._rprec

        ! Blade forces
        do m=1, mmend
            do n=1, nnend
               do q=1, qqend

                    b= bladePoints(m,n,q,:)
                    dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5

                    if (dist .le. projectradius) then
                    ! The value of the kernel. This is the actual smoothing function
                     force(1:2) = force(1:2) +  bladeForces(m,n,q,1:2) *exp(-(dist/epsilon)**2)
                    endif

                enddo
            enddo
        enddo
        force(1:2)=force(1:2)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5
    !~         if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel=exp(-(dist/nacelleEpsilon)**2.) / ( (nacelleEpsilon**3.)*(pi**1.5) )
                !write(*,*) 'kernel Value= ', kernel
                force(1:2) = force(1:2)+turbineArray(i) % nacelleForce(1:2) *  &
                             kernel *const2
    !~          integrateNacelleForce=integrateNacelleForce+force(1) * dx *dy * dz * z_i**3

    !~         endif
        endif


        bodyForceUV(1:2,c) = force(1:2)
    !~     if (abs(bodyForceUV(1,c)) .gt. 0) then
    !~                 write(*,*) 'bodyForceUV is: ', bodyForceUV(1,c)
    !~     endif
    enddo


    do c=1,forceFieldW(i) % c
        a= forceFieldW(i) %  location(1:3,c)
        force=0._rprec

        ! Blade forces
        do m=1,mmend
            do n=1,nnend
               do q=1,qqend

                    b= bladePoints(m,n,q,:)
                    dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5

                    if (dist .le. projectradius) then
                    ! The value of the kernel. This is the actual smoothing function
                    force(3) = force(3) +  bladeForces(m,n,q,3) &
                               *exp(-(dist/epsilon)**2)
                    endif

                enddo
            enddo
        enddo
        force(3)=force(3)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5
            if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel=exp(-(dist/nacelleEpsilon)**2.)/(nacelleepsilon**3.*pi**1.5)
                force(3) = force(3)+turbineArray(i) % nacelleForce(3) *           &
                           kernel *const2
            endif
        endif

        bodyForceW(3,c) = force(3)
    enddo

! The Spalart method uses the local velocity field.
! For this reason it needs to be done explicitly in this module
! and cannot be generally coded from the actuator_turbine_model module
elseif (turbineArray(i) % sampling == 'Spalart') then

    !~  call myClock % start()
    do c=1,forceFieldUV(i) % c
        a= forceFieldUV(i) %  location(1:3,c)
        force=0._rprec
        ! Indices for velocity field
        ii = forceFieldUV(i) % ijk(1,c)
        jj = forceFieldUV(i) % ijk(2,c)
        kk = forceFieldUV(i) % ijk(3,c)

        ! Blade forces
        do m=1, mmend
            do n=1, nnend
               do q=1, qqend

                    b= bladePoints(m,n,q,:)
                    dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5

                    if (dist .le. projectradius) then
                        ! The value of the kernel.
                        ! This is the actual smoothing function
                        ! Divide by velocity magnitude
                         force(1) = force(1) +  bladeForces(m,n,q,1) *         &
                                      exp(-(dist/epsilon)**2)                  &
                         / (turbineArray(i) % Vmag(m,n,q)) *                     &
                         ( u(ii,jj,kk)  * u_star +                             &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))

                         force(2) = force(2) +  bladeForces(m,n,q,2) *         &
                                      exp(-(dist/epsilon)**2)                  &
                         / turbineArray(i) % Vmag(m,n,q) *                     &
                         ( v(ii,jj,kk) * u_star +                              &
                         turbineArray(i) % bladeAlignedVectors(m,n,q,2,2) *    &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))

                    endif

                enddo
            enddo
        enddo
        force(1:2)=force(1:2)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5
    !~         if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-(dist/nacelleEpsilon)**2.) /                       &
                         ( (nacelleEpsilon**3.)*(pi**1.5) )
                !write(*,*) 'kernel Value= ', kernel
                force(1:2) = force(1:2)+turbineArray(i) % nacelleForce(1:2) *  &
                             kernel *const2
    !~          integrateNacelleForce=integrateNacelleForce+force(1) * dx *dy * dz * z_i**3

    !~         endif
        endif


        bodyForceUV(1:2,c) = force(1:2)
    enddo


    do c=1,forceFieldW(i) % c
        a= forceFieldW(i) %  location(1:3,c)
        force=0._rprec
        ! Indices for velocity field
        ii = forceFieldW(i) % ijk(1,c)
        jj = forceFieldW(i) % ijk(2,c)
        kk = forceFieldW(i) % ijk(3,c)

        ! Blade forces
        do m=1,mmend
            do n=1,nnend
               do q=1,qqend

                    b= bladePoints(m,n,q,:)
                    dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5

                    if (dist .le. projectradius) then
                        ! The value of the kernel.
                        ! This is the actual smoothing function
                        force(3) = force(3) +  bladeForces(m,n,q,3) *          &
                                   exp(-(dist/epsilon)**2)                     &
                         / turbineArray(i) % Vmag(m,n,q) *                     &
                         ( w(ii,jj,kk) * u_star +                              &
                         turbineArray(i) % bladeAlignedVectors(m,n,q,2,3) *    &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))
                    endif

                enddo
            enddo
        enddo
        force(3)=force(3)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist=((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)**0.5
            if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel=exp(-(dist/nacelleEpsilon)**2.)/(nacelleepsilon**3.*pi**1.5)
                force(3) = force(3)+turbineArray(i) % nacelleForce(3) *           &
                           kernel *const2
            endif
        endif

        bodyForceW(3,c) = force(3)
    enddo

endif

end subroutine atm_lesgo_convolute_force


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_apply_force()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will apply the blade force onto the CFD grid by using the convolution
! function in the ATM library
implicit none

integer :: c,m
integer :: i,j,k

do m=1, numberOfTurbines

    if (turbineArray(m) % operate) then
        ! Impose force field onto the flow field variables
        ! The forces are non-dimensionalized here as well
        do c=1,forceFieldUV(m) % c
            i=forceFieldUV(m) % ijk(1,c)
            j=forceFieldUV(m) % ijk(2,c)
            k=forceFieldUV(m) % ijk(3,c)

            fxa(i,j,k) = fxa(i,j,k) + forceFieldUV(m) % force(1,c)
            fya(i,j,k) = fya(i,j,k) + forceFieldUV(m) % force(2,c)

        enddo

        do c=1,forceFieldW(m) % c
            i=forceFieldW(m) % ijk(1,c)
            j=forceFieldW(m) % ijk(2,c)
            k=forceFieldW(m) % ijk(3,c)

            fza(i,j,k) = fza(i,j,k) + forceFieldW(m) % force(3,c)

        enddo
    endif
enddo

end subroutine atm_lesgo_apply_force


end module atm_lesgo_interface





