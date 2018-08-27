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
module actuator_turbine_model
!*******************************************************************************
! This module has the subroutines to provide all calculations for use in the 
! actuator turbine model (ATM)

! Imported modules
use atm_base ! Include basic types and precission of real numbers

use atm_input_util ! Utilities to read input files

implicit none

! Declare everything private except for subroutine which will be used
private 
public :: atm_initialize, numberOfTurbines,                                    &
          atm_computeBladeForce, atm_update,                                   &
          vector_add, vector_divide, vector_mag, distance,                     &
          atm_output, atm_process_output,                                      &
          atm_initialize_output, atm_computeNacelleForce, atm_write_restart,   &
          atm_compute_cl_correction

! The very crucial parameter pi
real(rprec), parameter :: pi=acos(-1._rprec) 

! These are used to do unit conversions
real(rprec) :: degRad = pi/180._rprec ! Degrees to radians conversion
real(rprec) :: rpmRadSec =  pi/30._rprec ! Set the revolutions/min to radians/s 

logical :: pastFirstTimeStep ! Establishes if we are at the first time step

! Subroutines for the actuator turbine model 
! All suboroutines names start with (atm_) 
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize()
! This subroutine initializes the ATM. It calls the subroutines in
! atm_input_util to read the input data and creates the initial geometry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer :: i
logical :: file_exists

pastFirstTimeStep=.false. ! The first time step not reached yet

write(*,*) 'Reading Actuator Turbine Model Input...'
call read_input_conf()  ! Read input data
write(*,*) 'Done Reading Actuator Turbine Model Input'
do i = 1,numberOfTurbines
    inquire(file = "./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/actuatorPoints", exist=file_exists)

    ! Creates the ATM points defining the geometry
    call atm_create_points(i) 
    ! This will create the first yaw alignment
    turbineArray(i) % deltaNacYaw = turbineArray(i) % nacYaw
    call atm_yawNacelle(i)

    if (file_exists .eqv. .true.) then
        write(*,*) 'Reading bladePoints from Previous Simulation'
        call atm_read_actuator_points(i)
    endif

    call atm_calculate_variables(i) ! Calculates variables depending on input

    inquire(file = "./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/restart", exist=file_exists)

    if (file_exists .eqv. .true.) then
        write(*,*) 'Reading Turbine Properties from Previous Simulation'
        call atm_read_restart(i)
    endif

end do

pastFirstTimeStep=.true. ! Past the first time step

end subroutine atm_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_read_actuator_points(i)
! This subroutine reads the location of the actuator points
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original location
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i ! Indicates the turbine number

integer :: j, m, n, q

j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=1, file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                  "/actuatorPoints", action='read')

do m=1, turbineModel(j) % numBl
    do n=1, turbineArray(i) %  numAnnulusSections
        do q=1, turbineArray(i) % numBladePoints
            read(1,*) turbineArray(i) % bladePoints(m,n,q,:)
        enddo
    enddo
enddo

close(1)

end subroutine atm_read_actuator_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_read_restart(i)
! This subroutine reads the rotor speed
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original omega
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i  ! Indicates the turbine number

! Open the file at the last line (append)
open( unit=1, file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/restart", action='read') !, position='append')

! Bring the pointer to the last line
!~ backspace 1

! Read past the first line
read(1,*)

! Read the restart variables
read(1,*) turbineArray(i) % rotSpeed
read(1,*) turbineArray(i) % torqueGen
read(1,*) turbineArray(i) % torqueRotor
read(1,*) turbineArray(i) % u_infinity
read(1,*) turbineArray(i) % induction_a
read(1,*) turbineArray(i) % PitchControlAngle
read(1,*) turbineArray(i) % IntSpeedError
read(1,*) turbineArray(i) % nacYaw
read(1,*) turbineArray(i) % rotorApex
read(1,*) turbineArray(i) % uvShaft 
close(1)

write(*,*) ' RotSpeed Value from previous simulation is ',                     &
                turbineArray(i) % rotSpeed
write(*,*) ' torqueGen Value from previous simulation is ',                    &
                turbineArray(i) % torqueGen
write(*,*) ' torqueRotor Value from previous simulation is ',                  &
                turbineArray(i) % torqueRotor
write(*,*) ' PitchControlAngle Value from previous simulation is ',            &
                turbineArray(i) % PitchControlAngle
write(*,*) ' IntSpeedError Value from previous simulation is ',                &
                turbineArray(i) % IntSpeedError
write(*,*) ' Yaw Value from previous simulation is ',                          &
                turbineArray(i) % nacYaw
write(*,*) ' Rotor Apex Value from previous simulation is ',                   &
                turbineArray(i) % rotorApex
write(*,*) ' uvShaft Value from previous simulation is ',                      &
                turbineArray(i) % uvShaft


end subroutine atm_read_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_write_restart(i)
! This subroutine reads the rotor speed
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original omega
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: pointsFile=787 ! File to write the actuator points
integer :: restartFile=21 ! File to write restart data
integer j, m,n,q ! counters

! Open the file 
open( unit=restartFile, file="./turbineOutput/"//                              &
            trim(turbineArray(i) % turbineName)//"/restart", status="replace")

write(restartFile,*) 'RotSpeed ', 'torqueGen ', 'torqueRotor ', 'u_infinity ', &
                     'induction_a ', 'PitchControlAngle ', 'IntSpeedError ',   &
                     'nacYaw ', 'rotorApex ', 'uvShaft'
! Store the rotSpeed value 
write(restartFile,*) turbineArray(i) % rotSpeed
write(restartFile,*) turbineArray(i) % torqueGen
write(restartFile,*) turbineArray(i) % torqueRotor
write(restartFile,*) turbineArray(i) % u_infinity
write(restartFile,*) turbineArray(i) % induction_a
write(restartFile,*) turbineArray(i) % PitchControlAngle
write(restartFile,*) turbineArray(i) % IntSpeedError
write(restartFile,*) turbineArray(i) % nacYaw
write(restartFile,*) turbineArray(i) % rotorApex
write(restartFile,*) turbineArray(i) % uvShaft
close(restartFile)

! Write the actuator points at every time-step regardless
j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=pointsFile, status="replace", file="./turbineOutput/"//              &
                      trim(turbineArray(i) % turbineName)//"/actuatorPoints")

do m=1, turbineModel(j) % numBl
    do n=1, turbineArray(i) %  numAnnulusSections
        do q=1, turbineArray(i) % numBladePoints
            ! A new file will be created each time-step with the proper
            ! location of the blades
            write(pointsFile,*) turbineArray(i) % bladePoints(m,n,q,:)
        enddo
    enddo
enddo

close(pointsFile)


end subroutine atm_write_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize_output()
! This subroutine initializes the output files for the ATM
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
logical :: file_exists
integer :: i

! Write to the screen output start
call atm_print_initialize()

do i = 1,numberOfTurbines

    inquire(file="./turbineOutput/"//                                   &
                     trim(turbineArray(i) % turbineName),EXIST=file_exists)

    if (file_exists .eqv. .false.) then

        ! Create turbineOutput directory
        call system("mkdir -vp turbineOutput/"//                               &
                     trim(turbineArray(i) % turbineName)) 

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/power") 
        write(1,*) 'time PowerRotor powerGen '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/thrust") 
        write(1,*) 'time thrust '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/RotSpeed") 
        write(1,*) 'time RotSpeed'
        close(1)
        
        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Yaw") 
        write(1,*) 'time deltaNacYaw NacYaw'
        close(1)
        
        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/lift")
        write(1,*) 'turbineNumber bladeNumber '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/drag")
        write(1,*) 'turbineNumber bladeNumber '
        close(1)
        
        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Cl")
        write(1,*) 'turbineNumber bladeNumber Cl'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Cd")
        write(1,*) 'turbineNumber bladeNumber Cd'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/alpha")
        write(1,*) 'turbineNumber bladeNumber alpha'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Vrel")
        write(1,*) 'turbineNumber bladeNumber Vrel'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Vaxial")
        write(1,*) 'turbineNumber bladeNumber Vaxial'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Vtangential")
        write(1,*) 'turbineNumber bladeNumber Vtangential'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/tangentialForce")
        write(1,*) 'turbineNumber bladeNumber tangentialForce'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/axialForce")
        write(1,*) 'turbineNumber bladeNumber axialForce'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/nacelle")
        write(1,*) 'time Velocity-no-correction Velocity-w-correction'
        close(1)

    endif
enddo

end subroutine atm_initialize_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_create_points(i)
! This subroutine generate the set of blade points for each turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
integer :: m ! Indicates the blade point number
integer :: n ! Indicates number of actuator section
integer :: k
real(rprec), dimension (3) :: root ! Location of rotor apex
real(rprec) :: beta ! Difference between coning angle and shaft tilt
real(rprec) :: dist ! Distance from each actuator point
integer,     pointer :: numBladePoints
integer,     pointer :: numBl
integer,     pointer :: numAnnulusSections
real(rprec), pointer :: NacYaw             
real(rprec),  pointer :: db(:)
real(rprec),  pointer :: bladePoints(:,:,:,:)
real(rprec),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: azimuth
real(rprec), pointer :: rotSpeed
real(rprec), pointer :: ShftTilt
real(rprec), pointer :: PreCone
real(rprec), pointer :: towerShaftIntersect(:)
real(rprec), pointer :: baseLocation(:)
real(rprec), pointer :: TowerHt
real(rprec), pointer :: Twr2Shft
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: uvShaftDir
real(rprec), pointer :: uvShaft(:)
real(rprec), pointer :: uvTower(:)
real(rprec), pointer :: TipRad
real(rprec), pointer :: HubRad
real(rprec), pointer :: annulusSectionAngle
real(rprec), pointer :: solidity(:,:,:)

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Variables to be used locally. They are stored in local variables within the 
! subroutine for easier code following. The values are then passed to the 
! proper type
numBladePoints => turbineArray(i) % numBladePoints
numBl=>turbineModel(j) % numBl
numAnnulusSections=>turbineArray(i) % numAnnulusSections

! Allocate variables depending on specific turbine properties and general
! turbine model properties
allocate(turbineArray(i) % db(numBladePoints))

allocate(turbineArray(i) % bladePoints(numBl, numAnnulusSections, &
         numBladePoints,3))
         
allocate(turbineArray(i) % bladeRadius(numBl,numAnnulusSections,numBladePoints))  

allocate(turbineArray(i) % solidity(numBl,numAnnulusSections,numBladePoints))  

! Assign Pointers turbineArray denpendent (i)
db=>turbineArray(i) % db
bladePoints=>turbineArray(i) % bladePoints
bladeRadius=>turbineArray(i) % bladeRadius
solidity=>turbineArray(i) % solidity
azimuth=>turbineArray(i) % azimuth
rotSpeed=>turbineArray(i) % rotSpeed
towerShaftIntersect=>turbineArray(i) % towerShaftIntersect
baseLocation=>turbineArray(i) % baseLocation
uvShaft=>turbineArray(i) % uvShaft
uvTower=>turbineArray(i) % uvTower
rotorApex=>turbineArray(i) % rotorApex
uvShaftDir=>turbineArray(i) % uvShaftDir
nacYaw=>turbineArray(i) % nacYaw
numAnnulusSections=>turbineArray(i) % numAnnulusSections
annulusSectionAngle=>turbineArray(i) % annulusSectionAngle


! Assign Pointers turbineModel (j)
ShftTilt=>turbineModel(j) % ShftTilt
preCone=>turbineModel(j) % preCone
TowerHt=>turbineModel(j) % TowerHt
Twr2Shft=> turbineModel(j) % Twr2Shft
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
TipRad=>turbineModel(j) % TipRad
HubRad=>turbineModel(j) % HubRad
PreCone=>turbineModel(j) %PreCone

!!-- Do all proper conversions for the required variables
! Convert nacelle yaw from compass directions to the standard convention
!~ call atm_compassToStandard(nacYaw)

! The nacelle Yaw is set to 0 deg in the streamwise direction
write(*,*) 'NacYaw is ', nacYaw
! Turbine specific
azimuth = degRad * azimuth
rotSpeed = rpmRadSec * rotSpeed
nacYaw = degRad * nacYaw

! Turbine model specific
shftTilt = degRad * shftTilt 
preCone = degRad * preCone

! Calculate tower shaft intersection and rotor apex locations. (The i-index is 
! at the turbine array level for each turbine and the j-index is for each type 
! of turbine--if all turbines are the same, j- is always 0.)  The rotor apex is
! not yet rotated for initial yaw that is done below.
towerShaftIntersect = turbineArray(i) % baseLocation
towerShaftIntersect(3) = towerShaftIntersect(3) + TowerHt + Twr2Shft
rotorApex = towerShaftIntersect
rotorApex(1) = rotorApex(1) +  (OverHang + UndSling) * cos(ShftTilt)
rotorApex(3) = rotorApex(3) +  (OverHang + UndSling) * sin(ShftTilt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Create Nacelle Point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
turbineArray(i) % nacelleLocation = rotorApex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Create the first set of actuator points                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the vector along the shaft pointing in the direction of the wind
uvShaftDir = OverHang / abs( OverHang )

! Define the vector along the shaft pointing in the direction of the wind
uvShaft = vector_add(rotorApex , - towerShaftIntersect)
uvShaft = vector_divide(uvShaft, vector_mag(uvShaft))
uvShaft = vector_multiply( uvShaft,uvShaftDir)

! Define vector aligned with the tower pointing from the ground to the nacelle
uvTower = vector_add(towerShaftIntersect, - baseLocation)
uvTower = vector_divide( uvTower, vector_mag(uvTower))

! Define thickness of each blade section
do k=1, numBladePoints
    db(k) = (TipRad - HubRad)/(numBladePoints)
enddo

! This creates the first set of points
do k=1, numBl
    root = rotorApex
    beta = PreCone - ShftTilt
    root(1)= root(1) + HubRad*sin(beta)
    root(3)= root(3) + HubRad*cos(beta)
!~     dist = HubRad
    dist = 0.
    
    ! Number of blade points for the first annular section
    do m=1, numBladePoints
        dist = dist + 0.5*(db(m))
        bladePoints(k,1,m,1) = root(1) + dist*sin(beta)
        bladePoints(k,1,m,2) = root(2)
        bladePoints(k,1,m,3) = root(3) + dist*cos(beta)
        do n=1,numAnnulusSections
!~             bladeRadius(k,n,m) = dist
            bladeRadius(k,n,m) = dist + HubRad
            solidity(k,n,m)=1./numAnnulusSections
        enddo
        dist = dist + 0.5*db(m)
    enddo

    ! If there are more than one blade create the points of other blades by
    ! rotating the points of the first blade
    if (k > 1) then
        do m=1, numBladePoints
            bladePoints(k,1,m,:)=rotatePoint(bladePoints(k,1,m,:), rotorApex,  &
            uvShaft,(360.0/NumBl)*(k-1)*degRad)
        enddo
    endif

    ! Rotate points for all the annular sections
    if (numAnnulusSections .lt. 2) cycle ! Cycle if only one section (ALM)
    do n=2, numAnnulusSections
        do m=1, numBladePoints
            bladePoints(k,n,m,:) =                                             &
            rotatePoint(bladePoints(k,1,m,:), rotorApex,                       &
            uvShaft,(annulusSectionAngle/(numAnnulusSections))*(n-1.)*degRad)
        enddo
    enddo
enddo

! Apply the first rotation
turbineArray(i) % deltaAzimuth = azimuth
call atm_rotateBlades(i)

end subroutine atm_create_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_update(i, dt, time)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
real(rprec), intent(in) :: dt                            ! Time step
real(rprec), intent(in) :: time                          ! Simulation time

! Rotate the blades
call atm_computeRotorSpeed(i,dt) 
call atm_rotateBlades(i)

call atm_control_yaw(i, time)

!~ if(pastFirstTimeStep) then
    ! Compute the lift correction for this case
!~     call atm_compute_cl_correction(i)
!~ endif

end subroutine atm_update

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_control_yaw(i, time)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
real(rprec), intent(in) :: time                          ! Simulation time

integer :: j                                             ! Turbine Type ID

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Will calculate the yaw angle  and yaw the nacelle (from degrees to radians)
if ( turbineModel(j) % YawControllerType == "timeYawTable" ) then
    turbineArray(i) % deltaNacYaw = interpolate(time,                          &
    turbineModel(j) % yaw_time(:), turbineModel(j) % yaw_angle(:)) * degRad -  &
    turbineArray(i) % NacYaw

    ! Yaw only if angle is greater than given tolerance
    if (abs(turbineArray(i) % deltaNacYaw) > 0.00000001) then
        call atm_yawNacelle(i)
    endif

!~     write(*,*) 'Delta Yaw is', turbineArray(i) % deltaNacYaw/degRad
!~     write(*,*) 'Nacelle Yaw is', turbineArray(i) % NacYaw/degRad
endif

end subroutine atm_control_yaw

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeRotorSpeed(i,dt)
! This subroutine rotates the turbine blades 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                     ! Turbine number
real(rprec), intent(in) :: dt                ! time step
integer :: j                                 ! Turbine type
!real(rprec) :: deltaAzimuth                  ! Angle of rotation

! Pointers to turbineArray(i)
real(rprec), pointer :: rotSpeed, torqueGen, torqueRotor, fluidDensity

! Pointers to turbineModel(j)
real(rprec), pointer :: GBRatio, CutInGenSpeed, RatedGenSpeed
real(rprec), pointer :: Region2StartGenSpeed, Region2EndGenSpeed
real(rprec), pointer :: CutInGenTorque,RateLimitGenTorque,RatedGenTorque
real(rprec), pointer :: KGen,TorqueControllerRelax, DriveTrainIner

! Other variables to be used
real(rprec) :: torqueGenOld, genSpeed, dGenSpeed, Region2StartGenTorque
real(rprec) :: torqueSlope, Region2EndGenTorque

! Pitch Controller values
real(rprec) :: GK, KP, KI, SpeedError
real(rprec), pointer :: IntSpeedError, PitchControlAngle

j=turbineArray(i) % turbineTypeID

rotSpeed=>turbineArray(i) % rotSpeed
torqueGen=>turbineArray(i) % torqueGen
torqueRotor => turbineArray(i) % torqueRotor
fluidDensity => turbineArray(i) % fluidDensity

GBRatio => turbineModel(j) % GBRatio
CutInGenSpeed => turbineModel(j) % CutInGenSpeed
CutInGenTorque => turbineModel(j) % CutInGenTorque
Region2StartGenSpeed => turbineModel(j) % Region2StartGenSpeed
KGen => turbineModel(j) % KGen
RatedGenSpeed => turbineModel(j) % RatedGenSpeed
Region2EndGenSpeed => turbineModel(j) % Region2EndGenSpeed
RatedGenTorque => turbineModel(j) % RatedGenTorque
RateLimitGenTorque => turbineModel(j) % RateLimitGenTorque
TorqueControllerRelax => turbineModel(j) % TorqueControllerRelax
DriveTrainIner => turbineModel(j) % DriveTrainIner

IntSpeedError => turbineArray(i) % IntSpeedError
PitchControlAngle => turbineArray(i) % PitchControlAngle

    ! No torque controller option
    if (turbineModel(j) % TorqueControllerType == "none") then

    elseif (turbineModel(j) % TorqueControllerType == "fiveRegion") then

        ! Get the generator speed.
        genSpeed = (rotSpeed/rpmRadSec)*GBRatio

        ! Save the generator torque from the last time step.
        torqueGenOld = torqueGen
            
        ! Region 1.
        if (genSpeed < CutInGenSpeed) then

            torqueGen = CutInGenTorque

        ! Region 1-1/2.
        elseif ((genSpeed >= CutInGenSpeed) .and.                              &
               (genSpeed < Region2StartGenSpeed)) then

        dGenSpeed = genSpeed - CutInGenSpeed
        Region2StartGenTorque = KGen * Region2StartGenSpeed *                  &
                                       Region2StartGenSpeed
        torqueSlope = (Region2StartGenTorque - CutInGenTorque) /               &
                      ( Region2StartGenSpeed - CutInGenSpeed )
        torqueGen = CutInGenTorque + torqueSlope*dGenSpeed

        ! Region 2.
        elseif ((genSpeed >= Region2StartGenSpeed) .and.                       &
                 (genSpeed < Region2EndGenSpeed)) then

                torqueGen = KGen * genSpeed * genSpeed

        ! Region 2-1/2.
        elseif ((genSpeed >= Region2EndGenSpeed) .and.                         &
                 (genSpeed < RatedGenSpeed)) then

                dGenSpeed = genSpeed - Region2EndGenSpeed
                Region2EndGenTorque = KGen * Region2EndGenSpeed *              &
                                             Region2EndGenSpeed
                torqueSlope = (RatedGenTorque - Region2EndGenTorque) /         &
                              ( RatedGenSpeed - Region2EndGenSpeed )
                torqueGen = Region2EndGenTorque + torqueSlope*dGenSpeed

        ! Region 3.
        elseif (genSpeed >= RatedGenSpeed) then

                torqueGen = RatedGenTorque
        endif

        ! Limit the change in generator torque if after first time step
        ! (otherwise it slowly ramps up from its zero initialized value--we
        ! want it to instantly be at its desired value on the first time
        ! step, but smoothly vary from there).
        if ((abs((torqueGen - torqueGenOld)/dt) > RateLimitGenTorque) &
              .and. (pastFirstTimeStep)) then

            if (torqueGen > torqueGenOld) then

                torqueGen = torqueGenOld + (RateLimitGenTorque * dt);

            elseif (torqueGen <= torqueGenOld) then

                torqueGen = torqueGenOld - (RateLimitGenTorque * dt);
            endif
        endif

        ! Update the rotor speed.
        rotSpeed = rotSpeed + TorqueControllerRelax * (dt/DriveTrainIner) *       &
                              (torqueRotor*fluidDensity - GBRatio*torqueGen)

        if (turbineModel(j) % PitchControllerType == "none") then
            ! Limit the rotor speed to be positive and such that the generator 
            !does not turn faster than rated.
            rotSpeed = max(0.0,rotSpeed)
            rotSpeed = min(rotSpeed,(RatedGenSpeed*rpmRadSec)/GBRatio)
        endif

    ! Torque control for fixed tip speed ratio
    ! Note that this current method does NOT support Coning in the rotor
    elseif (turbineModel(j) % TorqueControllerType == "fixedTSR") then

        if (pastFirstTimeStep) then
            ! Integrate the velocity along all actuator points
            call atm_integrate_u(i)   
    
            ! Match the rotor speed to a given TSR
            rotSpeed = turbineArray(i) % u_infinity_mean *     &
                       turbineArray(i) % TSR / turbineModel(j) % tipRad

            ! Important to get rid of negative values
            rotSpeed = max(0.0,rotSpeed)
        endif
    endif

    ! Pitch controllers (If there's no pitch controller, then don't do anything)
    if (turbineModel(j) % PitchControllerType == "gainScheduledPI") then

       ! Get the generator speed.
       genSpeed = (rotSpeed/rpmRadSec)*GBRatio

       ! Calculate the gain
       GK =  1.0/(1.0 + PitchControlAngle/turbineModel(j) % PitchControlAngleK)

       ! Calculate the Proportional and Integral terms
       KP = GK*turbineModel(j) % PitchControlKP0
       KI = GK*turbineModel(j) % PitchControlKI0

       ! Get speed error (generator in rpm) and update integral
       ! Integral is saturated to not push the angle beyond its limits
       SpeedError = genSpeed - RatedGenSpeed
       !write(*,*) 'Speed Error is: ', speedError

       IntSpeedError = IntSpeedError + SpeedError*dt
       IntSpeedError = min( max(IntSpeedError,                                 &
                       turbineModel(j) % PitchControlAngleMin/KI),             &
                       turbineModel(j) % PitchControlAngleMax/KI)

       ! Apply PI controller and saturate
       PitchControlAngle = KP*SpeedError + KI*IntSpeedError
       PitchControlAngle = min( max( PitchControlAngle,                        &
                           turbinemodel(j) % PitchControlAngleMin),            &
                           turbineModel(j) % PitchControlAngleMax)

    endif

    ! Compute the change in blade position at new rotor speed.
    turbineArray(i) % deltaAzimuth = rotSpeed * dt

end subroutine atm_computeRotorSpeed

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_rotateBlades(i)
! This subroutine rotates the turbine blades 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
!real(rprec), intent(in) :: dt                            ! time step
integer :: j                                 ! Turbine type
integer :: m, n, q                           ! Counters tu be used in do loops
real(rprec) :: deltaAzimuth, deltaAzimuthI   ! Angle of rotation
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: rotSpeed
real(rprec), pointer :: uvShaft(:)
real(rprec), pointer :: azimuth

j=turbineArray(i) % turbineTypeID

! Variables which are used by pointers
rotorApex=> turbineArray(i) % rotorApex
rotSpeed=>turbineArray(i) % rotSpeed
uvShaft=>turbineArray(i) % uvShaft
azimuth=>turbineArray(i) % azimuth

! Angle of rotation
deltaAzimuth = turbineArray(i) % deltaAzimuth

! Check the rotation direction first and set the local delta azimuth
! variable accordingly.
if (turbineArray(i) % rotationDir == "cw") then
    deltaAzimuthI = deltaAzimuth
else if (turbineArray(i) % rotationDir == "ccw") then
    deltaAzimuthI = -deltaAzimuth
endif

! Loop through all the points and rotate them accordingly
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
            turbineArray(i) %   bladePoints(m,n,q,:)=rotatePoint(              &
            turbineArray(i) % bladePoints(m,n,q,:), rotorApex, uvShaft,        &
            deltaAzimuthI)
        enddo
    enddo
enddo

if (pastFirstTimeStep) then
    azimuth = azimuth + deltaAzimuth;
        if (azimuth .ge. 2.0 * pi) then
            azimuth =azimuth - 2.0 *pi;
        endif
endif

end subroutine atm_rotateBlades

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compute_cl_correction(i)
! This subroutine compute the correction for the Cl
! It needs to be run only once in the initialization
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i             ! Turbine number
integer :: j                         ! Turbine type
integer :: m, n, q, k                ! Counters tu be used in do loops
real(rprec) :: a,b,c                 ! Correction coefficients
real(rprec) :: chord                 ! The chord value at the tip
!~ real(rprec) :: chord_r               ! The chord value at the root
real(rprec) :: r                     ! The radial distance from the tip
!~ real(rprec) :: r_r                   ! The radial distance from the root
real(rprec) :: eps_opt               ! The optimal epsilon
real(rprec) :: eps_s                 ! The epsilon value in the simulation
real(rprec) :: dG                    ! The finite difference G
real(rprec) :: dz                    ! The difference in distance between points
real(rprec) :: up_s                  ! The velocity perturbation for simulation
real(rprec) :: up_o                  ! The velocity perturbation for optimal
real(rprec) :: dup                   ! The difference in velocity perturbation
real(rprec) :: dWt                   ! The second order correction for tip
real(rprec) :: dWr                    ! The second order correction for tip

! Terms for correction equation
!~ real(rprec) :: term0, term1_tip, term2_tip, term1_root, term2_root 
!~ real(rprec) :: term1_tip_d, term2_tip_d, term1_root_d, term2_root_d 

! Constants for the tip vortex solution
a=0.029
b=-2./3.
c=0.357

! Simulation epsilon
eps_s = turbineArray(i) % epsilon

! The turbine type number
j=turbineArray(i) % turbineTypeID

! Firs compute the function G
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
            ! Compute G
!~             turbineArray(i) % G(m,n,q) = 1./2. * turbineArray(i) % Cl(m,n,q) * &
!~                 turbineArray(i) % chord(m,n,q) *                               &
!~                 turbineArray(i) % Vmag(m,n,q)**2
            turbineArray(i) % G(m,n,q) = 1./2. *  turbineArray(i) % Cl_b(m,n,q) * &
                turbineArray(i) % chord(m,n,q) *                               &
                turbineArray(i) % Vmag(m,n,q)**2

            turbineArray(i) % epsilon_opt(m,n,q) =                             &
                turbineArray(i) % chord(m,n,q) * turbineArray(i) % optimalEpsilonChord   

        enddo
    enddo
enddo

! Firs compute the function G
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
            ! Compute the difference in G using finite difference
            ! The tip and root are a special case
            if (q.eq.1) then
                turbineArray(i) % dG(m,n,q) = turbineArray(i) % G(m,n,1)
            elseif (q.eq.turbineArray(i) % numBladePoints) then
                turbineArray(i) % dG(m,n,q) =                                  &
                    - turbineArray(i) % G(m,n,turbineArray(i) % numBladePoints)
            else
                ! The central finite difference in G
                turbineArray(i) %  dG(m,n,q) = (turbineArray(i) % G(m,n,q+1) - &
                        turbineArray(i) % G(m,n,q-1))/2
            endif
!~     write(*,*) 'dG of q=', q, turbineArray(i) %  dG(m,n,q)

            ! First zero out the Cl correction
            turbineArray(i) % cl_correction(m, n, q) = 0.
        enddo
    enddo
enddo

! Apply the correction to the ALM points
! First compute the function G
do n=1, turbineArray(i) % numAnnulusSections
    do m=1, turbineModel(j) % numBl
        ! This is the first loop over all actuator points
        do q=1, turbineArray(i) % numBladePoints
            !!
            !!  Compute the contribution from the tip
            !!
            chord = turbineArray(i) % chord(m,n,                               &
                            turbineArray(i) % numBladePoints)

            ! Distance from the tip 
            r = abs(turbineArray(i) % bladeRadius(m,n,q) -                     &
                        turbineModel(j) % TipRad)

            ! The optimal epsilon value in meters
            eps_opt = turbineArray(i) % epsilon_opt(m,n,q)
            dWt =    1./2. *  a *                                              &
                    turbineArray(i) % Cl_b(m,n,turbineArray(i) % numBladePoints)*&
                        ! First term
                        ((eps_opt/chord)**(1.+b) * (chord / r)**2 *            &
                            (1. - exp(-c*abs(r/eps_opt)**3)) -                 &
                        ! Second term
                        (eps_s/chord)**(1.+b) * (chord / r)**2 *               &
                            (1. - exp(-c*abs(r/eps_s)**3)))

            ! Distance from the tip 
            chord = turbineArray(i) % chord(m,n,1)
            r = abs(turbineArray(i) % bladeRadius(m,n,q) -                     &
                        turbineModel(j) % HubRad)

            ! The difference in velocity from the root
            dWr =   1./2. *  a *                                               &
                    turbineArray(i) % Cl_b(m,n,1) *                              &
                        ! First term
                        ((eps_opt/chord)**(1.+b) * (chord / r)**2 *            &
                            (1. - exp(-c*abs(r/eps_opt)**3)) -                 &
                        ! Second term
                        (eps_s/chord)**(1.+b) * (chord / r)**2 *               &
                            (1. - exp(-c*abs(r/eps_s)**3)))

            turbineArray(i) % cl_correction(m, n, q) =                         &
                turbineArray(i) % cl_correction(m, n, q) + 2. * pi * (dWr + dWt)

            ! Now loop through all actuator points
            do k=1, turbineArray(i) % numBladePoints

                ! Compute dup only if not at the same actuator point
                if (k == q) then
                   dup = 0.
                else
                    ! Define the dG
                    dG = turbineArray(i) % dG(m,n,k)
                    ! Compute the difference in distance
                    dz = turbineArray(i) % bladeRadius(m,n,q) -                &
                    turbineArray(i) % bladeRadius(m,n,k)
                    ! Perturbation velocity for simulation epsilon
                    up_s = -dG / (4. * pi * dz) * (1. - exp(-(dz/eps_s)**2))
                    ! Perturbation velocity for optimal epsilon
                    up_o = -dG / (4. * pi * dz) * (1. - exp(-(dz/eps_opt)**2))
                    ! The difference in up
                    dup = up_o - up_s
!~                         write(*,*) 'dz',dz,'du', dup, up_o, up_s

                endif

                ! Additive correction
                turbineArray(i) % cl_correction(m, n, q) =                     &
                    turbineArray(i) % cl_correction(m, n, q) +                 &
                    2. * pi * dup / turbineArray(i) % Vmag(m,n,q)**2
            enddo
!~     write(*,*) 'dz',dz,'Cl Correction', q,' is:', turbineArray(i) % cl_correction(m, n, q)
!~     write(*,*) 'Epsilon', eps_s, eps_opt, chord
        enddo
    enddo
enddo

!~ ! Correction for the tip
!~ do q=1, turbineArray(i) % numBladePoints
!~     do n=1, turbineArray(i) % numAnnulusSections
!~         do m=1, turbineModel(j) % numBl

!~             ! The chord
 !!!           chord = turbineArray(i) % chord(m,n,q)
!~             chord = turbineArray(i) % chord(m,n,                               &
!~                             turbineArray(i) % numBladePoints)
!~             chord_r = turbineArray(i) % chord(m,n,1)

!~             ! Compute the optimal epsilon
!~             turbineArray(i) % epsilon_opt(m,n,q) =                             &
!~                 chord * turbineArray(i) % optimalEpsilonChord

!~             ! The optimal epsilon value in meters
!~             eps_opt = turbineArray(i) % epsilon_opt(m,n,q)

!~             ! Distance from the tip 
!~             r = abs(turbineArray(i) % bladeRadius(m,n,q)                       &
!~                         -                                                      &
!~                         turbineModel(j) % TipRad)
!~             r_r = abs(turbineArray(i) % bladeRadius(m,n,q)                     &
!~                         -                                                      &
!~                         turbineModel(j) % HubRad)

!~             ! The correction eta for both tip and root
!~             ! The correction needs to be multiplied by both radii^2
!~             term0 = 2. * (r/chord)**2 * (r_r/chord_r)**2

!~             ! Numerator terms
!~             term1_tip = (r_r/chord_r)**2 * (r/chord)/2 *                       &
!~                             (1. - exp(-r**2/eps_opt**2))
!~             term2_tip = (r_r/chord_r)**2 * 2. * pi * a *                       &
!~                             (eps_opt/chord)**(1.+b) *                          &
!~                             (1. - exp(-c*abs(r/eps_opt)**3))
!~             term1_root = (r/chord)**2 * (r_r/chord_r)/2 *                      &
!~                              (1. - exp(-r_r**2/eps_opt**2))
!~             term2_root = (r/chord)**2 * 2. * pi * a *                          &
!~                              (eps_opt/chord_r)**(1.+b) *                       &
!~                                  (1. - exp(-c*abs(r_r/eps_opt)**3))

!~             ! Denominator terms
!~             term1_tip_d = (r_r/chord_r)**2 * (r/chord)/2 *                     &
!~                               (1. - exp(-r**2/eps_s**2)) 
!~             term2_tip_d = (r_r/chord_r)**2 * 2. * pi * a *                     &
!~                               (eps_s/chord)**(1.+b) *                          &
!~                                 (1. - exp(-c*abs(r/eps_s)**3))
!~             term1_root_d = (r/chord)**2 * (r_r/chord_r)/2 *                    &
!~                                (1. - exp(-r_r**2/eps_s**2))
!~             term2_root_d = (r/chord)**2 * 2. * pi * a *                        &
!~                                (eps_s/chord_r)**(1.+b) *                       &
!~                                 (1. - exp(-c*abs(r_r/eps_s)**3))

!~             ! Correction for the root
!~             if (turbineArray(i) % rootALMCorrection .eqv. .true.)  then

!~                 turbineArray(i) % cl_correction(m, n, q) =                     &
!~                      (term0 -term1_tip   + term2_tip                           &
!~                      -term1_root   + term2_root )/                             &
!~                     (term0 -term1_tip_d + term2_tip_d                          & 
!~                     -term1_root_d + term2_root_d )
!~             else
!~                 ! The correction eta
!~                 turbineArray(i) % cl_correction(m, n, q) =                     &
!~                      (term0 -term1_tip   + term2_tip )/                        &
!~                     (term0 -term1_tip_d + term2_tip_d )

!~             endif
!~         enddo
!~     enddo
!~ enddo

end subroutine atm_compute_cl_correction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_calculate_variables(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Calculates the variables of the model that need information from the input
! files. It runs after reading input information.
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
integer :: m, n, q ! Looping indices
integer,     pointer :: NumSec  ! Number of sections in lookup table
real(rprec),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: projectionRadius, projectionRadiusNacelle
real(rprec), pointer :: sphereRadius
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: TipRad
real(rprec), pointer :: PreCone

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Number of sections in lookup table
NumSec => turbineModel(j) % NumSec

! Pointers dependent on turbineArray (i)
bladeRadius=>turbineArray(i) % bladeRadius
projectionRadius=>turbineArray(i) % projectionRadius
projectionRadiusNacelle=>turbineArray(i) % projectionRadiusNacelle
sphereRadius=>turbineArray(i) % sphereRadius

! Pointers dependent on turbineType (j)
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
TipRad=>turbineModel(j) % TipRad
PreCone=>turbineModel(j) %PreCone

! First compute the radius of the force projection (to the radius where the 
! projection is only 0.001 its maximum value - this seems to recover 99.9% of 
! the total forces when integrated
projectionRadius= turbineArray(i) % epsilon * sqrt(log(1.0/0.001))
projectionRadiusNacelle= turbineArray(i) % nacelleEpsilon*sqrt(log(1.0/0.001))

sphereRadius=sqrt(((OverHang + UndSling) + TipRad*sin(PreCone))**2 &
+ (TipRad*cos(PreCone))**2) + projectionRadius


! Compute the optimum value of epsilon for each blade section
do m=1, turbineModel(j) % numBl
    do n=1, turbineArray(i) %  numAnnulusSections
        do q=1, turbineArray(i) % numBladePoints

            ! Interpolate quantities through section
            turbineArray(i) % twistAng(m,n,q) =                                &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % twist(1:NumSec) )
          
            turbineArray(i) % chord(m,n,q) =                                   &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % chord(1:NumSec) )
            
            turbineArray(i) % sectionType(m,n,q) =                             &
                                   interpolate_i(bladeRadius(m,n,q),           &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % sectionType(1:NumSec))
        enddo
    enddo
enddo

!~ ! Compute the lift correction for this case
!~ call atm_compute_cl_correction(i)

end subroutine atm_calculate_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeBladeForce(i,m,n,q,U_local)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the wind vectors by projecting the velocity 
! onto the transformed coordinates system
implicit none

integer, intent(in) :: i,m,n,q
! i - turbineTypeArray
! n - numAnnulusSections
! q - numBladePoints
! m - numBl
real(rprec), intent(in) :: U_local(3)    ! The local velocity at this point

! Local variables
integer :: j,k ! Use to identify turbine type (j) and length of airoilTypes (k)
integer :: sectionType_i ! The type of airfoil
real(rprec) :: twistAng_i, chord_i, windAng_i, db_i, sigma, base_alpha
!real(rprec) :: solidity_i
real(rprec), dimension(3) :: dragVector, liftVector

! Pointers to be used
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: bladeAlignedVectors(:,:,:,:,:)
real(rprec), pointer :: windVectors(:,:,:,:)
real(rprec),  pointer :: bladePoints(:,:,:,:)
real(rprec), pointer :: rotSpeed
integer,     pointer :: numSec     
real(rprec),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: PreCone
real(rprec), pointer :: solidity(:,:,:),cl(:,:,:),cd(:,:,:),alpha(:,:,:)
real(rprec), pointer :: Vmag(:,:,:)

! Identifier for the turbine type
j= turbineArray(i) % turbineTypeID

! Pointers to trubineArray (i)
rotorApex => turbineArray(i) % rotorApex
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
windVectors => turbineArray(i) % windVectors
bladePoints => turbineArray(i) % bladePoints
rotSpeed => turbineArray(i) % rotSpeed
solidity=> turbineArray(i) % solidity
bladeRadius => turbineArray(i) % bladeRadius
cd => turbineArray(i) % cd       ! Drag coefficient
cl => turbineArray(i) % cl       ! Lift coefficient
alpha => turbineArray(i) % alpha ! Angle of attack
Vmag => turbineArray(i) % Vmag ! Velocity magnitude

!turbineTypeID => turbineArray(i) % turbineTypeID
NumSec => turbineModel(j) % NumSec

! Pointers for turbineModel (j)
PreCone => turbineModel(j) % PreCone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This will compute the vectors defining the local coordinate 
! system of the actuator point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define vector in z'
! If clockwise rotating, this vector points along the blade toward the tip.
! If counter-clockwise rotating, this vector points along the blade towards 
! the root.
if (turbineArray(i) % rotationDir == "cw")  then
    bladeAlignedVectors(m,n,q,3,:) =      &
                                     vector_add(bladePoints(m,n,q,:),-rotorApex)
elseif (turbineArray(i) % rotationDir == "ccw") then
    bladeAlignedVectors(m,n,q,3,:) =      &
                                     vector_add(-bladePoints(m,n,q,:),rotorApex)
endif

bladeAlignedVectors(m,n,q,3,:) =  &
                        vector_divide(bladeAlignedVectors(m,n,q,3,:),   &
                        vector_mag(bladeAlignedVectors(m,n,q,3,:)) )

! Define vector in y'
bladeAlignedVectors(m,n,q,2,:) = cross_product(bladeAlignedVectors(m,n,q,3,:), &
                                 turbineArray(i) % uvShaft)

bladeAlignedVectors(m,n,q,2,:) = vector_divide(bladeAlignedVectors(m,n,q,2,:), &
                                 vector_mag(bladeAlignedVectors(m,n,q,2,:)))

! Define vector in x'
bladeAlignedVectors(m,n,q,1,:) = cross_product(bladeAlignedVectors(m,n,q,2,:), &
                                 bladeAlignedVectors(m,n,q,3,:))

bladeAlignedVectors(m,n,q,1,:) = vector_divide(bladeAlignedVectors(m,n,q,1,:), &
                                 vector_mag(bladeAlignedVectors(m,n,q,1,:)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
! This concludes the definition of the local coordinate system


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now put the velocity in that cell into blade-oriented coordinates and add on 
! the velocity due to blade rotation.
windVectors(m,n,q,1) = dot_product(bladeAlignedVectors(m,n,q,1,:) , U_local)
windVectors(m,n,q,2) = dot_product(bladeAlignedVectors(m,n,q,2,:), U_local) + &
                      (rotSpeed * bladeRadius(m,n,q) * cos(PreCone))
windVectors(m,n,q,3) = dot_product(bladeAlignedVectors(m,n,q,3,:), U_local)

! Interpolated quantities through section
twistAng_i = turbineArray(i) % twistAng(m,n,q)
chord_i = turbineArray(i) % chord(m,n,q)
sectionType_i = turbineArray(i) % sectionType(m,n,q)

! Velocity magnitude
Vmag(m,n,q)=sqrt( windVectors(m,n,q,1)**2+windVectors(m,n,q,2)**2 )

! Angle between wind vector components
windAng_i = atan2( windVectors(m,n,q,1), windVectors(m,n,q,2) ) /degRad

! Local angle of attack
alpha(m,n,q) = windAng_i - twistAng_i - turbineArray(i) % Pitch - &
               turbineArray(i) % PitchControlAngle

! Total number of entries in lists of AOA, cl and cd
k = turbineModel(j) % airfoilType(sectionType_i) % n

! Lift coefficient
cl(m,n,q)= interpolate(alpha(m,n,q),                                           &
                 turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k),      &
                 turbineModel(j) % airfoilType(sectionType_i) % cl(1:k) )

base_alpha = 90. + alpha(m,n,q) - windAng_i
! Lift coefficient base (zero velocity in tangential direction)
turbineArray(i) % cl_b(m,n,q)= interpolate(base_alpha,                         &
                 turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k),      &
                 turbineModel(j) % airfoilType(sectionType_i) % cl(1:k) )

! Correct the lift coefficient
if (turbineArray(i) % tipALMCorrection .eqv. .true.)  then
    cl(m,n,q) = cl(m,n,q) + turbineArray(i) % cl_correction(m, n, q)
!~     cl(m,n,q) = cl(m,n,q) * turbineArray(i) % cl_correction(m, n, q)
endif

! Drag coefficient
cd(m,n,q)= interpolate(alpha(m,n,q),                                           &
                 turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k),      &
                turbineModel(j) % airfoilType(sectionType_i) % cd(1:k) )

db_i = turbineArray(i) % db(q) 

! Lift force
turbineArray(i) % lift(m,n,q) = 0.5_rprec * cl(m,n,q) * (Vmag(m,n,q)**2) *     &
                                chord_i * db_i * solidity(m,n,q)

! Drag force
turbineArray(i) % drag(m,n,q) = 0.5_rprec * cd(m,n,q) * (Vmag(m,n,q)**2) *     &
                                chord_i * db_i * solidity(m,n,q)

! This vector projects the drag onto the local coordinate system
dragVector = bladeAlignedVectors(m,n,q,1,:)*windVectors(m,n,q,1) +             &
             bladeAlignedVectors(m,n,q,2,:)*windVectors(m,n,q,2)

dragVector = vector_divide(dragVector,vector_mag(dragVector) )

! Lift vector
liftVector = cross_product(dragVector,bladeAlignedVectors(m,n,q,3,:) )
liftVector = liftVector/vector_mag(liftVector)

! Apply the lift and drag as vectors
liftVector = -turbineArray(i) % lift(m,n,q) * liftVector;
dragVector = -turbineArray(i) % drag(m,n,q) * dragVector;

! The blade force is the total lift and drag vectors 
turbineArray(i) % bladeForces(m,n,q,:) = vector_add(liftVector, dragVector)

! Find the component of the blade element force/density in the axial 
! (along the shaft) direction.
turbineArray(i) % axialForce(m,n,q) = dot_product(                           &
        -turbineArray(i) % bladeForces(m,n,q,:), turbineArray(i) % uvShaft)

! Find the component of the blade element force/density in the tangential 
! (torque-creating) direction.
turbineArray(i) % tangentialForce(m,n,q) = dot_product(                      &
       turbineArray(i) % bladeForces(m,n,q,:), bladeAlignedVectors(m,n,q,2,:))

! Change this back to radians
windAng_i = windAng_i * degRad
! The solidity
sigma = chord_i * turbineModel(j) % NumBl/ (2.*pi * bladeRadius(m,n,q) )

! Calculate the induction factor
turbineArray(i) % induction_a(m,n,q) = 1. / ( 4. * sin(windAng_i)**2 /  &
                (sigma * ( Cl(m,n,q) * cos(windAng_i) +   &
                Cd(m,n,q) * sin(windAng_i))) + 1.)

!~ write(*,*) 'Induction ', turbineArray(i) % induction_a(m,n,q)
! Calculate u infinity
turbineArray(i) % u_infinity(m,n,q) = windVectors(m,n,q,1) !/    &
!~                              (1. - turbineArray(i) % induction_a(m,n,q))
!~ turbineArray(i) % u_infinity(m,n,q) = Vmag(m,n,q) * sin(windAng_i) / &
!~ (1-turbineArray(i) % induction_a(m,n,q))


!~             turbineArray(i) % u_infinity = turbineArray(i) % u_infinity  +     &
!~                              windVectors(m,n,q,1) /    &
!~                              (1. - turbineArray(i) % induction_a(m,n,q))
        
! Calculate output quantities based on each point
call atm_process_output(i,m,n,q)

end subroutine atm_computeBladeForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeNacelleForce(i,U_local)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the force from the nacelle
implicit none

integer, intent(in) :: i ! i - turbineTypeArray
real(rprec), intent(in), dimension(3) :: U_local ! Velocity input

integer :: j ! j - turbineModel
real(rprec) :: V ! Velocity projected 
real(rprec), dimension(3) :: nacelleAlignedVector ! Nacelle vector
real(rprec) :: area, drag

! Identifier for the turbine type
j= turbineArray(i) % turbineTypeID

area = pi * turbineModel(j) % hubRad **2 

nacelleAlignedVector = turbineArray(i) % uvShaft

! Velocity projected in the direction of the nacelle
V = dot_product( nacelleAlignedVector , U_local)
!~     write(*,*) 'Nacelle Velocity before correction is: ', V
!~     write(*,*) 'nacelleAlignedVector is: ', nacelleAlignedVector

! The sampled velocity (uncorrected)
turbineArray(i) % VelNacelle_sampled = V

! Apply the velocity correction
V = V / (1. - .25/ pi * turbineArray(i) % nacelleCd * area /                   &
                turbineArray(i) % nacelleEpsilon**2 )
!~ write(*,*) 'Nacelle Velocity after correction is: ', V

! The velocity (corrected)
turbineArray(i) % VelNacelle_corrected = V

if (V .ge. 0.) then
    ! Drag force
    drag = 0.5_rprec * turbineArray(i) % nacelleCd * (V**2.) * area

    ! Drag Vector
    turbineArray(i) % nacelleForce = - drag * nacelleAlignedVector
!~     write(*,*) 'Nacelle Cd= ', turbineArray(i) % nacelleCd
!~     write(*,*) 'Nacelle Force is: ', turbineArray(i) % nacelleForce
endif

end subroutine atm_computeNacelleForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_integrate_u(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the induction factor a
! for each actuator point
implicit none

integer, intent(in) :: i ! i - turbineTypeArray
real(rprec), pointer, dimension(:,:,:) :: u_infinity, induction_a
real(rprec), pointer :: u_infinity_mean

induction_a => turbineArray(i) % induction_a
u_infinity_mean => turbineArray(i) % u_infinity_mean
u_infinity => turbineArray(i) % u_infinity

u_infinity_mean = sum(u_infinity) / size(u_infinity) / &
(1. - sum(induction_a) / size(induction_a))

!~ write(*,*) "U infinity is", u_infinity_mean, size(u_infinity)

end subroutine atm_integrate_u


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_yawNacelle(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine yaws the nacelle according to the yaw angle
implicit none

integer, intent(in) :: i
integer :: j 
integer :: m,n,q
! Perform rotation for the turbine.
j=turbineArray(i) % turbineTypeID 

! Rotate the rotor apex first.
turbineArray(i) % rotorApex = rotatePoint(turbineArray(i) % rotorApex,         &
                              turbineArray(i) % towerShaftIntersect,           &
                              turbineArray(i) % uvTower,                       &
                              turbineArray(i) % deltaNacYaw)

! Recompute the shaft unit vector since the shaft has rotated.
turbineArray(i) % uvShaft =                                                    &
                             turbineArray(i) % rotorApex -                     &
                             turbineArray(i) % towerShaftIntersect

turbineArray(i) % uvShaft = vector_divide(turbineArray(i) % uvShaft,           &
                            vector_mag(turbineArray(i) % uvShaft)) 
turbineArray(i) % uvShaft = vector_multiply( turbineArray(i) % uvShaft,        &
                                             turbineArray(i) % uvShaftDir )

! Rotate turbine blades, blade by blade, point by point.
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) %  numAnnulusSections
        do m=1, turbineModel(j) % numBl
            turbineArray(i) % bladePoints(m,n,q,:) =                           &
            rotatePoint( turbineArray(i) % bladePoints(m,n,q,:),               &
            turbineArray(i) % towerShaftIntersect,                             &
            turbineArray(i) % uvTower,                                         &
            turbineArray(i) % deltaNacYaw )                                    
        enddo
    enddo
enddo

! Compute the new yaw angle and make sure it isn't bigger than 2*pi.
if (pastFirstTimeStep) then
    turbineArray(i) % nacYaw = turbineArray(i) % nacYaw +                      &
                               turbineArray(i) % deltaNacYaw
    if (turbineArray(i) % nacYaw .ge. 2.0 * pi) then
        turbineArray(i) % nacYaw = turbineArray(i) % nacYaw - 2.0 * pi
    endif
endif

end subroutine atm_yawNacelle

!~ !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!~ subroutine atm_compassToStandard(dir)
!~ ! This function converts nacelle yaw from compass directions to the standard
!~ ! convention of 0 degrees on the + x axis with positive degrees
!~ ! in the counter-clockwise direction.
!~ !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!~ implicit none
!~ real(rprec), intent(inout) :: dir
!~ dir = dir + 180.0
!~ if (dir .ge. 360.0) then
!~     dir = dir - 360.0
!~ endif
!~ dir = 90.0 - dir
!~ if (dir < 0.0) then
!~     dir = dir + 360.0
!~ endif
!~  
!~ end subroutine atm_compassToStandard

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_output(i, jt_total, time)
! This subroutine will calculate and write the output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: jt_total ! Number of iteration fed in from solver
real(rprec), intent(in) :: time  ! time from simulation
integer, intent(in) :: i  ! The turbine number
integer :: j, m
integer :: powerFile=11, rotSpeedFile=12, bladeFile=13, liftFile=14, dragFile=15
integer :: ClFile=16, CdFile=17, alphaFile=18, VrelFile=19
integer :: VaxialFile=20, VtangentialFile=21, pitchFile=22, thrustFile=23
integer :: tangentialForceFile=24, axialForceFile=25, yawfile=26, nacelleFile=27

! Output only if the number of intervals is right
if ( mod(jt_total-1, outputInterval) == 0) then

    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!  Writing Actuator Turbine Model output  !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    
    j=turbineArray(i) % turbineTypeID ! The turbine type ID

    ! File for power output
    open(unit=powerFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/power")

    ! File for thrust
    open(unit=thrustFile,position="append",                                    &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/thrust")

    ! File for rotor speed
    open(unit=RotSpeedFile,position="append",                                  &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/RotSpeed")

    ! File for yaw
    open(unit=YawFile,position="append",                                  &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Yaw")

    ! File for blade output
    open(unit=bladeFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/blade")

    open(unit=liftFile,position="append",                                      &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/lift")

    open(unit=dragFile,position="append",                                      &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/drag")
    
    open(unit=ClFile,position="append",                                        &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cl")

    open(unit=CdFile,position="append",                                        &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cd")

    open(unit=alphaFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/alpha")

    open(unit=VrelFile,position="append",                                      &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vrel")

    open(unit=VaxialFile,position="append",                                    &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vaxial")

    open(unit=VtangentialFile,position="append",                               &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vtangential")

    open(unit=tangentialForceFile,position="append",                           &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/tangentialForce")

    open(unit=axialForceFile,position="append",                                &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/axialForce")

    open(unit=pitchFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/pitch")

    open(unit=nacelleFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/nacelle")

    call atm_compute_power(i)
    write(powerFile,*) time, turbineArray(i) % powerRotor,                     &
                       turbineArray(i) % powerGen
    write(thrustFile,*) time, turbineArray(i) % thrust
    write(RotSpeedFile,*) time, turbineArray(i) % RotSpeed
    write(pitchFile,*) time, turbineArray(i) % PitchControlAngle,              &
                       turbineArray(i) % IntSpeedError
    write(YawFile,*) time, turbineArray(i) % deltaNacYaw,                      &
                           turbineArray(i) % NacYaw
    write(nacelleFile,*) time, turbineArray(i) % VelNacelle_sampled,           &
                           turbineArray(i) % VelNacelle_corrected

    ! Will write only the first actuator section of the blade
    do m=1, turbineModel(j) % numBl
        write(bladeFile,*) i, m, turbineArray(i) % bladeRadius(m,1,:)
        write(liftFile,*) i, m, turbineArray(i) % lift(m,1,:)/                 &
                                turbineArray(i) % db(:)
        write(dragFile,*) i, m, turbineArray(i) % drag(m,1,:)/                 &
                                turbineArray(i) % db(:)
        write(ClFile,*) i, m, turbineArray(i) % cl(m,1,:)
        write(CdFile,*) i, m, turbineArray(i) % cd(m,1,:)
        write(alphaFile,*) i, m, turbineArray(i) % alpha(m,1,:)
        write(VrelFile,*) i, m, turbineArray(i) % Vmag(m,1,:)
        write(VaxialFile,*) i, m, turbineArray(i) % windVectors(m,1,:,1)
        write(VtangentialFile,*) i, m, turbineArray(i) %                       &
                                       windVectors(m,1,:,2)
        write(tangentialForceFile,*) i, m, turbineArray(i) %                   &
                                            tangentialForce(m,1,:)
        write(axialForceFile,*) i, m, turbineArray(i) % axialForce(m,1,:)

    enddo
    
        ! Write blade points 
!~         call atm_write_blade_points(i,jt_total)

    ! Close all the files 
    close(powerFile)
    close(thrustFile)
    close(rotSpeedFile)
    close(bladeFile)
    close(liftFile)
    close(dragFile)
    close(ClFile)
    close(CdFile)
    close(alphaFile)
    close(VrelFile)
    close(VaxialFile)
    close(VtangentialFile)
    close(pitchFile)
    close(tangentialForceFile)
    close(axialForceFile)
    close(yawFile)
    close(nacelleFile)

    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!  Done Writing Actuator Turbine Model output  !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

endif

end subroutine atm_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compute_power(i)
! This subroutine will calculate the total power of the turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j

j = turbineArray(i) % turbineTypeID

turbineArray(i) % powerRotor = turbineArray(i) % torqueRotor *                 &
    turbineArray(i) % rotSpeed * turbineArray(i) % fluidDensity
if (turbineModel(j) % TorqueControllerType == "fiveRegion") then
    turbineArray(i) % powerGen = turbineArray(i) % torqueGen *                 &
        turbineArray(i) % rotSpeed * turbineModel(j) % GBRatio
else
    turbineArray(i) % powerGen = turbineArray(i) % powerRotor
endif

write(*,*) 'Turbine ',i,' (Aerodynamic, Generator) Power is: ',                &
    turbineArray(i) % powerRotor, turbineArray(i) % powerGen

end subroutine atm_compute_power

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_write_blade_points(i,time_counter)
! This subroutine writes the position of all the blades at each time step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i, time_counter
integer :: m, n, q, j

j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=231, file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//  &
               '/blades'//trim(int2str(time_counter))//".vtk")

! Write the points to the blade file
do m=1, turbineModel(j) % numBl

    do n=1, turbineArray(i) %  numAnnulusSections

        do q=1, turbineArray(i) % numBladePoints

            write(231,*) turbineArray(i) % bladePoints(m,n,q,:)

        enddo

    enddo

enddo

close(231)

end subroutine atm_write_blade_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_process_output(i,m,n,q)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will process the output for an individual actuator point
! Important quantities such as power, thrust, ect are calculated here
implicit none

integer, intent(in) :: i,m,n,q
integer :: j
! Identify the turbine type
j=turbineArray(i) % turbineTypeID

turbineArray(i) % axialForce(m,n,q) = dot_product(                             &
-turbineArray(i) % bladeForces(m,n,q,:), turbineArray(i) % uvShaft)

! Find the component of the blade element force/density in the 
! tangential (torque-creating) direction.
turbineArray(i) % tangentialForce(m,n,q) = dot_product(                        &
                  turbineArray(i) % bladeForces(m,n,q,:) ,                     &
                  turbineArray(i) % bladeAlignedVectors(m,n,q,2,:))

! Add this blade element's contribution to thrust to the total 
! turbine thrust.
turbineArray(i) % thrust = turbineArray(i) % thrust +                          &
                           turbineArray(i) % axialForce(m,n,q)

! Add this blade element's contribution to aerodynamic torque to 
! the total turbine aerodynamic torque.
turbineArray(i) % torqueRotor = turbineArray(i) % torqueRotor +                &
                                turbineArray(i) % tangentialForce(m,n,q) *     &
                                turbineArray(i) % bladeRadius(m,n,q) *         &
                                cos(turbineModel(j) % PreCone)


end subroutine

!~ !-------------------------------------------------------------------------------
!~ function atm_convoluteForce(i,m,n,q,xyz)
!~ !-------------------------------------------------------------------------------
!~ ! This subroutine will convolute the body forces onto a point xyz
!~ integer, intent(in) :: i,m,n,q
!~ ! i - turbineTypeArray
!~ ! n - numAnnulusSections
!~ ! q - numBladePoints
!~ ! m - numBl
!~ real(rprec), intent(in) :: xyz(3)    ! Point onto which to convloute the force 
!~ real(rprec) :: Force(3)   ! The blade force to be convoluted
!~ real(rprec) :: dis                ! Distance onto which convolute the force
!~ real(rprec) :: atm_convoluteForce(3)    ! The local velocity at this point
!~ real(rprec) :: kernel                ! Gaussian dsitribution value
!~ 
!~ ! Distance from the point of the force to the point where it is being convoluted
!~ dis=distance(xyz,turbineArray(i) % bladepoints(m,n,q,:))
!~ 
!~ ! The force which is being convoluted
!~ Force=turbineArray(i) % bladeForces(m,n,q,:)
!~ 
!~ ! The value of the kernel. This is the actual smoothing function
!~ kernel=exp(-(dis/turbineArray(i) % epsilon)**2._rprec) /                       &
!~ ((turbineArray(i) % epsilon**3._rprec)*(pi**1.5_rprec))
!~ 
!~ ! The force times the kernel will give the force/unitVolume
!~ atm_convoluteForce = Force * kernel
!~ 
!~ return
!~ end function atm_convoluteForce



end module actuator_turbine_model

















