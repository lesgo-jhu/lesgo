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
public :: atm_initialize, numberOfTurbines,                                 &
          atm_computeBladeForce, atm_update,                                &
          vector_add, vector_divide, vector_mag, distance,                  &
          atm_convoluteForce, atm_output, atm_process_output,               &
          atm_initialize_output, atm_computeNacelleForce

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

call read_input_conf()  ! Read input data

do i = 1,numberOfTurbines
    inquire(file = "./turbineOutput/turbine"//trim(int2str(i))//              &
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

open(unit=1,         &
     file="./turbineOutput/turbine"//trim(int2str(i))//"/actuatorPoints", action='read')

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
subroutine atm_initialize_output()
! This subroutine initializes the output files for the ATM
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
logical :: file_exists
integer :: i

! Write to the screen output start
call atm_print_initialize()

do i = 1,numberOfTurbines

    inquire(file='./turbineOutput/turbine'//trim(int2str(i)),EXIST=file_exists)

    if (file_exists .eqv. .false.) then

        ! Create turbineOutput directory
        call system("mkdir -vp turbineOutput/turbine"//trim(int2str(i))) 

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/power") 
        write(1,*) 'turbineNumber Power'
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/lift")
        write(1,*) 'turbineNumber bladeNumber '
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/drag")
        write(1,*) 'turbineNumber bladeNumber '
        close(1)
        
        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/Cl")
        write(1,*) 'turbineNumber bladeNumber Cl'
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/Cd")
        write(1,*) 'turbineNumber bladeNumber Cd'
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/alpha")
        write(1,*) 'turbineNumber bladeNumber alpha'
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/Vrel")
        write(1,*) 'turbineNumber bladeNumber Vrel'
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/Vaxial")
        write(1,*) 'turbineNumber bladeNumber Vaxial'
        close(1)

        open(unit=1, file="./turbineOutput/turbine"//trim(int2str(i))//"/Vtangential")
        write(1,*) 'turbineNumber bladeNumber Vtangential'
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
call atm_compassToStandard(nacYaw)

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
    dist = HubRad
    
    ! Number of blade points for the first annular section
    do m=1, numBladePoints
        dist = dist + 0.5*(db(m))
        bladePoints(k,1,m,1) = root(1) + dist*sin(beta)
        bladePoints(k,1,m,2) = root(2)
        bladePoints(k,1,m,3) = root(3) + dist*cos(beta)
        do n=1,numAnnulusSections
            bladeRadius(k,n,m) = dist
            solidity(k,n,m)=1./numAnnulusSections
        enddo
        dist = dist + 0.5*db(m)
    enddo

    ! If there are more than one blade create the points of other blades by
    ! rotating the points of the first blade
    if (k > 1) then
        do m=1, numBladePoints
            bladePoints(k,1,m,:)=rotatePoint(bladePoints(k,1,m,:), rotorApex, &
            uvShaft,(360.0/NumBl)*(k-1)*degRad)
        enddo
    endif

    ! Rotate points for all the annular sections
    if (numAnnulusSections .lt. 2) cycle ! Cycle if only one section (ALM)
    do n=2, numAnnulusSections
        do m=1, numBladePoints
            bladePoints(k,n,m,:) =                                       &
            rotatePoint(bladePoints(k,1,m,:), rotorApex,                 &
            uvShaft,(annulusSectionAngle/(numAnnulusSections))*(n-1.)*degRad)
        enddo
    enddo
enddo

end subroutine atm_create_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_update(dt)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i                                 ! Turbine number
real(rprec), intent(in) :: dt                ! Time step

! Loop through all turbines and rotate the blades
do i = 1, numberOfTurbines
    call atm_rotateBlades(i,dt)
end do

end subroutine atm_update

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_rotateBlades(i,dt)
! This subroutine rotates the turbine blades 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
real(rprec), intent(in) :: dt                            ! time step
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
deltaAzimuth = rotSpeed * dt

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
subroutine atm_calculate_variables(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Calculates the variables of the model that need information from the input
! files. It runs after reading input information.
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
real(rprec), pointer :: projectionRadius, projectionRadiusNacelle
real(rprec), pointer :: sphereRadius
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: TipRad
real(rprec), pointer :: PreCone

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Pointers dependent on turbineArray (i)
projectionRadius=>turbineArray(i) % projectionRadius
projectionRadiusNacelle=>turbineArray(i) % projectionRadiusNacelle
sphereRadius=>turbineArray(i) % sphereRadius

! Pointers dependent on turbineType (j)
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
TipRad=>turbineModel(j) % TipRad
PreCone=>turbineModel(j) %PreCone

! First compute the radius of the force projection (to the radius where the 
! projection is only 0.0001 its maximum value - this seems to recover 99.99% of 
! the total forces when integrated
projectionRadius= turbineArray(i) % epsilon * sqrt(log(1.0/0.0001))
projectionRadiusNacelle= turbineArray(i) % nacelleEpsilon*sqrt(log(1.0/0.0001))

sphereRadius=sqrt(((OverHang + UndSling) + TipRad*sin(PreCone))**2 &
+ (TipRad*cos(PreCone))**2) + projectionRadius

! Calculate the smearing value epsilon for the nacelle
! It will be the hub radius unless the grid is coarser, and thus
! the epsilon value will be used
if (turbineArray(i) % nacelle) then
    if ( turbineArray(i) % epsilon > 2 * turbineModel(j) % hubRad ) then
        turbineArray(i) % nacelleEpsilon = turbineArray(i) % epsilon
    else
        turbineArray(i) % nacelleEpsilon = 2 * turbineModel(j) % hubRad
    endif
endif

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
real(rprec) :: twistAng_i, chord_i, windAng_i, db_i
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

! Interpolate quantities through section
twistAng_i = interpolate(bladeRadius(m,n,q),                                   &
                       turbineModel(j) % radius(1:NumSec),   &
                       turbineModel(j) % twist(1:NumSec) )

chord_i = interpolate(bladeRadius(m,n,q),                                      &
                       turbineModel(j) % radius(1:NumSec),   &
                       turbineModel(j) % chord(1:NumSec) )

sectionType_i = interpolate_i(bladeRadius(m,n,q),                              &
                       turbineModel(j) % radius(1:NumSec),   &
                       turbineModel(j) % sectionType(1:NumSec))

! Velocity magnitude
Vmag(m,n,q)=sqrt( windVectors(m,n,q,1)**2+windVectors(m,n,q,2)**2 )

! Angle between wind vector components
windAng_i = atan2( windVectors(m,n,q,1), windVectors(m,n,q,2) ) /degRad

! Local angle of attack
alpha(m,n,q) = windAng_i - twistAng_i - turbineArray(i) % Pitch

! Total number of entries in lists of AOA, cl and cd
k = turbineModel(j) % airfoilType(sectionType_i) % n

! Lift coefficient
cl(m,n,q)= interpolate(alpha(m,n,q),                                           &
                 turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k),      &
                 turbineModel(j) % airfoilType(sectionType_i) % cl(1:k) )

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

if (V > 0.) then
    ! Drag force
    drag = 0.5_rprec * turbineArray(i) % nacelleCd * (V**2) * area
    
    ! Drag Vector
    turbineArray(i) % nacelleForce = - drag * nacelleAlignedVector

write(*,*) 'Nacelle V is: ', V
endif

end subroutine atm_computeNacelleForce

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
                              turbineArray(i) % deltaNacYaw);

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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compassToStandard(dir)
! This function converts nacelle yaw from compass directions to the standard
! convention of 0 degrees on the + x axis with positive degrees
! in the counter-clockwise direction.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
real(rprec), intent(inout) :: dir
dir = dir + 180.0
if (dir .ge. 360.0) then
    dir = dir - 360.0
endif
dir = 90.0 - dir
if (dir < 0.0) then
    dir = dir + 360.0
endif
 
end subroutine atm_compassToStandard

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_output(jt_total)
! This subroutine will calculate the output of the model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: jt_total ! Number of iteration fed in from solver
integer :: i, j, m, n, q
integer :: powerFile=11, bladeFile=12, liftFile=13, dragFile=14
integer :: ClFile=15, CdFile=16, alphaFile=17, VrelFile=18
integer :: VaxialFile=19, VtangentialFile=20
integer :: pointsFile=787 ! File to write the actuator points

! Output only if the number of intervals is right
if ( mod(jt_total-1, outputInterval) == 0) then
        
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!  Writing Actuator Turbine Model output  !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    do i=1,numberOfTurbines

        j=turbineArray(i) % turbineTypeID ! The turbine type ID

        ! Files for power output
        open(unit=powerFile,position="append",                                 &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/power")
    
        ! File for blade output
        open(unit=bladeFile,position="append",                                 &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/blade")

        open(unit=liftFile,position="append",                                  &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/lift")

        open(unit=dragFile,position="append",                                  &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/drag")
        
        open(unit=ClFile,position="append",                                    &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/Cl")

        open(unit=CdFile,position="append",                                    &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/Cd")

        open(unit=alphaFile,position="append",                                 &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/alpha")

        open(unit=VrelFile,position="append",                                  &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/Vrel")

        open(unit=VaxialFile,position="append",                                &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/Vaxial")

        open(unit=VtangentialFile,position="append",                           &
        file="./turbineOutput/turbine"//trim(int2str(i))//"/Vtangential")

        call atm_compute_power(i)
        write(powerFile,*) i, turbineArray(i) % powerRotor

        ! Will write only the first actuator section of the blade
        do m=1, turbineModel(j) % numBl
            write(bladeFile,*) i, m, turbineArray(i) % bladeRadius(m,1,:)
            write(liftFile,*) i, m, turbineArray(i) % lift(m,1,:)/   &
                                    turbineArray(i) % db(:)
            write(dragFile,*) i, m, turbineArray(i) % drag(m,1,:)/   &
                                    turbineArray(i) % db(:)
            write(ClFile,*) i, m, turbineArray(i) % cl(m,1,:)
            write(CdFile,*) i, m, turbineArray(i) % cd(m,1,:)
            write(alphaFile,*) i, m, turbineArray(i) % alpha(m,1,:)
            write(VrelFile,*) i, m, turbineArray(i) % Vmag(m,1,:)
            write(VaxialFile,*) i, m, turbineArray(i) % windVectors(m,1,:,1)
            write(VtangentialFile,*) i, m, turbineArray(i) %                   &
                                           windVectors(m,1,:,2)

        enddo
    
        ! Write blade points 
        call atm_write_blade_points(i,jt_total)

    enddo

    ! Close all the files 
    close(powerFile)
    close(bladeFile)
    close(liftFile)
    close(dragFile)
    close(ClFile)
    close(CdFile)
    close(alphaFile)
    close(VrelFile)
    close(VaxialFile)
    close(VtangentialFile)

    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!  Done Writing Actuator Turbine Model output  !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

endif

! Write the actuator points at every time-step regardless
do i=1,numberOfTurbines

    j=turbineArray(i) % turbineTypeID ! The turbine type ID

    open(unit=pointsFile,status="replace",                                     &
         file="./turbineOutput/turbine"//trim(int2str(i))//"/actuatorPoints")

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

enddo

end subroutine atm_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compute_power(i)
! This subroutine will calculate the total power of the turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i

turbineArray(i) % powerRotor = turbineArray(i) % torqueRotor *                 &
                               turbineArray(i) % rotSpeed

write(*,*) 'Turbine ',i,' Power is: ', turbineArray(i) % powerRotor

end subroutine atm_compute_power

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_write_blade_points(i,time_counter)
! This subroutine writes the position of all the blades at each time step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i, time_counter
integer :: m, n, q, j

j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=231, file="./turbineOutput/turbine"//trim(int2str(i))//'/blades'     &
                     //trim(int2str(time_counter))//".vtk")

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

!-------------------------------------------------------------------------------
function atm_convoluteForce(i,m,n,q,xyz)
!-------------------------------------------------------------------------------
! This subroutine will convolute the body forces onto a point xyz
integer, intent(in) :: i,m,n,q
! i - turbineTypeArray
! n - numAnnulusSections
! q - numBladePoints
! m - numBl
real(rprec), intent(in) :: xyz(3)    ! Point onto which to convloute the force 
real(rprec) :: Force(3)   ! The blade force to be convoluted
real(rprec) :: dis                ! Distance onto which convolute the force
real(rprec) :: atm_convoluteForce(3)    ! The local velocity at this point
real(rprec) :: kernel                ! Gaussian dsitribution value

! Distance from the point of the force to the point where it is being convoluted
dis=distance(xyz,turbineArray(i) % bladepoints(m,n,q,:))

! The force which is being convoluted
Force=turbineArray(i) % bladeForces(m,n,q,:)

! The value of the kernel. This is the actual smoothing function
kernel=exp(-(dis/turbineArray(i) % epsilon)**2._rprec) /                       &
((turbineArray(i) % epsilon**3._rprec)*(pi**1.5_rprec))

! The force times the kernel will give the force/unitVolume
atm_convoluteForce = Force * kernel

return
end function atm_convoluteForce



end module actuator_turbine_model

















