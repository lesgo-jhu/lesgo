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
public :: atm_initialize, numberOfTurbines,                    &
          atm_computeBladeForce, atm_update,                                &
          vector_add, vector_divide, vector_mag, distance, atm_convoluteForce

! These are used to do unit conversions
real(rprec), parameter :: pi=acos(-1.) !pi= 3.141592653589793238462643383279502884197169399375
real(rprec) :: degRad = pi/180. ! Degrees to radians conversion
real(rprec) :: rpmRadSec =  pi/30. ! Set the revolutions/min to radians/s 

!integer :: i, j, k ! Counters

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
pastFirstTimeStep=.false. ! The first time step not reached yet

call read_input_conf()  ! Read input data

do i = 1,numberOfTurbines
    call atm_create_points(i)   ! Creates the ATM points defining the geometry

    ! This will create the first yaw alignment
    turbineArray(i) % deltaNacYaw = turbineArray(i) % nacYaw
    call atm_yawNacelle(i)
write(*,*) 'turbineArray(i) % deltaNacYaw = ',turbineArray(i) % deltaNacYaw
!write(*,*) 'Points in blade ',i,' = ', turbineArray(i) % bladePoints

    call atm_calculate_variables(i) ! Calculates variables depending on input
end do
pastFirstTimeStep=.true. ! Past the first time step

end subroutine atm_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_create_points(i)
! This subroutine generate the set of blade points for each turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
type(real(rprec)),  pointer :: db(:)
type(real(rprec)),  pointer :: bladePoints(:,:,:,:)
type(real(rprec)),  pointer :: bladeRadius(:,:,:)
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
call compassToStandard(nacYaw)
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
!write(*,*) 'turbineArray(i) % baseLocation = ', turbineArray(i) % baseLocation
towerShaftIntersect(3) = towerShaftIntersect(3) + TowerHt + Twr2Shft
!write(*,*) 'towerShaftIntersect(3) = ', towerShaftIntersect(3)
rotorApex = towerShaftIntersect
!write(*,*) 'rotorApex = ', rotorApex

rotorApex(1) = rotorApex(1) +  (OverHang + UndSling) * cos(ShftTilt)
rotorApex(3) = rotorApex(3) +  (OverHang + UndSling) * sin(ShftTilt)
!write(*,*) 'rotorApex = ', rotorApex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Create the first set of actuator points                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the vector along the shaft pointing in the direction of the wind
uvShaftDir = OverHang / abs( OverHang )
!write(*,*) 'uvShaftDir = ', uvShaftDir

!write(*,*) 'rotorApex = ', rotorApex
!write(*,*) 'towerShaftIntersect = ', towerShaftIntersect

! Define the vector along the shaft pointing in the direction of the wind                               
uvShaft = vector_add(rotorApex , - towerShaftIntersect)
uvShaft = vector_divide(uvShaft, vector_mag(uvShaft))
uvShaft = vector_multiply( uvShaft,uvShaftDir)
!write(*,*) 'uvShaft = ', uvShaft

!write(*,*) 'uvShaft = ',uvShaft
! Define vector aligned with the tower pointing from the ground to the nacelle
uvTower = vector_add(towerShaftIntersect, - baseLocation)
uvTower = vector_divide( uvTower, vector_mag(uvTower))
!write(*,*) 'uvTower = ', uvTower

! Define thickness of each blade section
!write(*,*) 'TipRad = ',TipRad
!write(*,*)'HubRad = ',HubRad
!write(*,*)'numBladePoints = ', numBladePoints
do k=1, numBladePoints
    db(k) = (TipRad - HubRad)/(numBladePoints)
!write(*,*) 'db(k) = ', db(k)
enddo

! This creates the first set of points
do k=1, numBl
    root = rotorApex
    beta = PreCone - ShftTilt
    root(1)= root(1) + HubRad*sin(beta)
    root(3)= root(3) + HubRad*cos(beta)
    dist = HubRad
!write(*,*) 'beta, dis = ', beta, dist
!write(*,*) 'root = ',root
    ! Number of blade points for the first annular section
    do m=1, numBladePoints
        dist = dist + 0.5*(db(m))                                                  ! Bug in OpenFOAM Code
        bladePoints(k,1,m,1) = root(1) + dist*sin(beta)
        bladePoints(k,1,m,2) = root(2)
        bladePoints(k,1,m,3) = root(3) + dist*cos(beta)
        do n=1,numAnnulusSections
            bladeRadius(k,n,m) = dist
            solidity(k,n,m)=1./numAnnulusSections
        enddo
        dist = dist + 0.5*db(m)
!write(*,*) 'bladePoints(k,1,m,:) = ', bladePoints(k,1,m,:)
!write(*,*) 'bladeRadius(k,1,m) = ',  bladeRadius(k,1,m)
    enddo

    if (k > 1) then
        do m=1, numBladePoints
            bladePoints(k,1,m,:)=rotatePoint(bladePoints(k,1,m,:), rotorApex, &
            uvShaft,(360.0/NumBl)*(k-1)*degRad)
!write(*,*) 'bladePoints(k,1,m,:) = ', bladePoints(k,1,m,:)
!write(*,*) '(360.0/NumBl)*k*degRad = ',(360.0/NumBl)*(k-1)*degRad
!write(*,*) 'bladeRadius(k,1,m) = ',  bladeRadius(k,1,m)
        enddo
    endif

    ! Rotate points for all the annular sections
!    if (numAnnulusSections .lt. 2) cycle ! Cycle if only one annular section

    if (numAnnulusSections .lt. 2) cycle
    do n=2, numAnnulusSections
        do m=1, numBladePoints
            bladePoints(k,n,m,:) =                                       &
            rotatePoint(bladePoints(k,1,m,:), rotorApex,                 &
            uvShaft,(annulusSectionAngle/(numAnnulusSections))*(n-1)*degRad)
        enddo
    enddo
enddo

open(unit=787,file='./output/points')
do k=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
            write(787,*) bladePoints(m,n,k,:)/400.
        enddo
    enddo
enddo
close(787)

end subroutine atm_create_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_update(dt)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer :: i                                 ! Turbine number
real(rprec), intent(in) :: dt                ! Time step
!write(*,*) 'dt = ', dt

do i = 1, numberOfTurbines
    call atm_rotateBlades(i,dt)              ! Rotate the blades of each turbine
end do

end subroutine atm_update

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_rotateBlades(i,dt)
! This subroutine rotates the turbine blades 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!write(*,*) 'dt = ', dt
!write(*,*) 'turbineArray(i) % rotationDir = ' ,turbineArray(i) % rotationDir
!write(*,*) 'rotSpeed = ', rotSpeed
!write(*,*) 'deltaAzimuth = ', deltaAzimuth
!write(*,*) 'deltaAzimuthI = ', deltaAzimuthI

do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
!write(*,*) 'turbineArray(i) % bladePoints(m,n,q,:) before= ', turbineArray(i) % bladePoints(m,n,q,:)
!write(*,*) 'rotorApex = ',rotorApex
!write(*,*) ' uvShaft = ',uvShaft
!write(*,*) ' deltaAzimuthI = ',deltaAzimuthI
            turbineArray(i) %   bladePoints(m,n,q,:)=rotatePoint(              &
            turbineArray(i) % bladePoints(m,n,q,:), rotorApex, uvShaft,        &
            deltaAzimuthI)
!write(*,*) 'turbineArray(i) % bladePoints(m,n,q,:) after= ', turbineArray(i) % bladePoints(m,n,q,:)
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
integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
real(rprec), pointer :: projectionRadius
real(rprec), pointer :: sphereRadius
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: TipRad
real(rprec), pointer :: PreCone

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Pointers dependent on turbineArray (i)
projectionRadius=>turbineArray(i) % projectionRadius
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

sphereRadius=sqrt(((OverHang + UndSling) + TipRad*sin(PreCone))**2 &
+ (TipRad*cos(PreCone))**2) + projectionRadius

end subroutine atm_calculate_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeBladeForce(i,m,n,q,U_local)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the wind vectors by projecting the velocity 
! onto the transformed coordinates system
integer, intent(in) :: i,m,n,q
! i - turbineTypeArray
! n - numAnnulusSections
! q - numBladePoints
! m - numBl
real(rprec), intent(in) :: U_local(3)    ! The local velocity at this point

! Local variables
integer :: j,k ! Use to identify turbine type (j) and length of airoilTypes (k)
integer :: sectionType_i ! The type of airfoil
real(rprec) :: cl_i, cd_i, twistAng_i, chord_i, Vmag_i, windAng_i, alpha_i, db_i
!real(rprec) :: solidity_i
real(rprec), dimension(3) :: dragVector, liftVector

! Pointers to be used
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: bladeAlignedVectors(:,:,:,:,:)
real(rprec), pointer :: windVectors(:,:,:,:)
type(real(rprec)),  pointer :: bladePoints(:,:,:,:)
real(rprec), pointer :: rotSpeed
integer,     pointer :: numSec     
type(real(rprec)),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: PreCone
real(rprec), pointer :: solidity(:,:,:)

j= turbineArray(i) % turbineTypeID

! Pointers to trubineArray (i)
rotorApex => turbineArray(i) % rotorApex
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
windVectors => turbineArray(i) % windVectors
bladePoints => turbineArray(i) % bladePoints
rotSpeed => turbineArray(i) % rotSpeed
solidity=> turbineArray(i) % solidity
! Pointers for turbineModel (j)
!turbineTypeID => turbineArray(i) % turbineTypeID
NumSec => turbineModel(j) % NumSec
bladeRadius => turbineArray(i) % bladeRadius
!cd(m,n,q) => turbineArray(i) % cd(m,n,q)    ! Drag coefficient
!cl(m,n,q) => turbineArray(i) % cl(m,n,q)    ! Lift coefficient
!alpha(m,n,q) => turbineArray(i) % alpha(m,n,q) ! Anlge of attack

PreCone => turbineModel(j) % PreCone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This will compute the vectors defining the local coordinate 
! system of the actuator point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define vector in z'
! If clockwise rotating, this vector points along the blade toward the tip.
! If counter-clockwise rotating, this vector points along the blade towards 
! the root.

!write(*,*) ' turbineArray(i) % rotationDir = ', turbineArray(i) % rotationDir
if (turbineArray(i) % rotationDir == "cw")  then
    bladeAlignedVectors(m,n,q,3,:) =      &
                                     vector_add(bladePoints(m,n,q,:),-rotorApex)
elseif (turbineArray(i) % rotationDir == "ccw") then
    bladeAlignedVectors(m,n,q,3,:) =      &
                                     vector_add(-bladePoints(m,n,q,:),rotorApex)
endif

!write(*,*) 'bladePoints(m,n,q,:) = ', bladePoints(m,n,q,:)
!write(*,*)'rotorApex = ',rotorApex
!write(*,*) 'bladeAlignedVectors(m,n,q,3,:) 2= ' , bladeAlignedVectors(m,n,q,3,:)

bladeAlignedVectors(m,n,q,3,:) =  &
                        vector_divide(bladeAlignedVectors(m,n,q,3,:),   &
                        vector_mag(bladeAlignedVectors(m,n,q,3,:)) )

!write(*,*) 'bladeAlignedVectors(m,n,q,3,:) = ' , bladeAlignedVectors(m,n,q,3,:)

! Define vector in y'
bladeAlignedVectors(m,n,q,2,:) = cross_product(bladeAlignedVectors(m,n,q,3,:), &
                                 turbineArray(i) % uvShaft)
!write(*,*) 'cross_product(bladeAlignedVectors(m,n,q,3,:),  turbineArray(i) % uvShaft)',cross_product(bladeAlignedVectors(m,n,q,3,:), turbineArray(i) % uvShaft)
!write(*,*) 'size is ', size(bladeAlignedVectors)

!write(*,*) 'bladeAlignedVectors(m,n,q,3,:) = ', bladeAlignedVectors(m,n,q,3,:)
!write(*,*) 'uvShaft = ', turbineArray(i) % uvShaft
!write(*,*) 'bladeAlignedVectors(m,n,q,2,:) = ' , bladeAlignedVectors(m,n,q,2,:)

bladeAlignedVectors(m,n,q,2,:) = vector_divide(bladeAlignedVectors(m,n,q,2,:), &
                                 vector_mag(bladeAlignedVectors(m,n,q,2,:)))

!write(*,*) 'bladeAlignedVectors(m,n,q,2,:) = ' , bladeAlignedVectors(m,n,q,2,:)

! Define vector in x'
bladeAlignedVectors(m,n,q,1,:) = cross_product(bladeAlignedVectors(m,n,q,2,:), &
                                 bladeAlignedVectors(m,n,q,3,:))
!write(*,*) 'bladeAlignedVectors(m,n,q,1,:) = ', bladeAlignedVectors(m,n,q,1,:)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bladeAlignedVectors(m,n,q,1,:) = vector_divide(bladeAlignedVectors(m,n,q,1,:), &
                                 vector_mag(bladeAlignedVectors(m,n,q,1,:)))
!write(*,*) bladeAlignedVectors(m,n,q,1,:)
! This concludes the definition of the local corrdinate system


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now put the velocity in that cell into blade-oriented coordinates and add on 
! the velocity due to blade rotation.
windVectors(m,n,q,1) = dot_product(bladeAlignedVectors(m,n,q,1,:) , U_local)
windVectors(m,n,q,2) = dot_product(bladeAlignedVectors(m,n,q,2,:), U_local) + &
                      (rotSpeed * bladeRadius(m,n,q) * cos(PreCone))
windVectors(m,n,q,3) = dot_product(bladeAlignedVectors(m,n,q,3,:), U_local)
!write(*,*) 'rotSpeed ', rotSpeed
!write(*,*) 'bladeRadius(m,n,q) ',   bladeRadius(m,n,q)
!write(*,*) 'Precone is ',Precone
!write(*,*) 'windVectors(m,n,q,:) = ', windVectors(m,n,q,:)

! Interpolate quantities through section
twistAng_i = interpolate(bladeRadius(m,n,q),                                   &
                       turbineModel(j) % radius(1:NumSec),   &
                       turbineModel(j) % twist(1:NumSec) )   
chord_i = interpolate(bladeRadius(m,n,q),                                      &
                       turbineModel(j) % radius(1:NumSec),   &
                       turbineModel(j) % chord(1:NumSec) )

sectionType_i = interpolate_i(bladeRadius(m,n,q),                              &     ! Problem here with turbineModel(i) % sectionType(1:turbineModel(i)% NumSec
                       turbineModel(j) % radius(1:NumSec),   &
                       turbineModel(j) % sectionType(1:NumSec))
!write(*,*) 'bladeRadius(m,n,q) = ', bladeRadius(m,n,q)
!write(*,*) 'turbineModel(j) % radius(1:NumSec) = ', turbineModel(j) % radius(1:NumSec)
!write(*,*) 'turbineModel(j) % sectionType(1:NumSec) = ', turbineModel(j) % sectionType(1:NumSec)
!write(*,*) 'sectionType_i = ', sectionType_i
! Velocity magnitude
Vmag_i=sqrt( windVectors(m,n,q,1)**2+windVectors(m,n,q,2)**2 )

! Angle between wind vector components
windAng_i = atan2( windVectors(m,n,q,1), windVectors(m,n,q,2) ) /degRad

! Local angle of attack
alpha_i=windAng_i-twistAng_i - turbineArray(i) % Pitch

!write(*,*) 'Error YES Here'
!write(*,*) sectionType_i

! Total number of entries in lists of AOA, cl and cd
k = turbineModel(j) % airfoilType(sectionType_i) % n

!write(*,*) 'Error NOT Here'


! Lift coefficient
cl_i= interpolate(alpha_i,                                                     &
                 turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k),      &
                 turbineModel(j) % airfoilType(sectionType_i) % cl(1:k) )


!write(*,*) 'turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k) = ',turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k)
!write(*,*) 'turbineModel(j) % airfoilType(sectionType_i) % cl(1:k) = ',turbineModel(j) % airfoilType(sectionType_i) % cl(1:k)
!write(*,*) 'Cl = ',cl_i
! Drag coefficient
cd_i= interpolate(alpha_i,                                                  &
                 turbineModel(j) % airfoilType(sectionType_i) % AOA(1:k),      &
                turbineModel(j) % airfoilType(sectionType_i) % cd(1:k) )
!write(*,*) 'Cd = ',cd_i

db_i = turbineArray(i) % db(q) 

! Lift force
turbineArray(i) % lift(m,n,q) = 0.5 * cl_i * Vmag_i**2 * chord_i * db_i * solidity(m,n,q)
!write(*,*) 'turbineArray(i) % lift(m,n,q) = ',turbineArray(i) % lift(m,n,q)

! Drag force
turbineArray(i) % drag(m,n,q) = 0.5 * cd_i * Vmag_i**2 * chord_i * db_i * solidity(m,n,q)

dragVector = bladeAlignedVectors(m,n,q,1,:)*windVectors(m,n,q,1) +  &
             bladeAlignedVectors(m,n,q,2,:)*windVectors(m,n,q,2)

dragVector = vector_divide(dragVector,vector_mag(dragVector) )

! Lift vector
liftVector = cross_product(dragVector,bladeAlignedVectors(m,n,q,3,:) )
liftVector = vector_divide(liftVector,vector_mag(liftVector))

liftVector = -turbineArray(i) % lift(m,n,q) * liftVector;
dragVector = -turbineArray(i) % drag(m,n,q) * dragVector;

turbineArray(i) % bladeForces(m,n,q,:) = vector_add(liftVector, dragVector)
turbineArray(i) % bladeForcesDummy(m,n,q,:) = turbineArray(i) % bladeForces(m,n,q,:)
!write(*,*) 'turbineArray(i) % bladeForces(m,n,q,:) = ', turbineArray(i) % bladeForces(m,n,q,:)


end subroutine atm_computeBladeForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_yawNacelle(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine yaws the nacelle according to the yaw angle
integer, intent(in) :: i
integer :: j 
integer :: m,n,q
! Perform rotation for the turbine.
j=turbineArray(i) % turbineTypeID 

! Rotate the rotor apex first.
turbineArray(i) % rotorApex = rotatePoint(turbineArray(i) % rotorApex,       &
                              turbineArray(i) % towerShaftIntersect,         &
                              turbineArray(i) % uvTower,                     &
                              turbineArray(i) % deltaNacYaw);

! Recompute the shaft unit vector since the shaft has rotated.
turbineArray(i) % uvShaft =                                                &
                             turbineArray(i) % rotorApex -                 &
                             turbineArray(i) % towerShaftIntersect
	
turbineArray(i) % uvShaft = vector_divide(turbineArray(i) % uvShaft,       &
                            vector_mag(turbineArray(i) % uvShaft)) 
turbineArray(i) % uvShaft = vector_multiply( turbineArray(i) % uvShaft,    &
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
function atm_convoluteForce(i,m,n,q,xyz)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the wind vectors by projecting the velocity 
! onto the transformed coordinates system
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

dis=distance(xyz,turbineArray(i) % bladepoints(m,n,q,:))
Force=turbineArray(i) % bladeForces(m,n,q,:)
kernel=exp(-(dis/turbineArray(i) % epsilon)**2.) / &
((turbineArray(i) % epsilon**3.)*(pi**1.5))
atm_convoluteForce = Force * kernel

!if (dis .le. 10.) then
!    write(*,*) 'dis= ', dis, 'kernel = ', kernel
!    write(*,*) 'epsilon = ', turbineArray(i) % epsilon
!endif
return
end function atm_convoluteForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine compassToStandard(dir)
! This function converts nacelle yaw from compass directions to the standard
! convention of 0 degrees on the + x axis with positive degrees
! in the counter-clockwise direction.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec), intent(inout) :: dir
dir = dir + 180.0
if (dir .ge. 360.0) then
    dir = dir - 360.0
endif
dir = 90.0 - dir
if (dir < 0.0) then
    dir = dir + 360.0
endif
 
end subroutine compassToStandard

!-------------------------------------------------------------------------------
function interpolate(xp,x,y)
! This function interpolates xp from x and y 
!-------------------------------------------------------------------------------
real(rprec), dimension(:), intent(in) :: x,y
real(rprec), intent(in) ::  xp
real(rprec) :: xa, xb, ya, yb
integer :: i,p
real(rprec) :: interpolate
p=size(x)

    do i=2,p
        if ( ( xp .ge. x(i-1) ) .and. ( xp .le. x(i) ) ) then
            xa=x(i-1)
            xb=x(i)
            ya=y(i-1)
            yb=y(i)
            interpolate = ya + (yb-ya) * (xp-xa) / (xb-xa) 
        endif
    enddo

return
end function interpolate

!-------------------------------------------------------------------------------
integer function interpolate_i(xp,x,y)
! This function interpolates xp from x and y 
!-------------------------------------------------------------------------------
real(rprec), dimension(:), intent(in) :: x
integer, dimension(:), intent(in) :: y
real(rprec), intent(in) ::  xp
real(rprec) :: xa,xb,ya,yb
integer :: i,p
p=size(x)

if (xp .lt. x(1)) then
    interpolate_i=y(1) 
    else if (xp .gt. x(p)) then
        interpolate_i=y(p)
    else
        do i=2,p
            if ( ( xp .ge. x(i-1) ) .and. ( xp .le. x(i) ) ) then
                xa=x(i-1)
                xb=x(i)
                ya=real(y(i-1),rprec)
                yb=real(y(i),rprec)
                interpolate_i=nint( ya + (yb-ya) * (xp-xa) / (xb-xa) )
!write(*,*) 'Interpolation = ', nint( ya + (yb-ya) * (xp-xa) / (yb-ya) )
                endif
        enddo
endif
!write(*,*) 'Value of y = ', ya, yb
!write(*,*) 'Value of x = ', xa, xb
!write(*,*) 'Value of xp = ', xp
!write(*,*) 'Value of p = ', p
!write(*,*) 'Interpolated Value = ', interpolate_i
return
end function interpolate_i

!-------------------------------------------------------------------------------
function vector_add(a,b)
! This function adds 2 vectors (arrays real(rprec), dimension(3))
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a,b
real(rprec), dimension(3) :: vector_add
vector_add(1)=a(1)+b(1)
vector_add(2)=a(2)+b(2)
vector_add(3)=a(3)+b(3)
return
end function vector_add

!-------------------------------------------------------------------------------
function vector_divide(a,b)
! This function divides one vector (array real(rprec), dimension(3) by a number)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec), intent(in) :: b
real(rprec), dimension(3) :: vector_divide
vector_divide(1)=a(1)/b
vector_divide(2)=a(2)/b
vector_divide(3)=a(3)/b
return
end function vector_divide

!-------------------------------------------------------------------------------
function vector_multiply(a,b)
! This function multiplies one vector (array real(rprec), dimension(3) by
! a real(rprec) number)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec), intent(in) :: b
real(rprec), dimension(3) :: vector_multiply
vector_multiply(1)=a(1)*b
vector_multiply(2)=a(2)*b
vector_multiply(3)=a(3)*b
return
end function vector_multiply

!-------------------------------------------------------------------------------
function vector_mag(a)
! This function calculates the magnitude of a vector
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec) :: vector_mag
vector_mag=abs(sqrt(a(1)**2+a(2)**2+a(3)**2))
return
end function vector_mag

!-------------------------------------------------------------------------------
function rotatePoint(point_in, rotationPoint, axis, angle)
! This function performs rotation of a point with respect to an axis or rotation
! and a certain angle
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: point_in
real(rprec), dimension(3), intent(in) :: rotationPoint
real(rprec), dimension(3), intent(in) :: axis
real(rprec), intent(in) :: angle
real(rprec), dimension(3,3) :: RM ! Rotation Matrix tensor
real(rprec), dimension(3) :: rotatePoint, point

point=point_in

RM(1,1) = axis(1)**2 + (1.0 - axis(1)**2) * cos(angle)
RM(1,2) = axis(1) * axis(2) * (1.0 - cos(angle)) - axis(3) * sin(angle)
RM(1,3) = axis(1) * axis(3) * (1.0 - cos(angle)) + axis(2) * sin(angle)
RM(2,1) = axis(1) * axis(2) * (1.0 - cos(angle)) + axis(3) * sin(angle)
RM(2,2) = axis(2)**2 + (1.0 - axis(2)**2) * cos(angle)
RM(2,3) = axis(2) * axis(3) * (1.0 - cos(angle)) - axis(1) * sin(angle)
RM(3,1) = axis(1) * axis(3) * (1.0 - cos(angle)) - axis(2) * sin(angle)
RM(3,2) = axis(2) * axis(3) * (1.0 - cos(angle)) + axis(1) * sin(angle)
RM(3,3) = axis(3)**2 + (1.0 - axis(3)**2) * cos(angle)

! Rotation matrices make a rotation about the origin, so need to subtract 
! rotation point off the point to be rotated
point=vector_add(point,-rotationPoint)

! Perform rotation (multiplication matrix and vector)
point=matrix_vector(RM,point)

! Return the rotated point to its new location relative to the rotation point
rotatePoint = point + rotationPoint

return 
end function rotatePoint

!-------------------------------------------------------------------------------
function matrix_vector(RM,point)
! This function multiplies a matrix and a vector
!-------------------------------------------------------------------------------
real(rprec), dimension(3,3), intent(in) :: RM ! Matrix
real(rprec), dimension(3), intent(in) :: point ! vector point
real(rprec), dimension(3) :: matrix_vector
! Perform rotation
matrix_vector(1)=RM(1,1)*point(1)+RM(1,2)*point(2)+RM(1,3)*point(3)
matrix_vector(2)=RM(2,1)*point(1)+RM(2,2)*point(2)+RM(2,3)*point(3)
matrix_vector(3)=RM(3,1)*point(1)+RM(3,2)*point(2)+RM(3,3)*point(3)
return
end function matrix_vector

!-------------------------------------------------------------------------------
function cross_product(u,v)
! This function calculates the cross product of 2 vectors
!-------------------------------------------------------------------------------
real(rprec), intent(in) :: u(3),v(3)
real(rprec) :: cross_product(3)
cross_product(1) = u(2)*v(3)-u(3)*v(2)
cross_product(2) = u(3)*v(1)-u(1)*v(3)
cross_product(3) = u(1)*v(2)-u(2)*v(1)
return
end function cross_product

!-------------------------------------------------------------------------------
function distance(a,b)
! This function calculates the distance between a(1,2,3) and b(1,2,3)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a,b
real(rprec) :: distance
distance=abs(sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2))
return
end function distance


end module actuator_turbine_model

















