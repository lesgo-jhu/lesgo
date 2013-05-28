!*******************************************************************************
module actuator_turbine_model
!*******************************************************************************
! This module has the subroutines to provide all calculations for use in the 
! actuator turbine model (ATM)

! Imported modules
use atm_input_util ! Utilities to read input files

implicit none

! Declare everything private except for subroutine which will be used
private 
public :: atm_initialize, atm_forcing, numberOfTurbines

! These are used to do unit conversions
real(rprec), parameter :: pi= 3.141592653589793238462643383279502884197169399375
real(rprec) :: degRad = pi/180. ! Degrees to radians conversion
real(rprec) :: rpmRadSec =  pi/30. ! Set the revolutions/min to radians/s 

integer :: i, j, k ! Counters

logical :: pastFirstTimeStep ! Establishes if we are at the first time step

! Pointers to be used from the turbineArray and Turbinemodel modules
! It is very important to have the pointers pointing into the right variable
! in each subroutine
type(real(rprec)),  pointer :: db(:)
type(real(rprec)),  pointer :: bladePoints(:,:,:,:)
type(real(rprec)),  pointer :: bladeRadius(:,:,:)
integer,     pointer :: numBladePoints
integer,     pointer :: numBl
integer,     pointer :: numAnnulusSections
real(rprec), pointer :: NacYaw             
real(rprec), pointer :: azimuth
real(rprec), pointer :: rotSpeed
real(rprec), pointer :: ShftTilt
real(rprec), pointer :: towerShaftIntersect(:)
real(rprec), pointer :: baseLocation(:)
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: uvShaft(:)
real(rprec), pointer :: uvTower(:)
real(rprec), pointer :: TowerHt
real(rprec), pointer :: Twr2Shft
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: uvShaftDir
real(rprec), pointer :: TipRad
real(rprec), pointer :: HubRad
real(rprec), pointer :: PreCone
real(rprec), pointer :: projectionRadius  ! Radius up to which forces are spread
real(rprec), pointer :: sphereRadius ! Radius of the sphere of forces 

! Subroutines for the actuator turbine model 
! All suboroutines names start with (atm_) 
contains 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize()
! This subroutine initializes the ATM. It calls the subroutines in
! atm_input_util to read the input data and creates the initial geometry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
pastFirstTimeStep=.false. ! The first time step not reached yet

call read_input_conf()  ! Read input data

do i = 1,numberOfTurbines
    call atm_create_points(i)   ! Creates the ATM points defining the geometry
    call calculate_variables(i) ! Calculates variables depending on input
end do
pastFirstTimeStep=.true. ! Past the first time step

end subroutine atm_initialize


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_update(dt)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) :: dt                            ! Time step

do i = 1, numberOfTurbines
    call atm_rotateBlades(dt,i)              ! Rotate the blades of each turbine
end do

end subroutine atm_update


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_rotateBlades(dt,i)
! This subroutine rotates the turbine blades 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer :: i                                 ! Turbine number
integer :: j                                 ! Turbine type
integer :: k, n, m                           ! Counters tu be used in do loops
real(rprec) :: dt                            ! time step
real(rprec) :: deltaAzimuth, deltaAzimuthI   ! Angle of rotation


! Variables which are used by pointers
rotorApex=> turbineArray(i) % rotorApex
j=turbineArray(i) % turbineTypeID
rotSpeed=>turbineArray(i) % rotSpeed
uvShaft=>turbineArray(i) % uvShaft
azimuth=>turbineArray(i) % azimuth

! Angle of rotation
deltaAzimuth = rotSpeed * dt;

! Check the rotation direction first and set the local delta azimuth
! variable accordingly.
if (turbineArray(i) % rotationDir == "cw") then
	   deltaAzimuthI = deltaAzimuth
else if (turbineArray(i) % rotationDir == "ccw") then
    deltaAzimuthI =-deltaAzimuth
end if

do m=1, turbineArray(i) % numAnnulusSections
    do n=1, turbineArray(i) % numBladePoints
        do k=1, turbineModel(j) % numBl
            turbineArray(i) %   bladePoints(k,n,m,:)=rotatePoint(              &
            turbineArray(i) % bladePoints(k,n,m,:), rotorApex, uvShaft,        &
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
subroutine atm_create_points(i)
! This subroutine generate the set of blade points for each turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
integer :: m ! Indicates the blade point number
real(rprec), dimension (3) :: root ! Location of rotor apex
real(rprec) :: beta ! Difference between coning angle and shaft tilt
real(rprec) :: dist ! Distance from each actuator point
! Width of the actuator section

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Variables to be used locally. They are stored in local variables within the 
! subroutine for easier code following. The values are then passed to the 
! proper type
numBladePoints => turbineArray(i) % numBladePoints
numBl=>turbineModel(j) % numBl
numAnnulusSections=>turbineArray(i) % numAnnulusSections


nacYaw=>turbineArray(i) % nacYaw

! Allocate variables depending on specific turbine properties and general
! turbine model properties
allocate(turbineArray(i) % db(numBladePoints))

allocate(turbineArray(i) % bladePoints(numBl, numAnnulusSections, &
         numBladePoints,3))
         
allocate(turbineArray(i) % bladeRadius(numBl,numAnnulusSections,numBladePoints))  

! Assign Pointers
db=>turbineArray(i) % db
bladePoints=>turbineArray(i) % bladePoints
bladeRadius=>turbineArray(i) % bladeRadius
azimuth=>turbineArray(i) % azimuth
rotSpeed=>turbineArray(i) % rotSpeed
ShftTilt=>turbineModel(j) % ShftTilt
preCone=>turbineModel(j) % preCone
towerShaftIntersect=>turbineArray(i) % towerShaftIntersect
baseLocation=>turbineArray(i) % baseLocation
TowerHt=>turbineModel(j) % TowerHt
Twr2Shft=> turbineModel(j) % Twr2Shft
rotorApex=>turbineArray(i) % rotorApex
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
uvShaftDir=>turbineArray(i) % uvShaftDir
uvShaft=>turbineArray(i) % uvShaft
uvTower=>turbineArray(i) % uvTower
TipRad=>turbineModel(j) % TipRad
HubRad=>turbineModel(j) % HubRad
PreCone=>turbineModel(j) %PreCone

!!-- Do all proper conversions for the required variables
! Convert nacelle yaw from compass directions to the standard convention
call compassToStandard(nacYaw)
! Turbine specific
azimuth = degRad*(azimuth)
rotSpeed = rpmRadSec * rotSpeed
nacYaw =degRad * nacYaw
! Turbine model specific
shftTilt = degRad *  shftTilt 
preCone =degRad * preCone

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
!                  Create the first set of actuator points                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the vector along the shaft pointing in the direction of the wind
uvShaftDir = OverHang / abs( OverHang )
! Define the vector along the shaft pointing in the direction of the wind                               
uvShaft = vector_add(rotorApex , - towerShaftIntersect)
uvShaft = vector_divide(uvShaft, vector_mag(uvShaft))
! Define vector aligned with the tower pointing from the ground to the nacelle
uvTower = vector_add(towerShaftIntersect, - baseLocation)
uvTower = vector_divide( uvTower, vector_mag(uvTower))
! Define thickness of each blade section
do k=1, numBladePoints
    db(k) = (TipRad- HubRad)/(numBladePoints)
enddo

root = rotorApex
beta = PreCone - ShftTilt

! This creates the first set of points
do k=1, numBl
    root(1)= root(1) + HubRad*sin(beta)
    root(3)= root(3) + HubRad*cos(beta)
    dist = HubRad
    do m=1, numBladePoints
        dist = dist + 0.5*(db(k))
        bladePoints(k,1,m,1) = root(1) + dist*sin(beta)
        bladePoints(k,1,m,2) = root(2)
        bladePoints(k,1,m,3) = root(3) + dist*cos(beta);
        bladeRadius(k,1,m) = dist;
        dist = dist + 0.5*(db(k))
    enddo
    if (k > 1) then
        do m=1, numBladePoints
            bladePoints(k,1,m,:)=rotatePoint(bladePoints(k,1,m,:), rotorApex, &
            uvShaft,(360.0/NumBl)*k*degRad)
        enddo
    endif
    
enddo


end subroutine atm_create_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine calculate_variables(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Calculates the variables of the model that need information from the input
! files. It runs after reading input information.
integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Declare the required pointers
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
projectionRadius=>turbineArray(i) % projectionRadius
sphereRadius=>turbineArray(i) % sphereRadius
TipRad=>turbineModel(j) % TipRad
PreCone=>turbineModel(j) %PreCone

! First compute the radius of the force projection (to the radius where the 
! projection is only 0.0001 its maximum value - this seems to recover 99.99% of 
! the total forces when integrated
projectionRadius= turbineArray(i) % epsilon * sqrt(log(1.0/0.0001))

sphereRadius=sqrt(((OverHang + UndSling) + TipRad*sin(PreCone))**2 &
+ (TipRad*cos(PreCone))**2) + projectionRadius

end subroutine calculate_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_forcing()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine atm_forcing

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
real(rprec) :: b
real(rprec), dimension(3) :: vector_divide
vector_divide(1)=a(1)/b
vector_divide(2)=a(2)/b
vector_divide(3)=a(3)/b
return
end function vector_divide

!-------------------------------------------------------------------------------
function vector_mag(a)
! This function calculates the magnitude of a vector
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a
real(rprec) :: vector_mag
vector_mag=abs(sqrt(a(1)**2+a(2)**2+a(3)**2))
return
end function vector_mag

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
return 
end subroutine compassToStandard

!-------------------------------------------------------------------------------
function rotatePoint(point, rotationPoint, axis, angle)
! This function performs rotation of a point with respect to an axis or rotation
! and a certain angle
!-------------------------------------------------------------------------------
real(rprec), dimension(3) :: point
real(rprec), dimension(3) :: rotationPoint
real(rprec), dimension(3) :: axis
real(rprec) :: angle
real(rprec), dimension(3,3) :: RM ! Rotation Matrix tensor
real(rprec), dimension(3) :: rotatePoint

RM(1,1) = sqrt(axis(1)) + (1.0 - sqrt(axis(1))) * cos(angle); 
RM(1,2) = axis(1) * axis(2) * (1.0 - cos(angle)) - axis(3) * sin(angle); 
RM(1,3) = axis(1) * axis(3) * (1.0 - cos(angle)) + axis(2) * sin(angle);
RM(2,1) = axis(1) * axis(2) * (1.0 - cos(angle)) + axis(3) * sin(angle); 
RM(2,2) = sqrt(axis(2)) + (1.0 - sqrt(axis(2))) * cos(angle);
RM(2,3) = axis(2) * axis(3) * (1.0 - cos(angle)) - axis(1) * sin(angle);
RM(3,1) = axis(1) * axis(3) * (1.0 - cos(angle)) - axis(2) * sin(angle);
RM(3,2) = axis(2) * axis(3) * (1.0 - cos(angle)) + axis(1) * sin(angle);
RM(3,3) = sqrt(axis(3)) + (1.0 - sqrt(axis(3))) * cos(angle);

! Rotation matrices make a rotation about the origin, so need to subtract 
! rotation point off the point to be rotated
point=vector_add(point,-rotationPoint)

! Perform rotation
rotatePoint(1)=RM(1,1)*point(1)+RM(1,2)*point(2)+RM(1,3)*point(3)
rotatePoint(2)=RM(2,1)*point(1)+RM(2,2)*point(2)+RM(2,3)*point(3)
rotatePoint(3)=RM(3,1)*point(1)+RM(3,2)*point(2)+RM(3,3)*point(3)

! Return the rotated point to its new location relative to the rotation point
rotatePoint=rotatePoint+rotationPoint

return 
end function rotatePoint




end module actuator_turbine_model

















