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
module atm_input_util
!*******************************************************************************
! This module reads the input files for the actuator turbine model module (ATM)

! Module for dynamic allocation variables
use atm_base

implicit none

! The variables for the ATM are defined here
integer :: numberOfTurbines
integer :: outputInterval

! This type will store the necessary variables for -each- turbine 
! To declare: type(turbineArray_t), allocatable, dimension(:) :: turbineArray
! To access the variables do: turbineArray(n) % variableNeeded  
type turbineArray_t 
    ! Variables that are read from input files
    character(128) :: turbineName ! Name of turbine ('turbine1')
    character(128) :: turbineType ! Name of turbine type ('NREL5MWRef')
    real(rprec), dimension(3) :: baseLocation ! Location of the base (0 0 0)
    integer :: numBladePoints ! Number of points along each blade
    character(128) :: bladeUpdateType ! 
    real(rprec) :: epsilon ! Width of the smearing Gaussian function (5.0)
    character(128) :: rotationDir ! Direction of rotation ('cw')
    real(rprec) :: Azimuth           
    real(rprec) :: RotSpeed           ! Speed of the rotor (rpm)
    real(rprec) :: Pitch              
    real(rprec) :: NacYaw             
    real(rprec) :: fluidDensity      
    integer :: numAnnulusSections ! Number of annulus sections on each blade
    real(rprec) :: AnnulusSectionAngle ! Number of annulus sections on each blade
    real(rprec) :: deltaNacYaw
    ! Not read variables
    real(rprec) :: thrust ! Total turbine thrust
    real(rprec) :: torqueRotor ! Rotor torque
    real(rprec) :: torqueGen ! Generator torque
    real(rprec) :: powerRotor ! Rotor Power

    integer :: turbineTypeID ! Identifies the type of turbine   
    
    !!-- Important geometry data 
    ! Collection of all the actuator points (blade, annular section, point, 3)
    type(real(rprec)), allocatable, dimension(:,:,:,:) :: bladePoints
    ! The solidity at each actuator section  
    type(real(rprec)), allocatable, dimension(:,:,:) :: solidity     
    ! Collection of radius of each point (different because of coning)
    type(real(rprec)), allocatable, dimension(:,:,:) :: bladeRadius
    ! Forces on each actuator point (blade, annular section, point, 3)
    type(real(rprec)), allocatable, dimension(:,:,:,:) :: bladeForces
    ! This dummy variable is used for MPI purposes of doing MPI_SUM on
    ! the forces for all processors, since forces are only computed at
    ! processors which contain the points
    type(real(rprec)), allocatable, dimension(:,:,:,:) :: bladeForcesDummy
    ! Vectors at each actuator point defining the local reference frame
    ! (blade, annular section, point, 3, 3) (three vectors)
    type(real(rprec)), allocatable, dimension(:,:,:,:,:) :: bladeAlignedVectors
    ! The wind U projected onto the bladeAlignedVectors plus rotational speed
    ! (blade, annular section, point, 3, 3) (three vectors)
    type(real(rprec)), allocatable, dimension(:,:,:,:) :: windVectors
    ! Angle of attack at each each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: alpha
    ! Velocity magnitud at each each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: Vmag
    ! Lift coefficient at each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: Cl
    ! Drag coeficient at each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: Cd
    ! Lift at each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: lift
    ! Drag at each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: drag
    ! Axial force at each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: axialForce
    ! Tangential force at each actuator point
    type(real(rprec)), allocatable, dimension(:,:,:) :: tangentialForce


    ! An indicator of shaft direction.  The convention is that when viewed
    ! from upwind, the rotor turns clockwise for positive rotation angles,
    ! regardless of if it is an upwind or downwind turbine.  uvShaft is
    ! found by subtracting the rotor apex location from the tower shaft
    ! intersection point.  This vector switches direciton depending on 
    ! if the turbine is upwind or downwind, so this uvShaftDir multiplier
    ! makes the vector consistent no matter what kind of turbine
     real(rprec) :: uvShaftDir

    ! Define the vector along the shaft pointing in the direction of the wind
    real(rprec), dimension(3) :: uvShaft
  
    ! List of locations of the rotor apex relative to the origin (m)
    real(rprec), dimension(3) :: rotorApex
  
    ! List of locations of the intersection of the tower axis and the shaft 
    ! centerline relative to the origin (m).
    real(rprec), dimension(3) :: towerShaftIntersect
  
     ! Unit vector pointing along the tower (axis of yaw).
    real(rprec), dimension(3) :: uvTower

    ! Width of the actuator section
    type(real(rprec)), allocatable, dimension(:) :: db

    ! Sphere radius which defines a sphere from the center of the rotor and
    ! identifies the volume onto which forces are applied
    real(rprec) :: projectionRadius
    real(rprec) :: sphereRadius

end type turbineArray_t

! This type stores all the airfoils and their AOA, Cd, Cl values
type airfoilType_t
    character(128) :: airfoilName          ! The type of Airfoil ('Cylinder1')
    integer :: n                           ! Number of data points
    ! The maximum number of points is chosen to be 150. If airfoil data has 
    ! more than this then this number should be modified!
    type(real(rprec)), dimension(150) :: AOA    ! Angle of Attack
    type(real(rprec)), dimension(150) :: Cd     ! Drag coefficient
    type(real(rprec)), dimension(150) :: Cl     ! Lift coefficient
    type(real(rprec)), dimension(150) :: Cm     ! Moment coefficient

end type airfoilType_t

type turbineModel_t
    character(128) :: turbineType ! The type of turbine ('NREL5MWRef')
    integer :: NumBl          ! Number of turbine blades
    integer :: NumSec         ! Number of sections
    real(rprec) :: TipRad ! Radius from the root to the tip of the blade
    real(rprec) :: HubRad ! Radius of the hub
    real(rprec) :: UndSling !
    real(rprec) :: OverHang
    real(rprec) :: TowerHt
    real(rprec) :: Twr2Shft
    real(rprec) :: ShftTilt
    real(rprec) :: PreCone      
    real(rprec) :: GBRatio  
    real(rprec) :: GenIner
    real(rprec) :: HubIner
    real(rprec) :: BladeIner
    real(rprec) :: DriveTrainIner
    real(rprec) :: TorqueControllerType
    real(rprec) :: CutInGenSpeed
    real(rprec) :: RatedGenSpeed
    real(rprec) :: Region2StartGenSpeed
    real(rprec) :: Region2EndGenSpeed
    real(rprec) :: CutInGenTorque
    real(rprec) :: RatedGenTorque
    real(rprec) :: RateLimitGenTorque
    real(rprec) :: KGen
    real(rprec) :: TorqueControllerRelax

    ! Blade section quantities
    real(rprec), dimension(25) :: chord, twist, radius
    integer, dimension(25) :: sectionType

    ! The airfoil type properties ( includes AOA, Cl, and Cd) Attempt 1
    type(airfoilType_t), allocatable, dimension(:) :: airfoilType
!    ! The airfoil type properties ( includes AOA, Cl, and Cd) Attempt 2
!    type(DynamicList_), pointer :: allAirfoils

end type turbineModel_t

! Declare turbine array variable
type(turbineArray_t), allocatable, dimension(:) , target :: turbineArray

! Declare turbine model variable (stores information for turbine models)
type(turbineModel_t), allocatable, dimension(:), target :: turbineModel

! Name of the utility used  
character (*), parameter :: mod_name = 'atm_input_util'
character (*), parameter :: input_conf = './inputATM/turbineArrayProperties'
character (*), parameter :: comment = '!'
character (*), parameter :: block_entry = '{' ! The start of a block
character (*), parameter :: block_exit = '}' ! The end of a block 
character (*), parameter :: equal = '='
character (*), parameter :: esyntax = 'syntax error at line'
character (*), parameter :: array_entry = '(' ! The start of an array
character (*), parameter :: array_exit  = ')' ! The end of an array

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Thresh hold for evaluating differences in floating point values.
!real(rprec), parameter :: thresh = 1.0e-6_rprec

! Variables used to read lines in file
integer :: block_entry_pos ! Determines if the line is the start of a block 
integer :: block_exit_pos ! Determines if the line is the end of a block
integer :: array_entry_pos ! Determines if the line is the start of a block 
integer :: array_exit_pos ! Determines if the line is the end of a block
integer :: equal_pos ! Determines if there is an equal sign
integer :: ios 
logical :: exst ! Used to determine existence of a file

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_input_conf()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character (*), parameter :: sub_name = mod_name // '.read_input_conf'
integer :: n = 0 ! Counter for the wind turbines
integer :: lun =1 ! Reference number for input file
integer :: line ! Counts the current line in a file
character (128) :: buff ! Stored the read line
integer, pointer :: numBladePoints, numAnnulusSections

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

! Open file
if (exst) then
    ! Open the input file
    open (lun, file=input_conf, action='read')
else
    ! Error for non existing file
    call error ('file ' // input_conf // ' does not exist')
endif

! Read the file line by line *Counter starts at 0 and modified inside subroutine
line = 0
do
! Read line by line (lun=file number) 
    call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                   array_entry_pos, array_exit_pos, equal_pos, ios )
                   
    if (ios /= 0) exit ! Exit if reached end of file

    ! This will read the numberOfTurbines integer
        if( buff(1:16) == 'numberOfTurbines' ) then
            read(buff(17:), *) numberOfTurbines
            write(*,*) 'numberOfTurbines is: ', numberOfTurbines
            ! Allocate space for the wind turbine variables
            allocate(turbineArray(numberOfTurbines))
            cycle
        else if( buff(1:14) == 'outputInterval' ) then
            read(buff(15:), *) outputInterval
            write(*,*)  'outputInterval is: ', outputInterval
            cycle
        endif

    if (block_entry_pos /= 0) then ! This will start reading turbine block
        n = n + 1     ! Increment turbine counter
        if (n .gt. numberOfTurbines) exit
        ! Read the name of the turbine
        read(buff(1:index(buff, block_entry)-1), *) turbineArray(n) &
        % turbineName
    endif  
    if (block_entry_pos == 0) then ! This will start reading turbine block
        if( buff(1:11) == 'turbineType' ) then
            read(buff(12:), *) turbineArray(n) % turbineType
            write(*,*)  'turbineType is: ', turbineArray(n) % turbineType
        endif
        if( buff(1:12) == 'baseLocation' ) then
            read(buff(13:), *) turbineArray(n) % baseLocation  
            write(*,*)  'baseLocation is: ', turbineArray(n) % baseLocation      
        endif        
        if( buff(1:14) == 'numBladePoints' ) then
            read(buff(15:), *) turbineArray(n) % numBladePoints
            write(*,*)  'numBladePoints is: ', turbineArray(n) % numBladePoints
            ! Allocation depending on the number of blade points
        endif
        if( buff(1:7) == 'epsilon' ) then
            read(buff(8:), *) turbineArray(n) % epsilon
            write(*,*)  'epsilon is: ', turbineArray(n) % epsilon
        endif
        if( buff(1:11) == 'rotationDir' ) then
            read(buff(12:), *) turbineArray(n) % rotationDir
            write(*,*)  'rotationDir is: ', turbineArray(n) % rotationDir
        endif
        if( buff(1:7) == 'Azimuth' ) then
            read(buff(8:), *) turbineArray(n) % Azimuth
            write(*,*)  'Azimuth is: ', turbineArray(n) % Azimuth
        endif
        if( buff(1:8) == 'RotSpeed' ) then
            read(buff(9:), *) turbineArray(n) % RotSpeed
            write(*,*)  'RotSpeed is: ', turbineArray(n) % RotSpeed
        endif
        if( buff(1:5) == 'Pitch' ) then
            read(buff(6:), *) turbineArray(n) % Pitch
            write(*,*)  'Pitch is: ', turbineArray(n) % Pitch
        endif
        if( buff(1:6) == 'NacYaw' ) then
            read(buff(7:), *) turbineArray(n) % NacYaw
            write(*,*)  'NacYaw is: ', turbineArray(n) % NacYaw
        endif
        if( buff(1:12) == 'fluidDensity' ) then
            read(buff(13:), *) turbineArray(n) % fluidDensity
            write(*,*)  'fluidDensity is: ', turbineArray(n) % fluidDensity
        endif
        if( buff(1:18) == 'numAnnulusSections' ) then
            read(buff(19:), *) turbineArray(n) % numAnnulusSections
            write(*,*)  'numAnnulusSections is: ', &
                         turbineArray(n) % numAnnulusSections
        endif      
        if( buff(1:19) == 'annulusSectionAngle' ) then
            read(buff(20:), *) turbineArray(n) % annulusSectionAngle
            write(*,*)  'annulusSectionAngle is: ', &
                         turbineArray(n) % annulusSectionAngle
        endif   
    endif        
end do

if( .not. allocated( turbineArray ) ) then 
    write(*,*) 'Did not allocate memory for turbineArray' 
stop 
endif
close (lun)

call read_turbine_model_variables ()

end subroutine read_input_conf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_turbine_model_variables ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
integer :: i,j ! counter
integer :: numTurbinesDistinct ! Number of different turbine types
character(128) :: currentTurbineType ! Will store turbineType in loop
character(128) :: input_turbine
integer :: lun =19  ! Reference number for input file
integer :: line ! Counts the current line in a file
character (128) :: buff ! Stored the read line
integer:: numAirfoils ! Number of distinct airfoils
integer :: k, p ! Used to loop through aifoil types and character counter
integer :: numAnnulusSections, numBladePoints, numBl, numSec

! Name of all the airfoils types (max 20) If more than this increase the number
character(128), dimension(20) :: airfoils 

! Initialize variables for the loop
! Will find the number of distincet turbines to allocate memory for turbineModel
numTurbinesDistinct=1
currentTurbineType=turbineArray(1) % turbineType
! Find how many turbine types there are
do i=1,numberOfTurbines
    if (turbineArray(i) % turbineType .ne. currentTurbineType) then
        numTurbinesDistinct = numTurbinesDistinct+1
        currentTurbineType=turbineArray(i) % turbineType
    endif
    turbineArray(i) % turbineTypeID=numTurbinesDistinct
enddo

! Allocate space for turbine model variables
allocate(turbineModel(numTurbinesDistinct))

! This will store the turbine types on each turbine model ("NREL5MW")
numTurbinesDistinct=1
currentTurbineType=turbineArray(1) % turbineType
turbineModel(numTurbinesDistinct) % turbineType = turbineArray(1) % turbineType
do i=1,numberOfTurbines
    if (turbineArray(i) % turbineType .ne. currentTurbineType) then
    numTurbinesDistinct = numTurbinesDistinct+1
    turbineModel(numTurbinesDistinct) % turbineType = &
    turbineArray(i) % turbineType
    endif
enddo

! Read the input properties for each turbine type
do i = 1, numTurbinesDistinct
    input_turbine = './inputATM/' // turbineModel(i) % turbineType 
    ! Check that the configuration file exists
    inquire (file=input_turbine, exist=exst)
    ! Open file
    if (exst) then
        ! Open the input file
        open (lun, file=input_turbine, action='read')
    else
        ! Error for non existing file
        call error ('file ' // input_turbine // ' does not exist')
    endif

    write(*,*) 'Reading Turbine Model Properties for: ', &
                turbineModel(i) % turbineType
    ! Read the file line by line - starts at 0 and modified inside subroutine
    line = 0
    do
    ! Read line by line (lun=file number) 
        call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                       array_entry_pos, array_exit_pos, equal_pos, ios )
        if (ios /= 0) exit
        ! This will all the input variables
        if( buff(1:5) == 'NumBl' ) then
            read(buff(6:), *) turbineModel(i) % NumBl
            write(*,*) 'NumBl is: ', turbineModel(i) % NumBl
        endif
        if( buff(1:6) == 'TipRad' ) then
            read(buff(7:), *) turbineModel(i) % TipRad
            write(*,*) 'TipRad is: ', turbineModel(i) % TipRad
        endif
        if( buff(1:6) == 'HubRad' ) then
            read(buff(7:), *) turbineModel(i) % HubRad
            write(*,*) 'HubRad is: ', turbineModel(i) % HubRad
        endif
        if( buff(1:8) == 'UndSling' ) then
            read(buff(9:), *) turbineModel(i) % UndSling
            write(*,*) 'UndSling is: ', turbineModel(i) % UndSling
        endif
        if( buff(1:8) == 'OverHang' ) then
            read(buff(9:), *) turbineModel(i) % OverHang
            write(*,*) 'OverHang is: ', turbineModel(i) % OverHang
        endif
        if( buff(1:7) == 'TowerHt' ) then
            read(buff(8:), *) turbineModel(i) % TowerHt
            write(*,*) 'TowerHt is: ', turbineModel(i) % TowerHt
        endif
        if( buff(1:8) == 'Twr2Shft' ) then
            read(buff(9:), *) turbineModel(i) % Twr2Shft
            write(*,*) 'Twr2Shft is: ', turbineModel(i) % Twr2Shft
        endif
        if( buff(1:8) == 'ShftTilt' ) then
            read(buff(9:), *) turbineModel(i) % ShftTilt
            write(*,*) 'ShftTilt is: ', turbineModel(i) % ShftTilt
        endif
        if( buff(1:7) == 'PreCone' ) then
            read(buff(8:), *) turbineModel(i) % PreCone
            write(*,*) 'PreCone is: ', turbineModel(i) % PreCone
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the airfoils
        if ( buff(1:8) == 'Airfoils' ) then ! Start reading airfoil block
            numAirfoils=0 ! Conuter for the number of distince airfoils
            array_entry_pos=0 ! If 'Airfoils(' then make array_entry_pos zero 

            do while (array_entry_pos == 0) 
                call readline( lun, line, buff, block_entry_pos,              &
                block_exit_pos, array_entry_pos, array_exit_pos, equal_pos, ios)
                if (ios /= 0) exit        ! exit if end of file reached
                
                if (array_entry_pos /= 0) then
                    call readline( lun, line, buff, block_entry_pos,          &
                                   block_exit_pos, array_entry_pos,           &
                                   array_exit_pos, equal_pos, ios)
                   if (ios /= 0) exit        ! exit if end of file reached
                endif

                if (array_exit_pos /= 0) exit     ! exit if end of file reached
                numAirfoils = numAirfoils + 1     ! Increment airfoil counter
                airfoils(numAirfoils)=buff ! Stores the name of the airfoil
            enddo
            
            ! Allocate the airfoilTypes
            allocate(turbineModel(i) % airfoilType(numAirfoils))

            ! Loop through all the airfoil types and look-up lift and drag
            do k=1,numAirfoils
                call eat_whitespace (airfoils(k)) ! Eliminate white space
                p=len(trim(airfoils(k)))-1        ! Length without last element
                ! Airfoil type (2:p) is used to eliminate the commas ""
                turbineModel(i) % airfoilType % airfoilName =  airfoils(k)(2:p)
                ! Read each airfoil accordingly
                call read_airfoil( turbineModel(i) % airfoilType(k) ) 
            enddo
        endif  ! End of Airfoil loop

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the blade properties
        if ( buff(1:9) .eq. 'BladeData' ) then ! Start reading blade data block
                NumSec = 0
                do while (array_exit_pos .eq. 0) 
                    read (lun, '(a)', iostat=ios) buff ! Read the comment line
                    if (scan(buff,'(' ) .ne. 0) cycle
                    if (scan(buff,'!' ) .ne. 0) cycle
                    if (scan(buff,')' ) .ne. 0) exit

                    ! Number of sections will account for all blade sections
                    NumSec=NumSec+1
                    turbineModel(i) % NumSec = NumSec

                    ! Read in radius, chord, twist, type
                    read(buff,*) turbineModel(i) % radius(NumSec),           &
                    turbineModel(i) % chord(NumSec),   &
                    turbineModel(i) % twist(NumSec),   &
                    turbineModel(i) % sectionType(NumSec)

                    ! Add one to airfoil identifier. List starts at 0, now 
                    ! it will start at 1
                    turbineModel(i) % sectionType( NumSec )  =               &
                    turbineModel(i) % sectionType( NumSec ) + 1
                enddo
        endif
    enddo
    close (lun)      
enddo

! Allocate variables inside turbineArray
do  i=1,numberOfTurbines
numBladePoints = turbineArray(i) % numBladePoints
numAnnulusSections = turbineArray(i) % numAnnulusSections
j=turbineArray(i) % turbineTypeID
numBl=turbineModel(j) % numBl

    allocate(turbineArray(i) % bladeForces(numBl,          &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % bladeForcesDummy(numBl,          &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % bladeAlignedVectors(numBl,  &
             numAnnulusSections, numBladePoints,3,3) )
    allocate(turbineArray(i) % windVectors(numBl,          &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % alpha(numBl,                &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Vmag(numBl,                 &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Cl(numBl,                   &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Cd(numBl,                   &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % lift(numBl,                 &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % drag(numBl,                 &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % axialForce(numBl,           &
             numAnnulusSections, numBladePoints))
    allocate(turbineArray(i) % tangentialForce(numBl,      &
             numAnnulusSections, numBladePoints) )

enddo

end subroutine read_turbine_model_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_airfoil( airfoilType )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine reads the angle of attack, lift and drag for a specific 
! airfoil type
integer :: lun=17 ! File identifier
integer :: q ! Counter
character(128) :: input_airfoil
type(airfoilType_t), intent(inout) :: airfoilType
real(rprec) :: AOA, Cd, Cl, Cm
! Name of the input file to read
input_airfoil= './inputATM/AeroData/' //trim (airfoilType % airfoilName)  &
               // '.dat'
write (*,*) input_airfoil
! Open airfoil input file
open (lun, file=input_airfoil, action='read')

! Skip the first 14 lines of the Fortran input file 
do q=1,14
    read(lun,*)
enddo

q=0 ! Initialize the counter back to 1 for the first element in the list

AOA=-181.
do while (AOA .lt. 180.00)
    q=q+1
    read(lun,*) AOA, Cl , Cd, Cm 
    airfoilType % AOA(q) = AOA
    airfoilType % Cd(q) = Cd
    airfoilType % Cl(q) = Cl
    airfoilType % Cm(q) = Cm
enddo
airfoilType % n = q

close(lun)

end subroutine read_airfoil

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine readline(lun, line, buff, block_entry_pos, block_exit_pos, & 
                    array_entry_pos, array_exit_pos, equal_pos, ios )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! This subroutine reads the specified line and determines the attributes
! of the contents of the line.
!
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line

character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, block_exit_pos, equal_pos, ios, &
                        array_entry_pos, array_exit_pos

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do     
    line = line + 1
    read (lun, '(a)', iostat=ios) buff

    if (ios /= 0) exit
    
    call eat_whitespace (buff)

    if (verify (buff, ' ') == 0) cycle  !--drop blank lines
  
    if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

    block_entry_pos = index( buff, block_entry )
    block_exit_pos  = index( buff, block_exit )
    array_entry_pos = index( buff, array_entry )
    array_exit_pos  = index( buff, array_exit )
    equal_pos       = index( buff, equal )
    exit
enddo 

return
end subroutine readline

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine eat_whitespace (buff, whtspc)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! eats leading and intermediate whitespace, fill trailing space with
! blanks
!
implicit none

character (*), intent (inout) :: buff
character (*), intent (in), optional :: whtspc 
character (*), parameter :: whtspc_default = achar (9) // achar (32)
                            !--add more characters here if needed
character (1), parameter :: fill_char = ' '
character (1) :: tmp (len (buff))
character (1) :: fill (len (buff))

fill = fill_char
tmp = transfer (buff, tmp)

if (present (whtspc)) then
  tmp = pack (tmp, scan (tmp, whtspc) == 0, fill)
else
  tmp = pack (tmp, scan (tmp, whtspc_default) == 0, fill)
end if

buff = transfer (tmp, buff)

end subroutine eat_whitespace

end module atm_input_util































