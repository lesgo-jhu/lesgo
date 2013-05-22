!*******************************************************************************
module atm_input_util
!*******************************************************************************
! This module reads the input files for the actuator turbine model module (ATM)
implicit none

! Precision of real numbers
integer, parameter :: rprec = kind (1.d0)

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
    
    ! Not read variables
    integer :: turbineTypeID ! Identifies the type of turbine   
    
    !!-- Important geometry data 
    ! Collection of all the actuator points (blade, annular section, point, 3)
    type(real(rprec)), allocatable, dimension(:,:,:,:) :: bladePoints       
    ! Collection of radius of each point (different because of coning)
    type(real(rprec)), allocatable, dimension(:,:,:) :: bladeRadius
    
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

type turbineModel_t
    character(128) :: turbineType ! The type of turbine ('NREL5MWRef')
    integer :: NumBl ! Number of turbine blades
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

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Thresh hold for evaluating differences in floating point values.
!real(rprec), parameter :: thresh = 1.0e-6_rprec

! Variables used to read lines in file
integer :: block_entry_pos ! Determines if the line is the start of a block 
integer :: block_exit_pos ! Determines if the line is the end of a block
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
                   equal_pos, ios )
                   
    if (ios /= 0) exit
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
integer :: lun =1 ! Reference number for input file
integer :: line ! Counts the current line in a file
character (128) :: buff ! Stored the read line

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

! This will store the turbine types on each turbine model
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
    ! Read the file line by line *Counter starts at 0 and modified inside subroutine
    line = 0
    do
    ! Read line by line (lun=file number) 
        call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                       equal_pos, ios )
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
    enddo                                  
close (lun)      
enddo

end subroutine read_turbine_model_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine readline(lun, line, buff, block_entry_pos, block_exit_pos, & 
                    equal_pos, ios )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! This subroutine reads the specified line and determines the attributes
! of the contents of the line.
!
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line

character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, block_exit_pos, equal_pos, ios

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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine error (msg)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
character (*), intent (in) :: msg

write (*, '(1x,a)') '*****ERROR*****'
write (*, '(1x,a)') 'In atm_input_util:'
write (*, '(1x,a)') trim (msg)
write (*, '(1x,a)') '***************'
write (*, '(1x,a)') 'Program aborted'

stop
end subroutine error









end module atm_input_util































