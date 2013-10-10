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
! actuator airfoil model (atm)

! Imported modules
use atm_base ! Include basic types and precission of real numbers

implicit none

! Declare everything private except for subroutine which will be used
private 
public :: atm_initialize, projectionRadius, epsilon, bladeForces, outputInterval,  &
          cd, cl, alpha, lift, drag, windVectors, bladePoints, numOfPoints,  &
          bladeScalarDummy, bladeVectorDummy, Vmag, atm_convoluteForce,  &
          atm_computeBladeForce, atm_output, atm_initialize_output

! The very crucial parameter pi
real(rprec), parameter :: pi=acos(-1._rprec)
real(rprec) :: degRad = pi/180._rprec ! Degrees to radians conversion

! Variables made public
integer, parameter :: numOfPoints=100, lt_n = 3 
integer, parameter :: N=100, outputInterval=20
real(rprec) :: epsilon, projectionRadius, twistAng, chord, db
real(rprec) :: cl(N), cd(N),cm(N), alpha(N), lift(N), drag(N)
real(rprec) :: bladePoints(N,3), bladeForces(N,3), bladeRadius(N)
real(rprec) :: bladeAlignedVectors(N,3), windVectors(N,3), Vmag(N)
real(rprec) :: bladeVectorDummy(N,3), bladeScalarDummy(N)

real(rprec) :: AOA(lt_N), cd_table(lt_N), cl_table(lt_N), cm_table(lt_N)

logical :: pastFirstTimeStep ! Establishes if we are at the first time step

! Subroutines for the actuator turbine model 
! All suboroutines names start with (atm_) 
contains 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize()
! This subroutine initializes the atm. It calls the subroutines in
! atm_input_util to read the input data and creates the initial geometry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer :: i
real(rprec) :: read_AOA, read_Cl, read_Cd, read_Cm
! Initialize all variables
epsilon=5.5
projectionRadius= epsilon * sqrt(log(1.0/0.0001))
twistAng = 13.
chord = 4.

! Read in the airfoil variables
open (1, file='Cylinder1.dat', action='read')
! Skip the first 14 lines of the Fortran input file 
do i=1,14
    read(1,*)
enddo
do i=1,lt_n
    read(1,*) read_AOA, read_Cl, read_Cd, read_Cm
    AOA(i)=read_AOA
    Cl(i)=read_Cl
    Cd(i)=read_Cd
    Cm(i)=read_Cm
enddo
write(*,*) AOA
write(*,*) Cd
write(*,*) Cl

close(1)

end subroutine atm_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize_output()
! This subroutine initializes the output files for the atm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
logical :: file_exists

inquire(file='./airfoilOutput',EXIST=file_exists)

if (file_exists .eqv. .false.) then

    ! Create turbineOutput directory
    call system("mkdir -vp airfoilOutput")

    open(unit=1, file="./airfoilOutput/lift")
    write(1,*) 'turbineNumber bladeNumber '
    close(1)

    open(unit=1, file="./airfoilOutput/drag")
    write(1,*) 'turbineNumber bladeNumber '
    close(1)
    
    open(unit=1, file="./airfoilOutput/Cl")
    write(1,*) 'turbineNumber bladeNumber Cl'
    close(1)

    open(unit=1, file="./airfoilOutput/Cd")
    write(1,*) 'turbineNumber bladeNumber Cd'
    close(1)

    open(unit=1, file="./airfoilOutput/alpha")
    write(1,*) 'turbineNumber bladeNumber alpha'
    close(1)

    open(unit=1, file="./airfoilOutput/Vrel")
    write(1,*) 'turbineNumber bladeNumber Vrel'
    close(1)

    open(unit=1, file="./airfoilOutput/Vaxial")
    write(1,*) 'turbineNumber bladeNumber Vaxial'
    close(1)

    open(unit=1, file="./airfoilOutput/Vtangential")
    write(1,*) 'turbineNumber bladeNumber Vtangential'
    close(1)

endif

end subroutine atm_initialize_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_create_points()
! This subroutine generates the set of blade points for each turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: m ! Indicates the blade point number
real(rprec), dimension (3) :: root ! Location of rotor apex

root(1)=0.5
root(2)=0.5
root(3)=-2*pi
db=2*pi/numOfPoints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Create the set of actuator points                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bladePoints(1,:)=root
do m=2, numOfPoints
    bladePoints(m,1)=root(1)
    bladePoints(m,3)=root(3)
    bladePoints(m,2)=bladePoints(m,2)+db
    bladeRadius(m)=bladePoints(m,2)
enddo

end subroutine atm_create_points


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeBladeForce(i,U_local)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the wind vectors by projecting the velocity 
! onto the transformed coordinates system
implicit none

integer, intent(in) :: i

real(rprec), intent(in) :: U_local(3)    ! The local velocity at this point

! Local variables
real(rprec) :: windAng_i
real(rprec), dimension(3) :: dragVector, liftVector

! Pointers to trubineArray (i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This will compute the vectors defining the local coordinate 
! system of the actuator point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define vector in z'
bladeAlignedVectors(3,:) =  (/ 0, 0, 1 /)

! Define vector in y'
bladeAlignedVectors(2,:) = (/ 0, 1, 0 /)

! Define vector in x'
bladeAlignedVectors(1,:) = (/ 1, 0, 0 /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
! This concludes the definition of the local coordinate system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now put the velocity in that cell into blade-oriented coordinates and add on 
! the velocity due to blade rotation.
windVectors(i,1) = dot_product(bladeAlignedVectors(1,:) , U_local)
windVectors(i,2) = dot_product(bladeAlignedVectors(2,:), U_local) 
windVectors(i,3) = dot_product(bladeAlignedVectors(3,:), U_local)

! Velocity magnitude
Vmag(i)=sqrt( windVectors(i,1)**2+windVectors(i,2)**2 )

! Angle between wind vector components
windAng_i = atan2( windVectors(i,1), windVectors(i,2) ) /degRad

! Local angle of attack
alpha(i) = windAng_i - twistAng

! Lift coefficient
cl(i)= interpolate(alpha(i), AOA, cl )

! Drag coefficient
cd(i)= interpolate(alpha(i), AOA, cd )

! Lift force
lift(i) = 0.5_rprec * cl(i) * (Vmag(i)**2) * chord * db

! Drag force
drag(i) = 0.5_rprec * cd(i) * (Vmag(i)**2) * chord * db

! This vector projects the drag onto the local coordinate system
dragVector = bladeAlignedVectors(1,:)*windVectors(i,1) +  &
             bladeAlignedVectors(2,:)*windVectors(i,2)

dragVector = vector_divide(dragVector,vector_mag(dragVector) )

! Lift vector
liftVector = cross_product(dragVector,bladeAlignedVectors(3,:) )
liftVector = liftVector/vector_mag(liftVector)

! Apply the lift and drag as vectors
liftVector = -lift(i) * liftVector
dragVector = -drag(i) * dragVector

! The blade force is the total lift and drag vectors 
bladeForces(i,:) = vector_add(liftVector, dragVector)

! Calculate output quantities based on each point
!call atm_process_output(i,m,n,q)

end subroutine atm_computeBladeForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_output(jt_total, outputInterval )
! This subroutine will calculate the output of the model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

! Number of iteration fed in from solver
integer, intent(in) :: jt_total, outputInterval 

integer :: powerFile=11, bladeFile=12, liftFile=13, dragFile=14
integer :: ClFile=15, CdFile=16, alphaFile=17, VrelFile=18
integer :: VaxialFile=19, VtangentialFile=20

! Output only if the number of intervals is right
if ( mod(jt_total-1, outputInterval) == 0) then
        
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!  Writing Actuator Turbine Model output  !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

        ! File for blade output
        open(unit=bladeFile,position="append",                                 &
        file="./airfoilOutput/blade")

        open(unit=liftFile,position="append",                                  &
        file="./airfoilOutput/lift")

        open(unit=dragFile,position="append",                                  &
        file="./airfoilOutput/drag")
        
        open(unit=ClFile,position="append",                                    &
        file="./airfoilOutput/Cl")

        open(unit=CdFile,position="append",                                    &
        file="./airfoilOutput/Cd")

        open(unit=alphaFile,position="append",                                 &
        file="./airfoilOutput/alpha")

        open(unit=VrelFile,position="append",                                 &
        file="./airfoilOutput/Vrel")

        open(unit=VaxialFile,position="append",                                 &
        file="./airfoilOutput/Vaxial")

        open(unit=VtangentialFile,position="append",                                 &
        file="./airfoilOutput/Vtangential")

        ! Will write only the first actuator section of the blade
            write(bladeFile,*) bladeRadius(:)
            write(liftFile,*) lift(:)/ db
            write(dragFile,*) drag(:)/ db
            write(ClFile,*) cl(:)
            write(CdFile,*) cd(:)
            write(alphaFile,*) alpha(:)
            write(VrelFile,*) Vmag(:)
            write(VaxialFile,*) windVectors(:,1)
            write(VtangentialFile,*) windVectors(:,2)
    
        ! Write blade points 
!        call atm_write_blade_points(i,jt_total)

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
    write(*,*) '!  Done Writing Actuator Airfoil Model output  !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

endif

end subroutine atm_output

!-------------------------------------------------------------------------------
function atm_convoluteForce(i,xyz)
!-------------------------------------------------------------------------------
! This subroutine will convolute the body forces onto a point xyz
integer, intent(in) :: i
! i - number of sector
real(rprec), intent(in) :: xyz(3)    ! Point onto which to convloute the force 
real(rprec) :: Force(3)   ! The blade force to be convoluted
real(rprec) :: dis                ! Distance onto which convolute the force
real(rprec) :: atm_convoluteForce(3)    ! The local velocity at this point
real(rprec) :: kernel                ! Gaussian dsitribution value

! Distance from the point of the force to the point where it is being convoluted
dis=distance(xyz,bladePoints(i,:))

! The force which is being convoluted
Force= bladeForces(i,:)

! The value of the kernel. This is the actual smoothing function
kernel=exp(-(dis/epsilon)**2._rprec) / &
((epsilon**3._rprec)*(pi**1.5_rprec))

! The force times the kernel will give the force/unitVolume
atm_convoluteForce = Force * kernel

return
end function atm_convoluteForce



end module actuator_turbine_model

















