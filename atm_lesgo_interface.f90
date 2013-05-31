!*******************************************************************************
module atm_lesgo_interface
!*******************************************************************************
! This module interfaces actuator turbine module with lesgo
! It is a lesgo specific module, unlike the atm module

! Remember to always dimensinoalize the variables from LESGO
! Length is non-dimensionalized by z_i

! Lesgo data used regarding the grid
use param, only : dt ,nx,ny,nz,jzmin,jzmax,dx,dy,dz,coord,lbz,nproc, z_i 
! nx, ny, nz - nodes in every direction
! z_i - non-dimensionalizing length
! dt - time-step 

use grid_defs, only: grid 

! Actuator Turbine Model module
use actuator_turbine_model
use atm_input_util, only : rprec, turbineArray

implicit none

private
public atm_lesgo_initialize

! This is a dynamic list that stores all the points in the domain with a body 
! force due to the turbines. 
type bodyForce_t
    ! c0 counter for all the points in the domain with forces
    ! i,j,k stores the index for the point in the domain
    integer c0,i,j,k
    real(rprec), dimension(3) :: force ! Force vector
    type(bodyForce_t), pointer :: next
    type(bodyForce_t), pointer :: previous
end type bodyForce_t



! This is the body force vector which has the forces by the atm
type(bodyForce_t), allocatable, dimension(:) :: bodyForce

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

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_initialize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model

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
                end if
            end do 
        end do
    end do
end do
write(*,*) c0

end subroutine atm_lesgo_initialize


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_update ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
! Calls function within the actuator_turbine_model module

call atm_forcing ()

end subroutine atm_lesgo_update


!-------------------------------------------------------------------------------
function distance(a,b)
! This function calculates the distance between (a,b,c) and (d,e,f)
!-------------------------------------------------------------------------------
real(rprec), dimension(3), intent(in) :: a,b
real(rprec) :: distance
distance=sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
return
end function distance

end module atm_lesgo_interface


