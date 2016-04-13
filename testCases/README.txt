!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  
!
!  This directory contains a set of test-cases for lesgo
!  The test cases are defined below
!
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            actuator_turbine_model                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The case has 2 aligned wind turbines under uniform inlow
! The cases uses the fringe region to create a uniform inflow velocity
!
! executable - lesgo-mpi-ATM
! The input for this case is in the directory './inputATM'
!
! ./inputATM/turbineArrayProperties - a file specifying the location and 
! properties of each turbine
!
! NREL5MW - A reference file with all the specs for the turbine type  (In this 
! case it is the NREL 5 MW)
!
! ./inputATM/Airfoils - This directory contains the files with 
! lift and drag coefficients

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               LES_Lagrangian                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A standard Large-Eddy Simulation of an Atmospheric Boundary Layer
! executable - lesgo-mpi
!
! The Number of grid points is 32 x 32 x 32
! The dimensions are 2 pi x 2 pi x 1
! The bottom boundary condition is an equilibrium wall model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         LES_Integral_Wall_Model                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A standard Large-Eddy Simulation of an Atmospheric Boundary Layer
! executable - lesgo-mpi
!
! The Number of grid points is 32 x 32 x 32
! The dimensions are 2 pi x 2 pi x 1
! The bottom boundary condition is an integral wall model (Xiang Yang)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                  precursor                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a concurrent-precursor simulation
! There are 2 domains
! red - this is used to develop the inflow in a perdiodic channel
! blue -  this domain contains the simulation of interest with a devloping ABL
! The case setup is for a single wind farm from the Chamorro experiment
! Reference: Wind Turbine Large-Eddy Simulations on Very Coarse Grid Resolutions using an Actuator Line Model
!            LAM Tossas, RJAM Stevens, C Meneveau 2016
! This is commonly used for simulating finite-wind farms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                turbines (ADM)                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A perdiodic simulation with a wind farm with an Actuator Disk Model
! The thrust coefficient and wind farm characteristics are specified in 
! lesgo.conf
