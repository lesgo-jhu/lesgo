!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

subroutine interpolag_Ssim ()
! This subroutine takes the arrays F_LM and F_MM from the previous  
!   timestep and essentially moves the values around to follow the 
!   corresponding particles. The (x, y, z) value at the current 
!   timestep will be the (x-u*dt, y-v*dt, z-w*dt) value at the 
!   previous timestep.  Since particle motion does not conform to
!   the grid, an interpolation will be required.  Variables should 
!   be on the w-grid.

! This subroutine assumes that dt and cs_count are chosen such that
!   the Lagrangian CFL in the z-direction will never exceed 1.  If the
!   Lag. CFL in the x-direction is less than one this should generally
!   be satisfied.

use types,only:rprec
use param
use sgs_param, only: F_LM, F_MM, lagran_dt
!$if ($DYN_TN)
!use sgs_param, only: F_ee2, F_deedt2, ee_past
!$endif
use messages
use sim_param,only:u,v,w
use grid_defs,only:grid 
use functions, only:trilinear_interp
$if ($MPI)
use mpi_defs, only:mpi_sync_real_array,MPI_SYNC_DOWNUP
$endif
use cfl_util, only : get_max_cfl
implicit none

$if ($VERBOSE)
character (*), parameter :: sub_name = 'interpolag_Ssim'
$endif
$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

real(rprec), dimension(3) :: xyz_past
real(rprec), dimension(ld,ny,lbz:nz) :: tempF_LM, tempF_MM
!$if ($DYN_TN)
!real(rprec), dimension(ld,ny,lbz:nz) :: tempF_ee2, tempF_deedt2, tempee_past
!$endif
integer :: i,j,k,kmin

real (rprec) :: lcfl

real(rprec), pointer, dimension(:) :: x,y,z

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

! Perform (backwards) Lagrangian interpolation
    ! F_* arrays should be synced at this point (for MPI)

    ! Create dummy arrays so information will not be overwritten during interpolation
        tempF_LM = F_LM
        tempF_MM = F_MM
!        $if ($DYN_TN)
!        tempF_ee2 = F_ee2
!        tempF_deedt2 = F_deedt2
!        tempee_past = ee_past  
!        $endif 
    
        ! Loop over domain (within proc): for each, calc xyz_past then trilinear_interp
        ! Variables x,y,z, F_LM, F_MM, etc are on w-grid
        ! Interpolation out of top/bottom of domain is not permitted.
        ! Note: x,y,z values are only good for k=1:nz-1 within each proc
            if ( coord.eq.0 ) then
                kmin = 2 
                ! At the bottom-most level (at the wall) the velocities are zero.
                ! Since there is no movement the values of F_LM, F_MM, etc should
                !   not change and no interpolation is necessary.
            else
                kmin = 1
            endif
        ! Intermediate levels
            do k=kmin,nz-1
            do j=1,ny
            do i=1,nx
                ! Determine position at previous timestep (u,v interp to w-grid) 
                xyz_past(1) = x(i) - 0.5_rprec*(u(i,j,k-1)+u(i,j,k))*lagran_dt
                xyz_past(2) = y(j) - 0.5_rprec*(v(i,j,k-1)+v(i,j,k))*lagran_dt
                xyz_past(3) = z(k) - w(i,j,k)*lagran_dt
                
                ! Interpolate   
                F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,lbz:nz),lbz,xyz_past)
                F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                $if ($DYN_TN)
!                F_ee2(i,j,k) = trilinear_interp(tempF_ee2(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                F_deedt2(i,j,k) = trilinear_interp(tempF_deedt2(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                ee_past(i,j,k) = trilinear_interp(tempee_past(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                $endif 
            enddo
            enddo
            enddo               
        ! Top-most level should not allow negative w
            $if ($MPI)
            if (coord.eq.nproc-1) then
            $endif
                k = nz
                do j=1,ny
                do i=1,nx
                    ! Determine position at previous timestep (u,v interp to w-grid)
                    xyz_past(1) = x(i) - 0.5_rprec*(u(i,j,k-1)+u(i,j,k))*lagran_dt
                    xyz_past(2) = y(j) - 0.5_rprec*(v(i,j,k-1)+v(i,j,k))*lagran_dt  
                    xyz_past(3) = z(k) - max(0.0_rprec,w(i,j,k))*lagran_dt
                    
                    ! Interpolate
                    F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,lbz:nz),lbz,xyz_past)
                    F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                    $if ($DYN_TN)
!                    F_ee2(i,j,k) = trilinear_interp(tempF_ee2(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                    F_deedt2(i,j,k) = trilinear_interp(tempF_deedt2(1:nx,1:ny,lbz:nz),lbz,xyz_past)
!                    ee_past(i,j,k) = trilinear_interp(tempee_past(1:nx,1:ny,lbz:nz),lbz,xyz_past)    
!                    $endif
                enddo
                enddo    
            $if ($MPI)
            endif    
            $endif
            
         ! Share new data between overlapping nodes
         $if ($MPI)
            call mpi_sync_real_array( F_LM, 0, MPI_SYNC_DOWNUP )  
            call mpi_sync_real_array( F_MM, 0, MPI_SYNC_DOWNUP )   
!            $if ($DYN_TN)
!            call mpi_sync_real_array( F_ee2, 0, MPI_SYNC_DOWNUP )
!            call mpi_sync_real_array( F_deedt2, 0, MPI_SYNC_DOWNUP )
!            call mpi_sync_real_array( ee_past, 0, MPI_SYNC_DOWNUP )
!            $endif 
        $endif   

    $if ($DEBUG)
    if (DEBUG) write (*, *) 'interpolag_Ssim: after interpolation'
    $endif

! Compute the Lagrangian CFL number and print to screen
!   Note: this is only in the x-direction... not good for complex geometry cases
    if (mod (jt_total, lag_cfl_count) .eq. 0) then
        lcfl = get_max_cfl()
        lcfl = lcfl*lagran_dt/dt  
        $if($MPI)
            if(coord.eq.0) print*, 'Lagrangian CFL condition= ', lcfl
        $else
            print*, 'Lagrangian CFL condition= ', lcfl
        $endif
    endif
    
$if ($VERBOSE)
call exit_sub (sub_name)
$endif

nullify(x,y,z)

end subroutine interpolag_Ssim
