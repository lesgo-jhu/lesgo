subroutine interpolag_Sdep()
! This subroutine takes the arrays F_{LM,MM,QN,NN} from the previous  
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
use sgsmodule

use grid_defs,only:grid_t 
use functions, only:trilinear_interp
$if ($MPI)
use mpi_defs, only:mpi_sync_real_array,MPI_SYNC_DOWNUP
$endif
implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rprec), dimension(3) :: xyz_past
real(rprec), dimension(ld,ny,$lbz:nz) :: tempF_LM, tempF_MM, tempF_QN, tempF_NN
$if ($DYN_TN)
real(rprec), dimension(ld,ny,$lbz:nz) :: tempF_ee2, tempF_deedt2, tempee_past
$endif
integer :: i,j,k,kmin,jz

real (rprec) :: lcfl

real(rprec), pointer, dimension(:) :: x,y,z

real(kind=rprec), save, dimension(nx+2,ny+2,nz+2) :: Beta_t

!---------------------------------------------------------------------
$if ($VERBOSE)
write (*, *) 'started interpolag_Sdep'
$endif

nullify(x,y,z)
x => grid_t % x
y => grid_t % y
z => grid_t % z

! Perform (backwards) Lagrangian interpolation

    ! Create dummy arrays so information will not be overwritten during interpolation
        tempF_LM(:,:,1:nz) = F_LM(:,:,1:nz)
        tempF_MM(:,:,1:nz) = F_MM(:,:,1:nz)  
        tempF_QN(:,:,1:nz) = F_QN(:,:,1:nz)
        tempF_NN(:,:,1:nz) = F_NN(:,:,1:nz)          
        $if ($DYN_TN)
        tempF_ee2 = F_ee2
        tempF_deedt2 = F_deedt2
        tempee_past = ee_past  
        $endif 

        ! Loop over domain (within proc): for each, calc xyz_past then trilinear_interp
        ! Variables x,y,z_lag, F_LM, F_MM, F_QN, F_NN, etc are on w-grid
        ! Interpolation out of top/bottom of domain is not permitted.
        ! Note: x,y,z_lag values are only good for k=1:nz-1 within each proc
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord.eq.0)) then
                kmin = 2                    
                ! At the bottom-most level (at the wall) the velocities are zero.
                ! Since there is no movement the values of F_LM, F_MM, etc should
                !   not change and no interpolation is necessary.
                x_lag(:,:,1) = 0.0_rprec
                y_lag(:,:,1) = 0.0_rprec
                z_lag(:,:,1) = 0.0_rprec                
            else
                kmin = 1
            endif
        ! Intermediate levels
            do k=kmin,nz-1
            do j=1,ny
            do i=1,nx
                ! Determine position at previous timestep
                xyz_past(1) = x(i) - x_lag(i,j,k)
                xyz_past(2) = y(j) - y_lag(i,j,k)
                xyz_past(3) = z(k) - z_lag(i,j,k)                  

                ! Interpolate   
                F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                F_QN(i,j,k) = trilinear_interp(tempF_QN(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                F_NN(i,j,k) = trilinear_interp(tempF_NN(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)                          
                $if ($DYN_TN)
                F_ee2(i,j,k) = trilinear_interp(tempF_ee2(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                F_deedt2(i,j,k) = trilinear_interp(tempF_deedt2(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                ee_past(i,j,k) = trilinear_interp(tempee_past(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                $endif 
            enddo
            enddo
            enddo               
        ! Top-most level should not allow negative w
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord.eq.nproc-1)) then
                k = nz
                do j=1,ny
                do i=1,nx
                    ! Determine position at previous timestep
                    xyz_past(1) = x(i) - x_lag(i,j,k)
                    xyz_past(2) = y(j) - y_lag(i,j,k)   
                                              
                    if (z_lag(i,j,k).le.0) then    ! trilinear_interp won't work (boo)
                        ! this is cheating, but it should work for now...
                        xyz_past(3) = z(k) - 1.0e-10_rprec
                    else
                        xyz_past(3) = z(k) - z_lag(i,j,k)
                    endif
                    
                    ! Interpolate
                    F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_QN(i,j,k) = trilinear_interp(tempF_QN(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_NN(i,j,k) = trilinear_interp(tempF_NN(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)                      
                    $if ($DYN_TN)
                    F_ee2(i,j,k) = trilinear_interp(tempF_ee2(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_deedt2(i,j,k) = trilinear_interp(tempF_deedt2(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    ee_past(i,j,k) = trilinear_interp(tempee_past(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)    
                    $endif
                enddo
                enddo    
            endif          
        
         ! Share new data between overlapping nodes
         $if ($MPI)
            call mpi_sync_real_array( F_LM, MPI_SYNC_DOWNUP )  
            call mpi_sync_real_array( F_MM, MPI_SYNC_DOWNUP )   
            call mpi_sync_real_array( F_QN, MPI_SYNC_DOWNUP )  
            call mpi_sync_real_array( F_NN, MPI_SYNC_DOWNUP )              
            $if ($DYN_TN)
            call mpi_sync_real_array( F_ee2, MPI_SYNC_DOWNUP )
            call mpi_sync_real_array( F_deedt2, MPI_SYNC_DOWNUP )
            call mpi_sync_real_array( ee_past, MPI_SYNC_DOWNUP )
            $endif 
        $endif   

! Compute the Lagrangian CFL number and print to screen
!   Note: this is only in the x-direction... not good for complex geometry cases
    if (mod (jt, cfl_count) .eq. 0) then
        lcfl = 0._rprec
        do jz = 1, nz  
            lcfl = max ( lcfl,  maxval (abs (x_lag(1:nx, :, jz)))/dx )
        enddo
        print*, 'Lagrangian CFL condition= ', lcfl
    endif

! Reset Lagrangian u/v/w variables for use during next set of cs_counts
    x_lag = 0._rprec;    y_lag = 0._rprec;    z_lag = 0._rprec  
    lagran_dt = 0._rprec
        
$if ($VERBOSE)
write (*, *) 'finished interpolag_Sdep'
$endif

nullify(x,y,z)

end subroutine interpolag_Sdep
