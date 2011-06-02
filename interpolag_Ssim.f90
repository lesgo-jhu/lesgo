subroutine interpolag_Ssim ()
! This subroutine takes the arrays F_LM and F_MM from the previous  
!   timestep and essentially moves the values around to follow the 
!   corresponding particles. The (x, y, z) value at the current 
!   timestep will be the (x-u*dt, y-v*dt, z-w*dt) value at the 
!   previous timestep.  Since particle motion does not conform to
!   the grid, an interpolation will be required.  Variables should 
!   be on the w-grid.

use types,only:rprec
use param
use sgsmodule
use messages
use grid_defs,only:grid_t 
use functions, only:trilinear_interp
$if ($MPI)
use mpi_defs, only:mpi_sync_real_array,MPI_SYNC_DOWNUP
$endif
implicit none

character (*), parameter :: sub_name = 'interpolag_Ssim'
$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rprec), dimension(3) :: xyz_past
real(rprec), dimension(ld,ny,$lbz:nz) :: tempF_LM, tempF_MM
$if ($DYN_TN)
real(rprec), dimension(ld,ny,$lbz:nz) :: tempF_ee2, tempF_deedt2, tempee_past
$endif
integer :: i,j,k,kmin,jz

real (rprec) :: lcfl

real(rprec), pointer, dimension(:) :: x,y,z

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

nullify(x,y,z)
x => grid_t % x
y => grid_t % y
z => grid_t % z

! Perform (backwards) Lagrangian interpolation

    ! Create dummy arrays so information will not be overwritten during interpolation
        tempF_LM(:,:,1:nz) = F_LM(:,:,1:nz)
        tempF_MM(:,:,1:nz) = F_MM(:,:,1:nz)  
        $if ($DYN_TN)
        tempF_ee2 = F_ee2
        tempF_deedt2 = F_deedt2
        tempee_past = ee_past  
        $endif 
    
        ! Loop over domain (within proc): for each, calc xyz_past then trilinear_interp
        !   Note: we need to put u_lag and v_lag on w-nodes (same as Cs, F_LM, F_MM, etc)
        !  Also, do not allow interpolation out of top/bottom of domain.
        
        ! Bottom-most level for Cs, F_LM, F_MM, etc is on uvp-nodes
        !   This requires special treatment for k=1,2 since F_*(i,j,k=1) is closer 
        !   to F_*(i,j,k=2) than the usual spacing, dz.  
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord.eq.0)) then
                kmin = 3 
                do k=1,2                     
                do j=1,ny   
                do i=1,nx
                    ! Determine position at previous timestep
                    if (k.eq.1) then    ! UVP-node; If w is positive, set w=0            
                        xyz_past(1) = x(i) - u_lag(i,j,k)*dt
                        xyz_past(2) = y(j) - v_lag(i,j,k)*dt
                        xyz_past(3) = z(k) - 0.5_rprec*min(w_lag(i,j,k+1),0._rprec)*dt                  
                    else                ! W-node; If w is positive, multiply by two
                        xyz_past(1) = x(i) - 0.5*(u_lag(i,j,k)+u_lag(i,j,k-1))*dt
                        xyz_past(2) = y(j) - 0.5*(v_lag(i,j,k)+v_lag(i,j,k-1))*dt                    
                        xyz_past(3) = z(k) - (w_lag(i,j,k) + max(w_lag(i,j,k),0._rprec))*dt 
                    endif
               
                    ! Interpolate
                    F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    $if ($DYN_TN)
                    F_ee2(i,j,k) = trilinear_interp(tempF_ee2(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_deedt2(i,j,k) = trilinear_interp(tempF_deedt2(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    ee_past(i,j,k) = trilinear_interp(tempee_past(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    $endif 
                enddo
                enddo
                enddo
            else
                kmin = 1
            endif
        ! Intermediate levels on w-nodes
            do k=kmin,nz-1
            do j=1,ny
            do i=1,nx
                ! Determine position at previous timestep
                xyz_past(1) = x(i) - 0.5*(u_lag(i,j,k)+u_lag(i,j,k-1))*dt
                xyz_past(2) = y(j) - 0.5*(v_lag(i,j,k)+v_lag(i,j,k-1))*dt
                xyz_past(3) = z(k) - w_lag(i,j,k)*dt  ! minus dz/2 (since these values are on w-grid)
                                             ! but this dz/2 would be added back during interpolation                     

                ! Interpolate   
                F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
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
                    xyz_past(1) = x(i) - 0.5*(u_lag(i,j,k)+u_lag(i,j,k-1))*dt
                    xyz_past(2) = y(j) - 0.5*(v_lag(i,j,k)+v_lag(i,j,k-1))*dt
                                              
                    if (w_lag(i,j,k).le.0) then    ! trilinear_interp won't work (boo)
                        ! this is cheating, but it should work for now...
                        xyz_past(3) = z(k) - 1.0e-10_rprec*dt
                    else
                        xyz_past(3) = z(k) - w_lag(i,j,k)*dt
                    endif
                    
                    ! Interpolate
                    F_LM(i,j,k) = trilinear_interp(tempF_LM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
                    F_MM(i,j,k) = trilinear_interp(tempF_MM(1:nx,1:ny,$lbz:nz),$lbz,xyz_past)
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
            $if ($DYN_TN)
            call mpi_sync_real_array( F_ee2, MPI_SYNC_DOWNUP )
            call mpi_sync_real_array( F_deedt2, MPI_SYNC_DOWNUP )
            call mpi_sync_real_array( ee_past, MPI_SYNC_DOWNUP )
            $endif 
        $endif   

    $if ($DEBUG)
    if (DEBUG) write (*, *) 'interpolag_Ssim: after interpolation'
    $endif

! Compute the Lagrangian CFL number and print to screen
!   Note: this is only in the x-direction... not good for complex geometry cases
    if (mod (jt, cfl_count) .eq. 0) then
        lcfl = 0._rprec
        do jz = 1, nz
            lcfl = max ( lcfl,  maxval (abs (u_lag(1:nx, :, jz))) )*dt/dx
        enddo
        print*, 'Lagrangian CFL condition= ', lcfl
    endif

! Reset Lagrangian u/v/w variables for use during next set of cs_counts
    u_lag = 0._rprec;    v_lag = 0._rprec;    w_lag = 0._rprec

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

nullify(x,y,z)

end subroutine interpolag_Ssim
