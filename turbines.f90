!!
!!  Copyright (C) 2010-2016  Johns Hopkins University
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

module turbines
! This module contains all of the subroutines associated with drag-disk turbines:

use param
use turbines_base
use grid_m
use messages
use string_util

implicit none

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize

real(rprec) :: T_avg_dim_file
real(rprec), dimension(:), allocatable :: z_tot

character (100) :: string1
real(rprec) :: eps !epsilon used for disk velocity time-averaging

integer :: i,j,k,i2,j2,k2,i3,j3,i4,j4,b,l,s,ssx,ssy,ssz, p
integer :: imax,jmax,kmax,count_i,count_n,icp,jcp,kcp
integer :: min_i,max_i,min_j,max_j,min_k,max_k,cut
integer :: k_start, k_end
logical :: exst, opn

logical :: turbine_in_proc=.false.      !init, do not change this

real(rprec), pointer, dimension(:) :: buffer_array
real(rprec) :: buffer
logical :: buffer_logical
integer, dimension(:), allocatable :: turbine_in_proc_array
integer :: turbine_in_proc_cnt = 0
integer, dimension(:), allocatable :: file_id,file_id2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_init()
!  This subroutine creates the 'turbine' folder and starts the turbine forcing output files.
!  It also creates the indicator function (Gaussian-filtered from binary locations - in or out)
!  and sets values for turbine type (node locations, etc)
use open_file_fid_mod
implicit none

real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.turbines_init'
integer :: fid
real(rprec) :: delta2

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

! Allocate and initialize
allocate(turbine_in_proc_array(nproc-1))
allocate(z_tot(nz_tot))
allocate(file_id(nloc))
allocate(file_id2(nloc))
turbine_in_proc_array = 0

nullify(buffer_array)
allocate(buffer_array(nloc))

!Create turbine directory
call system("mkdir -vp turbine") 

!z_tot for total domain (since z is local to the processor)
do k=1,nz_tot
    z_tot(k) = (k - 0.5_rprec) * dz
enddo

!Compute a lookup table object for the indicator function 
delta2 = alpha**2 * (dx**2 + dy**2 + dz**2)
call turb_ind_func%init(delta2, thk_all, dia_all, max( max(nx, ny), nz) )

!find turbine nodes - including filtered ind, n_hat, num_nodes, and nodes for each turbine
!each processor finds turbines in its domain
call turbines_nodes

if (turbine_cumulative_time) then
    if (coord == 0) then
        string1 = path // 'turbine/turbine_u_d_T.dat'
        inquire (file=string1, exist=exst)
        if (exst) then
            write(*,*) 'Reading from file turbine_u_d_T.dat'
            fid = open_file_fid( string1, 'rewind', 'formatted' )
            do i=1,nloc
                read(fid,*) wind_farm%turbine(i)%u_d_T    
            enddo    
            read(fid,*) T_avg_dim_file
            if (T_avg_dim_file /= T_avg_dim) then
                write(*,*) 'Time-averaging window does not match value in turbine_u_d_T.dat'
            endif
            close (fid)
        else  
            write (*, *) 'File ', trim(string1), ' not found'
            write (*, *) 'Assuming u_d_T = -1. for all turbines'
            do k=1,nloc
                wind_farm%turbine(k)%u_d_T = -1.
            enddo
        endif                                    
    endif
else
    write (*, *) 'Assuming u_d_T = -1 for all turbines'
    do k=1,nloc
        wind_farm%turbine(k)%u_d_T = -1.
    enddo    
endif

! Set all Ct_prime to reference in input file
do k=1,nloc
    wind_farm%turbine(k)%Ct_prime = Ct_prime
enddo 
   
if (coord .eq. nproc-1) then
    string1=path // 'output/vel_top_of_domain.dat'
    fid = open_file_fid( string1, 'rewind', 'formatted' )
    write(fid,*) 'total_time','u_HI'
    close(fid)
endif

! Generate the files for the turbine forcing output
do s=1,nloc
   if(coord==0) then
   call string_splice( string1, path // 'turbine/turbine_', s, '_forcing.dat' )
   file_id(s) = open_file_fid( string1, 'append', 'formatted' )
   endif
enddo

! Generate the files for the turbine velocity output
do s=1,nloc
#ifdef PPMPI
   kcp = nint(wind_farm%turbine(s)%height/dz + 0.5)
   k_start =  1+coord*(nz-1)
   k_end = nz-1+coord*(nz-1)

   if (kcp>=k_start .and. kcp<=k_end) then
       call string_splice( string1, path // 'turbine/turbine_', s, '_velcenter.dat' )
       call string_concat (string1, '.c', coord)
       file_id2(s) = open_file_fid( string1, 'append', 'formatted' )
   endif
#else
   call string_splice( string1, path // 'turbine/turbine_', s, '_velcenter.dat' )
   file_id2(s) = open_file_fid( string1, 'append', 'formatted' )
#endif
enddo

nullify(x,y,z)

end subroutine turbines_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_nodes
!This subroutine locates nodes for each turbine and builds the arrays: ind, n_hat, num_nodes, and nodes
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_nodes'

real(rprec) :: R_t,rx,ry,rz,r,r_norm,r_disk

real(rprec), pointer :: p_xloc => null(), p_yloc=> null(), p_height=> null()
real(rprec), pointer :: p_dia => null(), p_thk=> null(), p_theta1=> null(), p_theta2=> null()
real(rprec), pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3=> null() 

real(rprec), pointer, dimension(:) :: x, y, z

real(rprec) :: filt
real(rprec), dimension(:), allocatable :: sumA, turbine_vol

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

turbine_in_proc = .false.

allocate(sumA(nloc))
allocate(turbine_vol(nloc))
sumA = 0

do s=1,nloc
    
    count_n = 0    !used for counting nodes for each turbine
    count_i = 1    !index count - used for writing to array "nodes"

    !set pointers
    p_xloc => wind_farm%turbine(s)%xloc     
    p_yloc => wind_farm%turbine(s)%yloc  
    p_height => wind_farm%turbine(s)%height 
    p_dia => wind_farm%turbine(s)%dia 
    p_thk => wind_farm%turbine(s)%thk
    p_theta1 => wind_farm%turbine(s)%theta1
    p_theta2 => wind_farm%turbine(s)%theta2
    p_nhat1 => wind_farm%turbine(s)%nhat(1)
    p_nhat2 => wind_farm%turbine(s)%nhat(2)
    p_nhat3 => wind_farm%turbine(s)%nhat(3)

    !identify "search area"
    imax = p_dia/dx + 2
    jmax = p_dia/dy + 2
    kmax = p_dia/dz + 2

    !determine unit normal vector for each turbine	
    p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat3 = sin(pi*p_theta2/180.)

    !determine nearest (i,j,k) to turbine center
    icp = nint(p_xloc/dx)
    jcp = nint(p_yloc/dy)
    kcp = nint(p_height/dz + 0.5)

    !determine limits for checking i,j,k
    !due to spectral BCs, i and j may be < 1 or > nx,ny
    !the mod function accounts for this when these values are used
    min_i = icp-imax
    max_i = icp+imax
    min_j = jcp-jmax
    max_j = jcp+jmax
    min_k = max((kcp-kmax),1)
    max_k = min((kcp+kmax),nz_tot)
    wind_farm%turbine(s)%nodes_max(1) = min_i
    wind_farm%turbine(s)%nodes_max(2) = max_i
    wind_farm%turbine(s)%nodes_max(3) = min_j
    wind_farm%turbine(s)%nodes_max(4) = max_j
    wind_farm%turbine(s)%nodes_max(5) = min_k
    wind_farm%turbine(s)%nodes_max(6) = max_k      

    !check neighboring grid points	
    !update num_nodes, nodes, and ind for this turbine
    !and split domain between processors
    !z(nz) and z(1) of neighboring coords match so each coord gets (local) 1 to nz-1
    wind_farm%turbine(s)%ind = 0.
    wind_farm%turbine(s)%nodes = 0.
    wind_farm%turbine(s)%num_nodes = 0.
    count_n = 0
    count_i = 1
#ifdef PPMPI
    k_start = max(1+coord*(nz-1), min_k)
    k_end = min(nz-1+coord*(nz-1), max_k)
#else
    k_start = 1
    k_end = nz
#endif
    
    do k=k_start,k_end  !global k     
        do j=min_j,max_j
            do i=min_i,max_i
                !vector from center point to this node is (rx,ry,rz) with length r
                if (i<1) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (x(i2)-L_x) - p_xloc
                elseif (i>nx) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (L_x+x(i2)) - p_xloc
                else
                    i2 = i
                    rx = x(i) - p_xloc 
                endif            
                if (j<1) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (y(j2)-L_y) - p_yloc                
                elseif (j>ny) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (L_y+y(j2)) - p_yloc
                else
                    j2 = j
                    ry = y(j) - p_yloc 
                endif                      
                rz = z_tot(k) - p_height 
                r = sqrt(rx*rx + ry*ry + rz*rz)
                !length projected onto unit normal for this turbine
                r_norm = abs(rx*p_nhat1 + ry*p_nhat2 + rz*p_nhat3)
                !(remaining) length projected onto turbine disk
                r_disk = sqrt(r*r - r_norm*r_norm)
                ! get the filter value
                filt = turb_ind_func%val(r_disk, r_norm)
                if ( filt > filter_cutoff ) then
                    wind_farm%turbine(s)%ind(count_i) = filt
                    wind_farm%turbine(s)%nodes(count_i,1) = i2
                    wind_farm%turbine(s)%nodes(count_i,2) = j2
                    wind_farm%turbine(s)%nodes(count_i,3) = k - coord*(nz-1) !local k
                    count_n = count_n + 1
                    count_i = count_i + 1
                    turbine_in_proc = .true.
                    sumA(s) = sumA(s) + filt * dx * dy * dz
                endif
           enddo
       enddo
    enddo
    wind_farm%turbine(s)%num_nodes = count_n
    
    ! Calculate turbine volume
    turbine_vol(s) = pi/4. * (wind_farm%turbine(s)%dia)**2 * wind_farm%turbine(s)%thk
    
enddo

! Sum the indicator function across all processors if using MPI
#ifdef PPMPI 
buffer_array = sumA
call MPI_Allreduce(buffer_array, sumA, nloc, MPI_rprec, MPI_SUM, comm, ierr)
#endif

! Normalize the indicator function
do s = 1, nloc
    wind_farm%turbine(s)%ind = wind_farm%turbine(s)%ind(:) * turbine_vol(s) / sumA(s)
end do

!each processor sends its value of turbine_in_proc
!if false, disk-avg velocity will not be sent (since it will always be 0.)
#ifdef PPMPI
turbine_in_proc_cnt = 0
if (coord == 0) then
    do i=1,nproc-1
        call MPI_recv(buffer_logical, 1, MPI_logical, i, 2, comm, status, ierr )
        if (buffer_logical) then
            turbine_in_proc_cnt = turbine_in_proc_cnt + 1
            turbine_in_proc_array(turbine_in_proc_cnt) = i
        endif
    enddo
else
    call MPI_send(turbine_in_proc, 1, MPI_logical, 0, 2, comm, ierr )
endif
#endif

! Cleanup
deallocate(sumA)
deallocate(turbine_vol)
nullify(x,y,z)

end subroutine turbines_nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_forcing()
! This subroutine applies the drag-disk forcing
use sim_param, only : u,v,w, fxa,fya,fza
use functions, only : linear_interp, interp_to_uv_grid
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_forcing'

real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec), pointer :: p_Ct_prime => null()

real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vels, disk_force
real(rprec), allocatable, dimension(:,:,:) :: w_uv ! Richard: This 3D matrix can relatively easy be prevented
real(rprec), pointer, dimension(:) :: y,z

nullify(y,z)
y => grid % y
z => grid % z

allocate(w_uv(ld,ny,lbz:nz))

#ifdef PPMPI
call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)     !syncing intermediate w-velocities!
#endif

w_uv = interp_to_uv_grid(w, lbz)

! Do interpolation for dynamically changing parameters
do s = 1, nloc
    if (dyn_theta1) wind_farm%turbine(s)%theta1 =                       &
        linear_interp(theta1_time, theta1_arr(s,:), total_time_dim)
    if (dyn_theta2) wind_farm%turbine(s)%theta2 =                       &
        linear_interp(theta2_time, theta2_arr(s,:), total_time_dim)
    if (dyn_Ct_prime) wind_farm%turbine(s)%Ct_prime =                   &
        linear_interp(Ct_prime_time, Ct_prime_arr(s,:), total_time_dim)        
end do

! Recompute the turbine position if theta1 or theta2 can change
if (dyn_theta1 .or. dyn_theta2) call turbines_nodes

!Each processor calculates the weighted disk-averaged velocity
disk_avg_vels = 0.
if (turbine_in_proc) then

    !for each turbine:        
    do s=1,nloc      
         
        !set pointers
        p_u_d => wind_farm%turbine(s)%u_d   

        !calculate total disk-averaged velocity for each turbine (current,instantaneous)    
        !u_d is the velocity in the normal direction	  
        !weighted average using "ind"            
        p_u_d = 0.
        do l=1,wind_farm%turbine(s)%num_nodes   
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            p_u_d = p_u_d + (wind_farm%turbine(s)%nhat(1)*u(i2,j2,k2) &
                           + wind_farm%turbine(s)%nhat(2)*v(i2,j2,k2) &
                           + wind_farm%turbine(s)%nhat(3)*w_uv(i2,j2,k2)) &
                * wind_farm%turbine(s)%ind(l)        
         enddo
         
        ! write center disk velocity to file                   
        if (modulo (jt_total, tbase) == 0) then
            icp = nint(wind_farm%turbine(s)%xloc/dx)
            jcp = nint(wind_farm%turbine(s)%yloc/dy)
            kcp = nint(wind_farm%turbine(s)%height/dz + 0.5)
#ifdef PPMPI
            k_start =  1+coord*(nz-1)
            k_end = nz-1+coord*(nz-1)
#else
            k_start = 1
            k_end = nz
#endif
            if (kcp>=k_start .and. kcp<=k_end) then
                kcp=kcp-k_start+1
                write( file_id2(s), *) total_time_dim, u(icp,jcp,kcp), w_uv(icp,jcp,kcp)
            endif
        endif
        
        !write this value to the array (which will be sent to coord 0)
        disk_avg_vels(s) = p_u_d
        
    enddo
        
endif        

!send the disk-avg values to coord==0
#ifdef PPMPI 
call mpi_barrier (comm,ierr)

if (coord == 0) then
    do i=1,turbine_in_proc_cnt
        j = turbine_in_proc_array(i)
        buffer_array = 0.
        call MPI_recv( buffer_array, nloc, MPI_rprec, j, 3, comm, status, ierr )
        disk_avg_vels = disk_avg_vels + buffer_array
    enddo
elseif (turbine_in_proc) then
    call MPI_send( disk_avg_vels, nloc, MPI_rprec, 0, 3, comm, ierr )
endif                          
#endif

!Coord==0 takes that info and calculates total disk force, then sends it back
if (coord == 0) then           
    !update epsilon for the new timestep (for cfl_dt)
    if (T_avg_dim > 0.) then
        eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
    else
        eps = 1.
    endif

    !for each turbine:        
    do s=1,nloc            
         
        !set pointers
        p_u_d => wind_farm%turbine(s)%u_d   
        p_u_d_T => wind_farm%turbine(s)%u_d_T   
        p_f_n => wind_farm%turbine(s)%f_n                
        p_Ct_prime => wind_farm%turbine(s)%Ct_prime  
        
        !volume correction:
        !since sum of ind is turbine volume/(dx*dy*dz) (not exactly 1.)
        p_u_d = disk_avg_vels(s) * wind_farm%turbine(s)%vol_c

        !add this current value to the "running average" (first order relaxation)
        p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

        !calculate total thrust force for each turbine  (per unit mass)
        !force is normal to the surface (calc from u_d_T, normal to surface)
        p_f_n = -0.5*p_Ct_prime*abs(p_u_d_T)*p_u_d_T/wind_farm%turbine(s)%thk       
        !write values to file                   
            if (modulo (jt_total, tbase) == 0) then
               write( file_id(s), *) total_time_dim, p_u_d, p_u_d_T, p_f_n, &
                   -0.5*p_Ct_prime*(p_u_d_T*p_u_d_T*p_u_d_T)*pi/(4.*sx*sy)
            endif 
        !write force to array that will be transferred via MPI    
        disk_force(s) = p_f_n
        
    enddo                   
endif

!send total disk force to the necessary procs (with turbine_in_proc==.true.)
#ifdef PPMPI
if (coord == 0) then  
    do i=1,turbine_in_proc_cnt
        j = turbine_in_proc_array(i)
        call MPI_send( disk_force, nloc, MPI_rprec, j, 5, comm, ierr )
    enddo              
elseif (turbine_in_proc) then
    call MPI_recv( disk_force, nloc, MPI_rprec, 0, 5, comm, status, ierr )
endif          
#endif 

!apply forcing to each node
if (turbine_in_proc) then
    do s=1,nloc        
        do l=1,wind_farm%turbine(s)%num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            ind2 = wind_farm%turbine(s)%ind(l)
            fxa(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(1)*ind2
            fya(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(2)*ind2
            fza(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(3)*ind2
        enddo
    enddo
endif    

deallocate(w_uv)

!spatially average velocity at the top of the domain and write to file
if (coord .eq. nproc-1) then
    string1=path // 'output/vel_top_of_domain.dat'
    open(unit=1,file=string1,status='unknown',form='formatted',action='write',position='append')
    write(1,*) total_time, sum(u(:,:,nz-1))/(nx*ny)
    close(1)
endif

nullify(y,z)

end subroutine turbines_forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_finalize ()
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_finalize'

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
call turbines_checkpoint
    
!deallocate
deallocate(wind_farm%turbine) 
deallocate(buffer_array)
    
end subroutine turbines_finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_checkpoint ()
use open_file_fid_mod
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_checkpoint'
integer :: fid

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
if (coord == 0) then  
    string1 = path // 'turbine/turbine_u_d_T.dat'
    fid = open_file_fid( string1, 'rewind', 'formatted' )
    do i=1,nloc
        write(fid,*) wind_farm%turbine(i)%u_d_T
    enddo           
    write(fid,*) T_avg_dim
    close (fid)
endif
    
end subroutine turbines_checkpoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbine_vel_init(zo_high)
!  called from ic.f90 if initu, lbc_mom==1, S_FLAG are all false.
!  this accounts for the turbines when creating the initial velocity profile.

use param, only: zo
implicit none
character (*), parameter :: sub_name = mod_name // '.turbine_vel_init'

real(rprec), intent(inout) :: zo_high
real(rprec) :: cft, nu_w, exp_KE, induction_factor, Ct_noprime

! Convert Ct' to Ct
! a = Ct'/(4+Ct'), Ct = 4a(1-a)
induction_factor = Ct_prime / (4._rprec + Ct_prime)
Ct_noprime = 4*(induction_factor) * (1 - induction_factor)

! friction coefficient, cft
cft = pi*Ct_noprime/(4.*sx*sy)
       
!wake viscosity
nu_w = 28.*sqrt(0.5*cft)

!turbine friction height, Calaf, Phys. Fluids 22, 2010
zo_high = height_all*(1.+0.5*dia_all/height_all)**(nu_w/(1.+nu_w))* &
  exp(-1.*(0.5*cft/(vonk**2) + (log(height_all/zo* &
  (1.-0.5*dia_all/height_all)**(nu_w/(1.+nu_w))) )**(-2) )**(-0.5) )

exp_KE =  0.5*(log(0.45/zo_high)/0.4)**2
  
if(.false.) then
    write(*,*) 'sx,sy,cft: ',sx,sy,cft
    write(*,*) 'nu_w: ',nu_w    
    write(*,*) 'zo_high: ',zo_high
    write(*,*) 'approx expected KE: ', exp_KE
endif
end subroutine turbine_vel_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module turbines
