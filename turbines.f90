!!
!!  Copyright (C) 2010-2013  Johns Hopkins University
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
! This module currently contains all of the subroutines associated with drag-disk turbines:
!  (1) turbines_init
!      This routine creates the 'turbine' folder and starts the turbine forcing output files.
!      It also creates the indicator function (Gaussian-filtered from binary locations - in or out)
!      and sets values for turbine type (node locations, etc)
!  (2) turbines_forcing
!      Applies the drag-disk forcing
use param
use turbines_base
use grid_defs, only: grid 
use messages
use string_util

implicit none
include 'tecryte.h'

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize

real(rprec) :: Ct_prime_05
real(rprec) :: T_avg_dim_file
real(rprec), dimension(:), allocatable :: z_tot

character (100) :: string1
!real(rprec), dimension(:,:,:), allocatable :: large_node_array    !used for visualizing node locations. Richard: ! removed large 3D array to limit memory use
!real(rprec), dimension(:,:,:), allocatable :: large_node_array_filtered ! Richard: ! removed large 3D array to limit memory use
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
character (*), parameter :: mod_name = 'turbines'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_init()

use open_file_fid_mod
implicit none
real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.turbines_init'

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

! Allocate and initialize
!allocate(large_node_array(nx,ny,nz_tot)) ! removed large 3D array to limit memory use
!allocate(large_node_array_filtered(nx,ny,nz_tot)) ! removed large 3D array to limit memory use
allocate(turbine_in_proc_array(nproc-1))
allocate(z_tot(nz_tot))
allocate(file_id(nloc))
allocate(file_id2(nloc))
turbine_in_proc_array = 0

nullify(buffer_array)
allocate(buffer_array(nloc))

!new variables for optimization:
    Ct_prime_05 = -0.5*Ct_prime

!z_tot for total domain (since z is local to the processor)
    do k=1,nz_tot
        z_tot(k) = (k - 0.5_rprec) * dz
    enddo

!find turbine nodes - including unfiltered ind, n_hat, num_nodes, and nodes for each turbine
!each processor finds turbines in the entire domain
!    large_node_array = 0.
!    call turbines_nodes(large_node_array)
    call turbines_nodes                  ! removed large 3D array to limit memory use
    
!    if (coord == 0) then
!      !to write the node locations to file
!      string1 = path // 'turbine/nodes_unfiltered.dat'
!      call write_tecplot_header_ND(string1,'rewind', 4, (/nx+1, ny+1, nz_tot/), '"x", "y", "z", "nodes_unfiltered"', numtostr(0,1), 1)
!      call write_real_data_3D(string1, 'append','formatted', 1, nx, ny, nz_tot, (/large_node_array/), 4, x,y,z_tot)
!    endif

!1.smooth/filter indicator function                     
!2.associate new nodes with turbines                               
!3.normalize such that each turbine's ind integrates to turbine volume
!4.split domain between processors 
    call turbines_filter_ind()
    
    if (turbine_cumulative_time) then
        if (coord == 0) then
            string1 = path // 'turbine/turbine_u_d_T.dat'
            inquire (file=string1, exist=exst)
            if (exst) then
                write(*,*) 'Reading from file turbine_u_d_T.dat'
                inquire (unit=1, opened=opn)
                if (opn) call error (sub_name, 'unit 1 already open, mark1')                
                open (1, file=string1)
                do i=1,nloc
                    read(1,*) wind_farm%turbine(i)%u_d_T    
                enddo    
                read(1,*) T_avg_dim_file
                if (T_avg_dim_file /= T_avg_dim) then
                    write(*,*) 'Time-averaging window does not match value in turbine_u_d_T.dat'
                endif
                close (1)
            else  
                write (*, *) 'File ', trim(string1), ' not found'
                write (*, *) 'Assuming u_d_T = -7. for all turbines'
                do k=1,nloc
                    wind_farm%turbine(k)%u_d_T = -7.
                enddo
            endif                                         
        endif
    else
        write (*, *) 'Assuming u_d_T = -7 for all turbines'
        do k=1,nloc
            wind_farm%turbine(k)%u_d_T = -7.
        enddo    
    endif
   
if (coord .eq. nproc-1) then
    string1=path // 'output/vel_top_of_domain.dat'
    open(unit=1,file=string1,status='unknown',form='formatted',action='write',position='rewind')
    write(1,*) 'total_time','u_HI'
    close(1)
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

$if ($MPI)
   kcp = nint(wind_farm%turbine(s)%height/dz + 0.5)
   k_start =  1+coord*(nz-1)
   k_end = nz-1+coord*(nz-1)

   if (kcp>=k_start .and. kcp<=k_end) then
       call string_splice( string1, path // 'turbine/turbine_', s, '_velcenter.dat' )
       call string_concat (string1, '.c', coord)
       file_id2(s) = open_file_fid( string1, 'append', 'formatted' )
   endif

$else
   call string_splice( string1, path // 'turbine/turbine_', s, '_velcenter.dat' )
   file_id2(s) = open_file_fid( string1, 'append', 'formatted' )
$endif

enddo

nullify(x,y,z)

end subroutine turbines_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine turbines_nodes(array) ! removed large 3D array to limit memory use
subroutine turbines_nodes
!This subroutine locates nodes for each turbine and builds the arrays: ind, n_hat, num_nodes, and nodes
implicit none
real(rprec) :: R_t,rx,ry,rz,r,r_norm,r_disk
!real(rprec), dimension(nx,ny,nz_tot) :: array ! removed large 3D array to limit memory use
real(rprec), pointer :: p_xloc => null(), p_yloc=> null(), p_height=> null()
real(rprec), pointer :: p_dia => null(), p_thk=> null(), p_theta1=> null(), p_theta2=> null()
real(rprec), pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3=> null() 
real(rprec), pointer, dimension(:) :: x, y, z

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

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
    R_t = p_dia/2.
            !if (coord == 0) then
            !    if (verbose) then
            !      write(*,*) '     rad:',R_t
            !      write(*,*) '     dx:',dx
            !      write(*,*) '     dy:',dy
            !      write(*,*) '     dz:',dz
            !      write(*,*) '     thk:',p_thk
            !    endif
            !endif
    imax = R_t/dx + 2
    jmax = R_t/dy + 2
    kmax = R_t/dz + 2

    !determine unit normal vector for each turbine	
    p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat3 = sin(pi*p_theta2/180.)
            !if (coord == 0) then
            !    if (verbose) then
            !      write(*,*) '     n_hat', p_nhat1, p_nhat2, p_nhat3
            !    endif
            !endif

    !determine nearest (i,j,k) to turbine center
    icp = nint(p_xloc/dx)
    jcp = nint(p_yloc/dy)
    kcp = nint(p_height/dz + 0.5)
            !if (coord == 0) then
            !    if (verbose) then
            !      write(*,*) '     turbine center (i,j,k)',icp,jcp,kcp
            !      write(*,*) '     turbine center (x,y,z)',x(icp),y(jcp),z_tot(kcp)
            !    endif
            !endif

    !determine limits for checking i,j,k
        !due to spectral BCs, i and j may be < 1 or > nx,ny
        !the mod function accounts for this when these values are used
    min_i = icp-imax
    max_i = icp+imax
    min_j = jcp-jmax
    max_j = jcp+jmax
    min_k = max((kcp-kmax),1)
    max_k = min((kcp+kmax),nz_tot)
            !if (coord == 0) then
            !    if (verbose) then
            !      write(*,*) '     i limits, range', min_i, max_i, dx*(max_i-min_i)
            !      write(*,*) '     j limits, range', min_j, max_j, dy*(max_j-min_j)
            !      write(*,*) '     k limits, range', min_k, max_k, dz*(max_k-min_k)
            !    endif
            !endif
            wind_farm%turbine(s)%nodes_max(1) = min_i
            wind_farm%turbine(s)%nodes_max(2) = max_i
            wind_farm%turbine(s)%nodes_max(3) = min_j
            wind_farm%turbine(s)%nodes_max(4) = max_j
            wind_farm%turbine(s)%nodes_max(5) = min_k
            wind_farm%turbine(s)%nodes_max(6) = max_k            

    !check neighboring grid points	
    do k=min_k,max_k
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
                !if r_disk<R_t and r_norm within thk/2 from turbine -- this node is part of the turbine
                if ( (r_disk .LE. R_t) .AND. (r_norm .LE. p_thk/2.) ) then
                    !if (coord == 0) then
                    !    if (verbose) then
                    !      write(*,*) '     FOUND NODE', i,j,k
                    !    endif
                    !endif
!                    array(i2,j2,k) = 1.  ! removed large 3D array to limit memory use
                    wind_farm%turbine(s)%ind(count_i) = 1. 
                    wind_farm%turbine(s)%nodes(count_i,1) = i2
                    wind_farm%turbine(s)%nodes(count_i,2) = j2
                    wind_farm%turbine(s)%nodes(count_i,3) = k   !global k (might be out of this proc's range)
                    count_n = count_n + 1
                    count_i = count_i + 1
                endif
           enddo
       enddo
    enddo
    wind_farm%turbine(s)%num_nodes = count_n

    if (coord == 0) then
        write(*,*) '     Turbine #',s,'has',count_n,'unfiltered nodes in entire domain'
    endif
    
enddo

nullify(x,y,z)
end subroutine turbines_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_filter_ind()
! This subroutine takes ind and nodes for each turbine and filters according to
! alpha and ifilter from wind_farm
!       1.smooth/filter indicator function                                  CHANGE IND
!       2.normalize such that each turbine's ind integrates to 1.           CHANGE IND
!       3.associate new nodes with turbines                                 CHANGE NODES, NUM_NODES       

implicit none
!real(rprec), dimension(nx,ny,nz_tot) :: out_a, g, g_shift, fg
! Richard: Commented out g_shift. This large array is not used  fg can easily be replaced by single double
real(rprec), dimension(nx,ny,nz_tot) :: out_a, g
real(rprec), dimension(nx,ny,nz_tot) :: temp_array
!real(rprec), dimension(nx,ny,nz) :: temp_array_2 ! Richard: Commented out to save i/o operations and memory
real(rprec) :: sumG,delta2,r2,sumA
real(rprec) :: turbine_vol
real(rprec) :: fg ! Removed the 3D matrix as it is only used as a temporary dummy variable

real(rprec), pointer, dimension(:) :: x,y,z

!logical :: verbose = .false.
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

!create convolution function, centered at (nx/2,ny/2,(nz_tot-1)/2) and normalized
!if(wind_farm%ifilter==2) then		!2-> Gaussian
delta2 = alpha**2 * (dx**2 + dy**2 + dz**2)
      do k=1,nz_tot
        do j=1,ny
          do i=1,nx
            r2 = ((real(i)-nx/2.)*dx)**2 + ((real(j)-ny/2.)*dy)**2 + ((real(k)-(nz_tot-1)/2.)*dz)**2
            g(i,j,k) = sqrt(6./(pi*delta2))*6./(pi*delta2)*exp(-6.*r2/delta2)
          enddo
        enddo
      enddo
!endif

!normalize the convolution function
sumG = sum(g(:,:,:))*dx*dy*dz
g = g/sumG

!display the convolution function
    !if (coord == 0) then
    !    if(.false.) then
    !      write(*,*) 'Convolution function'
    !      write(*,*) g
    !      write(*,*) 'integral of g(i,j,k): ',sumG
    !    endif       
    !endif

    !to write the data to file, centered at (i,j,k=(nz_tot-1)/2)
!    if (coord == 0) then    
!        i=nx/2
!        j=ny/2
!        do k2=1,nz_tot
!          do j2=1,ny
!            do i2=1,nx
!            g_shift(i2,j2,k2) = g( mod(i2-i+nx/2+nx-1,nx)+1 , mod(j2-j+ny/2+ny-1,ny)+1, k2)
!            enddo
!          enddo
!        enddo

        !if (.false.) then
        !    fname0 = path // 'turbine/convolution_function.dat'
        !    call write_tecplot_header_ND(fname0,'rewind', 4, (/nx,ny,nz_tot/), '"x","y","z","g"', convtostr(1,1), 1)
        !    call write_real_data_3D(fname0, 'append', 'formatted', 1, nx, ny, nz_tot, (/g_shift/), 0, x, y, z_tot)

        !    if (verbose) then
        !        write(*,*) 'Convolution function written to Tecplot file.'
        !    endif
        !endif
!    endif

!filter indicator function for each turbine
do b=1,nloc
    
    !if (coord == 0) then
    !    if (verbose) then
    !        write(*,*) 'Filtering turbine Number ',b
    !    endif
    !endif

    !create the input array (nx,ny,nz_tot) from a list of included nodes
        temp_array = 0.
        do l=1,wind_farm%turbine(b)%num_nodes
            i2 = wind_farm%turbine(b)%nodes(l,1)
            j2 = wind_farm%turbine(b)%nodes(l,2)
            k2 = wind_farm%turbine(b)%nodes(l,3)
            temp_array(i2,j2,k2) = wind_farm%turbine(b)%ind(l)
        enddo

    !perform convolution on temp_array --> out_a    
        out_a=0.

        min_i = wind_farm%turbine(b)%nodes_max(1) 
        max_i = wind_farm%turbine(b)%nodes_max(2) 
        min_j = wind_farm%turbine(b)%nodes_max(3) 
        max_j = wind_farm%turbine(b)%nodes_max(4) 
        min_k = wind_farm%turbine(b)%nodes_max(5)
        max_k = wind_farm%turbine(b)%nodes_max(6) 
        cut = trunc   

        !if (coord == 0) then
        !    if (verbose) then
        !        write(*,*) 'search over: ',min_i-cut,max_i+cut,min_j-cut,max_j+cut,min_k-cut,max_k+cut
        !    endif
        !endif

        !convolution computed for points (i4,j4,k)
        !only compute for nodes near the turbine (defined by cut aka trunc)
        do k=max(min_k-cut,1),min(max_k+cut,nz_tot)    
        do j=(min_j-cut),(max_j+cut)
        do i=(min_i-cut),(max_i+cut)
        
            i4 = mod(i+nx-1,nx)+1       !since values may be out 1-nx,1-ny domain (spectral BCs)
            j4 = mod(j+ny-1,ny)+1              
        
          !for each (i4,j4,k), center convolution function on that point and 'integrate' 
          !relative coords are (ssx,ssy,ssz). absolute coords of other/surrounding points are (i2,j2,k2)
          !only need to consider other/surrounding points near (i4,j4,k) since conv. function is compact
          do k2=max(k-trunc,1),min(k+trunc,nz_tot)     !currently using truncated Gaussian
          do j2=j-trunc,j+trunc
          do i2=i-trunc,i+trunc

            i3 = mod(i2+nx-1,nx)+1      !since values may be out 1-nx,1-ny domain (spectral BCs)
            j3 = mod(j2+ny-1,ny)+1             
          
            ssx = mod(i2-i+nx/2+nx-1,nx)+1
            ssy = mod(j2-j+ny/2+ny-1,ny)+1       
            ssz = k2-k+(nz_tot-1)/2       !since no spectral BCs in z-direction
                             
            if( ssz < 1) then
                !fg(i2,j2,k2) = 0.
                fg=0.
                write(*,*) 'See turbines.f90, ssz < 1'                    
            elseif( ssz > nz_tot ) then
                !fg(i2,j2,k2) = 0.
                fg=0.
                write(*,*) 'See turbines.f90, ssz > nz_tot'                    
            else
                !fg(i3,j3,k2) = temp_array(i3,j3,k2)*g(ssx,ssy,ssz)
                !out_a(i4,j4,k) = out_a(i4,j4,k) + fg(i3,j3,k2)*dx*dy*dz
                fg = temp_array(i3,j3,k2)*g(ssx,ssy,ssz)
                out_a(i4,j4,k) = out_a(i4,j4,k) + fg*dx*dy*dz
            endif    
        
        enddo
       enddo
      enddo
     enddo
    enddo
        enddo

        !if (coord == 0) then
        !    if (verbose) then
        !        write(*,*) 'Convolution complete for turbine ',b
        !    endif
        !endif

    !normalize this "indicator function" such that it integrates to turbine volume
    sumA = 0.
    do k=1,nz_tot
     do j=1,ny
      do i=1,nx
            if (out_a(i,j,k) < filter_cutoff) then
            out_a(i,j,k) = 0.     !don't want to include too many nodes (truncated Gaussian?)
            else
                sumA = sumA + out_a(i,j,k)*dx*dy*dz
            endif            
      enddo
     enddo
    enddo
    turbine_vol = pi/4. * (wind_farm%turbine(b)%dia)**2 * wind_farm%turbine(b)%thk
    out_a = turbine_vol/sumA*out_a

!    if (coord == 0) then
!        large_node_array_filtered = large_node_array_filtered + out_a ! removed large 3D array to limit memory use
!    endif

    !update num_nodes, nodes, and ind for this turbine
    !and split domain between processors
    !z(nz) and z(1) of neighboring coords match so each coord gets (local) 1 to nz-1
    wind_farm%turbine(b)%ind = 0.
    wind_farm%turbine(b)%nodes = 0.
    wind_farm%turbine(b)%num_nodes = 0.
    count_n = 0
    count_i = 1
    $if ($MPI)
        k_start = 1+coord*(nz-1)
        k_end = nz-1+coord*(nz-1)
    $else
        k_start = 1
        k_end = nz
    $endif

    do k=k_start,k_end  !global k     
     do j=1,ny
      do i=1,nx
       if (out_a(i,j,k) > filter_cutoff) then
        wind_farm%turbine(b)%ind(count_i) = out_a(i,j,k)
        wind_farm%turbine(b)%nodes(count_i,1) = i
        wind_farm%turbine(b)%nodes(count_i,2) = j
        wind_farm%turbine(b)%nodes(count_i,3) = k - coord*(nz-1)   !local k
        count_n = count_n + 1
        count_i = count_i + 1
        turbine_in_proc = .true.                    
       endif
      enddo
     enddo
    enddo
    wind_farm%turbine(b)%num_nodes = count_n
    
    if (count_n > 0) then
        $if ($MPI)

            call string_splice( string1, 'Turbine number ', b,' has ', count_n,' filtered nodes in coord ', coord )
            write(*,*) trim(string1)
            
        $else

            call string_splice( string1, 'Turbine number ',b,' has ',count_n,' filtered nodes' )
            write(*,*) trim(string1)

        $endif
    endif

enddo

!    !test to make sure domain is divided correctly:
!        temp_array_2 = 0.
!        do b=1,nloc
!        do l=1,wind_farm%turbine(b)%num_nodes
!            i2 = wind_farm%turbine(b)%nodes(l,1)
!            j2 = wind_farm%turbine(b)%nodes(l,2)
!            k2 = wind_farm%turbine(b)%nodes(l,3)
!            temp_array_2(i2,j2,k2) = wind_farm%turbine(b)%ind(l)
!        enddo   
!        enddo
!        !write to file with .dat.c* extension
!        call string_splice( string1, path // 'turbine/nodes_filtered_c.dat' // '.c', coord )
!        
!        call write_tecplot_header_ND(string1,'rewind', 4, (/nx,ny,nz/), '"x","y","z","nodes_filtered_c"', numtostr(1,1), 1)
!        call write_real_data_3D(string1, 'append', 'formatted', 1, nx, ny, nz, (/temp_array_2/), 0, x, y, z(1:nz))      

!    if (coord == 0) then
!        string1 = path // 'turbine/nodes_filtered.dat'
!        call write_tecplot_header_ND(string1,'rewind', 4, (/nx,ny,nz_tot/), '"x","y","z","nodes_filtered"', numtostr(1,1), 1)
!        call write_real_data_3D(string1, 'append', 'formatted', 1, nx, ny, nz_tot, (/large_node_array_filtered/), 0, x, y, z_tot)                       
!    endif

!each processor sends its value of turbine_in_proc
!if false, disk-avg velocity will not be sent (since it will always be 0.)
!############################################## 2
$if ($MPI)
    if (coord == 0) then
        if (turbine_in_proc) then
            write(*,*),'Coord 0 has turbine nodes' 
        endif
        do i=1,nproc-1
            call MPI_recv( buffer_logical, 1, MPI_logical, i, 2, comm, status, ierr )

            if (buffer_logical) then
                call string_splice( string1, 'Coord', i,' has turbine nodes')
                write(*,*) trim(string1)
                turbine_in_proc_cnt = turbine_in_proc_cnt + 1
                turbine_in_proc_array(turbine_in_proc_cnt) = i
            endif
        enddo
    else
        call MPI_send( turbine_in_proc, 1, MPI_logical, 0, 2, comm, ierr )
    endif
$endif
nullify(x,y,z)
!##############################################
end subroutine turbines_filter_ind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_forcing()
use sim_param, only: u,v,w, fxa
use functions, only: interp_to_uv_grid

implicit none
real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vels, disk_force
real(rprec), allocatable, dimension(:,:,:) :: w_uv ! Richard: This 3D matrix can relatively easy be prevented
real(rprec), pointer, dimension(:) :: y,z

nullify(y,z)
y => grid % y
z => grid % z

allocate(w_uv(ld,ny,lbz:nz))

$if ($MPI)
    call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)     !syncing intermediate w-velocities!
$endif

w_uv = interp_to_uv_grid(w, lbz)

disk_avg_vels = 0.

!Each processor calculates the weighted disk-averaged velocity
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
                !write center disk velocity to file                   
                if (modulo (jt_total, tbase) == 0) then
                icp = nint(wind_farm%turbine(s)%xloc/dx)
                jcp = nint(wind_farm%turbine(s)%yloc/dy)
                kcp = nint(wind_farm%turbine(s)%height/dz + 0.5)
                $if ($MPI)
                k_start =  1+coord*(nz-1)
                k_end = nz-1+coord*(nz-1)
                $else
                k_start = 1
                k_end = nz
                $endif
                if (kcp>=k_start .and. kcp<=k_end) then
                kcp=kcp-k_start+1
                call write_real_data( file_id2(s), 'formatted', 3, (/ total_time_dim, u(icp,jcp,kcp), w_uv(icp,jcp,kcp) /))
                endif
                endif
        !write this value to the array (which will be sent to coord 0)
            disk_avg_vels(s) = p_u_d
            
        enddo
        
endif        

!send the disk-avg values to coord==0
$if ($MPI) 

    call mpi_barrier (comm,ierr)

    !############################################## 3
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
    !##############################################              
$endif

!Coord==0 takes that info and calculates total disk force, then sends it back
if (coord == 0) then           
    !update epsilon for the new timestep (for cfl_dt) 
        eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)

    !for each turbine:        
        do s=1,nloc            
             
        !set pointers
            p_u_d => wind_farm%turbine(s)%u_d   
            p_u_d_T => wind_farm%turbine(s)%u_d_T   
            p_f_n => wind_farm%turbine(s)%f_n                  
            
        !volume correction:
        !since sum of ind is turbine volume/(dx*dy*dz) (not exactly 1.)
            p_u_d = disk_avg_vels(s) * wind_farm%turbine(s)%vol_c
    
        !add this current value to the "running average" (first order relaxation)
            p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

        !calculate total thrust force for each turbine  (per unit mass)
        !force is normal to the surface (calc from u_d_T, normal to surface)
            p_f_n = Ct_prime_05*abs(p_u_d_T)*p_u_d_T/wind_farm%turbine(s)%thk       
            !write values to file                   
                if (modulo (jt_total, tbase) == 0) then
                call write_real_data( file_id(s), 'formatted', 5, (/total_time_dim, p_u_d, p_u_d_T, p_f_n, Ct_prime_05*(p_u_d_T*p_u_d_T*p_u_d_T)*pi/(4.*sx*sy) /))
                endif 
            !write force to array that will be transferred via MPI    
            disk_force(s) = p_f_n
            
        enddo           
        
endif

!send total disk force to the necessary procs (with turbine_in_proc==.true.)
$if ($MPI)
    !############################################## 5  
        if (coord == 0) then
          
            do i=1,turbine_in_proc_cnt
                j = turbine_in_proc_array(i)
                call MPI_send( disk_force, nloc, MPI_rprec, j, 5, comm, ierr )
            enddo              
                        
        elseif (turbine_in_proc) then
            call MPI_recv( disk_force, nloc, MPI_rprec, 0, 5, comm, status, ierr )
        endif     
    !##############################################     
$endif 
    
!apply forcing to each node
if (turbine_in_proc) then

    do s=1,nloc     
            
        do l=1,wind_farm%turbine(s)%num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            ind2 = wind_farm%turbine(s)%ind(l)
            fxa(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(1)*ind2 
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
    if (coord == 0) then  
        string1 = path // 'turbine/turbine_u_d_T.dat'    
        inquire (unit=1, opened=opn)
        if (opn) call error (sub_name, 'unit 1 already open, mark6')        
        open (unit=1,file = string1, status='unknown',form='formatted', action='write',position='rewind')
        do i=1,nloc
            write(1,*) wind_farm%turbine(i)%u_d_T
        enddo           
        write(1,*) T_avg_dim
        close (1)
    endif
    
!deallocate
deallocate(wind_farm%turbine) 
deallocate(buffer_array)
    
end subroutine turbines_finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbine_vel_init(zo_high)
!  called from ic.f90 if initu, dns_bc, S_FLAG are all false.
!  this accounts for the turbines when creating the initial velocity profile.

use param, only: zo
implicit none
real(rprec), intent(inout) :: zo_high
real(rprec) :: cft,nu_w,exp_KE

!friction coefficient, cft
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
end module turbines
