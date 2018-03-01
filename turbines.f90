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
!!public :: turbines_forcing_fourier
public :: turbines_RNL2
public :: turbines_RNL3  ! seems to work
public :: turbines_RNL4  ! experimental
public :: turbines_RNL5  ! trying w/ disk averaging
public :: disk_avg_kx

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

character(*), parameter :: param_dat = path // 'turb_param.dat'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine disk_avg_kx()
use sim_param, only: u
use derivatives, only: dft_direct_back_2d_n_yonlyC
implicit none

real(rprec), dimension(ld) :: buff
real(rprec), dimension(nloc, ld) :: disk_avgs, disk_buff
real(rprec), dimension(ld, ny, lbz:nz) :: vel_temp
real(rprec), dimension(ld, ny) :: uhat
integer :: jz

disk_avgs(:,:) = 0._rprec
disk_buff(:,:) = 0._rprec
buff(:) = 0._rprec
vel_temp(:,:,:) = 0._rprec
uhat(:,:) = 0._rprec

if (turbine_in_proc) then
   vel_temp(:,:,:) = 0._rprec
   uhat(:,:) = 0._rprec
   do jz=1,nz
      uhat(:,:) = u(:,:,jz)
      call dft_direct_back_2d_n_yonlyC( uhat(:,:) )
      vel_temp(:,:,jz) = uhat(:,:)
   enddo
endif

if (turbine_in_proc) then
   ! for each turbine
   do s=1,nloc
      buff(:) = 0._rprec
      do l=1,wind_farm%turbine(s)%num_nodes
         j2 = wind_farm%turbine(s)%nodes(l,2)
         k2 = wind_farm%turbine(s)%nodes(l,3)
         buff(:) = buff(:) + vel_temp(:, j2, k2)
      enddo
      disk_avgs(s,:) = buff(:) / real(wind_farm%turbine(s)%num_nodes,rprec)
   enddo
endif

call mpi_barrier(comm, ierr)
if (coord == 0) then
   do i=1,turbine_in_proc_cnt
      j = turbine_in_proc_array(i)
      disk_buff = 0.
      call MPI_recv( disk_buff, nloc*ld, MPI_rprec, j, 3, comm, status, ierr )
      disk_avgs(:,:) = disk_avgs(:,:) + disk_buff(:,:)
   enddo
elseif (turbine_in_proc) then
   call MPI_send( disk_avgs(:,:), nloc*ld, MPI_rprec, 0, 3, comm, ierr )
endif
call mpi_barrier(comm, ierr)

!!$if (coord == 0) then
!!$   do s=1,nloc
!!$      print*, 'sss', s, disk_avgs(s,:)
!!$   enddo
!!$endif

call mpi_barrier(comm, ierr)
if (coord == 0) then

   do i=1,turbine_in_proc_cnt
      j = turbine_in_proc_array(i)
      call MPI_send( disk_avgs(:,:), nloc*ld, MPI_rprec, j, 5, comm, ierr )
   enddo
elseif (turbine_in_proc) then

   call MPI_recv( disk_buff(:,:), nloc*ld, MPI_rprec, 0, 5, comm, status, ierr )
   disk_avgs(:,:) = disk_buff(:,:)

endif

do s=1,nloc
   wind_farm%turbine(s)%u_d_kx(:) = disk_avgs(s,:)
enddo

end subroutine disk_avg_kx

subroutine turbines_init()

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

if (fourier) then   !!jb
   do s=1,nloc
      allocate( wind_farm % turbine(s) % u_d_kx(ld) )
   enddo
endif
!call disk_avg_kx()   !!jb

!new variables for optimization:
    Ct_prime_05 = -0.5*Ct_prime

!z_tot for total domain (since z is local to the processor)
    do k=1,nz_tot
        z_tot(k) = (k - 0.5_rprec) * dz
    enddo

    call place_turbines    !!jb

!find turbine nodes - including unfiltered ind, n_hat, num_nodes, and nodes for each turbine
!each processor finds turbines in the entire domain
!    large_node_array = 0.
!    call turbines_nodes(large_node_array)
    if (fourier) then
       call turbines_nodes_fourier
    else
       call turbines_nodes            ! removed large 3D array to limit memory use
    endif
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

    print*, 'coord, turbine_in_proc: ', coord, turbine_in_proc
    if (coord == 1) then
    print*, '>>before >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
    do s=1,nloc
       do l=1,wind_farm % turbine(s) % num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            k = wind_farm%turbine(s)%ind(l)
            print*,'s,l,i,j,k,ind',s, l, i2, j2, k2, k
       enddo
    enddo
    endif

    if (fourier) then
       !call turbines_filter_ind_fourier()
       call turbines_filter_ind2()
    else
       call turbines_filter_ind()
    endif
    
    print*, 'coord, turbine_in_proc: ', coord, turbine_in_proc
    if (coord == 1) then
    print*, '>>after >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
    do s=1,nloc
       do l=1,wind_farm % turbine(s) % num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            k = wind_farm%turbine(s)%ind(l)
            print*,'s,l,i,j,k,ind',s, l, i2, j2, k2, k
       enddo
    enddo
    endif
  


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
   file_id(s) = open_file( string1, 'append', 'formatted' )
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
       file_id2(s) = open_file( string1, 'append', 'formatted' )
   endif

$else
   call string_splice( string1, path // 'turbine/turbine_', s, '_velcenter.dat' )
   file_id2(s) = open_file( string1, 'append', 'formatted' )
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

subroutine turbines_nodes_fourier
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
                    i2 = mod(i+nxp-1,nxp)+1
                    rx = (x(i2)-L_x) - p_xloc
                elseif (i>nxp) then
                    i2 = mod(i+nxp-1,nxp)+1
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
end subroutine turbines_nodes_fourier

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
    if (coord == 0) then
       write(*,*) 'Convolution function'
       !write(*,*) g
       write(*,*) 'integral of g(i,j,k): ',sumG
    endif

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
subroutine turbines_filter_ind2()
! This subroutine takes ind and nodes for each turbine and filters according to
! alpha and ifilter from wind_farm
!       1.smooth/filter indicator function                                  CHANGE IND
!       2.normalize such that each turbine's ind integrates to 1.           CHANGE IND
!       3.associate new nodes with turbines                                 CHANGE NODES, NUM_NODES       

implicit none
!real(rprec), dimension(nx,ny,nz_tot) :: out_a, g, g_shift, fg
! Richard: Commented out g_shift. This large array is not used  fg can easily be replaced by single double
real(rprec), dimension(ny,nz_tot) :: out_a, g
real(rprec), dimension(ny,nz_tot) :: temp_array
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
delta2 = alpha**2 * (dy**2 + dz**2)
      do k=1,nz_tot
        do j=1,ny
            r2 = ((real(j)-ny/2.)*dy)**2 + ((real(k)-(nz_tot-1)/2.)*dz)**2
            g(j,k) = 6./(pi*delta2)*exp(-6.*r2/delta2)
        enddo
      enddo
!endif

!normalize the convolution function
sumG = sum(g(:,:))*dy*dz
g = g/sumG

    !display the convolution function
    if (coord == 0) then
       write(*,*) '2D Convolution function'
       !write(*,*) g
       write(*,*) 'integral of g(i,j,k): ',sumG
    endif

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
            j2 = wind_farm%turbine(b)%nodes(l,2)
            k2 = wind_farm%turbine(b)%nodes(l,3)
            temp_array(j2,k2) = wind_farm%turbine(b)%ind(l)
        enddo

    !perform convolution on temp_array --> out_a    
        out_a=0.

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
        
            j4 = mod(j+ny-1,ny)+1              
        
          !for each (i4,j4,k), center convolution function on that point and 'integrate' 
          !relative coords are (ssx,ssy,ssz). absolute coords of other/surrounding points are (i2,j2,k2)
          !only need to consider other/surrounding points near (i4,j4,k) since conv. function is compact
          do k2=max(k-trunc,1),min(k+trunc,nz_tot)     !currently using truncated Gaussian
          do j2=j-trunc,j+trunc
            j3 = mod(j2+ny-1,ny)+1             
          
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
                fg = temp_array(j3,k2)*g(ssy,ssz)
                out_a(j4,k) = out_a(j4,k) + fg*dy*dz
            endif    
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
            if (out_a(j,k) < filter_cutoff) then
            out_a(j,k) = 0.     !don't want to include too many nodes (truncated Gaussian?)
            else
                sumA = sumA + out_a(j,k)*dy*dz
            endif
     enddo
    enddo
    turbine_vol = pi/4. * (wind_farm%turbine(b)%dia)**2 !* wind_farm%turbine(b)%thk
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
       if (out_a(j,k) > filter_cutoff) then
        wind_farm%turbine(b)%ind(count_i) = out_a(j,k)
        !wind_farm%turbine(b)%nodes(count_i,1) = i
        wind_farm%turbine(b)%nodes(count_i,2) = j
        wind_farm%turbine(b)%nodes(count_i,3) = k - coord*(nz-1)   !local k
        count_n = count_n + 1
        count_i = count_i + 1
        turbine_in_proc = .true.                    
       endif
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
end subroutine turbines_filter_ind2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_filter_ind_fourier()
! This subroutine takes ind and nodes for each turbine and filters according to
! alpha and ifilter from wind_farm
!       1.smooth/filter indicator function                                  CHANGE IND
!       2.normalize such that each turbine's ind integrates to 1.           CHANGE IND
!       3.associate new nodes with turbines                                 CHANGE NODES, NUM_NODES       

implicit none
!real(rprec), dimension(nx,ny,nz_tot) :: out_a, g, g_shift, fg
! Richard: Commented out g_shift. This large array is not used  fg can easily be replaced by single double
real(rprec), dimension(nxp,ny,nz_tot) :: out_a, g
real(rprec), dimension(nxp,ny,nz_tot) :: temp_array
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
          do i=1,nxp
            r2 = ((real(i)-nxp/2.)*dx)**2 + ((real(j)-ny/2.)*dy)**2 + ((real(k)-(nz_tot-1)/2.)*dz)**2
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
        
            i4 = mod(i+nxp-1,nxp)+1       !since values may be out 1-nx,1-ny domain (spectral BCs)
            j4 = mod(j+ny-1,ny)+1              
        
          !for each (i4,j4,k), center convolution function on that point and 'integrate' 
          !relative coords are (ssx,ssy,ssz). absolute coords of other/surrounding points are (i2,j2,k2)
          !only need to consider other/surrounding points near (i4,j4,k) since conv. function is compact
          do k2=max(k-trunc,1),min(k+trunc,nz_tot)     !currently using truncated Gaussian
          do j2=j-trunc,j+trunc
          do i2=i-trunc,i+trunc

            i3 = mod(i2+nxp-1,nxp)+1      !since values may be out 1-nx,1-ny domain (spectral BCs)
            j3 = mod(j2+ny-1,ny)+1             
          
            ssx = mod(i2-i+nxp/2+nxp-1,nxp)+1
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
      do i=1,nxp
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
      do i=1,nxp
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
end subroutine turbines_filter_ind_fourier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_forcing()
use sim_param, only: u,v,w, fxa
use functions, only: interp_to_uv_grid
use derivatives, only: dft_direct_back_2d_n_yonlyC, dft_direct_forw_2d_n_yonlyC

implicit none
real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vels, disk_force
real(rprec), allocatable, dimension(:,:,:) :: w_uv ! Richard: This 3D matrix can relatively easy be prevented
real(rprec), pointer, dimension(:) :: y,z
integer :: jz


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
  if (fourier) then
  do jz=1,nz
    call dft_direct_back_2d_n_yonlyC( u(:,:,jz) )
  enddo
  endif


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
                if (fourier) then
                p_u_d = p_u_d + (wind_farm%turbine(s)%nhat(1)*u(1,j2,k2))* wind_farm%turbine(s)%ind(l)  !! just kx=0 part! (purely real)
                else
                p_u_d = p_u_d + (wind_farm%turbine(s)%nhat(1)*u(i2,j2,k2) &
                               + wind_farm%turbine(s)%nhat(2)*v(i2,j2,k2) &
                               + wind_farm%turbine(s)%nhat(3)*w_uv(i2,j2,k2)) &
                    * wind_farm%turbine(s)%ind(l)        
                endif
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
            disk_avg_vels(s) = p_u_d !/ real(wind_farm%turbine(s)%num_nodes, rprec)
            
         enddo

  if (fourier) then
  do jz=1,nz
    call dft_direct_forw_2d_n_yonlyC( u(:,:,jz) )
  enddo
  endif
        
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
            if (fourier) then
               p_u_d = disk_avg_vels(s) * wind_farm%turbine(s)%vol_c
            else
               p_u_d = disk_avg_vels(s) * wind_farm%turbine(s)%vol_c
            endif
        !add this current value to the "running average" (first order relaxation)
            p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

        !calculate total thrust force for each turbine  (per unit mass)
        !force is normal to the surface (calc from u_d_T, normal to surface)
            if (fourier) then
               p_f_n = Ct_prime_05*abs(p_u_d_T)*p_u_d_T !/wind_farm%turbine(s)%thk       
            else
               p_f_n = Ct_prime_05*abs(p_u_d_T)*p_u_d_T/wind_farm%turbine(s)%thk       
            endif

            !write values to file                   
                if (modulo (jt_total, tbase) == 0) then
                !call write_real_data( file_id(s), 'formatted', 5, (/total_time_dim, p_u_d, p_u_d_T, p_f_n, Ct_prime_05*(p_u_d_T*p_u_d_T*p_u_d_T)*pi/(4.*sx*sy) /))
                call write_real_data( file_id(s), 'formatted', 5, (/total_time_dim, p_u_d, p_u_d_T, p_f_n, (p_u_d_T*p_u_d_T*p_u_d_T) /))   !!jb
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
       !print*,'s,nhat,thk: ',s,wind_farm%turbine(s)%nhat(1),wind_farm%turbine(s)%thk,wind_farm%turbine(s)%vol_c
       if (fourier) then
         do l=1,wind_farm%turbine(s)%num_nodes
            !i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            ind2 = wind_farm%turbine(s)%ind(l)
            !print*,'s, l, ind2: ', s, l, ind2
            fxa(1,j2,k2) = fxa(1,j2,k2) + disk_force(s)*wind_farm%turbine(s)%nhat(1)*ind2
         enddo
       else
         do l=1,wind_farm%turbine(s)%num_nodes
            i2 = wind_farm%turbine(s)%nodes(l,1)
            j2 = wind_farm%turbine(s)%nodes(l,2)
            k2 = wind_farm%turbine(s)%nodes(l,3)
            ind2 = wind_farm%turbine(s)%ind(l)
            fxa(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(1)*ind2 
         enddo
       endif

    enddo

  if (fourier) then
     fxa(:,:,:) = fxa(:,:,:) / L_x * real(num_x, rprec)
     do jz=1,nz
        call dft_direct_forw_2d_n_yonlyC( fxa(:,:,jz) )
     enddo  
  endif

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

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$subroutine turbines_forcing_fourier()
!!$use sim_param, only: u,  v,  w,  fxa
!!$use sim_param, only: uF, vF, wF, fxaF
!!$use functions, only: interp_to_uv_grid
!!$use derivatives, only: phys2waveF_pr
!!$
!!$implicit none
!!$real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
!!$real(rprec) :: ind2
!!$real(rprec), dimension(nloc) :: disk_avg_vels, disk_force
!!$real(rprec), allocatable, dimension(:,:,:) :: w_uv ! Richard: This 3D matrix can relatively easy be prevented
!!$real(rprec), pointer, dimension(:) :: y,z
!!$
!!$nullify(y,z)
!!$y => grid % y
!!$z => grid % z
!!$
!!$allocate(w_uv(nxp+2,ny,lbz:nz))
!!$
!!$$if ($MPI)
!!$    call mpi_sync_real_array(wF, 0, MPI_SYNC_DOWNUP)     !syncing intermediate w-velocities!
!!$$endif
!!$
!!$w_uv = interp_to_uv_grid(wF, lbz)
!!$
!!$disk_avg_vels = 0.
!!$
!!$!Each processor calculates the weighted disk-averaged velocity
!!$if (turbine_in_proc) then
!!$
!!$    !for each turbine:        
!!$        do s=1,nloc      
!!$             
!!$        !set pointers
!!$            p_u_d => wind_farm%turbine(s)%u_d   
!!$
!!$        !calculate total disk-averaged velocity for each turbine (current,instantaneous)    
!!$        !u_d is the velocity in the normal direction	  
!!$        !weighted average using "ind"            
!!$            p_u_d = 0.
!!$            do l=1,wind_farm%turbine(s)%num_nodes   
!!$                i2 = wind_farm%turbine(s)%nodes(l,1)
!!$                j2 = wind_farm%turbine(s)%nodes(l,2)
!!$                k2 = wind_farm%turbine(s)%nodes(l,3)
!!$                p_u_d = p_u_d + (wind_farm%turbine(s)%nhat(1)*uF(i2,j2,k2) &
!!$                               + wind_farm%turbine(s)%nhat(2)*vF(i2,j2,k2) &
!!$                               + wind_farm%turbine(s)%nhat(3)*w_uv(i2,j2,k2)) &
!!$                    * wind_farm%turbine(s)%ind(l)        
!!$             enddo
!!$                !write center disk velocity to file                   
!!$                if (modulo (jt_total, tbase) == 0) then
!!$                icp = nint(wind_farm%turbine(s)%xloc/dx)
!!$                jcp = nint(wind_farm%turbine(s)%yloc/dy)
!!$                kcp = nint(wind_farm%turbine(s)%height/dz + 0.5)
!!$                $if ($MPI)
!!$                k_start =  1+coord*(nz-1)
!!$                k_end = nz-1+coord*(nz-1)
!!$                $else
!!$                k_start = 1
!!$                k_end = nz
!!$                $endif
!!$                if (kcp>=k_start .and. kcp<=k_end) then
!!$                kcp=kcp-k_start+1
!!$                call write_real_data( file_id2(s), 'formatted', 3, (/ total_time_dim, uF(icp,jcp,kcp), w_uv(icp,jcp,kcp) /))
!!$                endif
!!$                endif
!!$        !write this value to the array (which will be sent to coord 0)
!!$            disk_avg_vels(s) = p_u_d
!!$            
!!$        enddo
!!$        
!!$endif        
!!$
!!$!send the disk-avg values to coord==0
!!$$if ($MPI) 
!!$
!!$    call mpi_barrier (comm,ierr)
!!$
!!$    !############################################## 3
!!$    if (coord == 0) then
!!$        do i=1,turbine_in_proc_cnt
!!$            j = turbine_in_proc_array(i)
!!$            buffer_array = 0.
!!$            call MPI_recv( buffer_array, nloc, MPI_rprec, j, 3, comm, status, ierr )
!!$            disk_avg_vels = disk_avg_vels + buffer_array
!!$            
!!$        enddo              
!!$                
!!$    elseif (turbine_in_proc) then
!!$        call MPI_send( disk_avg_vels, nloc, MPI_rprec, 0, 3, comm, ierr )
!!$    endif            
!!$    !##############################################              
!!$$endif
!!$
!!$!Coord==0 takes that info and calculates total disk force, then sends it back
!!$if (coord == 0) then           
!!$    !update epsilon for the new timestep (for cfl_dt) 
!!$        eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
!!$
!!$    !for each turbine:        
!!$        do s=1,nloc            
!!$             
!!$        !set pointers
!!$            p_u_d => wind_farm%turbine(s)%u_d   
!!$            p_u_d_T => wind_farm%turbine(s)%u_d_T   
!!$            p_f_n => wind_farm%turbine(s)%f_n                  
!!$            
!!$        !volume correction:
!!$        !since sum of ind is turbine volume/(dx*dy*dz) (not exactly 1.)
!!$            p_u_d = disk_avg_vels(s) * wind_farm%turbine(s)%vol_c
!!$    
!!$        !add this current value to the "running average" (first order relaxation)
!!$            p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d
!!$
!!$        !calculate total thrust force for each turbine  (per unit mass)
!!$        !force is normal to the surface (calc from u_d_T, normal to surface)
!!$            p_f_n = Ct_prime_05*abs(p_u_d_T)*p_u_d_T/wind_farm%turbine(s)%thk       
!!$            !write values to file                   
!!$                if (modulo (jt_total, tbase) == 0) then
!!$                call write_real_data( file_id(s), 'formatted', 5, (/total_time_dim, p_u_d, p_u_d_T, p_f_n, Ct_prime_05*(p_u_d_T*p_u_d_T*p_u_d_T)*pi/(4.*sx*sy) /))
!!$                endif 
!!$            !write force to array that will be transferred via MPI    
!!$            disk_force(s) = p_f_n
!!$            
!!$        enddo           
!!$        
!!$endif
!!$
!!$!send total disk force to the necessary procs (with turbine_in_proc==.true.)
!!$$if ($MPI)
!!$    !############################################## 5  
!!$        if (coord == 0) then
!!$          
!!$            do i=1,turbine_in_proc_cnt
!!$                j = turbine_in_proc_array(i)
!!$                call MPI_send( disk_force, nloc, MPI_rprec, j, 5, comm, ierr )
!!$            enddo              
!!$                        
!!$        elseif (turbine_in_proc) then
!!$            call MPI_recv( disk_force, nloc, MPI_rprec, 0, 5, comm, status, ierr )
!!$        endif     
!!$    !##############################################     
!!$$endif 
!!$    
!!$!apply forcing to each node
!!$if (turbine_in_proc) then
!!$
!!$    do s=1,nloc     
!!$            
!!$        do l=1,wind_farm%turbine(s)%num_nodes
!!$            i2 = wind_farm%turbine(s)%nodes(l,1)
!!$            j2 = wind_farm%turbine(s)%nodes(l,2)
!!$            k2 = wind_farm%turbine(s)%nodes(l,3)
!!$            ind2 = wind_farm%turbine(s)%ind(l)
!!$            fxaF(i2,j2,k2) = disk_force(s)*wind_farm%turbine(s)%nhat(1)*ind2 
!!$        enddo
!!$
!!$    enddo
!!$    
!!$endif   
!!$
!!$call phys2waveF_pr(fxaF, fxa) 
!!$
!!$deallocate(w_uv)
!!$
!!$!spatially average velocity at the top of the domain and write to file
!!$if (coord .eq. nproc-1) then
!!$    string1=path // 'output/vel_top_of_domain.dat'
!!$    open(unit=1,file=string1,status='unknown',form='formatted',action='write',position='append')
!!$    write(1,*) total_time, sum(u(:,:,nz-1))/(nxp*ny)
!!$    close(1)
!!$endif
!!$
!!$nullify(y,z)
!!$
!!$end subroutine turbines_forcing_fourier
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_RNL()
use param, only: ld_big, ny2, nxp, ny, L_x
use sim_param, only: u, fxa
use derivatives, only: convolve, convolve2
use derivatives, only: dft_direct_back_2d_n_yonlyC_big
use derivatives, only: dft_direct_forw_2d_n_yonlyC_big
use derivatives, only: dft_direct_back_2d_n_yonlyC
use derivatives, only: dft_direct_forw_2d_n_yonlyC
use messages
use fft
!use functions, only: interp_to_uv_grid

implicit none
!real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
!real(rprec) :: ind2
!real(rprec), dimension(nloc) :: disk_avg_vels, disk_force
!real(rprec), allocatable, dimension(:,:,:) :: w_uv ! Richard: This 3D matrix can relatively easy be prevented
real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nxp) :: xp
real(rprec), dimension(ld_big, ny2) :: fhat_big
real(rprec), dimension(ld_big, ny2) :: uhat_big
real(rprec), dimension(ld_big, ny2) :: fxa_big

real(rprec), dimension(ld, ny) :: fhat
real(rprec), dimension(ld, ny) :: uhat

complex(rprec) :: imag = (0._rprec, 1._rprec)
integer :: jx,jy,jz
real(rprec) :: r1,r2,r3,r4,r5,r6,r7,rad

rad = 0.5_rprec * dia_all

fhat_big = 0._rprec
uhat_big = 0._rprec
fxa_big = 0._rprec

fhat = 0._rprec
uhat = 0._rprec

nullify(y,z)
y => grid % y  ! non-dimensional, 0 to L_y
z => grid % z  ! goes from 0 to 1+, coord is taken into account

!do jx=1,nxp !+1
!  xp(jx) = (jx - 1)*dx    !! dx = L_x / nxp
!enddo

do jx=1,num_x
  xp(jx) = (jx - 1)*L_x/real(num_x,rprec)
enddo

!!$print*, 'Y >>'
!!$write(*,*) coord, y(:)
!!$call mpi_barrier(comm, ierr)
!!$print*, 'Z >>'
!!$write(*,*) coord, z(:)

!!$print*, 'rad =', rad
!!$do jz = 1, nz
!!$do jy = 1, ny
!!$   r = sqrt( abs(y(jy)-y(ny/2))**2 + abs(z(jz)-height_all)**2 )
!!$   if ( r .le. rad ) then
!!$      write(*,*) 'jz, jy: ', jz,jy
!!$   endif
!!$enddo
!!$enddo

fxa(:,:,:) = 0._rprec
do jz = 1, nz

!call padd( uhat_big(:,:), u(:,:,jz) )
!call dft_direct_back_2d_n_yonlyC_big( uhat_big(:,:) )

uhat(:,:) = u(:,:,jz)
call dft_direct_back_2d_n_yonlyC( uhat(:,:) )

!fhat_big(:,:) = convolve( uhat_big(:,:), uhat_big(:,:) )
fhat(:,:) = convolve2( uhat(:,:), uhat(:,:) )

!! aligned
!do jy = 1, ny !ny2
!   r1 = sqrt( abs(y(jy)-0.261799)**2 + abs(z(jz)-height_all)**2 )
!   r2 = sqrt( abs(y(jy)-0.785398)**2 + abs(z(jz)-height_all)**2 )
!   r3 = sqrt( abs(y(jy)-1.309)**2 + abs(z(jz)-height_all)**2 )
!   r4 = sqrt( abs(y(jy)-1.8326)**2 + abs(z(jz)-height_all)**2 )
!   r5 = sqrt( abs(y(jy)-2.35619)**2 + abs(z(jz)-height_all)**2 )
!   r6 = sqrt( abs(y(jy)-2.87979)**2 + abs(z(jz)-height_all)**2 )
!   if ( r1 .le. rad .or. r2 .le. rad .or. r3 .le. rad .or. r4 .le. rad .or. r5 .le. rad .or. r6 .le. rad) then
!   do jx=1,4  !4
!   !fxa_big(1:ld,jy) = fxa_big(1:ld,jy) + fhat_big(1:ld,jy)*exp( imag*kx(1:ld,jy)*xp(jx) )
!   fxa(1:ld,jy,jz) = fxa(1:ld,jy,jz) + fhat(1:ld,jy)*exp( imag*kx(1:ld,jy)*xp(jx) )
!   enddo
!   endif
!enddo

!! staggered
do jy = 1, ny !ny2
   r1 = sqrt( abs(y(jy)-0.261799)**2 + abs(z(jz)-height_all)**2 )
   r2 = sqrt( abs(y(jy)-0.785398)**2 + abs(z(jz)-height_all)**2 )
   r3 = sqrt( abs(y(jy)-1.309)**2 + abs(z(jz)-height_all)**2 )
   r4 = sqrt( abs(y(jy)-1.8326)**2 + abs(z(jz)-height_all)**2 )
   r5 = sqrt( abs(y(jy)-2.35619)**2 + abs(z(jz)-height_all)**2 )
   r6 = sqrt( abs(y(jy)-2.87979)**2 + abs(z(jz)-height_all)**2 )
   if ( r1 .le. rad .or. r2 .le. rad .or. r3 .le. rad .or. r4 .le. rad .or. r5 .le. rad .or. r6 .le. rad) then
   do jx=1,num_x  !4
   !fxa_big(1:ld,jy) = fxa_big(1:ld,jy) + fhat_big(1:ld,jy)*exp(imag*kx(1:ld,jy)*xp(jx) )
   fxa(1:ld,jy,jz) = fxa(1:ld,jy,jz) + fhat(1:ld,jy)*exp(imag*kx(1:ld,jy)*xp(jx) )
   enddo
   endif
enddo

!do jy = 1, ny !ny2
!   r1 = sqrt( abs(y(jy)-0.5235987)**2 + abs(z(jz)-height_all)**2 )
!   r2 = sqrt( abs(y(jy)-1.04719755)**2 + abs(z(jz)-height_all)**2 )
!   r3 = sqrt( abs(y(jy)-1.5707963)**2 + abs(z(jz)-height_all)**2 )
!   r4 = sqrt( abs(y(jy)-2.0943951)**2 + abs(z(jz)-height_all)**2 )
!   r5 = sqrt( abs(y(jy)-2.6179938)**2 + abs(z(jz)-height_all)**2 )
!   r6 = sqrt( abs(y(jy)-3.14159265)**2 + abs(z(jz)-height_all)**2 )
!   r7 = sqrt( abs(y(jy)-0.0)**2 + abs(z(jz)-height_all)**2 )
!   if ( r1 .le. rad .or. r2 .le. rad .or. r3 .le. rad .or. r4 .le. rad) then
!   do jx=1,2  !4
!   !fxa_big(1:ld,jy) = fxa_big(1:ld,jy) + fhat_big(1:ld,jy)*exp(imag*kx(1:ld,jy)*xp(jx) )
!   fxa(1:ld,jy,jz) = fxa(1:ld,jy,jz) + fhat(1:ld,jy)*exp(imag*kx(1:ld,jy)*xp(jx) )
!   enddo
!   endif
!   if ( r5 .le. rad .or. r6 .le. rad .or. r7 .le. rad) then
!   do jx=1,2  !4
!   !fxa_big(1:ld,jy) = fxa_big(1:ld,jy) + fhat_big(1:ld,jy)*exp(imag*kx(1:ld,jy)*xp(jx) )
!   fxa(1:ld,jy,jz) = fxa(1:ld,jy,jz) + fhat(1:ld,jy)*exp(imag*kx(1:ld,jy)*xp(jx) )
!   enddo
!   endif
!enddo

!call dft_direct_forw_2d_n_yonlyC_big( fxa_big(:,:) )
call dft_direct_forw_2d_n_yonlyC( fxa(:,:,jz) )
!call unpadd( fxa(:,:,jz), fxa_big(:,:) )
enddo

!! eliminate nonzero kx modes
fxa(3:,:,:) = 0._rprec
!fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0/L_x,rprec)   !C_Tprime = 4/3
fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0,rprec)   !C_Tprime = 4/3
!fxa(:,:,23:) = 0._rprec

end subroutine turbines_RNL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_RNL2()
use param, only: ld_big, ny2, nxp, ny, L_x
use sim_param, only: u, fxa
use derivatives, only: convolve, convolve2
use derivatives, only: dft_direct_back_2d_n_yonlyC_big
use derivatives, only: dft_direct_forw_2d_n_yonlyC_big
use derivatives, only: dft_direct_back_2d_n_yonlyC
use derivatives, only: dft_direct_forw_2d_n_yonlyC
use messages
use fft

implicit none

real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nxp) :: xp
real(rprec), dimension(ld_big, ny2) :: fhat_big
real(rprec), dimension(ld_big, ny2) :: uhat_big
real(rprec), dimension(ld_big, ny2) :: fxa_big

real(rprec), dimension(ld, ny) :: fhat
real(rprec), dimension(ld, ny) :: uhat

complex(rprec) :: imag = (0._rprec, 1._rprec)
integer :: jx,jy,jz
real(rprec) :: r1,r2,r3,r4,r5,r6,r7,rad

rad = 0.5_rprec * dia_all

fhat_big = 0._rprec
uhat_big = 0._rprec
fxa_big = 0._rprec

fhat = 0._rprec
uhat = 0._rprec

nullify(y,z)
y => grid % y  ! non-dimensional, 0 to L_y
z => grid % z  ! goes from 0 to 1+, coord is taken into account

!do jx=1,nxp !+1
!  xp(jx) = (jx - 1)*dx    !! dx = L_x / nxp
!enddo

do jx=1,num_x
  xp(jx) = (jx - 1)*L_x/real(num_x,rprec)
enddo

fxa(:,:,:) = 0._rprec
do jz = 1, nz

uhat(:,:) = u(:,:,jz)
call dft_direct_back_2d_n_yonlyC( uhat(:,:) )

fhat(:,:) = convolve2( uhat(:,:), uhat(:,:) )

do jy = 1, ny
   r1 = sqrt( abs(y(jy)-0.261799)**2 + abs(z(jz)-height_all)**2 ) ! aligned
   r2 = sqrt( abs(y(jy)-0.785398)**2 + abs(z(jz)-height_all)**2 )
   r3 = sqrt( abs(y(jy)-1.309)**2 + abs(z(jz)-height_all)**2 )
   r4 = sqrt( abs(y(jy)-1.8326)**2 + abs(z(jz)-height_all)**2 )
   r5 = sqrt( abs(y(jy)-2.35619)**2 + abs(z(jz)-height_all)**2 )
   r6 = sqrt( abs(y(jy)-2.87979)**2 + abs(z(jz)-height_all)**2 )
   !r1 = sqrt( abs(y(jy)-0.261799)**2 + abs(z(jz)-height_all)**2 ) ! staggered
   !r2 = sqrt( abs(y(jy)-0.785398)**2 + abs(z(jz)-height_all)**2 )
   !r3 = sqrt( abs(y(jy)-1.309)**2 + abs(z(jz)-height_all)**2 )
   !r4 = sqrt( abs(y(jy)-1.8326)**2 + abs(z(jz)-height_all)**2 )
   !r5 = sqrt( abs(y(jy)-2.35619)**2 + abs(z(jz)-height_all)**2 )
   !r6 = sqrt( abs(y(jy)-2.87979)**2 + abs(z(jz)-height_all)**2 )
   if ( r1 .le. rad .or. r2 .le. rad .or. r3 .le. rad .or. r4 .le. rad .or. r5 .le. rad .or. r6 .le. rad) then
       fxa(1:ld,jy,jz) = fhat(1:ld,jy) !* real(num_x, rprec)
   endif
enddo

call dft_direct_forw_2d_n_yonlyC( fxa(:,:,jz) )
enddo

!! eliminate nonzero kx modes
!fxa(3:,:,:) = 0._rprec
!fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0,rprec)   !C_Tprime = 4/3
fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0*num_x/L_x,rprec)   !C_Tprime = 4/3

end subroutine turbines_RNL2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_RNL3()
use param, only: ld_big, ny2, nxp, ny, L_x
use sim_param, only: u, fxa
use turbines_base, only: ind
use derivatives, only: convolve, convolve2
use derivatives, only: dft_direct_back_2d_n_yonlyC_big
use derivatives, only: dft_direct_forw_2d_n_yonlyC_big
use derivatives, only: dft_direct_back_2d_n_yonlyC
use derivatives, only: dft_direct_forw_2d_n_yonlyC
use messages
use fft

implicit none

real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nxp) :: xp
real(rprec), dimension(ld_big, ny2) :: fhat_big
real(rprec), dimension(ld_big, ny2) :: uhat_big
real(rprec), dimension(ld_big, ny2) :: fxa_big

real(rprec), dimension(ld, ny) :: fhat
real(rprec), dimension(ld, ny) :: uhat

complex(rprec) :: imag = (0._rprec, 1._rprec)
integer :: jx,jy,jz
real(rprec) :: r1,r2,r3,r4,r5,r6,r7,rad

rad = 0.5_rprec * dia_all

fhat_big = 0._rprec
uhat_big = 0._rprec
fxa_big = 0._rprec

fhat = 0._rprec
uhat = 0._rprec

nullify(y,z)
y => grid % y  ! non-dimensional, 0 to L_y
z => grid % z  ! goes from 0 to 1+, coord is taken into account

!do jx=1,nxp !+1
!  xp(jx) = (jx - 1)*dx    !! dx = L_x / nxp
!enddo

do jx=1,num_x
  xp(jx) = (jx - 1)*L_x/real(num_x,rprec)
enddo

fxa(:,:,:) = 0._rprec
do jz = 1, nz
  uhat(:,:) = u(:,:,jz)
  call dft_direct_back_2d_n_yonlyC( uhat(:,:) )
  fhat(:,:) = convolve2( uhat(:,:), uhat(:,:) )
  fxa(1:ld, :, jz) = fhat(1:ld, :)
enddo

do jx = 1, ld
  fxa(jx, 1:ny, 1:nz) = fxa(jx, 1:ny, 1:nz) * ind(1:ny, 1:nz) !*6._rprec
enddo

do jz = 1, nz
  call dft_direct_forw_2d_n_yonlyC( fxa(:,:,jz) )
enddo

!! eliminate nonzero kx modes
!fxa(3:,:,:) = 0._rprec
fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0/L_x,rprec)   !C_Tprime = 4/3

!if (coord .eq. 0) then
!do jz=1,nz
!write(*,*) jz, fxa(1,1:ny, jz)
!enddo
!endif

end subroutine turbines_RNL3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_RNL4()
use param, only: ld_big, ny2, nxp, ny, L_x
use sim_param, only: u, fxa !, pow
use turbines_base, only: ind
use derivatives, only: convolve, convolve2
use derivatives, only: dft_direct_back_2d_n_yonlyC_big
use derivatives, only: dft_direct_forw_2d_n_yonlyC_big
use derivatives, only: dft_direct_back_2d_n_yonlyC
use derivatives, only: dft_direct_forw_2d_n_yonlyC
use messages
use fft

implicit none

real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nxp) :: xp
real(rprec), dimension(ld_big, ny2) :: fhat_big
real(rprec), dimension(ld_big, ny2) :: uhat_big
real(rprec), dimension(ld_big, ny2) :: fxa_big

real(rprec), dimension(ld, ny) :: fhat
real(rprec), dimension(ld, ny) :: uhat

complex(rprec) :: imag = (0._rprec, 1._rprec)
integer :: jx,jy,jz
real(rprec) :: r1,r2,r3,r4,r5,r6,r7,rad

rad = 0.5_rprec * dia_all

fhat_big = 0._rprec
uhat_big = 0._rprec
fxa_big = 0._rprec

fhat = 0._rprec
uhat = 0._rprec

nullify(y,z)
y => grid % y  ! non-dimensional, 0 to L_y
z => grid % z  ! goes from 0 to 1+, coord is taken into account

!do jx=1,nxp !+1
!  xp(jx) = (jx - 1)*dx    !! dx = L_x / nxp
!enddo

do jx=1,num_x
  xp(jx) = (jx - 1)*L_x/real(num_x,rprec)
enddo

fxa(:,:,:) = 0._rprec
do jz = 1, nz
  uhat(:,:) = u(:,:,jz)
  call dft_direct_back_2d_n_yonlyC( uhat(:,:) )
  fhat(:,:) = convolve2( uhat(:,:), uhat(:,:) )
  fxa(1:ld, :, jz) = fhat(1:ld, :)
enddo

do jx = 1, ld
  fxa(jx, 1:ny, 1:nz) = fxa(jx, 1:ny, 1:nz) * ind(1:ny, 1:nz) !*3._rprec
enddo

!! eliminate nonzero kx modes
!fxa(3:,:,:) = 0._rprec
fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0/L_x,rprec)   !C_Tprime = 4/3

do jz = 1, nz
   fhat(:,:) = convolve2( fxa(1:ld,:,jz), uhat(:,:) )   !! power calc
   !pow(:,:,jz) = fhat(:,:)
enddo

do jz = 1, nz
  call dft_direct_forw_2d_n_yonlyC( fxa(:,:,jz) )
enddo


!if (coord .eq. 0) then
!do jz=1,nz
!write(*,*) jz, fxa(1,1:ny, jz)
!enddo
!endif

end subroutine turbines_RNL4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_RNL5()
use param, only: ld_big, ny2, nxp, ny, L_x, lbz, nz
use sim_param, only: u, fxa
use turbines_base, only: ind
use derivatives, only: convolve, convolve2
use derivatives, only: dft_direct_back_2d_n_yonlyC_big
use derivatives, only: dft_direct_forw_2d_n_yonlyC_big
use derivatives, only: dft_direct_back_2d_n_yonlyC
use derivatives, only: dft_direct_forw_2d_n_yonlyC
use derivatives, only: wave2phys
use messages
use fft

implicit none

real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nxp) :: xp
real(rprec), dimension(ld_big, ny2) :: fhat_big
real(rprec), dimension(ld_big, ny2) :: uhat_big
real(rprec), dimension(ld_big, ny2) :: fxa_big

real(rprec), dimension(ld, ny) :: fhat
real(rprec), dimension(ld, ny) :: uhat

real(rprec), dimension(ld,ny,lbz:nz) :: temp

complex(rprec) :: imag = (0._rprec, 1._rprec)
integer :: jx,jy,jz
real(rprec) :: r1,r2,r3,r4,r5,r6,r7,rad

rad = 0.5_rprec * dia_all

fhat_big = 0._rprec
uhat_big = 0._rprec
fxa_big = 0._rprec

fhat = 0._rprec
uhat = 0._rprec

nullify(y,z)
y => grid % y  ! non-dimensional, 0 to L_y
z => grid % z  ! goes from 0 to 1+, coord is taken into account

!do jx=1,nxp !+1
!  xp(jx) = (jx - 1)*dx    !! dx = L_x / nxp
!enddo

do jx=1,num_x
  xp(jx) = (jx - 1)*L_x/real(num_x,rprec)
enddo

call disk_avg_kx()

fxa(:,:,:) = 0._rprec

if (turbine_in_proc) then
   do s=1,nloc
      temp(:,:,:) = 0._rprec    !! temp vel field
      do l=1,wind_farm%turbine(s)%num_nodes
         jy = wind_farm%turbine(s)%nodes(l,2)
         jz = wind_farm%turbine(s)%nodes(l,3)
         temp(1:ld,jy,jz) = wind_farm % turbine(s) % u_d_kx(1:ld)
      enddo
      do jz=1,nz
         uhat(:,:) = temp(:,:,jz)
         uhat(3:,:) = 0._rprec
         fhat(:,:) = convolve2( uhat(:,:), uhat(:,:) )
         fxa(1:ld, :, jz) = fxa(1:ld,:,jz) + fhat(1:ld, :)
      enddo
   enddo
endif

do jz = 1, nz
  call dft_direct_forw_2d_n_yonlyC( fxa(:,:,jz) )
enddo

!! eliminate nonzero kx modes
!fxa(3:,:,:) = 0._rprec
!C_Tprime = 4/3
fxa(:,:,:) = fxa(:,:,:) * real(-1.0/2.0*4.0/3.0/L_x,rprec) !* 0.1_rprec   

if (coord == 0) then
do s=1,nloc
!print*, coord, s, '>>', wind_farm % turbine(s) % u_d_kx(:) 
do jy=1,ny
   temp(:,jy,1) = wind_farm % turbine(s) % u_d_kx(:) 
enddo

!!$temp(3:,:,:) = 0._rprec
!!$if ( s == 1 ) then
!!$   print*, 'a1', temp(:,1,1)
!!$   fhat(:,:) = convolve2( temp(:,:,1), temp(:,:,1) )
!!$   print*, 'b1', fhat(:,1)
!!$   fhat(:,:) = convolve2(   fhat(:,:), temp(:,:,1) )
!!$   print*, 'c1', fhat(:,1)
!!$endif

call dft_direct_forw_2d_n_yonlyC( fhat(:,:) )
do jz=1,nz
temp(:,:,jz) = fhat(:,:)
enddo
temp(3:,:,:) = 0._rprec
call wave2phys( temp )
!print*,'>>> s: ', s
!print*, 'aaa', temp(:,1,1)
!print*, 'bbb', temp(:,2,2)

enddo
endif




end subroutine turbines_RNL5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************************
subroutine place_turbines
!*******************************************************************************
!
! This subroutine places the turbines on the domain. It also sets the values for
! each individual turbine. After the subroutine is called, the following values 
! are set for each turbine in wind_farm: xloc, yloc, height, dia, thk, theta1, 
! theta2, and Ct_prime.
!
use param, only: pi, z_i
use open_file_fid_mod
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

real(rprec) :: sxx, syy, shift_base, const
real(rprec) :: dummy, dummy2
logical :: exst
integer :: fid

! Read parameters from file if needed
if (read_param) then
    ! Check if file exists and open
    inquire (file = param_dat, exist = exst)
    if (.not. exst) then
        call error (sub_name, 'file ' // param_dat // 'does not exist')
    end if

    ! Check that there are enough lines from which to read data
    nloc = count_lines(param_dat)
    !if (nloc < num_x*num_y) then
    !    nloc = num_x*num_y
    !    call error(sub_name, param_dat // 'must have num_x*num_y lines')
    !else if (nloc > num_x*num_y) then
    !    call warn(sub_name, param_dat // ' has more than num_x*num_y lines. '  &
    !              // 'Only reading first num_x*num_y lines')
    !end if
    write(*,*) "Number of turbines (nloc):", nloc

    ! Read from parameters file, which should be in this format:
    ! xloc [meters], yloc [meters], height [meters], dia [meters], thk [meters], 
    ! theta1 [degrees], theta2 [degrees], Ct_prime [-]
    write(*,*) "Reading from", param_dat
    fid = open_file_fid(param_dat, 'rewind', 'formatted')
    do k = 1, nloc
        read(fid,*) wind_farm%turbine(k)%xloc, wind_farm%turbine(k)%yloc,      &
            wind_farm%turbine(k)%height, wind_farm%turbine(k)%dia,             &
            wind_farm%turbine(k)%thk, wind_farm%turbine(k)%theta1,             &
            wind_farm%turbine(k)%theta2, wind_farm%turbine(k)%Ct_prime
    end do
    close(fid)
    
    ! Make lengths dimensionless
    do k = 1, nloc
        wind_farm%turbine(k)%xloc = wind_farm%turbine(k)%xloc / z_i
        wind_farm%turbine(k)%yloc = wind_farm%turbine(k)%yloc / z_i
        wind_farm%turbine(k)%height = wind_farm%turbine(k)%height / z_i
        wind_farm%turbine(k)%dia = wind_farm%turbine(k)%dia / z_i
        wind_farm%turbine(k)%thk = wind_farm%turbine(k)%thk / z_i
    end do
else
    ! Set values for each turbine based on values in input file
    wind_farm%turbine(:)%height = height_all
    wind_farm%turbine(:)%dia = dia_all
    wind_farm%turbine(:)%thk = thk_all
    wind_farm%turbine(:)%theta1 = theta1_all
    wind_farm%turbine(:)%theta2 = theta2_all
    wind_farm%turbine(:)%Ct_prime = Ct_prime

    ! Set baseline locations (evenly spaced, not staggered aka aligned)
    k = 1
    sxx = sx * dia_all  ! x-spacing with units to match those of L_x
    syy = sy * dia_all  ! y-spacing
    do i = 1,num_x
        do j = 1,num_y
            wind_farm%turbine(k)%xloc = sxx*real(2*i-1)/2
            wind_farm%turbine(k)%yloc = syy*real(2*j-1)/2
            k = k + 1
        end do
    end do

 
    ! Place turbines based on orientation flag
    ! This will shift the placement relative to the baseline locations abive
    select case (orientation)
        ! Evenly-spaced, not staggered
        case (1)
    
        ! Evenly-spaced, horizontally staggered only
        ! Shift each row according to stag_perc    
        case (2)
            do i = 2, num_x
                do k = 1+num_y*(i-1), num_y*i
                    shift_base = syy * stag_perc/100.
                    wind_farm%turbine(k)%yloc = mod( wind_farm%turbine(k)%yloc &
                                                    + (i-1)*shift_base , L_y )
                end do
            end do
 
        ! Evenly-spaced, only vertically staggered (by rows)
        case (3)
            ! Make even rows taller
            do i = 2, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.+stag_perc/100.)
                end do
            end do
            ! Make odd rows shorter
            do i = 1, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.-stag_perc/100.)
                end do
            end do
       
        ! Evenly-spaced, only vertically staggered, checkerboard pattern
        case (4)
            k = 1
            do i = 1, num_x 
                do j = 1, num_y
                    ! this should alternate between 1, -1
                    const = 2.*mod(real(i+j),2.)-1.  
                    wind_farm%turbine(k)%height = height_all                   &
                                                  *(1.+const*stag_perc/100.)
                    k = k + 1
                end do
            end do

        ! Aligned, but shifted forward for efficient use of simulation space 
        ! during CPS runs
        case (5)
        ! Shift in spanwise direction: Note that stag_perc is now used
            k=1
            dummy=stag_perc                                                    &
                  *(wind_farm%turbine(2)%yloc - wind_farm%turbine(1)%yloc)
            do i = 1, num_x
                do j = 1, num_y
                    dummy2=dummy*(i-1)         
                    wind_farm%turbine(k)%yloc=mod( wind_farm%turbine(k)%yloc   &
                                                  + dummy2,L_y)
                    k=k+1
                end do
            end do      

        case default
            call error (sub_name, 'invalid orientation')
        
    end select
end if
end subroutine place_turbines

!*******************************************************************************
function count_lines(fname) result(N)
!*******************************************************************************
!
! This function counts the number of lines in a file
!
use open_file_fid_mod
use messages
use param, only : CHAR_BUFF_LENGTH
implicit none
character(*), intent(in) :: fname
logical :: exst
integer :: fid, ios
integer :: N

character(*), parameter :: sub_name = mod_name // '.count_lines'

! Check if file exists and open
inquire (file = trim(fname), exist = exst)
if (.not. exst) then
    call error (sub_name, 'file ' // trim(fname) // 'does not exist')
end if
fid = open_file_fid(trim(fname), 'rewind', 'formatted')

! count number of lines and close
ios = 0
N = 0
do
    read(fid, *, IOstat = ios)
    if (ios /= 0) exit
    N = N + 1
end do

! Close file
close(fid)

end function count_lines

end module turbines
