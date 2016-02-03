module nacelle
!! Module that introduces nacelle in actuator line simulations.
!! The function 'nacelle_init' should be called before time integration starts. This function looks up where the nacelle is
!! placed and builds up an indicator function that is used to distribute the forces later. The relevant parameters are given
!! at the beginning of the function 'nacelle_init'.
!! The function 'nacelle_forcing' should be called during time integration. Note that the code uses the variable fxa to transfer
!! the forces to lesgo. This variable must thus be made available in the code.
use param, only: L_x,L_y,L_z,z_i,dx,dy,dz,nx,ny,nz,nz_tot,pi,comm,ierr,coord,MPI_RPREC,MPI_SUM,dt_dim,jt_total,path,total_time_dim
use grid_defs, only: grid 
use types, only : rprec
!use turbines_base, only: wind_farm,alpha,trunc,Ct_prime,filter_cutoff,T_avg_dim,tbase,turbine_cumulative_time,nacelle_fraction,nloc
use string_util

use atm_input_util, only : turbineArray, turbineModel, numberOfTurbines

implicit none
!include 'tecryte.h'
real(rprec), dimension(:), allocatable :: z_tot                  !! vertical grid in LESGOi

!!! For each nacelle we have the following properties
!! First store in temporary array to make sure that we can allocate the right memory size for these arrays, which are unknow in advance
real(rprec), dimension(:,:,:), allocatable :: nacelle_nodes_temp !! The computational nodes located in the nacelle 
real(rprec), dimension(:,:), allocatable :: nacelle_ind_temp     !! indicator function of the nacelle
real(rprec), dimension(:,:,:), allocatable :: nacelle_nodes      !! The computational nodes located in the nacelle 
real(rprec), dimension(:,:), allocatable :: nacelle_ind          !! indicaton function of the nacelle
real(rprec), dimension(:), allocatable :: nacelle_xloc           !! streamwise location
real(rprec), dimension(:), allocatable :: nacelle_yloc           !! spanwise location
real(rprec), dimension(:), allocatable :: nacelle_height         !! height of nacelle center
real(rprec), dimension(:), allocatable :: nacelle_dia            !! nacelle diameter
real(rprec), dimension(:), allocatable :: nacelle_r              !! nacelle radius
real(rprec), dimension(:), allocatable :: nacelle_thick          !! nacelle tickness (in streamwise direction)
real(rprec), dimension(:), allocatable :: nacelle_vol            !! volume of nacelle disk
real(rprec), dimension(:), allocatable :: nacelle_theta          !! angle of nacelle
real(rprec), dimension(:), allocatable :: nacelle_Cd             !! Drag coefficient of the nacelle
real(rprec), dimension(:), allocatable :: nacelle_nhat1          !! x1 direction
real(rprec), dimension(:), allocatable :: nacelle_nhat2          !! x2 direction
real(rprec), dimension(:), allocatable :: nacelle_nhat3          !! x3 direction
real(rprec), dimension(:), allocatable :: nacelle_u_avg          !! time averaged nacelle velocity
integer    , dimension(:), allocatable :: nacelle_num_nodes      !! Indicates the actual counted number of grid points in the nacelle area
integer    , dimension(:), allocatable :: nacelle_file_id        !! Used for output files
integer    , dimension(:), allocatable :: nacelle_min_i,nacelle_max_i,nacelle_min_j !! counters
integer    , dimension(:), allocatable :: nacelle_max_j,nacelle_min_k,nacelle_max_k !! counters
!! End define for each nacelle seperately
integer :: nacelle_trunc !! fixed parameter
real(rprec) :: nacelle_alpha !! indicate howmuch Gaussian kernel is smoothed

real(rprec) ::   filter_cutoff = 1e-2                            !! Filtering out small forces

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nacelle_init()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
character (100) :: string1
logical :: exst
integer :: k,ss
allocate(z_tot(nz_tot))

!z_tot for total domain (since z is local to the processor)
do k=1,nz_tot
   z_tot(k) = (k - 0.5_rprec) * dz
enddo

allocate(nacelle_xloc(numberOfTurbines))
allocate(nacelle_yloc(numberOfTurbines))
allocate(nacelle_height(numberOfTurbines))
allocate(nacelle_dia(numberOfTurbines))
allocate(nacelle_r(numberOfTurbines))
allocate(nacelle_thick(numberOfTurbines))
allocate(nacelle_vol(numberOfTurbines))
allocate(nacelle_theta(numberOfTurbines))
allocate(nacelle_Cd(numberOfTurbines))
allocate(nacelle_nhat1(numberOfTurbines))
allocate(nacelle_nhat2(numberOfTurbines))
allocate(nacelle_nhat3(numberOfTurbines))
allocate(nacelle_u_avg(numberOfTurbines))
allocate(nacelle_num_nodes(numberOfTurbines))
allocate(nacelle_file_id(numberOfTurbines))
allocate(nacelle_min_i(numberOfTurbines))
allocate(nacelle_min_j(numberOfTurbines))
allocate(nacelle_min_k(numberOfTurbines))
allocate(nacelle_max_i(numberOfTurbines))
allocate(nacelle_max_j(numberOfTurbines))
allocate(nacelle_max_k(numberOfTurbines))

do ss=1,numberOfTurbines
!! Properties of the nacelle should be set here 
nacelle_Cd(ss)     =4._rprec                    ! Drag coefficient of nacelle
!nacelle_xloc(ss)   =wind_farm%turbine(ss)%xloc   ! streamwise location of nacelle (in non dimensional units) 
!nacelle_yloc(ss)   =wind_farm%turbine(ss)%yloc   ! spanwise location of nacelle (in non dimensional units)
!nacelle_height(ss) =wind_farm%turbine(ss)%height ! vertical height of nacelle center (in non dimensional units)
!nacelle_dia(ss)    =wind_farm%turbine(ss)%dia*nacelle_fraction    ! diameter of nacelle (in non dimensional units)
!nacelle_thick(ss)  =wind_farm%turbine(ss)%thk    ! tickness of nacelle (in non dimensional units)
!nacelle_theta(ss)  =wind_farm%turbine(ss)%theta1 ! angle of the nacelle with respect to streamwise direction
nacelle_xloc(ss)   =turbineArray(ss) % nacelleLocation(1) / z_i   ! streamwise location of nacelle (in non dimensional units) 
nacelle_yloc(ss)   =turbineArray(ss) % nacelleLocation(2) / z_i   ! spanwise location of nacelle (in non dimensional units)
nacelle_height(ss) =turbineModel(ss) % TowerHt / z_i ! vertical height of nacelle center (in non dimensional units)
nacelle_dia(ss)    =turbineModel(ss) % HubRad * 2. / z_i    ! diameter of nacelle (in non dimensional units)
nacelle_thick(ss)  = 1.01 * dx    ! tickness of nacelle (in non dimensional units)
nacelle_theta(ss)  = 0. ! angle of the nacelle with respect to streamwise direction
!nacelle_alpha  =alpha                       ! parameter indicating the smoothing of the force
!nacelle_trunc  =trunc                       ! parameter indicating over how many grid points the force is smoothed
nacelle_alpha  = 1.3                       ! parameter indicating the smoothing of the force
nacelle_trunc  = 3                       ! parameter indicating over how many grid points the force is smoothed

!nacelle_thick(ss)  = max(nacelle_thick(ss),dx) ! reset nacelle tickness when dx is larger than nacelle
nacelle_vol(ss)    = pi/4.*(nacelle_dia(ss))**2*nacelle_thick(ss) ! calculate the nacelle volume

if (coord==0) then
write(*,*) '------ nacelle input ------'
write(*,*) 'nacelle_Cd', nacelle_Cd(ss)
write(*,*) 'nacelle_xloc', nacelle_xloc(ss)   
write(*,*) 'nacelle_yloc', nacelle_yloc(ss)
write(*,*) 'nacelle_height', nacelle_height(ss)
write(*,*) 'nacelle_dia', nacelle_dia(ss)
write(*,*) 'nacelle_thick', nacelle_thick(ss)
write(*,*) 'nacelle_theta', nacelle_theta(ss)
write(*,*) 'nacelle_alpha', nacelle_alpha
write(*,*) 'nacelle_trunc', nacelle_trunc
endif
enddo

!! Finds the grid cells that are located in the nacelle area
call nacelle_nodes_func

!! Makes an indicate function for the nacelle based on the poinst in the nacelle area. The indicator function is
!! used to calculate the forces that need to be applied in the nacelle area. The magnitude of the indicate function
!! indicates how much of the force is applied in a particular gridpoint. At the moment use a Gaussian smoothing kernel.
call nacelle_filter_ind()

!if (turbine_cumulative_time) then
!     string1 = path // 'nacelle/nacelle_u_d_T.dat'
!     inquire (file=string1,exist=exst)
!     if (exst) then
!        write(*,*) 'Reading from file nacelle_u_d_T.dat'
!        open(1,file=string1)
!        do ss=1,nloc
!        read(1,*) nacelle_u_avg(ss)
!        enddo
!        close(1)
!     else
!     do ss=1,nloc
!     nacelle_u_avg(ss) = -7._rprec
!     enddo
!     endif   
!else
!do ss=1,nloc
!nacelle_u_avg(ss) = -7._rprec
!enddo
!endif

! Generate the files for the turbine forcing output
!if(coord==0) then
!do ss=1,nloc
!call string_splice( string1, path // 'nacelle/nacelle_', ss, '_forcing.dat' )
!nacelle_file_id(ss) = open_file( string1, 'append', 'formatted' )
!enddo
!endif

end subroutine nacelle_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nacelle_nodes_func
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(rprec),pointer,dimension(:)::x,y,z
real(rprec) :: rx,ry,rz,r,r_norm,r_disk
integer :: i,j,k,i2,j2
integer :: imax,jmax,kmax,count_i,count_n,icp,jcp,kcp,ss,buffer

!! Set the grid variables
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

do ss=1,numberOfTurbines
nacelle_num_nodes(ss) = 0  !number of nodes that are influenced by nacelle
count_n = 0            !used for counting nodes for each nacelle 
count_i = 1            !index count - used for writing to array "nodes"

!!identify "search area"
nacelle_r(ss)  = nacelle_dia(ss)/2.0d0 !! set the radius of the nacelle
imax       = nacelle_r(ss)/dx + 2  !! number of streamwise points to be searched for nacelle radius (+2 is safety) 
jmax       = nacelle_r(ss)/dy + 2  !! number of spanwise points to be searched for nacelle radius (+2 is safety)
kmax       = nacelle_r(ss)/dz + 2  !! number of vertical points to be searched for nacelle radius (+2 is safety)

!!determine unit normal vector with respect to nacelle
nacelle_nhat1(ss) = -cos(pi*nacelle_theta(ss)/180.)
nacelle_nhat2(ss) = -sin(pi*nacelle_theta(ss)/180.)
nacelle_nhat3(ss) =  0.0d0

!!determine nearest (i,j,k) to center nacelle
icp = nint(nacelle_xloc(ss)/dx)
jcp = nint(nacelle_yloc(ss)/dy)
kcp = nint(nacelle_height(ss)/dz + 0.5)

!!determine limits for checking i,j,k
nacelle_min_i(ss) = icp-imax
nacelle_max_i(ss) = icp+imax
nacelle_min_j(ss) = jcp-jmax
nacelle_max_j(ss) = jcp+jmax
nacelle_min_k(ss) = max((kcp-kmax),1)
nacelle_max_k(ss) = min((kcp+kmax),nz_tot)
    
!!check grid points around the nacelle center.	
do k=nacelle_min_k(ss),nacelle_max_k(ss)
   do j=nacelle_min_j(ss),nacelle_max_j(ss)
      do i=nacelle_min_i(ss),nacelle_max_i(ss)
      !!vector from center point to this node is (rx,ry,rz) with length r
      !!if statements are used to handle periodic boundary conditions
      if (i<1) then
      i2 = mod(i+nx-1,nx)+1
      rx = (x(i2)-L_x) - nacelle_xloc(ss)
      elseif (i>nx) then
      i2 = mod(i+nx-1,nx)+1
      rx = (L_x+x(i2)) - nacelle_xloc(ss)
      else
      i2 = i
      rx = x(i) - nacelle_xloc(ss)
      endif            
      if (j<1) then
      j2 = mod(j+ny-1,ny)+1
      ry = (y(j2)-L_y) - nacelle_yloc(ss)                
      elseif (j>ny) then
      j2 = mod(j+ny-1,ny)+1
      ry = (L_y+y(j2)) - nacelle_yloc(ss)
      else
      j2 = j
      ry = y(j) - nacelle_yloc(ss)
      endif                      
      rz = z_tot(k) - nacelle_height(ss)
      r = sqrt(rx*rx + ry*ry + rz*rz)
      !length projected onto unit normal for nacelle
      r_norm = abs(rx*nacelle_nhat1(ss) + ry*nacelle_nhat2(ss) + rz*nacelle_nhat3(ss))
      !(remaining) length projected onto nacelle
      r_disk = sqrt(r*r - r_norm*r_norm)
      !if r_disk<R_t and r_norm within thk/2 from nacelle than in nacelle
      if ( (abs(ry) .LE. 0.5*dy) .AND. (abs(rx) .LE. 0.5*dx ) .AND. abs(rz) .LE. nacelle_r(ss)) then ! RS makes sure there is always a tower
!      if ( (r_disk .LE. nacelle_r(ss)) .AND. (r_norm .LE. nacelle_thick(ss)/2.) ) then
!         nacelle_ind(count_i)     = 1._rprec 
!         nacelle_nodes(count_i,1) = i2
!         nacelle_nodes(count_i,2) = j2
!         nacelle_nodes(count_i,3) = k   !global k (might be out of this proc's range)
         count_n = count_n + 1
         count_i = count_i + 1
      endif
      enddo ! end i
   enddo ! end j 
enddo ! end k

!! allocate here
nacelle_num_nodes(ss) = count_n !! set num_nodes equal to the number of grid points found in the nacelle area
enddo

buffer=maxval(nacelle_num_nodes)
allocate(nacelle_ind_temp(buffer,numberOfTurbines))
allocate(nacelle_nodes_temp(buffer,3,numberOfTurbines))

do ss=1,numberOfTurbines
count_n = 0            !used for counting nodes for each nacelle 
count_i = 1            !index count - used for writing to array "nodes"

!!check grid points around the nacelle center.  
do k=nacelle_min_k(ss),nacelle_max_k(ss)
   do j=nacelle_min_j(ss),nacelle_max_j(ss)
      do i=nacelle_min_i(ss),nacelle_max_i(ss)
      !!vector from center point to this node is (rx,ry,rz) with length r
      !!if statements are used to handle periodic boundary conditions
      if (i<1) then
      i2 = mod(i+nx-1,nx)+1
      rx = (x(i2)-L_x) - nacelle_xloc(ss)
      elseif (i>nx) then
      i2 = mod(i+nx-1,nx)+1
      rx = (L_x+x(i2)) - nacelle_xloc(ss)
      else
      i2 = i
      rx = x(i) - nacelle_xloc(ss)
      endif
      if (j<1) then
      j2 = mod(j+ny-1,ny)+1
      ry = (y(j2)-L_y) - nacelle_yloc(ss)
      elseif (j>ny) then
      j2 = mod(j+ny-1,ny)+1
      ry = (L_y+y(j2)) - nacelle_yloc(ss)
      else
      j2 = j
      ry = y(j) - nacelle_yloc(ss)
      endif
      rz = z_tot(k) - nacelle_height(ss)
      r = sqrt(rx*rx + ry*ry + rz*rz)
      !length projected onto unit normal for nacelle
      r_norm = abs(rx*nacelle_nhat1(ss) + ry*nacelle_nhat2(ss) + rz*nacelle_nhat3(ss))
      !(remaining) length projected onto nacelle
      r_disk = sqrt(r*r - r_norm*r_norm)
      !if r_disk<R_t and r_norm within thk/2 from nacelle than in nacelle
      if ( (abs(ry) .LE. 0.5*dy) .AND. (abs(rx) .LE. 0.5*dx ) .AND. abs(rz) .LE. nacelle_r(ss)) then ! RS makes sure there is always a tower
      !if ( (r_disk .LE. nacelle_r(ss)) .AND. (r_norm .LE. nacelle_thick(ss)/2.) ) then
         nacelle_ind_temp(count_i,ss)     = 1._rprec
         nacelle_nodes_temp(count_i,1,ss) = i2
         nacelle_nodes_temp(count_i,2,ss) = j2
         nacelle_nodes_temp(count_i,3,ss) = k   !global k (might be out of this proc's range)
         count_n = count_n + 1
         count_i = count_i + 1
      endif
      enddo ! end i
   enddo ! end j 
enddo ! end k

if (coord==0) then
write(*,*) 'nacelle_num_nodes',coord,nacelle_num_nodes(ss)
endif
!! Check the conditions in the two loops above when this message is shown
if (count_n.ne.nacelle_num_nodes(ss)) write(*,*) 'something is going wrong above 1' 
enddo

nullify(x,y,z)

end subroutine nacelle_nodes_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nacelle_filter_ind()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(rprec), dimension(:,:,:,:), allocatable :: out_a
real(rprec), dimension(:,:,:), allocatable :: g,temp_array
real(rprec) :: sumA,sumG,delta2,r2,dxdydz
integer :: i,j,k,i2,j2,k2,i3,j3,i4,j4,l,ssx,ssy,ssz
integer :: count_i,count_n
integer :: k_start,k_end,ss,nz_nacelle,buffer

dxdydz=dx*dy*dz
nz_nacelle=nz_tot*(maxval(nacelle_height)+maxval(nacelle_dia))+2*nacelle_trunc+3
if (coord==0) write(*,*) 'nz_nacelle',nz_nacelle

!! Allocation of temporary arrays. 
!! It should be possible to optimize memory usage here. It is only used at the startup of the simulation
allocate(out_a(nx,ny,nz_nacelle,numberOfTurbines))
allocate(temp_array(nx,ny,nz_nacelle))

!! Gaussian smoothing kernel, centered around middle grid point. Matrix made smaller to limit memory use. Can probably be limited more
allocate(g(nx/4:nx/4*3,ny/4:ny/4*3,nz_tot/4:nz_tot/4*3))
delta2 = nacelle_alpha**2 * (dx**2 + dy**2 + dz**2)
do k=nz_tot/4,nz_tot/4*3
   do j=ny/4,ny/4*3
      do i=nx/4,nx/4*3
      r2 = ((real(i)-nx/2.)*dx)**2 + ((real(j)-ny/2.)*dy)**2 + ((real(k)-(nz_tot-1)/2.)*dz)**2
      g(i,j,k) = sqrt(6./(pi*delta2))*6./(pi*delta2)*exp(-6.*r2/delta2)
      enddo
   enddo
enddo

!!normalize the convolution function
sumG = sum(g(:,:,:))*dxdydz
g = g/sumG

!! perform convolution on temp_array to out_a    
out_a=0.

do ss=1,numberOfTurbines   
!! create the input array (nx,ny,nz_tot) from a list of included nodes
temp_array = 0.
do l=1,nacelle_num_nodes(ss)
   i2                   = nacelle_nodes_temp(l,1,ss)
   j2                   = nacelle_nodes_temp(l,2,ss)
   k2                   = nacelle_nodes_temp(l,3,ss)
   temp_array(i2,j2,k2) = nacelle_ind_temp(l,ss)
enddo

!!convolution computed for points (i4,j4,k)
!!only compute for nodes near the nacelle (possible influence region indicated by trunc)
do k=max(nacelle_min_k(ss)-nacelle_trunc,1),min(nacelle_max_k(ss)+nacelle_trunc,nz_tot)
do j=(nacelle_min_j(ss)-nacelle_trunc),(nacelle_max_j(ss)+nacelle_trunc)
do i=(nacelle_min_i(ss)-nacelle_trunc),(nacelle_max_i(ss)+nacelle_trunc)
        
   i4 = mod(i+nx-1,nx)+1       !since values may be out 1-nx,1-ny domain (spectral BCs)
   j4 = mod(j+ny-1,ny)+1              
        
   !!for each (i4,j4,k), center convolution function on that point and 'integrate' 
   !!relative coords are (ssx,ssy,ssz). absolute coords of other/surrounding points are (i2,j2,k2)
   !!only need to consider other/surrounding points near (i4,j4,k) since conv. function is compact

   do k2=max(k-nacelle_trunc,1),min(k+nacelle_trunc,nz_tot)     !currently using truncated Gaussian
   do j2=j-nacelle_trunc,j+nacelle_trunc
   do i2=i-nacelle_trunc,i+nacelle_trunc

      i3 = mod(i2+nx-1,nx)+1      !since values may be out 1-nx,1-ny domain (spectral BCs)
      j3 = mod(j2+ny-1,ny)+1             
          
      !! Look up position with respect to center of smoothing kernel
      ssx = mod(i2-i+nx/2+nx-1,nx)+1
      ssy = mod(j2-j+ny/2+ny-1,ny)+1       
      ssz = k2-k+(nz_tot-1)/2       !since no spectral BCs in z-direction
                            
      out_a(i4,j4,k,ss) = out_a(i4,j4,k,ss) + temp_array(i3,j3,k2)*g(ssx,ssy,ssz)*dxdydz
   enddo
   enddo
   enddo
enddo
enddo
enddo

! To make sure we use the right value for the normalization
do k=1,nz_nacelle  !global k
do j=1,ny
do i=1,nx
   if (out_a(i,j,k,ss) < filter_cutoff) then
      out_a(i,j,k,ss)=0._rprec
   endif
enddo
enddo
enddo

!!normalize this "indicator function" 
sumA = sum(out_a(1:nx,1:ny,1:nz_nacelle,ss))*dxdydz
out_a(1:nx,1:ny,1:nz_nacelle,ss)= nacelle_vol(ss)/sumA*out_a(1:nx,1:ny,1:nz_nacelle,ss)

!!select the nodes that belong to this core 
$if ($MPI)
   k_start = 1+coord*(nz-1)
   k_end = nz-1+coord*(nz-1)
$else
   k_start = 1
   k_end = nz
$endif

count_n=0 !! reset counter
count_i=1 !! reset counter
do k=k_start,min(k_end,nz_nacelle)  !global k
do j=1,ny
do i=1,nx
   !! store all points in which the indicator function has some value. This is the area of the nacelle where the forces should be applied
   if (out_a(i,j,k,ss) > filter_cutoff) then
!      nacelle_ind  (count_i)   = out_a(i,j,k,ss)
!      nacelle_nodes(count_i,1) = i
!      nacelle_nodes(count_i,2) = j
!      nacelle_nodes(count_i,3) = k - coord*(nz-1)  !local k
      count_n = count_n + 1
      count_i = count_i + 1
   endif
enddo
enddo
enddo

nacelle_num_nodes(ss)=count_n
enddo

!! allocate here
buffer=maxval(nacelle_num_nodes)
allocate(nacelle_ind(buffer,numberOfTurbines))
allocate(nacelle_nodes(buffer,3,numberOfTurbines))

do ss=1,numberOfTurbines
count_n=0 !! reset counter
count_i=1 !! reset counter
do k=k_start,min(k_end,nz_nacelle)  !global k
do j=1,ny
do i=1,nx
   !! store all points in which the indicator function has some value. This is the area of the nacelle where the forces should be applied
   if (out_a(i,j,k,ss) > filter_cutoff) then
      nacelle_ind  (count_i,ss)   = out_a(i,j,k,ss)
      nacelle_nodes(count_i,1,ss) = i
      nacelle_nodes(count_i,2,ss) = j
      nacelle_nodes(count_i,3,ss) = k - coord*(nz-1)  !local k
      count_n = count_n + 1
      count_i = count_i + 1
   endif
enddo
enddo
enddo

if (nacelle_num_nodes(ss)>0) then
write(*,*) 'nacelle_num_nodes',coord,nacelle_num_nodes(ss)
endif
!! Check the conditions in the two loops above when this message is shown
if (count_n.ne.nacelle_num_nodes(ss)) write(*,*) 'something is going wrong above 2' 
enddo

deallocate(out_a,g,temp_array)
deallocate(nacelle_ind_temp,nacelle_nodes_temp)

end subroutine nacelle_filter_ind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nacelle_forcing()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use sim_param, only: u,fxa
implicit none
real(rprec)::nacelle_avg_vels,nacelle_avg_vels2,nacelle_force,u_nacelle
integer :: i2,j2,k2,l,ss

do ss=1,numberOfTurbines
nacelle_avg_vels = 0.0d0 !! average velocity in the nacelle area
nacelle_force    = 0.0d0 !! force that should be applied in the nacelle area

!! Each core calculates the wheighted average of the velocity of grid points in the nacelle area
!! For cores without nodes num_nodes=0, and this loop is not computed 
do l=1,nacelle_num_nodes(ss)  
   i2 = nacelle_nodes(l,1,ss)
   j2 = nacelle_nodes(l,2,ss)
   k2 = nacelle_nodes(l,3,ss)
   nacelle_avg_vels = nacelle_avg_vels-u(i2,j2,k2)*nacelle_ind(l,ss)        
enddo

! Determine the sum over all processors to get the total weighted sum
$if ($MPI) 
call mpi_allreduce(nacelle_avg_vels,nacelle_avg_vels2, 1, MPI_RPREC, MPI_SUM, comm, ierr)
$endif

!volume correction
u_nacelle = nacelle_avg_vels2 * dx*dy*dz/nacelle_vol(ss) 
!calculate total thrust force for each nacelle  (per unit mass)
nacelle_force = 0.5_rprec*nacelle_Cd(ss)*abs(u_nacelle)*u_nacelle/nacelle_thick(ss)

!if (coord==0) then
!if (modulo (jt_total, tbase) == 0) then
!call write_real_data( nacelle_file_id(ss), 'formatted', 4, (/total_time_dim, u_nacelle, nacelle_u_avg(ss), nacelle_force/))
!endif
!endif 

!! Each core applies the forces to nodes in its domain based on the force calculated above and the fraction
!! of the force that should be applied per point (indicator function) 
!! For cores without nodes num_nodes=0, and this loop is not computed 
do l=1,nacelle_num_nodes(ss)
   i2 = nacelle_nodes(l,1,ss)
   j2 = nacelle_nodes(l,2,ss)
   k2 = nacelle_nodes(l,3,ss)
   fxa(i2,j2,k2) = fxa(i2,j2,k2)+nacelle_force*nacelle_ind(l,ss)
enddo
enddo

end subroutine nacelle_forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine nacelle_finalize ()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!implicit none
!character (100) :: string1
!integer :: ss

!write disk-averaged velocity to file along with T_avg_dim
!if (coord == 0) then
!    string1 = path // 'nacelle/nacelle_u_d_T.dat'
!    open (unit=1,file = string1, status='unknown',form='formatted', action='write',position='rewind')
!    do ss=1,numberOfTurbines
!    write(1,*) nacelle_u_avg(ss)
!    enddo
!    close (1)
!endif
!
!end subroutine nacelle_finalize

end module nacelle
