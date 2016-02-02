module tower
!! Module that introduces tower in actuator line simulations.
!! The function 'tower_init' should be called before time integration starts. This function looks up where the tower is
!! placed and builds up an indicator function that is used to distribute the forces later. The relevant parameters are given
!! at the beginning of the function 'tower_init'.
!! The function 'tower_forcing' should be called during time integration. Note that the code uses the variable fxa to transfer
!! the forces to lesgo. This variable must thus be made available in the code.
use param, only: L_x,L_y,L_z,z_i,dx,dy,dz,nx,ny,nz,nz_tot,pi,comm,ierr,coord,MPI_RPREC,MPI_SUM,dt_dim,jt_total,path,total_time_dim
use grid_defs, only: grid 
use types, only : rprec
!use turbines_base, only: wind_farm,alpha,trunc,Ct_prime,filter_cutoff,T_avg_dim,tbase,turbine_cumulative_time,tower_fraction,nloc
use nacelle, only: nacelle_dia 
!use string_util
use atm_input_util, only : turbineArray, turbineModel, numberOfTurbines

implicit none
!include 'tecryte.h'
real(rprec), dimension(:), allocatable :: z_tot                  !! vertical grid in LESGOi

!!! For each tower we have the following properties
!! First store in temporary array to make sure that we can allocate the right memory size for these arrays, which are unknow in advance
real(rprec), dimension(:,:,:), allocatable :: tower_nodes_temp !! The computational nodes located in the tower 
real(rprec), dimension(:,:), allocatable :: tower_ind_temp     !! indicator function of the tower 
real(rprec), dimension(:,:,:), allocatable :: tower_nodes      !! The computational nodes located in the tower 
real(rprec), dimension(:,:), allocatable :: tower_ind          !! indicaton function of the tower
real(rprec), dimension(:), allocatable :: tower_xloc           !! streamwise location
real(rprec), dimension(:), allocatable :: tower_yloc           !! spanwise location
real(rprec), dimension(:), allocatable :: tower_height         !! height of tower center
real(rprec), dimension(:), allocatable :: tower_dia            !! tower diameter
real(rprec), dimension(:), allocatable :: tower_r              !! tower radius
real(rprec), dimension(:), allocatable :: tower_thick          !! tower tickness (in streamwise direction)
real(rprec), dimension(:), allocatable :: tower_vol            !! volume of tower cylinder
real(rprec), dimension(:), allocatable :: tower_theta          !! angle of tower
real(rprec), dimension(:), allocatable :: tower_Cd             !! Drag coefficient of the tower
real(rprec), dimension(:), allocatable :: tower_nhat1          !! x1 direction
real(rprec), dimension(:), allocatable :: tower_nhat2          !! x2 direction
real(rprec), dimension(:), allocatable :: tower_nhat3          !! x3 direction
real(rprec), dimension(:), allocatable :: tower_u_avg          !! time averaged tower velocity
integer    , dimension(:), allocatable :: tower_num_nodes      !! Indicates the actual counted number of grid points in the tower area
integer    , dimension(:), allocatable :: tower_file_id        !! Used for output files
integer    , dimension(:), allocatable :: tower_min_i,tower_max_i,tower_min_j !! counters
integer    , dimension(:), allocatable :: tower_max_j,tower_min_k,tower_max_k !! counters
!! End define for each tower seperately
integer :: tower_trunc !! fixed parameter
real(rprec) :: tower_alpha !! indicate howmuch Gaussian kernel is smoothed

real(rprec) ::   filter_cutoff = 1e-2                          !! Filtering out small forces
real(rprec) ::   tower_fraction = 0.1/3.                       !! Tower diameter 
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tower_init()
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

allocate(tower_xloc(numberOfTurbines))
allocate(tower_yloc(numberOfTurbines))
allocate(tower_height(numberOfTurbines))
allocate(tower_dia(numberOfTurbines))
allocate(tower_r(numberOfTurbines))
allocate(tower_thick(numberOfTurbines))
allocate(tower_vol(numberOfTurbines))
allocate(tower_theta(numberOfTurbines))
allocate(tower_Cd(numberOfTurbines))
allocate(tower_nhat1(numberOfTurbines))
allocate(tower_nhat2(numberOfTurbines))
allocate(tower_nhat3(numberOfTurbines))
allocate(tower_u_avg(numberOfTurbines))
allocate(tower_num_nodes(numberOfTurbines))
allocate(tower_file_id(numberOfTurbines))
allocate(tower_min_i(numberOfTurbines))
allocate(tower_min_j(numberOfTurbines))
allocate(tower_min_k(numberOfTurbines))
allocate(tower_max_i(numberOfTurbines))
allocate(tower_max_j(numberOfTurbines))
allocate(tower_max_k(numberOfTurbines))

do ss=1,numberOfTurbines
!! Properties of the tower should be set here 
tower_Cd(ss)     =0.6863_rprec ! RS tower update
!tower_xloc(ss)   =wind_farm%turbine(ss)%xloc   ! streamwise location of tower (in non dimensional units) 
!tower_yloc(ss)   =wind_farm%turbine(ss)%yloc   ! spanwise location of tower (in non dimensional units)
!tower_height(ss) =wind_farm%turbine(ss)%height ! vertical height of tower center (in non dimensional units)
!tower_dia(ss)    =wind_farm%turbine(ss)%dia*tower_fraction    ! diameter of tower (in non dimensional units)
!tower_thick  =wind_farm%turbine(1)%thk    ! tickness of tower (in non dimensional units)
!tower_theta(ss)  =wind_farm%turbine(ss)%theta1 ! angle of the tower with respect to streamwise direction
!tower_alpha  =alpha                       ! parameter indicating the smoothing of the force
!tower_trunc  =trunc                       ! parameter indicating over how many grid points the force is smoothed
tower_xloc(ss)   = turbineArray(ss) % baseLocation(1) / z_i   ! streamwise location of tower (in non dimensional units) 
tower_yloc(ss)   = turbineArray(ss) % baseLocation(2) / z_i   ! spanwise location of tower (in non dimensional units)
tower_height(ss) = turbineModel(ss) % TowerHt / z_i ! vertical height of tower center (in non dimensional units)
tower_dia(ss)    = turbineModel(ss) % TipRad * 2. / z_i * tower_fraction    ! diameter of tower (in non dimensional units)
tower_thick(ss)  = 1.01*dx    ! tickness of tower (in non dimensional units)
tower_theta  = 0. ! angle of the tower with respect to streamwise direction
tower_alpha  =1.3                       ! parameter indicating the smoothing of the force
!tower_trunc  =2                           ! parameter indicating over how many !grid points the force is smoothed
tower_trunc  =3                       ! parameter indicating over how many grid points the force is smoothed

!tower_thick  = max(tower_thick,dx) ! set tower tickness when dx is larger than tower
tower_vol(ss)     = pi/4.*(tower_dia(ss))**2*(tower_height(ss)-nacelle_dia(ss)) !RS tower update

if (coord==0) then
write(*,*) '------ tower input ------'
write(*,*) 'tower_Cd', tower_Cd(ss)     
write(*,*) 'tower_xloc', tower_xloc(ss)   
write(*,*) 'tower_yloc', tower_yloc(ss)
write(*,*) 'tower_height', tower_height(ss)
write(*,*) 'tower_dia', tower_dia(ss)
!write(*,*) 'tower_thick', tower_thick  
write(*,*) 'tower_theta', tower_theta(ss)
write(*,*) 'tower_alpha', tower_alpha
write(*,*) 'tower_trunc', tower_trunc
endif
enddo

!! Finds the grid cells that are located in the tower area
call tower_nodes_func

!! Makes an indicate function for the tower based on the poinst in the tower area. The indicator function is
!! used to calculate the forces that need to be applied in the tower area. The magnitude of the indicate function
!! indicates how much of the force is applied in a particular gridpoint. At the moment use a Gaussian smoothing kernel.
call tower_filter_ind()

!if (turbine_cumulative_time) then
!     string1 = path // 'tower/tower_u_d_T.dat'
!     inquire (file=string1,exist=exst)
!     if (exst) then
!        write(*,*) 'Reading from file tower_u_d_T.dat'
!        open(1,file=string1)
!        do ss=1,nloc
!        read(1,*) tower_u_avg(ss)
!        enddo
!        close(1)
!     else
!     do ss=1,nloc
!     tower_u_avg(ss) = -7._rprec
!     enddo
!     endif   
!else
!do ss=1,nloc
!tower_u_avg(ss) = -7._rprec
!enddo
!endif

! Generate the files for the turbine forcing output
!if(coord==0) then
!do ss=1,nloc
!call string_splice( string1, path // 'tower/tower_', ss, '_forcing.dat' )
!tower_file_id(ss) = open_file( string1, 'append', 'formatted' )
!enddo
!endif

end subroutine tower_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tower_nodes_func
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(rprec),pointer,dimension(:)::x,y,z
real(rprec) :: rx,ry,rz,r,r_disk
integer :: i,j,k,i2,j2
integer :: imax,jmax,kmax,count_i,count_n,icp,jcp,kcp,kcor,ss,buffer

!! Set the grid variables
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

do ss=1,numberOfTurbines
tower_num_nodes(ss) = 0  !number of nodes that are influenced by tower
count_n = 0            !used for counting nodes for each tower 
count_i = 1            !index count - used for writing to array "nodes"

!!identify "search area"
tower_r(ss)  = tower_dia(ss)/2.0d0 !! set the radius of the tower
imax       = tower_r(ss)/dx + 2  !! number of streamwise points to be searched for tower radius (+2 is safety) 
jmax       = tower_r(ss)/dy + 2  !! number of spanwise points to be searched for tower radius (+2 is safety)
kmax       = tower_r(ss)/dz + 2  !! number of vertical points to be searched for tower radius (+2 is safety)
kcor       = nacelle_dia(ss)/dz !RS tower update

!!determine unit normal vector with respect to tower
tower_nhat1(ss) = -cos(pi*tower_theta(ss)/180.)
tower_nhat2(ss) = -sin(pi*tower_theta(ss)/180.)
tower_nhat3(ss) =  0.0d0

!!determine nearest (i,j,k) to center tower
icp = nint(tower_xloc(ss)/dx)
jcp = nint(tower_yloc(ss)/dy)
kcp = nint(tower_height(ss)/dz + 0.5)

!!determine limits for checking i,j,k
tower_min_i(ss) = icp-imax
tower_max_i(ss) = icp+imax
tower_min_j(ss) = jcp-jmax
tower_max_j(ss) = jcp+jmax
tower_min_k(ss) = max((kcp-kmax),1)
tower_max_k(ss) = min((kcp+kmax),nz_tot)
tower_max_k(ss) = kcp-kcor !RS tower update

!!check grid points around the tower center.	
!do k=tower_min_k,tower_max_k
do k=1,tower_max_k(ss) !RS tower update
   do j=tower_min_j(ss),tower_max_j(ss)
      do i=tower_min_i(ss),tower_max_i(ss)
      !!vector from center point to this node is (rx,ry,rz) with length r
      !!if statements are used to handle periodic boundary conditions
      if (i<1) then
      i2 = mod(i+nx-1,nx)+1
      rx = (x(i2)-L_x) - tower_xloc(ss)
      elseif (i>nx) then
      i2 = mod(i+nx-1,nx)+1
      rx = (L_x+x(i2)) - tower_xloc(ss)
      else
      i2 = i
      rx = x(i) - tower_xloc(ss)
      endif            
      if (j<1) then
      j2 = mod(j+ny-1,ny)+1
      ry = (y(j2)-L_y) - tower_yloc(ss)
      elseif (j>ny) then
      j2 = mod(j+ny-1,ny)+1
      ry = (L_y+y(j2)) - tower_yloc(ss)
      else
      j2 = j
      ry = y(j) - tower_yloc(ss)
      endif
      !rz = z_tot(k) - tower_height 
      rz=0._rprec
      r = sqrt(rx*rx + ry*ry + rz*rz)
      !length projected onto unit normal for tower
      !r_norm = abs(rx*tower_nhat1 + ry*tower_nhat2 + rz*tower_nhat3)
      !(remaining) length projected onto tower
      !r_disk = sqrt(r*r - r_norm*r_norm)
      r_disk = sqrt(r*r)
      !if r_disk<R_t and r_norm within thk/2 from tower than in tower
     ! if ( (r_disk .LE. tower_r) .AND. (r_norm .LE. tower_thick/2.) ) then
!      if ( (r_disk .LE. tower_r) ) then
      if ( (abs(ry) .LE. 0.5*dy) .AND. (abs(rx) .LE. 0.5*dx )) then ! RS makes sure there is always a tower
!          write(*,*) ry,tower_r(ss),rx,0.5*dx,i,j,k
!         tower_ind(count_i)     = 1._rprec 
!         tower_nodes(count_i,1) = i2
!         tower_nodes(count_i,2) = j2
!         tower_nodes(count_i,3) = k   !global k (might be out of this proc's range)
         count_n = count_n + 1
         count_i = count_i + 1
      endif
      enddo ! end i
   enddo ! end j 
enddo ! end k

tower_num_nodes(ss) = count_n !! set num_nodes equal to the number of grid points found in the tower area

! If these conditions are true the code above needs to be changed. Now assumes
! we should place the tower in the closest grid points
!write(*,*) maxval(tower_r), 0.5*dx
if (tower_r(ss)>0.5*dx) write(*,*) 'tower warning dx'
if (tower_r(ss)>0.5*dy) write(*,*) 'tower warning dy'

enddo

!! allocate here
buffer=maxval(tower_num_nodes)
allocate(tower_ind_temp(buffer,numberOfTurbines))
allocate(tower_nodes_temp(buffer,3,numberOfTurbines))

do ss=1,numberOfTurbines
count_n = 0            !used for counting nodes for each tower 
count_i = 1            !index count - used for writing to array "nodes"

!!check grid points around the tower center.  
!do k=tower_min_k,tower_max_k
do k=1,tower_max_k(ss) !RS tower update
   do j=tower_min_j(ss),tower_max_j(ss)
      do i=tower_min_i(ss),tower_max_i(ss)
      !!vector from center point to this node is (rx,ry,rz) with length r
      !!if statements are used to handle periodic boundary conditions
      if (i<1) then
      i2 = mod(i+nx-1,nx)+1
      rx = (x(i2)-L_x) - tower_xloc(ss)
      elseif (i>nx) then
      i2 = mod(i+nx-1,nx)+1
      rx = (L_x+x(i2)) - tower_xloc(ss)
      else
      i2 = i
      rx = x(i) - tower_xloc(ss)
      endif
      if (j<1) then
      j2 = mod(j+ny-1,ny)+1
      ry = (y(j2)-L_y) - tower_yloc(ss)
      elseif (j>ny) then
      j2 = mod(j+ny-1,ny)+1
      ry = (L_y+y(j2)) - tower_yloc(ss)
      else
      j2 = j
      ry = y(j) - tower_yloc(ss)
      endif
      !rz = z_tot(k) - tower_height 
      rz=0._rprec
      r = sqrt(rx*rx + ry*ry + rz*rz)
      !length projected onto unit normal for tower
      !r_norm = abs(rx*tower_nhat1 + ry*tower_nhat2 + rz*tower_nhat3)
      !(remaining) length projected onto tower
      !r_disk = sqrt(r*r - r_norm*r_norm)
      r_disk = sqrt(r*r)
      !if r_disk<R_t and r_norm within thk/2 from tower than in tower
     ! if ( (r_disk .LE. tower_r) .AND. (r_norm .LE. tower_thick/2.) ) then
      !if ( (r_disk .LE. tower_r) ) then
      if ( (abs(ry) .LE. 0.5*dy) .AND. (abs(rx) .LE. 0.5*dx )) then ! RS makes sure there is always a tower
         tower_ind_temp(count_i,ss)     = 1._rprec
         tower_nodes_temp(count_i,1,ss) = i2
         tower_nodes_temp(count_i,2,ss) = j2
         tower_nodes_temp(count_i,3,ss) = k   !global k (might be out of this proc's range)
         count_n = count_n + 1
         count_i = count_i + 1
      endif
      enddo ! end i
   enddo ! end j 
enddo ! end k

if (coord==0) then
write(*,*) 'tower_num_nodes',coord,tower_num_nodes(ss)
endif
!! Check the conditions in the two loops above when this message is shown
if (count_n.ne.tower_num_nodes(ss)) write(*,*) 'something is going wrong above 1' 

enddo

nullify(x,y,z)

end subroutine tower_nodes_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tower_filter_ind()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(rprec), dimension(:,:,:,:), allocatable :: out_a
real(rprec), dimension(:,:,:), allocatable :: g,temp_array
real(rprec) :: sumA,sumG,delta2,r2,dxdydz
integer :: i,j,k,i2,j2,k2,i3,j3,i4,j4,l,ssx,ssy,ssz
integer :: count_i,count_n
integer :: k_start,k_end,ss,buffer,nz_tower

dxdydz=dx*dy*dz
nz_tower=nz_tot*(maxval(tower_height)+maxval(tower_dia))+2*tower_trunc+1
if (coord==0) write(*,*) 'nz_tower',nz_tower

!! Allocation of temporary arrays. 
!! It should be possible to optimize memory usage here. It is only used at the startup of the simulation
allocate(out_a(nx,ny,nz_tower,numberOfTurbines))
allocate(temp_array(nx,ny,nz_tower))

!! Gaussian smoothing kernel, centered around middle grid point. Matrix made smaller to limit memory use. Can probably be limited more
allocate(g(nx/4:nx/4*3,ny/4:ny/4*3,nz_tot/4:nz_tot/4*3))
delta2 = tower_alpha**2 * (dx**2 + dy**2 + dz**2)
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
do l=1,tower_num_nodes(ss)
   i2                   = tower_nodes_temp(l,1,ss)
   j2                   = tower_nodes_temp(l,2,ss)
   k2                   = tower_nodes_temp(l,3,ss)
   temp_array(i2,j2,k2) = tower_ind_temp(l,ss)
enddo

!!convolution computed for points (i4,j4,k)
!!only compute for nodes near the tower (possible influence region indicated by trunc)
!do k=max(tower_min_k-tower_trunc,1),min(tower_max_k+tower_trunc,nz_tot)
do k=1,min(tower_max_k(ss)+tower_trunc,nz_tot) !RS tower update
do j=(tower_min_j(ss)-tower_trunc),(tower_max_j(ss)+tower_trunc)
do i=(tower_min_i(ss)-tower_trunc),(tower_max_i(ss)+tower_trunc)
        
   i4 = mod(i+nx-1,nx)+1       !since values may be out 1-nx,1-ny domain (spectral BCs)
   j4 = mod(j+ny-1,ny)+1              
        
   !!for each (i4,j4,k), center convolution function on that point and 'integrate' 
   !!relative coords are (ssx,ssy,ssz). absolute coords of other/surrounding points are (i2,j2,k2)
   !!only need to consider other/surrounding points near (i4,j4,k) since conv. function is compact

   do k2=max(k-tower_trunc,1),min(k+tower_trunc,nz_tot)     !currently using truncated Gaussian
   do j2=j-tower_trunc,j+tower_trunc
   do i2=i-tower_trunc,i+tower_trunc

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
do k=1,nz_tower  !global k
do j=1,ny
do i=1,nx
   if (out_a(i,j,k,ss) < filter_cutoff) then
      out_a(i,j,k,ss)=0._rprec
   endif
enddo
enddo
enddo

!!normalize this "indicator function" 
sumA = sum(out_a(1:nx,1:ny,1:nz_tower,ss))*dxdydz
out_a(1:nx,1:ny,1:nz_tower,ss) = tower_vol(ss)/sumA*out_a(1:nx,1:ny,1:nz_tower,ss)

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
do k=k_start,min(k_end,nz_tower)  !global k
do j=1,ny
do i=1,nx
   !! store all points in which the indicator function has some value. This is the area of the tower where the forces should be applied
   if (out_a(i,j,k,ss) > filter_cutoff/4._rprec) then ! !RS tower update
!      tower_ind  (count_i)   = out_a(i,j,k,ss)
!      tower_nodes(count_i,1) = i
!      tower_nodes(count_i,2) = j
!      tower_nodes(count_i,3) = k - coord*(nz-1)  !local k
!      if (coord==0) write(*,*) 'tower',i,j,k
      count_n = count_n + 1
      count_i = count_i + 1
   endif !RS tower update
enddo
enddo
enddo

tower_num_nodes(ss)=count_n
enddo

!! allocate here
buffer=maxval(tower_num_nodes)
allocate(tower_ind(buffer,numberOfTurbines))
allocate(tower_nodes(buffer,3,numberOfTurbines))

do ss=1,numberOfTurbines
count_n=0 !! reset counter
count_i=1 !! reset counter
do k=k_start,min(k_end,nz_tower)  !global k
do j=1,ny
do i=1,nx
   !! store all points in which the indicator function has some value. This is the area of the tower where the forces should be applied
   if (out_a(i,j,k,ss) > filter_cutoff/4._rprec) then !RS update
      tower_ind  (count_i,ss)   = out_a(i,j,k,ss)
      tower_nodes(count_i,1,ss) = i
      tower_nodes(count_i,2,ss) = j
      tower_nodes(count_i,3,ss) = k - coord*(nz-1)  !local k
      count_n = count_n + 1
      count_i = count_i + 1
   endif
enddo
enddo
enddo

if (tower_num_nodes(ss)>0) then
write(*,*) 'tower_num_nodes',coord,tower_num_nodes(ss)
endif
!! Check the conditions in the two loops above when this message is shown
if (count_n.ne.tower_num_nodes(ss)) write(*,*) 'something is going wrong above 2' 
enddo 
 
deallocate(out_a,g,temp_array)
deallocate(tower_ind_temp,tower_nodes_temp)

end subroutine tower_filter_ind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tower_forcing()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use sim_param, only: u,fxa
implicit none
real(rprec)::tower_avg_vels,tower_avg_vels2,tower_force,u_tower,eps,T_avg_dim
integer :: i2,j2,k2,l,ss

T_avg_dim=0.02

do ss=1,numberOfTurbines
tower_avg_vels = 0.0d0 !! average velocity in the tower area
tower_force    = 0.0d0 !! force that should be applied in the tower area

!! Each core calculates the wheighted average of the velocity of grid points in the tower area
!! For cores without nodes num_nodes=0, and this loop is not computed 
do l=1,tower_num_nodes(ss)
   i2 = tower_nodes(l,1,ss)
   j2 = tower_nodes(l,2,ss)
   k2 = tower_nodes(l,3,ss)
   tower_avg_vels = tower_avg_vels-u(i2,j2,k2)*tower_ind(l,ss)        
enddo

! Determine the sum over all processors to get the total weighted sum
$if ($MPI) 
call mpi_allreduce(tower_avg_vels,tower_avg_vels2, 1, MPI_RPREC, MPI_SUM, comm, ierr)
$endif

!update epsilon for the new timestep (for cfl_dt) 
eps = (dt_dim / T_avg_dim) / (1._rprec + dt_dim / T_avg_dim)
!volume correction
u_tower = tower_avg_vels2 * dx*dy*dz/tower_vol(ss)
!add this current value to the "running average" (first order relaxation)
tower_u_avg(ss) = (1.-eps)*tower_u_avg(ss) + eps*u_tower

!calculate total thrust force for each tower  (per unit mass)
tower_force = 0.5_rprec*tower_Cd(ss)*abs(tower_u_avg(ss))*tower_u_avg(ss)/(pi*tower_dia(ss)/4._rprec) !RS tower update

!if (coord==0) then
!if (modulo (jt_total, tbase) == 0) then
!call write_real_data( tower_file_id(ss), 'formatted', 4, (/total_time_dim, u_tower, tower_u_avg(ss), tower_force/))
!endif
!endif 

!! Each core applies the forces to nodes in its domain based on the force calculated above and the fraction
!! of the force that should be applied per point (indicator function) 
!! or cores without nodes num_nodes=0, and this loop is not computed 
do l=1,tower_num_nodes(ss)
   i2 = tower_nodes(l,1,ss)
   j2 = tower_nodes(l,2,ss)
   k2 = tower_nodes(l,3,ss)
   fxa(i2,j2,k2) = fxa(i2,j2,k2)+tower_force*tower_ind(l,ss) !RS tower update
enddo

enddo

end subroutine tower_forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine tower_finalize ()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!implicit none
!character (100) :: string1
!integer :: ss
!!write disk-averaged velocity to file along with T_avg_dim
!if (coord == 0) then
!    string1 = path // 'tower/tower_u_d_T.dat'
!    open (unit=1,file = string1, status='unknown',form='formatted', action='write',position='rewind')
!    do ss=1,nloc
!    write(1,*) tower_u_avg(ss)
!    enddo
!    close (1)
!endif
!
!end subroutine tower_finalize

end module tower
