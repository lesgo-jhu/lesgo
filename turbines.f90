module turbines
use types,only:rprec
use param
use stat_defs, only:wind_farm_t
use grid_defs, only:x,y,z
use io

implicit none

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize

integer :: nloc 
integer :: num_x,num_y
real(rprec) :: height_all,dia_all,thk_all,theta1_all,theta2_all
real(rprec) :: Ct_prime,Ct_noprime   					    !thrust coefficient
real(rprec) :: T_avg_dim, T_avg_dim_file
real(rprec), dimension(nz_tot) :: z_tot
real(rprec) :: sx,sy

character (64) :: fname0, fname, fname3, fname4, var_list, temp
real(rprec), dimension(nx,ny,nz_tot) :: large_node_array    !used for visualizing node locations
real(rprec), dimension(nx,ny,nz_tot) :: large_node_array_filtered

real(rprec) :: eps											!epsilon used for disk velocity time-averaging

integer :: i,j,k,i2,j2,k2,b,l,s,nn,ssx,ssy,ssz
integer :: imax,jmax,kmax,count_i,count_n,icp,jcp,kcp
integer :: min_i,max_i,min_j,max_j,min_k,max_k,cut
integer :: k_start, k_end
real(rprec) :: a0,a1,a2,a3,a4
character (4) :: string1, string2, string3
logical :: exst

logical :: turbine_in_proc=.false.      !false if there is no turbine (partial or whole) in this processor

real(rprec), pointer, dimension(:) :: buffer_array
real :: buffer
logical :: buffer_logical
integer, dimension(nproc-1) :: turbine_in_proc_array = 0
integer :: turbine_in_proc_cnt = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_init()
implicit none

!##############################  SET BY USER  ############################################
!set turbine parameters
!turbines are numbered as follows:
!   #1 = turbine nearest (x,y)=(0,0)
!   #2 = next turbine in the x-direction, etc.

    num_x = 4               !number of turbines in x-direction
    num_y = 6               !number of turbines in y-direction  
	nloc = num_x*num_y      !number of turbines (locations) 
	allocate(wind_farm_t%turbine_t(nloc))
    allocate(buffer_array(nloc))

    !x,y-locations
        !for evenly-spaced turbines, not staggered
        if(.true.) then
            
            k=1
            do j=1,num_y
                do i=1,num_x
                    wind_farm_t%turbine_t(k)%xloc = L_x*real(2*i-1)/real(2*num_x)
                    wind_farm_t%turbine_t(k)%yloc = L_y*real(2*j-1)/real(2*num_y)
                    k = k + 1
                enddo
            enddo
        !for other orientation
        else
            !** set individual locations here
        endif

    !height, diameter, and thickness
        !same values for all
        if(.true.) then
            height_all = 100.       !turbine height, dimensional [m]
            dia_all = 100.	        !turbine diameter, dimensional [m]
            thk_all = 10.	        !turbine disk thickness, dimensional [m]
            
            !non-dimensionalize values by z_i
            height_all = height_all/z_i
            dia_all = dia_all/z_i
            thk_all = thk_all/z_i
            thk_all = max(thk_all,dx*1.01)	
            
            do k=1,nloc
                wind_farm_t%turbine_t(k)%height = height_all
                wind_farm_t%turbine_t(k)%dia = dia_all
                wind_farm_t%turbine_t(k)%thk = thk_all
            enddo
        !for varying height and/or diameter and/or thickness
        else
            !** set individual values here
        endif          

    !orientation (angles)
        !same values for all
        if(.true.) then
            theta1_all = 0.     !angle CCW(from above) from -x direction [degrees]
            theta2_all = 0.     !angle above the horizontal, from -x dir [degrees]
            
            do k=1,nloc
                wind_farm_t%turbine_t(k)%theta1 = theta1_all
                wind_farm_t%turbine_t(k)%theta2 = theta2_all
            enddo
        !for varying angles
        else
            !** set individual values here
        endif     

    !filtering operation
        wind_farm_t%ifilter = 2			    !Filter type: 0-> none 1->cut off 2->Gaussian 3->Top-hat	
        wind_farm_t%alpha = 1.5			    !filter size is alpha*(grid spacing)
        wind_farm_t%trunc = 3               !truncated Gaussian - how many grid points in any direction
        wind_farm_t%filter_cutoff = 1e-2    !ind only includes values above this cutoff

    !other
	    Ct_prime = 1.33		!thrust coefficient
        Ct_noprime = 0.75   !a=1/4
        T_avg_dim = 600.     !time-averaging 'window' for one-sided exp. weighting (seconds)
        
        sx = L_x/(num_x*dia_all)        !spacing in x-dir, multiple of DIA            
        sy = L_y/(num_y*dia_all)        !spacing in y-dir, multiple of DIA
!#########################################################################################

!z_tot for total domain (since z is local to the processor)
    do k=1,nz_tot
        z_tot(k) = (k - 0.5_rprec) * dz
    enddo

!find turbine nodes - including unfiltered ind, n_hat, num_nodes, and nodes for each turbine
!each processor finds turbines in the entire domain
    large_node_array = 0.
	call turbines_nodes(large_node_array)

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
	!to write the node locations to file
	  fname0 = 'nodes_unfiltered.dat'
	  call write_tecplot_header_ND(fname0,'rewind', 4, (/nx+1, ny+1, nz_tot/), '"x", "y", "z", "nodes_unfiltered"', 0, 1)
	  call write_real_data_3D(fname0, 'append','formatted', 1, nx, ny, nz_tot, (/large_node_array/), 4, x,y,z_tot)
    endif

!1.smooth/filter indicator function                     
!2.associate new nodes with turbines                               
!3.normalize such that each turbine's ind integrates to turbine volume
!4.split domain between processors 
    call turbines_filter_ind()
    
!set variables for time-averaging velocity 
    eps = dt_dim/T_avg_dim / (1. + dt_dim/T_avg_dim)
    
    if (cumulative_time) then
        if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
            fname4 = 'turbine_u_d_T.dat'
            inquire (file=fname4, exist=exst)
            if (exst) then
                write(*,*) 'Reading from file turbine_u_d_T.dat'
                open (1, file=fname4)
                do i=1,nloc
                    read(1,*) wind_farm_t%turbine_t(i)%u_d_T    
                enddo    
                read(1,*) T_avg_dim_file
                if (T_avg_dim_file /= T_avg_dim) then
                    write(*,*) 'Time-averaging window does not match value in turbine_u_d_T.dat'
                endif
                close (1)
            else  
                write (*, *) 'File ', fname4, ' not found'
                write (*, *) 'Assuming u_d_T = -7. for all turbines'
                do k=1,nloc
                    wind_farm_t%turbine_t(k)%u_d_T = -7.
                enddo
            end if                 
        endif
    else
        write (*, *) 'Assuming u_d_T = -7 for all turbines'
        do k=1,nloc
            wind_farm_t%turbine_t(k)%u_d_T = -7.
        enddo    
    endif
    
!create files to store turbine forcing data
    if (cumulative_time) then
        if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
            var_list = '"t (s)", "u_d", "u_d_T", "f_n", "P"'
            call write_tecplot_header_xyline('turbine_x1y1_forcing.dat','rewind', var_list)     ! 1 --turbine number
            call write_tecplot_header_xyline('turbine_x2y1_forcing.dat','rewind', var_list)     ! 2
            call write_tecplot_header_xyline('turbine_x3y1_forcing.dat','rewind', var_list)     ! 3
            call write_tecplot_header_xyline('turbine_x4y1_forcing.dat','rewind', var_list)     ! 4
            
            call write_tecplot_header_xyline('turbine_x1y2_forcing.dat','rewind', var_list)     ! 5
            call write_tecplot_header_xyline('turbine_x1y3_forcing.dat','rewind', var_list)     ! 9
            call write_tecplot_header_xyline('turbine_x1y4_forcing.dat','rewind', var_list)     ! 13
            call write_tecplot_header_xyline('turbine_x1y5_forcing.dat','rewind', var_list)     ! 17
            call write_tecplot_header_xyline('turbine_x1y6_forcing.dat','rewind', var_list)     ! 21
        endif
    endif
	
end subroutine turbines_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_nodes(array)
!This subroutine locates nodes for each turbine and builds the arrays: ind, n_hat, num_nodes, and nodes
implicit none

real(rprec) :: R_t,rx,ry,rz,r,r_norm,r_disk
real(rprec), dimension(nx,ny,nz_tot) :: array

real(rprec), pointer :: p_dia => null(), p_thk=> null(), p_theta1=> null(), p_theta2=> null()
real(rprec), pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3=> null() 
real(rprec), pointer :: p_xloc => null(), p_yloc=> null(), p_height=> null()

logical :: verbose = .false.

do s=1,nloc
    
    count_n = 0		!used for counting nodes for each turbine
    count_i = 1		!index count - used for writing to array "nodes"

    !set pointers
        p_xloc => wind_farm_t%turbine_t(s)%xloc     
        p_yloc => wind_farm_t%turbine_t(s)%yloc  
        p_height => wind_farm_t%turbine_t(s)%height 
        p_dia => wind_farm_t%turbine_t(s)%dia 
        p_thk => wind_farm_t%turbine_t(s)%thk
        p_theta1 => wind_farm_t%turbine_t(s)%theta1
        p_theta2 => wind_farm_t%turbine_t(s)%theta2
        p_nhat1 => wind_farm_t%turbine_t(s)%nhat(1)
        p_nhat2 => wind_farm_t%turbine_t(s)%nhat(2)
        p_nhat3 => wind_farm_t%turbine_t(s)%nhat(3)

    !identify "search area"
    	R_t = p_dia/2.
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                if (verbose) then
                  write(*,*) '     rad:',R_t
                  write(*,*) '     dx:',dx
                  write(*,*) '     dy:',dy
                  write(*,*) '     dz:',dz
                  write(*,*) '     thk:',p_thk
                endif
            endif
    	imax = R_t/dx + 2
    	jmax = R_t/dy + 2
    	kmax = R_t/dz + 2

    !determine unit normal vector for each turbine	
    	p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    	p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    	p_nhat3 = sin(pi*p_theta2/180.)
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                if (verbose) then
                  write(*,*) '     n_hat', p_nhat1, p_nhat2, p_nhat3
                endif
            endif

    !determine nearest (i,j,k) to turbine center
    	icp = nint(p_xloc/dx)
    	jcp = nint(p_yloc/dy)
    	kcp = nint(p_height/dz + 0.5)
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                if (verbose) then
                  write(*,*) '     turbine center (i,j,k)',icp,jcp,kcp
                  write(*,*) '     turbine center (x,y,z)',x(icp),y(jcp),z_tot(kcp)
                endif
            endif

    !determine limits for checking i,j,k
    	min_i = max((icp-imax),1)
    	max_i = min((icp+imax),nx)
    	min_j = max((jcp-jmax),1)
    	max_j = min((jcp+jmax),ny)
    	min_k = max((kcp-kmax),1)
    	max_k = min((kcp+kmax),nz_tot)
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                if (verbose) then
                  write(*,*) '     i limits, range', min_i, max_i, dx*(max_i-min_i)
                  write(*,*) '     j limits, range', min_j, max_j, dy*(max_j-min_j)
                  write(*,*) '     k limits, range', min_k, max_k, dz*(max_k-min_k)
                endif
            endif
            wind_farm_t%turbine_t(s)%nodes_max(1) = min_i
            wind_farm_t%turbine_t(s)%nodes_max(2) = max_i
            wind_farm_t%turbine_t(s)%nodes_max(3) = min_j
            wind_farm_t%turbine_t(s)%nodes_max(4) = max_j
            wind_farm_t%turbine_t(s)%nodes_max(5) = min_k
            wind_farm_t%turbine_t(s)%nodes_max(6) = max_k            

    !check neighboring grid points	
	do k=min_k,max_k
		do j=min_j,max_j
			do i=min_i,max_i
			!vector from center point to this node is (rx,ry,rz) with length r
				rx = x(i) - p_xloc
				ry = y(j) - p_yloc
				rz = z_tot(k) - p_height 
				r = sqrt(rx*rx + ry*ry + rz*rz)
			!length projected onto unit normal for this turbine
				r_norm = abs(rx*p_nhat1 + ry*p_nhat2 + rz*p_nhat3)
			!(remaining) length projected onto turbine disk
				r_disk = sqrt(r*r - r_norm*r_norm)
			!if r_disk<R_t and r_norm within thk/2 from turbine -- this node is part of the turbine
				if ( (r_disk .LE. R_t) .AND. (r_norm .LE. p_thk/2.) ) then
                    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                        if (verbose) then
                          write(*,*) '     FOUND NODE', i,j,k
                        endif
                    endif
					array(i,j,k) = 1.
                    wind_farm_t%turbine_t(s)%ind(count_i) = 1. 
					wind_farm_t%turbine_t(s)%nodes(count_i,1) = i
					wind_farm_t%turbine_t(s)%nodes(count_i,2) = j
					wind_farm_t%turbine_t(s)%nodes(count_i,3) = k   !global k (might be out of this proc's range)
                    count_n = count_n + 1
					count_i = count_i + 1
				endif
			enddo
		enddo
	enddo
	wind_farm_t%turbine_t(s)%num_nodes = count_n

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        write(*,*) '     Turbine #',s,'has',count_n,'unfiltered nodes in entire domain'
    endif
    
enddo
		
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

real(rprec), dimension(nx,ny,nz_tot) :: out_a, g, g_shift, fg
real(rprec), dimension(nx,ny,nz_tot) :: temp_array
real(rprec), dimension(nx,ny,nz) :: temp_array_2
real(rprec) :: sumG,delta2,r2,sumA
real(rprec) :: turbine_vol

logical :: verbose = .false.

!create convolution function, centered at (nx/2,ny/2,(nz_tot-1)/2) and normalized
	if(wind_farm_t%ifilter==0) then		!0-> none
	  !
	elseif(wind_farm_t%ifilter==1) then  	!1-> cutoff/sharp spectral
	  !
	elseif(wind_farm_t%ifilter==2) then		!2-> Gaussian
	  delta2 = wind_farm_t%alpha**2 * (dx**2 + dy**2 + dz**2)
      do k=1,nz_tot
    	  do j=1,ny
    	    do i=1,nx
    	      r2 = ((real(i)-nx/2.)*dx)**2 + ((real(j)-ny/2.)*dy)**2 + ((real(k)-(nz_tot-1)/2.)*dz)**2
    	      g(i,j,k) = sqrt(6./(pi*delta2))*6./(pi*delta2)*exp(-6.*r2/delta2)
    	    enddo
    	  enddo
      enddo
	elseif(wind_farm_t%ifilter==3) then		!3-> box/tophat
	  !
	elseif(wind_farm_t%ifilter==4) then
	  !g = ??    insert any function here
	endif
	
!normalize the convolution function
	sumG = 0.
    do k=1,nz_tot
	  do j=1,ny-1		!since j=1 and j=ny are the same point
	    do i=1,nx-1		!since i=1 and i=nx are the same point
	      sumG = sumG + g(i,j,k)*dx*dy*dz
	    enddo
	  enddo
    enddo
	g = g/sumG

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        if (verbose) then
            write(*,*) 'Convolution function created.'
        endif
    endif

!display the convolution function
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        if(.false.) then
          write(*,*) 'Convolution function'
          write(*,*) g
          write(*,*) 'integral of g(i,j,k): ',sumG
        endif       
    endif

	!to write the data to file, centered at (i,j,k=(nz_tot-1)/2)
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then    
        i=nx/2
        j=ny/2
        do k2=1,nz_tot
          do j2=1,ny
            do i2=1,nx
            g_shift(i2,j2,k2) = g( mod(i2-i+nx/2+(nx-1)-1,(nx-1))+1 , mod(j2-j+ny/2+(ny-1)-1,(ny-1))+1, k2)
            enddo
          enddo
        enddo

        if (.false.) then
            fname0 = 'convolution_function.dat'
            call write_tecplot_header_ND(fname0,'rewind', 4, (/nx+1,ny+1,nz_tot/), '"x","y","z","g"', 1, 1)
            call write_real_data_3D(fname0, 'append', 'formatted', 1, nx, ny, nz_tot, (/g_shift/), 4, x, y, z_tot)

            if (verbose) then
                write(*,*) 'Convolution function written to Tecplot file.'
            endif
        endif
    endif

!filter indicator function for each turbine
do b=1,nloc
    
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        if (verbose) then
            write(*,*) 'Filtering turbine Number ',b
        endif
    endif

    !create the input array (nx,ny,nz_tot) from a list of included nodes
        temp_array = 0.
        do l=1,wind_farm_t%turbine_t(b)%num_nodes
            i2 = wind_farm_t%turbine_t(b)%nodes(l,1)
            j2 = wind_farm_t%turbine_t(b)%nodes(l,2)
            k2 = wind_farm_t%turbine_t(b)%nodes(l,3)	
            temp_array(i2,j2,k2) = wind_farm_t%turbine_t(b)%ind(l)
        enddo

    !perform convolution on temp_array --> out_a    
        out_a=0.

        min_i = wind_farm_t%turbine_t(b)%nodes_max(1) 
        max_i = wind_farm_t%turbine_t(b)%nodes_max(2) 
        min_j = wind_farm_t%turbine_t(b)%nodes_max(3) 
        max_j = wind_farm_t%turbine_t(b)%nodes_max(4) 
        min_k = wind_farm_t%turbine_t(b)%nodes_max(5)
        max_k = wind_farm_t%turbine_t(b)%nodes_max(6) 
        cut = wind_farm_t%trunc   

        if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
            if (verbose) then
                write(*,*) 'search over: ',min_i-cut,max_i+cut,min_j-cut,max_j+cut,min_k-cut,max_k+cut
            endif
        endif

        do k=max(min_k-cut,1),min(max_k+cut,nz_tot)    !only compute for nodes near the turbine
    	do j=max(min_j-cut,1),min(max_j+cut,ny)
    	do i=max(min_i-cut,1),min(max_i+cut,nx)
		
          do k2=k-wind_farm_t%trunc,k+wind_farm_t%trunc     !currently using truncated Gaussian
    	  do j2=j-wind_farm_t%trunc,j+wind_farm_t%trunc
    	  do i2=i-wind_farm_t%trunc,i+wind_farm_t%trunc

            if ( (i2>0).and.(j2>0).and.(k2>0).and.(i2<=nx).and.(j2<=ny).and.(k2<=nz_tot) ) then
                ssx = mod(i2-i+nx/2+(nx-1)-1,(nx-1))+1
    	        ssy = mod(j2-j+ny/2+(ny-1)-1,(ny-1))+1       
                ssz = k2-k+(nz_tot-1)/2       !since no spectral BCs in z-direction
            
                if (ssx < 1) then
                    fg(i2,j2,k2) = 0.
                elseif (ssx > nx) then
                    fg(i2,j2,k2) = 0.
                elseif( ssy < 1 ) then
                    fg(i2,j2,k2) = 0.
                elseif( ssy > ny) then
                    fg(i2,j2,k2) = 0.
                elseif( ssz < 1) then
                    fg(i2,j2,k2) = 0.
                elseif( ssz > nz_tot ) then
                    fg(i2,j2,k2) = 0.
                else
                    fg(i2,j2,k2) = temp_array(i2,j2,k2)*g(ssx,ssy,ssz)
                    out_a(i,j,k) = out_a(i,j,k) + fg(i2,j2,k2)*dx*dy*dz
                endif    

            endif	
    	    
    	  enddo
    	  enddo
          enddo
    	enddo
    	enddo
        enddo

        if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
            if (verbose) then
                write(*,*) 'Convolution complete for turbine ',b
            endif
        endif

    !normalize this "indicator function" such that it integrates to turbine volume
	sumA = 0.
    do k=1,nz_tot
	  do j=1,ny-1		!since j=1 and j=ny are the same point
	    do i=1,nx-1		!since i=1 and i=nx are the same point
            if (out_a(i,j,k) < wind_farm_t%filter_cutoff) then
	            out_a(i,j,k) = 0.     !don't want to include too many nodes (truncated Gaussian?)
            else
                sumA = sumA + out_a(i,j,k)*dx*dy*dz
            endif            
	    enddo
	  enddo
    enddo
    turbine_vol = pi/4. * (wind_farm_t%turbine_t(b)%dia)**2 * wind_farm_t%turbine_t(b)%thk
	out_a = turbine_vol/sumA*out_a

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        large_node_array_filtered = large_node_array_filtered + out_a
    endif

    !update num_nodes, nodes, and ind for this turbine
    !and split domain between processors
    !z(nz) and z(1) of neighboring coords match so each coord gets (local) 1 to nz-1
    wind_farm_t%turbine_t(b)%ind = 0.
    wind_farm_t%turbine_t(b)%nodes = 0.
    wind_farm_t%turbine_t(b)%num_nodes = 0.
    count_n = 0
    count_i = 1
    if (.not. USE_MPI) then
        k_start = 1
        k_end = nz
    else !this is the global k. searching over local k=1,nz-1
        k_start = 1+coord*(nz-1)
        k_end = nz-1+coord*(nz-1)
    endif
	do k=k_start,k_end  !global k     
		do j=1,ny
			do i=1,nx
				if (out_a(i,j,k) > wind_farm_t%filter_cutoff) then
                    wind_farm_t%turbine_t(b)%ind(count_i) = out_a(i,j,k)		
					wind_farm_t%turbine_t(b)%nodes(count_i,1) = i
					wind_farm_t%turbine_t(b)%nodes(count_i,2) = j
					wind_farm_t%turbine_t(b)%nodes(count_i,3) = k - coord*(nz-1)   !local k
                    count_n = count_n + 1
					count_i = count_i + 1
                    turbine_in_proc = .true.                    
				endif
			enddo
		enddo
	enddo
	wind_farm_t%turbine_t(b)%num_nodes = count_n
    
    if (count_n > 0) then
        if (.not. USE_MPI) then
            write (string1, '(i3)') b
            string1 = trim(adjustl(string1))
            write (string2, '(i4)') count_n
            string2 = trim(adjustl(string2))            
          
            write(*,*) 'Turbine number ',string1,' has ',string2,' filtered nodes' 
        else
            write (string1, '(i3)') b
            string1 = trim(adjustl(string1))
            write (string2, '(i4)') count_n
            string2 = trim(adjustl(string2)) 
            write (string3, '(i3)') coord
            string3 = trim(adjustl(string3))             

            write(*,*) 'Turbine number ',string1,' has ',string2,' filtered nodes in coord ', string3
        endif
    endif

enddo

    !test to make sure domain is divided correctly:
    if (.false.) then
        temp_array_2 = 0.
        do b=1,nloc
        do l=1,wind_farm_t%turbine_t(b)%num_nodes
            i2 = wind_farm_t%turbine_t(b)%nodes(l,1)
            j2 = wind_farm_t%turbine_t(b)%nodes(l,2)
            k2 = wind_farm_t%turbine_t(b)%nodes(l,3)	
            temp_array_2(i2,j2,k2) = wind_farm_t%turbine_t(b)%ind(l)
        enddo   
        enddo
        !write to file with .dat.c* extension
            fname3 = 'nodes_filtered_c.dat'
            write (temp, '(".c",i0)') coord
            fname3 = trim (fname3) // temp
            call write_tecplot_header_ND(fname3,'rewind', 4, (/nx+1,ny+1,nz/), '"x","y","z","nodes_filtered_c"', 1, 1)
            call write_real_data_3D(fname3, 'append', 'formatted', 1, nx, ny, nz, (/temp_array_2/), 4, x, y, z(1:nz))      
    endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    fname3 = 'nodes_filtered.dat'
    call write_tecplot_header_ND(fname3,'rewind', 4, (/nx+1,ny+1,nz_tot/), '"x","y","z","nodes_filtered"', 1, 1)
    call write_real_data_3D(fname3, 'append', 'formatted', 1, nx, ny, nz_tot, (/large_node_array_filtered/), 4, x, y, z_tot)                       
endif

!each processor sends its value of turbine_in_proc
!if false, disk-avg velocity will not be sent (since it will always be 0.)
!############################################## 2
$if ($MPI)
    if (rank == 0) then
        do i=1,nproc-1
            call MPI_recv( buffer_logical, 1, MPI_logical, i, 2, MPI_COMM_WORLD, status, ierr )
            if (buffer_logical==.true.) then
                write (string3, '(i3)') i
                string3 = trim(adjustl(string3))       
                print*,'I am coord 0 and coord ',string3,' has turbine nodes'            
                turbine_in_proc_cnt = turbine_in_proc_cnt + 1
                turbine_in_proc_array(turbine_in_proc_cnt) = i
            endif
        enddo
    else
        call MPI_send( turbine_in_proc, 1, MPI_logical, 0, 2, MPI_COMM_WORLD, ierr )
    endif
$endif
!##############################################

end subroutine turbines_filter_ind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_forcing()
use sim_param, only: u,v,w
use immersedbc, only: fx,fy,fz
$if ($MPI)
    use mpi_defs, only: mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

real(rprec), pointer :: p_u_d => null(), p_nhat1=> null(), p_nhat2=> null(), p_nhat3=> null()
real(rprec), pointer :: p_u_d_T => null(), p_dia => null(), p_thk=> null(), p_f_n => null()
integer, pointer :: p_num_nodes=> null()

real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vels, disk_force

integer :: w_uv_tag_turbines = -1

$if ($MPI)
    call mpi_sync_real_array(w, MPI_SYNC_DOWNUP)     !syncing intermediate w-velocities!
$endif
call interp_to_uv_grid(w, w_uv, w_uv_tag_turbines)

disk_avg_vels = 0.

!Each processor calculates the weighted disk-averaged velocity
if (turbine_in_proc == .true.) then

    !for each turbine:        
        do s=1,nloc      
             
        !set pointers
            p_u_d => wind_farm_t%turbine_t(s)%u_d   
            p_num_nodes => wind_farm_t%turbine_t(s)%num_nodes
            p_nhat1 => wind_farm_t%turbine_t(s)%nhat(1)
            p_nhat2 => wind_farm_t%turbine_t(s)%nhat(2)
            p_nhat3 => wind_farm_t%turbine_t(s)%nhat(3)

        !calculate total disk-averaged velocity for each turbine (current,instantaneous)    
        !u_d is the velocity in the normal direction	  
        !weighted average using "ind"            
            p_u_d = 0.
            do l=1,p_num_nodes      
                i2 = wind_farm_t%turbine_t(s)%nodes(l,1)
                j2 = wind_farm_t%turbine_t(s)%nodes(l,2)
                k2 = wind_farm_t%turbine_t(s)%nodes(l,3)	
                p_u_d = p_u_d + (p_nhat1*u(i2,j2,k2) + p_nhat2*v(i2,j2,k2) + p_nhat3*w_uv(i2,j2,k2)) &
                    * wind_farm_t%turbine_t(s)%ind(l)        
            enddo              

        !write this value to the array (which will be sent to coord 0)
            disk_avg_vels(s) = p_u_d
            
        enddo
        
endif        

!send the disk-avg values to coord==0
$if ($MPI) 
    !############################################## 3
    if (rank == 0) then
        do i=1,turbine_in_proc_cnt
            j = turbine_in_proc_array(i)
            buffer_array = 0.
            call MPI_recv( buffer_array, nloc, MPI_rprec, j, 3, MPI_COMM_WORLD, status, ierr )
            disk_avg_vels = disk_avg_vels + buffer_array
            
        enddo              
                
    elseif (turbine_in_proc == .true.) then
        call MPI_send( disk_avg_vels, nloc, MPI_rprec, 0, 3, MPI_COMM_WORLD, ierr )
    endif            
    !##############################################              
$endif

!Coord==0 takes that info and calculates total disk force, then sends it back
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then           

    !for each turbine:        
        do s=1,nloc            
             
        !set pointers
            p_u_d => wind_farm_t%turbine_t(s)%u_d   
            p_u_d_T => wind_farm_t%turbine_t(s)%u_d_T   
            p_f_n => wind_farm_t%turbine_t(s)%f_n               
            p_dia => wind_farm_t%turbine_t(s)%dia 
            p_thk => wind_farm_t%turbine_t(s)%thk               
            
        !correction:
        !since sum of ind is turbine volume/(dx*dy*dz) (not exactly 1.)
            p_u_d = disk_avg_vels(s)
            p_u_d = p_u_d *dx*dy*dz/(pi/4.*(wind_farm_t%turbine_t(s)%dia)**2 * wind_farm_t%turbine_t(s)%thk)
    
        !add this current value to the "running average" (first order relaxation)
            p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

        !calculate total thrust force for each turbine  (per unit mass)
        !force is normal to the surface (calc from u_d_T, normal to surface)
            p_f_n = -0.5*Ct_prime*abs(p_u_d_T)*p_u_d_T/p_thk        

            !write values to file 
                a0 = jt_total*dt_dim
                a1 = p_u_d
                a2 = p_u_d_T
                a3 = p_f_n
                a4 = -0.5*Ct_prime*(a2*a2*a2)*pi/(4.*sx*sy)  !power density

                if (s==1) then
                    call write_real_data('turbine_x1y1_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/))  
                elseif (s==2) then
                    call write_real_data('turbine_x2y1_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/))          
                elseif (s==3) then
                    call write_real_data('turbine_x3y1_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/)) 
                elseif (s==4) then
                    call write_real_data('turbine_x4y1_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/)) 
                elseif (s==5) then
                    call write_real_data('turbine_x1y2_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/)) 
                elseif (s==9) then
                    call write_real_data('turbine_x1y3_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/)) 
                elseif (s==13) then
                    call write_real_data('turbine_x1y4_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/)) 
                elseif (s==17) then
                    call write_real_data('turbine_x1y5_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/))
                elseif (s==21) then
                    call write_real_data('turbine_x1y6_forcing.dat', 'append', 'formatted', 5, (/a0,a1,a2,a3,a4/))                     
                endif          
                
            !write force to array that will be transferred via MPI    
            disk_force(s) = p_f_n
         
        enddo    
        
endif

!send total disk force to the necessary procs (with turbine_in_proc==.true.)
$if ($MPI)
    !############################################## 4  
        if (rank == 0) then
          
            do i=1,turbine_in_proc_cnt
                j = turbine_in_proc_array(i)
                call MPI_send( disk_force, nloc, MPI_rprec, j, 4, MPI_COMM_WORLD, ierr )
            enddo              
                        
        elseif (turbine_in_proc == .true.) then
            call MPI_recv( disk_force, nloc, MPI_rprec, 0, 4, MPI_COMM_WORLD, status, ierr )
        endif     
    !##############################################     
$endif 
    
!apply forcing to each node
if (turbine_in_proc == .true.) then

    !initialize forces to zero (only local values?)
    fx = 0.
    fy = 0.
    fz = 0.    

    do s=1,nloc
    
        !set pointers
            p_num_nodes => wind_farm_t%turbine_t(s)%num_nodes
            p_nhat1 => wind_farm_t%turbine_t(s)%nhat(1)
            p_nhat2 => wind_farm_t%turbine_t(s)%nhat(2)
            p_nhat3 => wind_farm_t%turbine_t(s)%nhat(3)    
            
        do l=1,p_num_nodes
            i2 = wind_farm_t%turbine_t(s)%nodes(l,1)
            j2 = wind_farm_t%turbine_t(s)%nodes(l,2)
            k2 = wind_farm_t%turbine_t(s)%nodes(l,3)
            ind2 = wind_farm_t%turbine_t(s)%ind(l)			
            fx(i2,j2,k2) = disk_force(s)*p_nhat1*ind2                            
            fy(i2,j2,k2) = disk_force(s)*p_nhat2*ind2   
            !fz(i2,j2,k2) = disk_force(s)*p_nhat3*ind2   !<< different points than fx,fy... check this
            fz(i2,j2,k2) = 0.5*disk_force(s)*p_nhat3*ind2
            fz(i2,j2,k2+1) = 0.5*disk_force(s)*p_nhat3*ind2
        enddo

    enddo
    
endif    


end subroutine turbines_forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_finalize()
implicit none
!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    fname4 = 'turbine_u_d_T.dat'
    open (unit = 2,file = fname4, status='unknown',form='formatted', action='write',position='rewind')
    do i=1,nloc
        write(2,*) wind_farm_t%turbine_t(i)%u_d_T    
    enddo    
    write(2,*) T_avg_dim
endif
    
end subroutine turbines_finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbine_vel_init(zo_high)
!  called from ic.f90 if initu, dns_bc, S_FLAG are all false.
!  this accounts for the turbines when creating the initial velocity profile.

use bottombc, only: zo_avg
implicit none

real(rprec), intent(inout) :: zo_high
real(rprec) :: cft,nu_w,exp_KE

!friction coefficient, cft
    cft = pi*Ct_noprime/(4.*sx*sy)
       
!wake viscosity
    nu_w = 28.*sqrt(0.5*cft)

!turbine friction height, Calaf, Phys. Fluids 22, 2010
    zo_high = height_all*(1.+0.5*dia_all/height_all)**(nu_w/(1.+nu_w))* &
      exp(-1.*(0.5*cft/(vonk**2) + (log(height_all/zo_avg* &
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
