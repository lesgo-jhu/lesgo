module turbines
use types,only:rprec
use param, only: nx,ny,nz,pi,L_x,L_y,L_z,dx,dy,dz,ld,jt_total,dt_dim
use stat_defs, only:wind_farm_t
use grid_defs, only:x,y,z
use io

implicit none

save
private

public :: turbines_init, turbines_forcing  

integer :: nloc
integer :: num_x,num_y
real :: space_x,space_y
real :: height_all,dia_all,thk_all,theta1_all,theta2_all
real :: Ct									        !thrust coefficient

character (64) :: fname0, fname, fname3, var_list
real(rprec), dimension(nx,ny,nz) :: large_node_array       !used for visualizing node locations

real :: eps											!epsilon used for disk velocity time-averaging

integer :: i,j,k,i2,j2,k2,b,l,s,nn,sx,sy,sz
integer :: imax,jmax,kmax,count_i,count_n,icp,jcp,kcp
integer :: min_i,max_i,min_j,max_j,min_k,max_k,cut
real(rprec) :: a0,a1,a2,a3,a4,a5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_init()
!PROBLEMS WITH MPI - MAKE SURE TURBINES ARE COMPLETELY WITHIN ONE PROC'S DOMAIN

!##############################  SET BY USER  ############################################
!set turbine parameters
!turbines are numbered as follows:
!   #1 = turbine nearest (x,y)=(0,0)
!   #2 = next turbine in the x-direction, etc.

	nloc = 4      		!number of turbines (locations) 
	allocate(wind_farm_t%turbine_t(nloc))

    !x,y-locations
        !for evenly-spaced turbines, not staggered
        if(.true.) then
            num_x = 2           !number of turbines in x-direction
            !space_x = 7.        !spacing in x-dir, multiple of DIA
            num_y = nloc/num_x  !number of turbines in y-direction            
            !space_y = 5.        !spacing in y-dir, multiple of DIA
            
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
            height_all = 100./1000.             !turbine height    [m]
            dia_all = 100./1000.	            !turbine diameter  [m] 
            thk_all = max(10./1000.,dx*1.01)	!turbine disk thickness [m]
            
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
	    Ct = 1.33		!thrust coefficient
!#########################################################################################

!find turbine nodes - including unfiltered ind, n_hat, num_nodes, and nodes for each turbine
	large_node_array = 0.
	call turbines_nodes(large_node_array)

	!to write the node locations to file
	  fname0 = 'node_loc_turbine.dat'
	  call write_tecplot_header_ND(fname0,'rewind', 4, (/nx, ny, nz/), '"x", "y", "z", "nodes"', 0, 1)
	  call write_real_data_3D(fname0, 'append','formatted', 1, nx, ny, nz, (/large_node_array/),x,y,z)

!1.smooth/filter indicator function                     
!2.associate new nodes with turbines                               
!3.normalize such that each turbine's ind integrates to turbine volume
    call turbines_filter_ind()

!set variables for time-averaging velocity 
    do k=1,nloc
	    wind_farm_t%turbine_t(k)%u_d_T = 0.
        wind_farm_t%turbine_t(k)%u_d_flag = 1
    enddo
	eps = 0.05 / (1 + 0.05)
	!time average over T(sec) = dt_dim/0.05 so that eps~0.05
	!T_min = dt_dim/0.05*1./60.	
    
    var_list = '"t (s)", "u_d", "u_d_T", "f_n"'
    call write_tecplot_header_xyline('turbine1_forcing.dat','rewind', var_list)   
    call write_tecplot_header_xyline('turbine2_forcing.dat','rewind', var_list) 
    call write_tecplot_header_xyline('turbine3_forcing.dat','rewind', var_list) 
	
end subroutine turbines_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_nodes(array)
!This subroutine locates nodes for each turbine and builds the arrays: ind, n_hat, num_nodes, and nodes
	
real :: R_t,rx,ry,rz,r,r_norm,r_disk
real(rprec), dimension(nx,ny,nz) :: array

real, pointer :: p_dia => null(), p_thk=> null(), p_theta1=> null(), p_theta2=> null()
real, pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3=> null() 
real, pointer :: p_xloc => null(), p_yloc=> null(), p_height=> null()

logical :: verbose
verbose = .true.

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
    		if (verbose) then
    		  write(*,*) 'R_t,dx,dy,dz,dthk',R_t,dx,dy,dz,p_thk
    		endif
    	imax = R_t/dx + 1
    	jmax = R_t/dy + 1
    	kmax = R_t/dz + 1

    !determine unit normal vector for each turbine	
    	p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    	p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    	p_nhat3 = sin(pi*p_theta2/180.)
    		if (verbose) then
    		  write(*,*) 'n_hat', p_nhat1, p_nhat2, p_nhat3
    		endif

    !determine nearest (i,j,k) to turbine center
    	icp = nint(p_xloc/dx)
    	jcp = nint(p_yloc/dy)
    	kcp = nint(p_height/dz + 0.5)
    		if (verbose) then
    		  write(*,*) 'turbine center',icp,jcp,kcp
    		  write(*,*) 'turbine center',x(icp),y(jcp),z(kcp)
    		endif

    !determine limits for checking i,j,k
    	min_i = max((icp-imax),1)
    	max_i = min((icp+imax),nx)
    	min_j = max((jcp-jmax),1)
    	max_j = min((jcp+jmax),ny)
    	min_k = max((kcp-kmax),1)
    	max_k = min((kcp+kmax),nz)
    		if (verbose) then
    		  write(*,*) 'i limits', min_i, max_i, dx*(max_i-min_i)
    		  write(*,*) 'j limits', min_j, max_j, dy*(max_j-min_j)
    		  write(*,*) 'k limits', min_k, max_k, dz*(max_k-min_k)
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
				rz = z(k) - p_height
				r = sqrt(rx*rx + ry*ry + rz*rz)
			!length projected onto unit normal for this turbine
				r_norm = abs(rx*p_nhat1 + ry*p_nhat2 + rz*p_nhat3)
			!(remaining) length projected onto turbine disk
				r_disk = sqrt(r*r - r_norm*r_norm)
			!if r_disk<R_t and r_norm within thk/2 from turbine -- this node is part of the turbine
				if ( (r_disk .LE. R_t) .AND. (r_norm .LE. p_thk/2.) ) then
					if (verbose) then
					  write(*,*) 'FOUND NODE', i,j,k
					endif
					array(i,j,k) = 1.
                    wind_farm_t%turbine_t(s)%ind(count_i) = 1.		
					wind_farm_t%turbine_t(s)%nodes(count_i,1) = i
					wind_farm_t%turbine_t(s)%nodes(count_i,2) = j
					wind_farm_t%turbine_t(s)%nodes(count_i,3) = k
                    count_n = count_n + 1
					count_i = count_i + 1
				endif
			enddo
		enddo
	enddo
	wind_farm_t%turbine_t(s)%num_nodes = count_n

    if (verbose) then
        write(*,*) 'Turbine ',s,': ',count_n,' nodes'
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

real(rprec), dimension(nx,ny,nz) :: out_a, g, g_shift, fg
real, dimension(nx,ny,nz) :: temp_array
real :: sumG,delta2,r2,sumA
real :: turbine_vol

logical :: verbose
verbose = .true.

!create convolution function, centered at (nx/2,ny/2,nz/2) and normalized
	if(wind_farm_t%ifilter==0) then		!0-> none
	  !
	elseif(wind_farm_t%ifilter==1) then  	!1-> cutoff/sharp spectral
	  !
	elseif(wind_farm_t%ifilter==2) then		!2-> Gaussian
	  delta2 = wind_farm_t%alpha**2 * (dx**2 + dy**2 + dz**2)
      do k=1,nz
    	  do j=1,ny
    	    do i=1,nx
    	      r2 = ((real(i)-nx/2.)*dx)**2 + ((real(j)-ny/2.)*dy)**2 + ((real(k)-nz/2.)*dz)**2
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
    do k=1,nz-1         !since k=1 and k=nz are the same point
	  do j=1,ny-1		!since j=1 and j=ny are the same point
	    do i=1,nx-1		!since i=1 and i=nx are the same point
	      sumG = sumG + g(i,j,k)*dx*dy*dz
	    enddo
	  enddo
    enddo
	g = g/sumG

    if (verbose) then
        write(*,*) 'Convolution function created.'
    endif

!display the convolution function
	if(.false.) then
	  write(*,*) 'Convolution function'
	  write(*,*) g
	endif
	write(*,*) 'integral of g(i,j,k): ',sumG

	!to write the data to file, centered at (i,j,k=nz/2)
	i=nx/2
	j=ny/2
    do k2=1,nz
	  do j2=1,ny
	    do i2=1,nx
	    g_shift(i2,j2,k2) = g( mod(i2-i+nx/2+(nx-1)-1,(nx-1))+1 , mod(j2-j+ny/2+(ny-1)-1,(ny-1))+1, k2)
	    enddo
	  enddo
    enddo
	  fname0 = 'convolution_function.dat'
	  call write_tecplot_header_ND(fname0,'rewind', 4, (/nx,ny,nz/), '"x","y","z","g"', 1, 1)
	  call write_real_data_3D(fname0, 'append', 'formatted', 1, nx, ny, nz, (/g_shift/),x,y,z)

    if (verbose) then
        write(*,*) 'Convolution function written to Tecplot file.'
    endif

!filter indicator function for each turbine
do b=1,nloc

    if (verbose) then
        write(*,*) 'Turbine Number ',b
    endif

    !create the input array (nx,ny,nz) from a list of included nodes
        temp_array = 0.
        do l=1,wind_farm_t%turbine_t(b)%num_nodes
            i2 = wind_farm_t%turbine_t(b)%nodes(l,1)
            j2 = wind_farm_t%turbine_t(b)%nodes(l,2)
            k2 = wind_farm_t%turbine_t(b)%nodes(l,3)	
            temp_array(i2,j2,k2) = wind_farm_t%turbine_t(b)%ind(l)

            !if (verbose) then
                !write(*,*) 'Writing node to temp_array ',i2,j2,k2
            !endif
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

        if (verbose) then
            write(*,*) 'search over: ',min_i-cut,max_i+cut,min_j-cut,max_j+cut,min_k-cut,max_k+cut
        endif

        do k=max(min_k-cut,1),min(max_k+cut,nz)    !only compute for nodes near the turbine
    	do j=max(min_j-cut,1),min(max_j+cut,ny)
    	do i=max(min_i-cut,1),min(max_i+cut,nx)
		
          do k2=k-wind_farm_t%trunc,k+wind_farm_t%trunc     !currently using truncated Gaussian
    	  do j2=j-wind_farm_t%trunc,j+wind_farm_t%trunc
    	  do i2=i-wind_farm_t%trunc,i+wind_farm_t%trunc

            if ( (i2>0).and.(j2>0).and.(k2>0).and.(i2<=nx).and.(j2<=ny).and.(k2<=nz) ) then
                sx = mod(i2-i+nx/2+(nx-1)-1,(nx-1))+1
    	        sy = mod(j2-j+ny/2+(ny-1)-1,(ny-1))+1       
                sz = k2-k+nz/2       !since no spectral BCs in z-direction
            
                if (sx < 1) then
                    fg(i2,j2,k2) = 0.
                elseif (sx > nx) then
                    fg(i2,j2,k2) = 0.
                elseif( sy < 1 ) then
                    fg(i2,j2,k2) = 0.
                elseif( sy > ny) then
                    fg(i2,j2,k2) = 0.
                elseif( sz < 1) then
                    fg(i2,j2,k2) = 0.
                elseif( sz > nz ) then
                    fg(i2,j2,k2) = 0.
                else
                    fg(i2,j2,k2) = temp_array(i2,j2,k2)*g(sx,sy,sz)
                    out_a(i,j,k) = out_a(i,j,k) + fg(i2,j2,k2)*dx*dy*dz
                endif    

            endif	
    	    
    	  enddo
    	  enddo
          enddo
    	enddo
    	enddo
        enddo

        if (verbose) then
            write(*,*) 'Convolution complete for turbine ',b
        endif

    !normalize this "indicator function" such that it integrates to turbine volume
	sumA = 0.
    do k=1,nz-1         !since k=1 and k=nz are the same point
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

    if (verbose) then
        write(*,*) 'sumA,turbine_vol = ',sumA,turbine_vol
            if (b==1) then
        	  fname3 = 'convolution_out.dat'
        	  call write_tecplot_header_ND(fname3,'rewind', 4, (/nx,ny,nz/), '"x","y","z","out"', 1, 1)
        	  call write_real_data_3D(fname3, 'append', 'formatted', 1, nx, ny, nz, (/out_a/),x,y,z)
            endif
    endif

    !update num_nodes, nodes, and ind for this turbine
    wind_farm_t%turbine_t(b)%ind = 0.
    wind_farm_t%turbine_t(b)%nodes = 0.
    wind_farm_t%turbine_t(b)%num_nodes = 0.
    count_n = 0
    count_i = 1
	do k=1,nz
		do j=1,ny
			do i=1,nx
				if (out_a(i,j,k) > wind_farm_t%filter_cutoff) then
                    wind_farm_t%turbine_t(b)%ind(count_i) = out_a(i,j,k)		
					wind_farm_t%turbine_t(b)%nodes(count_i,1) = i
					wind_farm_t%turbine_t(b)%nodes(count_i,2) = j
					wind_farm_t%turbine_t(b)%nodes(count_i,3) = k
                    count_n = count_n + 1
					count_i = count_i + 1
				endif
			enddo
		enddo
	enddo
	wind_farm_t%turbine_t(b)%num_nodes = count_n
    if (verbose) then
        write(*,*) 'Turbine number ',b,' has ',count_n,' nodes' 
    endif

enddo

end subroutine turbines_filter_ind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine turbines_forcing()
use sim_param, only: u,v,w
use immersedbc, only: fx,fy,fz

real, pointer :: p_u_d => null(), p_nhat1=> null(), p_nhat2=> null(), p_nhat3=> null()
real, pointer :: p_u_d_T => null(), p_dia => null(), p_thk=> null(), p_f_n => null()
integer, pointer :: p_u_d_flag=> null(), p_num_nodes=> null()

real :: ind2

!for each turbine:        
    do s=1,nloc            
         
    !set pointers
        p_u_d => wind_farm_t%turbine_t(s)%u_d   
        p_u_d_T => wind_farm_t%turbine_t(s)%u_d_T   
        p_num_nodes => wind_farm_t%turbine_t(s)%num_nodes
        p_nhat1 => wind_farm_t%turbine_t(s)%nhat(1)
        p_nhat2 => wind_farm_t%turbine_t(s)%nhat(2)
        p_nhat3 => wind_farm_t%turbine_t(s)%nhat(3)
        p_u_d_flag => wind_farm_t%turbine_t(s)%u_d_flag
        p_f_n => wind_farm_t%turbine_t(s)%f_n               
        p_dia => wind_farm_t%turbine_t(s)%dia 
        p_thk => wind_farm_t%turbine_t(s)%thk

    !calculate total disk-averaged velocity for each turbine (current,instantaneous)    
    !u_d and u_d_T are velocities in the normal direction	 
    !for calculating average velocity, nodes are weighted equally    
        call interp_to_uv_grid(w, w_uv, w_uv_tag) 
        p_u_d = 0.
        do l=1,p_num_nodes
            i2 = wind_farm_t%turbine_t(s)%nodes(l,1)
            j2 = wind_farm_t%turbine_t(s)%nodes(l,2)
            k2 = wind_farm_t%turbine_t(s)%nodes(l,3)	
            p_u_d = p_u_d + (p_nhat1*u(i2,j2,k2) + p_nhat2*v(i2,j2,k2) + p_nhat3*w_uv(i2,j2,k2))/p_num_nodes
        enddo  

    !add this current value to the "running average" (first order relaxation)
        if (p_u_d_flag) then
            p_u_d_T = p_u_d
            p_u_d_flag = 0
        else
            p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d
        endif

    !calculate total thrust force for each turbine  (per unit mass)
    !force is normal to the surface (calc from u_d_T, normal to surface)
        p_f_n = -0.5*Ct*abs(p_u_d_T)*p_u_d_T   !/p_thk            !<<<<<< try using p_u_d instead of p_u_d_T

        if (s<4) then
            a0 = jt_total*dt_dim
            a1 = p_u_d
            a2 = p_u_d_T
            a3 = p_f_n
        endif
        if (s==1) then
            call write_real_data('turbine1_forcing.dat', 'append', 4, (/a0,a1,a2,a3/))  
        elseif (s==2) then
            call write_real_data('turbine2_forcing.dat', 'append', 4, (/a0,a1,a2,a3/))          
        elseif (s==3) then
            call write_real_data('turbine3_forcing.dat', 'append', 4, (/a0,a1,a2,a3/))    
        endif
        
    !apply forcing to each node
        do l=1,p_num_nodes
            i2 = wind_farm_t%turbine_t(s)%nodes(l,1)
            j2 = wind_farm_t%turbine_t(s)%nodes(l,2)
            k2 = wind_farm_t%turbine_t(s)%nodes(l,3)
            ind2 = wind_farm_t%turbine_t(s)%ind(l)			
            fx(i2,j2,k2) = p_f_n*p_nhat1*ind2                            
            fy(i2,j2,k2) = p_f_n*p_nhat2*ind2   
            fz(i2,j2,k2) = p_f_n*p_nhat3*ind2   !<< different points than fx,fy... check this
            !fz(i2,j2,k2) = 0.5*p_f_n*p_nhat3*ind2
            !fz(i2,j2,k2+1) = 0.5*p_f_n*p_nhat3*ind2
            
            !if (s==1) then
            !    a0 = jt_total*dt_dim
            !    a4 = l
            !    a5 = p_f_n*p_nhat1*ind2             
            !endif
        enddo

    enddo
	
    if(.false.) then        !need to make this occur only at last time step
    !to write the data to file
        fname = 'force_turbine.dat'
        call write_tecplot_header_ND(fname,'rewind', 4, (/nx, ny, nz/), '"x", "y", "z", "force_x"', 1, 1)
        call write_real_data_3D(fname, 'append','formatted', 1, nx, ny, nz, (/fx(1:nx,:,:)/),x(1:nx),y,z)
    endif

end subroutine turbines_forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module turbines
