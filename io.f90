module io
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, write_inflow_file, path,  &
                  USE_MPI, coord, rank, nproc, jt_total, total_time, total_time_dim
use param, only : cumulative_time, fcumulative_time
use sim_param, only : w, dudz, dvdz
use messages
use strmod
use sgsmodule,only:Cs_opt2

implicit none

save
private

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

!!$public openfiles,output_loop,output_final,                   &
!!$     inflow_write, avg_stats
public jt_total, openfiles, inflow_read, inflow_write, output_loop, output_final
public mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2
public w_uv, dudz_uv, w_uv_tag, dudz_uv_tag, interp_to_uv_grid, stats_init
!public write_tecplot_header_xyline, write_tecplot_header_ND
!public write_real_data, write_real_data_1D, write_real_data_2D, write_real_data_3D

character (*), parameter :: mod_name = 'io'

integer,parameter::jx_pls=1,jx_ple=1,width=1
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width+1
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz):: &
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2

real(rprec), dimension(ld,ny,$lbz:nz) :: w_uv, dudz_uv, dvdz_uv

integer :: w_uv_tag, dudz_uv_tag, dvdz_uv_tag
                   
!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine interp_to_uv_grid(var,var_uv,tag)
!**********************************************************************
!  This function interpolates the array var, which resides on the w-grid,
!  onto the uv-grid variable var_uv using linear interpolation. It is 
!  important to note that message passing is required for MPI cases and 
!  all processors must call this routine. If this subroutine is call from a 
!  only a subset of the total processors, the code will hang due to the usage
!  of the syncronous send/recv functions and certain processors waiting
!  to recv data but it never gets there.
!
!  NOTE: It is assumed that the size of var and var_uv are the same as the
!  coord (processor) domain and that k=nz-1 (k=0) and k=1 (k=nz) are overlap
!  nodes - no interpolation is performed for k=0 and k=nz
!
use types, only : rprec
use param,only : nz,ld,jt
use sim_param, only : w, dudz
$if ($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

real(rprec), dimension(:,:,:), intent(IN) :: var
real(rprec), dimension(:,:,:), intent(OUT) :: var_uv
integer, intent(INOUT) :: tag
!real(rprec), dimension(2) :: var
integer :: lbx,ubx,lby,uby,lbz,ubz
integer :: i,j,k

character (*), parameter :: sub_name = mod_name // '.interp_to_uv_grid'

if(tag == jt) then
$if ($VERBOSE)
  call mesg(sub_name, 'Interpolation already performed for current time step')
$endif
  return
endif
  
lbx=lbound(var,1); ubx=ubound(var,1)
lby=lbound(var,2); uby=ubound(var,2)
lbz=lbound(var,3); ubz=ubound(var,3)

do k=lbz+1,ubz-1
  do j=lby,uby
    do i=lbx,ubx
      var_uv(i,j,k) = 0.5 * (var(i,j,k+1) + var(i,j,k))
    enddo
  enddo
enddo

$if ($MPI)

!  Take care of top "physical" boundary
if(coord == nproc - 1) var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

!  Sync all overlapping data
call mpi_sync_real_array( var_uv, MPI_SYNC_DOWNUP )

$else

!  Take care of top "physical" boundary
var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

$endif
  
tag = jt ! Set identifying tag 

return 

!!$if($MPI)
!deallocate(buf)
!!$endif

end subroutine interp_to_uv_grid

!**********************************************************************
subroutine openfiles()
!**********************************************************************
$if($CFL_DT)
use param, only : dt, cfl_f
$endif
use sim_param,only:path
implicit none

logical :: exst

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    
    $if($CFL_DT)
    read(1, *) jt_total, total_time, total_time_dim, dt, cfl_f
    $else
    read (1, *) jt_total, total_time, total_time_dim
    $endif
    
    close (1)
  else  !--assume this is the first run on cumulative time
    if( .not. USE_MPI .or. ( USE_MPI .and. coord == 0 ) ) then
      write (*, *) '--> Assuming jt_total = 0, total_time = 0., total_time_dim = 0.'
    endif
    jt_total = 0
    total_time = 0.
    total_time_dim = 0._rprec
  end if

end if

!!--see also energy () subr. for flushing this file
!if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
!  open(13, file=path//'output/check_ke.out', position='append')
!end if


end subroutine openfiles

!**********************************************************************
subroutine output_loop(jt)
!**********************************************************************
!
!  This subroutine is called every time step and acts as a driver for 
!  computing statistics and outputing instantaneous data. No actual
!  calculations are performed here.
!
use param, only : tavg_calc, tavg_nstart, tavg_nend
use param, only : spectra_calc, spectra_nstart, spectra_nend
use param, only : point_calc, point_nstart, point_nend, point_nskip
use param, only : domain_calc, domain_nstart, domain_nend, domain_nskip
use param, only : xplane_calc, xplane_nstart, xplane_nend, xplane_nskip
use param, only : yplane_calc, yplane_nstart, yplane_nend, yplane_nskip
use param, only : zplane_calc, zplane_nstart, zplane_nend, zplane_nskip

implicit none
integer,intent(in)::jt

!  Determine if time summations are to be calculated
if(tavg_calc) then
!  Check if we are in the time interval for running summations
  if(jt >= tavg_nstart .and. jt <= tavg_nend) then
    if(jt == tavg_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Starting running time summation from ', tavg_nstart, ' to ', tavg_nend
        write(*,*) '-------------------------------'   
      endif

      call tavg_init()

    endif
!  Compute running summations
    call tavg_compute ()
  endif 
  
endif

if( spectra_calc ) then
  !  Check if we are in the time interval for running summations
  if(jt >= spectra_nstart .and. jt <= spectra_nend) then
    if(jt == spectra_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Starting running spectra calculations from ', spectra_nstart, ' to ', spectra_nend
        write(*,*) '-------------------------------'
      endif

      call spectra_init()

    endif
!  Compute running summations
    call spectra_compute ()
  endif

endif


!  Determine if instantaneous point velocities are to be recorded
if(point_calc) then
  if(jt >= point_nstart .and. jt <= point_nend .and. ( jt == point_nstart .or. mod(jt,point_nskip)==0) ) then
    if(jt == point_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then   
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous point velocities from ', point_nstart, ' to ', point_nend
        write(*,"(1a,i9)") 'Iteration skip:', point_nskip
        write(*,*) '-------------------------------'
      endif
    endif
    call inst_write(1)
  endif
endif
  
!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
  if(jt >= domain_nstart .and. jt <= domain_nend .and. ( jt == domain_nstart .or. mod(jt,domain_nskip)==0) ) then
    if(jt == domain_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous domain velocities from ', domain_nstart, ' to ', domain_nend
        write(*,"(1a,i9)") 'Iteration skip:', domain_nskip
        write(*,*) '-------------------------------'
      endif

    endif
    call inst_write(2)
  endif
endif 

!  Determine if instantaneous x-plane velocities are to be recorded
if(xplane_calc) then
  if(jt >= xplane_nstart .and. jt <= xplane_nend .and. ( jt == xplane_nstart .or. mod(jt,xplane_nskip)==0) ) then
    if(jt == xplane_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous x-plane velocities from ', xplane_nstart, ' to ', xplane_nend
        write(*,"(1a,i9)") 'Iteration skip:', xplane_nskip
        write(*,*) '-------------------------------'
      endif

    endif
    call inst_write(3)
  endif
endif  

!  Determine if instantaneous y-plane velocities are to be recorded
if(yplane_calc) then
  if(jt >= yplane_nstart .and. jt <= yplane_nend .and. ( jt == yplane_nstart .or. mod(jt,yplane_nskip)==0) ) then
    if(jt == yplane_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous y-plane velocities from ', yplane_nstart, ' to ', yplane_nend
        write(*,"(1a,i9)") 'Iteration skip:', yplane_nskip
        write(*,*) '-------------------------------'
      endif
    endif
    call inst_write(4)
  endif
endif    

!  Determine if instantaneous z-plane velocities are to be recorded
if(zplane_calc) then
  if(jt >= zplane_nstart .and. jt <= zplane_nend .and. ( jt == zplane_nstart .or. mod(jt,zplane_nskip)==0) ) then
    if(jt == zplane_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous z-plane velocities from ', zplane_nstart, ' to ', zplane_nend
        write(*,"(1a,i9)") 'Iteration skip:', zplane_nskip
        write(*,*) '-------------------------------'
      endif
    endif
    call inst_write(5)
  endif
endif


return
end subroutine output_loop



!**********************************************************************
subroutine inst_write(itype)
!**********************************************************************
!  This subroutine writes the instantaneous values
!  at specified i,j,k locations
use functions, only : linear_interp, trilinear_interp
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use grid_defs, only : x,y,z,zw
use sim_param, only : u,v,w,dudx,dvdy,dwdz
use stat_defs, only : xplane_t, yplane_t, zplane_t, point_t
$if($MPI)
use mpi
use param, only : ld, ny, nz, MPI_RPREC, down, up, comm, status, ierr
$endif

$if($LVLSET)
use level_set_base, only : phi
use immersedbc, only : fx, fy, fz, fxa, fya, fza
$endif

use param, only : jt_total, dt_dim, nx, ny, nz,dx,dy,dz,z_i,L_x,L_y,L_z,coord
implicit none

include 'tecio.h'

integer, intent(IN) :: itype

character(25) :: cl, ct
character (64) :: fname
$if($MPI)
character(64) :: temp
$endif
character(128) :: var_list
integer :: n, i, j, k, nvars

real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi, fx_tot, fy_tot, fz_tot

$if($DEBUG)
real(rprec), allocatable, dimension(:,:,:) :: divvel
$endif

$if($PGI)
real(rprec), allocatable, dimension(:) :: u_inter
$endif

! real(rprec) :: dnx, dny, dnz

!  Write point data; assumes files have been opened properly
!  in stats_init

!  Make sure w has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)

!$if($MPI)
!!  Sync fx; can not use mpi_sync_real_array since its not allocated from 0 -> nz
!call mpi_sendrecv (fx(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!  fx(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!  comm, status, ierr)
!$endif

if(itype==1) then

  do n=1,point_nloc

    !  For parallel runs check if data is on correct proc
    $if ($MPI)
    if(point_t(n) % coord == coord) then
    $endif

    call write_real_data(point_t(n) % fname, 'append', 'formatted', 4, (/ total_time, &
      trilinear_interp(u(1:nx,1:ny,1:nz), 1, point_loc(n)%xyz), &
      trilinear_interp(v(1:nx,1:ny,1:nz), 1, point_loc(n)%xyz), &
      trilinear_interp(w_uv(1:nx,1:ny,1:nz), 1, point_loc(n)%xyz) /))


    $if ($MPI)
    endif
    $endif

  enddo
!  Instantaneous write for entire domain
elseif(itype==2) then

!  Convert total iteration time to string
  write(ct,*) jt_total
!  Open file which to write global data
  write (fname,*) 'output/vel.', trim(adjustl(ct)),'.dat'
  fname = trim(adjustl(fname))

  $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
  $endif

  $if($LVLSET)
  !write(7,*) 'variables = "x", "y", "z", "u", "v", "w", "phi"';
  var_list = '"x", "y", "z", "u", "v", "w", "phi"'
  nvars = 7
  $else
  !write(7,*) 'variables = "x", "y", "z", "u", "v", "w"';
  var_list = '"x", "y", "z", "u", "v", "w"'
  nvars = 6
  $endif
  
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), var_list, coord, 2, real(total_time,4))
  
  $if($LVLSET)
  call write_real_data_3D(fname, 'append', 'formatted', 3, nx, ny, nz, &
    (/ u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w_uv(1:nx,1:ny,1:nz) /), & 
    4, x, y, z(1:nz))
  call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny, nz, &
    (/ phi(1:nx,1:ny,1:nz) /), 4)
  
  $else
    $if($PGI)
    allocate(u_inter(nx*ny*nz*3))
    icount=0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          icount=icount+1
          u_inter(icount) = u(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=1,ny
        do i=1,nx
          icount=icount+1
          u_inter(icount) = v(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=1,ny
        do i=1,nx
          icount=icount+1
          u_inter(icount) = w(i,j,k)
        enddo
      enddo
    enddo
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
      u_inter, 4, x, y, z(1:nz))
    deallocate(u_inter)
    $else
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
      (/ u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w_uv(1:nx,1:ny,1:nz) /), &
      4, x, y, z(1:nz))
    $endif
  $endif

  !  Output Instantaneous Force Field for RNS Simulations
  !  Still need to put fz on uv grid may need a better way
  $if($LVLSET)
    $if($MPI)
    !  Sync fx; can't use mpi_sync_real_array since its not allocated from 0 -> nz
    call mpi_sendrecv (fx(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
                       fx(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
                       comm, status, ierr)
    call mpi_sendrecv (fy(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
                       fy(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
                       comm, status, ierr)
    call mpi_sendrecv (fz(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
                       fz(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
                       comm, status, ierr)
    call mpi_sendrecv (fxa(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
                       fxa(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
                       comm, status, ierr)
    call mpi_sendrecv (fya(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
                       fya(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
                       comm, status, ierr)
    call mpi_sendrecv (fza(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
                       fza(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
                       comm, status, ierr)
    $endif

    !  Sum both the induced and applied forces
    allocate(fx_tot(nx,ny,nz), fy_tot(nx,ny,nz), fz_tot(nx,ny,nz))
    fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
    fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
    fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)

    !  Open file which to write global data
    write (fname,*) 'output/force.', trim(adjustl(ct)),'.dat'
    fname = trim(adjustl(fname))

    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    var_list = '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>", "phi"'
    nvars = 7
    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), var_list, coord, 2, real(total_time,4))
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx, ny,nz, &
      (/ fx_tot, fy_tot, fz_tot /), 4, x, y, z(1:nz))
    call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny,nz, &
      (/ phi(1:nx,1:ny,1:nz) /), 4)

    deallocate(fx_tot, fy_tot, fz_tot)
  
  $endif

  $if($DEBUG)
  if(DEBUG) then
  !  Output divergence of velocity field
  allocate(divvel(nx,ny,nz))
  divvel=dudx(1:nx,1:ny,1:nz)+dvdy(1:nx,1:ny,1:nz)+dwdz(1:nx,1:ny,1:nz)

  !  Open file which to write global data
  write (fname,*) 'output/divvel.', trim(adjustl(ct)),'.dat'
  fname = trim(adjustl(fname))

  $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
  $endif

  $if($LVLSET)
  var_list = '"x", "y", "z", "divvel", "phi"'
  nvars = 5
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), var_list, coord, 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 2, nx, ny,nz, &
  (/ divvel, phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
  $else
  var_list = '"x", "y", "z", "divvel"'
  nvars = 4
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), var_list, coord, 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny,nz, &
  (/ divvel /), 4, x, y, z(1:nz))
  $endif

  deallocate(divvel)
  endif
  $endif


!  Write instantaneous x-plane values
elseif(itype==3) then

  allocate(ui(1,ny,nz), vi(1,ny,nz), wi(1,ny,nz))

!  Loop over all xplane locations
  do i=1,xplane_nloc

    write(cl,'(F9.4)') xplane_loc(i)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/vel.x-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    $if ($MPI)
!  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny+1, Nz /), &
      '"x", "y", "z", "u", "v", "w"', coord, 2, real(total_time,4))  
  
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(u(xplane_t(i) % istart,j,k), &
          u(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)
        vi(1,j,k) = linear_interp(v(xplane_t(i) % istart,j,k), &
          v(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)
        wi(1,j,k) = linear_interp(w_uv(xplane_t(i) % istart,j,k), &
          w_uv(xplane_t(i) % istart+1,j,k), dx, &
          xplane_t(i) % ldiff)
      enddo
    enddo

    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))     

    $if($LVLSET)

    write(fname,*) 'output/force.x-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    $if ($MPI)
    !  For MPI implementation
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
    $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny+1, Nz/), &
      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', coord, 2, real(total_time,4))

    !  Sum both induced forces, f{x,y,z}, and applied forces, f{x,y,z}a
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(fx(xplane_t(i) % istart,j,k), &
          fx(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff) + &
          linear_interp(fxa(xplane_t(i) % istart,j,k), &
          fxa(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)

        vi(1,j,k) = linear_interp(fy(xplane_t(i) % istart,j,k), &
          fy(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff) + &
          linear_interp(fya(xplane_t(i) % istart,j,k), &
          fya(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)

        wi(1,j,k) = linear_interp(fz(xplane_t(i) % istart,j,k), &
          fz(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff) + &
          linear_interp(fza(xplane_t(i) % istart,j,k), &
          fza(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)
      enddo
    enddo


    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))


    $endif
    
  enddo   
  
  deallocate(ui,vi,wi)

  
!  Write instantaneous y-plane values
elseif(itype==4) then
  
  allocate(ui(nx,1,nz), vi(nx,1,nz), wi(nx,1,nz))
  
!  Loop over all yplane locations
  do j=1,yplane_nloc

    write(cl,'(F9.4)') yplane_loc(j)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/vel.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    $if ($MPI)
!  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, 1, Nz/), &
      '"x", "y", "z", "u", "v", "w"', coord, 2, real(total_time,4)) 
    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(u(i,yplane_t(j) % istart,k), &
          u(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)
        vi(i,1,k) = linear_interp(v(i,yplane_t(j) % istart,k), &
          v(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)
        wi(i,1,k) = linear_interp(w_uv(i,yplane_t(j) % istart,k), &
          w_uv(i,yplane_t(j) % istart+1,k), dy, &
          yplane_t(j) % ldiff)
          
      enddo
    enddo
    
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))    
  
  $if($LVLSET)
  
    write(fname,*) 'output/force.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    $if ($MPI)
!  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
  
    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, 1, Nz/), &
      '"x", "y", "z", "fx", "fy", "fz"', coord, 2, real(total_time,4))  
  
    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(fx(i,yplane_t(j) % istart,k), &
          fx(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff) + &
          linear_interp(fxa(i,yplane_t(j) % istart,k), &
          fxa(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)

        vi(i,1,k) = linear_interp(fy(i,yplane_t(j) % istart,k), &
          fy(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff) + &
          linear_interp(fya(i,yplane_t(j) % istart,k), &
          fya(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)

        wi(i,1,k) = linear_interp(fz(i,yplane_t(j) % istart,k), &
          fz(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff) + &
          linear_interp(fza(i,yplane_t(j) % istart,k), &
          fza(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)

      enddo
    enddo
    
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))       
    
    $endif

  enddo  

  deallocate(ui,vi,wi)
  
!  Write instantaneous z-plane values
elseif(itype==5) then

  allocate(ui(nx,ny,1), vi(nx,ny,1), wi(nx,ny,1))

!  Loop over all zplane locations
  do k=1,zplane_nloc
  
    $if ($MPI)
    if(zplane_t(k) % coord == coord) then
    $endif

    write(cl,'(F9.4)') zplane_loc(k)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/vel.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, Ny+1, 1/), &
      '"x", "y", "z", "u", "v", "w"', coord, 2, real(total_time,4)) 

    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(u(i,j,zplane_t(k) % istart),u(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff)
        vi(i,j,1) = linear_interp(v(i,j,zplane_t(k) % istart),v(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff)
        wi(i,j,1) = linear_interp(w_uv(i,j,zplane_t(k) % istart), &
          w_uv(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff)

      enddo
    enddo

    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
    (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /))   
    
    $if($LVLSET)
    
    write(fname,*) 'output/force.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, Ny+1, 1/), &
      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', coord, 2, real(total_time,4))

    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(fx(i,j,zplane_t(k) % istart),fx(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff) + &
          linear_interp(fxa(i,j,zplane_t(k) % istart),fxa(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff)
        vi(i,j,1) = linear_interp(fy(i,j,zplane_t(k) % istart),fy(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff) + &
          linear_interp(fya(i,j,zplane_t(k) % istart),fya(i,j,zplane_t(k) % istart+1), &
          dz, zplane_t(k) % ldiff) 
        wi(i,j,1) = linear_interp(fz(i,j,zplane_t(k) % istart), &
          fz(i,j,zplane_t(k) % istart+1), dz, zplane_t(k) % ldiff) + &
          linear_interp(fza(i,j,zplane_t(k) % istart), &
          fza(i,j,zplane_t(k) % istart+1), dz, zplane_t(k) % ldiff)

      enddo
    enddo
 
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
    (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /) )      
    
    $endif

    $if ($MPI)
    endif
    $endif

  enddo  
  
  deallocate(ui,vi,wi)

else
  write(*,*) 'Error: itype not specified properly to inst_write!'
  stop
endif
return
end subroutine inst_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (lun)
use param, only : nz, S_FLAG
use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta
use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
use scalars_module, only : RHS_T
implicit none

integer, intent (in) :: lun

!---------------------------------------------------------------------

if (S_FLAG) then !WITH SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
              Cs_opt2, F_LM, F_MM, F_QN, F_NN, theta, RHS_T
else ! No SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
              Cs_opt2, F_LM, F_MM, F_QN, F_NN
end if

end subroutine checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by
!  inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final(jt, lun_opt)
use stat_defs, only : tavg_t, point_t
use param, only : tavg_calc, point_calc, point_nloc, spectra_calc

$if($CFL_DT)
use param, only : dt, cfl
$endif

implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun
integer, parameter :: lun_default = 11
integer::i,fid
integer :: lun

logical :: opn

!---------------------------------------------------------------------

if (present (lun_opt)) then
  lun = lun_opt
else
  lun = lun_default
end if

inquire (unit=lun, opened=opn)

if (.not. opn) then
  write (*, *) 'output_final: lun=', lun, ' is not open'
  stop
end if

rewind (lun)

call checkpoint (lun)

close (lun)


!  Update total_time.dat after simulation
if ((cumulative_time) .and. (lun == lun_default)) then
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)

    $if($CFL_DT)
    write(1, *) jt_total, total_time, total_time_dim, dt, cfl
    $else
    write(1, *) jt_total, total_time, total_time_dim
    $endif

    close(1)

  end if
end if

!  Check if average quantities are to be recorded
if(tavg_calc) call tavg_finalize()

!  Check if spectra is to be computed
if(spectra_calc) call spectra_finalize()
 
!  Close instantaneous velocity files
if(point_calc) then

  do i=1,point_nloc
    $if ($MPI)
    if(point_t(i)%coord == coord) then
    $endif

    fid=3000*i
    close(fid)

    $if ($MPI)
    endif
    $endif

  enddo

endif

return
end subroutine output_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_write ()
use param, only : jt_total, jt_start_write, buff_end,  &
                  read_inflow_file, write_inflow_file
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_write'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: field_file = 'output/inflow.vel.out'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80
integer, parameter :: field_lun = 81

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif
character (64) :: fname

integer, save :: rec = 0
integer :: nrec
integer :: iolen
integer :: iend, iend_w

logical, save :: initialized = .false.
logical :: opn, exst

!---------------------------------------------------------------------

!--option check
if ( read_inflow_file .and. write_inflow_file ) then
  write (*, *) sub // ': cannot have read_inflow_file and write_inflow_file'
  stop
end if

!--check consistency with inflow_cond
iend = floor (buff_end * nx + 1._rprec)
iend_w = modulo (iend - 1, nx) + 1

if (.not. initialized) then

  inquire ( unit=lun, exist=exst, opened=opn )
  if ( .not. exst ) then
    write (*, *) sub // ': lun = ', lun, ' does not exist'
    stop
  end if
  if (opn) then
    write (*, *) sub // ': lun = ', lun, ' is already open'
    stop
  end if

  if ( USE_MPI ) then
      write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
  else
      write ( fname, '(a)' ) inflow_file
  end if
  
  inquire ( file=fname, exist=exst, opened=opn )
  if (exst .and. opn) then
    write (*, *) sub // ': file = ', trim (fname), ' is already open'
    stop
  end if
  
  !--figure out the record length
  inquire (iolength=iolen) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)

  !--always add to existing inflow_file
  !--inflow_file records always start at 1
  if ( exst ) then

      !--figure out the number of records already in file
      call len_da_file (fname, iolen, nrec)

      write (*, *) sub // ': #records in ' // trim (fname) // '= ', nrec

      rec = nrec

  else

      rec = 0

  end if
  
  !--using direct-access file to allow implementation of 'inflow recycling'
  !  more easily
  !--may want to put in some simple checks on ny, nz
  open (unit=lun, file=fname, access='direct', action='write',  &
        recl=iolen)

  initialized = .true.

end if

if (jt_total == jt_start_write) then  !--write entire flow field out
  inquire (unit=field_lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    open (unit=field_lun, file=field_file, form='unformatted')
    call output_final (jt_total, field_lun)
  else
    write (*, *) sub // ': problem opening ' // field_file
    stop
  end if
end if

if (jt_total >= jt_start_write) then
  rec = rec + 1
  write (unit=lun, rec=rec) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
  $if ($DEBUG)
  if ( DEBUG ) write (*, *) sub // ': wrote record ', rec
  $endif
end if

end subroutine inflow_write

!**********************************************************************
subroutine inflow_read ()
!**********************************************************************
use param, only : ny, nz, pi, nsteps, jt_total, buff_end
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_read'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: debug_file = 'inflow_read_debug.dat'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80  !--inflow_write for now
integer, parameter :: lun_DEBUG = 88

integer, parameter :: l_blend = 300  !--length of blending zone (recycling)
                                     !--should correspond to integral scale
                                     !--this is number of t-steps
logical, parameter :: recycle = .false.

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
character (32) :: fmt
$endif
character (64) :: fname

!--check for consistency with sim_param here
!--could define a fortran integer lbz in sim_param, and make it visible
!  here, however, this may have complications elsewhere where the name lbz
!  is used.
$if ( $MPI )
    $define $lbz 0
$else
    $define $lbz 1
$endif

$if($DEBUG)
integer :: jy, jz
$endif

integer :: iend, iend_w
integer :: iolen
integer, save :: rec
integer, save :: nrec
integer :: recp

$if($DEBUG)
logical, save :: init_DEBUG = .false.
$endif
logical, save :: initialized = .false.
logical :: exst, opn

real (rprec) :: wgt

real (rprec) :: u_tmp(ny, $lbz:nz), v_tmp(ny, $lbz:nz), w_tmp(ny, $lbz:nz)

!---------------------------------------------------------------------

iend = floor ( buff_end * nx + 1.0_rprec )
iend_w = modulo ( iend - 1, nx ) + 1

if ( .not. initialized ) then

    inquire ( unit=lun, exist=exst, opened=opn )
    if ( .not. exst ) then
        write (*, *) sub // ': lun = ', lun, ' does not exist'
        stop
    end if
    if ( opn ) then
        write (*, *) sub // ': lun = ', lun, ' is already open'
        stop
    end if

    if ( USE_MPI ) then
        write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
    else
        write ( fname, '(a)' ) inflow_file
    end if
    
    inquire ( file=fname, exist=exst, opened=opn )
    if ( exst ) then
        if ( opn ) then
            write (*, *) sub // ': file = ', fname, ' is already open'
            stop
        end if
    else
        write (*, *) sub // ': file = ', fname, ' does not exist'
        stop
    end if

    !--can only reach this point if exst and .not. opn
  
    !--figure out the record length
    inquire ( iolength=iolen ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
  
    !--figure out the number of records
    call len_da_file ( fname, iolen, nrec )

    write (*, *) sub // ': number of records = ', nrec

    if ( recycle ) then
        !--check minimum length requirement
        !  checks that there are some points that will be non-blended
        
        if ( 2 * (l_blend - 1) > nrec ) then
            write (*, *) sub // ': ', fname, 'is too short to recycle'
            stop
        end if
    end if

    open ( unit=lun, file=fname, access='direct', action='read',  &
           recl=iolen )

    !--file always starts a record 1, but in continued runs, we may need to
    !  access a record that is not 1 to begin with
    !--actually, with wrap-around of records, it means the reading can start
    !  at any point in the file and be OK
    !--intended use: jt_total = 1 here at start of set of runs reading
    !  from the inflow_file, so the first record read will be record 1
    rec = jt_total - 1

    initialized = .true.

end if

rec = rec + 1

if ( recycle ) then
    rec = modulo ( rec - 1, nrec - l_blend + 1 ) + 1
else
    rec = modulo ( rec - 1, nrec ) + 1
end if

read ( unit=lun, rec=rec ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)

$if ($DEBUG)
if ( DEBUG ) write (*, *) sub // ' : read record ', rec
$endif
    
if ( recycle ) then

    if ( rec < l_blend ) then

        recp = nrec - l_blend + 1 + rec

        wgt = 0.5_rprec * ( 1.0_rprec -                              &
                            cos ( pi * real (rec, rprec) / l_blend ) )
            !--wgt = 0+ when rec = 1
            !  wgt = 1- when rec = l_blend

        read ( unit=lun, rec=recp ) u_tmp, v_tmp, w_tmp

        u(iend_w, :, :) = wgt * u(iend_w, :, :) + (1.0_rprec - wgt) * u_tmp
        v(iend_w, :, :) = wgt * v(iend_w, :, :) + (1.0_rprec - wgt) * v_tmp
        w(iend_w, :, :) = wgt * w(iend_w, :, :) + (1.0_rprec - wgt) * w_tmp

    end if

end if

$if ($DEBUG)
if ( DEBUG ) then  !--write out slices as an ascii time series

    if ( .not. init_DEBUG ) then

        inquire ( unit=lun_DEBUG, exist=exst, opened=opn )

        if ( exst .and. (.not. opn) ) then
        
            if ( USE_MPI ) then
                open ( unit=lun_DEBUG, file=debug_file // MPI_suffix )
            else
                open ( unit=lun_DEBUG, file=debug_file )
            end if
        
            write ( lun_DEBUG, '(a)' ) 'variables = "y" "z" "t" "u" "v" "w"'
            write ( lun_DEBUG, '(3(a,i0))' ) 'zone, f=point, i= ', ny,  &
                                             ', j= ', nz,               &
                                             ', k= ', nsteps

        else
            write (*, *) sub // ': problem opening debug file'
            stop
        end if

        init_DEBUG = .true.

    end if

    fmt = '(3(1x,i0),3(1x,es12.5))'
    do jz = 1, nz
        do jy = 1, ny
            write ( lun_DEBUG, fmt ) jy, jz, jt_total, u(iend_w, jy, jz),  &
                                     v(iend_w, jy, jz), w(iend_w, jy, jz)
        end do
    end do

end if
$endif

end subroutine inflow_read

!**********************************************************************
subroutine len_da_file(fname, lenrec, length)
!**********************************************************************
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
character (1) :: cdummy
integer :: lunit, nlo, nhi, mid, kode
logical :: exists, open
!
! find a free unit on which to open the file
!
do lunit = 99, 1, -1
  !--units to skip (compiler dependent)
  select case (lunit)
    case (5:6)
      !--do nothing
    case default
      inquire(unit=lunit, exist=exists, opened=open)
      if(exists .and. .not. open) exit
  end select
end do
open(unit=lunit, file=fname, access="direct", recl=lenrec, iostat=kode)
if(kode /= 0) then
  print *, 'error in len_da_file: ', trim(fname), ' does not exist'
  return
end if
!
! expansion phase
!
mid = 1
do
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode /= 0) exit
  mid = 2 * mid
end do
!
! length is between mid/2 and mid, do binary search to refine
!
nlo = mid/2
nhi = mid
do while(nhi - nlo > 1)
  mid = (nlo + nhi) / 2
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode == 0) then
     nlo = mid
  else
     nhi = mid
  end if
end do
length = nlo
close(unit=lunit)
return
end subroutine len_da_file

!**********************************************************************
subroutine stats_init ()
!***************************************************************
!  This subroutine allocates the memory for arrays
!  used for statistical calculations 

use param, only : L_x,L_y,L_z,dx,dy,dz,nx,ny,nz,nsteps,coord,nproc
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use param, only : spectra_calc, spectra_nloc, spectra_loc
use param, only : tavg_calc
use grid_defs
use functions, only : cell_indx
use messages
use stat_defs, only : point_t, xplane_t, yplane_t, zplane_t
use stat_defs, only : tavg_t, tavg_zplane_t, spectra_t
use stat_defs, only : type_set
implicit none

include 'tecio.h'

!character(120) :: cx,cy,cz
character(120) :: var_list
integer :: fid, i,j,k

logical :: exst

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

if( tavg_calc ) then

  allocate(tavg_t(nx,ny,nz))
  allocate(tavg_zplane_t(nz))
  
  !  Initialize arrays
  !tavg_t(:,:,:) = tavg(0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, 0._rprec)

  do k=1,Nz
    do j=1, Ny
      do i=1, Nx
        call type_set( tavg_t(i,j,k), 0._rprec )
      enddo
    enddo

    call type_set( tavg_zplane_t(k), 0._rprec )

  enddo
    
  !tavg_zplane_t(:) = tavg(0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, &
  !  0._rprec, 0._rprec, 0._rprec, 0._rprec, 0._rprec)

endif

if(spectra_calc) then

  allocate(spectra_t(spectra_nloc))
!  Initialize 
  spectra_t(:) % istart = -1
  spectra_t(:) % ldiff = 0._rprec
  spectra_t(:) % coord = -1

!  Compute istart and ldiff
  do k=1,spectra_nloc

    !  Initialize sub-arrays
    spectra_t(k) % power = 0._rprec

    $if ($MPI)
    if(spectra_loc(k) >= z(1) .and. spectra_loc(k) < z(nz)) then
      spectra_t(k) % coord = coord
      spectra_t(k) % istart = cell_indx('k',dz,spectra_loc(k))
      spectra_t(k) % ldiff = spectra_loc(k) - z(spectra_t(k) % istart)
    endif
    $else
    spectra_t(k) % coord = 0
    spectra_t(k) % istart = cell_indx('k',dz,spectra_loc(k))
    spectra_t(k) % ldiff = spectra_loc(k) - z(spectra_t(k) % istart)
    $endif

  enddo

endif

! Initialize information for x-planar stats/data
if(xplane_calc) then

  allocate(xplane_t(xplane_nloc))

  xplane_t(:) % istart = -1
  xplane_t(:) % ldiff = 0.
  
!  Compute istart and ldiff
  do i=1,xplane_nloc
    xplane_t(i) % istart = cell_indx('i', dx, xplane_loc(i))
    xplane_t(i) % ldiff = xplane_loc(i) - x(xplane_t(i) % istart)   
  enddo
    
endif

! Initialize information for y-planar stats/data
if(yplane_calc) then

  allocate(yplane_t(yplane_nloc))

  yplane_t(:) % istart = -1
  yplane_t(:) % ldiff = 0.
  
!  Compute istart and ldiff
  do j=1,yplane_nloc
    yplane_t(j) % istart = cell_indx('j', dy, yplane_loc(j))
    yplane_t(j) % ldiff = yplane_loc(j) - y(yplane_t(j) % istart)
  enddo
    
endif

! Initialize information for z-planar stats/data
if(zplane_calc) then

  allocate(zplane_t(zplane_nloc))

!  Initialize 
  zplane_t(:) % istart = -1
  zplane_t(:) % ldiff = 0. 
  zplane_t(:) % coord=-1 
  
!  Compute istart and ldiff
  do k=1,zplane_nloc

    $if ($MPI)
    if(zplane_loc(k) >= z(1) .and. zplane_loc(k) < z(nz)) then
      zplane_t(k) % coord = coord
      zplane_t(k) % istart = cell_indx('k',dz,zplane_loc(k))
      zplane_t(k) % ldiff = zplane_loc(k) - z(zplane_t(k) % istart)
    endif
    $else
    zplane_t(k) % coord = 0
    zplane_t(k) % istart = cell_indx('k',dz,zplane_loc(k))
    zplane_t(k) % ldiff = zplane_loc(k) - z(zplane_t(k) % istart)
    $endif

  enddo  
  
endif


!  Open files for instantaneous writing
if(point_calc) then

  allocate(point_t(point_nloc))

  !  Intialize the coord values (-1 shouldn't be used as coord so initialize to this)
  point_t % coord=-1

  do i=1,point_nloc
!  Find the processor in which this point lives
  $if ($MPI)
    if(point_loc(i)%xyz(3) >= z(1) .and. point_loc(i)%xyz(3) < z(nz)) then
      point_t(i) % coord = coord
  
      point_t(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
      point_t(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
      point_t(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))
  
      point_t(i) % xdiff = point_loc(i)%xyz(1) - x(point_t(i) % istart)
      point_t(i) % ydiff = point_loc(i)%xyz(2) - y(point_t(i) % jstart)
      point_t(i) % zdiff = point_loc(i)%xyz(3) - z(point_t(i) % kstart)

      fid=3000*i
  
      !  Can't concatenate an empty string
      point_t(i) % fname=''
      call strcat(point_t(i) % fname,'output/vel.x-')
      call strcat(point_t(i) % fname, point_loc(i)%xyz(1))
      call strcat(point_t(i) % fname,'.y-')
      call strcat(point_t(i) % fname,point_loc(i)%xyz(2))
      call strcat(point_t(i) % fname,'.z-')
      call strcat(point_t(i) % fname,point_loc(i)%xyz(3))
      call strcat(point_t(i) % fname,'.dat')
   
      !  Add tecplot header if file does not exist
      inquire (file=point_t(i) % fname, exist=exst)
      if (.not. exst) then
        var_list = '"t", "u", "v", "w"'
        call write_tecplot_header_xyline(point_t(i) % fname, 'rewind', var_list)
      endif 
  
    endif
  $else
    point_t(i) % coord = 0
    point_t(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
    point_t(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
    point_t(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))
    point_t(i) % xdiff = point_loc(i)%xyz(1) - x(point_t(i) % istart) 
    point_t(i) % ydiff = point_loc(i)%xyz(2) - y(point_t(i) % jstart) 
    point_t(i) % zdiff = point_loc(i)%xyz(3) - z(point_t(i) % kstart)
    !write(cx,'(F9.4)') point_t%xyz(1,i)
    !write(cy,'(F9.4)') point_t%xyz(2,i)
    !write(cz,'(F9.4)') point_t%xyz(3,i)

    fid=3000*i
    point_t(i) % fname=''
    call strcat(point_t(i) % fname,'output/vel.x-')
    call strcat(point_t(i) % fname, point_loc(i)%xyz(1))
    call strcat(point_t(i) % fname,'.y-')
    call strcat(point_t(i) % fname,point_loc(i)%xyz(2))
    call strcat(point_t(i) % fname,'.z-')
    call strcat(point_t(i) % fname,point_loc(i)%xyz(3))
    call strcat(point_t(i) % fname,'.dat')

    !  Add tecplot header if file does not exist
    inquire (file=point_t(i) % fname, exist=exst)
    if (.not. exst) then
      var_list = '"t (s)", "u", "v", "w"'
      call write_tecplot_header_xyline(point_t(i) % fname, 'rewind', var_list)
    endif 

  $endif
  
  enddo
endif

return
end subroutine stats_init

!**********************************************************************
subroutine tavg_init()
!**********************************************************************

!  Load tavg.out files

use param, only : coord, dt, USE_MPI, Nx, Ny, Nz
use messages
use stat_defs, only : tavg_t, tavg_total_time, operator(.MUL.)
use param, only : tavg_nstart, tavg_nend
implicit none

character (*), parameter :: sub_name = mod_name // '.tavg_init'
character (*), parameter :: ftavg_in = 'tavg.out'
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst
integer :: i,j,k

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') ftavg_in, MPI_suffix, coord
$else
fname = trim(adjustl(ftavg_in))
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  !  Nothing to read in
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,*) ' '
    write(*,*)'No previous time averaged data - starting from scratch.'
  endif

  tavg_total_time = 0._rprec
 
  return
endif

$if ($READ_BIG_ENDIAN)
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted', convert='big_endian')
$elseif ($READ_LITTLE_ENDIAN)
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted', convert='little_endian')  
$else
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted')
$endif

read (1) tavg_total_time
read (1) tavg_t

close(1)

! Now initialize all quantities for summation
do k=1, Nz
  do j=1, Ny
    do i=1, Nx
      tavg_t(i,j,k) = tavg_t(i,j,k) .MUL. tavg_total_time
    enddo
  enddo
enddo


!tavg_t % u = tavg_t % u * tavg_total_time
!tavg_t % v = tavg_t % v * tavg_total_time
!tavg_t % w = tavg_t % w * tavg_total_time
!tavg_t % u2 = tavg_t % u2 * tavg_total_time
!tavg_t % v2 = tavg_t % v2 * tavg_total_time
!tavg_t % w2 = tavg_t % w2 * tavg_total_time
!tavg_t % uw = tavg_t % uw * tavg_total_time
!tavg_t % vw = tavg_t % vw * tavg_total_time
!tavg_t % uv = tavg_t % uv * tavg_total_time
!tavg_t % dudz = tavg_t % dudz * tavg_total_time
!tavg_t % dvdz = tavg_t % dvdz * tavg_total_time
!
!tavg_t % txx = tavg_t % txx * tavg_total_time
!tavg_t % txy = tavg_t % txy * tavg_total_time
!tavg_t % tyy = tavg_t % tyy * tavg_total_time
!tavg_t % txz = tavg_t % txz * tavg_total_time
!tavg_t % tyz = tavg_t % tyz * tavg_total_time
!tavg_t % tzz = tavg_t % tzz * tavg_total_time

!tavg_t % fx = tavg_t % fx * tavg_total_time
!tavg_t % fy = tavg_t % fy * tavg_total_time
!tavg_t % fz = tavg_t % fz * tavg_total_time

!tavg_t % cs_opt2 = tavg_t % cs_opt2 * tavg_total_time

return
end subroutine tavg_init

!***************************************************************
subroutine tavg_compute()
!***************************************************************
!  This subroutine collects the stats for each flow 
!  variable quantity
use types, only : rprec
use stat_defs, only : tavg_t, tavg_zplane_t, tavg_total_time
use param, only : nx,ny,nz, dt
use sim_param, only : u,v,w, dudz, dvdz, txx, txy, tyy, txz, tyz, tzz
use immersedbc, only : fx, fy, fz, fxa, fya, fza

implicit none

!use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
integer :: i,j,k
real(rprec) :: u_p, v_p, w_p, dudz_p, dvdz_p

!  Make sure w stuff has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)
call interp_to_uv_grid(dudz, dudz_uv, dudz_uv_tag)
call interp_to_uv_grid(dvdz, dvdz_uv, dvdz_uv_tag)

do k=1,nz  
  do j=1,ny
    do i=1,nx
   
      !  Being cache friendly
      u_p = u(i,j,k)
      v_p = v(i,j,k) 
      w_p = w_uv(i,j,k)
      dudz_p = dudz_uv(i,j,k) 
      dvdz_p = dvdz_uv(i,j,k) 
          
      tavg_t(i,j,k)%u = tavg_t(i,j,k)%u + u_p * dt
      tavg_t(i,j,k)%v = tavg_t(i,j,k)%v + v_p * dt
      tavg_t(i,j,k)%w = tavg_t(i,j,k)%w + w_p * dt
      tavg_t(i,j,k)%u2 = tavg_t(i,j,k)%u2 + u_p * u_p * dt
      tavg_t(i,j,k)%v2 = tavg_t(i,j,k)%v2 + v_p * v_p * dt
      tavg_t(i,j,k)%w2 = tavg_t(i,j,k)%w2 + w_p * w_p * dt
      tavg_t(i,j,k)%uw = tavg_t(i,j,k)%uw+ u_p * w_p * dt
      tavg_t(i,j,k)%vw = tavg_t(i,j,k)%vw + v_p * w_p * dt
      tavg_t(i,j,k)%uv = tavg_t(i,j,k)%uv + u_p * v_p * dt
      tavg_t(i,j,k)%dudz = tavg_t(i,j,k)%dudz + dudz_p * dt
      tavg_t(i,j,k)%dvdz = tavg_t(i,j,k)%dvdz + dvdz_p * dt
      
      tavg_t(i,j,k)%txx = tavg_t(i,j,k)%txx + txx(i,j,k) * dt
      tavg_t(i,j,k)%txy = tavg_t(i,j,k)%txy + txy(i,j,k) * dt
      tavg_t(i,j,k)%tyy = tavg_t(i,j,k)%tyy + tyy(i,j,k) * dt
      tavg_t(i,j,k)%txz = tavg_t(i,j,k)%txz + txz(i,j,k) * dt
      tavg_t(i,j,k)%tyz = tavg_t(i,j,k)%tyz + tyz(i,j,k) * dt
      tavg_t(i,j,k)%tzz = tavg_t(i,j,k)%tzz + tzz(i,j,k) * dt
      
      tavg_t(i,j,k)%fx = tavg_t(i,j,k)%fx + (fx(i,j,k) + fxa(i,j,k)) * dt 
      tavg_t(i,j,k)%fy = tavg_t(i,j,k)%fy + (fy(i,j,k) + fya(i,j,k)) * dt 
      tavg_t(i,j,k)%fz = tavg_t(i,j,k)%fz + (fz(i,j,k) + fza(i,j,k)) * dt
 
      tavg_t(i,j,k)%cs_opt2 = tavg_t(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * dt 
      
    enddo
  enddo
enddo

! Update tavg_total_time for variable time stepping
tavg_total_time = tavg_total_time + dt

return

end subroutine tavg_compute

!**********************************************************************
subroutine tavg_finalize()
!**********************************************************************
use grid_defs, only : x,y,z
use stat_defs, only : tavg_t, tavg_zplane_t, tavg_total_time, tavg
use stat_defs, only : rs, rs_t, rs_zplane_t, cnpy_zplane_t 
use stat_defs, only : operator(.DIV.), operator(.MUL.)
use stat_defs, only :  operator(.ADD.), operator(.SUB.)
use stat_defs, only : type_set, type_zero_bogus
use stat_defs, only : rs_compute, cnpy_tavg_mul
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z, nz_tot

$if($MPI)
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use param, only : MPI_RPREC, rank_of_coord, comm, ierr
use stat_defs, only : rs, tavg
$endif

$if($LVLSET)
use level_set_base, only : phi
$endif

implicit none

include 'tecio.h'

character (*), parameter :: sub_name = mod_name // '.tavg_finalize'
character(64) :: fname_out, fname_vel, fname_vel2, fname_ddz, &
  fname_tau, fname_f, fname_rs, fname_cs
character(64) :: fname_vel_zplane, fname_vel2_zplane, &
  fname_ddz_zplane, fname_tau_zplane, fname_f_zplane, &
  fname_rs_zplane, fname_cs_zplane, fname_cnpy_zplane
  
integer :: i,j,k

$if($MPI)
character(64) :: temp
real(rprec), allocatable, dimension(:) :: z_tot
integer :: MPI_RS, MPI_CNPY, MPI_TSTATS
integer :: rs_type(1), rs_block(1), rs_disp(1)
integer :: cnpy_type(1), cnpy_block(1), cnpy_disp(1)
integer :: tavg_type(1), tavg_block(1), tavg_disp(1)

integer :: ip, kbuf_start, kbuf_end, ktot_start, ktot_end
integer, allocatable, dimension(:) :: gather_coord

type(rs), pointer, dimension(:) :: rs_zplane_tot_t
type(rs), pointer, dimension(:) :: rs_zplane_buf_t

type(rs), pointer, dimension(:) :: cnpy_zplane_tot_t
type(rs), pointer, dimension(:) :: cnpy_zplane_buf_t

type(tavg), pointer, dimension(:) :: tavg_zplane_tot_t
type(tavg), pointer, dimension(:) :: tavg_zplane_buf_t

real(rprec) :: fx_global, fy_global, fz_global
$endif

logical :: opn

type(rs) :: cnpy_avg_t
type(tavg) :: tavg_avg_t
real(rprec), parameter :: favg = real(nx*ny,kind=rprec)

allocate(rs_t(nx,ny,nz), rs_zplane_t(nz))
allocate(cnpy_zplane_t(nz))

$if($MPI)

rs_type = MPI_RPREC
rs_block = 6 ! Number of rs subtypes
rs_disp = 0

cnpy_type = MPI_RPREC
cnpy_block = 6 ! Number of cnpy subtypes
cnpy_disp = 0

tavg_type = MPI_RPREC
tavg_block = 21 ! Number of tavg subtypes
tavg_disp = 0

if(coord == 0) then

  !  Allocate space only on base processor for assembled z-plane data
  ! *_tot_t is the assembled data without the overlap nodes (the final stuff that is outputted)
  ! *_buf_t contains the overlap data and is used to recieve the z-plane data from all other nodes
  allocate(rs_zplane_tot_t(nz_tot), rs_zplane_buf_t(nz*nproc))
  allocate(cnpy_zplane_tot_t(nz_tot), cnpy_zplane_buf_t(nz*nproc))
  allocate(tavg_zplane_tot_t(nz_tot), tavg_zplane_buf_t(nz*nproc))
  
  allocate(z_tot(nz_tot))
  ! In order to ensure that *_tot_t is assembled correctly we make sure that the processor number 
  ! is consistent with the spatial location
  allocate(gather_coord(nproc))
    
  do k=1, nz_tot
  
    z_tot(k) = (dble(k) - 0.5_rprec) * dz
    
  enddo
  
elseif(coord == nproc - 1) then

  !  Zero bogus values
  do j=1, Ny
    do i=1,Nx
     call type_zero_bogus(tavg_t(i,j,nz))
    enddo
  enddo

endif

$else

  !  Zero bogus values
  do j=1, Ny
    do i=1,Nx
     call type_zero_bogus(tavg_t(i,j,nz))
    enddo
  enddo

$endif

! All processors need not do this, but that is ok
!  Set file names
fname_out = 'tavg.out'

fname_vel = 'output/vel_avg.dat'
fname_vel2 = 'output/vel2_avg.dat'
fname_ddz = 'output/ddz_avg.dat'
fname_tau = 'output/tau_avg.dat'
fname_f = 'output/force_avg.dat'
fname_rs = 'output/rs.dat'
fname_cs = 'output/cs_opt2.dat'

fname_vel_zplane = 'output/vel_zplane_avg.dat'
fname_vel2_zplane = 'output/vel2_zplane_avg.dat'
fname_ddz_zplane = 'output/ddz_zplane_avg.dat'
fname_tau_zplane = 'output/tau_zplane_avg.dat'
fname_f_zplane = 'output/force_zplane_avg.dat'
fname_rs_zplane = 'output/rs_zplane.dat'
fname_cnpy_zplane = 'output/cnpy_zplane.dat'
fname_cs_zplane = 'output/cs_opt2_zplane.dat'

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') coord
  fname_out = trim (fname_out) // temp
  
  fname_vel = trim (fname_vel) // temp
  fname_vel2 = trim (fname_vel2) // temp
  fname_ddz = trim (fname_ddz) // temp
  fname_tau = trim (fname_tau) // temp
  fname_f = trim (fname_f) // temp
  fname_rs = trim (fname_rs) // temp
  fname_cs = trim (fname_cs) // temp
   
$endif

!  Perform time averaging operation
!  tavg_t = tavg_t / tavg_total_time
do k=1,Nz
  do j=1, Ny
    do i=1, Nx
      tavg_t(i,j,k) = tavg_t(i,j,k) .DIV. tavg_total_time
    enddo
  enddo
enddo

!tavg_t % u = tavg_t % u / tavg_total_time
!tavg_t % v = tavg_t % v / tavg_total_time
!tavg_t % w = tavg_t % w / tavg_total_time
!tavg_t % u2 = tavg_t % u2 / tavg_total_time
!tavg_t % v2 = tavg_t % v2 / tavg_total_time
!tavg_t % w2 = tavg_t % w2 / tavg_total_time
!tavg_t % uw = tavg_t % uw / tavg_total_time
!tavg_t % vw = tavg_t % vw / tavg_total_time
!tavg_t % uv = tavg_t % uv / tavg_total_time
!tavg_t % dudz = tavg_t % dudz / tavg_total_time
!tavg_t % dvdz = tavg_t % dvdz / tavg_total_time

!tavg_t % txx = tavg_t % txx / tavg_total_time
!tavg_t % txy = tavg_t % txy / tavg_total_time
!tavg_t % tyy = tavg_t % tyy / tavg_total_time
!tavg_t % txz = tavg_t % txz / tavg_total_time
!tavg_t % tyz = tavg_t % tyz / tavg_total_time
!tavg_t % tzz = tavg_t % tzz / tavg_total_time

!tavg_t % fx = tavg_t % fx / tavg_total_time
!tavg_t % fy = tavg_t % fy / tavg_total_time
!tavg_t % fz = tavg_t % fz / tavg_total_time

!tavg_t % cs_opt2 = tavg_t % cs_opt2 / tavg_total_time

!  Average over z-planes
do k=1, nz
  
  !  Initialize to 0 for summations
  call type_set( tavg_zplane_t(k), 0._rprec )

  do j=1, Ny
    do i=1, Nx

      tavg_zplane_t(k) = tavg_zplane_t(k) .ADD. tavg_t(i,j,k)

    enddo
  enddo

  !  Divide by number of summation points 
  tavg_zplane_t(k) = tavg_zplane_t(k) .DIV. favg


  !tavg_zplane_t(k) % u = fa * sum( tavg_t(:,:,k) % u )
  !tavg_zplane_t(k) % v = fa * sum( tavg_t(:,:,k) % v )
  !tavg_zplane_t(k) % w = fa * sum( tavg_t(:,:,k) % w )
  !tavg_zplane_t(k) % u2 = fa * sum( tavg_t(:,:,k) % u2 )
  !tavg_zplane_t(k) % v2 = fa * sum( tavg_t(:,:,k) % v2 )
  !tavg_zplane_t(k) % w2 = fa * sum( tavg_t(:,:,k) % w2 )

  !tavg_zplane_t(k) % uw = fa * sum( tavg_t(:,:,k) % uw )
  !tavg_zplane_t(k) % vw = fa * sum( tavg_t(:,:,k) % vw )
  !tavg_zplane_t(k) % uv = fa * sum( tavg_t(:,:,k) % uv )

  !tavg_zplane_t(k) % dudz = fa * sum( tavg_t(:,:,k) % dudz )
  !tavg_zplane_t(k) % dvdz = fa * sum( tavg_t(:,:,k) % dvdz )
  
  !tavg_zplane_t(k) % txx = fa * sum( tavg_t(:,:,k) % txx )
  !tavg_zplane_t(k) % txy = fa * sum( tavg_t(:,:,k) % txy )
  !tavg_zplane_t(k) % tyy = fa * sum( tavg_t(:,:,k) % tyy )
  !tavg_zplane_t(k) % txz = fa * sum( tavg_t(:,:,k) % txz )
  !tavg_zplane_t(k) % tyz = fa * sum( tavg_t(:,:,k) % tyz )
  !tavg_zplane_t(k) % tzz = fa * sum( tavg_t(:,:,k) % tzz )

  !tavg_zplane_t(k) % fx = fa * sum( tavg_t(:,:,k) % fx )
  !tavg_zplane_t(k) % fy = fa * sum( tavg_t(:,:,k) % fy )
  !tavg_zplane_t(k) % fz = fa * sum( tavg_t(:,:,k) % fz )
  
  !tavg_zplane_t(k) % cs_opt2 = fa * sum( tavg_t(:,:,k) % cs_opt2 )
  
enddo

do k = 1, nz

  !  Initialize to 0
  call type_set( rs_zplane_t(k), 0._rprec)

  do j = 1, ny
    do i = 1, nx
    
      ! Compute the Reynolds stresses: bar(u_i * u_j) - bar(u_i) * bar(u_j)
      rs_t(i,j,k) = rs_compute( tavg_t(i,j,k) )

      !call rs_compute( rs_t(i,j,k), tavg_t(i,j,k) )
      !rs_t(i,j,k) % up2 = tavg_t(i,j,k) % u2 - tavg_t(i,j,k) % u * tavg_t(i,j,k) % u
      !rs_t(i,j,k) % vp2 = tavg_t(i,j,k) % v2 - tavg_t(i,j,k) % v * tavg_t(i,j,k) % v
      !rs_t(i,j,k) % wp2 = tavg_t(i,j,k) % w2 - tavg_t(i,j,k) % w * tavg_t(i,j,k) % w
      !rs_t(i,j,k) % upwp = tavg_t(i,j,k) % uw - tavg_t(i,j,k) % u * tavg_t(i,j,k) % w
      !rs_t(i,j,k) % vpwp = tavg_t(i,j,k) % vw - tavg_t(i,j,k) % v * tavg_t(i,j,k) % w
      !rs_t(i,j,k) % upvp = tavg_t(i,j,k) % uv - tavg_t(i,j,k) % u * tavg_t(i,j,k) % v

      rs_zplane_t(k) = rs_zplane_t(k) .ADD. rs_t(i,j,k) 

    enddo    
  enddo

  rs_zplane_t(k) = rs_zplane_t(k) .DIV. favg
  
  !  Compute the z-plane averaged Reynolds stresses: 
  !rs_zplane_t(k) % up2  = fa * sum( rs_t(:,:,k) % up2 )
  !rs_zplane_t(k) % vp2  = fa * sum( rs_t(:,:,k) % vp2 )
  !rs_zplane_t(k) % wp2  = fa * sum( rs_t(:,:,k) % wp2 )
  !rs_zplane_t(k) % upwp = fa * sum( rs_t(:,:,k) % upwp )
  !rs_zplane_t(k) % vpwp = fa * sum( rs_t(:,:,k) % vpwp )
  !rs_zplane_t(k) % upvp = fa * sum( rs_t(:,:,k) % upvp )
  
enddo


!  Compute the Canopy/Dispersive Stresses: <bar(u_i)*bar(u_j)>_xy - <bar(u_i)>_xy * <bar(u_j)>_xy
do k = 1, nz  

  ! Initialize to 0
  call type_set( cnpy_avg_t, 0._rprec)
  call type_set( tavg_avg_t, 0._rprec)

  do j=1, Ny
    do i=1, Nx

      cnpy_avg_t = cnpy_avg_t .ADD. cnpy_tavg_mul( tavg_t(i,j,k) )
      tavg_avg_t = tavg_avg_t .ADD. tavg_t(i,j,k)
      
    enddo
  enddo

  cnpy_avg_t = cnpy_avg_t .DIV. favg
  tavg_avg_t = tavg_avg_t .DIV. favg

  cnpy_zplane_t(k) = cnpy_avg_t .SUB. cnpy_tavg_mul( tavg_avg_t )

  !cnpy_zplane_t(k) % up2  = fa*sum(tavg_t(:,:,k)%u * tavg_t(:,:,k)%u) - (fa*sum( tavg_t(:,:,k)%u ))**2
  !cnpy_zplane_t(k) % vp2  = fa*sum(tavg_t(:,:,k)%v * tavg_t(:,:,k)%v) - (fa*sum( tavg_t(:,:,k)%v ))**2
  !cnpy_zplane_t(k) % wp2  = fa*sum(tavg_t(:,:,k)%w * tavg_t(:,:,k)%w) - (fa*sum( tavg_t(:,:,k)%w ))**2
  !cnpy_zplane_t(k) % upwp = fa*sum(tavg_t(:,:,k)%u * tavg_t(:,:,k)%w) - fa*sum( tavg_t(:,:,k)%u ) * fa*sum( tavg_t(:,:,k)%w )
  !cnpy_zplane_t(k) % vpwp = fa*sum(tavg_t(:,:,k)%v * tavg_t(:,:,k)%w) - fa*sum( tavg_t(:,:,k)%v ) * fa*sum( tavg_t(:,:,k)%w )
  !cnpy_zplane_t(k) % upvp = fa*sum(tavg_t(:,:,k)%u * tavg_t(:,:,k)%v) - fa*sum( tavg_t(:,:,k)%u ) * fa*sum( tavg_t(:,:,k)%v )
  
enddo


!  Write data to tavg.out
inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname_out, action='write', position='rewind', &
  form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname_out, action='write', position='rewind', &
  form='unformatted', convert='little_endian')
$else
open (1, file=fname_out, action='write', position='rewind', form='unformatted')
$endif

! write the entire structures
write (1) tavg_total_time
write (1) tavg_t          
close(1)

! ----- Write all the 3D data -----
! -- Work around for large data: only write 3 variables at a time. Since things are
! -- written in block format this can be done with mulitple calls to
! -- write_real_data_3D.
$if ($LVLSET)

call write_tecplot_header_ND(fname_vel, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<u>","<v>","<w>"', coord, 2)
!  write phi and x,y,z
call write_real_data_3D(fname_vel, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_vel, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % u, tavg_t % v, tavg_t % w /), 4)

$else

call write_tecplot_header_ND(fname_vel, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<u>","<v>","<w>"', coord, 2)
call write_real_data_3D(fname_vel, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % u, tavg_t % v, tavg_t % w /), &
  4, x, y, z(1:nz))

$endif

$if($LVLSET)

call write_tecplot_header_ND(fname_vel2, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
!  write phi and x,y,z
call write_real_data_3D(fname_vel2, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_vel2, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % u2, tavg_t % v2, tavg_t % w2 /), 4)
call write_real_data_3D(fname_vel2, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % uw, tavg_t % vw, tavg_t % uv /), 4)

$else

call write_tecplot_header_ND(fname_vel2, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
call write_real_data_3D(fname_vel2, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % u2, tavg_t % v2, tavg_t % w2 /), &
  4, x, y, z(1:nz)) 
call write_real_data_3D(fname_vel2, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % uw, tavg_t % vw, tavg_t % uv /), &
  4)

$endif

$if($LVLSET)  
call write_tecplot_header_ND(fname_ddz, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<dudz>","<dvdz>"', coord, 2)
!  write phi and x,y,z
call write_real_data_3D(fname_ddz, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_ddz, 'append', 'formatted', 2, nx, ny, nz, &
  (/ tavg_t % dudz, tavg_t % dvdz /), 4)

$else
call write_tecplot_header_ND(fname_ddz, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<dudz>","<dvdz>"', coord, 2)
call write_real_data_3D(fname_ddz, 'append', 'formatted', 2, nx, ny, nz, &
  (/ tavg_t % dudz, tavg_t % dvdz /), 4, x, y, z(1:nz))

$endif

$if($LVLSET)

call write_tecplot_header_ND(fname_tau, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', coord, 2)  
!  write phi and x,y,z
call write_real_data_3D(fname_tau, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_tau, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % txx, tavg_t % txy, tavg_t % tyy /), 4) 
call write_real_data_3D(fname_tau, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % txz, tavg_t % tyz, tavg_t % tzz /), 4)

$else

call write_tecplot_header_ND(fname_tau, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', coord, 2)
call write_real_data_3D(fname_tau, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % txx, tavg_t % txy, tavg_t % tyy /), &
  4, x, y, z(1:nz))
call write_real_data_3D(fname_tau, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % txz, tavg_t % tyz, tavg_t % tzz /), &
  4)

$endif  

$if($LVLSET)

call write_tecplot_header_ND(fname_f, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', coord, 2)
!  write phi and x,y,z
call write_real_data_3D(fname_f, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_f, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % fx, tavg_t % fy, tavg_t % fz /), 4)  

  $if($MPI)

  call mpi_allreduce(sum(tavg_t(1:nx,1:ny,1:nz-1)%fx), fx_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce(sum(tavg_t(1:nx,1:ny,1:nz-1)%fy), fy_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce(sum(tavg_t(1:nx,1:ny,1:nz-1)%fz), fz_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  
  fx_global = fx_global * dx * dy * dz
  fy_global = fy_global * dx * dy * dz
  fz_global = fz_global * dx * dy * dz

  if(coord == 0) then
    open(unit = 1, file = "output/force_total_avg.dat", status="unknown", position="rewind") 
    write(1,'(a,3e15.6)') '<fx>, <fy>, <fz> : ', fx_global, fy_global, fz_global
    close(1)
  endif

  $endif

    
 
  
$else

call write_tecplot_header_ND(fname_f, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', coord, 2)
call write_real_data_3D(fname_f, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % fx, tavg_t % fy, tavg_t % fz /), &
  4, x, y, z(1:nz))

$endif
  
$if($LVLSET)

call write_tecplot_header_ND(fname_rs, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
!  write phi and x,y,z
call write_real_data_3D(fname_rs, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_rs, 'append', 'formatted', 3, nx, ny, nz, &
  (/ rs_t % up2, rs_t%vp2, rs_t%wp2 /), 4)  
call write_real_data_3D(fname_rs, 'append', 'formatted', 3, nx, ny, nz, &
  (/ rs_t%upwp, rs_t%vpwp, rs_t%upvp /), 4) 

$else

call write_tecplot_header_ND(fname_rs, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)
call write_real_data_3D(fname_rs, 'append', 'formatted', 3, nx, ny, nz, &
  (/ rs_t % up2, rs_t%vp2, rs_t%wp2 /), &
  4, x, y, z(1:nz))
call write_real_data_3D(fname_rs, 'append', 'formatted', 3, nx, ny, nz, &
  (/ rs_t%upwp, rs_t%vpwp, rs_t%upvp /), &
  4)

$endif

$if($LVLSET)

call write_tecplot_header_ND(fname_cs, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<cs2>"', coord, 2)
!  write phi and x,y,z
call write_real_data_3D(fname_cs, 'append', 'formatted', 1, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
call write_real_data_3D(fname_cs, 'append', 'formatted', 1, nx, ny, nz, &
  (/ tavg_t % cs_opt2 /), 4)  

$else

call write_tecplot_header_ND(fname_cs, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<cs2>"', coord, 2)
call write_real_data_3D(fname_cs, 'append', 'formatted', 1, nx, ny, nz, &
  (/ tavg_t % cs_opt2 /), &
  4, x, y, z(1:nz))

$endif

! Construct zplane data 
$if($MPI)

  !  Create MPI type structures consistent with the derived types
  call MPI_TYPE_STRUCT(1, rs_block, rs_disp, rs_type, MPI_RS, ierr)
  Call MPI_Type_commit(MPI_RS,ierr)

  call MPI_TYPE_STRUCT(1, cnpy_block, cnpy_disp, cnpy_type, MPI_CNPY, ierr)
  Call MPI_Type_commit(MPI_CNPY,ierr)  

  call MPI_TYPE_STRUCT(1, tavg_block, tavg_disp, tavg_type, MPI_TSTATS, ierr)
  Call MPI_Type_commit(MPI_TSTATS,ierr)

  call mpi_gather( rs_zplane_t, nz, MPI_RS, rs_zplane_buf_t, nz, &
    MPI_RS, rank_of_coord(0), comm, ierr)
    
  call mpi_gather( cnpy_zplane_t, nz, MPI_CNPY, cnpy_zplane_buf_t, nz, &
    MPI_CNPY, rank_of_coord(0), comm, ierr)    
    
  call mpi_gather( tavg_zplane_t, nz, MPI_TSTATS, tavg_zplane_buf_t, nz, &
    MPI_TSTATS, rank_of_coord(0), comm, ierr)   
  
  call MPI_Type_free (MPI_RS, ierr)
  call MPI_Type_free (MPI_CNPY, ierr)
  call mpi_type_free (MPI_TSTATS, ierr)    
  
  !  Get the rank that was used for mpi_gather (ensure that assembly of {rs,tavg}_zplane_tot_t is
  !  done in increasing coord
  call mpi_gather( coord, 1, MPI_INTEGER, gather_coord, 1, &
  MPI_INTEGER, rank_of_coord(0), comm, ierr)
  
  if(coord == 0) then
  
  !  Need to remove overlapping nodes
  !! Set initial block of data  
  !  rs_zplane_tot_t(1:nz) = rs_zplane_buf_t(1:nz)
  !  cnpy_zplane_tot_t(1:nz) = cnpy_zplane_buf_t(1:nz)
  !  tavg_zplane_tot_t(1:nz) = tavg_zplane_buf_t(1:nz)
    
    do ip=1, nproc
      
      kbuf_start = gather_coord(ip)*nz + 1
      kbuf_end   = kbuf_start + nz - 1
      
      ktot_start = kbuf_start - gather_coord(ip)
      ktot_end = ktot_start + nz - 1
                
      rs_zplane_tot_t(ktot_start:ktot_end) = rs_zplane_buf_t(kbuf_start:kbuf_end)
      cnpy_zplane_tot_t(ktot_start:ktot_end) = cnpy_zplane_buf_t(kbuf_start:kbuf_end)
      tavg_zplane_tot_t(ktot_start:ktot_end) = tavg_zplane_buf_t(kbuf_start:kbuf_end)
      
    enddo   
    
    call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz_tot /), &
      '"z", "<u>","<v>","<w>"', coord, 2)
    call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot_t % u, tavg_zplane_tot_t % v, tavg_zplane_tot_t % w /), 0, z_tot)

    call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz_tot/), &
      '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
    call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot_t % u2, tavg_zplane_tot_t % v2, tavg_zplane_tot_t % w2, &
      tavg_zplane_tot_t % uw, tavg_zplane_tot_t % vw, tavg_zplane_tot_t % uv /), &
      0, z_tot) 
  
    call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz_tot/), &
      '"z", "<dudz>","<dvdz>"', coord, 2)
    call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, Nz_tot, &
      (/ tavg_zplane_tot_t % dudz, tavg_zplane_tot_t % dvdz /), 0, z_tot)
  
    call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', coord, 2)  
    call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot_t % txx, tavg_zplane_tot_t % txy, tavg_zplane_tot_t % tyy, &
      tavg_zplane_tot_t % txz, tavg_zplane_tot_t % tyz, tavg_zplane_tot_t % tzz /), 0, z_tot) 
  
    call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', coord, 2)
    call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot_t % fx, tavg_zplane_tot_t % fy, tavg_zplane_tot_t % fz /), 0, z_tot)  
  
    call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
    call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ rs_zplane_tot_t % up2, rs_zplane_tot_t%vp2, rs_zplane_tot_t%wp2, &
      rs_zplane_tot_t%upwp, rs_zplane_tot_t%vpwp, rs_zplane_tot_t%upvp /), 0, z_tot)    
 
    call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
    call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ cnpy_zplane_tot_t % up2, cnpy_zplane_tot_t%vp2, cnpy_zplane_tot_t%wp2, &
      cnpy_zplane_tot_t%upwp, cnpy_zplane_tot_t%vpwp, cnpy_zplane_tot_t%upvp /), 0, z_tot)         
      
    call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz_tot/), &
      '"z", "<cs2>"', coord, 2)
    call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, Nz_tot, &
      (/ tavg_zplane_tot_t % cs_opt2 /), 0, z_tot)        
      
    deallocate(tavg_zplane_tot_t, tavg_zplane_buf_t)
    deallocate(rs_zplane_tot_t, rs_zplane_buf_t)
    deallocate(cnpy_zplane_tot_t, cnpy_zplane_buf_t)
    deallocate(z_tot, gather_coord)
  
  endif

$else

call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz /), &
   '"z", "<u>","<v>","<w>"', coord, 2)  
call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane_t % u, tavg_zplane_t % v, tavg_zplane_t % w /), 0, z(1:nz))

call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz/), &
   '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane_t % u2, tavg_zplane_t % v2, tavg_zplane_t % w2, &
  tavg_zplane_t % uw, tavg_zplane_t % vw, tavg_zplane_t % uv /), 0, z(1:nz)) 
  
call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz/), &
   '"z", "<dudz>","<dvdz>"', coord, 2)
call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, nz, &
  (/ tavg_zplane_t % dudz, tavg_zplane_t % dvdz /), 0, z(1:nz))
  
call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', coord, 2)  
call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane_t % txx, tavg_zplane_t % txy, tavg_zplane_t % tyy, &
  tavg_zplane_t % txz, tavg_zplane_t % tyz, tavg_zplane_t % tzz /), 0, z(1:nz)) 
  
call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', coord, 2)
call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane_t % fx, tavg_zplane_t % fy, tavg_zplane_t % fz /), 0, z(1:nz))  
  
call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, nz, &
  (/ rs_zplane_t % up2, rs_zplane_t%vp2, rs_zplane_t%wp2, &
  rs_zplane_t%upwp, rs_zplane_t%vpwp, rs_zplane_t%upvp /), 0, z(1:nz))

call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, nz, &
  (/ cnpy_zplane_t % up2, cnpy_zplane_t%vp2, cnpy_zplane_t%wp2, &
  cnpy_zplane_t%upwp, cnpy_zplane_t%vpwp, cnpy_zplane_t%upvp /), 0, z(1:nz))  
  
call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz/), &
   '"z", "<cs2>"', coord, 2)
call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, nz, &
  (/ tavg_zplane_t % cs_opt2 /), 0, z(1:nz))    
  
$endif

$if ($MPI)
!  Sync data across all nodes for a subset of varibles which contain bogus initialization
!call mpi_sync_real_array(  tavg_t % txx, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array(  tavg_t % txy, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array(  tavg_t % tyy, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array(  tavg_t % txz, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array(  tavg_t % tyz, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array(  tavg_t % tzz, MPI_SYNC_DOWNUP )
$endif 

deallocate(tavg_t, tavg_zplane_t, rs_t, rs_zplane_t, cnpy_zplane_t)

return
end subroutine tavg_finalize

!**********************************************************************
subroutine spectra_init()
!**********************************************************************
use types, only : rprec
use messages
use param, only : coord, dt, USE_MPI, spectra_nloc, lh, nx
use stat_defs, only : spectra_t, spectra_total_time
implicit none

character (*), parameter :: sub_name = mod_name // '.spectra_init'
character (*), parameter :: fspectra_in = 'spectra.out'
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst
integer :: k

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') fspectra_in, MPI_suffix, coord
$else
fname = trim(adjustl(fspectra_in))
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  !  Nothing to read in
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,*) ' '
    write(*,*)'No previous spectra data - starting from scratch.'
  endif

  spectra_total_time = 0._rprec

  return
endif

$if ($READ_BIG_ENDIAN)
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted', convert='big_endian')
$elseif ($READ_LITTLE_ENDIAN)
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted', convert='little_endian')
$else
open (1, file=fname, action='read', position='rewind',  &
  form='unformatted')
$endif

read (1) spectra_total_time
do k=1, spectra_nloc
  read (1) spectra_t(k) % power
enddo

close(1)

!  Prepare for averaging
do k=1, spectra_nloc
  spectra_t(k) % power = spectra_t(k) % power * spectra_total_time
enddo


return
end subroutine spectra_init

!***************************************************************
subroutine spectra_compute()
!***************************************************************
use types, only : rprec
use sim_param, only : u
use param, only : Nx, Ny, dt, dz, lh
use param, only : spectra_nloc
use stat_defs, only : spectra_t, spectra_total_time
use functions, only : linear_interp
use fft, only : forw_spectra
implicit none

integer :: i, j, k
real(rprec), allocatable, dimension(:) :: ui, uhat, power

! Interpolation variable
allocate(ui(nx), uhat(nx), power(lh))

!  Loop over all spectra locations
do k=1, spectra_nloc
  
  $if ($MPI)
  if(spectra_t(k) % coord == coord) then
  $endif

  do j=1,Ny

    do i=1,Nx
        ! Interpolate to the specified plane
        ui(i) = linear_interp(u(i,j,spectra_t(k) % istart),u(i,j,spectra_t(k) % istart+1), &
                              dz, spectra_t(k) % ldiff)

    enddo

    ! 1) Compute uhat for the given j
    ui = ui - sum(ui) / Nx ! Remove the mean
    ! Compute FFT
    call rfftw_f77_one(forw_spectra, ui, uhat)
    !  Normalize
    uhat = uhat / Nx

    ! 2) Compute power spectra for given j
    power(1) = uhat(1)**2
    do i=2,lh-1
      power(i) = 2._rprec*(uhat(i)**2 + uhat(Nx-i+2)**2)
    enddo
    power(lh) = uhat(lh)**2 ! Nyquist

    ! Sum jth component 
    spectra_t(k) % power = spectra_t(k) % power + power

  enddo

  !  Update spectra total time
  spectra_total_time = spectra_total_time + dt

  $if ($MPI)
  endif
  $endif

enddo  
  
deallocate(ui, uhat, power)
return
end subroutine spectra_compute

!**********************************************************************
subroutine spectra_finalize()
!**********************************************************************
use types, only : rprec
use param, only : lh, spectra_nloc, spectra_loc
use messages
use fft, only : kx
use stat_defs, only : spectra_t, spectra_total_time
implicit none

include 'tecio.h'

character (*), parameter :: sub_name = mod_name // '.spectra_finalize'
character(25) :: cl
character (64) :: fname
character(64) :: fname_out

$if($MPI)
character(64) :: temp
$endif

integer :: i, k

logical :: opn

!  Set file names
fname_out = 'spectra.out'

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') coord
  fname_out = trim (fname_out) // temp  
$endif

!  Loop over all zplane locations
do k=1,spectra_nloc
  
  $if ($MPI)
  if(spectra_t(k) % coord == coord) then
  $endif

  ! Finalize averaging for power spectra (averaging over y and time)
  spectra_t(k) % power = (spectra_t(k) % power / Ny) / spectra_total_time

  !  Create unique file name
  write(cl,'(F9.4)') spectra_loc(k)
  !  Convert total iteration time to string
  write(fname,*) 'output/spectra.z-',trim(adjustl(cl)),'.dat'
  fname=trim(adjustl(fname))

  !  Omitting Nyquist from output
  call write_tecplot_header_ND(fname, 'rewind', 2, (/ lh-1/), &
    '"k", "E(k)"', k, 2 ) 
  call write_real_data_1D(fname, 'append', 'formatted', 1, lh-1, &
    (/ spectra_t(k) % power(1:lh-1) /), 0, (/ kx(1:lh-1,1) /))
    
  $if ($MPI)
  endif
  $endif

enddo  

!  Write data to tavg.out
inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname_out, action='write', position='rewind', &
  form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname_out, action='write', position='rewind', &
  form='unformatted', convert='little_endian')
$else
open (1, file=fname_out, action='write', position='rewind', form='unformatted')
$endif

! write the entire structures
write (1) spectra_total_time
do k=1, spectra_nloc
  write (1) spectra_t(k) % power
enddo
close(1)
  
deallocate(spectra_t)

return

end subroutine spectra_finalize

end module io
