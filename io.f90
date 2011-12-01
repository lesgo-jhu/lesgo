!///////////////////////////////////////////////////////////////////////////////
module io
!///////////////////////////////////////////////////////////////////////////////
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, path,  &
                  coord, rank, nproc, jt_total, total_time, &
                  total_time_dim, lbz, jzmin, jzmax
use param, only : cumulative_time, fcumulative_time
use sim_param, only : w, dudz, dvdz
use sgs_param,only:Cs_opt2
use string_util
use messages

implicit none

save
private

!!$public openfiles,output_loop,output_final,                   &
!!$     inflow_write, avg_stats
public jt_total, openfiles, output_loop, output_final

public stats_init

character (*), parameter :: mod_name = 'io'

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine openfiles()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : use_cfl_dt, dt, cfl_f
use sim_param,only:path

implicit none

logical :: exst

! Temporary values used to read time step and CFL from file
real(rprec) :: dt_r, cfl_r

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    
    read(1, *) jt_total, total_time, total_time_dim, dt_r, cfl_r
    
    close (1)
  else  !--assume this is the first run on cumulative time
    if( coord == 0 ) then
      write (*, *) '--> Assuming jt_total = 0, total_time = 0., total_time_dim = 0.'
    endif
    jt_total = 0
    total_time = 0.
    total_time_dim = 0._rprec
  end if

end if

! Update dynamic time stepping info if required; otherwise discard.
if( use_cfl_dt ) then
   dt = dt_r
   cfl_f = cfl_r
endif

return
end subroutine openfiles

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_loop(jt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      if (coord == 0) then
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
      if (coord == 0) then
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
      if (coord == 0) then   
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
      if (coord == 0) then        
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
      if (coord == 0) then        
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
      if (coord == 0) then        
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
      if (coord == 0) then        
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine inst_write(itype)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine is used to write all of the instantaneous data from
! lesgo to file. The types of data written are:
! 
!   points   : itype=1
!   domain   : itype=2
!   x-planes : itype=3
!   y-planes : itype=4
!   z-planes : itype=5
!
! For the points and planar data, this subroutine writes using the
! locations specfied from the param module. All of the data is written
! using the tecryte library. If additional instantenous values are
! desired to be written, they should be done so using this subroutine.
!
! REMARK: It may be desired to convert this subroutine to a driver and
! have it call appropriate subroutines to perform the writing, just to
! clean things up a bit
!
use functions, only : linear_interp, trilinear_interp, interp_to_uv_grid
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use grid_defs, only : grid_t
use sim_param, only : u,v,w,dudx,dvdy,dwdz
$if($DEBUG)
use sim_param, only : p, dpdx, dpdy, dpdz
use sim_param, only : RHSx, RHSy, RHSz
$endif
use stat_defs, only : xplane_t, yplane_t, zplane_t, point_t
$if($MPI)
use mpi
use param, only : ld, ny, nz, MPI_RPREC, down, up, comm, status, ierr
$endif

$if($LVLSET)
use level_set_base, only : phi
use sim_param, only : fx, fy, fz, fxa, fya, fza
$endif
use param, only : dx,dy,dz
use param, only : sgs_model
use sgs_param, only : F_LM,F_MM,F_QN,F_NN,beta,Cs_opt2,Nu_t
implicit none

include 'tecryte.h'      

integer, intent(IN) :: itype

character (*), parameter :: sub_name = mod_name // '.inst_write'

character(25) :: cl, ct
character (64) :: fname
$if($MPI)
character(64) :: temp
$endif
character(256) :: var_list
integer :: n, i, j, k, nvars

real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi
real(rprec), allocatable, dimension(:,:,:) :: w_uv

$if($LVLSET)
real(rprec), allocatable, dimension(:,:,:) :: fx_tot, fy_tot, fz_tot
$endif

$if($DEBUG)
real(rprec), allocatable, dimension(:,:,:) :: divvel
$endif

real(rprec), pointer, dimension(:) :: x,y,z,zw

$if($OUTPUT_EXTRA)
! Arrays used for outputing slices of LDSM variables
real(rprec), allocatable, dimension(:,:) :: F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s
real(rprec), allocatable, dimension(:,:,:) :: F_LM_uv,F_MM_uv,F_QN_uv,F_NN_uv,beta_uv,Cs_opt2_uv,Nu_t_uv
$endif

! Nullify pointers
nullify(x,y,z,zw)

! Set grid pointers
x => grid_t % x
y => grid_t % y
z => grid_t % z
zw => grid_t % zw

!  Allocate space for the interpolated w values
allocate(w_uv(nx,ny,lbz:nz))

!  Make sure w has been interpolated to uv-grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

$if($OUTPUT_EXTRA)
!  Allocate arrays and interpolate to uv grid for LDSM output
if( sgs_model == 4 .or. sgs_model == 5 ) then

  if( itype == 3 .or. itype == 4 .or. itype == 5 ) then

    allocate( F_LM_uv(nx,ny,nz), F_MM_uv(nx,ny,nz) )
    allocate( beta_uv(nx,ny,nz), Cs_opt2_uv(nx,ny,nz) )
    allocate( Nu_t_uv(nx,ny,nz) )

    F_LM_uv = interp_to_uv_grid( F_LM(1:nx,1:ny,1:nz), 1 )
    F_MM_uv = interp_to_uv_grid( F_MM(1:nx,1:ny,1:nz), 1 )
    beta_uv = interp_to_uv_grid( beta(1:nx,1:ny,1:nz), 1 )
    Cs_opt2_uv = interp_to_uv_grid( Cs_opt2(1:nx,1:ny,1:nz), 1 )
    Nu_t_uv = interp_to_uv_grid( Nu_t(1:nx,1:ny,1:nz), 1 )

    if( sgs_model == 5) then

      allocate( F_QN_uv(nx, ny, nz), F_NN_uv(nx,ny,nz) )

      F_QN_uv = interp_to_uv_grid( F_QN(1:nx,1:ny,1:nz), 1 )
      F_NN_uv = interp_to_uv_grid( F_NN(1:nx,1:ny,1:nz), 1 )

    endif

  endif

endif       
$endif

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

  !////////////////////////////////////////////
  !/// WRITE VELOCITY                       ///
  !////////////////////////////////////////////

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
  var_list = '"x", "y", "z", "u", "v", "w", "phi"'
  nvars = 7
  $else
  var_list = '"x", "y", "z", "u", "v", "w"'
  nvars = 6
  $endif
  
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))

  $if($LVLSET)
  call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
    (/ u(1:nx,1:ny,1:nz), &
    v(1:nx,1:ny,1:nz), &
    w_uv(1:nx,1:ny,1:nz), &
    phi(1:nx,1:ny,1:nz)/), & 
    4, x, y, z(1:nz))
  $else
  call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
    (/ u(1:nx,1:ny,1:nz), &
    v(1:nx,1:ny,1:nz), &
    w_uv(1:nx,1:ny,1:nz) /), &
    4, x, y, z(1:nz))
  $endif

  $if($MPI)
    ! Ensure that all processes finish before attempting to write 
    ! additional files. Otherwise it may flood the system with 
    ! too many I/O requests and crash the process 
    call mpi_barrier( comm, ierr )
  $endif

  !  Output instantaneous force field 
  $if($LVLSET)
    !////////////////////////////////////////////
    !/// WRITE FORCES                         ///
    !////////////////////////////////////////////

    ! Compute the total forces 
    call force_tot()

    !  Open file which to write global data
    write (fname,*) 'output/force.', trim(adjustl(ct)),'.dat'
    fname = trim(adjustl(fname))

    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    var_list = '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>", "phi"'
    nvars = 7

    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
    call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
      (/ fx_tot, &
      fy_tot, &
      fz_tot, &
      phi(1:nx,1:ny,1:nz) /), &
      4, x, y, z(1:nz))

    deallocate(fx_tot, fy_tot, fz_tot)

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

  $endif

  $if($DEBUG)

    !////////////////////////////////////////////
    !/// WRITE VELOCITY DIVERGENCE            ///
    !////////////////////////////////////////////


    !  Output divergence of velocity field
    allocate(divvel(nx,ny,nz))
    divvel = dudx(1:nx,1:ny,1:nz) + dvdy(1:nx,1:ny,1:nz) + dwdz(1:nx,1:ny,1:nz)

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
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 2, nx, ny,nz, &
           (/ divvel, phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
    $else
      var_list = '"x", "y", "z", "divvel"'
      nvars = 4
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny,nz, &
           (/ divvel /), 4, x, y, z(1:nz))
    $endif

    deallocate(divvel)

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

    !////////////////////////////////////////////
    !/// WRITE PRESSURE                       ///
    !////////////////////////////////////////////

    !  Open file which to write global data
    write (fname,*) 'output/pressure.', trim(adjustl(ct)),'.dat'
    fname = trim(adjustl(fname))
  
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    call pressure_sync()

    $if($LVLSET)
      var_list = '"x", "y", "z", "p", "dpdx", "dpdy", "dpdz", "phi"'
      nvars = 8
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 5, nx, ny,nz, &
           (/ p(1:nx,1:ny,1:nz), &
           dpdx(1:nx,1:ny,1:nz), &
           dpdy(1:nx,1:ny,1:nz), &
           interp_to_uv_grid(dpdz(1:nx,1:ny,1:nz),1), &
           phi(1:nx,1:ny,1:nz) /), &
           4, x, y, z(1:nz))
    $else
      var_list = '"x", "y", "z", "p", "dpdx", "dpdy", "dpdz"'
      nvars = 7
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
        (/ p(1:nx,1:ny,1:nz), &
        dpdx(1:nx,1:ny,1:nz), &
        dpdy(1:nx,1:ny,1:nz), &
        interp_to_uv_grid(dpdz(1:nx,1:ny,1:nz),1) /), &
        4, x, y, z(1:nz))
    $endif

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif
  
    !////////////////////////////////////////////
    !/// WRITE RHS                            ///
    !////////////////////////////////////////////

    !  Open file which to write global data
    write (fname,*) 'output/RHS.', trim(adjustl(ct)),'.dat'
    fname = trim(adjustl(fname))

    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif  
  
    call RHS_sync()

    $if($LVLSET)
      var_list = '"x", "y", "z", "RHSx", "RHSy", "RHSz", "phi"'
      nvars = 7
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
           (/ RHSx(1:nx,1:ny,1:nz), &
           RHSy(1:nx,1:ny,1:nz), &
           interp_to_uv_grid(RHSz(1:nx,1:ny,1:nz),1), &
           phi(1:nx,1:ny,1:nz) /), & 
           4, x, y, z(1:nz))
    $else
      var_list = '"x", "y", "z", "RHSx", "RHSy", "RHSz"'
      nvars = 6
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 3, nx, ny,nz, &
           (/ RHSx(1:nx,1:ny,1:nz), &
           RHSy(1:nx,1:ny,1:nz), &
           interp_to_uv_grid(RHSz(1:nx,1:ny,1:nz),1) /), &
           4, x, y, z(1:nz))
    $endif

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

  $endif

!  Write instantaneous x-plane values
elseif(itype==3) then

  allocate(ui(1,ny,nz), vi(1,ny,nz), wi(1,ny,nz))

  $if($LVLSET)
  call force_tot()
  $endif

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
      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4))  
  
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
      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', &
      numtostr(coord,6), 2, real(total_time,4))

    !  Sum both induced forces, f{x,y,z}, and applied forces, f{x,y,z}a
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(fx_tot(xplane_t(i) % istart,j,k), &
             fx_tot(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)

        vi(1,j,k) = linear_interp(fy_tot(xplane_t(i) % istart,j,k), &
             fy_tot(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)

        wi(1,j,k) = linear_interp(fz_tot(xplane_t(i) % istart,j,k), &
             fz_tot(xplane_t(i) % istart+1,j,k), dx, xplane_t(i) % ldiff)

      enddo
    enddo


    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))


    $endif

    $if($OUTPUT_EXTRA)

    !////////////////////////////////////////////
    !/// WRITE LDSM                           ///
    !////////////////////////////////////////////

    if( sgs_model == 4 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do j=1,Ny
!
          F_LM_s(j,k) = linear_interp(F_LM_uv(xplane_t(i) % istart,j,k), &
               F_LM_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          F_MM_s(j,k) = linear_interp(F_MM_uv(xplane_t(i) % istart,j,k), &
               F_MM_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          beta_s(j,k) = linear_interp(beta_uv(xplane_t(i) % istart,j,k), &
               beta_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          Cs_opt2_s(j,k) = linear_interp(Cs_opt2_uv(xplane_t(i) % istart,j,k), &
               Cs_opt2_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          Nu_t_s(j,k) = linear_interp(Nu_t_uv(xplane_t(i) % istart,j,k), &
               Nu_t_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)

        enddo
      enddo

      write(fname,*) 'output/ldsm.x-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
      fname=trim(adjustl(fname))

      $if ($MPI)
      !  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
      $endif      

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 8, (/ 1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 5, 1,ny,nz, &
           (/ F_LM_s, &
           F_MM_s, &
           beta_s, &
           Cs_opt2_s, &
           Nu_t_s /), &
           2, (/ xplane_loc(i) /), y, z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

    elseif( sgs_model == 5 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(F_QN_s(nx,nz),F_NN_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do j=1,Ny
!
          F_LM_s(j,k) = linear_interp(F_LM_uv(xplane_t(i) % istart,j,k), &
               F_LM_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          F_MM_s(j,k) = linear_interp(F_MM_uv(xplane_t(i) % istart,j,k), &
               F_MM_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          F_QN_s(j,k) = linear_interp(F_QN_uv(xplane_t(i) % istart,j,k), &
               F_QN_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)  
          F_NN_s(j,k) = linear_interp(F_NN_uv(xplane_t(i) % istart,j,k), &
               F_NN_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)                                         
          beta_s(j,k) = linear_interp(beta_uv(xplane_t(i) % istart,j,k), &
               beta_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          Cs_opt2_s(j,k) = linear_interp(Cs_opt2_uv(xplane_t(i) % istart,j,k), &
               Cs_opt2_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)
          Nu_t_s(j,k) = linear_interp(Nu_t_uv(xplane_t(i) % istart,j,k), &
               Nu_t_uv(xplane_t(i) % istart+1,j,k), &
               dx, xplane_t(i) % ldiff)

        enddo
      enddo

      write(fname,*) 'output/ldsm.x-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
      fname=trim(adjustl(fname))

      $if ($MPI)
      !  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
      $endif      

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 10, (/ 1, Ny+1, Nz/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 7, 1,ny,nz, &
        (/ F_LM_s, &
        F_MM_s, &
        F_QN_s, &
        F_NN_s, &
        beta_s, &
        Cs_opt2_s, &
        Nu_t_s /), &
        2, (/ xplane_loc(i) /), y, z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s)

    endif   

    $endif    
    
  enddo   
  
  deallocate(ui,vi,wi)

  $if($LVLSET)
  deallocate( fx_tot, fy_tot, fz_tot )
  $endif

!  Write instantaneous y-plane values
elseif(itype==4) then
  
  allocate(ui(nx,1,nz), vi(nx,1,nz), wi(nx,1,nz))

  $if($LVLSET)
  call force_tot()
  $endif
  
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
      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4)) 
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
      '"x", "y", "z", "fx", "fy", "fz"', numtostr(coord,6), 2, real(total_time,4))  
  
    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(fx_tot(i,yplane_t(j) % istart,k), &
             fx_tot(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)
        vi(i,1,k) = linear_interp(fy_tot(i,yplane_t(j) % istart,k), &
             fy_tot(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)
        wi(i,1,k) = linear_interp(fz_tot(i,yplane_t(j) % istart,k), &
             fz_tot(i,yplane_t(j) % istart+1,k), dy, yplane_t(j) % ldiff)

      enddo
    enddo
    
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))       
    
    $endif

    $if($OUTPUT_EXTRA)

    !////////////////////////////////////////////
    !/// WRITE LDSM                           ///
    !////////////////////////////////////////////

    if( sgs_model == 4 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do i=1,Nx
!
          F_LM_s(i,k) = linear_interp(F_LM_uv(i,yplane_t(j) % istart,k), &
               F_LM_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          F_MM_s(i,k) = linear_interp(F_MM_uv(i,yplane_t(j) % istart,k), &
               F_MM_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          beta_s(i,k) = linear_interp(beta_uv(i,yplane_t(j) % istart,k), &
               beta_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          Cs_opt2_s(i,k) = linear_interp(Cs_opt2_uv(i,yplane_t(j) % istart,k), &
               Cs_opt2_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          Nu_t_s(i,k) = linear_interp(Nu_t_uv(i,yplane_t(j) % istart,k), &
               Nu_t_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)

        enddo
      enddo

      write(fname,*) 'output/ldsm.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
      fname=trim(adjustl(fname))

      $if ($MPI)
      !  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
      $endif      

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 8, (/ Nx+1, 1, Nz/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 5, nx,1,nz, &
        (/ F_LM_s, &
        F_MM_s, &
        beta_s, &
        Cs_opt2_s, &
        Nu_t_s /), &
        1, x, (/ yplane_loc(j) /), z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

    elseif( sgs_model == 5 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(F_QN_s(nx,nz),F_NN_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do i=1,Nx

          F_LM_s(i,k) = linear_interp(F_LM_uv(i,yplane_t(j) % istart,k), &
               F_LM_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          F_MM_s(i,k) = linear_interp(F_MM_uv(i,yplane_t(j) % istart,k), &
               F_MM_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          F_QN_s(i,k) = linear_interp(F_QN_uv(i,yplane_t(j) % istart,k), &
               F_QN_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          F_NN_s(i,k) = linear_interp(F_NN_uv(i,yplane_t(j) % istart,k), &
               F_NN_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)                                      
          beta_s(i,k) = linear_interp(beta_uv(i,yplane_t(j) % istart,k), &
               beta_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          Cs_opt2_s(i,k) = linear_interp(Cs_opt2_uv(i,yplane_t(j) % istart,k), &
               Cs_opt2_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)
          Nu_t_s(i,k) = linear_interp(Nu_t_uv(i,yplane_t(j) % istart,k), &
               Nu_t_uv(i,yplane_t(j) % istart+1,k), &
               dy, yplane_t(j) % ldiff)

        enddo
      enddo

      write(fname,*) 'output/ldsm.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
      fname=trim(adjustl(fname))

      $if ($MPI)
      !  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
      $endif      

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ',  "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 10, (/ Nx+1, 1, Nz/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 7, nx,1,nz, &
        (/ F_LM_s, &
        F_MM_s, &
        F_QN_s, &
        F_NN_s, &
        beta_s, &
        Cs_opt2_s, &
        Nu_t_s /), &
        1, x, (/ yplane_loc(j) /), z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s)

    endif   

    $endif

  enddo  

  deallocate(ui,vi,wi)

  $if($LVLSET)
  deallocate(fx_tot, fy_tot, fz_tot)
  $endif
  
!  Write instantaneous z-plane values
elseif(itype==5) then

  allocate(ui(nx,ny,1), vi(nx,ny,1), wi(nx,ny,1))

  $if($LVLSET)
  call force_tot()
  $endif

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
      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4)) 

    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(u(i,j,zplane_t(k) % istart), &
             u(i,j,zplane_t(k) % istart+1), &
             dz, zplane_t(k) % ldiff)
        vi(i,j,1) = linear_interp(v(i,j,zplane_t(k) % istart), &
             v(i,j,zplane_t(k) % istart+1), &
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
      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', &
      numtostr(coord,6), 2, real(total_time,4))

    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(fx_tot(i,j,zplane_t(k) % istart), &
             fx_tot(i,j,zplane_t(k) % istart+1), &
             dz, zplane_t(k) % ldiff) 
        vi(i,j,1) = linear_interp(fy_tot(i,j,zplane_t(k) % istart), &
             fy_tot(i,j,zplane_t(k) % istart+1), &
             dz, zplane_t(k) % ldiff) 
        wi(i,j,1) = linear_interp(fz_tot(i,j,zplane_t(k) % istart), &
             fz_tot(i,j,zplane_t(k) % istart+1), &
             dz, zplane_t(k) % ldiff)

      enddo
    enddo
 
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
         (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /) )      
    
    $endif

    $if($OUTPUT_EXTRA)

    !////////////////////////////////////////////
    !/// WRITE LDSM                           ///
    !////////////////////////////////////////////

    if( sgs_model == 4 ) then

      allocate(F_LM_s(nx,ny),F_MM_s(nx,ny))
      allocate(beta_s(nx,ny),Cs_opt2_s(nx,ny))
      allocate(Nu_t_s(nx,ny))

      do j=1,Ny
        do i=1,Nx
!
          F_LM_s(i,j) = linear_interp(F_LM_uv(i,j,zplane_t(k) % istart), &
               F_LM_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)
          F_MM_s(i,j) = linear_interp(F_MM_uv(i,j,zplane_t(k) % istart), &
               F_MM_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)
          beta_s(i,j) = linear_interp(beta_uv(i,j,zplane_t(k) % istart), &
               beta_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)
          Cs_opt2_s(i,j) = linear_interp(Cs_opt2_uv(i,j,zplane_t(k) % istart), &
               Cs_opt2_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)
          Nu_t_s(i,j) = linear_interp(Nu_t_uv(i,j,zplane_t(k) % istart), &
               Nu_t_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)

        enddo
      enddo

      write(fname,*) 'output/ldsm.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
      fname=trim(adjustl(fname))

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 8, (/ Nx+1, Ny+1, 1/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))

      call write_real_data_3D(fname, 'append', 'formatted', 5, nx,ny,1, &
        (/ F_LM_s, &
        F_MM_s, &
        beta_s, &
        Cs_opt2_s, &
        Nu_t_s /), &
        4, x, y, (/ zplane_loc(k) /) )

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

     elseif( sgs_model == 5 ) then
 
      allocate(F_LM_s(nx,ny),F_MM_s(nx,ny))
      allocate(F_QN_s(nx,ny),F_NN_s(nx,ny))
      allocate(beta_s(nx,ny),Cs_opt2_s(nx,ny))
      allocate(Nu_t_s(nx,ny))

      do j=1,Ny
        do i=1,Nx
!
          F_LM_s(i,j) = linear_interp(F_LM_uv(i,j,zplane_t(k) % istart), &
               F_LM_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff) 
          F_MM_s(i,j) = linear_interp(F_MM_uv(i,j,zplane_t(k) % istart), &
               F_MM_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff) 
          F_QN_s(i,j) = linear_interp(F_QN_uv(i,j,zplane_t(k) % istart), &
               F_QN_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)
          F_NN_s(i,j) = linear_interp(F_NN_uv(i,j,zplane_t(k) % istart), &
               F_NN_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)
          beta_s(i,j) = linear_interp(beta_uv(i,j,zplane_t(k) % istart), &
               beta_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)            
          Cs_opt2_s(i,j) = linear_interp(Cs_opt2_uv(i,j,zplane_t(k) % istart), &
               Cs_opt2_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)                                         
          Nu_t_s(i,j) = linear_interp(Nu_t_uv(i,j,zplane_t(k) % istart), &
               Nu_t_uv(i,j,zplane_t(k) % istart+1), &
               dz, zplane_t(k) % ldiff)

        enddo
      enddo      

      write(fname,*) 'output/ldsm.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
      fname=trim(adjustl(fname))

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 10, (/ Nx+1, Ny+1, 1/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))

      call write_real_data_3D(fname, 'append', 'formatted', 7, nx,ny,1, &
        (/ F_LM_s, &
        F_MM_s, &
        F_QN_s, &
        F_NN_s, &
        beta_s, &
        Cs_opt2_s, &
        Nu_t_s /), &
        4, x, y, (/ zplane_loc(k) /) )         

      deallocate(F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s)

    endif
    $endif

    $if ($MPI)
    endif
    $endif

  enddo  
  
  deallocate(ui,vi,wi)
  
  $if($LVLSET)
  deallocate(fx_tot, fy_tot, fz_tot)
  $endif

else
  write(*,*) 'Error: itype not specified properly to inst_write!'
  stop
endif

deallocate(w_uv)
nullify(x,y,z,zw)

$if($LVLSET or $DEBUG)
contains
$endif

$if($LVLSET)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine force_tot()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$if($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
$endif
implicit none

! Zero bogus values
fx(:,:,nz) = 0._rprec
fy(:,:,nz) = 0._rprec
fz(:,:,nz) = 0._rprec

!  Sum both the induced and applied forces
allocate(fx_tot(nx,ny,nz), fy_tot(nx,ny,nz), fz_tot(nx,ny,nz))
fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)

$if($MPI)
!  Sync forces
! call mpi_sendrecv (fx_tot(:,:,1), nx*ny, MPI_RPREC, down, 1,  &
!                    fx_tot(:,:,nz), nx*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (fy_tot(:,:,1), nx*ny, MPI_RPREC, down, 1,  &
!                    fy_tot(:,:,nz), nx*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (fz_tot(:,:,1), nx*ny, MPI_RPREC, down, 1,  &
!                    fz_tot(:,:,nz), nx*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
call mpi_sync_real_array( fx_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fy_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fz_tot, 1, MPI_SYNC_DOWN )

$endif

! Put fz_tot on uv-grid
fz_tot(1:nx,1:ny,1:nz) = interp_to_uv_grid( fz_tot(1:nx,1:ny,1:nz), 1 )

return
end subroutine force_tot
$endif

$if($DEBUG)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pressure_sync()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
use param, only : ld
implicit none

! Reset bogus values
p(:,:,nz) = p(:,:,nz-1)
dpdx(:,:,nz) = dpdx(:,:,nz-1)
dpdy(:,:,nz) = dpdy(:,:,nz-1)
dpdz(:,:,nz) = dpdz(:,:,nz-1)

$if($MPI)
!  Sync pressure
! call mpi_sendrecv (p(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    p(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (dpdx(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    dpdx(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (dpdy(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    dpdy(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (dpdz(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    dpdz(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)                   
call mpi_sync_real_array( p, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( dpdx, 1 , MPI_SYNC_DOWN )
call mpi_sync_real_array( dpdy, 1 , MPI_SYNC_DOWN )
call mpi_sync_real_array( dpdz, 1 , MPI_SYNC_DOWN )

$endif

return
end subroutine pressure_sync

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine RHS_sync()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : ld
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
implicit none

! Reset bogus values
RHSx(:,:,nz) = RHSx(:,:,nz-1)
RHSy(:,:,nz) = RHSy(:,:,nz-1)
RHSz(:,:,nz) = RHSz(:,:,nz-1)

$if($MPI)
!  Sync RHS
! call mpi_sendrecv (RHSx(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    RHSx(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (RHSy(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    RHSy(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)
! call mpi_sendrecv (RHSz(:,:,1), ld*ny, MPI_RPREC, down, 1,  &
!                    RHSz(:,:,nz), ld*ny, MPI_RPREC, up, 1,   &
!                    comm, status, ierr)                   
call mpi_sync_real_array( RHSx, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( RHSy, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( RHSz, 0 , MPI_SYNC_DOWN )

$endif

return
end subroutine RHS_sync
$endif

end subroutine inst_write

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkpoint (lun)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use param, only : nz
use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
$if ($DYN_TN)
use sgs_param, only:F_ee2,F_deedt2,ee_past
$endif
implicit none

integer, intent (in) :: lun

!---------------------------------------------------------------------

    $if ($DYN_TN) 
      write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                  Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),     &
                  F_QN(:,:,1:nz), F_NN(:,:,1:nz),                         &
                  F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)                
    $else
      write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                  Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),     &
                  F_QN(:,:,1:nz), F_NN(:,:,1:nz)
    $endif

end subroutine checkpoint


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_final(jt, lun_opt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by
!  inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use stat_defs, only : tavg_t, point_t
use param, only : tavg_calc, point_calc, point_nloc, spectra_calc
use param, only : dt
use param, only : use_cfl_dt, cfl
use cfl_mod, only : get_max_cfl

implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun
integer, parameter :: lun_default = 11
integer::i,fid
integer :: lun

logical :: opn

real(rprec) :: cfl_w

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

! Set the current cfl to a temporary (write) value based whether CFL is
! specified or must be computed
if( use_cfl_dt ) then
   cfl_w = cfl
else
   cfl_w = get_max_cfl()
endif

!  Update total_time.dat after simulation
!if ((cumulative_time) .and. (lun == lun_default)) then
  if (coord == 0) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)

    write(1, *) jt_total, total_time, total_time_dim, dt, cfl_w

    close(1)

  end if
!end if

!  Check if average quantities are to be recorded
if(tavg_calc) call tavg_finalize()

!  Check if spectra is to be computed
if(spectra_calc) call spectra_finalize()
 
return
end subroutine output_final

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine len_da_file(fname, lenrec, length)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine stats_init ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine allocates the memory for arrays used for statistical
!  calculations
use param, only : L_x,L_y,L_z,dx,dy,dz,nx,ny,nz,nsteps,coord,nproc,lbz,lh
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use param, only : spectra_calc, spectra_nloc, spectra_loc
use param, only : tavg_calc
use grid_defs, only : grid_t
use functions, only : cell_indx
use stat_defs, only : point_t, xplane_t, yplane_t, zplane_t
use stat_defs, only : tavg_t, tavg_zplane_t, spectra_t
use stat_defs, only : type_set
implicit none

include 'tecryte.h'

!character(120) :: cx,cy,cz
character(120) :: var_list
integer :: fid, i,j,k

logical :: exst

real(rprec), pointer, dimension(:) :: x,y,z

nullify(x,y,z)

x => grid_t % x
y => grid_t % y
z => grid_t % z

if( tavg_calc ) then

  allocate(tavg_t(nx,ny,lbz:nz))
  allocate(tavg_zplane_t(nz))

  ! Initialize the derived types tavg_t and tavg_zplane_t  
  do k=1,Nz
    do j=1, Ny
      do i=1, Nx
        call type_set( tavg_t(i,j,k), 0._rprec )
      enddo
    enddo

    call type_set( tavg_zplane_t(k), 0._rprec )

  enddo

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
    allocate( spectra_t(k)%power(lh) )
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

      !  Can't concatenate an empty string
      point_t(i) % fname=''
      call string_concat(point_t(i) % fname,'output/vel.x-')
      call string_concat(point_t(i) % fname, point_loc(i)%xyz(1))
      call string_concat(point_t(i) % fname,'.y-')
      call string_concat(point_t(i) % fname,point_loc(i)%xyz(2))
      call string_concat(point_t(i) % fname,'.z-')
      call string_concat(point_t(i) % fname,point_loc(i)%xyz(3))
      call string_concat(point_t(i) % fname,'.dat')
   
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

    point_t(i) % fname=''
    call string_concat(point_t(i) % fname,'output/vel.x-')
    call string_concat(point_t(i) % fname, point_loc(i)%xyz(1))
    call string_concat(point_t(i) % fname,'.y-')
    call string_concat(point_t(i) % fname,point_loc(i)%xyz(2))
    call string_concat(point_t(i) % fname,'.z-')
    call string_concat(point_t(i) % fname,point_loc(i)%xyz(3))
    call string_concat(point_t(i) % fname,'.dat')

    !  Add tecplot header if file does not exist
    inquire (file=point_t(i) % fname, exist=exst)
    if (.not. exst) then
      var_list = '"t (s)", "u", "v", "w"'
      call write_tecplot_header_xyline(point_t(i) % fname, 'rewind', var_list)
    endif 

  $endif
  
  enddo
endif

nullify(x,y,z)

return
end subroutine stats_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Load tavg.out files
use param, only : coord, dt, Nx, Ny, Nz
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
  if (coord == 0) then
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

return
end subroutine tavg_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_compute()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine collects the stats for each flow 
!  variable quantity
use types, only : rprec
use stat_defs, only : tavg_t, tavg_zplane_t, tavg_total_time
use param, only : nx,ny,nz,dt,lbz,jzmin,jzmax
use sim_param, only : u,v,w, dudz, dvdz, txx, txy, tyy, txz, tyz, tzz
use sim_param, only : fx, fy, fz, fxa, fya, fza
use functions, only : interp_to_w_grid

implicit none

!use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
integer :: i,j,k
real(rprec) :: u_p, v_p, w_p
real(rprec), allocatable, dimension(:,:,:) :: u_w, v_w

allocate(u_w(nx,ny,lbz:nz),v_w(nx,ny,lbz:nz))

!  Interpolate velocities to w-grid
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( v(1:nx,1:ny,lbz:nz), lbz )

do k=jzmin,jzmax  
  do j=1,ny
    do i=1,nx
   
      u_p = u_w(i,j,k)
      v_p = v_w(i,j,k) 
      w_p = w(i,j,k)
          
      ! === uv-grid variables ===
      tavg_t(i,j,k)%txx = tavg_t(i,j,k)%txx + txx(i,j,k) * dt
      tavg_t(i,j,k)%txy = tavg_t(i,j,k)%txy + txy(i,j,k) * dt
      tavg_t(i,j,k)%tyy = tavg_t(i,j,k)%tyy + tyy(i,j,k) * dt
      tavg_t(i,j,k)%tzz = tavg_t(i,j,k)%tzz + tzz(i,j,k) * dt

      ! === w-grid variables === 
      tavg_t(i,j,k)%u = tavg_t(i,j,k)%u + u_p * dt                    
      tavg_t(i,j,k)%v = tavg_t(i,j,k)%v + v_p * dt                         
      tavg_t(i,j,k)%w = tavg_t(i,j,k)%w + w_p * dt
      tavg_t(i,j,k)%u2 = tavg_t(i,j,k)%u2 + u_p * u_p * dt
      tavg_t(i,j,k)%v2 = tavg_t(i,j,k)%v2 + v_p * v_p * dt
      tavg_t(i,j,k)%w2 = tavg_t(i,j,k)%w2 + w_p * w_p * dt
      tavg_t(i,j,k)%uw = tavg_t(i,j,k)%uw+ u_p * w_p * dt
      tavg_t(i,j,k)%vw = tavg_t(i,j,k)%vw + v_p * w_p * dt
      tavg_t(i,j,k)%uv = tavg_t(i,j,k)%uv + u_p * v_p * dt

      tavg_t(i,j,k)%dudz = tavg_t(i,j,k)%dudz + dudz(i,j,k) * dt
      tavg_t(i,j,k)%dvdz = tavg_t(i,j,k)%dvdz + dvdz(i,j,k) * dt

      tavg_t(i,j,k)%txz = tavg_t(i,j,k)%txz + txz(i,j,k) * dt
      tavg_t(i,j,k)%tyz = tavg_t(i,j,k)%tyz + tyz(i,j,k) * dt
      
    enddo
  enddo
enddo

do k=1,jzmax        ! these do not have a k=0 level (needed by coord==0)
                    ! this shoud be fixed in the near future...
  do j=1,ny
    do i=1,nx
 
      ! === uv-grid variables === 
      ! Includes both induced (IBM) and applied (RNS, turbines, etc.) forces
      tavg_t(i,j,k)%fx = tavg_t(i,j,k)%fx + (fx(i,j,k) + fxa(i,j,k)) * dt 
      tavg_t(i,j,k)%fy = tavg_t(i,j,k)%fy + (fy(i,j,k) + fya(i,j,k)) * dt 

      ! === w-grid variables === 
      tavg_t(i,j,k)%cs_opt2 = tavg_t(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * dt 

      ! Includes both induced (IBM) and applied (RNS, turbines, etc.) forces
      tavg_t(i,j,k)%fz = tavg_t(i,j,k)%fz + (fz(i,j,k) + fza(i,j,k)) * dt
      
    enddo
  enddo
enddo


deallocate( u_w, v_w )

! Update tavg_total_time for variable time stepping
tavg_total_time = tavg_total_time + dt

return

end subroutine tavg_compute

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_finalize()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use grid_defs, only : grid_t !x,y,z
use stat_defs, only : tavg_t, tavg_zplane_t, tavg_total_time, tavg
use stat_defs, only : rs, rs_t, rs_zplane_t, cnpy_zplane_t 
use stat_defs, only : operator(.DIV.), operator(.MUL.)
use stat_defs, only :  operator(.ADD.), operator(.SUB.)
use stat_defs, only : type_set, type_zero_bogus
use stat_defs, only : tavg_interp_to_uv_grid, tavg_interp_to_w_grid
use stat_defs, only : rs_compute, cnpy_tavg_mul
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z, nz_tot
$if($MPI)
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use param, only : MPI_RPREC, coord_of_rank, &
     rank_of_coord, comm, &
     ierr, down, up, status
use stat_defs, only : rs, tavg
$endif
$if($LVLSET)
use level_set_base, only : phi
$endif

implicit none

include 'tecryte.h'

character (*), parameter :: sub_name = mod_name // '.tavg_finalize'
character(64) :: fname_out, fname_vel, &
     fname_vel2, fname_ddz, &
     fname_tau, fname_f, &
     fname_rs, fname_cs

character(64) :: fname_vel_zplane, fname_vel2_zplane, &
  fname_ddz_zplane, fname_tau_zplane, fname_f_zplane, &
  fname_rs_zplane, fname_cs_zplane, fname_cnpy_zplane
  
integer :: i,j,k

$if($MPI)
character(64) :: temp

integer :: MPI_RS, MPI_CNPY, MPI_TAVG
integer :: rs_type(1), rs_block(1), rs_disp(1)
integer :: cnpy_type(1), cnpy_block(1), cnpy_disp(1)
integer :: tavg_type(1), tavg_block(1), tavg_disp(1)

! Definitions for reconstructing z-planar averaged data
integer :: sendsize
integer, allocatable, dimension(:) :: recvsize, recvstride
real(rprec), allocatable, dimension(:) :: z_tot, zw_tot
type(rs), allocatable, dimension(:) :: rs_zplane_tot_t
type(rs), allocatable, dimension(:) :: cnpy_zplane_tot_t
type(tavg), allocatable, dimension(:) :: tavg_zplane_tot_t

$endif

logical :: opn

type(rs) :: cnpy_avg_t
type(tavg) :: tavg_avg_t

real(rprec) :: favg

$if($LVLSET)
real(rprec) :: fx_global, fy_global, fz_global
$endif

real(rprec), pointer, dimension(:) :: x,y,z,zw

favg = real(nx*ny,kind=rprec)

nullify(x,y,z,zw)

allocate(rs_t(nx,ny,lbz:nz), rs_zplane_t(nz))
allocate(cnpy_zplane_t(nz))

x => grid_t % x
y => grid_t % y
z => grid_t % z
zw => grid_t % zw

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

! Zero bogus values
call type_zero_bogus( tavg_t(:,:,nz) )

$if($MPI)
! Build MPI types for derived types
rs_type = MPI_RPREC
rs_block = 6 ! Number of rs subtypes
rs_disp = 0

cnpy_type = MPI_RPREC
cnpy_block = 6 ! Number of cnpy subtypes
cnpy_disp = 0

tavg_type = MPI_RPREC
tavg_block = 21 ! Number of tavg subtypes
tavg_disp = 0

!  Create MPI type structures consistent with the derived types
call MPI_TYPE_STRUCT(1, rs_block, rs_disp, rs_type, MPI_RS, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_RS:', ierr)
Call MPI_Type_commit(MPI_RS,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_RS:', ierr)


call MPI_TYPE_STRUCT(1, cnpy_block, cnpy_disp, cnpy_type, MPI_CNPY, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_CNPY:', ierr)
Call MPI_Type_commit(MPI_CNPY,ierr)  
if(ierr /= 0) call error(sub_name,'Error in committing MPI_CNPY:', ierr)
   
call MPI_TYPE_STRUCT(1, tavg_block, tavg_disp, tavg_type, MPI_TAVG, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_TAVG:', ierr)
Call MPI_Type_commit(MPI_TAVG,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_TAVG:', ierr)

!  Allocate space only on base processor for assembled z-plane data
! *_tot_t is the assembled data without the overlap nodes (the final stuff that is outputted)
! All processes need to allocate the space even though it is not directly used
allocate(z_tot(nz_tot))
allocate(zw_tot(nz_tot))
allocate(rs_zplane_tot_t(nz_tot))
allocate(cnpy_zplane_tot_t(nz_tot))
allocate(tavg_zplane_tot_t(nz_tot))

$endif

!  Perform time averaging operation
!  tavg_t = tavg_t / tavg_total_time
do k=jzmin,jzmax
  do j=1, Ny
    do i=1, Nx
      tavg_t(i,j,k) = tavg_t(i,j,k) .DIV. tavg_total_time
    enddo
  enddo
enddo

!  Sync entire tavg_t structure
$if($MPI)
call mpi_sendrecv (tavg_t(:,:,1), nx*ny, MPI_TAVG, down, 1,  &
                   tavg_t(:,:,nz), nx*ny, MPI_TAVG, up, 1,   &
                   comm, status, ierr)
call mpi_sendrecv (tavg_t(:,:,nz-1), nx*ny, MPI_TAVG, up, 2,  &
                   tavg_t(:,:,0), nx*ny, MPI_TAVG, down, 2,   &
                   comm, status, ierr)
$endif

! Anything with velocity is on the w-grid so these
!   values for coord==0 at k=1 should be zeroed (otherwise bogus)
if (jzmin==0) then  !coord==0 only
  tavg_t(:,:,0:1)%u = 0.0_rprec
  tavg_t(:,:,0:1)%v = 0.0_rprec
  tavg_t(:,:,0:1)%w = 0.0_rprec
  tavg_t(:,:,0:1)%u2 = 0.0_rprec
  tavg_t(:,:,0:1)%v2 = 0.0_rprec
  tavg_t(:,:,0:1)%w2 = 0.0_rprec
  tavg_t(:,:,0:1)%uw = 0.0_rprec
  tavg_t(:,:,0:1)%vw = 0.0_rprec
  tavg_t(:,:,0:1)%uv = 0.0_rprec
endif

! Interpolate between grids where necessary
!tavg_t = tavg_interp_to_uv_grid( tavg_t )
tavg_t = tavg_interp_to_w_grid( tavg_t )

! Anything with velocity is on the w-grid so these
!   values for coord==0 at k=1 should be zeroed (otherwise bogus)
if (jzmin==0) then  !coord==0 only
  tavg_t(:,:,0:1)%fx = 0.0_rprec 
  tavg_t(:,:,0:1)%fy = 0.0_rprec   

  tavg_t(:,:,0:1)%txx = 0.0_rprec    
  tavg_t(:,:,0:1)%txy = 0.0_rprec    
  tavg_t(:,:,0:1)%tyy = 0.0_rprec    
  tavg_t(:,:,0:1)%tzz = 0.0_rprec
endif

!  Average over z-planes
do k=1,nz
  
  !  Initialize to 0 for summations
  call type_set( tavg_zplane_t(k), 0._rprec )

  do j=1, Ny
    do i=1, Nx

      tavg_zplane_t(k) = tavg_zplane_t(k) .ADD. tavg_t(i,j,k)

    enddo
  enddo

  !  Divide by number of summation points 
  tavg_zplane_t(k) = tavg_zplane_t(k) .DIV. favg

enddo

! Compute the Reynolds stresses: bar(u_i * u_j) - bar(u_i) * bar(u_j)
rs_t = rs_compute( tavg_t , lbz)

! Compute planar averaged Reynolds stress
do k = 1, nz

  !  Initialize to 0
  call type_set( rs_zplane_t(k), 0._rprec)

  do j = 1, ny
    do i = 1, nx
      rs_zplane_t(k) = rs_zplane_t(k) .ADD. rs_t(i,j,k) 
    enddo    
  enddo

  rs_zplane_t(k) = rs_zplane_t(k) .DIV. favg
  
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


! ----- Write all the 3D data -----
! -- Work around for large data: only write 3 variables at a time. Since things are
! -- written in block format this can be done with mulitple calls to
! -- write_real_data_3D.
$if ($LVLSET)

call write_tecplot_header_ND(fname_vel, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<u>","<v>","<w>"', numtostr(coord, 6), 2)
!  write phi and x,y,z
call write_real_data_3D(fname_vel, 'append', 'formatted', 4, nx, ny, nz, &
     (/ phi(1:nx,1:ny,1:nz), &
     tavg_t(:,:,1:nz) % u, &
     tavg_t(:,:,1:nz) % v, &
     tavg_t(:,:,1:nz) % w /), &
     4, x, y, zw(1:nz))
$else

call write_tecplot_header_ND(fname_vel, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<u>","<v>","<w>"', numtostr(coord, 6), 2)
call write_real_data_3D(fname_vel, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t(:,:,1:nz) % u, &
  tavg_t(:,:,1:nz) % v, &
  tavg_t(:,:,1:nz) % w /), &
  4, x, y, zw(1:nz))

$endif

$if($LVLSET)

call write_tecplot_header_ND(fname_vel2, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
   numtostr(coord,6), 2)
!  write phi and x,y,z
call write_real_data_3D(fname_vel2, 'append', 'formatted', 7, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz), &
  tavg_t(:,:,1:nz) % u2, &
  tavg_t(:,:,1:nz) % v2, &
  tavg_t(:,:,1:nz) % w2, &
  tavg_t(:,:,1:nz) % uw, &
  tavg_t(:,:,1:nz) % vw, &
  tavg_t(:,:,1:nz) % uv /), &
  4, x, y, zw(1:nz))

$else

call write_tecplot_header_ND(fname_vel2, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
   numtostr(coord,6), 2)
call write_real_data_3D(fname_vel2, 'append', 'formatted', 6, nx, ny, nz, &
  (/ tavg_t(:,:,1:nz) % u2, &
  tavg_t(:,:,1:nz) % v2, &
  tavg_t(:,:,1:nz) % w2, &
  tavg_t(:,:,1:nz) % uw, &
  tavg_t(:,:,1:nz) % vw, &
  tavg_t(:,:,1:nz) % uv /), &
  4, x, y, zw(1:nz))

$endif

$if($LVLSET)  

call write_tecplot_header_ND(fname_ddz, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
!  write phi and x,y,z
call write_real_data_3D(fname_ddz, 'append', 'formatted', 3, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz), &
  tavg_t(:,:,1:nz) % dudz, &
  tavg_t(:,:,1:nz) % dvdz /), &
  4, x, y, zw(1:nz))

$else

call write_tecplot_header_ND(fname_ddz, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
call write_real_data_3D(fname_ddz, 'append', 'formatted', 2, nx, ny, nz, &
  (/ tavg_t(:,:,1:nz) % dudz, &
  tavg_t(:,:,1:nz) % dvdz /), &
  4, x, y, zw(1:nz))

$endif

$if($LVLSET)

call write_tecplot_header_ND(fname_tau, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
   numtostr(coord,6), 2)  
!  write phi and x,y,z
call write_real_data_3D(fname_tau, 'append', 'formatted', 7, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz), &
  tavg_t(:,:,1:nz) % txx, &
  tavg_t(:,:,1:nz) % txy, &
  tavg_t(:,:,1:nz) % tyy, &
  tavg_t(:,:,1:nz) % txz, &
  tavg_t(:,:,1:nz) % tyz, &
  tavg_t(:,:,1:nz) % tzz /), &
  4, x, y, zw(1:nz))

$else

call write_tecplot_header_ND(fname_tau, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
   numtostr(coord,6), 2)
call write_real_data_3D(fname_tau, 'append', 'formatted', 6, nx, ny, nz, &
  (/ tavg_t(:,:,1:nz) % txx, &
  tavg_t(:,:,1:nz) % txy, &
  tavg_t(:,:,1:nz) % tyy, &
  tavg_t(:,:,1:nz) % txz, &
  tavg_t(:,:,1:nz) % tyz, &
  tavg_t(:,:,1:nz) % tzz /), &
  4, x, y, zw(1:nz))

$endif  

$if($LVLSET)

call write_tecplot_header_ND(fname_f, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', &
   numtostr(coord,6), 2)
!  write phi and x,y,z
call write_real_data_3D(fname_f, 'append', 'formatted', 4, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz), &
  tavg_t(:,:,1:nz) % fx, &
  tavg_t(:,:,1:nz) % fy, &
  tavg_t(:,:,1:nz) % fz /), &
  4, x, y, zw(1:nz))

  $if($MPI)

  call mpi_allreduce(sum(tavg_t(1:nx,1:ny,1:nz-1)%fx), fx_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce(sum(tavg_t(1:nx,1:ny,1:nz-1)%fy), fy_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  call mpi_allreduce(sum(tavg_t(1:nx,1:ny,1:nz-1)%fz), fz_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)

  $else

  fx_global = sum(tavg_t(1:nx,1:ny,1:nz-1)%fx)
  fy_global = sum(tavg_t(1:nx,1:ny,1:nz-1)%fy)
  fz_global = sum(tavg_t(1:nx,1:ny,1:nz-1)%fz)

  $endif
  
  if (coord == 0) then
    open(unit = 1, file = "output/force_total_avg.dat", status="unknown", position="rewind") 
    write(1,'(a,3e15.6)') '<fx>, <fy>, <fz> : ', fx_global, fy_global, fz_global
    close(1)
  endif

$else

call write_tecplot_header_ND(fname_f, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', &
   numtostr(coord,6), 2)
call write_real_data_3D(fname_f, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t(:,:,1:nz) % fx, &
  tavg_t(:,:,1:nz) % fy, &
  tavg_t(:,:,1:nz) % fz /), &
  4, x, y, zw(1:nz))

$endif
  
$if($LVLSET)

call write_tecplot_header_ND(fname_rs, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', &
   numtostr(coord,6), 2)  
!  write phi and x,y,z
call write_real_data_3D(fname_rs, 'append', 'formatted', 7, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz), &
  rs_t(:,:,1:nz)%up2, &
  rs_t(:,:,1:nz)%vp2, &
  rs_t(:,:,1:nz)%wp2, &
  rs_t(:,:,1:nz)%upwp, &
  rs_t(:,:,1:nz)%vpwp, &
  rs_t(:,:,1:nz)%upvp /), &
  4, x, y, zw(1:nz))

$else

call write_tecplot_header_ND(fname_rs, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', &
   numtostr(coord,6), 2)
call write_real_data_3D(fname_rs, 'append', 'formatted', 6, nx, ny, nz, &
  (/ rs_t(:,:,1:nz)%up2, &
  rs_t(:,:,1:nz)%vp2, &
  rs_t(:,:,1:nz)%wp2, &
  rs_t(:,:,1:nz)%upwp, &
  rs_t(:,:,1:nz)%vpwp, &
  rs_t(:,:,1:nz)%upvp /), &
  4, x, y, zw(1:nz))

$endif

$if($LVLSET)

call write_tecplot_header_ND(fname_cs, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "phi", "<cs2>"', numtostr(coord,6), 2)
!  write phi and x,y,z
call write_real_data_3D(fname_cs, 'append', 'formatted', 2, nx, ny, nz, &
  (/ phi(1:nx,1:ny,1:nz), &
  tavg_t(:,:,1:nz) % cs_opt2 /), &
  4, x, y, zw(1:nz))

$else

call write_tecplot_header_ND(fname_cs, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<cs2>"', numtostr(coord,6), 2)
call write_real_data_3D(fname_cs, 'append', 'formatted', 1, nx, ny, nz, &
  (/ tavg_t(:,:,1:nz)% cs_opt2 /), &
  4, x, y, zw(1:nz))

$endif

! Reconstruct z-plane data 
$if($MPI)

  allocate(recvsize(nproc),recvstride(nproc))

  sendsize=nz-1
  recvsize=nz-1
  recvstride=coord_of_rank*(nz-1)

  ! Set very bottom values
  if( coord == 0 ) then

     z_tot(1) = z(1)
     zw_tot(1) = zw(1)
     rs_zplane_tot_t(1) = rs_zplane_t(1)
     cnpy_zplane_tot_t(1) = cnpy_zplane_t(1)
     tavg_zplane_tot_t(1) = tavg_zplane_t(1)

  endif

  call mpi_gatherv( z(2), sendsize, MPI_RPREC, &
       z_tot(2), recvsize, recvstride, MPI_RPREC, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( zw(2), sendsize, MPI_RPREC, &
       zw_tot(2), recvsize, recvstride, MPI_RPREC, &
       rank_of_coord(0), comm, ierr)
  
  call mpi_gatherv( rs_zplane_t(2), sendsize, MPI_RS, &
       rs_zplane_tot_t(2), recvsize, recvstride, MPI_RS, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( cnpy_zplane_t(2), sendsize, MPI_CNPY, &
       cnpy_zplane_tot_t(2), recvsize, recvstride, MPI_CNPY, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( tavg_zplane_t(2), sendsize, MPI_TAVG, &
       tavg_zplane_tot_t(2), recvsize, recvstride, MPI_TAVG, &
       rank_of_coord(0), comm, ierr)

  deallocate(recvsize, recvstride)

  call MPI_Type_free (MPI_RS, ierr)
  call MPI_Type_free (MPI_CNPY, ierr)
  call mpi_type_free (MPI_TAVG, ierr)    

  ! Write reconstructed data only if bottom processor
  if( coord == 0 ) then
  
    call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz_tot /), &
      '"z", "<u>","<v>","<w>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot_t % u, tavg_zplane_tot_t % v, tavg_zplane_tot_t % w /), 0, zw_tot)

    call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz_tot/), &
      '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
      numtostr(coord,6), 2)
    call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot_t % u2, tavg_zplane_tot_t % v2, tavg_zplane_tot_t % w2, &
      tavg_zplane_tot_t % uw, tavg_zplane_tot_t % vw, tavg_zplane_tot_t % uv /), &
      0, zw_tot) 
  
    call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz_tot/), &
      '"z", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, Nz_tot, &
      (/ tavg_zplane_tot_t % dudz, tavg_zplane_tot_t % dvdz /), 0, zw_tot)
  
    call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
      numtostr(coord,6), 2)  
    call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot_t % txx, tavg_zplane_tot_t % txy, tavg_zplane_tot_t % tyy, &
      tavg_zplane_tot_t % txz, tavg_zplane_tot_t % tyz, tavg_zplane_tot_t % tzz /), 0, zw_tot) 
  
    call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot_t % fx, tavg_zplane_tot_t % fy, tavg_zplane_tot_t % fz /), 0, zw_tot)  
  
    call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ rs_zplane_tot_t % up2, rs_zplane_tot_t%vp2, rs_zplane_tot_t%wp2, &
      rs_zplane_tot_t%upwp, rs_zplane_tot_t%vpwp, rs_zplane_tot_t%upvp /), 0, zw_tot)    
 
    call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ cnpy_zplane_tot_t % up2, cnpy_zplane_tot_t%vp2, cnpy_zplane_tot_t%wp2, &
      cnpy_zplane_tot_t%upwp, cnpy_zplane_tot_t%vpwp, cnpy_zplane_tot_t%upvp /), 0, zw_tot)         
      
    call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz_tot/), &
      '"z", "<cs2>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, Nz_tot, &
      (/ tavg_zplane_tot_t % cs_opt2 /), 0, zw_tot)        

    deallocate(z_tot, zw_tot)      
    deallocate(tavg_zplane_tot_t)
    deallocate(rs_zplane_tot_t)
    deallocate(cnpy_zplane_tot_t)
  
  endif

$else

call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz /), &
   '"z", "<u>","<v>","<w>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane_t % u, tavg_zplane_t % v, tavg_zplane_t % w /), 0, zw(1:nz))

call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz/), &
   '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
   numtostr(coord,6), 2)
call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane_t % u2, tavg_zplane_t % v2, tavg_zplane_t % w2, &
  tavg_zplane_t % uw, tavg_zplane_t % vw, tavg_zplane_t % uv /), 0, zw(1:nz)) 
  
call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz/), &
   '"z", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, nz, &
  (/ tavg_zplane_t % dudz, tavg_zplane_t % dvdz /), 0, zw(1:nz))
  
call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
   numtostr(coord,6), 2)  
call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane_t % txx, tavg_zplane_t % txy, tavg_zplane_t % tyy, &
  tavg_zplane_t % txz, tavg_zplane_t % tyz, tavg_zplane_t % tzz /), 0, zw(1:nz)) 
  
call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', numtostr(coord,6), 2)
call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane_t % fx, tavg_zplane_t % fy, tavg_zplane_t % fz /), 0, zw(1:nz))  
  
call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, nz, &
  (/ rs_zplane_t % up2, rs_zplane_t%vp2, rs_zplane_t%wp2, &
  rs_zplane_t%upwp, rs_zplane_t%vpwp, rs_zplane_t%upvp /), 0, zw(1:nz))

call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, nz, &
  (/ cnpy_zplane_t % up2, cnpy_zplane_t%vp2, cnpy_zplane_t%wp2, &
  cnpy_zplane_t%upwp, cnpy_zplane_t%vpwp, cnpy_zplane_t%upvp /), 0, zw(1:nz))  
  
call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz/), &
   '"z", "<cs2>"', numtostr(coord,6), 2)
call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, nz, &
  (/ tavg_zplane_t % cs_opt2 /), 0, zw(1:nz))    
  
$endif

deallocate(tavg_t, tavg_zplane_t, rs_t, rs_zplane_t, cnpy_zplane_t)

nullify(x,y,z,zw)

return
end subroutine tavg_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use messages
use param, only : coord, dt, spectra_nloc, lh, nx
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
  if (coord == 0) then
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_compute()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      power(i) = uhat(i)**2 + uhat(Nx-i+2)**2
    enddo
    power(lh) = uhat(lh)**2 ! Nyquist

    ! Sum jth component and normalize
    spectra_t(k) % power = spectra_t(k) % power + power / nx

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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_finalize()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param, only : lh, spectra_nloc, spectra_loc
use fft, only : kx
use stat_defs, only : spectra_t, spectra_total_time
implicit none

include 'tecryte.h'

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
    '"k", "E(k)"', numtostr(k, 6), 2 ) 
  call write_real_data_1D(fname, 'append', 'formatted', 1, lh-1, &
    (/ spectra_t(k) % power(1:lh-1) /), 0, (/ kx(1:lh-1,1) /))
    
  $if ($MPI)
  endif
  $endif

enddo  

!  Write data to spectra.out
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

! write the accumulated time and power
write (1) spectra_total_time
do k=1, spectra_nloc
  write (1) spectra_t(k) % power
enddo
close(1)
  
deallocate(spectra_t)

return

end subroutine spectra_finalize

end module io
