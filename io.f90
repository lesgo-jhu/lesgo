!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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

!///////////////////////////////////////////////////////////////////////////////
module io
!///////////////////////////////////////////////////////////////////////////////
use types,only:rprec
use param, only : ld,nx,ny,nz,nz_tot,path,coord,rank,nproc,jt_total
use param, only : total_time,total_time_dim,lbz,jzmin,jzmax,cumulative_time,fcumulative_time
use sim_param,only:w,dudz,dvdz
use sgs_param,only:Cs_opt2
use string_util
use messages
$if ($MPI)
use mpi
$endif
implicit none
save
private

public jt_total, openfiles, energy, output_loop,output_final,output_init

character (*), parameter :: mod_name = 'io'

! Output file id's (see README for assigned values)
!integer :: ke_fid
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine openfiles()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : use_cfl_dt, dt, cfl_f
implicit none
include 'tecryte.h'

logical :: exst

! Temporary values used to read time step and CFL from file
real(rprec) :: dt_r, cfl_r

! Open output files using tecryte library
! Kinetic energy (check_ke.dat)
!ke_fid = open_file( path // 'output/check_ke.dat', 'append', 'formatted' )

!dt_r=0.0_rprec
if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then

    open (1, file=fcumulative_time)   
    read(1, *) jt_total, total_time, total_time_dim, dt_r, cfl_r
    close (1)
 
    ! Update dynamic time stepping info if required; otherwise discard.
    if( use_cfl_dt ) then
      dt = dt_r
      cfl_f = cfl_r
    endif
   
  else  !--assume this is the first run on cumulative time
    if( coord == 0 ) then
      write (*, *) '--> Assuming jt_total = 0, total_time = 0.0'
    endif
    jt_total = 0
    total_time = 0.0_rprec
    total_time_dim = 0.0_rprec
  end if

end if

!! Update dynamic time stepping info if required; otherwise discard.
!if( use_cfl_dt ) then
!   dt = dt_r
!   cfl_f = cfl_r
!endif

return
end subroutine openfiles

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine energy (ke)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types,only:rprec
use param
use sim_param,only:u,v,w
implicit none
integer :: jx, jy, jz, nan_count

real(rprec)::KE,temp_w
$if ($MPI)
  real (rprec) :: ke_global
$endif

! Initialize variables
ke=0._rprec
do jz=1,nz-1
    do jy=1,ny
        do jx=1,nx
            
            temp_w = 0.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
            ke = ke + (u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)
                        
        end do
    end do
end do

! Perform spatial averaging
ke = ke*0.5_rprec/(nx*ny*(nz-1))

$if ($MPI)
   call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
   if (rank == 0) then  !--note its rank here, not coord
   open(2,file=path // 'output/check_ke.dat',status='unknown',form='formatted',position='append')
   ke = ke_global/nproc
   write(2,*) total_time,ke
   close(2)
   end if
$else
open(2,file=path // 'output/check_ke.dat',status='unknown',form='formatted',position='append')
write(2,*) total_time,ke
close(2)
$endif

end subroutine energy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_loop()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This subroutine is called every time step and acts as a driver for 
!  computing statistics and outputing instantaneous data. No actual
!  calculations are performed here.
!
use param, only : jt_total, dt
use param, only : checkpoint_data, checkpoint_nskip
use param, only : tavg_calc, tavg_nstart, tavg_nend, tavg_nskip
use param, only : domain_calc, domain_nstart, domain_nend, domain_nskip
use param, only : xplane_calc, xplane_nstart, xplane_nend, xplane_nskip
use param, only : yplane_calc, yplane_nstart, yplane_nend, yplane_nskip
use param, only : zplane_calc, zplane_nstart, zplane_nend, zplane_nskip
use stat_defs, only : tavg_initialized,tavg_dt
$if (not $FFTW3)
use stat_defs, only: spectra_initialized,spectra_dt
$endif
implicit none

! Determine if we are to checkpoint intermediate times
if( checkpoint_data ) then
   ! Now check if data should be checkpointed this time step
   ! Don't checkpoint for last time step since this is done in output_final
!   if ( jt_total < nsteps .and. modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
   if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
endif 

!  Determine if time summations are to be calculated
if (tavg_calc) then

  ! Are we between the start and stop timesteps?
  if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then

    ! Every timestep (between nstart and nend), add to tavg_dt
    tavg_dt = tavg_dt + dt

    ! Are we at the beginning or a multiple of nstart?
    if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then

      ! Check if we have initialized tavg
      if (.not.tavg_initialized) then
        if (coord == 0) then
          write(*,*) '-------------------------------'   
          write(*,"(1a,i9,1a,i9)") 'Starting running time summation from ', tavg_nstart, ' to ', tavg_nend
          write(*,*) '-------------------------------'   
        endif  ! coord==0
      
        call tavg_init()
      else
        call tavg_compute ()
      endif  ! init
     
    endif  ! mod of nskip

  endif  ! between nstart and nend

endif  ! tavg_calc
     
!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
  if(jt_total >= domain_nstart .and. jt_total <= domain_nend .and. ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
    if(jt_total == domain_nstart) then
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
  if(jt_total >= xplane_nstart .and. jt_total <= xplane_nend .and. ( mod(jt_total-xplane_nstart,xplane_nskip)==0) ) then
    if(jt_total == xplane_nstart) then
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
  if(jt_total >= yplane_nstart .and. jt_total <= yplane_nend .and. ( mod(jt_total-yplane_nstart,yplane_nskip)==0) ) then
    if(jt_total == yplane_nstart) then
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
  if(jt_total >= zplane_nstart .and. jt_total <= zplane_nend .and. ( mod(jt_total-zplane_nstart,zplane_nskip)==0) ) then
    if(jt_total == zplane_nstart) then
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
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use param, only : dx,dy,dz
use grid_defs, only : grid
use sim_param, only : u,v,w
use stat_defs, only : xplane, yplane, zplane
$if($MPI)
use param, only :ny,nz,comm,ierr
$endif
use sgs_param, only : Nu_t
implicit none
integer, intent(IN) :: itype

character (128) :: fname
integer :: n, i, j, k

real(rprec), allocatable, dimension(:,:,:) :: w_uv

!  Allocate space for the interpolated w values
allocate(w_uv(nx,ny,lbz:nz))

!  Make sure w has been interpolated to uv-grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

if(itype==1) then

!  Instantaneous write for entire domain
elseif(itype==2) then

  !////////////////////////////////////////////
  !/// WRITE VELOCITY                       ///
  !////////////////////////////////////////////

  $if( $BINARY )
  call string_splice( fname, path // 'output/binary_vel.', jt_total,'.dat')
  $endif

  $if ($MPI)
  call string_concat( fname, '.c', coord )
  $endif
 
  $if($BINARY)
  open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*(nz-1)*rprec/4)
  write(13,rec=1) u(:nx,:ny,1:nz-1)
  write(13,rec=2) v(:nx,:ny,1:nz-1)
  write(13,rec=3) w_uv(:nx,:ny,1:nz-1)
  write(13,rec=4) cs_opt2(:nx,:ny,1:nz-1)
  write(13,rec=5) Nu_t(:nx,:ny,1:nz-1) 
  close(13)
  $endif

  $if($MPI)
    ! Ensure that all processes finish before attempting to write 
    ! additional files. Otherwise it may flood the system with 
    ! too many I/O requests and crash the process 
    call mpi_barrier( comm, ierr )
  $endif

!  Write instantaneous x-plane values
elseif(itype==3) then

    if (coord>=nproc/4-1 .and. coord<=nproc/4*3+1) then
    $if ($BINARY)
    call string_splice(fname,path //'output/binary_vel_x_',jt_total,'.dat')
    $else
    call string_splice(fname,path //'output/vel_x_',jt_total,'.dat')
    $endif

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

! 0.5_rprec * (var(:,:,2:ubz) + var(:,:,1:ubz-1))
    $if ($BINARY)
    open(unit=13,file=fname,form='unformatted',convert='big_endian',access='direct',recl=1*ny*nz*rprec/4)
    write(13,rec=1)    u(nx*3/8,1:ny,1:nz)
    write(13,rec=2)    v(nx*3/8,1:ny,1:nz)
    write(13,rec=3) w_uv(nx*3/8,1:ny,1:nz)
    write(13,rec=4)    u(nx*4/8,1:ny,1:nz)
    write(13,rec=5)    v(nx*4/8,1:ny,1:nz)
    write(13,rec=6) w_uv(nx*4/8,1:ny,1:nz)
    write(13,rec=7)    u(nx*5/8,1:ny,1:nz)
    write(13,rec=8)    v(nx*5/8,1:ny,1:nz)
    write(13,rec=9) w_uv(nx*5/8,1:ny,1:nz)
    write(13,rec=10)   u(nx*6/8,1:ny,1:nz)
    write(13,rec=11)   v(nx*6/8,1:ny,1:nz)
    write(13,rec=12)w_uv(nx*6/8,1:ny,1:nz)
    close(13)
    $endif
    endif


!  Write instantaneous y-plane values
elseif(itype==4) then
    
    if (coord>=nproc/4-1 .and. coord<=nproc/4*3+1) then 
    $if ($BINARY)
    call string_splice(fname,path //'output/binary_vel_y_',jt_total,'.dat')
    $else
    call string_splice(fname,path //'output/vel_y_',jt_total,'.dat')
    $endif

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

! 0.5_rprec * (var(:,:,2:ubz) + var(:,:,1:ubz-1))
    $if ($BINARY)
    open(unit=13,file=fname,form='unformatted',convert='big_endian',access='direct',recl=nx*1*nz*rprec/4)
    write(13,rec=1)    u(1:nx,ny/2,1:nz)
    write(13,rec=2)    v(1:nx,ny/2,1:nz)
    write(13,rec=3) w_uv(1:nx,ny/2,1:nz)
    close(13) 
    $endif
    endif
  
!  Write instantaneous z-plane values
elseif(itype==5) then

!  Loop over all zplane locations
  do k=1,nz
 
    if (k+coord*(nz-1)==(nz_tot-1)/2) then
 
    call string_splice( fname, path // 'output/binary_vel_z_', jt_total, '.dat')

    open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*1*rprec/4)
    write(13,rec=1) u(1:nx,1:ny,k)
    write(13,rec=2) v(1:nx,1:ny,k)
    write(13,rec=3) w_uv(1:nx,1:ny,k)
    close(13)
   
    endif  
    
  enddo  
  
else
  write(*,*) 'Error: itype not specified properly to inst_write!'
  stop
endif

deallocate(w_uv)
end subroutine inst_write

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkpoint ()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : nz, checkpoint_file, tavg_calc
$if($MPI)
use param, only : comm,ierr
$endif
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
use param, only : jt_total, total_time, total_time_dim, dt,use_cfl_dt, cfl
use cfl_util, only : get_max_cfl
use stat_defs, only : tavg_initialized
use string_util, only : string_concat
implicit none

character(64) :: fname
real(rprec) :: cfl_w

fname = checkpoint_file
$if ($MPI)
call string_concat( fname, '.c', coord )
$endif

!  Open vel.out (lun_default in io) for final output
$if ($WRITE_BIG_ENDIAN)
open(11,file=fname,form='unformatted', convert='big_endian', status='unknown', position='rewind')
$elseif ($WRITE_LITTLE_ENDIAN)
open(11,file=fname,form='unformatted', convert='little_endian', status='unknown', position='rewind')
$else
open(11,file=fname,form='unformatted', status='unknown', position='rewind')
$endif

write (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
     RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
     Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),     &
     F_QN(:,:,1:nz), F_NN(:,:,1:nz)

! Close the file to ensure that the data is flushed and written to file
close(11)

$if($MPI)
call mpi_barrier( comm, ierr )
$endif

$if ($DYN_TN) 
! Write running average variables to file
  fname = path // 'dyn_tn.out'
  $if ($MPI)
  call string_concat( fname, '.c', coord)
  $endif
  
  $if ($WRITE_BIG_ENDIAN)
  open(unit=13,file=fname,form='unformatted',position='rewind',convert='big_endian')
  $elseif ($WRITE_LITTLE_ENDIAN)
  open(unit=13,file=fname,form='unformatted',position='rewind',convert='little_endian')
  $else
  open(unit=13,file=fname,form='unformatted',position='rewind')
  $endif

  write(13) F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)

  close(13)

  $if($MPI)
  call mpi_barrier( comm, ierr )
  $endif
  
$endif

! Checkpoint time averaging restart data
if( tavg_calc .and. tavg_initialized ) call tavg_checkpoint()

! Write time and current simulation state
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

return
end subroutine checkpoint


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_final()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use stat_defs, only : tavg_initialized
use param, only : tavg_calc
$if(not $FFTW3)
use stat_defs, only : spectra_initialized
use param, only : spectra_calc
$endif
implicit none

! Perform final checkpoing
call checkpoint()

!  Check if average quantities are to be recorded
if(tavg_calc .and. tavg_initialized ) call tavg_finalize()

$if(not $FFTW3)
!  Check if spectra is to be computed
if(spectra_calc .and. spectra_initialized ) call spectra_finalize()
$endif

return
end subroutine output_final

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_init ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine allocates the memory for arrays used for statistical
!  calculations
use param, only : dx,dy,dz,nx,ny,nz,lbz
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use param, only : tavg_calc
use grid_defs, only : grid
use functions, only : cell_indx
use stat_defs, only : xplane, yplane, zplane
use stat_defs, only : tavg, tavg_zplane
$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs
$endif
use stat_defs, only : type_set
implicit none

include 'tecryte.h'

character(120) :: var_list, fname
integer :: i,j,k

logical :: exst

real(rprec), pointer, dimension(:) :: x,y,z

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

if( tavg_calc ) then

  allocate(tavg(nx,ny,lbz:nz))
  allocate(tavg_zplane(nz))

  $if($OUTPUT_EXTRA)
  allocate(tavg_sgs(nx,ny,nz))
  $endif

  ! Initialize the derived types tavg and tavg_zplane  
  do k=1,Nz
    do j=1, Ny
      do i=1, Nx
        call type_set( tavg(i,j,k), 0._rprec )
        $if($OUTPUT_EXTRA)
        call type_set( tavg_sgs(i,j,k), 0._rprec )
        $endif        
      enddo
    enddo

    call type_set( tavg_zplane(k), 0._rprec )

  enddo

endif

! Initialize information for x-planar stats/data
if(xplane_calc) then

  allocate(xplane(xplane_nloc))

  xplane(:) % istart = -1
  xplane(:) % ldiff = 0.
  
!  Compute istart and ldiff
  do i=1,xplane_nloc
    xplane(i) % istart = cell_indx('i', dx, xplane_loc(i))
    xplane(i) % ldiff = xplane_loc(i) - x(xplane(i) % istart)   
  enddo
    
endif

! Initialize information for y-planar stats/data
if(yplane_calc) then

  allocate(yplane(yplane_nloc))

  yplane(:) % istart = -1
  yplane(:) % ldiff = 0.
  
!  Compute istart and ldiff
  do j=1,yplane_nloc
    yplane(j) % istart = cell_indx('j', dy, yplane_loc(j))
    yplane(j) % ldiff = yplane_loc(j) - y(yplane(j) % istart)
  enddo
    
endif

! Initialize information for z-planar stats/data
if(zplane_calc) then

  allocate(zplane(zplane_nloc))

!  Initialize 
  zplane(:) % istart = -1
  zplane(:) % ldiff = 0. 
  zplane(:) % coord=-1 
  
!  Compute istart and ldiff
  do k=1,zplane_nloc

    $if ($MPI)
    if(zplane_loc(k) >= z(1) .and. zplane_loc(k) < z(nz)) then
      zplane(k) % coord = coord
      zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
      zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
    endif
    $else
    zplane(k) % coord = 0
    zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
    zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
    $endif

  enddo  
  
endif

nullify(x,y,z)

return
end subroutine output_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Load tavg.out files
use messages
use stat_defs, only : tavg, tavg_total_time, tavg_dt, tavg_initialized
use stat_defs, only : operator(.MUL.)
$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs, tavg_total_time_sgs
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.tavg_init'
character (*), parameter :: ftavg_in = path // 'tavg.out'
$if($OUTPUT_EXTRA)
character (*), parameter :: ftavg_sgs_in = path // 'tavg_sgs.out'
$endif
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst
!integer :: i,j,k

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

fname = ftavg_in
$if ($MPI)
call string_concat( fname, MPI_suffix, coord )
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  !  Nothing to read in
  if (coord == 0) then
    write(*,*) ' '
    write(*,*)'No previous time averaged data - starting from scratch.'
  endif

    tavg_total_time = 0.0_rprec
    ! note: tavg was already initialized to zero in output_init routine
 
else

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
    read (1) tavg

    close(1)

endif

!------
$if($OUTPUT_EXTRA)

    fname = ftavg_sgs_in
    $if ($MPI)
    call string_concat( fname, MPI_suffix, coord )
    $endif

    inquire (file=fname, exist=exst)
    if (.not. exst) then
        !  Nothing to read in
        if(coord == 0) then
            write(*,*) ' '
            write(*,*)'No previous time averaged data (sgs) - starting from scratch.'
        endif

        tavg_total_time_sgs = 0._rprec  
    else
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

        read (1) tavg_total_time_sgs
        read (1) tavg_sgs

        close(1)    
    endif
    
$endif

! Initialize tavg_dt
tavg_dt = 0.0_rprec

! Set global switch that tavg as been initialized
tavg_initialized = .true. 

return
end subroutine tavg_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_compute()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine collects the stats for each flow 
!  variable quantity
use stat_defs, only : tavg,tavg_total_time,tavg_dt
$if($OUTPUT_EXTRA)
use param, only : sgs_model
use stat_defs, only : tavg_sgs, tavg_total_time_sgs
use sgs_param
$endif
use param, only : nx,ny,nz,lbz,jzmax
use sim_param, only : u,v,w, dudz, dvdz, txx, txy, tyy, txz, tyz, tzz
$if($TURBINES)
use sim_param, only : fxa
$endif
$if($LVLSET)
use sim_param, only : fx, fy, fz, fxa, fya, fza
$endif
use functions, only : interp_to_uv_grid

implicit none

integer :: i,j,k

real(rprec) :: u_p, v_p, w_p, w_p2
real(rprec), allocatable, dimension(:,:,:) :: w_uv
allocate(w_uv(nx,ny,lbz:nz))
w_uv(1:nx,1:ny,lbz:nz)= interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz )

$if($MPI)
k=0
do j=1,ny
   do i=1,nx

      u_p = u(i,j,k)
      v_p = v(i,j,k) 
      w_p = w(i,j,k)
      w_p2= w_uv(i,j,k)

      tavg(i,j,k)%u = tavg(i,j,k)%u + u_p * tavg_dt                    
      tavg(i,j,k)%v = tavg(i,j,k)%v + v_p * tavg_dt                         
      tavg(i,j,k)%w = tavg(i,j,k)%w + w_p * tavg_dt

      tavg(i,j,k)%u2 = tavg(i,j,k)%u2 + u_p * u_p * tavg_dt
      tavg(i,j,k)%v2 = tavg(i,j,k)%v2 + v_p * v_p * tavg_dt
      tavg(i,j,k)%w2 = tavg(i,j,k)%w2 + w_p * w_p * tavg_dt
      tavg(i,j,k)%uv = tavg(i,j,k)%uv + u_p * v_p * tavg_dt
      tavg(i,j,k)%uw = tavg(i,j,k)%uw + u_p * w_p2 * tavg_dt
      tavg(i,j,k)%vw = tavg(i,j,k)%vw + v_p * w_p2 * tavg_dt
$if ($TURBINES )
      tavg(i,j,k)%fx = tavg(i,j,k)%fx + fxa(i,j,k) * tavg_dt
$endif
      tavg(i,j,k)%cs_opt2 = tavg(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * tavg_dt

   enddo
enddo
$endif

do k=1,jzmax  
  do j=1,ny
    do i=1,nx
   
      u_p = u(i,j,k)
      v_p = v(i,j,k) 
      w_p = w(i,j,k)
      w_p2= w_uv(i,j,k)
    
      tavg(i,j,k)%u = tavg(i,j,k)%u + u_p * tavg_dt                    
      tavg(i,j,k)%v = tavg(i,j,k)%v + v_p * tavg_dt                         
      tavg(i,j,k)%w = tavg(i,j,k)%w + w_p * tavg_dt

      tavg(i,j,k)%u2 = tavg(i,j,k)%u2 + u_p * u_p * tavg_dt
      tavg(i,j,k)%v2 = tavg(i,j,k)%v2 + v_p * v_p * tavg_dt
      tavg(i,j,k)%w2 = tavg(i,j,k)%w2 + w_p * w_p * tavg_dt
      tavg(i,j,k)%uv = tavg(i,j,k)%uv + u_p * v_p * tavg_dt
      tavg(i,j,k)%uw = tavg(i,j,k)%uw + u_p * w_p2 * tavg_dt
      tavg(i,j,k)%vw = tavg(i,j,k)%vw + v_p * w_p2 * tavg_dt
$if ($TURBINES )      
      tavg(i,j,k)%fx = tavg(i,j,k)%fx + (             fxa(i,j,k)) * tavg_dt 
$endif
      tavg(i,j,k)%cs_opt2 = tavg(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * tavg_dt
      
    enddo
  enddo
enddo

$if ($OUTPUT_EXTRA)
do k=1,jzmax       
  do j=1,ny
    do i=1,nx
       ! === w-grid variables ===
       tavg_sgs(i,j,k)%Nu_t = tavg_sgs(i,j,k)%Nu_t + Nu_t(i,j,k) * tavg_dt
    enddo
 enddo
enddo

$endif

deallocate(w_uv)

! Update tavg_total_time for variable time stepping
tavg_total_time = tavg_total_time + tavg_dt
$if($OUTPUT_EXTRA)
tavg_total_time_sgs = tavg_total_time_sgs + tavg_dt
$endif

! Set tavg_dt back to zero for next increment
tavg_dt = 0.0_rprec

return

end subroutine tavg_compute

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_finalize()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use grid_defs, only : grid
use stat_defs, only : tavg_t, tavg_total_time, tavg
use stat_defs, only : rs_t, rs 
use stat_defs, only : operator(.DIV.), operator(.MUL.)
use stat_defs, only : operator(.ADD.), operator(.SUB.)
use stat_defs, only : rs_compute,tavg_interp_to_uv_grid
$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs, tavg_total_time_sgs
$endif
use param, only : ny,nz,nz_tot
$if($MPI)
use mpi_defs, only : mpi_sync_real_array,MPI_SYNC_DOWNUP
use param, only : ierr,comm
$endif

implicit none
character(64) :: fname_velb, &
     fname_vel2b, fname_ddzb, &
     fname_taub, fname_fb, &
     fname_rsb, fname_csb
$if($OUTPUT_EXTRA)  
character(64) :: fname_sgs_TnNu
$endif  
integer :: i,j,k

allocate(rs(nx,ny,lbz:nz))

fname_velb = path // 'output/binary_vel_avg.dat'
fname_vel2b = path // 'output/binary_vel2_avg.dat'
fname_ddzb = path // 'output/binary_ddz_avg.dat'
fname_taub = path // 'output/binary_tau_avg.dat'
fname_fb = path // 'output/binary_force_avg.dat'
fname_rsb = path // 'output/binary_rs.dat'
fname_csb = path // 'output/binary_cs_opt2.dat'

$if($OUTPUT_EXTRA)  
fname_sgs_TnNu = path // 'output/TnNu_avg.dat'
$endif  
  
$if ($MPI)
  !  For MPI implementation     
  call string_concat( fname_velb, '.c', coord)
  call string_concat( fname_vel2b, '.c', coord)
  call string_concat( fname_ddzb, '.c', coord)
  call string_concat( fname_taub, '.c', coord)
  call string_concat( fname_fb, '.c', coord)
  call string_concat( fname_rsb, '.c', coord)
  call string_concat( fname_csb, '.c', coord)
  
  $if($OUTPUT_EXTRA)  
  call string_concat( fname_sgs_TnNu, '.c', coord)
  $endif    
$endif

! Final checkpoint all restart data
call tavg_checkpoint()

$if($MPI)
call mpi_barrier( comm, ierr )
$endif

!  Perform time averaging operation
!  tavg = tavg / tavg_total_time
do k=jzmin,jzmax
  do j=1, Ny
    do i=1, Nx
      tavg(i,j,k) = tavg(i,j,k) .DIV. tavg_total_time
    enddo
  enddo
enddo

$if ($OUTPUT_EXTRA)
do k=1,jzmax
  do j=1, Ny
    do i=1, Nx
      tavg_sgs(i,j,k) = tavg_sgs(i,j,k) .DIV. tavg_total_time_sgs
    enddo
  enddo
enddo
$endif

!  Sync entire tavg structure
$if($MPI)
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%v, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%w, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%u2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%v2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%w2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%uw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%vw, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%uv, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%fx, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( tavg(1:nx,1:ny,lbz:nz)%cs_opt2, 0, MPI_SYNC_DOWNUP )
$if ($OUTPUT_EXTRA)
call mpi_sync_real_array( tavg_sgs(1:nx,1:ny,lbz:nz)%Nu_t, 0, MPI_SYNC_DOWNUP )
$endif
$endif

! ----- Write all the 3D data -----
open(unit=13,file=fname_velb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec/4)
write(13,rec=1) tavg(:nx,:ny,1:nz)%u
write(13,rec=2) tavg(:nx,:ny,1:nz)%v
write(13,rec=3) tavg(:nx,:ny,1:nz)%w
close(13)

open(unit=13,file=fname_vel2b,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec/4)
write(13,rec=1) tavg(:nx,:ny,1:nz)%u2
write(13,rec=2) tavg(:nx,:ny,1:nz)%v2
write(13,rec=3) tavg(:nx,:ny,1:nz)%w2
write(13,rec=4) tavg(:nx,:ny,1:nz)%uw
write(13,rec=5) tavg(:nx,:ny,1:nz)%vw
write(13,rec=6) tavg(:nx,:ny,1:nz)%uv
close(13)

open(unit=13,file=fname_fb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec/4)
write(13,rec=1) tavg(:nx,:ny,1:nz)%fx
close(13)

open(unit=13,file=fname_csb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec/4)
write(13,rec=1) tavg(:nx,:ny,1:nz)%cs_opt2 
close(13)

$if($OUTPUT_EXTRA)
open(unit=13,file=fname_sgs_TnNu,form='unformatted',convert='big_endian',access='direct',recl=nx*ny*nz*rprec/4)
write(13,rec=1) tavg_sgs(:nx,:ny,1:nz)%Nu_t
close(13)
$endif

! Do the Reynolds stress calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
tavg = tavg_interp_to_uv_grid( tavg )

rs = rs_compute( tavg , lbz)

open(unit=13,file=fname_rsb,form='unformatted',convert='big_endian',access='direct',recl=nx*ny*nz*rprec/4)
write(13,rec=1) rs(:nx,:ny,1:nz)%up2
write(13,rec=2) rs(:nx,:ny,1:nz)%vp2
write(13,rec=3) rs(:nx,:ny,1:nz)%wp2
write(13,rec=4) rs(:nx,:ny,1:nz)%upwp
write(13,rec=5) rs(:nx,:ny,1:nz)%vpwp
write(13,rec=6) rs(:nx,:ny,1:nz)%upvp
close(13)

deallocate(rs)

$if($MPI)
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
$endif

end subroutine tavg_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_checkpoint()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'tavg_finalize' at the end of the
! simulation.
!
use param, only : checkpoint_tavg_file
use stat_defs, only : tavg_total_time, tavg
$if($OUTPUT_EXTRA)
use param, only : checkpoint_tavg_sgs_file
use stat_defs, only : tavg_total_time_sgs, tavg_sgs
$endif
implicit none

character(*), parameter :: sub_name = mod_name // '.tavg_checkpoint'

character(64) :: fname
logical :: opn

fname = checkpoint_tavg_file
$if($MPI)
call string_concat( fname, '.c', coord)
$endif

!  Write data to tavg.out
inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, action='write', position='rewind', &
  form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, action='write', position='rewind', &
  form='unformatted', convert='little_endian')
$else
open (1, file=fname, action='write', position='rewind', form='unformatted')
$endif

! write the entire structures
write (1) tavg_total_time
write (1) tavg
close(1)

!----
$if($OUTPUT_EXTRA)
fname = checkpoint_tavg_sgs_file
  $if($MPI)
  call string_concat( fname, '.c', coord)
  $endif

  !  Write data to tavg_sgs.out
  $if ($WRITE_BIG_ENDIAN)
  open (1, file=fname, action='write', position='rewind', &
       form='unformatted', convert='big_endian')
  $elseif ($WRITE_LITTLE_ENDIAN)
  open (1, file=fname, action='write', position='rewind', &
       form='unformatted', convert='little_endian')
  $else
  open (1, file=fname, action='write', position='rewind', form='unformatted')
  $endif

  ! write the entire structures
  write (1) tavg_total_time_sgs
  write (1) tavg_sgs
  close(1)
$endif

return
end subroutine tavg_checkpoint

end module io
