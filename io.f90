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
$if ($MPI)
!include 'mpif.h'
$endif
save
private

public jt_total, openfiles, closefiles, energy, output_loop, output_final,output_init

character (*), parameter :: mod_name = 'io'

! Output file id's (see README for assigned values)
!~ integer :: ke_fid

! Where to start start and end with nz index.
! For coord==nproc-1 nz_end=0 else  nz_end=1
! For coord==0       nz_st=0 else   nz_st=1
integer :: nz_st, nz_end

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine openfiles()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : use_cfl_dt, dt, cfl_f
implicit none
!~ include 'tecryte.h'

logical :: exst

! Temporary values used to read time step and CFL from file
real(rprec) :: dt_r, cfl_r

! Open output files using tecryte library
! Kinetic energy (check_ke.dat)
!ke_fid = open_file( path // 'output/check_ke.dat', 'append', 'formatted' )

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then

    open (1, file=fcumulative_time)   
    read(1, *) jt_total, total_time, total_time_dim, dt_r, cfl_r
    close (1)
    
  else  !--assume this is the first run on cumulative time
    if( coord == 0 ) then
      write (*, *) '--> Assuming jt_total = 0, total_time = 0.0'
    endif
    jt_total = 0
    total_time = 0.0_rprec
    total_time_dim = 0.0_rprec
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
subroutine closefiles()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine is used to close all open files used by the main
! program. These files are opened by calling 'openfiles'.
implicit none

! Close kinetic energy file
!close( ke_fid )

return
end subroutine closefiles

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine energy (ke)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types,only:rprec
use param
use sim_param,only:u,v,w
use messages
$if ($XLF)
  use ieee_arithmetic  !--for NAN checking
$endif
implicit none

!~ include 'tecryte.h'

$if ($DEBUG)
character (*), parameter :: sub_name = 'energy'


integer, parameter :: NAN_MAX = 10
                      !--write this many NaN's before calling error (to aid
                      !  diagnosis of problem)
logical, parameter :: DEBUG = .true.
$endif
integer :: jx, jy, jz, nan_count

$if($DEBUG)
logical :: nan
$endif

real(kind=rprec)::KE,temp_w
$if ($MPI)
  real (rprec) :: ke_global
$endif

! Initialize variables
nan_count = 0
ke=0._rprec

z_loop: do jz=1,nz-1
    do jy=1,ny
        do jx=1,nx
            
            temp_w = 0.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
            ke = ke + (u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)
            
            $if ($DEBUG)
            if (DEBUG) then
                $if ($IFORT || $IFC)
                    nan = isnan (ke)
                $elsif ($XLF)
                    !--this is a bit verbose, should make into sub-program
                    if (ieee_support_datatype (ke)) then
                        if (ieee_support_nan (ke)) nan = ieee_is_nan (ke)
                    end if
                $endif
 
                if (nan) then
                    nan_count = nan_count + 1
                    write (*, *) 'NaN in ke at (jx, jy, jz) =', jx, jy, jz
                    write (*, *) 'jt_total = ', jt_total
                    write (*, *) 'u = ', u(jx, jy, jz)
                    write (*, *) 'v = ', v(jx, jy, jz)
                    write (*, *) 'w = ', w(jx, jy, jz)
                    if ( nan_count >= NAN_MAX ) exit z_loop
                end if
            end if
            $endif 
            
        end do
    end do
end do z_loop

! Perform spatial averaging
ke = ke*0.5_rprec/(nx*ny*(nz-1))

! Check if NaN's where found
$if ($DEBUG)
if ( nan_count > 0 ) call error (sub_name, 'NaN found')
$endif

!~ $if ($MPI)
!~   call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
!~   if (rank == 0) then  !--note its rank here, not coord
!~     ke = ke_global/nproc
!~     call write_real_data( ke_fid, 'formatted', 2, (/ total_time, ke /))
!~   end if
!~ $else
!~   call write_real_data( ke_fid, 'formatted', 2, (/ total_time, ke /))
!~ $endif

! Open KE file 
open(2,file='output/check_ke.dat',status='unknown',form='formatted',           &
     position='append')

$if ($MPI)

    call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
    if (rank == 0) then  !--note its rank here, not coord
        ke = ke_global/nproc
         write(2,*) total_time,ke
    end if

$else

write(2,*) total_time,ke

$endif

close(2)

end subroutine energy

$if($CGNS)
$if ($MPI)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine write_parallel_cgns ( file_name, nx, ny, nz, nz_tot, start_n,       &
                                     end_n, xin, yin, zin, num_fields,         &
                                     fieldNames, input )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Grid information
use grid_defs, only : grid

! Load variables to write
use param, only : total_time !,output_velocity,output_pressure

implicit none

! This subroutine writes parallel CGNS file output
include 'cgnslib_f.h'

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
character(*), intent(in) :: file_name  ! Name of file to be written
character(*), intent(in), dimension(:) :: fieldNames ! Name of fields we are writting
real(rprec), intent(in), dimension(:) :: input ! Data to be written
real(rprec), intent(in), dimension(:) :: xin,yin,zin ! Coordinates to write
integer, intent(in) :: start_n(:)  ! Where the total node counter starts nodes
integer, intent(in) :: end_n(:)  ! Where the total node counter ends nodes

integer :: fn          ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base        ! base number
integer :: zone        ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1        ! solution number
integer :: field     ! section number
integer :: sizes(3,3)    ! Sizes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: x,y,z

! The total number of nodes in this processor
nnodes=nx*ny*nz

! Create grid points
do k=1,nz
    do j=1,ny
        do i=1,nx
            x(i,j,k) = xin(i)
            y(i,j,k) = yin(j)
            z(i,j,k) = zin(k)
        enddo
    enddo
enddo

! Sizes, used to create zone
sizes(:,1) = (/nx,ny,nz_tot/)
sizes(:,2) = (/nx-1,ny-1,nz_tot-1/)
sizes(:,3) = (/0 , 0, 0/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
endif

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX',              &
                      (/nx,ny,nz/), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY',              &
                      (/nx,ny,nz/), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ',              &
                      (/nx,ny,nz/), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the coordinate data in parallel to the queue
call cgp_queue_set_f(1, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f
 
! This is done for the 3 dimensions x,y and z
! It writes the coordinates
call cgp_coord_write_data_f(fn, base, zone, 1,   &
                            start_n, end_n, x(1:nx,1:ny,1:nz), ier)    
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_data_f(fn, base, zone, 2,   &
                            start_n, end_n, y(1:nx,1:ny,1:nz), ier)   
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_data_f(fn, base, zone, 3,   &
                            start_n, end_n, z(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f
    
! Write out the queued coordinate data
call cgp_queue_flush_f(ier)
if (ier .ne. CG_OK) call cgp_error_exit_f
call cgp_queue_set_f(0, ier)

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
                           field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
                                input((i-1)*nnodes+1:(i)*nnodes), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

enddo

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

end subroutine write_parallel_cgns
$endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine write_serial_cgns ( file_name, nx, ny, nz, xin, yin, zin, num_fields,  &
                                     fieldNames, input )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Grid information
use grid_defs, only : grid

! Load variables to write
use param, only : total_time !,output_velocity,output_pressure

implicit none

! This subroutine writes parallel CGNS file output
include 'cgnslib_f.h'

integer, intent(in) :: nx, ny, nz, num_fields
character(*), intent(in) :: file_name  ! Name of file to be written
character(*), intent(in), dimension(:) :: fieldNames ! Name of fields we are writting
real(rprec), intent(in), dimension(:) :: input ! Data to be written
real(rprec), intent(in), dimension(:) :: xin,yin,zin ! Coordinates to write

integer :: fn          ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base        ! base number
integer :: zone        ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1        ! solution number
integer :: field     ! section number
integer :: sizes(3,3)    ! Sizes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: x,y,z

! The total number of nodes in this processor
nnodes=nx*ny*nz

! Create grid points
do k=1,nz
    do j=1,ny
        do i=1,nx
            x(i,j,k) = xin(i)
            y(i,j,k) = yin(j)
            z(i,j,k) = zin(k)
        enddo
    enddo
enddo

! Sizes, used to create zone
sizes(:,1) = (/nx,ny,nz/)
sizes(:,2) = (/nx-1,ny-1,nz-1/)
sizes(:,3) = (/0 , 0, 0/)

! Open CGNS file
call cg_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
endif

! Create data nodes for coordinates
call cg_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX',               &
                      x(1:nx,1:ny,1:nz), (/nx,ny,nz/), ier)
if (ier .ne. CG_OK) call cg_error_exit_f

call cg_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY',               &
                      y(1:nx,1:ny,1:nz), (/nx,ny,nz/), ier)
if (ier .ne. CG_OK) call cg_error_exit_f

call cg_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ',               &
                      z(1:nx,1:ny,1:nz), (/nx,ny,nz/), ier)
if (ier .ne. CG_OK) call cg_error_exit_f

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cg_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
                           input((i-1)*nnodes+1:(i)*nnodes), field, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

!~     call cg_field_write_data_f(fn, base, zone, sol, field,  &
!~                                 input((i-1)*nnodes+1:(i)*nnodes), ier)
!~     if (ier .ne. CG_OK) call cg_error_exit_f

enddo

! Close the file
call cg_close_f(fn, ier)
if (ier .ne. CG_OK) call cg_error_exit_f

end subroutine write_serial_cgns
! END OF  CGNS  PART
$endif 

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
$if (not $FFTW3)
use param, only : spectra_calc, spectra_nstart, spectra_nend, spectra_nskip
$endif
use param, only : point_calc, point_nstart, point_nend, point_nskip
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
     

!  Determine if spectra are to be calculated
$if (not $FFTW3)
if (spectra_calc) then

  ! Are we between the start and stop timesteps?
  if ((jt_total >= spectra_nstart).and.(jt_total <= spectra_nend)) then

    ! Every timestep (between nstart and nend), add to spectra_dt
    spectra_dt = spectra_dt + dt

    ! Are we at the beginning or a multiple of nstart?
    if ( mod(jt_total-spectra_nstart,spectra_nskip)==0 ) then

      ! Check if we have initialized spectra
      if (.not.spectra_initialized) then
        if (coord == 0) then
         write(*,*) '-------------------------------'
         write(*,"(1a,i9,1a,i9)") 'Starting running spectra calculations from ', spectra_nstart, ' to ', spectra_nend
         write(*,*) '-------------------------------'
        endif  ! coord==0
      
        call spectra_init()
      else
        call spectra_compute ()
      endif  ! init
     
    endif  ! mod of nskip

  endif  ! between nstart and nend

endif  ! spectra_calc
$endif

!  Determine if instantaneous point velocities are to be recorded
if(point_calc) then
  if(jt_total >= point_nstart .and. jt_total <= point_nend .and. ( mod(jt_total-point_nstart,point_nskip)==0) ) then
    if(jt_total == point_nstart) then
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
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use param, only : dx,dy,dz
use grid_defs, only : grid
use sim_param, only : u,v,w
$if($DEBUG)
use sim_param, only : p, dpdx, dpdy, dpdz,RHSx, RHSy, RHSz
$endif
use stat_defs, only : xplane, yplane, zplane, point
$if($MPI)
use param, only :ny,nz,comm,ierr
$endif
$if($LVLSET)
use level_set_base, only : phi
use sim_param, only : fx,fy,fz,fxa,fya,fza
$endif
implicit none

!include 'tecryte.h'      

integer, intent(IN) :: itype

$if($VERBOSE)
character (*), parameter :: sub_name = mod_name // '.inst_write'
$endif

character (64) :: fname
integer :: n, i, j, k
$if(not $BINARY)
character (64) :: var_list
integer :: nvars
$endif

$if($CGNS)
character (64) :: fname_cgns ! Name for CGNS output file
$endif

real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi,w_uv

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
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

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
    if(point(n) % coord == coord) then
    $endif

    ! Want to replace with write based on fid
!    call write_real_data(point(n) % fid, 'formatted', 4, (/ total_time, &
!         trilinear_interp(u(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz), &
!         trilinear_interp(v(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz), &
!         trilinear_interp(w_uv(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz) /))
    

    $if ($MPI)
    endif
    $endif

  enddo
!  Instantaneous write for entire domain
elseif(itype==2) then

  !////////////////////////////////////////////
  !/// WRITE VELOCITY                       ///
  !////////////////////////////////////////////

  $if( $BINARY )
  call string_splice( fname, path // 'output/binary_vel.', jt_total,'.dat')
  $else
!~   call string_splice( fname, path // 'output/vel.', jt_total, '.dat')
  $endif

  $if ($MPI)
      call string_concat( fname, '.c', coord )
  $endif

  ! Write CGNS Output
  $if ($CGNS and $MPI)
      call string_splice( fname_cgns, path //'output/output_', jt_total,'.cgns')

      call write_parallel_cgns(fname_cgns,nx,ny, nz - nz_end, nz_tot,           &
      (/ 1, 1,   (nz-1)*coord + 1 /),                                          &
      (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                              &
      x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                    &
      4, (/ 'VelocityX', 'VelocityY', 'VelocityZ', 'Pressure ' /),             &
      (/ u(1:nx,1:ny,1:(nz-nz_end)), v(1:nx,1:ny,1:(nz-nz_end)),                 &
         w_uv(1:nx,1:ny,1:(nz-nz_end)), u(1:nx,1:ny,1:(nz-nz_end)) /) )
  $endif
  
  $if($BINARY)
  open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
  write(13,rec=1) u(:nx,:ny,1:nz)
  write(13,rec=2) v(:nx,:ny,1:nz)
  write(13,rec=3) w_uv(:nx,:ny,1:nz)
  close(13)
  $else

    $if($LVLSET)
    var_list = '"x", "y", "z", "u", "v", "w", "phi"'
    nvars = 7
!    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
!    call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
!         (/ u(1:nx,1:ny,1:nz), &
!         v(1:nx,1:ny,1:nz), &
!         w_uv(1:nx,1:ny,1:nz), &
!         phi(1:nx,1:ny,1:nz)/), & 
!         4, x, y, z(1:nz))

    $else   
    var_list = '"x", "y", "z", "u", "v", "w"'
    nvars = 6
!    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
!    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
!         (/ u(1:nx,1:ny,1:nz), &
!         v(1:nx,1:ny,1:nz), &
!         w_uv(1:nx,1:ny,1:nz) /), &
!         4, x, y, z(1:nz))
    $endif
  
  $endif
  

  $if($MPI)
    ! Ensure that all processes finish before attempting to write 
    ! additional files. Otherwise it may flood the system with 
    ! too many I/O requests and crash the process 
    call mpi_barrier( comm, ierr )
  $endif

  !  Output instantaneous force field 
  $if(not $BINARY)
  $if($LVLSET)
    !////////////////////////////////////////////
    !/// WRITE FORCES                         ///
    !////////////////////////////////////////////

    ! Compute the total forces 
    call force_tot()

    !  Open file which to write global data
!~     call string_splice( fname, path // 'output/force.', jt_total, '.dat')

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

!~     var_list = '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>", "phi"'
!~     nvars = 7

!    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
!    call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
!                            (/ fx_tot, fy_tot, fz_tot, &
!                            phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))

    deallocate(fx_tot, fy_tot, fz_tot)

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

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
    $if( $BINARY )
    call string_splice( fname, path // 'output/binary_divvel.', jt_total,'.dat')
    $else
!~     call string_splice( fname, path // 'output/divvel.', jt_total, '.dat')
    $endif

    $if ($MPI)
      call string_concat( fname, '.c', coord )
    $endif

    $if($BINARY)
    open(unit=13,file=fname2,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
    write(13,rec=1) divvel(:nx,:ny,1:nz)
    close(13)
     
    $else
  
      $if($LVLSET)
!~       var_list = '"x", "y", "z", "divvel", "phi"'
!~       nvars = 5
!      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
!      call write_real_data_3D(fname, 'append', 'formatted', 2, nx, ny,nz, &
!           (/ divvel, phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
      $else   
!~       var_list = '"x", "y", "z", "divvel"'
!~       nvars = 4
!      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
!      call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny,nz, &
!           (/ divvel /), 4, x, y, z(1:nz))
      $endif

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

    $if($BINARY)
    call string_splice( fname, path // 'output/binary_pressure.', jt_total,'.dat')
    $else
    call string_splice( fname, path // 'output/pressure.', jt_total, '.dat')
    $endif

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

    call pressure_sync()

    $if($BINARY)
    open(unit=13,file=fname2,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
    write(13,rec=1) p(:nx,:ny,1:nz)
    write(13,rec=2) dpdx(:nx,:ny,1:nz)
    write(13,rec=3) dpdy(:nx,:ny,1:nz)
    write(13,rec=4) interp_to_uv_grid(dpdz(:nx,:ny,1:nz),1)
    close(13)

    $else
    
      $if($LVLSET)

      var_list = '"x", "y", "z", "p", "dpdx", "dpdy", "dpdz", "phi"'
      nvars = 8
!      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
!      call write_real_data_3D(fname, 'append', 'formatted', 5, nx, ny,nz, &
!           (/ p(1:nx,1:ny,1:nz), &
!           dpdx(1:nx,1:ny,1:nz), &
!           dpdy(1:nx,1:ny,1:nz), &
!           interp_to_uv_grid(dpdz(1:nx,1:ny,1:nz),1), &
!           phi(1:nx,1:ny,1:nz) /), &
!           4, x, y, z(1:nz))

      $else

      var_list = '"x", "y", "z", "p", "dpdx", "dpdy", "dpdz"'
      nvars = 7
!      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
!      call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
!        (/ p(1:nx,1:ny,1:nz), &
!        dpdx(1:nx,1:ny,1:nz), &
!        dpdy(1:nx,1:ny,1:nz), &
!        interp_to_uv_grid(dpdz(1:nx,1:ny,1:nz),1) /), &
!        4, x, y, z(1:nz))

      $endif

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
    $if($BINARY)
    call string_splice( fname, path // 'output/binary_RHS.', jt_total,'.dat')
    $else
    call string_splice( fname, path // 'output/RHS.', jt_total, '.dat')
    $endif

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif
  
    call RHS_sync()

    $if($BINARY)
    open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
    write(13,rec=1) RHSx(:nx,:ny,1:nz)
    write(13,rec=2) RHSy(:nx,:ny,1:nz)
    write(13,rec=3) dpdy(:nx,:ny,1:nz)
    write(13,rec=4) interp_to_uv_grid(RHSz(:nx,:ny,1:nz),1)
    close(13)
    
    $else

      $if($LVLSET)
      
      var_list = '"x", "y", "z", "RHSx", "RHSy", "RHSz", "phi"'
      nvars = 7
!      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
!      call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
!           (/ RHSx(1:nx,1:ny,1:nz), &
!           RHSy(1:nx,1:ny,1:nz), &
!           interp_to_uv_grid(RHSz(1:nx,1:ny,1:nz),1), &
!           phi(1:nx,1:ny,1:nz) /), & 
!           4, x, y, z(1:nz))

      $else

      var_list = '"x", "y", "z", "RHSx", "RHSy", "RHSz"'
      nvars = 6
!      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
!      call write_real_data_3D(fname, 'append', 'formatted', 3, nx, ny,nz, &
!           (/ RHSx(1:nx,1:ny,1:nz), &
!           RHSy(1:nx,1:ny,1:nz), &
!           interp_to_uv_grid(RHSz(1:nx,1:ny,1:nz),1) /), &
!           4, x, y, z(1:nz))

      $endif

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

     call string_splice( fname, path // 'output/vel.x-', xplane_loc(i), '.', jt_total, '.dat')


!    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny+1, Nz /), &
!      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4))  
  
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(u(xplane(i) % istart,j,k), &
             u(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
        vi(1,j,k) = linear_interp(v(xplane(i) % istart,j,k), &
             v(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
        wi(1,j,k) = linear_interp(w_uv(xplane(i) % istart,j,k), &
             w_uv(xplane(i) % istart+1,j,k), dx, &
             xplane(i) % ldiff)
      enddo
    enddo
  
     $if ($MPI)
         call string_concat( fname, '.c', coord )
     $endif

    $if ($CGNS and $MPI)    
        call string_splice( fname_cgns, path // 'output/plane_x_plane',        &
                            xplane_loc(i),'_', jt_total, '.cgns')

        call write_parallel_cgns (fname_cgns,1,ny, nz - nz_end, nz_tot,   &
                                    (/ 1, 1,   (nz-1)*coord + 1 /),            &
                                    (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /), &
                                xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),  &
                          3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),      &
                          (/ ui(1,1:ny,1:(nz-nz_end)), vi(1,1:ny,1:(nz-nz_end)), &
                             wi(1,1:ny,1:(nz-nz_end)) /) )
    $endif
     
!    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
!      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))     

    $if($LVLSET)

    call string_splice( fname, path // 'output/force.x-', xplane_loc(i), '.', jt_total, '.dat')

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

!    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny+1, Nz/), &
!      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', &
!      numtostr(coord,6), 2, real(total_time,4))

    !  Sum both induced forces, f{x,y,z}, and applied forces, f{x,y,z}a
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(fx_tot(xplane(i) % istart,j,k), &
             fx_tot(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)

        vi(1,j,k) = linear_interp(fy_tot(xplane(i) % istart,j,k), &
             fy_tot(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)

        wi(1,j,k) = linear_interp(fz_tot(xplane(i) % istart,j,k), &
             fz_tot(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)

      enddo
    enddo


!    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
!      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))




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
          F_LM_s(j,k) = linear_interp(F_LM_uv(xplane(i) % istart,j,k), &
               F_LM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          F_MM_s(j,k) = linear_interp(F_MM_uv(xplane(i) % istart,j,k), &
               F_MM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          beta_s(j,k) = linear_interp(beta_uv(xplane(i) % istart,j,k), &
               beta_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Cs_opt2_s(j,k) = linear_interp(Cs_opt2_uv(xplane(i) % istart,j,k), &
               Cs_opt2_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Nu_t_s(j,k) = linear_interp(Nu_t_uv(xplane(i) % istart,j,k), &
               Nu_t_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.x-', xplane_loc(i), '.', jt_total, '.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

!      call write_tecplot_header_ND(fname, 'rewind', 8, (/ 1, Ny+1, Nz/), &
!           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

!      call write_real_data_3D(fname, 'append', 'formatted', 5, 1,ny,nz, &
!                              (/ F_LM_s, F_MM_s, beta_s, Cs_opt2_s, Nu_t_s /), &
!                              2, (/ xplane_loc(i) /), y, z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

    elseif( sgs_model == 5 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(F_QN_s(nx,nz),F_NN_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do j=1,Ny
!
          F_LM_s(j,k) = linear_interp(F_LM_uv(xplane(i) % istart,j,k), &
               F_LM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          F_MM_s(j,k) = linear_interp(F_MM_uv(xplane(i) % istart,j,k), &
               F_MM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          F_QN_s(j,k) = linear_interp(F_QN_uv(xplane(i) % istart,j,k), &
               F_QN_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)  
          F_NN_s(j,k) = linear_interp(F_NN_uv(xplane(i) % istart,j,k), &
               F_NN_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)                                         
          beta_s(j,k) = linear_interp(beta_uv(xplane(i) % istart,j,k), &
               beta_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Cs_opt2_s(j,k) = linear_interp(Cs_opt2_uv(xplane(i) % istart,j,k), &
               Cs_opt2_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Nu_t_s(j,k) = linear_interp(Nu_t_uv(xplane(i) % istart,j,k), &
               Nu_t_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)

        enddo
      enddo

      call string_splice( fname , path // 'output/ldsm.x-', xplane_loc(i), '.', jt_total, '.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

!      call write_tecplot_header_ND(fname, 'rewind', 10, (/ 1, Ny+1, Nz/), &
!        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

!      call write_real_data_3D(fname, 'append', 'formatted', 7, 1,ny,nz, &
!                             (/ F_LM_s, F_MM_s, F_QN_s, F_NN_s, beta_s, Cs_opt2_s, Nu_t_s /), &
!                             2, (/ xplane_loc(i) /), y, z(1:nz)) 

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

!    call string_splice( fname, path // 'output/vel.y-', yplane_loc(j), '.', jt_total, '.dat')


!    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, 1, Nz/), &
!      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4)) 

    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(u(i,yplane(j) % istart,k), &
             u(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        vi(i,1,k) = linear_interp(v(i,yplane(j) % istart,k), &
             v(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        wi(i,1,k) = linear_interp(w_uv(i,yplane(j) % istart,k), &
             w_uv(i,yplane(j) % istart+1,k), dy, &
             yplane(j) % ldiff)
        
      enddo
    enddo
  
    $if ($MPI)
        call string_concat( fname, '.c', coord )
    $endif    

    $if ($CGNS and $MPI)    
        call string_splice( fname_cgns, path // 'output/plane_y_plane',        &
                            yplane_loc(j),'_', jt_total, '.cgns')

        call write_parallel_cgns (fname_cgns,nx,1, nz - nz_end, nz_tot,   &
                                    (/ 1, 1,   (nz-1)*coord + 1 /),            &
                                    (/ nx, 1, (nz-1)*(coord+1) + 1 - nz_end /), &
                                x(1:nx) , yplane_loc(j:j) , z(1:(nz-nz_end) ),  &
                          3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),      &
                          (/ ui(1:nx,1,1:(nz-nz_end)), vi(1:nx,1,1:(nz-nz_end)), &
                             wi(1:nx,1,1:(nz-nz_end)) /) )
    $endif

!    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
!      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))    
  
    $if($LVLSET)

    call string_splice( fname, path // 'output/force.y-', yplane_loc(j), '.', jt_total, '.dat')

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

!    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, 1, Nz/), &
!      '"x", "y", "z", "fx", "fy", "fz"', numtostr(coord,6), 2, real(total_time,4))  
  
    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(fx_tot(i,yplane(j) % istart,k), &
             fx_tot(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        vi(i,1,k) = linear_interp(fy_tot(i,yplane(j) % istart,k), &
             fy_tot(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        wi(i,1,k) = linear_interp(fz_tot(i,yplane(j) % istart,k), &
             fz_tot(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)

      enddo
    enddo
    
!    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
!      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))       
    
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
          F_LM_s(i,k) = linear_interp(F_LM_uv(i,yplane(j) % istart,k), &
               F_LM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_MM_s(i,k) = linear_interp(F_MM_uv(i,yplane(j) % istart,k), &
               F_MM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          beta_s(i,k) = linear_interp(beta_uv(i,yplane(j) % istart,k), &
               beta_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Cs_opt2_s(i,k) = linear_interp(Cs_opt2_uv(i,yplane(j) % istart,k), &
               Cs_opt2_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Nu_t_s(i,k) = linear_interp(Nu_t_uv(i,yplane(j) % istart,k), &
               Nu_t_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.y-', yplane_loc(j), '.', jt_total, '.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

!      call write_tecplot_header_ND(fname, 'rewind', 8, (/ Nx+1, 1, Nz/), &
!        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

!      call write_real_data_3D(fname, 'append', 'formatted', 5, nx,1,nz, &
!                              (/ F_LM_s, F_MM_s, beta_s, Cs_opt2_s, Nu_t_s /), &
!                              1, x, (/ yplane_loc(j) /), z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

    elseif( sgs_model == 5 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(F_QN_s(nx,nz),F_NN_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do i=1,Nx

          F_LM_s(i,k) = linear_interp(F_LM_uv(i,yplane(j) % istart,k), &
               F_LM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_MM_s(i,k) = linear_interp(F_MM_uv(i,yplane(j) % istart,k), &
               F_MM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_QN_s(i,k) = linear_interp(F_QN_uv(i,yplane(j) % istart,k), &
               F_QN_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_NN_s(i,k) = linear_interp(F_NN_uv(i,yplane(j) % istart,k), &
               F_NN_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)                                      
          beta_s(i,k) = linear_interp(beta_uv(i,yplane(j) % istart,k), &
               beta_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Cs_opt2_s(i,k) = linear_interp(Cs_opt2_uv(i,yplane(j) % istart,k), &
               Cs_opt2_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Nu_t_s(i,k) = linear_interp(Nu_t_uv(i,yplane(j) % istart,k), &
               Nu_t_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.y-', yplane_loc(j), '.', jt_total,'.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ',  "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

!      call write_tecplot_header_ND(fname, 'rewind', 10, (/ Nx+1, 1, Nz/), &
!        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

!      call write_real_data_3D(fname, 'append', 'formatted', 7, nx,1,nz, &
!                              (/ F_LM_s, F_MM_s, F_QN_s, F_NN_s, beta_s, Cs_opt2_s, Nu_t_s /), &
!                              1, x, (/ yplane_loc(j) /), z(1:nz)) 

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
    if(zplane(k) % coord == coord) then
    $endif

    $if ($BINARY)
    call string_splice( fname, path // 'output/binary_vel.z-', zplane_loc(k), '.', jt_total, '.dat')
    $else
    call string_splice( fname, path // 'output/vel.z-', zplane_loc(k), '.', jt_total, '.dat')    
    $endif
    
    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(u(i,j,zplane(k) % istart), &
             u(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)
        vi(i,j,1) = linear_interp(v(i,j,zplane(k) % istart), &
             v(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)
        wi(i,j,1) = linear_interp(w_uv(i,j,zplane(k) % istart), &
             w_uv(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)

      enddo
    enddo

    $if ($CGNS)
        call string_splice( fname_cgns, path // 'output/plane_z_plane',        &
                            zplane_loc(k),'_', jt_total, '.cgns')
!~         call write_serial_cgns(fname_cgns,nx,ny,1,x,y,zplane_loc(k:k), 3,      &
!~         (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                           &
!~         (/ ui(1,1:ny,1:nz),vi(1,1:ny,1:nz),wi(1,1:ny,1:nz) /) )

        call write_serial_cgns ( fname_cgns, nx, ny,1,x,y,zplane_loc(k:k), 3,   &
                               (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),    &
                    (/ ui(1:nx,1:ny,1),vi(1:nx,1:ny,1),wi(1:nx,1:ny,1) /) )
    $endif
    
    $if ($BINARY)
    open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*1*rprec)
    write(13,rec=1) ui(1:nx,1:ny,1)
    write(13,rec=2) vi(1:nx,1:ny,1)
    write(13,rec=3) wi(1:nx,1:ny,1)
    close(13)
    
    $else
!    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, Ny+1, 1/), &
!      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4)) 
!    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
!    (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /))   

    $endif
    
    $if($LVLSET)

    call string_splice( fname, path // 'output/force.z-', zplane_loc(k), '.', jt_total, '.dat')

!    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, Ny+1, 1/), &
!      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', &
!      numtostr(coord,6), 2, real(total_time,4))

    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(fx_tot(i,j,zplane(k) % istart), &
             fx_tot(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff) 
        vi(i,j,1) = linear_interp(fy_tot(i,j,zplane(k) % istart), &
             fy_tot(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff) 
        wi(i,j,1) = linear_interp(fz_tot(i,j,zplane(k) % istart), &
             fz_tot(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)

      enddo
    enddo
 
!    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
!         (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /) )      
    
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
          F_LM_s(i,j) = linear_interp(F_LM_uv(i,j,zplane(k) % istart), &
               F_LM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          F_MM_s(i,j) = linear_interp(F_MM_uv(i,j,zplane(k) % istart), &
               F_MM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          beta_s(i,j) = linear_interp(beta_uv(i,j,zplane(k) % istart), &
               beta_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          Cs_opt2_s(i,j) = linear_interp(Cs_opt2_uv(i,j,zplane(k) % istart), &
               Cs_opt2_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          Nu_t_s(i,j) = linear_interp(Nu_t_uv(i,j,zplane(k) % istart), &
               Nu_t_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.z-', zplane_loc(k), '.', jt_total, '.dat')

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

!      call write_tecplot_header_ND(fname, 'rewind', 8, (/ Nx+1, Ny+1, 1/), &
!        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))

!      call write_real_data_3D(fname, 'append', 'formatted', 5, nx,ny,1, &
!                              (/ F_LM_s, F_MM_s, beta_s, Cs_opt2_s, Nu_t_s /), &
!                              4, x, y, (/ zplane_loc(k) /) )

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

     elseif( sgs_model == 5 ) then
 
      allocate(F_LM_s(nx,ny),F_MM_s(nx,ny))
      allocate(F_QN_s(nx,ny),F_NN_s(nx,ny))
      allocate(beta_s(nx,ny),Cs_opt2_s(nx,ny))
      allocate(Nu_t_s(nx,ny))

      do j=1,Ny
        do i=1,Nx
!
          F_LM_s(i,j) = linear_interp(F_LM_uv(i,j,zplane(k) % istart), &
               F_LM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff) 
          F_MM_s(i,j) = linear_interp(F_MM_uv(i,j,zplane(k) % istart), &
               F_MM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff) 
          F_QN_s(i,j) = linear_interp(F_QN_uv(i,j,zplane(k) % istart), &
               F_QN_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          F_NN_s(i,j) = linear_interp(F_NN_uv(i,j,zplane(k) % istart), &
               F_NN_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          beta_s(i,j) = linear_interp(beta_uv(i,j,zplane(k) % istart), &
               beta_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)            
          Cs_opt2_s(i,j) = linear_interp(Cs_opt2_uv(i,j,zplane(k) % istart), &
               Cs_opt2_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)                                         
          Nu_t_s(i,j) = linear_interp(Nu_t_uv(i,j,zplane(k) % istart), &
               Nu_t_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)

        enddo
      enddo      

      call string_splice( fname, path // 'output/ldsm.z-', zplane_loc(k), '.', jt_total, '.dat')

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

!      call write_tecplot_header_ND(fname, 'rewind', 10, (/ Nx+1, Ny+1, 1/), &
!        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))

!      call write_real_data_3D(fname, 'append', 'formatted', 7, nx,ny,1, &
!                              (/ F_LM_s, F_MM_s, F_QN_s, F_NN_s, beta_s, Cs_opt2_s, Nu_t_s /), &
!                              4, x, y, (/ zplane_loc(k) /) )         

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

! Richard: Might not be necessary to do this as the function only seems to be called when LVLSET is activated
$if($TURBINES and not $LVLSET)
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = 0._rprec
fz_tot = 0._rprec
$elseif($LVLSET)
fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)
$else
fx_tot = 0._rprec
fy_tot = 0._rprec
fz_tot = 0._rprec
$endif

$if($MPI)
!  Sync forces
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
call mpi_sync_real_array( RHSx, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( RHSy, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( RHSz, 0 , MPI_SYNC_DOWN )
$endif

return
end subroutine RHS_sync
$endif

end subroutine inst_write

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkpoint ()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : nz, checkpoint_file, tavg_calc
$if(not $FFTW3)
use param, only : spectra_calc
use stat_defs, only : spectra_initialized
$endif
$if($MPI)
use param, only : comm,ierr
$endif
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
$if($DYN_TN)
use sgs_param, only: F_ee2, F_deedt2, ee_past
$endif
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
$if(not $FFTW3)
! Checkpoint spectra restart data
if( spectra_calc .and. spectra_initialized ) call spectra_checkpoint()
$endif

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
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
$if(not $FFTW3)
use param, only : spectra_calc, spectra_nloc, spectra_loc
$endif
use param, only : tavg_calc
use grid_defs, only : grid
use functions, only : cell_indx
use stat_defs, only : point, xplane, yplane, zplane
use stat_defs, only : tavg, tavg_zplane
$if(not $FFTW3)
use stat_defs, only : spectra
$endif
$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs
$endif
use stat_defs, only : type_set
implicit none

!include 'tecryte.h'
!logical :: exst
!character(120) :: var_list
character(120) :: fname

integer :: i,j,k

real(rprec), pointer, dimension(:) :: x,y,z

! This adds one more element to the last processor (which contains an extra one)
! Processor nproc-1 has data from 1:nz
! Rest of processors have data from 1:nz-1
$if $MPI
    if ( coord == nproc-1 ) then
        nz_end=0
    else
        nz_end=1
    endif
$else
    nz_end=0    
$endif

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

$if(not $FFTW3)
if(spectra_calc) then

  allocate(spectra(spectra_nloc))
!  Initialize 
  spectra(:) % istart = -1
  spectra(:) % ldiff = 0._rprec
  spectra(:) % coord = -1

!  Compute istart and ldiff
  do k=1,spectra_nloc

    !  Initialize sub-arrays
    allocate( spectra(k)%power(lh) )
    spectra(k) % power = 0._rprec

    $if ($MPI)
    if(spectra_loc(k) >= z(1) .and. spectra_loc(k) < z(nz)) then
      spectra(k) % coord = coord
      spectra(k) % istart = cell_indx('k',dz,spectra_loc(k))
      spectra(k) % ldiff = spectra_loc(k) - z(spectra(k) % istart)
    endif
    $else
    spectra(k) % coord = 0
    spectra(k) % istart = cell_indx('k',dz,spectra_loc(k))
    spectra(k) % ldiff = spectra_loc(k) - z(spectra(k) % istart)
    $endif

  enddo

endif
$endif

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


!  Open files for instantaneous writing
if(point_calc) then

  allocate(point(point_nloc))

  !  Intialize the coord values (-1 shouldn't be used as coord so initialize to this)
  point % coord=-1
  point % fid = -1

  do i=1,point_nloc
    !  Find the processor in which this point lives
    $if ($MPI)
    if(point_loc(i)%xyz(3) >= z(1) .and. point_loc(i)%xyz(3) < z(nz)) then
    $endif

    point(i) % coord = coord
  
    point(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
    point(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
    point(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))
    
    point(i) % xdiff = point_loc(i)%xyz(1) - x(point(i) % istart)
    point(i) % ydiff = point_loc(i)%xyz(2) - y(point(i) % jstart)
    point(i) % zdiff = point_loc(i)%xyz(3) - z(point(i) % kstart)
    
    call string_splice(fname, path // 'output/vel.x-', point_loc(i)%xyz(1), &
         '.y-', point_loc(i)%xyz(2), &
         '.z-', point_loc(i)%xyz(3),'.dat')
    
    !  Add tecplot header if file does not exist
!~     inquire (file=fname, exist=exst)
!~     if (exst) then
!~        point(i) % fid = open_file( fname, 'append', 'formatted' )
!~     else
!~ 
!~        point(i) % fid = open_file( fname, 'rewind', 'formatted' )
!~        var_list = '"t", "u", "v", "w"'
!~        ! Compilation error
!~ !       call write_tecplot_header_xyline( point(i) % fid, var_list )
!~ 
!~     endif
    
    $if($MPI)
    endif
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
use stat_defs, only : tavg_sgs, tavg_total_time_sgs,tavg_time_stamp
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
use functions, only : interp_to_w_grid

implicit none

integer :: i,j,k

real(rprec) :: u_p, v_p, w_p
real(rprec), allocatable, dimension(:,:,:) :: u_w, v_w

allocate(u_w(nx,ny,lbz:nz),v_w(nx,ny,lbz:nz))

!  Interpolate velocities to w-grid
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( v(1:nx,1:ny,lbz:nz), lbz )

$if($MPI)
k=0
do j=1,ny
   do i=1,nx

      u_p = u_w(i,j,k)
      v_p = v_w(i,j,k) 
      w_p = w(i,j,k)

      ! === w-grid variables === 
      tavg(i,j,k)%u = tavg(i,j,k)%u + u_p * tavg_dt                    
      tavg(i,j,k)%v = tavg(i,j,k)%v + v_p * tavg_dt                         
      tavg(i,j,k)%w = tavg(i,j,k)%w + w_p * tavg_dt

      tavg(i,j,k)%u2 = tavg(i,j,k)%u2 + u_p * u_p * tavg_dt
      tavg(i,j,k)%v2 = tavg(i,j,k)%v2 + v_p * v_p * tavg_dt
      tavg(i,j,k)%w2 = tavg(i,j,k)%w2 + w_p * w_p * tavg_dt
      tavg(i,j,k)%uv = tavg(i,j,k)%uv + u_p * v_p * tavg_dt
      tavg(i,j,k)%uw = tavg(i,j,k)%uw + u_p * w_p * tavg_dt
      tavg(i,j,k)%vw = tavg(i,j,k)%vw + v_p * w_p * tavg_dt

      tavg(i,j,k)%dudz = tavg(i,j,k)%dudz + dudz(i,j,k) * tavg_dt
      tavg(i,j,k)%dvdz = tavg(i,j,k)%dvdz + dvdz(i,j,k) * tavg_dt

      ! === uv-grid variables ===
      tavg(i,j,k)%txx = tavg(i,j,k)%txx + txx(i,j,k) * tavg_dt
      tavg(i,j,k)%tyy = tavg(i,j,k)%tyy + tyy(i,j,k) * tavg_dt
      tavg(i,j,k)%tzz = tavg(i,j,k)%tzz + tzz(i,j,k) * tavg_dt
      tavg(i,j,k)%txy = tavg(i,j,k)%txy + txy(i,j,k) * tavg_dt
      ! === w-grid variables === 
      tavg(i,j,k)%txz = tavg(i,j,k)%txz + txz(i,j,k) * tavg_dt
      tavg(i,j,k)%tyz = tavg(i,j,k)%tyz + tyz(i,j,k) * tavg_dt

   enddo
enddo
$endif

do k=1,jzmax  
  do j=1,ny
    do i=1,nx
   
      u_p = u_w(i,j,k)
      v_p = v_w(i,j,k) 
      w_p = w(i,j,k)
          
      ! === w-grid variables === 
      tavg(i,j,k)%u = tavg(i,j,k)%u + u_p * tavg_dt                    
      tavg(i,j,k)%v = tavg(i,j,k)%v + v_p * tavg_dt                         
      tavg(i,j,k)%w = tavg(i,j,k)%w + w_p * tavg_dt

      tavg(i,j,k)%u2 = tavg(i,j,k)%u2 + u_p * u_p * tavg_dt
      tavg(i,j,k)%v2 = tavg(i,j,k)%v2 + v_p * v_p * tavg_dt
      tavg(i,j,k)%w2 = tavg(i,j,k)%w2 + w_p * w_p * tavg_dt
      tavg(i,j,k)%uv = tavg(i,j,k)%uv + u_p * v_p * tavg_dt
      tavg(i,j,k)%uw = tavg(i,j,k)%uw + u_p * w_p * tavg_dt
      tavg(i,j,k)%vw = tavg(i,j,k)%vw + v_p * w_p * tavg_dt

      tavg(i,j,k)%dudz = tavg(i,j,k)%dudz + dudz(i,j,k) * tavg_dt
      tavg(i,j,k)%dvdz = tavg(i,j,k)%dvdz + dvdz(i,j,k) * tavg_dt

      ! === uv-grid variables ===
      tavg(i,j,k)%txx = tavg(i,j,k)%txx + txx(i,j,k) * tavg_dt
      tavg(i,j,k)%tyy = tavg(i,j,k)%tyy + tyy(i,j,k) * tavg_dt
      tavg(i,j,k)%tzz = tavg(i,j,k)%tzz + tzz(i,j,k) * tavg_dt
      tavg(i,j,k)%txy = tavg(i,j,k)%txy + txy(i,j,k) * tavg_dt
      ! === w-grid variables === 
      tavg(i,j,k)%txz = tavg(i,j,k)%txz + txz(i,j,k) * tavg_dt
      tavg(i,j,k)%tyz = tavg(i,j,k)%tyz + tyz(i,j,k) * tavg_dt

$if ($TURBINES and not $LVLSET)      
      ! Includes both induced (IBM) and applied (RNS, turbines, etc.) forces 
      ! === uv-grid variables === 
      tavg(i,j,k)%fx = tavg(i,j,k)%fx + (             fxa(i,j,k)) * tavg_dt 
!      tavg(i,j,k)%fy = tavg(i,j,k)%fy + (fy(i,j,k)             ) * tavg_dt
      ! === w-grid variables === 
!      tavg(i,j,k)%fz = tavg(i,j,k)%fz + (fz(i,j,k)             ) * tavg_dt
$elseif ($LVLSET )
      ! Includes both induced (IBM) and applied (RNS, turbines, etc.) forces 
      ! === uv-grid variables === 
      tavg(i,j,k)%fx = tavg(i,j,k)%fx + (fx(i,j,k) + fxa(i,j,k)) * tavg_dt 
      tavg(i,j,k)%fy = tavg(i,j,k)%fy + (fy(i,j,k) + fya(i,j,k)) * tavg_dt
      ! === w-grid variables === 
      tavg(i,j,k)%fz = tavg(i,j,k)%fz + (fz(i,j,k) + fza(i,j,k)) * tavg_dt
$else
      ! Includes both induced (IBM) and applied (RNS, turbines, etc.) forces 
      ! === uv-grid variables === 
!      tavg(i,j,k)%fx = tavg(i,j,k)%fx + (fx(i,j,k)             ) * tavg_dt 
!      tavg(i,j,k)%fy = tavg(i,j,k)%fy + (fy(i,j,k)             ) * tavg_dt
      ! === w-grid variables === 
!      tavg(i,j,k)%fz = tavg(i,j,k)%fz + (fz(i,j,k)             ) * tavg_dt
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

if (sgs_model.gt.1) then
   do k=1,jzmax       
      do j=1,ny
         do i=1,nx
            ! === w-grid variables ===
            tavg_sgs(i,j,k)%ee_now = tavg_sgs(i,j,k)%ee_now + ee_now(i,j,k) * tavg_dt
         enddo
      enddo
   enddo
endif

if (sgs_model.gt.3) then
   do k=1,jzmax       
      do j=1,ny
         do i=1,nx
            ! === w-grid variables ===
            tavg_sgs(i,j,k)%Tn = tavg_sgs(i,j,k)%Tn + Tn_all(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_LM = tavg_sgs(i,j,k)%F_LM + F_LM(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_MM = tavg_sgs(i,j,k)%F_MM + F_MM(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_QN = tavg_sgs(i,j,k)%F_QN + F_QN(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_NN = tavg_sgs(i,j,k)%F_NN + F_NN(i,j,k) * tavg_dt
            
            $if ($DYN_TN)
            tavg_sgs(i,j,k)%F_ee2 = tavg_sgs(i,j,k)%F_ee2 + F_ee2(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_deedt2 = tavg_sgs(i,j,k)%F_deedt2 + F_deedt2(i,j,k) * tavg_dt
            $endif         
         enddo
      enddo
   enddo
endif
$endif

deallocate( u_w, v_w )

$if($BINARY)
! Added to compute the average of u on the grid points and prevent the 
! linear interpolation, which is not exactly valid due to log profile in boundary layer
do k=lbz,jzmax  
  do j=1,ny
    do i=1,nx
    u_p = u(i,j,k)
    enddo
  enddo
enddo
$endif

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
use grid_defs, only : grid !x,y,z
use stat_defs, only : tavg_t, tavg_zplane, tavg_total_time, tavg
use stat_defs, only : rs_t, rs, rs_zplane, cnpy_zplane 
use stat_defs, only : operator(.DIV.), operator(.MUL.)
use stat_defs, only :  operator(.ADD.), operator(.SUB.)
use stat_defs, only : type_set, type_zero_bogus
use stat_defs, only : tavg_interp_to_uv_grid, tavg_interp_to_w_grid
use stat_defs, only : rs_compute, cnpy_tavg_mul
$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs, tavg_total_time_sgs
$endif
use param, only : ny,nz,nz_tot
$if($MPI)
use mpi_defs, only : mpi_sync_real_array
use param, only : MPI_RPREC,coord_of_rank,rank_of_coord,comm,ierr,down,up,status
use stat_defs, only : rs_t, tavg_t
$endif
$if($LVLSET)
use level_set_base, only : phi
$endif

implicit none

!~ include 'tecryte.h'

character (*), parameter :: sub_name = mod_name // '.tavg_finalize'

! For .dat output
character(64) :: fname_vel, &
     fname_vel2, fname_ddz, &
     fname_tau, fname_f, &
     fname_rs, fname_cs, fname_u_vel_grid
     
! For bindary output
character(64) :: fname_velb, &
     fname_vel2b, fname_ddzb, &
     fname_taub, fname_fb, &
     fname_rsb, fname_csb, fname_u_vel_gridb

$if($CGNS)
! For CGNS
character(64) :: fname_vel_cgns, &
     fname_vel2_cgns, fname_ddz_cgns, &
     fname_tau_cgns, fname_f_cgns, &
     fname_rs_cgns, fname_cs_cgns, fname_u_vel_grid_cgns
$endif
     
character(64) :: fname_vel_zplane, fname_vel2_zplane, &
  fname_ddz_zplane, fname_tau_zplane, fname_f_zplane, &
  fname_rs_zplane, fname_cs_zplane, fname_cnpy_zplane
$if($OUTPUT_EXTRA)  
character(64) :: fname_sgs_TnNu, fname_sgs_Fsub
character(64) :: fname_sgs_ee
$endif  

integer :: i,j,k

$if($MPI)

integer :: MPI_RS, MPI_CNPY, MPI_TAVG
integer :: rs_type(1), rs_block(1), rs_disp(1)
integer :: cnpy_type(1), cnpy_block(1), cnpy_disp(1)
integer :: tavg_type(1), tavg_block(1), tavg_disp(1)

! Definitions for reconstructing z-planar averaged data
integer :: sendsize
integer, allocatable, dimension(:) :: recvsize, recvstride
real(rprec), allocatable, dimension(:) :: z_tot, zw_tot
type(rs_t), allocatable, dimension(:) :: rs_zplane_tot
type(rs_t), allocatable, dimension(:) :: cnpy_zplane_tot
type(tavg_t), allocatable, dimension(:) :: tavg_zplane_tot

$endif

type(rs_t) :: cnpy_avg
type(tavg_t) :: tavg_avg

real(rprec) :: favg

$if($LVLSET)
real(rprec) :: fx_global, fy_global, fz_global
$endif

real(rprec), pointer, dimension(:) :: x,y,z,zw

favg = real(nx*ny,kind=rprec)

nullify(x,y,z,zw)

allocate(rs(nx,ny,lbz:nz), rs_zplane(nz))
allocate(cnpy_zplane(nz))

x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! All processors need not do this, but that is ok
!  Set file names

! .dat
fname_vel = path // 'output/vel_avg.dat'
fname_vel2 = path // 'output/vel2_avg.dat'
fname_ddz = path // 'output/ddz_avg.dat'
fname_tau = path // 'output/tau_avg.dat'
fname_f = path // 'output/force_avg.dat'
fname_rs = path // 'output/rs.dat'
fname_cs = path // 'output/cs_opt2.dat'
fname_u_vel_grid = path // 'output/u_grid_vel.dat'

$if($CGNS)
! CGNS
fname_vel_cgns = path // 'output/vel_avg.cgns'
fname_vel2_cgns = path // 'output/vel2_avg.cgns'
fname_ddz_cgns = path // 'output/ddz_avg.cgns'
fname_tau_cgns = path // 'output/tau_avg.cgns'
fname_f_cgns = path // 'output/force_avg.cgns'
fname_rs_cgns = path // 'output/rs.cgns'
fname_cs_cgns = path // 'output/cs_opt2.cgns'
fname_u_vel_grid_cgns = path // 'output/u_grid_vel.cgns'
$endif

! Binary
fname_velb = path // 'output/binary_vel_avg.dat'
fname_vel2b = path // 'output/binary_vel2_avg.dat'
fname_ddzb = path // 'output/binary_ddz_avg.dat'
fname_taub = path // 'output/binary_tau_avg.dat'
fname_fb = path // 'output/binary_force_avg.dat'
fname_rsb = path // 'output/binary_rs.dat'
fname_csb = path // 'output/binary_cs_opt2.dat'
fname_u_vel_gridb = path // 'output/binary_u_grid_vel.dat'

fname_vel_zplane = path // 'output/vel_zplane_avg.dat'
fname_vel2_zplane = path // 'output/vel2_zplane_avg.dat'
fname_ddz_zplane = path // 'output/ddz_zplane_avg.dat'
fname_tau_zplane = path // 'output/tau_zplane_avg.dat'
fname_f_zplane = path // 'output/force_zplane_avg.dat'
fname_rs_zplane = path // 'output/rs_zplane.dat'
fname_cnpy_zplane = path // 'output/cnpy_zplane.dat'
fname_cs_zplane = path // 'output/cs_opt2_zplane.dat'

$if($OUTPUT_EXTRA)  
fname_sgs_TnNu = path // 'output/TnNu_avg.dat'
fname_sgs_Fsub = path // 'output/Fsub_avg.dat'
fname_sgs_ee = path // 'output/ee_avg.dat'
$endif  
  
$if ($MPI)
  $if($BINARY)
  !  For MPI implementation     
  call string_concat( fname_velb, '.c', coord)
  call string_concat( fname_vel2b, '.c', coord)
  call string_concat( fname_ddzb, '.c', coord)
  call string_concat( fname_taub, '.c', coord)
  call string_concat( fname_fb, '.c', coord)
  call string_concat( fname_rsb, '.c', coord)
  call string_concat( fname_csb, '.c', coord)
  call string_concat( fname_u_vel_gridb, '.c',coord)
  $else
  call string_concat( fname_vel, '.c', coord)
  call string_concat( fname_vel2, '.c', coord)
  call string_concat( fname_ddz, '.c', coord)
  call string_concat( fname_tau, '.c', coord)
  call string_concat( fname_f, '.c', coord)
  call string_concat( fname_rs, '.c', coord)
  call string_concat( fname_cs, '.c', coord)
  call string_concat( fname_u_vel_grid, '.c', coord)
  $endif
  
  $if($OUTPUT_EXTRA)  
  call string_concat( fname_sgs_TnNu, '.c', coord)
  call string_concat( fname_sgs_Fsub, '.c', coord)
  call string_concat( fname_sgs_ee, '.c', coord)
  $endif    
$endif

! Final checkpoint all restart data
call tavg_checkpoint()

$if($MPI)
call mpi_barrier( comm, ierr )
$endif

! Zero bogus values
call type_zero_bogus( tavg(:,:,nz) )

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
allocate(rs_zplane_tot(nz_tot))
allocate(cnpy_zplane_tot(nz_tot))
allocate(tavg_zplane_tot(nz_tot))

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
call mpi_sendrecv (tavg(:,:,1), nx*ny, MPI_TAVG, down, 1,  &
                   tavg(:,:,nz), nx*ny, MPI_TAVG, up, 1,   &
                   comm, status, ierr)
call mpi_sendrecv (tavg(:,:,nz-1), nx*ny, MPI_TAVG, up, 2,  &
                   tavg(:,:,0), nx*ny, MPI_TAVG, down, 2,   &
                   comm, status, ierr)
$endif

! Anything with velocity is on the w-grid so these
!   values for coord==0 at k=1 should be zeroed (otherwise bogus)
if (jzmin==0) then  !coord==0 only
  tavg(:,:,0:1)%u = 0.0_rprec
  tavg(:,:,0:1)%v = 0.0_rprec
  tavg(:,:,0:1)%w = 0.0_rprec
  tavg(:,:,0:1)%u2 = 0.0_rprec
  tavg(:,:,0:1)%v2 = 0.0_rprec
  tavg(:,:,0:1)%w2 = 0.0_rprec
  tavg(:,:,0:1)%uw = 0.0_rprec
  tavg(:,:,0:1)%vw = 0.0_rprec
  tavg(:,:,0:1)%uv = 0.0_rprec
endif

! Interpolate between grids where necessary
tavg = tavg_interp_to_w_grid( tavg )

! Anything with velocity is on the w-grid so these
!   values for coord==0 at k=1 should be zeroed (otherwise bogus)
if (jzmin==0) then  !coord==0 only
  tavg(:,:,0:1)%fx = 0.0_rprec 
  tavg(:,:,0:1)%fy = 0.0_rprec   

  tavg(:,:,0:1)%txx = 0.0_rprec    
  tavg(:,:,0:1)%txy = 0.0_rprec    
  tavg(:,:,0:1)%tyy = 0.0_rprec    
  tavg(:,:,0:1)%tzz = 0.0_rprec
endif

!  Average over z-planes
do k=1,nz
  
  !  Initialize to 0 for summations
  call type_set( tavg_zplane(k), 0._rprec )

  do j=1, Ny
    do i=1, Nx

      tavg_zplane(k) = tavg_zplane(k) .ADD. tavg(i,j,k)

    enddo
  enddo

  !  Divide by number of summation points 
  tavg_zplane(k) = tavg_zplane(k) .DIV. favg

enddo

! Compute the Reynolds stresses: bar(u_i * u_j) - bar(u_i) * bar(u_j)
rs = rs_compute( tavg , lbz)

! Compute planar averaged Reynolds stress
do k = 1, nz

  !  Initialize to 0
  call type_set( rs_zplane(k), 0._rprec)

  do j = 1, ny
    do i = 1, nx
      rs_zplane(k) = rs_zplane(k) .ADD. rs(i,j,k) 
    enddo    
  enddo

  rs_zplane(k) = rs_zplane(k) .DIV. favg
  
enddo


!  Compute the Canopy/Dispersive Stresses: <bar(u_i)*bar(u_j)>_xy - <bar(u_i)>_xy * <bar(u_j)>_xy
do k = 1, nz  

  ! Initialize to 0
  call type_set( cnpy_avg, 0._rprec)
  call type_set( tavg_avg, 0._rprec)

  do j=1, Ny
    do i=1, Nx
      cnpy_avg = cnpy_avg .ADD. cnpy_tavg_mul( tavg(i,j,k) )
      tavg_avg = tavg_avg .ADD. tavg(i,j,k)
    enddo
  enddo

  cnpy_avg = cnpy_avg .DIV. favg
  tavg_avg = tavg_avg .DIV. favg

  cnpy_zplane(k) = cnpy_avg .SUB. cnpy_tavg_mul( tavg_avg )

  !cnpy_zplane(k) % up2  = fa*sum(tavg(:,:,k)%u * tavg(:,:,k)%u) - (fa*sum( tavg(:,:,k)%u ))**2
  !cnpy_zplane(k) % vp2  = fa*sum(tavg(:,:,k)%v * tavg(:,:,k)%v) - (fa*sum( tavg(:,:,k)%v ))**2
  !cnpy_zplane(k) % wp2  = fa*sum(tavg(:,:,k)%w * tavg(:,:,k)%w) - (fa*sum( tavg(:,:,k)%w ))**2
  !cnpy_zplane(k) % upwp = fa*sum(tavg(:,:,k)%u * tavg(:,:,k)%w) - fa*sum( tavg(:,:,k)%u ) * fa*sum( tavg(:,:,k)%w )
  !cnpy_zplane(k) % vpwp = fa*sum(tavg(:,:,k)%v * tavg(:,:,k)%w) - fa*sum( tavg(:,:,k)%v ) * fa*sum( tavg(:,:,k)%w )
  !cnpy_zplane(k) % upvp = fa*sum(tavg(:,:,k)%u * tavg(:,:,k)%v) - fa*sum( tavg(:,:,k)%u ) * fa*sum( tavg(:,:,k)%v )
  
enddo

$if($MPI)
call mpi_barrier( comm, ierr )
$endif
$if($CGNS and $MPI)
    ! Write CGNS Data
    call write_parallel_cgns (fname_vel_cgns,nx,ny, nz - nz_end, nz_tot,        &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                      &
    3, (/ 'VelocityX', 'VelocityY', 'VelocityZ'/),                             &
    (/ tavg(1:nx,1:ny,1:nz- nz_end) % u,                                        &
    tavg(1:nx,1:ny,1:nz- nz_end) % v,                                           &
    tavg(1:nx,1:ny,1:nz- nz_end) % w /) )
    
    call write_parallel_cgns(fname_vel2_cgns,nx,ny,nz- nz_end,nz_tot,           &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                   &
    (/ 'Mean--uu', 'Mean--vv', 'Mean--ww','Mean--uw','Mean--vw','Mean--uv'/),  &
    (/ tavg(1:nx,1:ny,1:nz- nz_end) % u2,                                       &
    tavg(1:nx,1:ny,1:nz- nz_end) % v2,                                          &
    tavg(1:nx,1:ny,1:nz- nz_end) % w2,                                          &
    tavg(1:nx,1:ny,1:nz- nz_end) % uw,                                          &
    tavg(1:nx,1:ny,1:nz- nz_end) % vw,                                          &
    tavg(1:nx,1:ny,1:nz- nz_end) % uv /) )
    
    call write_parallel_cgns(fname_ddz_cgns,nx,ny,nz- nz_end,nz_tot,            &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 2, (/ 'dudz----', 'dvdz----'/),      &
    (/ tavg(1:nx,1:ny,1:nz- nz_end) % dudz,                                     &
    tavg(1:nx,1:ny,1:nz- nz_end) % dvdz /) )
    
    call write_parallel_cgns(fname_tau_cgns,nx,ny,nz- nz_end,nz_tot,            &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                   &
    (/ 'Tao--txx', 'Tao--txy', 'Tao--tyy','Tao--txz','Tao--tyz','Tao--tzz'/),  &
    (/ tavg(1:nx,1:ny,1:nz- nz_end) % txx,                                      &
    tavg(1:nx,1:ny,1:nz- nz_end) % txy,                                         &
    tavg(1:nx,1:ny,1:nz- nz_end) % tyy,                                         &
    tavg(1:nx,1:ny,1:nz- nz_end) % txz,                                         &
    tavg(1:nx,1:ny,1:nz- nz_end) % tyz,                                         &
    tavg(1:nx,1:ny,1:nz- nz_end) % tzz /) )  
      
    call write_parallel_cgns(fname_f_cgns,nx,ny,nz- nz_end,nz_tot,              &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 3,                                   &
    (/ 'bodyForX', 'bodyForY', 'bodyForZ'/),                                   &
    (/ tavg(1:nx,1:ny,1:nz- nz_end) % fx,                                       &
    tavg(1:nx,1:ny,1:nz- nz_end) % fy,                                          &
    tavg(1:nx,1:ny,1:nz- nz_end) % fz /) )
    
    call write_parallel_cgns(fname_rs_cgns,nx,ny,nz- nz_end,nz_tot,             &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                   &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ rs(1:nx,1:ny,1:nz) % up2,                                               &
    rs(1:nx,1:ny,1:nz- nz_end) % vp2,                                           &
    rs(1:nx,1:ny,1:nz- nz_end) % wp2,                                           &
    rs(1:nx,1:ny,1:nz- nz_end) % upwp,                                          &
    rs(1:nx,1:ny,1:nz- nz_end) % vpwp,                                          &
    rs(1:nx,1:ny,1:nz- nz_end) % upvp /) )
    
    call write_parallel_cgns(fname_cs_cgns,nx,ny,nz- nz_end,nz_tot,             &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                                   &
    (/ 'Cs Coeff'/),  (/ tavg(1:nx,1:ny,1:nz- nz_end) % cs_opt2 /) )
$endif

! ----- Write all the 3D data -----
$if($BINARY)
open(unit=13,file=fname_velb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%u
write(13,rec=2) tavg(:nx,:ny,1:nz)%v
write(13,rec=3) tavg(:nx,:ny,1:nz)%w
close(13)

open(unit=13,file=fname_vel2b,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%u2
write(13,rec=2) tavg(:nx,:ny,1:nz)%v2
write(13,rec=3) tavg(:nx,:ny,1:nz)%w2
write(13,rec=4) tavg(:nx,:ny,1:nz)%uw
write(13,rec=5) tavg(:nx,:ny,1:nz)%vw
write(13,rec=6) tavg(:nx,:ny,1:nz)%uv
close(13)

open(unit=13,file=fname_ddzb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%dudz
write(13,rec=2) tavg(:nx,:ny,1:nz)%dvdz
close(13)

open(unit=13,file=fname_taub,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%txx
write(13,rec=2) tavg(:nx,:ny,1:nz)%txy
write(13,rec=3) tavg(:nx,:ny,1:nz)%tyy
write(13,rec=4) tavg(:nx,:ny,1:nz)%txz
write(13,rec=5) tavg(:nx,:ny,1:nz)%tyz
write(13,rec=6) tavg(:nx,:ny,1:nz)%tzz
close(13)

open(unit=13,file=fname_fb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%fx
write(13,rec=2) tavg(:nx,:ny,1:nz)%fy
write(13,rec=3) tavg(:nx,:ny,1:nz)%fz
close(13)

open(unit=13,file=fname_rsb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) rs(:nx,:ny,1:nz)%up2 
write(13,rec=2) rs(:nx,:ny,1:nz)%vp2 
write(13,rec=3) rs(:nx,:ny,1:nz)%wp2 
write(13,rec=4) rs(:nx,:ny,1:nz)%upwp
write(13,rec=5) rs(:nx,:ny,1:nz)%vpwp
write(13,rec=6) rs(:nx,:ny,1:nz)%upvp
close(13)

open(unit=13,file=fname_csb,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tavg(:nx,:ny,1:nz)%cs_opt2 
close(13)

$else

  $if ($LVLSET)

!~   call write_tecplot_header_ND(fname_vel, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<u>","<v>","<w>"', numtostr(coord, 6), 2)
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_vel, 'append', 'formatted', 4, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        tavg(:,:,1:nz) % u, &
!~        tavg(:,:,1:nz) % v, &
!~        tavg(:,:,1:nz) % w /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_vel2, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
!~        numtostr(coord,6), 2)
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_vel2, 'append', 'formatted', 7, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        tavg(:,:,1:nz) % u2, &
!~        tavg(:,:,1:nz) % v2, &
!~        tavg(:,:,1:nz) % w2, &
!~        tavg(:,:,1:nz) % uw, &
!~        tavg(:,:,1:nz) % vw, &
!~        tavg(:,:,1:nz) % uv /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_ddz, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_ddz, 'append', 'formatted', 3, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        tavg(:,:,1:nz) % dudz, &
!~        tavg(:,:,1:nz) % dvdz /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_tau, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
!~        numtostr(coord,6), 2)  
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_tau, 'append', 'formatted', 7, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        tavg(:,:,1:nz) % txx, &
!~        tavg(:,:,1:nz) % txy, &
!~        tavg(:,:,1:nz) % tyy, &
!~        tavg(:,:,1:nz) % txz, &
!~        tavg(:,:,1:nz) % tyz, &
!~        tavg(:,:,1:nz) % tzz /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_f, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', &
!~        numtostr(coord,6), 2)
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_f, 'append', 'formatted', 4, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        tavg(:,:,1:nz) % fx, &
!~        tavg(:,:,1:nz) % fy, &
!~        tavg(:,:,1:nz) % fz /), &
!~        4, x, y, zw(1:nz))

    $if($MPI)

    call mpi_allreduce(sum(tavg(1:nx,1:ny,1:nz-1)%fx), fx_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
    call mpi_allreduce(sum(tavg(1:nx,1:ny,1:nz-1)%fy), fy_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
    call mpi_allreduce(sum(tavg(1:nx,1:ny,1:nz-1)%fz), fz_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)

    $else

    fx_global = sum(tavg(1:nx,1:ny,1:nz-1)%fx)
    fy_global = sum(tavg(1:nx,1:ny,1:nz-1)%fy)
    fz_global = sum(tavg(1:nx,1:ny,1:nz-1)%fz)
    
    $endif

!~   call write_tecplot_header_ND(fname_rs, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', &
!~        numtostr(coord,6), 2)  
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_rs, 'append', 'formatted', 7, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        rs(:,:,1:nz)%up2, &
!~        rs(:,:,1:nz)%vp2, &
!~        rs(:,:,1:nz)%wp2, &
!~        rs(:,:,1:nz)%upwp, &
!~        rs(:,:,1:nz)%vpwp, &
!~        rs(:,:,1:nz)%upvp /), &
!~        4, x, y, zw(1:nz))
!~   
!~   call write_tecplot_header_ND(fname_cs, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<cs2>"', numtostr(coord,6), 2)
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_cs, 'append', 'formatted', 2, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz), &
!~        tavg(:,:,1:nz) % cs_opt2 /), &
!~        4, x, y, zw(1:nz))

  $else

!~   call write_tecplot_header_ND(fname_vel, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<u>","<v>","<w>"', numtostr(coord, 6), 2)
!~   call write_real_data_3D(fname_vel, 'append', 'formatted', 3, nx, ny, nz, &
!~        (/ tavg(:,:,1:nz) % u, &
!~        tavg(:,:,1:nz) % v, &
!~        tavg(:,:,1:nz) % w /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_vel2, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
!~        numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_vel2, 'append', 'formatted', 6, nx, ny, nz, &
!~        (/ tavg(:,:,1:nz) % u2, &
!~        tavg(:,:,1:nz) % v2, &
!~        tavg(:,:,1:nz) % w2, &
!~        tavg(:,:,1:nz) % uw, &
!~        tavg(:,:,1:nz) % vw, &
!~        tavg(:,:,1:nz) % uv /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_ddz, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_ddz, 'append', 'formatted', 2, nx, ny, nz, &
!~        (/ tavg(:,:,1:nz) % dudz, &
!~        tavg(:,:,1:nz) % dvdz /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_tau, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
!~        numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_tau, 'append', 'formatted', 6, nx, ny, nz, &
!~        (/ tavg(:,:,1:nz) % txx, &
!~        tavg(:,:,1:nz) % txy, &
!~        tavg(:,:,1:nz) % tyy, &
!~        tavg(:,:,1:nz) % txz, &
!~        tavg(:,:,1:nz) % tyz, &
!~        tavg(:,:,1:nz) % tzz /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_f, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', &
!~        numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_f, 'append', 'formatted', 3, nx, ny, nz, &
!~        (/ tavg(:,:,1:nz) % fx, &
!~        tavg(:,:,1:nz) % fy, &
!~        tavg(:,:,1:nz) % fz /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_rs, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', &
!~        numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_rs, 'append', 'formatted', 6, nx, ny, nz, &
!~        (/ rs(:,:,1:nz)%up2, &
!~        rs(:,:,1:nz)%vp2, &
!~        rs(:,:,1:nz)%wp2, &
!~        rs(:,:,1:nz)%upwp, &
!~        rs(:,:,1:nz)%vpwp, &
!~        rs(:,:,1:nz)%upvp /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_cs, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<cs2>"', numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_cs, 'append', 'formatted', 1, nx, ny, nz, &
!~        (/ tavg(:,:,1:nz)% cs_opt2 /), &
!~        4, x, y, zw(1:nz))

  $endif

$endif

!----
$if($OUTPUT_EXTRA)

  $if($LVLSET)

!~   call write_tecplot_header_ND(fname_sgs_TnNu, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<Tn>", "<Nu_t>"', numtostr(coord,6), 2)
!~   !  write phi and x,y,z
!~   call write_real_data_3D(fname_sgs_TnNu, 'append', 'formatted', 1, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, zw(1:nz))
!~   call write_real_data_3D(fname_sgs_TnNu, 'append', 'formatted', 2, nx, ny, nz, &
!~        (/ tavg_sgs % Tn, tavg_sgs % Nu_t /), 4)  
!~   
!~   call write_tecplot_header_ND(fname_sgs_Fsub, 'rewind', 8, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "phi", "<F_LM>", "<F_MM>", "<F_QN>", "<F_NN>"', numtostr(coord,6), 2)
!~     !  write phi and x,y,z
!~   call write_real_data_3D(fname_sgs_Fsub, 'append', 'formatted', 1, nx, ny, nz, &
!~        (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, zw(1:nz))
!~   call write_real_data_3D(fname_sgs_Fsub, 'append', 'formatted', 4, nx, ny, nz, &
!~        (/ tavg_sgs % F_LM, tavg_sgs % F_MM, tavg_sgs % F_QN, tavg_sgs % F_NN /), 4)  
!~ 
  $else
!~   
!~   call write_tecplot_header_ND(fname_sgs_TnNu, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<Tn>", "<Nu_t>"', numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_sgs_TnNu, 'append', 'formatted', 2, nx, ny, nz, &
!~        (/ tavg_sgs % Tn, tavg_sgs % Nu_t /), &
!~        4, x, y, zw(1:nz))
!~ 
!~   call write_tecplot_header_ND(fname_sgs_Fsub, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
!~        '"x", "y", "z", "<F_LM>", "<F_MM>", "<F_QN>", "<F_NN>"', numtostr(coord,6), 2)
!~   call write_real_data_3D(fname_sgs_Fsub, 'append', 'formatted', 4, nx, ny, nz, &
!~        (/ tavg_sgs % F_LM, tavg_sgs % F_MM, tavg_sgs % F_QN, tavg_sgs % F_NN /), &
!~        4, x, y, zw(1:nz))

  $endif
  
    
  $if($DYN_TN)
    
    $if($LVLSET)

!~     call write_tecplot_header_ND(fname_sgs_ee, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
!~          '"x", "y", "z", "phi", "<ee_now>", "<F_ee2>", "<F_deedt2>"', numtostr(coord,6), 2)
!~     !  write phi and x,y,z
!~     call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 1, nx, ny, nz, &
!~          (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, zw(1:nz))
!~     call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 3, nx, ny, nz, &
!~          (/ tavg_sgs % ee_now, tavg_sgs % F_ee2, tavg_sgs % F_deedt2 /), 4)  
!~ 
    $else
!~     
!~     call write_tecplot_header_ND(fname_sgs_ee, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
!~          '"x", "y", "z", "<ee_now>", "<F_ee2>", "<F_deedt2>"', numtostr(coord,6), 2)
!~     call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 3, nx, ny, nz, &
!~          (/ tavg_sgs % ee_now, tavg_sgs % F_ee2, tavg_sgs % F_deedt2 /), &
!~          4, x, y, zw(1:nz))
    
    $endif        
    
  $else
  
    $if($LVLSET)

!~     call write_tecplot_header_ND(fname_sgs_ee, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
!~          '"x", "y", "z", "phi", "<ee_now>"', numtostr(coord,6), 2)
!~     !  write phi and x,y,z
!~     call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 1, nx, ny, nz, &
!~          (/ phi(1:nx,1:ny,1:nz) /), 4, x, y, zw(1:nz))
!~     call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 1, nx, ny, nz, &
!~          (/ tavg_sgs % ee_now /), 4)  
    
    $else

!~     call write_tecplot_header_ND(fname_sgs_ee, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
!~          '"x", "y", "z", "<ee_now>"', numtostr(coord,6), 2)
!~     call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 1, nx, ny, nz, &
!~          (/ tavg_sgs % ee_now /), &
!~          4, x, y, zw(1:nz))
    
    $endif        

  $endif
    
$endif
!----


! Construct zplane data 
$if($MPI)

  allocate(recvsize(nproc),recvstride(nproc))

  sendsize=nz-1
  recvsize=nz-1
  recvstride=coord_of_rank*(nz-1)

  ! Set very bottom values
  if( coord == 0 ) then

     z_tot(1) = z(1)
     zw_tot(1) = zw(1)
     rs_zplane_tot(1) = rs_zplane(1)
     cnpy_zplane_tot(1) = cnpy_zplane(1)
     tavg_zplane_tot(1) = tavg_zplane(1)

  endif

  call mpi_gatherv( z(2), sendsize, MPI_RPREC, &
       z_tot(2), recvsize, recvstride, MPI_RPREC, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( zw(2), sendsize, MPI_RPREC, &
       zw_tot(2), recvsize, recvstride, MPI_RPREC, &
       rank_of_coord(0), comm, ierr)
  
  call mpi_gatherv( rs_zplane(2), sendsize, MPI_RS, &
       rs_zplane_tot(2), recvsize, recvstride, MPI_RS, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( cnpy_zplane(2), sendsize, MPI_CNPY, &
       cnpy_zplane_tot(2), recvsize, recvstride, MPI_CNPY, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( tavg_zplane(2), sendsize, MPI_TAVG, &
       tavg_zplane_tot(2), recvsize, recvstride, MPI_TAVG, &
       rank_of_coord(0), comm, ierr)

  deallocate(recvsize, recvstride)

  call MPI_Type_free (MPI_RS, ierr)
  call MPI_Type_free (MPI_CNPY, ierr)
  call mpi_type_free (MPI_TAVG, ierr)    

  ! Write reconstructed data only if bottom processor
  if( coord == 0 ) then
  
!~     call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz_tot /), &
!~       '"z", "<u>","<v>","<w>"', numtostr(coord,6), 2)
!~     call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, Nz_tot, &
!~       (/ tavg_zplane_tot % u, tavg_zplane_tot % v, tavg_zplane_tot % w /), 0, zw_tot)
!~ 
!~     call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz_tot/), &
!~       '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
!~       numtostr(coord,6), 2)
!~     call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, Nz_tot, &
!~       (/ tavg_zplane_tot % u2, tavg_zplane_tot % v2, tavg_zplane_tot % w2, &
!~       tavg_zplane_tot % uw, tavg_zplane_tot % vw, tavg_zplane_tot % uv /), &
!~       0, zw_tot) 
!~   
!~     call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz_tot/), &
!~       '"z", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
!~     call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, Nz_tot, &
!~       (/ tavg_zplane_tot % dudz, tavg_zplane_tot % dvdz /), 0, zw_tot)
!~   
!~     call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz_tot/), &
!~       '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
!~       numtostr(coord,6), 2)  
!~     call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, Nz_tot, &
!~       (/ tavg_zplane_tot % txx, tavg_zplane_tot % txy, tavg_zplane_tot % tyy, &
!~       tavg_zplane_tot % txz, tavg_zplane_tot % tyz, tavg_zplane_tot % tzz /), 0, zw_tot) 
!~   
!~     call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz_tot/), &
!~       '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', numtostr(coord,6), 2)
!~     call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, Nz_tot, &
!~       (/ tavg_zplane_tot % fx, tavg_zplane_tot % fy, tavg_zplane_tot % fz /), 0, zw_tot)  
!~   
!~     call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz_tot/), &
!~       '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
!~     call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, Nz_tot, &
!~       (/ rs_zplane_tot % up2, rs_zplane_tot%vp2, rs_zplane_tot%wp2, &
!~       rs_zplane_tot%upwp, rs_zplane_tot%vpwp, rs_zplane_tot%upvp /), 0, zw_tot)    
!~  
!~     call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz_tot/), &
!~       '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
!~     call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, Nz_tot, &
!~       (/ cnpy_zplane_tot % up2, cnpy_zplane_tot%vp2, cnpy_zplane_tot%wp2, &
!~       cnpy_zplane_tot%upwp, cnpy_zplane_tot%vpwp, cnpy_zplane_tot%upvp /), 0, zw_tot)         
!~       
!~     call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz_tot/), &
!~       '"z", "<cs2>"', numtostr(coord,6), 2)
!~     call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, Nz_tot, &
!~       (/ tavg_zplane_tot % cs_opt2 /), 0, zw_tot)        

    deallocate(z_tot, zw_tot)      
    deallocate(tavg_zplane_tot)
    deallocate(rs_zplane_tot)
    deallocate(cnpy_zplane_tot)
  
  endif

$else

!~ call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz /), &
!~    '"z", "<u>","<v>","<w>"', numtostr(coord,6), 2)  
!~ call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, nz, &
!~   (/ tavg_zplane % u, tavg_zplane % v, tavg_zplane % w /), 0, zw(1:nz))
!~ 
!~ call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz/), &
!~    '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
!~    numtostr(coord,6), 2)
!~ call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, nz, &
!~   (/ tavg_zplane % u2, tavg_zplane % v2, tavg_zplane % w2, &
!~   tavg_zplane % uw, tavg_zplane % vw, tavg_zplane % uv /), 0, zw(1:nz)) 
!~   
!~ call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz/), &
!~    '"z", "<dudz>","<dvdz>"', numtostr(coord,6), 2)
!~ call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, nz, &
!~   (/ tavg_zplane % dudz, tavg_zplane % dvdz /), 0, zw(1:nz))
!~   
!~ call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz/), &
!~    '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
!~    numtostr(coord,6), 2)  
!~ call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, nz, &
!~   (/ tavg_zplane % txx, tavg_zplane % txy, tavg_zplane % tyy, &
!~   tavg_zplane % txz, tavg_zplane % tyz, tavg_zplane % tzz /), 0, zw(1:nz)) 
!~   
!~ call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz/), &
!~    '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', numtostr(coord,6), 2)
!~ call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, nz, &
!~   (/ tavg_zplane % fx, tavg_zplane % fy, tavg_zplane % fz /), 0, zw(1:nz))  
!~   
!~ call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz/), &
!~    '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
!~ call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, nz, &
!~   (/ rs_zplane % up2, rs_zplane%vp2, rs_zplane%wp2, &
!~   rs_zplane%upwp, rs_zplane%vpwp, rs_zplane%upvp /), 0, zw(1:nz))
!~ 
!~ call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz/), &
!~    '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
!~ call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, nz, &
!~   (/ cnpy_zplane % up2, cnpy_zplane%vp2, cnpy_zplane%wp2, &
!~   cnpy_zplane%upwp, cnpy_zplane%vpwp, cnpy_zplane%upvp /), 0, zw(1:nz))  
!~   
!~ call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz/), &
!~    '"z", "<cs2>"', numtostr(coord,6), 2)
!~ call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, nz, &
!~   (/ tavg_zplane % cs_opt2 /), 0, zw(1:nz))    
  
$endif


deallocate(tavg, tavg_zplane, rs, rs_zplane, cnpy_zplane)
$if($OUTPUT_EXTRA)
deallocate(tavg_sgs)
$endif

nullify(x,y,z,zw)

$if($MPI)
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
$endif

return
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

! RICHARD: SPECTRA NOT YET CONSIDERED WITH FFTW3
$if (not $FFTW3)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : checkpoint_spectra_file
use messages
use param, only : spectra_nloc
use stat_defs, only : spectra, spectra_total_time, spectra_initialized,spectra_dt
implicit none

character (*), parameter :: sub_name = mod_name // '.spectra_init'
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst
integer :: k

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

fname = checkpoint_spectra_file
$if ($MPI)
call string_concat( fname, MPI_suffix, coord )
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
   !  Nothing to read in
   if (coord == 0) then
      write(*,*) ' '
      write(*,*)'No previous spectra data - starting from scratch.'
   endif

   spectra_total_time = 0._rprec
   ! note: spectra already initialized during output_init

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

   read (1) spectra_total_time
   do k=1, spectra_nloc
      read (1) spectra(k) % power
   enddo

   close(1)

endif

! Initialize spectra_dt
spectra_dt = 0.0_rprec

! Set global switch that spectra as been initialized
spectra_initialized = .true.

! Set global switch that spectra as been initialized
spectra_initialized = .false.

return
end subroutine spectra_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_compute()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use sim_param, only : u
use param, only : nx,ny,dz,lh,spectra_nloc
use stat_defs, only : spectra, spectra_total_time, spectra_dt
use functions, only : linear_interp
implicit none

integer :: i, j, k
real(rprec), allocatable, dimension(:) :: ui, uhat, power

real(rprec) :: nxny

! Interpolation variable
allocate(ui(nx), uhat(nx), power(lh))


nxny = real(nx*ny,rprec)

!  Loop over all spectra locations
do k=1, spectra_nloc
  
  $if ($MPI)
  if(spectra(k) % coord == coord) then
  $endif

  do j=1,Ny

    do i=1,Nx
        ! Interpolate to the specified plane
        ui(i) = linear_interp(u(i,j,spectra(k) % istart),u(i,j,spectra(k) % istart+1), &
                              dz, spectra(k) % ldiff)

    enddo

    ! 1) Compute uhat for the given j
    ui = ui - sum(ui) / Nx ! Remove the mean
    ! Compute FFT
    $if ($FFTW3)
    write(*,*) 'spectra not calculated yet in FFTW3 mode'
    $else
    call rfftw_f77_one(forw_spectra, ui, uhat)
    $endif
    !  Normalize
    uhat = uhat / Nx

    ! 2) Compute power spectra for given j
    power(1) = uhat(1)**2
    do i=2,lh-1
      power(i) = uhat(i)**2 + uhat(Nx-i+2)**2
    enddo
    power(lh) = uhat(lh)**2 ! Nyquist

    ! Sum jth component, normalize, and include averaging over y 
    spectra(k) % power = spectra(k) % power + spectra_dt * power / nxny

  enddo

  $if ($MPI)
  endif
  $endif

enddo  

!  Update spectra total time
spectra_total_time = spectra_total_time + spectra_dt
  
deallocate(ui, uhat, power)

! Set spectra_dt back to zero for next use
spectra_dt = 0.0_rprec

return
end subroutine spectra_compute

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_finalize()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : lh, spectra_nloc, spectra_loc
use fft, only : kx
use stat_defs, only : spectra, spectra_total_time
implicit none

!~ include 'tecryte.h'

$if($VERBOSE)
character (*), parameter :: sub_name = mod_name // '.spectra_finalize'
$endif
character (64) :: fname

integer ::  k

! Perform final checkpoint 
call spectra_checkpoint()

!  Loop over all zplane locations
do k=1,spectra_nloc
  
  $if ($MPI)
  if(spectra(k) % coord == coord) then
  $endif

  ! Finalize averaging for power spectra (averaging over time)
  spectra(k) % power = spectra(k) % power / spectra_total_time

  !  Create unique file name
  call string_splice( fname, path // 'output/spectra.z-', spectra_loc(k),'.dat' )

  !  Omitting Nyquist from output
!~   call write_tecplot_header_ND(fname, 'rewind', 2, (/ lh-1/), &
!~     '"k", "E(k)"', numtostr(k, 6), 2 ) 
!~   call write_real_data_1D(fname, 'append', 'formatted', 1, lh-1, &
!~     (/ spectra(k) % power(1:lh-1) /), 0, (/ kx(1:lh-1,1) /))
!~     
  $if ($MPI)
  endif
  $endif

enddo  
  
deallocate(spectra)

return

end subroutine spectra_finalize

!*******************************************************************************
subroutine spectra_checkpoint()
!*******************************************************************************
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'spectra_finalize' at the end of the 
! simulation.    
! 
use param, only : checkpoint_spectra_file, spectra_nloc
$if($MPI)
use param, only : comm, ierr
$endif
use stat_defs, only : spectra, spectra_total_time

implicit none

character (*), parameter :: sub_name = mod_name // '.spectra_checkpoint'
character (64) :: fname

logical :: opn
integer :: k

fname = checkpoint_spectra_file
$if($MPI)
call string_concat( fname, '.c', coord)
$endif

!  Write data to spectra.out
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

! write the accumulated time and power
write (1) spectra_total_time
do k=1, spectra_nloc
  write (1) spectra(k) % power
enddo
close(1)

$if($MPI)
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
$endif

return
end subroutine spectra_checkpoint
$endif

end module io
