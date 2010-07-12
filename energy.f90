subroutine energy (ke)
use types,only:rprec
use param
use sim_param,only:u,v,w
use messages
$if ($XLF)
  use ieee_arithmetic  !--for NAN checking
$endif
implicit none

character (*), parameter :: sub_name = 'energy'

integer, parameter :: NAN_MAX = 10
                      !--write this many NAN's before calling error (to aid
                      !  diagnosis of problem)

$if ($DEBUG)
logical, parameter :: DEBUG = .true.
$endif

logical, parameter :: flush = .false.

integer :: jx, jy, jz
integer :: nan_count

$if($DEBUG)
logical :: nan
$endif

real(kind=rprec)::KE,denom,temp_w
$if ($MPI)
  real (rprec) :: ke_global
$endif

!---------------------------------------------------------------------

nan_count = 0

ke=0._rprec
denom=2._rprec*(nx*ny*(nz-1))

z_loop: do jz=1,nz-1
    do jy=1,ny
        do jx=1,nx
            !  Should perform trilinear interpolation on the cell
            temp_w=.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
            ke=ke+(u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)/denom
            !ke=ke+(u(jx,jy,jz)**2+v(jx,jy,jz)**2+interp_to_uv_grid('w',jx,jy,jz)**2)/denom
            
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
                    write (*, *) 'jt = ', jt
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

if ( nan_count > 0 ) call error (sub_name, 'NaN found')

$if ($MPI)

  call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  if (rank == 0) then  !--note its rank here, not coord
    ke = ke_global/nproc
    write (13, *) (jt_total-1) * dt, ke
  end if
  !if (rank == 0) ke = ke_global/nproc  !--its rank here, not coord

$else

  write (13, *) (jt_total-1) * dt, ke

$endif

if ( flush ) then
  if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
    close (13)
    open ( 13, file=path//'output/check_ke.out', position='append' )
  end if
end if

end subroutine energy
