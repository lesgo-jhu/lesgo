!**********************************************************************
program tsum_post
!**********************************************************************
use types, only : rprec
use mpi_defs
use grid_defs, only : x,y,z
use stat_defs, only : tavg_t, tsum_t, rs_t
use param, only : nx, ny, nz, USE_MPI

implicit none

logical, parameter :: rs_output=.true.
logical, parameter :: uvw_avg_output=.true.
integer, parameter :: iter_start=390, iter_stop=540, iter_skip=50 ! In thousands
character(50) :: ci,fname
character(50) :: fin, ftec,fdir
integer :: i,j,k,n
integer :: nf,ndirs
real(rprec) :: favg

call initialize_mpi()

ndirs = (iter_stop - iter_start)/iter_skip + 1 ! # of iteration sets

favg = 1._rprec/(ndirs * iter_skip * 1000._rprec ) ! 1/(total # of iterations)

!  Allocate and initialize
allocate(x(nx),y(ny),z(nz));

allocate(tavg_t%u(nx, ny, nz))
allocate(tavg_t%v(nx, ny, nz))
allocate(tavg_t%w(nx, ny, nz))
allocate(tavg_t%u2(nx, ny, nz))
allocate(tavg_t%v2(nx, ny, nz))
allocate(tavg_t%w2(nx, ny, nz))
allocate(tavg_t%uw(nx, ny, nz))
allocate(tavg_t%vw(nx, ny, nz))
allocate(tavg_t%uv(nx, ny, nz))
allocate(tavg_t%dudz(nx, ny, nz))

allocate(tsum_t%u(nx, ny, nz))
allocate(tsum_t%v(nx, ny, nz))
allocate(tsum_t%w(nx, ny, nz))
allocate(tsum_t%u2(nx, ny, nz))
allocate(tsum_t%v2(nx, ny, nz))
allocate(tsum_t%w2(nx, ny, nz))
allocate(tsum_t%uw(nx, ny, nz))
allocate(tsum_t%vw(nx, ny, nz))
allocate(tsum_t%uv(nx, ny, nz))
allocate(tsum_t%dudz(nx, ny, nz))

tavg_t%u=0.
tavg_t%v=0.
tavg_t%w=0.
tavg_t%u2=0.
tavg_t%v2=0.
tavg_t%w2=0.
tavg_t%uw=0.
tavg_t%vw=0.
tavg_t%uv=0.
tavg_t%dudz=0.

tsum_t%u=0.
tsum_t%v=0.
tsum_t%w=0.
tsum_t%u2=0.
tsum_t%v2=0.
tsum_t%w2=0.
tsum_t%uw=0.
tsum_t%vw=0.
tsum_t%uv=0.
tsum_t%dudz=0.

do nf=1,ndirs
  write(ci,'(i0)') iter_start + (nf - 1)*iter_skip
  fdir =  'output.' // trim(adjustl(ci)) // 'k'
  fname =  trim(adjustl(fdir)) // '/tsum.out'
 
  write(*,*) fname 

  call load_data(fname) !,x,y,z,u,v,w   

  tavg_t%u = tavg_t%u + favg*tsum_t%u
  tavg_t%v = tavg_t%v + favg*tsum_t%v
  tavg_t%w = tavg_t%w + favg*tsum_t%w
  tavg_t%u2 = tavg_t%u2 + favg*tsum_t%u2
  tavg_t%v2 = tavg_t%v2 + favg*tsum_t%v2
  tavg_t%w2 = tavg_t%w2 + favg*tsum_t%w2
  tavg_t%uw = tavg_t%uw + favg*tsum_t%uw
  tavg_t%vw = tavg_t%vw + favg*tsum_t%vw
  tavg_t%uv = tavg_t%uv + favg*tsum_t%uv
  tavg_t%dudz = tavg_t%dudz + favg*tsum_t%dudz

enddo

deallocate(tsum_t%u, &
tsum_t%v, &
tsum_t%w, &
tsum_t%u2, &
tsum_t%v2, &
tsum_t%w2, &
tsum_t%uw, &
tsum_t%vw, &
tsum_t%uv, &
tsum_t%dudz)

if(rs_compute) then
  allocate(rs_t%up2(nx, ny, nz))
  allocate(rs_t%vp2(nx, ny, nz))
  allocate(rs_t%wp2(nx, ny, nz))
  allocate(rs_t%upwp(nx, ny, nz))
  allocate(rs_t%vpwp(nx, ny, nz))
  allocate(rs_t%upvp(nx, ny, nz))
  rs_t%up2=0.
  rs_t%vp2=0.
  rs_t%wp2=0.
  rs_t%upwp=0.
  rs_t%vpwp=0.
  rs_t%upvp=0.
  do k=1,nz
    do j=1,ny
      do i=1,nx
        rs_t%up2(i,j,k)=tavg_t%u2(i,j,k) - tavg_t%u(i,j,k)*tavg_t%u(i,j,k)
        rs_t%vp2(i,j,k)=tavg_t%v2(i,j,k) - tavg_t%v(i,j,k)*tavg_t%v(i,j,k)
        rs_t%wp2(i,j,k)=tavg_t%w2(i,j,k) - tavg_t%w(i,j,k)*tavg_t%w(i,j,k)
        rs_t%upwp(i,j,k)=tavg_t%uw(i,j,k) - tavg_t%u(i,j,k)*tavg_t%w(i,j,k)
        rs_t%vpwp(i,j,k)=tavg_t%vw(i,j,k) - tavg_t%v(i,j,k)*tavg_t%w(i,j,k)
        rs_t%upvp(i,j,k)=tavg_t%uv(i,j,k) - tavg_t%u(i,j,k)*tavg_t%v(i,j,k)
      enddo
    enddo
  enddo
  write(ftec,*) 'rs-',iter_start,'k-',iter_stop,'k.dat'
  ftec = trim(adjustl(ftec))

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') mpirank
  ftec = trim (ftec) // temp
$endif

!  Create tecplot formatted velocity field file  
  open (unit = 7,file = ftec, status='unknown',form='formatted', &
    action='write',position='rewind')
  write(7,*) 'variables= "x", "y", "z", "up2", "vp2", "wp2", "upwp", "vpwp", "upvp"'
  write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
  write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''  
  do k=1,nz
    do j=1,ny
      do i=1,nx
  !  Write spatially averaged, temporally averaged quantities   
       write(7,*) x(i), y(j), z(k), rs_t%up2(i,j,k), rs_t%vp2(i,j,k), rs_t%wp2(i,j,k), &
         rs_t%upwp(i,j,k), rs_t%vpwp(i,j,k), rs_t%upvp(i,j,k)
      enddo
    enddo
  enddo
  close(7)

  deallocate(rs_t%up2, rs_t%vp2, rs_t%wp2, &
    rs_t%upwp, rs_t%vpwp, rs_t%upvp)

endif

if(uvw_avg_output) then

  write(ftec,*) 'uvw_avg-',iter_start,'k-',iter_stop,'k.dat'

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') mpirank
  ftec = trim (ftec) // temp
$endif

  open (unit = 7,file = ftec, status='unknown',form='formatted', &
    action='write',position='rewind')

  write(7,*) 'variables= "x", "y", "z", "<u>", "<v>", "<w>"'
  write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
  write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''	

  do k=1,nz
    do j=1,ny
      do i=1,nx
        write(7,*) x(i), y(j), z(k), tavg_t%u(i,j,k), tavg_t%v(i,j,k), tavg_t%w(i,j,k)
      enddo
    enddo
  enddo
  close(7)
endif


deallocate(tavg_t%u, &
tavg_t%v, &
tavg_t%w, &
tavg_t%u2, &
tavg_t%v2, &
tavg_t%w2, &
tavg_t%uw, &
tavg_t%vw, &
tavg_t%uv, &
tavg_t%dudz)

call finalize_mpi()

stop
end program tsum_post

!**********************************************************************
subroutine load_data(fname) 
!**********************************************************************
!
!  This subroutine loads binary formatted data files for both
!  serial and parallel cases
!
use stat_defs, only : tsum_t
use param, only : nx, ny, nz, USE_MPI

implicit none

character(*), intent(IN) :: fname
character(len=120) :: fname_coord
integer :: i,j,k

if(.not. USE_MPI) then
  write(*,"(1a,1a)") ' Processing File : ', fname

  open(unit = 7,file = fname, status='old',form='unformatted', &
    action='read',position='rewind')

  do k=1,nz
    do j=1,ny
      do i=1,nx
        read(7)  x(i), y(j), z(k), tsum_t%u(i,j,k), tsum_t%v(i,j,k), tsum_t%w(i,j,k), &
         tsum_t%u(i,j,k), tsum_t%v(i,j,k), tsum_t%w(i,j,k), tsum_t%u2(i,j,k), &
         tsum_t%v2(i,j,k), tsum_t%w2(i,j,k), tsum_t%uw(i,j,k), &
         tsum_t%vw(i,j,k), tsum_t%uv(i,j,k), tsum_t%dudz(i,j,k)
      enddo
    enddo
  enddo

  close(7)
else

  write (temp, '(".c",i0)') mpirank
  fname_coord = trim (fname) // temp

  write(*,"(1a,1a)") ' Processing File : ', fname_coord
 
  open(unit = 7,file = fname_coord, status='old',form='unformatted', &
    action='read',position='rewind') 

!  Read data from input data file
    do k=1,nz
      do j=1,ny
        do i=1,nx
          read(7,*)  x(i), y(j), z(k), tsum_t%u(i,j,k), tsum_t%v(i,j,k), &
            tsum_t%w(i,j,k), tsum_t%u(i,j,k), tsum_t%v(i,j,k), &
            tsum_t%w(i,j,k), tsum_t%u2(i,j,k), tsum_t%v2(i,j,k), &
            tsum_t%w2(i,j,k), tsum_t%uw(i,j,k), tsum_t%vw(i,j,k), &
            tsum_t%uv(i,j,k), tsum_t%dudz(i,j,k)
        enddo
      enddo
    enddo

    close(7)

  enddo
endif

return
end subroutine load_data