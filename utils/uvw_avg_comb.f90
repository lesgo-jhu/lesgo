module io_data
use types, only : rprec
implicit none
real(rprec), allocatable, dimension(:) :: x,y,z
real(rprec), allocatable, dimension(:,:,:) :: u,v,w, ua, va, wa
end module io_data

program uvw_avg_comb
use io_data
use param, only : nx, ny, nz, nz_tot, pi, nproc, USE_MPI

implicit none

integer, parameter :: iter_start=390, iter_stop=540, iter_skip=50
character(50) :: ci,fname
character(50) :: fin, ftec,fdir
integer :: i,j,k,n
integer :: nf,ndirs

allocate(x(Nx),y(Ny),z(Nz_tot));
allocate(u(Nx,Ny,Nz_tot),v(Nx,Ny,Nz_tot),w(Nx,Ny,Nz_tot))
allocate(ua(Nx,Ny,Nz_tot),va(Nx,Ny,Nz_tot),wa(Nx,Ny,Nz_tot))

ndirs = (iter_stop - iter_start)/iter_skip + 1

do nf=1,ndirs
  write(ci,'(i0)') iter_start + (nf - 1)*iter_skip
  fdir =  'output.' // trim(adjustl(ci)) // 'k'
  fname =  trim(adjustl(fdir)) // '/uvw_avg.dat'
 
  write(*,*) fname 

  call load_data(fname) !,x,y,z,u,v,w   
  ua = ua + u/ndirs
  va = va + v/ndirs
  wa = wa + w/ndirs
enddo

!ua = ua/ndirs
!va = va/ndirs
!wa = wa/ndirs

write(ftec,*) 'uvw_avg-',iter_start,'k-',iter_stop,'k.dat'
ftec = trim(adjustl(ftec))
!  Create tecplot formatted velocity field file  
open (unit = 2,file = ftec, status='unknown',form='formatted', &
  action='write',position='rewind')
write(2,*) 'variables = "x", "y", "z", "<u>", "<v>", "<w>"'; 
write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
  nf,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz_tot
write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''

do k=1,nz_tot
  do j=1,ny
    do i=1,nx
      write(2,*) x(i), y(j), z(k), ua(i,j,k), va(i,j,k), wa(i,j,k)
    enddo
  enddo
enddo
close(2)

stop
end program uvw_avg_comb

!**********************************************************************
subroutine load_data(fname) 
!**********************************************************************
!
!  This subroutine loads binary formatted data files for both
!  serial and parallel cases
!
use io_data
use param, only : nx, ny, nz, nz_tot, pi, nproc, USE_MPI
implicit none
character(*), intent(IN) :: fname
character(len=120) :: fname_coord, temp,temp1,temp2,temp3
integer :: itemp, itemp2, itemp3
integer :: i,j,k
integer :: coord, nz_start

!if(USE_MPI) then
!  Same way that it is computed in lesgo
!  nz_proc = (nz-1)/nproc + 1
!  nz_tot = (nz_proc - 1)*nproc + 1
!else
!  nz_tot = nz
!  write(*,*) 'not using mpi'
!  stop
!endif


if(.not. USE_MPI) then
  write(*,"(1a,1a)") ' Processing File : ', fname

  open(unit = 7,file = fname, status='old',form='formatted', &
    action='read',position='rewind')
  read(7,"(1a)") temp
  read(7,"(1a)") temp
  read(7,"(1a)") temp

  do k=1,nz
    do j=1,ny
      do i=1,nx
        read(7,*) x(i), y(j), z(k), u(i,j,k), v(i,j,k), w(i,j,k)
      enddo
    enddo
  enddo

  close(7)
else
!  Initialize global z-starting value
  nz_start = 0

  do coord=0, nproc-1
  
    write (temp, '(".c",i0)') coord
    fname_coord = trim (fname) // temp

    write(*,"(1a,1a)") ' Processing File : ', fname_coord
 
    open(unit = 7,file = fname_coord, status='old',form='formatted', &
      action='read',position='rewind')

    read(7,"(1a)") temp
    read(7,"(1a)") temp
    read(7,"(1a)") temp
   

!  Read data from input data file
    do k=1,nz
      do j=1,ny
        do i=1,nx
          read(7,*) x(i), y(j), z(k + nz_start), u(i,j,k + nz_start), v(i,j,k + nz_start), w(i,j,k + nz_start)
        enddo
      enddo
    enddo

    close(7)

!  Update global z-starting value; takes into account the overlap at the block boundaries
    nz_start = nz_start + nz - 1

  enddo
endif

return
end subroutine load_data
