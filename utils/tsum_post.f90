!**********************************************************************
program tsum_post
!**********************************************************************
use types, only : rprec
! $if ($MPI)
! use mpi_defs
! $endif
use grid_defs, only : x,y,z
use stat_defs, only : tavg_t, tsum_t, rs_t
use param, only : nx, ny, nz, nproc, USE_MPI
use cylinder_skew_base_ls, only : phi

implicit none

logical, parameter :: rs_output=.true.
logical, parameter :: uvw_avg_output=.true.
logical, parameter :: tecio=.true.
integer, parameter :: iter_start=100, iter_stop=1000, iter_skip=100 ! In thousands
character(50) :: ci,fname,temp,fiter_start, fiter_stop
character(50) :: ftec, fdir
integer :: i,j,k
integer :: nf,ndirs,np
real(rprec) :: favg, rcount
real(rprec), allocatable :: sum_z(:)

! $if ($MPI)
! call initialize_mpi()
! $endif

do np=0,nproc-1
  ndirs = (iter_stop - iter_start)/iter_skip + 1 ! # of iteration sets

  favg = 1._rprec/(ndirs * iter_skip * 1000._rprec ) ! 1/(total # of iterations)
  write(*,*) '1/favg : ', 1./favg
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

  $if ($LVLSET)
  $if ($MPI)
  allocate(phi(nx+2,ny,0:nz))
  $else
  allocate(phi(nx+2,ny,1:nz))
  $endif
  $endif

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


  $if ($LVLSET)
!  Load phi data
  call load_data('phi.out',np)
  $endif

  do nf=1,ndirs
    write(ci,'(i0)') iter_start + (nf - 1)*iter_skip
    fdir =  'output.' // trim(adjustl(ci)) // 'k'
    fname =  trim(adjustl(fdir)) // '/tsum.out'
 
    write(*,*) fname 

    call load_data(fname,np) !,x,y,z,u,v,w   

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
 
  if(rs_output) then
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
    write(fiter_start, '(i0)') iter_start
    write(fiter_stop, '(i0)') iter_stop
    write(ftec,*) 'rs-'//trim(fiter_start)//'k-'//trim(fiter_stop)//'k.dat'
    ftec = trim(adjustl(ftec))
 
    $if ($MPI)
!  For MPI implementation     
    write (temp, '(".c",i0)') np
    ftec = trim (ftec) // temp
    $endif

!  Create tecplot formatted velocity field file  
    open (unit = 7,file = ftec, status='unknown',form='formatted', &
      action='write',position='rewind')
  
    if(tecio) then
      write(7,*) 'variables= "x", "y", "z", "up2", "vp2", "wp2", "upwp", "vpwp", "upvp", "phi"'
      write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
        1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
      $if ($LVLSET)
      write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''  
      $else
      write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''  
      $endif
    endif

    do k=1,nz
      do j=1,ny
        do i=1,nx
  !  Write spatially averaged, temporally averaged quantities
          $if ($LVLSET)   
          write(7,*) x(i), y(j), z(k), rs_t%up2(i,j,k), rs_t%vp2(i,j,k), rs_t%wp2(i,j,k), &
            rs_t%upwp(i,j,k), rs_t%vpwp(i,j,k), rs_t%upvp(i,j,k), phi(i,j,k)
          $else
          write(7,*) x(i), y(j), z(k), rs_t%up2(i,j,k), rs_t%vp2(i,j,k), rs_t%wp2(i,j,k), &
            rs_t%upwp(i,j,k), rs_t%vpwp(i,j,k), rs_t%upvp(i,j,k)
          $endif
        enddo
      enddo
    enddo
    close(7)

    write(fiter_start, '(i0)') iter_start
    write(fiter_stop, '(i0)') iter_stop
    write(ftec,*) 'rs_z-'//trim(fiter_start)//'k-'//trim(fiter_stop)//'k.dat'
    ftec = trim(adjustl(ftec))
 
    $if ($MPI)
!  For MPI implementation
    write (temp, '(".c",i0)') np
    ftec = trim (ftec) // temp
    $endif

!  Create tecplot formatted velocity field file
    open (unit = 7,file = ftec, status='unknown',form='formatted', &
      action='write',position='rewind')
 
    if(tecio) then
      write(7,*) 'variables= "z", "up2", "vp2", "wp2", "upwp", "vpwp", "upvp"'
      write(7,"(1a,i3,1a,i3)") 'ZONE T="', &
        1,'", DATAPACKING=POINT, k=', Nz
      write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
    endif

!  Allocate mem for Reynolds stress summation
    allocate(sum_z(6))
    do k=1,nz
      sum_z = 0._rprec
      rcount = 0.
      do j=1,ny
        do i=1,nx
          if(phi(i,j,k) >= 0._rprec) then
            rcount = rcount + 1.
            sum_z(1) = sum_z(1) + rs_t%up2(i,j,k)
            sum_z(2) = sum_z(2) + rs_t%vp2(i,j,k)
            sum_z(3) = sum_z(3) + rs_t%wp2(i,j,k)
            sum_z(4) = sum_z(4) + rs_t%upwp(i,j,k)
            sum_z(5) = sum_z(5) + rs_t%vpwp(i,j,k)
            sum_z(6) = sum_z(6) + rs_t%upvp(i,j,k)
          endif
        enddo
      enddo
      if(rcount > 0.5) then
        sum_z = sum_z / rcount
      else
        sum_z = 0._rprec
      endif
  
  !  Write spatially averaged, temporally averaged quantities
      write(7,*) z(k), sum_z
    enddo
    close(7)

    deallocate(rs_t%up2, rs_t%vp2, rs_t%wp2, &
      rs_t%upwp, rs_t%vpwp, rs_t%upvp)
    deallocate(sum_z)
  endif

  if(uvw_avg_output) then
    write(fiter_start, '(i0)') iter_start
    write(fiter_stop, '(i0)') iter_stop
    write(ftec,*) 'uvw_avg-'//trim(fiter_start)//'k-'//trim(fiter_stop)//'k.dat'
    ftec = trim(adjustl(ftec))
    $if ($MPI)
!  For MPI implementation     
    write (temp, '(".c",i0)') np
    ftec = trim (ftec) // temp
    $endif

    open (unit = 7,file = ftec, status='unknown',form='formatted', &
      action='write',position='rewind')

    if(tecio) then
      $if ($LVLSET)
      write(7,*) 'variables= "x", "y", "z", "<u>", "<v>", "<w>", "phi"'
      $else
      write(7,*) 'variables= "x", "y", "z", "<u>", "<v>", "<w>"'
      $endif
      write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
        1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
      $if ($LVLSET)
      write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''	
      $else
      write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''	
      $endif
    endif

    do k=1,nz
      do j=1,ny
        do i=1,nx
          $if ($LVLSET)
          write(7,*) x(i), y(j), z(k), tavg_t%u(i,j,k), tavg_t%v(i,j,k), tavg_t%w(i,j,k), phi(i,j,k)
          $else
          write(7,*) x(i), y(j), z(k), tavg_t%u(i,j,k), tavg_t%v(i,j,k), tavg_t%w(i,j,k)
          $endif         
        enddo
      enddo
    enddo
    close(7)

!  Perform z plane averaging
    write(fiter_start, '(i0)') iter_start
    write(fiter_stop, '(i0)') iter_stop
    write(ftec,*) 'uvw_avg_z-'//trim(fiter_start)//'k-'//trim(fiter_stop)//'k.dat'
    ftec = trim(adjustl(ftec))
    $if ($MPI)
!  For MPI implementation
    write (temp, '(".c",i0)') np
    ftec=trim(ftec) // temp
    $endif

    open (unit = 7,file = ftec, status='unknown',form='formatted', &
      action='write',position='rewind')

    if(tecio) then
      write(7,*) 'variables= "z", "<u>", "<v>", "<w>"'
      write(7,"(1a,i9,1a,i3)") 'ZONE T="', &
        1,'", DATAPACKING=POINT, k=', Nz
      write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE)')//''
    endif
!  Allocate mem for Velocity average stress summation
    allocate(sum_z(3))
    do k=1,nz
      sum_z = 0._rprec
      rcount = 0.
      do j=1,ny
        do i=1,nx
          if(phi(i,j,k) >= 0._rprec) then
            rcount = rcount + 1.
            sum_z(1) = sum_z(1) + tavg_t%u(i,j,k)
            sum_z(2) = sum_z(2) + tavg_t%v(i,j,k)
            sum_z(3) = sum_z(3) + tavg_t%w(i,j,k)
          endif
        enddo
      enddo
!  Make sure there is at least 1 point to be averaged 
      if(rcount > 0.5) then
        sum_z = sum_z / rcount
      else
        sum_z = 0._rprec
      endif

  !  Write spatially averaged, temporally averaged quantities
      write(7,*) z(k), sum_z

    enddo
    close(7)
    deallocate(sum_z)

  endif

  deallocate(x,y,z)
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
  

! $if ($MPI)
! call finalize_mpi()
! $endif

enddo

stop
end program tsum_post

!**********************************************************************
subroutine load_data(fbase,np) 
!**********************************************************************
!
!  This subroutine loads binary formatted data files for both
!  serial and parallel cases
!
! $if ($MPI)
! use mpi_defs, only : mpirank
! $endif
use grid_defs, only : x, y, z
use stat_defs, only : tsum_t
use param, only : nx, ny, nz, USE_MPI
use cylinder_skew_base_ls, only : phi

implicit none

character(*), intent(IN) :: fbase
integer, intent(IN) :: np
character(50) :: fname
character(15) :: temp
integer :: i,j,k

fname = trim(fbase)

$if ($MPI)
  write (temp, '(".c",i0)') np
  fname = trim (fname) // temp
$endif


write(*,"(1a,1a)") ' Processing File : ', fname
 
open(unit = 7,file = fname, status='old',form='unformatted', &
  action='read',position='rewind') 

if(fbase  == 'phi.out') then

!  Read binary data for lesgo
  open (unit=7, file=fname, status='old', form='unformatted', &
    action='read', position='rewind')

  read(7) phi
  close (7)

else
!  Read data from input data file
  do k=1,nz
    do j=1,ny
      do i=1,nx
        read(7)  x(i), y(j), z(k), tsum_t%u(i,j,k), tsum_t%v(i,j,k), &
          tsum_t%w(i,j,k), tsum_t%u(i,j,k), tsum_t%v(i,j,k), &
          tsum_t%w(i,j,k), tsum_t%u2(i,j,k), tsum_t%v2(i,j,k), &
          tsum_t%w2(i,j,k), tsum_t%uw(i,j,k), tsum_t%vw(i,j,k), &
          tsum_t%uv(i,j,k), tsum_t%dudz(i,j,k)
      enddo
    enddo
  enddo
  close(7)
endif


return
end subroutine load_data
