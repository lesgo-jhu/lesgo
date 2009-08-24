module convert_types
implicit none

integer, parameter :: rp = kind (1.d0)
end module convert_types


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--interpolate initial conditions-in physical space
!--for now, handle MPI files by reading all files one large array
!  will change later if it is a problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert_endian
use convert_types
implicit none

character (*), parameter :: path='./'
character (*), parameter :: fbase = path // 'vel.out'
character (*), parameter :: MPI_suffix = '.c'  !--must be proc number

logical, parameter :: DEBUG = .true.
logical, parameter :: MPI = .false.

integer, parameter :: np = 1  !--must be 1 when MPI_s is false
integer, parameter :: iendian = 1 ! 1 - little to big endian; 2 - big to little endian

!--MPI: these are the total sizes (include all processes)
integer, parameter :: nx = 64, ny = 64, nz = 65

character (64) :: fmt
character (128) :: fname
character (64) :: read_endian, write_endian

integer :: ip
integer :: lbz, ubz
integer :: i, j, k

!--save is here, b/c xlf likes to make these automatic!
real (rp), save, dimension (nx+2, ny, nz) :: u, v, w,       &
                                                     Rx, Ry , Rz,   &
                                                     cs, FLM, FMM,  &
                                                     FQN, FNN
real (rp) :: ke1 (ny, nz)  !--ke summed in xdir
real (rp) :: ke

!---------------------------------------------------------------------

if ((.not. MPI) .and. (np /= 1)) then
  write (*, *) 'when MPI is false, must have np = 1'
  stop
end if

if(iendian == 1) then
  read_endian = 'little_endian'
  write_endian = 'big_endian'
elseif(iendian == 2) then
  read_endian = 'little_endian'
  write_endian = 'big_endian' 
else
  write(*,*) 'Error: incorrect endian specification.'
  stop
endif

read_endian = trim(adjustl(read_endian))
write_endian = trim(adjustl(write_endian))

write (*, '(1x,a)') 'Going to interpolate "' // fbase // '"'

fmt = '(1x,5(a,i0))'
write (*, fmt) 'start size is (', nx, ' X ', ny,' X ', nz, ') / ',  &
               np, ' process(es), kind = ', rp

!--can save ram by doing one array at a time, but this requires changing
!  the output format to the following
!  write (1) array1
!  write (1) array2
!  etc.
if (.not. MPI) then 

  open (1, file=fbase, form='unformatted', convert=read_endian)
  read (1) u, v, w, Rx, Ry, Rz, cs, FLM, FMM, FQN, FNN
  close (1)

  !--save a copy
  !--need 'system' command to copy fbase to fsave more efficiently?
  open (1, file=fbase, form='unformatted', convert=write_endian)
  write (1) u, v, w, Rx, Ry, Rz, cs, FLM, FMM, FQN, FNN
  close (1)

else

  do ip = 0, np-1

    write (fname, '(a,a,i0)') trim (fbase), MPI_suffix, ip
    open (1, file=fname, form='unformatted')

    lbz = ip * (nz-1)/np + 1
    ubz = lbz + (nz-1)/np   
    read (1) u(:, :, lbz:ubz), v(:, :, lbz:ubz), w(:, :, lbz:ubz),     &
             Rx(:, :, lbz:ubz), Ry(:, :, lbz:ubz), Rz(:, :, lbz:ubz),  &
             cs(:, :, lbz:ubz), FLM(:, :, lbz:ubz),                      &
             FMM(:, :, lbz:ubz), FQN(:, :, lbz:ubz), FNN(:, :, lbz:ubz)

    close (1)

    !--save a copy
    write (fname, '(a,a,i0)') trim (fbase), MPI_suffix, ip
    open (1, file=fname, form='unformatted',convert='big_endian')

    lbz = ip * (nz-1)/np + 1
    ubz = lbz + (nz-1)/np   

    write(1) u(:, :, lbz:ubz), v(:, :, lbz:ubz), w(:, :, lbz:ubz),     &
             Rx(:, :, lbz:ubz), Ry(:, :, lbz:ubz), Rz(:, :, lbz:ubz),  &
             cs(:, :, lbz:ubz), FLM(:, :, lbz:ubz),                      &
             FMM(:, :, lbz:ubz), FQN(:, :, lbz:ubz), FNN(:, :, lbz:ubz)
    close (1)
   
  end do

end if

if (DEBUG) then

  ke1 = 0._rp
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ke1(j, k) = ke1(j, k) + 0.5_rp * (u(i, j, k)**2 +     &
                                                v(i, j, k)**2 +     &
                                                w(i, j, k)**2)
      end do
    end do
  end do

  ke = 0._rp
  do k = 1, nz
    do j = 1, ny
      ke = ke + ke1(j, k)
    end do
  end do
  ke = ke / (nx * ny * nz)

  write (*, *) 'ke = ', ke

end if

stop
end program convert_endian
