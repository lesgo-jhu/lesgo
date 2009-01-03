!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interp_types
integer, parameter :: rp_s = kind (1.d0)
integer, parameter :: rp_b = kind (1.d0)
end module interp_types

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--interpolate initial conditions-in physical space
!--for now, handle MPI files by reading all files one large array
!  will change later if it is a problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program interp
use interp_types
implicit none

character (*), parameter :: path='./'
character (*), parameter :: fbase = path // 'vel.out'
character (*), parameter :: fsave = path // 'save.vel.out'
character (*), parameter :: MPI_suffix = '.c'  !--must be proc number

logical, parameter :: DEBUG = .true.
logical, parameter :: MPI_s = .false.
logical, parameter :: MPI_b = .false.

integer, parameter :: np_s = 1  !--must be 1 when MPI_s is false
integer, parameter :: np_b = 1  !--must to 1 when MPI_b is false

!--MPI: these are the total sizes (include all processes)
integer, parameter :: nx_s = 64, ny_s = 64, nz_s = 65
integer, parameter :: nx_b = 32, ny_b = 32 , nz_b = 33

character (64) :: fmt
character (128) :: fname

integer :: ip
integer :: lbz, ubz
integer :: i, j, k

!--save is here, b/c xlf likes to make these automatic!
real (rp_s), save, dimension (nx_s+2, ny_s, nz_s) :: u_s, v_s, w_s,       &
                                                     Rx_s, Ry_s , Rz_s,   &
                                                     cs_s, FLM_s, FMM_s,  &
                                                     FQN_s, FNN_s
real (rp_b), save, dimension (nx_b+2, ny_b, nz_b) :: u_b, v_b, w_b,       &
                                                     Rx_b, Ry_b ,Rz_b,    &
                                                     cs_b, FLM_b, FMM_b,  &
                                                     FQN_b, FNN_b
real (rp_s) :: ke1_s (ny_s, nz_s)  !--ke summed in xdir
real (rp_b) :: ke1_b (ny_b, nz_b)
real (rp_s) :: ke_s
real (rp_b) :: ke_b

!---------------------------------------------------------------------

if ((.not. MPI_s) .and. (np_s /= 1)) then
  write (*, *) 'when MPI_s is false, must have np_s = 1'
  stop
end if

if ((.not. MPI_b) .and. (np_b /= 1)) then
  write (*, *) 'when MPI_b is false, must have np_b = 1'
  stop
end if

write (*, '(1x,a)') 'Going to interpolate "' // fbase // '"'

fmt = '(1x,5(a,i0))'
write (*, fmt) 'start size is (', nx_s, ' X ', ny_s,' X ', nz_s, ') / ',  &
               np_s, ' process(es), kind = ', rp_s
write (*, fmt) 'new size is (', nx_b, ' X ', ny_b,' X ', nz_b, ') / ',  &
               np_b, ' process(es), kind = ', rp_b

!--can save ram by doing one array at a time, but this requires changing
!  the output format to the following
!  write (1) array1
!  write (1) array2
!  etc.
if (.not. MPI_s) then 

  open (1, file=fbase, form='unformatted')
  read (1) u_s, v_s, w_s, Rx_s, Ry_s, Rz_s, cs_s, FLM_s, FMM_s, FQN_s, FNN_s
  close (1)

  !--save a copy
  !--need 'system' command to copy fbase to fsave more efficiently?
  open (1, file=fsave, form='unformatted')
  write (1) u_s, v_s, w_s, Rx_s, Ry_s, Rz_s, cs_s, FLM_s, FMM_s, FQN_s, FNN_s
  close (1)

else

  do ip = 0, np_s-1

    write (fname, '(a,a,i0)') trim (fbase), MPI_suffix, ip
    open (1, file=fname, form='unformatted')

    lbz = ip * (nz_s-1)/np_s + 1
    ubz = lbz + (nz_s-1)/np_s   
    read (1) u_s(:, :, lbz:ubz), v_s(:, :, lbz:ubz), w_s(:, :, lbz:ubz),     &
             Rx_s(:, :, lbz:ubz), Ry_s(:, :, lbz:ubz), Rz_s(:, :, lbz:ubz),  &
             cs_s(:, :, lbz:ubz), FLM_s(:, :, lbz:ubz),                      &
             FMM_s(:, :, lbz:ubz), FQN_s(:, :, lbz:ubz), FNN_s(:, :, lbz:ubz)

    close (1)

    !--save a copy
    write (fname, '(a,a,i0)') trim (fsave), MPI_suffix, ip
    open (1, file=fname, form='unformatted')

    lbz = ip * (nz_s-1)/np_s + 1
    ubz = lbz + (nz_s-1)/np_s   
    write (1) u_s(:, :, lbz:ubz), v_s(:, :, lbz:ubz), w_s(:, :, lbz:ubz),     &
              Rx_s(:, :, lbz:ubz), Ry_s(:, :, lbz:ubz), Rz_s(:, :, lbz:ubz),  &
              cs_s(:, :, lbz:ubz), FLM_s(:, :, lbz:ubz),                      &
              FMM_s(:, :, lbz:ubz), FQN_s(:, :, lbz:ubz), FNN_s(:, :, lbz:ubz)

    close (1)
   
  end do

end if

if (DEBUG) then

  ke1_s = 0._rp_s
  do k = 1, nz_s
    do j = 1, ny_s
      do i = 1, nx_s
        ke1_s(j, k) = ke1_s(j, k) + 0.5_rp_s * (u_s(i, j, k)**2 +     &
                                                v_s(i, j, k)**2 +     &
                                                w_s(i, j, k)**2)
      end do
    end do
  end do

  ke_s = 0._rp_s
  do k = 1, nz_s
    do j = 1, ny_s
      ke_s = ke_s + ke1_s(j, k)
    end do
  end do
  ke_s = ke_s / (nx_s * ny_s * nz_s)

  write (*, *) 'ke_s = ', ke_s

end if

call interpolate (nx_s, ny_s, nz_s, u_s, nx_b, ny_b, nz_b, u_b)
call interpolate (nx_s, ny_s, nz_s, v_s, nx_b, ny_b, nz_b, v_b)
call interpolate (nx_s, ny_s, nz_s, w_s, nx_b, ny_b, nz_b, w_b)
call interpolate (nx_s, ny_s, nz_s, Rx_s, nx_b, ny_b, nz_b, Rx_b)
call interpolate (nx_s, ny_s, nz_s, Ry_s, nx_b, ny_b, nz_b, Ry_b)
call interpolate (nx_s, ny_s, nz_s, Rz_s, nx_b, ny_b, nz_b, Rz_b)
call interpolate (nx_s, ny_s, nz_s, cs_s, nx_b, ny_b, nz_b, cs_b)
call interpolate (nx_s, ny_s, nz_s, FLM_s, nx_b, ny_b, nz_b, FLM_b)
call interpolate (nx_s, ny_s, nz_s, FMM_s, nx_b, ny_b, nz_b, FMM_b)
call interpolate (nx_s, ny_s, nz_s, FQN_s, nx_b, ny_b, nz_b, FQN_b)
call interpolate (nx_s, ny_s, nz_s, FNN_s, nx_b, ny_b, nz_b, FNN_b)

if (.not. MPI_b) then

  open (1, file=fbase, form='unformatted')
  write (1) u_b, v_b, w_b, Rx_b, Ry_b, Rz_b, cs_b, FLM_b, FMM_b, FQN_b, FNN_b
  close (1)

else

  do ip = 0, np_b-1

    write (fname, '(a,a,i0)') trim (fbase), MPI_suffix, ip
    open (1, file=fname, form='unformatted')

    lbz = ip * (nz_b-1)/np_b + 1
    ubz = lbz + (nz_b-1)/np_b  !--corresponds to nz in main code  
    write (1) u_b(:, :, lbz:ubz), v_b(:, :, lbz:ubz), w_b(:, :, lbz:ubz),     &
              Rx_b(:, :, lbz:ubz), Ry_b(:, :, lbz:ubz), Rz_b(:, :, lbz:ubz),  &
              cs_b(:, :, lbz:ubz), FLM_b(:, :, lbz:ubz),                      &
              FMM_b(:, :, lbz:ubz), FQN_b(:, :, lbz:ubz), FNN_b(:, :, lbz:ubz)

    close (1)

  end do

end if

if (DEBUG) then

  ke1_b = 0._rp_b
  do k = 1, nz_b
    do j = 1, ny_b
      do i = 1, nx_b
        ke1_b(j, k) = ke1_b(j, k) + 0.5_rp_b * (u_b(i, j, k)**2 +     &
                                                v_b(i, j, k)**2 +     &
                                                w_b(i, j, k)**2)
      end do
    end do
  end do

  ke_b = 0._rp_b
  do k = 1, nz_b
    do j = 1, ny_b
      ke_b = ke_b + ke1_b(j, k)
    end do
  end do
  ke_b = ke_b / (nx_b * ny_b * nz_b)

  write (*, *) 'ke_b = ', ke_b

  write (*, *) 'Rx_s@nz_s-1 = ', Rx_s(nx_s, ny_s, nz_s-1)
  write (*, *) 'Rx_s@nz_s = ', Rx_s(nx_s, ny_s, nz_s)
  write (*, *) 'Rx_b@nz_b-1 = ', Rx_b(nx_b, ny_b, nz_b-1)

end if

end program interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--if run into stack size problems, then make a module and make 
!  nx, ny, nz visible
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interpolate (nx_s, ny_s, nz_s, u_s, nx_b, ny_b, nz_b, u_b)
use interp_types
implicit none

integer, intent (in) :: nx_s, ny_s, nz_s, nx_b, ny_b, nz_b
real (rp_s), dimension (nx_s+2, ny_s, nz_s), intent (in) :: u_s
real (rp_b), dimension (nx_b+2, ny_b, nz_b), intent (out) :: u_b
real (rp_b) :: rx, ry, rz, P1, P2, P3, P4, P5, P6
real (rp_b) :: dx, dy, dz

integer :: jx, jy, jz, i, j, k, i_wrap, j_wrap

rx = real (nx_s, rp_b) / real (nx_b, rp_b)
ry = real (ny_s, rp_b) / real (ny_b, rp_b)
rz = real (nz_s-1, rp_b) / real(nz_b-1, rp_b)

!do jz = 1, nz_b-1  ! we already know the last point!
!$omp parallel do  &
!$omp private(i, j, k, i_wrap, j_wrap, dx, dy, dz, P1, P2, P3, P4, P5, P6)
do jz = 1, nz_b

  if (jz >= floor ((nz_s - 2) / rz) + 1) then
    !--do not use pts near top because of interpolations involving
    !  garbage (BOGUS) at top
    u_b(1:nx_b, 1:ny_b, jz) = u_b(1:nx_b, 1:ny_b, jz - 1)

  else

    k = floor ((jz-1) * rz) + 1
    dz = (jz-1) * rz - (k-1)

    do jy = 1, ny_b

      j = floor ((jy-1) * ry) + 1
      j_wrap = modulo (j+1-1, ny_s) + 1
      dy = (jy-1) * ry - (j-1)

      do jx = 1, nx_b

        i = floor ((jx-1) * rx) + 1
        i_wrap = modulo (i+1-1, nx_s) + 1
        dx = (jx-1) * rx - (i-1)

        P1 = (1._rp_b - dx) * u_s(i, j, k) + dx * u_s(i_wrap, j, k)
        P2 = (1._rp_b - dx) * u_s(i, j, k+1) + dx * u_s(i_wrap, j, k+1)
        P3 = (1._rp_b - dx) * u_s(i, j_wrap, k) + dx * u_s(i_wrap, j_wrap, k)
        P4 = (1._rp_b - dx) * u_s(i, j_wrap, k+1) +  &
             dx * u_s(i_wrap, j_wrap, k+1)
        P5 = (1._rp_b - dy) * P1 + dy * P3 
        P6 = (1._rp_b - dy) * P2 + dy * P4

        u_b(jx, jy, jz) = (1._rp_b - dz) * P5 + dz * P6
        
      end do

    end do

  end if

end do
!$omp end parallel do

! still have to do nz_b level
!do jy = 1, ny_b
!  j = floor ((jy-1) * ry) + 1
!  j_wrap = modulo (j+1-1, ny_s) + 1
!  dy = (jy-1) * ry - (j-1)
!  do jx = 1, nx_b
!    i = 1 + floor ((jx-1) * rx)
!    i_wrap = modulo (i+1-1, nx_s) + 1
!    dx = (jx-1) * rx - (i-1)
!
!    P1 = (1._rp_b - dx) * u_s(i, j, nz_s) + dx *u_s(i_wrap, j, nz_s)
!    P2 = (1._rp_b - dx) * u_s(i, j_wrap, nz_s) + dx * u_s(i_wrap, j_wrap, nz_s)
!
!    u_b(jx, jy, nz_b) = (1._rp_b - dy) * P1 + dy * P2
!
!  end do
!
!end do

u_b(nx_b+1:nx_b+2,:,:) = 0._rp_b

end subroutine interpolate
