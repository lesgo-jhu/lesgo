program qpost_slice
implicit none

character (*), parameter :: digfmt = 'i6.6'  !--time format in filenames
character (*), parameter :: fnamefmt = '(a,' // digfmt // ',a)'

character (*), parameter :: fmt = '(8(1x,es17.10))'

integer, parameter :: rp = kind (1.d0)

integer,parameter:: nx=64,ny=64,nz=65
                    !--these do not take skips into account

integer, parameter :: nx_skip = nx / 4
integer, parameter :: ny_skip = ny / 4
integer, parameter :: nz_skip = (nz - 1) / 4

integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
integer,parameter:: lh=nx/2+1,ld=2*lh

logical, parameter :: use_rms = .false.

!--merge output files by separate MPI processes
logical, parameter :: merge_MPI = .false.
integer, parameter :: np = 1
character (*), parameter :: MPI_suffix = '.c'  !--c for coord, then a number

real (rp),parameter::pi=3.1415926535897932384626433_rp
real (rp),parameter::L_x=1._rp, L_y= 1._rp
real (rp),parameter::z_i=1._rp,L_z=1._rp * z_i !1._rp - 1._rp / nz
real (rp),parameter::dz=L_z/z_i/(nz-1)
real (rp),parameter::dx=L_x/nx,dy=L_y/ny
real (rp), parameter :: BOGUS = -1234567890._rp

character (128) :: fname, fbase, favgout
character (128) :: MPI_fname(0:np-1)

integer :: ip
integer :: jx,jy,jz,jt
integer :: ios
integer :: vsize(3)
integer :: nstart, nstop, nstep
integer :: nfiles
integer :: lbz, ubz

logical :: exst

real (rp), allocatable, dimension (:, :, :) :: u, v, w
real (rp), allocatable, dimension (:, :, :) :: uavg, vavg, wavg
!--this rms part is experimental
real (rp), allocatable, dimension (:, :, :) :: urms, vrms, wrms
real (rp), allocatable :: coord(:)
real (rp), save :: phi(ld, ny, 0:nz)  !--level set function
                               !--0-level not used for non-MPI
real (rp) :: x, y, z
real (rp) :: u2, ww

!---------------------------------------------------------------------

!--read in level-set function, if it is there, else put 0
write (fname, '(a)') '../phi.out'

if (merge_MPI) then

  do ip = 0, np-1
  
    write (MPI_fname(ip), '(a,a,i0)') trim (fname), MPI_suffix, ip
    
    inquire (file=MPI_fname(ip), exist=exst)
    
    if (exst) then
      write (*, *) 'attempting to read phi for ip = ', ip
      open (1, file=MPI_fname(ip), form='unformatted')
      lbz = ip * (nz-1) / np  !--these will overlap
      ubz = lbz + (nz-1) / np + 1
      read (1) phi(:, :, lbz:ubz)
      close (1)
    else
      write (*, *) 'all or part of phi missing, setting phi to zero'
      phi = 0._rp
      exit
    end if

  end do

else

  inquire (file=fname, exist=exst)
  if (exst) then
    write (*, *) 'attempting to read phi'
    write (*, *) 'expected size=', size (phi, 1), size (phi, 2),  &
                 size (phi, 3) - 1

    open (1, file=fname, form='unformatted')
    read (1) phi(:, :, 1:nz)
    close (1)

    phi(:, :, 0) = BOGUS

  else
    write (*, *) 'phi missing, setting phi to zero'
    phi = 0._rp
  end if

end if

!--read nstart, nstop, nstep from command line
!--best usage would probably be e.g. $> echo "1000 2000 20" | ./qpost  
!--or $> echo "1000 2000 20" > file; ./qpost < file
write (*, *) 'Enter nstart, nstop, nstep:'

read (*, *) nstart, nstop, nstep

write (*, *) 'nstart = ', nstart
write (*, *) 'nstop = ', nstop
write (*, *) 'nstep = ', nstep

!--do some basic checks on nstart, nstop, nskip
if (nstart > nstop) then
  write (*, *) 'nstart > nstop, nothing to do'
  stop
end if

!--process x-slices
write (*, *) 'processing x-slices...'

allocate (u(nx/nx_skip, ny, nz), v(nx/nx_skip, ny, nz),  &
          w(nx/nx_skip, ny, nz))
allocate (uavg(nx/nx_skip, ny, nz), vavg(nx/nx_skip, ny, nz),  &
          wavg(nx/nx_skip, ny, nz))
if (use_rms) then
  allocate (urms(nx/nx_skip, ny, nz), vrms(nx/nx_skip, ny, nz),  &
            wrms(nx/nx_skip, ny, nz))
end if

vsize = (/ nx/nx_skip, ny, nz /)

write (fbase, '(a,i0,a)') 'vel.xskip', nx_skip, '.'
write (favgout, '(a,i0,a)') 'vel-avg.xskip', nx_skip, '.dat'

call process_slices (fbase, favgout, vsize)

if (use_rms) deallocate (wrms, vrms, urms)
deallocate (wavg, vavg, uavg, w, v, u)

!--process y-slices
write (*, *) 'processing y-slices...'

allocate (u(nx, ny/ny_skip, nz), v(nx, ny/ny_skip, nz),  &
          w(nx, ny/ny_skip, nz))
allocate (uavg(nx, ny/ny_skip, nz), vavg(nx, ny/ny_skip, nz),  &
          wavg(nx, ny/ny_skip, nz))
if (use_rms) then
  allocate (urms(nx, ny/ny_skip, nz), vrms(nx, ny/ny_skip, nz),  &
            wrms(nx, ny/ny_skip, nz))
end if

vsize = (/ nx, ny/ny_skip, nz /)

write (fbase, '(a,i0,a)') 'vel.yskip', ny_skip, '.'
write (favgout, '(a,i0,a)') 'vel-avg.yskip', ny_skip, '.dat'

call process_slices (fbase, favgout, vsize)

if (use_rms) deallocate (wrms, vrms, urms)
deallocate (wavg, vavg, uavg, w, v, u)

!--process z-slices
!--this is different from x,y slices since the data for a z-slice is
!  contained all in one file--it is not split across files
write (*, *) 'processing z-slices...'

allocate (u(nx, ny, (nz-1)/nz_skip), v(nx, ny_skip, (nz-1)/nz_skip),  &
          w(nx, ny, (nz-1)/nz_skip))
allocate (uavg(nx, ny, (nz-1)/nz_skip), vavg(nx, ny_skip, (nz-1)/nz_skip),  &
          wavg(nx, ny, (nz-1)/nz_skip))
if (use_rms) then
  allocate (urms(nx, ny, (nz-1)/nz_skip), vrms(nx, ny_skip, (nz-1)/nz_skip),  &
            wrms(nx, ny, (nz-1)/nz_skip))
end if

vsize = (/ nx, ny, (nz-1)/nz_skip /)

write (fbase, '(a,i0,a)') 'vel.zskip', nz_skip, '.'
write (favgout, '(a,i0,a)') 'vel-avg.zskip', nz_skip, '.dat'

call process_slices (fbase, favgout, vsize, .true.)

if (use_rms) deallocate (wrms, vrms, urms)
deallocate (wavg, vavg, uavg, w, v, u)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this doesn't work when zslices is false, i.e. if it false, then it
!  must not be present.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_slices (fbase, favgout, vsize, zslices)
implicit none

character (*), intent (in) :: fbase, favgout
logical, intent (in), optional :: zslices

integer, intent (in) :: vsize(3)

integer :: i, j, k

!---------------------------------------------------------------------

uavg = 0._rp
vavg = 0._rp
wavg = 0._rp

if (use_rms) then
  !--these will holds mean of u^2, v^2, w^2 at first
  urms = 0._rp
  vrms = 0._rp
  wrms = 0._rp
end if

do jt=nstart,nstop,nstep

  write (fname, fnamefmt) trim (fbase), jt, '.out'

  if (merge_MPI) then

    do ip = 0, np-1
      write (MPI_fname(ip), '(a,a,i0)') trim (fname), MPI_suffix, ip
    end do

    do ip = 0, np-1
   
      write (*, *) 'ip =', ip
    
      if (present (zslices)) then
        if (zslices) then
        
          !--kglobal = 1 + m * nz_skip: find m's corresponding to this ip
          lbz = 1 + ceiling ( real (ip * (nz-1) / np, rp) / nz_skip )
          ubz = 1 + floor (real ( ((ip + 1) * (nz-1) / np) - 1, rp ) / nz_skip)
          
          if (lbz <= ubz) then

            inquire (file=MPI_fname(ip), exist=exst)
            if (.not. exst) then
              write (*, '(1x,a)') 'file ' // trim (MPI_fname(ip)) //  &
                                  ' does not exist'
              write (*, '(1x,a)') 'stopping here (no averaging performed)'
              stop
            end if

            write (*, *) 'trying to open ' // trim (MPI_fname(ip))
            write (*, *) 'lbz,ubz =', lbz, ubz

            open (1, file=MPI_fname(ip), form='unformatted')
            read (1) u(:, :, lbz:ubz), v(:, :, lbz:ubz), w(:, :, lbz:ubz)
            close (1)
           
            write (*, *) 'reading complete'

          end if

        end if
      else

        inquire (file=MPI_fname(ip), exist=exst)
        if (.not. exst) then
          write (*, '(1x,a)') 'file ' // trim (MPI_fname(ip)) //  &
                                ' does not exist'
          write (*, '(1x,a)') 'stopping here (no averaging performed)'
          stop
        end if

        open (1, file=MPI_fname(ip), form='unformatted')
      
        !--note the small overlap: gives us an opportunity to check the
        !  overlapping is OK
        !--problem here  with pressure: 0 layer???
        lbz = ip * (nz-1) / np + 1
        ubz = lbz + (nz-1) / np
        read (1) u(:, :, lbz:ubz), v(:, :, lbz:ubz), w(:, :, lbz:ubz)

        close (1)
        
      end if
      
    end do

  else

    inquire (file=fname, exist=exst)
    if (.not. exst) then
      write (*, '(1x,a)') 'file ' // trim (fname) // ' does not exist'
      write (*, '(1x,a)') 'stopping here (no averaging performed)'
      stop
    end if

    open(1,file=fname,form='unformatted')
    read (1) u, v, w
    close(1)

  end if

  uavg = uavg + u
  vavg = vavg + v
  wavg = wavg + w

  if (use_rms) then
    !--calc mean of u^2, v^2, w^2
    urms = urms + u**2
    vrms = vrms + v**2
    wrms = wrms + w**2
  end if

  write (fname, fnamefmt) trim (fbase), jt, '.dat'
  write (*, '(a)') trim (fname)

  call write_vel (fname, vsize, u, v, w, phi, zslices)

end do

nfiles = (nstop - nstart + nstep) / nstep  !--watch truncation here
write (*, *) 'nfiles = ', nfiles

!--normalize all averages
uavg = uavg / nfiles
vavg = vavg / nfiles
wavg = wavg / nfiles

call write_vel (favgout, vsize, uavg, vavg, wavg, phi, zslices)

if (use_rms) then
  !--averages of u^2, v^2, w^2
  urms = urms / nfiles
  vrms = vrms / nfiles
  wrms = wrms / nfiles

  !--calculate rms as u' = sqrt ( E(u^2) - (E(u))^2 )
  do k = 1, vsize(3)
    do j = 1, vsize(2)
      do i = 1, vsize(1)
 
        urms(i, j, k) = sqrt (urms(i, j, k) - uavg(i, j, k)**2)
        vrms(i, j, k) = sqrt (vrms(i, j, k) - vavg(i, j, k)**2)
        wrms(i, j, k) = sqrt (wrms(i, j, k) - wavg(i, j, k)**2)

      end do
    end do
  end do

  call write_vel ('rms.' // favgout, vsize, urms, vrms, wrms, phi, zslices)

end if  

end subroutine process_slices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_vel (fname, isize, u, v, w, phi, zslices)
implicit none

character (*), intent (in) :: fname
integer, intent (in) :: isize(3)
real (rp), intent (in) :: u(isize(1), isize(2), isize(3)),  &
                          v(isize(1), isize(2), isize(3)),  &
                          w(isize(1), isize(2), isize(3))
real (rp), intent (in) :: phi(ld, ny, 0:nz)

logical, intent (in), optional :: zslices

integer :: jx, jy, jz
integer:: jjx, jjy
integer :: skip(3)
integer :: s
integer :: jjz

real (rp) :: x, y, z
real (rp) :: ww
real (rp) :: u2

!---------------------------------------------------------------------

s = 1

if (present (zslices)) then
  if (zslices) s = 0
end if

open (1, file=fname)
write(1, '(a)') 'variables="x" "y" "z" "u" "v" "w" "|u|" "phi"'
write(1, *) 'zone f=point, i=', isize(1), ', j=', isize(2),  &
            ', k=', isize(3)-s

!skip = (/ nx/isize(1), ny/isize(2), 1 /)  !--watch truncation here
skip = (/ nx/isize(1), ny/isize(2), (nz-1)/(isize(3)-s) /)
          !--watch truncation here
     
!--write velocity vectors interpolated to uvp-nodes
!--unless these are z-slices (cannot interpolate)
do jz = 1, isize(3)-s

  if (present (zslices)) then
    write (*, *) 'writing jz =', jz
  end if

  do jy = 1, isize(2)
    do jx = 1, isize(1)

      x=(jx-1)*dx * skip(1)
      y=(jy-1)*dy * skip(2)
      if (present (zslices)) then
        if (zslices) then
          z = (jz-1) * dz * skip(3) + 0.5_rp * dz
        else
          z=(jz-0.5_rp)*dz * skip(3)
        end if
      else
        z=(jz-0.5_rp)*dz * skip(3)
      end if

      if (present (zslices)) then
        if (zslices) then
          ww = w(jx, jy, jz)
        else
          ww = 0.5_rp * (w(jx, jy, jz) + w(jx, jy, jz+1))
        end if
      else
        ww = 0.5_rp * (w(jx, jy, jz) + w(jx, jy, jz+1))
      end if

      u2 = u(jx, jy, jz)**2 + v(jx, jy, jz)**2 + ww**2

      jjx = skip(1) * (jx - 1) + 1
      jjy = skip(2) * (jy - 1) + 1
      jjz = skip(3) * (jz - 1) + 1

      write (1, fmt) x, y, z, u(jx, jy, jz), v(jx, jy, jz), ww,  &
                     sqrt (u2), phi(jjx, jjy, jjz)

    end do
  end do
end do

close(1)

end subroutine write_vel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program qpost_slice
