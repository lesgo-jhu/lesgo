!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--averages data in resolvedCD-XXXXXX.dat and totalCD-XXXXXX.dat files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program mean_CD
implicit none

character (*), parameter :: ofile_suffix = 'avg-CD.dat'

integer, parameter :: wp = kind (1.d0)
integer, parameter :: ndat_max = 256
integer, parameter :: nicol = 2
integer, parameter :: nrcol = 4
integer, parameter :: gen_max = 12

logical, parameter :: sort_into_zones = .true.

character (128) :: fname, ofmt, prefix, ofile

integer :: i
integer :: jt
integer :: ios
integer :: nstart, nstop, nstep
integer :: nfiles
integer :: ndat
integer :: icol(nicol, ndat_max)
integer :: g
integer :: ng_count(0:gen_max)

logical :: ext
logical, save :: firstpass = .true.

real (wp) :: rcol(nrcol, ndat_max)
real (wp) :: avg(nrcol, ndat_max)
real (wp) :: mean(0:gen_max), rms(0:gen_max)  !--only for 1 col (CD)
    
!---------------------------------------------------------------------

write (*, *) 'Enter file prefix (e.g. resolved or total):'
read (*, *) prefix

!--read nstart, nstop, nstep from command line
!--best usage would probably be e.g. $> echo "1000 2000 20" | ./trees_post_ls
!--or $> echo "1000 2000 20" > file; ./trees_post_ls < file
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

avg = 0._wp

do jt = nstart, nstop, nstep

    write (fname, '(a,i6.6,a)') trim (prefix) // 'CD-', jt, '.dat'

    inquire (file=fname, exist=ext)
    if (.not. ext) then
      write (*, *) 'error: file ', fname, 'does not exist'
      stop
    end if

    open (1, file=fname, action='read', position='rewind')

    if (firstpass) then  !--determine how many lines
    
      ndat = 0
      
      do
        read (1, *, iostat=ios)
        if (ios /= 0) exit
        ndat=ndat+1
      end do

      write (*, *) 'ndat = ', ndat

      if (ndat > ndat_max) then
        write (*, *) 'error: ndat > ndat_max'
        stop
      end if
     
      rewind (1)  !--prepare to read in data

      firstpass = .false.

    end if
    
    do i = 1, ndat
      read (1, *) icol(:, i), rcol(:, i)
    end do

    close (1)
    
    avg = avg + rcol

end do

nfiles = (nstop - nstart + nstep) / nstep  !--watch truncation here
write (*, *) 'nfiles = ', nfiles

!--normalize avg
avg = avg / nfiles

!--output
write (ofmt, '(2(a,i0),a)') '(', nicol, '(i0,1x),', nrcol, '(es13.6,1x))'
ofile = trim (prefix) // '.' // ofile_suffix
open (1, file=ofile, action='write', position='rewind')

if (sort_into_zones) then
  !--sort data based on column 2
  !--this is ugly

  !--also calculate mean, rms based on generation
  mean = 0.0_wp
  rms = 0.0_wp
  ng_count = 0
  
  do g = minval (icol(2, :)), maxval (icol(2, :))
    
    if (g > gen_max .or. g < 0) then
      write (*, *) 'check data file, g out of range 0..gen_max'
      stop
    end if
    
    if (count (icol(2, :) == g) == 0) cycle  !--nothing to write here

    write (1, '(a)') 'zone, f=point'  !--tecplot zone separator

    do i = 1, ndat
      if (icol(2, i) == g ) then
        write (1, ofmt) icol(:, i), avg(:, i)
        mean(g) = mean(g) + avg(1, i)  !--avg(1, i) holds CD
        rms(g) = rms(g) + avg(1, i)**2  !--tmp storage for sum(CD^2)
        ng_count(g) = ng_count(g) + 1
      end if
    end do
    
  end do

  where (ng_count /= 0)
    mean = mean / ng_count
    rms = sqrt ( max (0.0_wp, rms / ng_count - mean**2) )
  !elsewhere
  !  mean = 0.0_wp
  !  rms = 0.0_wp  
  end where

  !--write mean, rms
  open (21, file=trim (prefix) // '.sorted-mean-rms.' // ofile_suffix)
  write (21, '(a)') 'zone, f=point'
  do g = 0, gen_max
    if ( ng_count(g) /= 0 ) then
      write (21, '(i0,2(1x,es13.6))') g, mean(g), rms(g)
    end if
  end do
  close (21)
    
else

  do i = 1, ndat
    write (1, ofmt) icol(:, i), avg(:, i)
  end do

end if

close (1)

end program mean_CD
