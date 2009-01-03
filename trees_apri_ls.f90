!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--stand alone program to calculate a-priori CD from tree_data files
!--supposed to be minimally dependent on param so do not have to
!  recompile for every different simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program trees_apri_ls
use types, rp => rprec
!use param, only : BOGUS, dt
use trees_base_ls
use trees_setup_ls, only : fill_tree_array
use trees_io_ls, only : read_ta_data
use trees_fmodel_ls, only : fcoeff_d, fcoeff_dls, write_dls_line,  &
                            fcoeff_nba, write_nba_line
use messages
implicit none

character (*), parameter :: sub_name = 'trees_apri_ls'
character (*), parameter :: fout_dls = 'dls.dat'
character (*), parameter :: fout_nba = 'nba.dat'

integer, parameter :: lunout_dls = 21
integer, parameter :: lunout_nba = 22

logical, parameter :: fcoeff_d_analysis = .false.
logical, parameter :: fcoeff_dls_analysis = .true.
logical, parameter :: fcoeff_nba_analysis = .true.

character (64) :: fmt
character (256) :: fname

integer :: nstart, nstop, nstep
integer :: gen, ngen
integer :: jt
integer :: ios

logical :: exst

real (rp) :: dt
real (rp), allocatable :: CD(:, :)

!---------------------------------------------------------------------

!call grid_initialize ()  !--is this really needed here?
call fill_tree_array ()
tree_array_initialized = .true.

ngen = tree_array(n_tree) % n_gen

open (lunout_dls, file=fout_dls)
open (lunout_nba, file=fout_nba)

if (fcoeff_d_analysis) then
  allocate (CD(nzone, 0:ngen))
end if

write (*, *) 'Enter nstart, nstop, nstep:'

read (*, *) nstart, nstop, nstep

write (*, *) 'nstart = ', nstart

!--do some basic checks on nstart, nstop, nskip
if (nstart > nstop) then
  write (*, *) 'nstart > nstop, nothing to do'
  stop
end if
write (*, *) 'nstop = ', nstop

do jt=nstart,nstop,nstep

  write (*, *) 'jt =', jt
  
  write (fname, '(a,i6.6,a)') 'tree_data.', jt, '.dat'
  inquire (file=fname, exist=exst)
  
  if (.not. exst) then
    write (fname, '(a,i6.6,a)') 'tree_data.', jt, '.dat.c0'
    inquire (file=fname, exist=exst)
    if (.not. exst) then
      write (*, *) 'tree_data file does not exist, jt =', jt
      stop
    end if
  end if
  
  !write (*, *) 'Enter file name (tree_data.XXXXXX.dat[.c0]): '
  ! read (*, *) fname
  ! write (*, *) 'file name = ', trim (fname)

  !--prompt user for time index
  !write (*, *) 'Enter time index (X''s in tree_data.XXXXXX.dat[.c0]): '
  !read (*, *) jt

  !write (fname, '(a,i6.6,a)',iostat=ios) 'tree_data.', jt, '.dat'
  !if (ios /= 0) call error (sub_name, 'i/o problem writing fname')

  !inquire (file=fname, exist=exst)
  !if (.not. exst) call error (sub_name, 'file ' // trim (fname) //   &
  !                            ' nonexistant')

  !!--extract jt from fname
  !read (fname, '("tree_data.",i6)') jt
  !write (*, *) 'jt = ', jt

  call read_ta_data (fname)

  if (fcoeff_d_analysis) then

    !--prompt user for dt
    write (*, *) 'Enter dt: '
    read (*, *) dt
    write (*, *) 'dt = ', dt
  
    do gen = 0, ngen
      call fcoeff_d (gen, nzone, CD(:, gen))
    end do

    write (fmt, '(a,i0,a,i0,a)') '(1x,es13.6,1x,', ngen + 1, '(', nzone,  &
                                 '(es13.6,1x),1x))'
                                                   !--gen starts at 0
    write (*, fmt) (jt - 1) * dt, CD

  end if

  if (fcoeff_dls_analysis) then

    !write (stdout, *) 'jt =', jt

    !--drag, lift, side analysis
    !do gen = 0 , ngen-1  !--no need to do gen
    do gen = 1, 1
      call fcoeff_dls (gen)  !--sets internal tree_variables
      call write_dls_line (jt, gen, lunout_dls)
    end do

  end if

  if (fcoeff_nba_analysis) then

    !write (stdout, *) 'jt =', jt

    !--normal, binormal, axial analysis
    !do gen = 0 , ngen-1  !--no need to do gen
    do gen = 1, 1
      call fcoeff_nba (gen)  !--sets internal tree_variables
      call write_nba_line (jt, gen, lunout_nba)
    end do

  end if

end do

end program trees_apri_ls
