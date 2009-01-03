!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--stand alone program to postprocess tree_data files
!--supposed to be minimally dependent on param so do not have to
!  recompile for every different simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program trees_post_ls
use types, rp => rprec
!use param, only : BOGUS, dt
use trees_base_ls
use trees_setup_ls, only : fill_tree_array
use trees_io_ls, only : read_ta_data
use trees_post_mod_ls, only : mean_ftot, write_drag_lift, write_CD,  &
                              itime, first_time
use messages
implicit none

character (*), parameter :: sub_name = 'trees_post_ls'
character (*), parameter :: fmt = '(i0,1x,3(es13.6,1x),i0)'

character (256) :: fname, ofname

integer :: gen, ngen
integer :: jt
integer :: ios
integer :: nstart, nstop, nstep
integer :: nfiles
integer :: zone
integer, allocatable :: navg(:, :)  !--number of branches used in average
integer :: gen_dl_min, gen_dl_max

logical :: exst
logical :: first

real (rp) :: dt
real (rp), allocatable :: ftot(:, :, :)  !--mean total force (at one time)
real (rp), allocatable :: avg_ftot(:, :, :)  !--time avg total force

!---------------------------------------------------------------------

!call grid_initialize ()  !--is this really needed here?
call fill_tree_array ()
tree_array_initialized = .true.

ngen = tree_array(n_tree) % n_gen

!allocate (ftot(nd, nzone, 0:ngen))
!allocate (navg(nzone, 0:ngen))
!allocate (avg_ftot(nd, nzone, 0:ngen))
!
!ftot = 0._rp
!avg_ftot = 0._rp

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

!write (*, *) 'Enter min., max. gen for drag-lift analysis:'
!read (*, *) gen_dl_min, gen_dl_max
!
!if ((gen_dl_min < 0) .or. (gen_dl_min > gen_dl_max) .or.  &
!    (gen_dl_max > ngen)) then
!  write (*, *) 'invalid generation range'
!  stop
!end if

do jt = nstart, nstop, nstep

  if (jt == nstart) then
    first_time = .true.
  else
    first_time = .false.
  end if
 
  itime = jt
 
  write (fname, '(a,i6.6,a)') 'tree_data.', jt, '.dat' 

  inquire (file=fname, exist=exst)
  if (.not. exst) then
    write (fname, '(a,i6.6,a)') 'tree_data.', jt, '.dat.c0' 
  end if
  
  inquire (file=fname, exist=exst)
  if (.not. exst) then
    call error (sub_name, 'file ' // trim (fname) // ' nonexistant')
  end if
  
  call read_ta_data (fname)

  !--analyze ftot
  !do gen = 0, ngen
  !  call mean_ftot (gen, nzone, ftot(:, :, gen), navg(:, gen))
  !end do

  !avg_ftot = avg_ftot + ftot
  
  !--analyze lift, drag
  !do gen = gen_dl_min, gen_dl_max
  !  call write_drag_lift (gen)
  !end do

  !--if want averages, then need to average the files with another program
  !  e.g. perl script
  write (ofname, '(a,i6.6,a)') 'resolvedCD-', jt, '.dat'
  call write_CD (ofname, 'resolved')
 
  write (ofname, '(a,i6.6,a)') 'totalCD-', jt, '.dat'
  call write_CD (ofname, 'total')

end do

nfiles = (nstop - nstart + nstep) / nstep  !--watch truncation here
write (*, *) 'nfiles = ', nfiles

!--finalize ftot stuff
!avg_ftot = avg_ftot / nfiles

!--output
!--write navg here so force totals can be recovered
!do gen = 0, ngen
!
!  write (fname, '(a,i0,a)') 'ftot-gen', gen, '.dat'
!  
!  open (1, file=fname, action='write', position='rewind')
!
!  do zone = 1, nzone
!    write (1, fmt) zone, avg_ftot(:, zone, gen), navg(zone, gen)
!  end do
!  
!  close (1)
!  
!end do

end program trees_post_ls
