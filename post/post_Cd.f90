!--command line arg stuff: intel compiler needs -Vaxlib to compile
!--ifc 7.1 also needs -nbs, because it has retarded defaults
program post_Cd
implicit none

character (*), parameter :: prog_name = 'post_Cd'
!character (*), parameter :: file_in = 'Cd_dynamic.dat'
!character (*), parameter :: file_out = 'Cd_dynamic_plot.dat'
character (*), parameter :: TeX_table_out = 'table.tex'

integer, parameter :: rp = kind (1.e0)
integer, parameter :: jt_per_call = 5
integer, parameter :: n_zone = 1

integer, parameter :: BUFF_LEN=128

!real (rp), parameter :: dt = 0.0002_rp

character (32) :: fmt
character (BUFF_LEN) :: buff
character (128) :: file_in
character (128) :: file_out  ! if not outfile specified, then append .plot

integer, external :: iargc

integer :: status
integer :: i
integer :: num_args
integer :: counter
integer :: n_calls
integer :: iost
integer :: z

logical :: f_flag, o_flag, t_flag

real (rp) :: Cd(n_zone), Cd_avg(n_zone), Cd_rms(n_zone)
real (rp) :: t
real (rp) :: dt

!---------------------------------------------------------------------

! process command line options
f_flag = .false.
o_flag = .false.
t_flag = .false.

num_args = iargc ()

write (*, *) prog_name // ': num_args = ', num_args

if (num_args <= 1) then
  write (*, *) 'usage: post_Cd -f<file> -t<dt>'
  stop
end if

call getarg (0, buff)
!call getarg (0, buff, status)
write (*, *) 'i, buff = ', 0, buff
!write (*, *) 'status = ', status

! skip command itself (i=0)
do i = 1, num_args

  call getarg (i, buff)
  !call getarg (i, buff, status)

  write (*, *) 'i, buff = ', i, buff
  !write (*, *) 'status = ', status

  !if (status < 0) then
  !  write (*, *) 'error processing command line arguments'
  !  stop
  !else if (status > BUFF_LEN) then
  !  write (*, '(a,i0,a)') ' argument ', i, ' is too long for buffer'
  !  write (*, *) 'buff = ', buff
  !  stop
  !end if

  select case (buff(1:2))
    case ('-t')

      if (t_flag) then
        write (*, *) 'error: more than one -t'
        stop
      end if
   
      t_flag = .true.
      read (buff(3:), *) dt
      !read (buff(3:status), *) dt
      
    case ('-f')
      
      if (f_flag) then
        write (*, *) 'error: more than one -f'
        stop
      end if
      
      f_flag = .true.
      read (buff(3:), *) file_in
      !read (buff(3:status), *) file_in

    case ('-o')

      if (o_flag) then
        write (*, *) 'error: more than one -o'
        stop
      end if

      o_flag = .true.
      read (buff(3:), *) file_out
      !read (buff(3:status), *) file_out

    case default
      write (*, *) 'error: invalid arguments'
      stop
  end select
  
end do

if ((.not. t_flag) .or. (.not. f_flag)) then
  write (*, *) 'need to specify -f & -t options'
  stop
end if

if (.not. o_flag) then
  file_out = trim (file_in) // '.plot'
  write (*, *) 'file_out = ', file_out
end if

! provide idiot info
write (*, *) prog_name // ': n_zone = ', n_zone
write (*, *) prog_name // ': jt_per_call = ', jt_per_call
write (*, *) prog_name // ': dt = ', dt

open (1, file = trim (file_in), status='old')

! calculate the avg
Cd_avg = 0._rp
counter = 0

avg_loop: do

  read (1, *, iostat=iost) n_calls, Cd
  if (iost /= 0) exit avg_loop

  Cd_avg = Cd_avg + Cd
  counter = counter + 1

end do avg_loop

if (counter > 0) then

  Cd_avg = Cd_avg / counter

  fmt = '(1x, a, i0, a, es12.5)'
  do z = 1, n_zone
    write (*, fmt) prog_name // ': Cd_avg(zone=', z, ') = ', Cd_avg(z)
  end do

  write (*, *) prog_name // ': counter = ', counter

else
  write (*, *) prog_name // ': error, counter = ', counter
  stop
end if

! calculate the rms
rewind (1)

Cd_rms = 0._rp

rms_loop: do

  read (1, *, iostat=iost) n_calls, Cd
  if (iost /= 0) exit rms_loop

  Cd_rms = Cd_rms + (Cd - Cd_avg)**2

end do rms_loop

! we already know counter > 0 from above
Cd_rms = sqrt (Cd_rms / counter)

do z = 1, n_zone
  write (*, fmt) prog_name // ': Cd_rms(zone=', z, ') = ', Cd_rms(z)
end do

! now output a file suitable for plotting: scale the time to be simulations
! nondimensional time units, instead of time steps.
! actually this is just a quickie--we'll add this to the permanent code
! so this part is no longer required

rewind (1)
open (21, file = trim (file_out))

write (fmt, *) '(', n_zone + 1, '(1x, es12.5))'

rescale_loop: do

  read (1, *, iostat = iost) n_calls, Cd
  if (iost /= 0) exit rescale_loop

  t = n_calls * (dt * jt_per_call)
  write (21, fmt) t, Cd

end do rescale_loop

close (21)
close (1)

! write latex table:
! --------------------------------
! | zone | C_D^{avg} | C_D^{rms} |
! --------------------------------
! |   1  |   XXX     |    xxx    |
! ...
! --------------------------------
open (21, file = TeX_table_out)

write (21, '(a)') '\begin{tabular}{|c|c|c|}'
write (21, '(a)') '  \hline'
write (21, '(a)') '  zone & $ \overline{C}_D $ & $ C_D'' $ \\ \hline'

fmt = '(a,i0,2(a,es12.5),a)'
do z = 1, n_zone
  write (21, fmt) '  ', z, ' & ', Cd_avg(z), ' & ', Cd_rms(z), ' \\ \hline'
end do

write (21, '(a)') '\end{tabular}'

close (21)

end program post_Cd
