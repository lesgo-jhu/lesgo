program debug_merge
implicit none

character (*), parameter :: conf_file = 'debug_merge.conf'

integer, parameter :: np = 4
integer, parameter :: BUF_SIZE = 256

character (BUF_SIZE) :: line
character (BUF_SIZE) :: base_file, path
character (32) :: suffix

integer :: overlap
integer :: mx, my, mz, mmz, mz_tot
integer :: i, j, k, k_tot
integer :: ip
integer :: begin

!---------------------------------------------------------------------

open (1, file=conf_file)
read (1, '(a)') base_file
read (1, '(a)') path
read (1, *) overlap
close (1)

write (*, *) 'base_file =', trim (base_file)
write (*, *) 'path =', trim (path)
write (*, *) 'overlap =', overlap

open (1, file=base_file)
read (1, '(a)') line
if (line(1:10) == '#debug_mod') then
  read (line(11:), *) mx, my, mz_tot
else
  write (*, *) 'invalid base_file format'
  stop
end if
close (1)

!--set mz using np and mz_tot
mz = (mz_tot - 1) / np + 1

open (1, file=trim (path) // trim (base_file) // ".MPI")
              !--this will be the merged file
write (1, '(a)') '# produced by debug_merge'

k_tot = 0
do ip = 0, np-1
  write (suffix, '(".MPI.c",i0)') ip
  open (21, file=trim (path) // trim (base_file) // trim (suffix))
  read (21, *)  !--this skips the #debug_mod line (could use to check np)

  if (ip == np-1) then
    mmz = mz
  else
    mmz = mz - overlap
  end if

  do k = 1, mmz
    k_tot = k_tot + 1  !--this accumulates over separate files
    do j = 1, my
      do i = 1, mx
        read (21, '(a)') line
        begin = scan (line, ':') 
        write (1, '(3(i0,1x),":",a)') i, j, k_tot, trim (line(begin+1:))
      end do
    end do
  end do

  close (21)

end do
end program debug_merge
