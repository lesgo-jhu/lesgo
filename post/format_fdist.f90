program format_fdist
implicit none

integer, parameter :: rp = kind (1.d0)

character (128) :: fname

integer :: i, j, k
integer :: z
integer :: nzone
integer :: wksp_size(3)

logical :: exst

real (rp), allocatable :: fdist(:, :, :, :)

!---------------------------------------------------------------------

write (*, *) 'Enter size of fdist array (4 numbers, including zones): '
read (*, *) wksp_size, nzone

if (any (wksp_size <= 0) .or. (nzone <= 0)) then
  write (*, *) 'Error: array dimensions must be positive'
  stop
end if

allocate ( fdist(wksp_size(1), wksp_size(2), wksp_size(3), nzone) )

write (*, *) 'Enter fdist file name to be converted: '
read (*, *) fname

inquire (file=fname, exist=exst)
if (.not. exst) then
  write (*, *) 'Error: ' // trim (fname) // ' does not exist, exiting'
end if

open (1, file=fname, action='read', form='unformatted', position='rewind')
read (1) fdist
close (1)

open (1, file=trim (fname) // '.dat', action='write', position='rewind')
write (1, '(a)') 'variables = "i" "j" "k" "fdist"'

do z = 1, nzone

  write (1, '(3(a,i0))') 'zone, f=point, i=', wksp_size(1),     &
                         ',j=', wksp_size(2), ',k=', wksp_size(3)
                   
  do k = 1, wksp_size(3)
    do j = 1, wksp_size(2)
      do i = 1, wksp_size(1)
      
        write (1, '(3(i0,1x),es13.6)') i, j, k, fdist(i, j, k, z)

      end do
    end do
  end do

end do

close (1)

end program format_fdist
