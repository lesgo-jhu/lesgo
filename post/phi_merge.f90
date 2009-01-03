program phi_merge
implicit none

integer, parameter :: rp = kind (1.d0)

integer, parameter :: nx=64, ny=64, nz=65
integer, parameter :: lh=nx/2+1, ld=2*lh

integer, parameter :: np = 32  !--number of files to merge
character (*), parameter :: MPI_suffix = '.c'  !--c for coord, then a number

character (128) :: fname
character (128) :: MPI_fname(0:np-1)

integer :: ip
integer :: i, j, k
integer :: lbz, ubz

logical :: exst

real (rp) :: phi(ld, ny, 0:nz)  !--level set function

!---------------------------------------------------------------------

!--read in level-set function, if it is there, else put 0
write (fname, '(a)') 'phi.out'

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
    write (*, *) 'Error: all or part of phi missing'
    stop
  end if

end do

open (1, file='phi.dat', action='write', position='rewind')
write (1, '(a)') 'variables = "i" "j" "k" "phi"'
write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny, ', k=', nz

do k = 1, nz
  do j = 1, ny
    do i = 1, nx
    
      write (1, '(3(i0,1x),es13.6)') i, j, k, phi(i, j, k)
      
    end do
  end do
end do

close (1)

end program phi_merge
