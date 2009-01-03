!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--merges phi and brindex
!--this is mainly just for debugging new version of trees_pre_ls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program merge_phi
use types, rp => rprec
use param, only : nx, ny, nztot => nz, BOGUS
use trees_base_ls, only : grid_initialize, pt_of_grid
implicit none

character (*), parameter :: fphi = 'phi.out'
character (*), parameter :: fbrindex = 'brindex.out'
character (*), parameter :: fphi_ascii = 'phi.dat'
character (*), parameter :: fbrindex_ascii = 'brindex.dat'

character (*), parameter :: MPI_suffix = '.c'
integer, parameter :: np = 4  !--number of files to merge

integer, parameter :: nz = (nztot - 1) / np + 1  !--local nz

character (128) :: fphi_in_MPI
character (128) :: fbrindex_in_MPI

integer :: ip
integer :: i, j, k, ktot
integer :: lbz, ubz
integer :: brindex(nx+2, ny, nztot)

real (rp) :: phi(nx+2, ny, 0:nztot)  !--nx+2 to be compatible with other stuff
                           !--0-level added for MPI stuff
real (rp) :: x, y, z

!---------------------------------------------------------------------

!--read in phi, brindex chunks
do ip = 0, np-1

  write (fphi_in_MPI, '(a,a,i0)') trim (fphi), MPI_suffix, ip
  open (1, file=fphi_in_MPI, form='unformatted')

  !--note some overlap here for local 0, nz levels (less MPI comms later)
  lbz = ip * (nz - 1)  !--0 level (local)
  ubz = lbz + nz       !--nz level (local)

  read (1) phi(:, :, lbz:ubz)  !--each file gets 0:nz_local
      
  close (1)

  write (fbrindex_in_MPI, '(a,a,i0)') trim (fbrindex), MPI_suffix, ip
  open (1, file=fbrindex_in_MPI, form='unformatted')

  !--overlap is different from above: brindex is only 1:nz_local-1
  lbz = ip * (nz - 1) + 1  !--1 level (local)
  ubz = lbz + (nz - 2)     !--nz-1 level (local)

  read (1) brindex(:, :, lbz:ubz)

  close (1)

end do

!--slight problem here: nz level was never written to the MPI_split file,
!  so it cannot be read in
brindex(:, :, nztot) = -1

!--write out as single files, both binary and ascii
!--need pt_of_grid etc to match non-split output
call grid_initialize ()

open (1, file=fphi, form='unformatted')
write (1) phi(:, :, 1:nztot)  !--all k = 1:nz are valid here
close (1)

open (1, file=fbrindex, form='unformatted')
write (1) brindex(:, :, 1:nztot)
close (1)

open (1, file=fphi_ascii)
write (1, '(a)') 'variables = "x" "y" "z" "phi"'
write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny, ', k=', nztot

do k = 1, nztot
  do j = 1, ny
    do i = 1, nx

      x = pt_of_grid (i, 1, 1)
      y = pt_of_grid (j, 2, 1)
      z = pt_of_grid (k, 3, 1)

      write (1, '(4(es12.5,1x))') x, y, z, phi(i, j, k)
      !write (1, '(3(i0,1x),es12.5)') i, j, k, phi(i, j, k)

    end do
  end do
end do

close (1)

open (1, file=fbrindex_ascii)
write (1, '(a)') 'variables = "x" "y" "z" "brindex"'
write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny, ', k=', nztot

do k = 1, nztot
  do j = 1, ny
    do i = 1, nx

      x = pt_of_grid (i, 1, 1)
      y = pt_of_grid (j, 2, 1)
      z = pt_of_grid (k, 3, 1)

      write (1, '(3(es12.5,1x),i0)') x, y, z, brindex(i, j, k)
      !write (1, '(4(i0,1x))') i, j, k, brindex(i, j, k)
      
    end do
  end do
end do

close (1)
    
end program merge_phi
