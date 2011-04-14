program trees_pre_ls
use types, only : rprec
use param, only : nx, ny, nz, BOGUS, nproc
use trees_base_ls, only : grid_initialize, pt_of_grid
use trees_setup_ls, only : fill_tree_array, sdistfcn_tree_array
use trees_io_ls, only : draw_tree_array
use trees_global_fmask_ls, only : calc_global_fmask_ta
implicit none

!character (*), parameter :: ftrees_conf = 'trees.conf'
character (*), parameter :: fdraw_out = 'draw_tree_array.dat'
character (*), parameter :: fphi_out = 'phi.dat'
character (*), parameter :: fphi_raw_out = 'phi.out'
character (*), parameter :: fbrindex_out = 'brindex.dat'
character (*), parameter :: fbrindex_raw_out = 'brindex.out'

logical, parameter :: do_write_ascii = .true.
logical, parameter :: do_calc_global_fmask = .true.

!--may choose to connect np with nproc in params and
!  MPI_split with $MPI
character (*), parameter :: MPI_suffix = '.c'
integer, parameter :: np = nproc
logical, parameter :: MPI_split = .false.

!integer, parameter :: nz = (nztot - 1) / np +1  !--local nz

character (128) :: fphi_out_MPI, fphi_raw_out_MPI
character (128) :: fbrindex_out_MPI, fbrindex_raw_out_MPI

integer :: ip, ipmin, ipmax
integer :: i, j, k, ktot
integer :: lbz, ubz
integer, allocatable :: brindex(:, :, :)

real (rprec), allocatable :: phi(:, :, :)
real (rprec) :: x, y, z

!---------------------------------------------------------------------

if (MPI_split) then

!   !--prompt user for chunks to process
!   write (*, '(1x,a,2(i0,a))') 'Enter starting chunk (', 0, '..', np-1, '):'
!   read (*, *) ipmin
!   write (*, *) 'read ipmin = ', ipmin
!   write (*, '(1x,a,2(i0,a))') 'Enter ending chunk (', ipmin, '..', np-1, '):'
!   read (*, *) ipmax
!   write (*, *) 'read ipmax = ', ipmax

  ipmin = 0
  ipmax = np - 1

  if  ((ipmin < 0) .or. (ipmin > np - 1)) then
    write (*, *) 'invalid ipmin = ', ipmin
    stop
  else if ((ipmax < ipmin) .or. (ipmax > np - 1)) then
    write (*, *) 'invalid ipmax = ', ipmax
    stop
  end if

else

  if (np /= 1) then
    write (*, *) 'np must be 1 for when MPI_split is false'
    stop
  end if

  ipmin = 0
  ipmax = np - 1
  
end if

!write (*, *) 'reading ' // ftrees_conf
!call read_trees_conf ()

call grid_initialize ()
call fill_tree_array ()
call draw_tree_array (fdraw_out)

allocate ( phi( nx+2, ny, 0:nz ) )
allocate ( brindex( nx+2, ny, nz ) )
    !--nx+2 to be compatible with other stuff
    !--0-level added for MPI stuff

do ip = ipmin, ipmax

  call sdistfcn_tree_array (ip, (/ nx+2, ny, nz /), phi, brindex)

  if (ip == 0) phi(:, :, 0) = BOGUS

  write (*, '(1x,a,i0,a,i0)') 'chunk ', ip, ' of ', np-1

  !write (*, *) 'np = ', np 
  !write (*, *) 'trees_pre: size (phi)= ', size (phi, 1), size (phi, 2),  &
  !             size (phi, 3) - 1  !--(-1) since ignore 0-level

  !--not sure whether to make phi(nz_global) = BOGUS or not, for now
  !  it is left as a valid value
  if (MPI_split) then

    write (fphi_raw_out_MPI, '(a,a,i0)') trim (fphi_raw_out), MPI_suffix, ip
    open (1, file=fphi_raw_out_MPI, form='unformatted')

    !--note some overlap here for local 0, nz levels (less MPI comms later)
    !lbz = ip * (nz - 1) / np       !--0 level (local)
    !ubz = lbz + (nz - 1) / np + 1  !--nz level (local)

    write (1) phi(:, :, 0:nz)  !--each file gets 0:nz_local
      
    close (1)

    !--branch index
    write (fbrindex_raw_out_MPI, '(a,a,i0)') trim (fbrindex_raw_out),  &
                                             MPI_suffix, ip
    open (1, file=fbrindex_raw_out_MPI, form='unformatted')

    !--overlap is different from above: brindex is only 1:nz_local-1
    !lbz = ip * (nz - 1) / np + 1   !--1 level (local)
    !ubz = lbz + (nz - 1) / np - 1  !--nz-1 level (local)

    write (1) brindex(:, :, 1:nz-1)

    close (1)

    if (do_write_ascii) then  !--ascii-files for checking
      !--phi
      write (fphi_out_MPI, '(a,a,i0)') trim (fphi_out), MPI_suffix, ip
      open (1, file=fphi_out_MPI)
      write (1, '(a)') 'variables = "x" "y" "z" "phi"'
      write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny,  &
                             ', k=', nz + 1

      !lbz = ip * (nz - 1) / np       !--0-level
      !ubz = lbz + (nz - 1) / np + 1  !--nz-level

      do k = 0, nz

        ktot = ip * (nz - 1) + k

        do j = 1, ny
          do i = 1, nx

            x = pt_of_grid (i, 1, 1)
            y = pt_of_grid (j, 2, 1)
            z = pt_of_grid (ktot, 3, 1)
      
            !--when ip == 0, this will write out BOGUS for k = 0
            write (1, '(4(1x,es12.5))') x, y, z, phi(i, j, k)
         
          end do
        end do
      end do
  
      close (1)

      !--brindex
      write (fbrindex_out_MPI, '(a,a,i0)') trim (fbrindex_out),  &
                                           MPI_suffix, ip
      open (1, file=fbrindex_out_MPI)
      write (1, '(a)') 'variables = "x" "y" "z" "brindex"'
      write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny,  &
                             ', k=', nz - 1

      !lbz = ip * (nz - 1) / np + 1      !--1 level
      !ubz = lbz + (nz - 1) / np - 1  !--nz-1 level

      do k = 1, nz-1

        ktot = ip * (nz - 1) + k

        do j = 1, ny
          do i = 1, nx
        
            x = pt_of_grid (i, 1, 1)
            y = pt_of_grid (j, 2, 1)
            z = pt_of_grid (ktot, 3, 1)
      
            !--when ip == 0, this will write out BOGUS for k = 0
            write (1, '(3(es12.5,1x),i0)') x, y, z, brindex(i, j, k)
         
          end do
        end do
      end do
  
      close (1)

    end if

  else  !--no MPI_split (ip = 1)

    open (1, file=fphi_raw_out, form='unformatted')
    write (1) phi(:, :, 1:nz)  !--all k = 1:nz are valid here
    close (1)

    open (1, file=fbrindex_raw_out, form='unformatted')
    write (1) brindex(:, :, 1:nz)
    close (1)

    if (do_write_ascii) then
      !--ascii-file
      open (1, file=fphi_out)
      write (1, '(a)') 'variables = "x" "y" "z" "phi"'
      write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny, ', k=', nz

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            x = pt_of_grid (i, 1, 1)
            y = pt_of_grid (j, 2, 1)
            z = pt_of_grid (k, 3, 1)
      
            write (1, '(4(es12.5,1x))') x, y, z, phi(i, j, k)

          end do
        end do
      end do

      close (1)

      open (1, file=fbrindex_out)
      write (1, '(a)') 'variables = "x" "y" "z" "brindex"'
      write (1, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny, ', k=', nz

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            x = pt_of_grid (i, 1, 1)
            y = pt_of_grid (j, 2, 1)
            z = pt_of_grid (k, 3, 1)
      
            write (1, '(3(es12.5,1x),i0)') x, y, z, brindex(i, j, k)

          end do
        end do
      end do

      close (1)
    
    end if

  end if

end do

deallocate ( brindex )
deallocate ( phi )

if (do_calc_global_fmask) then  !--global fmask initialization
    call calc_global_fmask_ta ( np, MPI_split )
        !--this also writes to file, which is why we pass MPI stuff
end if

end program trees_pre_ls
