!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--This version read vel and fxyz files and trees.conf instead of
!  relying upon tree_data files.  This allows to change the definition
!  of the velscale, for example.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program trees_full_apri_ls
use trees_ls, only : trees_ls_calc_init, apriori, phi, brindex,  &
                     brindex_initialized
use param, only : nz, jt, jt_total
use sim_param, only : u, v, w
use immersedbc, only : fx, fy, fz
use messages
implicit none

character (*), parameter :: sub_name = 'trees_full_apri_ls'
character (*), parameter :: path = 'output/'  !--path to vel, fxyz files

character (*), parameter :: MPI_suffix = '.c'
integer, parameter :: np = 1
logical, parameter :: merge_MPI = .false.

character (*), parameter :: brindex_fname = 'brindex.out'
character (*), parameter :: phi_fname = 'phi.out'

character (64) :: MPI_fname(np)
character (64) :: vel_fname
character (64) :: fxyz_fname
character (64) :: sn

integer :: ip
integer :: lbz2, ubz

logical :: exst

!---------------------------------------------------------------------

!--problem! simulation not writing txz, so the effect of wall stress
!  will not be included in the globalCD calc here
call sim_param_init ( 'u, v, w, txz' )

call mesg ( sub_name, 'reading phi' )

if ( merge_MPI ) then

    do ip = 0, np-1
        write (MPI_fname(ip), '(a,a,i0)') trim (phi_fname), MPI_suffix, ip
    end do

    do ip = 0, np-1
    
        inquire (file=MPI_fname(ip), exist=exst)
        if ( exst ) then
            call mesg ( sub_name, 'reading ' // trim ( MPI_fname(ip) ) )
        else
            write (*, '(1x,a)') 'file ' // trim (MPI_fname(ip)) //  &
                                ' does not exist'
            write (*, '(1x,a)') 'stopping here (no averaging performed)'
            stop
        end if

        open (1, file=MPI_fname(ip), form='unformatted')

        !--overlap is 0:nz_local for brindex
        lbz2 = ip * (nz-1) / np       !--0 level (local)
        ubz = lbz2 + (nz-1) / np + 1  !--nz level (local)

        if ( ip == 0 ) lbz2 = 1
            !--no 0-level for non-MPI compilations
            !--yes this is confusing, and may change this

        read (1) phi(:, :, lbz2:ubz)

        close (1)

    end do

else

    inquire (file=phi_fname, exist=exst)
    if ( exst ) then
        call mesg ( sub_name, 'reading ' // trim ( phi_fname ) )
    else
        write (*, '(1x,a)') 'file ' // trim (phi_fname) // ' does not exist'
        write (*, '(1x,a)') 'stopping here (no averaging performed)'
        stop
    end if

    open(1,file=phi_fname,form='unformatted')
    read (1) phi(:, :, 1:nz)  !--make sure consistent with trees_pre_ls
    close(1)

end if

!--need to initialize brindex manually here when merging MPI files
!--this must be done before trees_ls_calc_init, or else brindex_init will
!  be called (which cannot merge MPI files)
call mesg ( sub_name, 'reading brindex' )

if ( merge_MPI ) then

    do ip = 0, np-1
        write (MPI_fname(ip), '(a,a,i0)') trim (brindex_fname), MPI_suffix, ip
    end do

    do ip = 0, np-1
        inquire (file=MPI_fname(ip), exist=exst)
        if (.not. exst) then
            write (*, '(1x,a)') 'file ' // trim (MPI_fname(ip)) //  &
                                ' does not exist'
            write (*, '(1x,a)') 'stopping here (no averaging performed)'
            stop
        end if

        open (1, file=MPI_fname(ip), form='unformatted')

        !--overlap is 1:nz_local-1 for brindex
        lbz2 = ip * (nz-1) / np + 1   !--1 level (local)
        ubz = lbz2 + (nz-1) / np - 1  !--nz-1 level (local)
    
        read (1) brindex(:, :, lbz2:ubz)

        close (1)

    end do

else

    inquire (file=brindex_fname, exist=exst)
    if (.not. exst) then
        write (*, '(1x,a)') 'file ' // trim (brindex_fname) // ' does not exist'
        write (*, '(1x,a)') 'stopping here (no averaging performed)'
        stop
    end if

    open(1,file=brindex_fname,form='unformatted')
    read (1) brindex  !--make sure consistent with trees_pre_ls
    close(1)

end if

brindex_initialized = .true.
    !--to prevent call to brindex_init in trees_calc_init, which cannot
    !  merge MPI files

!--this is roughly where the loop will begin
write ( *, * ) 'enter jt:'
read ( *, * ) jt

jt_total = jt  !--some internals in apriori use jt_total

!--read velocity file(s)
!--make the following into a subroutine taking u, v, w as args?

write ( vel_fname, '(a,a,i6.6,a)' ) path, 'vel', jt, '.out'

if (merge_MPI) then

    do ip = 0, np-1
        write (MPI_fname(ip), '(a,a,i0)') trim (vel_fname), MPI_suffix, ip
    end do

    do ip = 0, np-1
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
            
        lbz2 = ip * (nz-1) / np + 1  !--1 level (local)
        ubz = lbz2 + (nz-1) / np     !--nz level (local)
    
        read (1) u(:, :, lbz2:ubz), v(:, :, lbz2:ubz), w(:, :, lbz2:ubz)
            !--do not read remaining contents of file
            
        close (1)

    end do

else

    inquire (file=vel_fname, exist=exst)
    if (.not. exst) then
        write (*, '(1x,a)') 'file ' // trim (vel_fname) // ' does not exist'
        write (*, '(1x,a)') 'stopping here (no averaging performed)'
        stop
    end if

    open(1,file=vel_fname,form='unformatted')
    read (1) u, v, w
        !--do not read remaining contents of file
    close(1)

end if
    
!--read fxyz files(s)
!--make a suroutine out of this fx, fy, fz as args?
write ( fxyz_fname, '(a,a,i6.6,a)' ) path, 'fxyz.', jt, '.dat'

if (merge_MPI) then

    do ip = 0, np-1
        write (MPI_fname(ip), '(a,a,i0)') trim (fxyz_fname), MPI_suffix, ip
    end do

    do ip = 0, np-1
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
        lbz2 = ip * (nz-1) / np + 1  !--1 level (local)
        ubz = lbz2 + (nz-1) / np     !--nz level (local)
        read (1) fx(:, :, lbz2:ubz), fy(:, :, lbz2:ubz), fz(:, :, lbz2:ubz)
        close (1)

    end do

else

    inquire (file=fxyz_fname, exist=exst)
    if (.not. exst) then
        write (*, '(1x,a)') 'file ' // trim (fxyz_fname) // ' does not exist'
        write (*, '(1x,a)') 'stopping here (no averaging performed)'
        stop
    end if

    open(1,file=fxyz_fname,form='unformatted')
    read (1) fx, fy, fz
    close(1)

end if


!--read trees.conf, set up tree "environment"
call trees_ls_calc_init ()

call apriori ()

end program trees_full_apri_ls
