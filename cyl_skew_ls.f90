!**********************************************************************
module cylinder_skew_ls
!**********************************************************************
use types, only : rprec
use param
use cylinder_skew_base_ls
implicit none

save
private

public :: cylinder_skew_init_ls, cylinder_skew_CD_ls

contains

!**********************************************************************
subroutine cylinder_skew_init_ls()
!**********************************************************************

implicit none

character(64) :: fname, temp
integer :: i,j,k,ng

!  Open file which to write global data
fname = path // 'cylinder_skew_ls_gen.out'
$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

!  Read in cylinder_skew_gen.dat file
open (unit = 2,file = fname, status='old',form='formatted', &
  action='read',position='rewind')

allocate(igen(ngen))
allocate(kbottom_inside(ngen))
allocate(kbottom(ngen))
allocate(dz_bottom(ngen))
allocate(ktop_inside(ngen))
allocate(ktop(ngen))
allocate(dz_top(ngen))
allocate(itype(nx+2,ny,nz))

do ng=1,ngen
  read(2,*) igen(ng), kbottom_inside(ng), kbottom(ng), &
    dz_bottom(ng), ktop_inside(ng), ktop(ng), &
    dz_top(ng)
enddo
close(2)

!  Check 1st generation only for need ground association
if(igen(1) /= -1) then
  !  Open file which to write global data
  fname = path // 'cylinder_skew_ls_point.out'
  $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
  $endif

  !  Read in cylinder_skew_gen.dat file
  open (unit = 2,file = fname, status='old',form='formatted', &
    action='read',position='rewind')
  do k=1,nz
    do j = 1,ny
      do i = 1,nx+2
        read(2,*) itype(i,j,k)
      enddo
    enddo
  enddo
   
  close(2)
endif

return
end subroutine cylinder_skew_init_ls

!**********************************************************************
subroutine cylinder_skew_CD_ls ()
!**********************************************************************
use immersedbc, only : fx
use sim_param, only : u
use io, only : jt_total
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.cylinder_skew_ls_CD'
character (*), parameter :: fCD_out = 'output/cylinder_skew_ls_CD.dat'
character(64) :: fname, temp

integer, parameter :: lun = 991  !--keep open between calls
integer, parameter :: n_calc_CD = 10  !--# t-steps between updates

real (rprec), parameter :: Ap = 1._rprec !--projected area

logical, save, dimension(10) :: file_init=.false. !  May want to change this to allocatable to
                                                  !  match the generation number
logical :: opn, exst

real (rprec) :: CD
real (rprec) :: Uinf   !--velocity scale used in calculation of CD
real (rprec) :: fD     !--drag, lift force
real (rprec) :: Uinf_global

integer :: i,j,k,ng
integer :: kstart, kend

real(rprec) :: dz_start, dz_end
real(rprec) :: dz_p

!---------------------------------------------------------------------

if (modulo (jt, n_calc_CD) /= 0) return  !--do nothing

!  The the global Uinf from the inlet plane; average for proc
Uinf = sum (u(1, :, 1:nz-1)) / (ny * (nz - 1))  !--measure at inflow plane

$if ($MPI)

  !  Sum Uinf from all procs and redestribute
  call mpi_allreduce(Uinf, Uinf_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
  !  Average over all procs; assuming distribution is even
  Uinf_global = Uinf_global / nproc

$else

  Uinf_global = Uinf

$endif

do ng=1,ngen
  if(igen(ng) /= -1) then
    
    fD=0.

    !  Check if bottom is in proc domain
    if(kbottom_inside(ng) == 1) then
      kstart = kbottom(ng)
      dz_start = dz_bottom(ng)
    else
      kstart = 1
      dz_start = dz
    endif

    if(ktop_inside(ng) == 1) then
      kend = ktop(ng)
      dz_end = dz_top(ng)
    else
      kend = nz-1 ! -1 to avoid interprocessor overlap
      dz_end = dz
    endif

     !--(-) since want force ON cylinder
     !--dx*dy*dz is since force is per cell (unit volume)
     !--may want to restrict this sum to points with phi < 0.
    do k=kstart,kend
      if(k==kstart) then
        dz_p = dz_start
      elseif(k==kend) then
        dz_p = dz_end
      else
        dz_p = dz
      endif
      if(ng==1) then !  Want to check with ground association
        do j=1,ny
          do i=1,nx
            if(itype(i,j,k) /= 0) fD = fD - fx(i,j,k) * dx * dy * dz_p
          enddo
        enddo
      else
        fD = fD - sum(fx(1:nx, :, k)) * dx * dy * dz_p
      endif
    enddo

    CD = fD / (0.5_rprec * Ap * Uinf_global**2)

    inquire (lun, exist=exst, opened=opn)

    if (.not. exst) call error (sub_name, 'unit', lun, ' nonexistant')
    if (opn) call error (sub_name, 'unit', lun, ' is already open')

    !  Open file which to write global data
    fname = trim(fCD_out)
    write(temp,'(".g",i0)') ng
    fname = trim (fname) // temp

    $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
    $endif

    open (lun, file=fname, position='append')

    if (.not. file_init(ng)) then  !--set up file for output

    !--write a header
      write (lun, '(a,es12.5)') '# Ap = ', Ap
      write (lun, '(a,es12.5)') '# Gen thickness = ', dz_start+(kend - kstart - 2)*dz+dz_end
      write (lun, '(a)') '# t, CD, fD, Uinf' 

      file_init(ng) = .true.

    end if


    !--output to file
    write (lun, '(4(es12.5,1x))') (jt_total * dt), CD, fD, Uinf_global

    close (lun)  !--only do this to force a flush

  end if

enddo

return
end subroutine cylinder_skew_CD_ls


end module cylinder_skew_ls
