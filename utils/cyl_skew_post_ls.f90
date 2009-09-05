program cylinder_skew_post
use types, only : rprec
use cylinder_skew_base_ls, only : skew_angle, d, l, scale_fact, ntrunk,ngen,nproc
implicit none

!  in thousands
integer, parameter :: iter_start=10; 
integer, parameter :: iter_step=10;
integer, parameter :: iter_end=10; 
integer, parameter :: niter=(iter_end - iter_start)/iter_step + 1

character(200) :: fdir, fname, temp
integer :: iter, iter_count, n,ng, np, nsamples, nsamples_tot, nstart
real(rprec) :: Ap, CD_avg, fD_avg, Uinf_avg, Ap_tot, CD_tot, fD_tot, Uinf_tot
real(rprec), dimension(:), allocatable :: CD, fD, Uinf
real(rprec), dimension(:,:), allocatable :: dat
logical :: exst

!  Check that all directories are present
do iter=iter_start,iter_end,iter_step
  fdir = 'output.'
  write (temp, '(i0)') iter
  fdir = trim (fdir) // temp
  fdir = trim(fdir) // 'k'
  fdir = trim(fdir)
  inquire(file=fdir, exist=exst)
  if(.not. exst) then
    write(*,*) 'Directory ', trim(fdir), ' not found. Please correct iter settings!'
  !  stop
  endif
enddo

Ap_tot = 0._rprec
CD_tot = 0._rprec
fD_tot = 0._rprec

!  Open output file
fname ='cylinder_skew_ls_CD.dat'
open (unit = 11,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')

!  Loop over each generation
do ng=1,ngen
  Ap = (scale_fact**2)**(ng-1) * d * (l*cos(skew_angle)) * (ntrunk**ng) ! Projected area for generation

  iter_count=0
  do iter=iter_start,iter_end,iter_step
    iter_count = iter_count + 1
    fdir = 'output.'
    write (temp, '(i0)') iter
    fdir = trim (fdir) // temp
    fdir = trim(fdir) // 'k'
    
    write(*,*) 'Checking for files in directory : ', fdir
    do np=0,nproc-1 ! Sum over all processors

      fname = trim(fdir) // '/cylinder_skew_ls_CD.dat.g'
      write (temp, '(i0)') ng
      fname = trim (fname) // temp
      fname = trim(fname) // '.c'
      write (temp, '(i0)') np
      fname = trim (fname) // temp
      write(*,*) 'Checking if file exists : ', fname
      
      inquire (file=fname, exist=exst)
      if(exst) then ! Generation is associated with proc
        write(*,*) 'File exists : ', fname
        call load_data()
        nsamples = size(dat,2); ! size(dat,2) should be the same for all files
        nstart = nsamples*(iter_count-1) !  Keep data in time order
        write(*,*) 'nsamples : ', nsamples
        write(*,*) 'nstart   : ', nstart
 
        if(.not. allocated(CD)) then
          write(*,*) 'Allocating CD and fD'
          nsamples_tot = niter*nsamples
          allocate(CD(nsamples_tot))
          allocate(fD(nsamples_tot))
          allocate(Uinf(nsamples_tot))
          CD=0._rprec;
          fD=0._rprec;
          Uinf=0._rprec
        endif
        
        !  Check that nsamples in current file is the same
        if(niter*nsamples /= ubound(CD,1)) then
          write(*,*) 
          write(*,*) 'Error: mismatch in data size in ', fname, '!'
          stop
        endif

        do n=1,nsamples
          CD(nstart + n) = CD(nstart + n) + dat(2,n) ! Sum proc contributions
          fD(nstart + n) = fD(nstart + n) + dat(3,n) ! Sum proc contributions
          Uinf(nstart + n) = dat(4,n) ! Uinf for all procs is the same; just over writting
        enddo
        write(*,*) 'Summed proc', np, ' contribution for gen ', ng
        deallocate(dat)
      endif        
    enddo
  enddo
  CD_avg = sum(CD)/nsamples_tot/Ap
  fD_avg = sum(fD)/nsamples_tot
  Uinf_avg = sum(Uinf)/nsamples_tot
  write(11,*) 'Gen, Ap, CD, fD, Uinf : ', ng, Ap, CD_avg, fD_avg, Uinf_avg

  CD_tot = CD_tot + Ap*CD_avg
  fD_tot = fD_tot + fD_avg
  Ap_tot = Ap_tot + Ap

  if(ng==1) then
  !  Open output file
    fname ='cylinder_skew_ls_CD_inst.dat'
    open (unit = 12,file = fname, status='unknown',form='formatted', &
      action='write',position='rewind')
    do n=1,nsamples_tot
      write(12,*) CD(n)/Ap
    enddo
    close(12)
  endif

  deallocate(CD)
  deallocate(fD)
  deallocate(Uinf)

enddo
CD_tot = CD_tot/Ap_tot
!  Just using the last Uinf_avg; it is a global quantity so Uinf_avg is the same for all generations
write(11,*) 'Gen_tot, Ap_tot, CD_tot, fD_tot, Uinf_tot : ',ngen, Ap_tot, CD_tot, fD_tot, Uinf_avg
close(11)

stop

contains

!**********************************************************************
subroutine load_data()
!**********************************************************************
implicit none

integer :: ios, n, nlines

open (unit = 10,file = fname, status='old',form='formatted', &
  action='read',position='rewind')

nlines = -3 ! First three lines are auxilary data
do
  read(10,'(a)',iostat=ios,end=100)
  nlines=nlines+1
enddo

100 close(10)
allocate(dat(4,nlines))

open(unit = 10,file=fname,action='read',status='old',form='formatted', &
  position='rewind')

read(10,*) temp
read(10,*) temp
read(10,*) temp

do n=1,nlines
  read(10,fmt=*,iostat=ios) dat(1,n), dat(2,n), dat(3,n), dat(4,n)
enddo
close(10)

return
end subroutine load_data

end program cylinder_skew_post
