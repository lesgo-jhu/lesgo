program cylinder_skew_post
use types, only : rprec
use cylinder_skew_base_ls, only : skew_angle, d, l, scale_fact, ntrunk,ngen,nproc
implicit none

!  in thousands
integer, parameter :: iter_start=390; 
integer, parameter :: iter_step=50;
integer, parameter :: iter_end=540; 
integer, parameter :: niter=(iter_end - iter_start)/iter_step + 1

character(200) :: fdir, fname, temp
integer :: iter, iter_count, n,ng, np, nsamples, nstart
real(rprec) :: Ap,Ap_tot,CD_tot,fD_tot
real(rprec), dimension(:), allocatable :: CD,fD
real(rprec), dimension(:,:), allocatable :: dat
logical :: exst

!  Check that all directories are present
do iter=iter_start,iter_end,iter_step
  fdir = 'output.'
  write (temp, '(i0)') iter
  fdir = trim (fdir) // temp
  fdir = trim(fdir) // 'k'
  inquire(file=fdir, exist=exst)
  if(.not. exst) then
    write(*,*) 'Directory', fname, ' not found. Please correct iter settings!'
    stop
  endif
enddo

Ap_tot = 0._rprec
CD_tot = 0._rprec
fD_tot = 0._rprec

!  Open output file
fname ='cylinder_skew_CD.dat'
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

      fname = trim(fdir) // '/cylinder_skew_CD.dat.g'
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
        if(.not. allocated(CD)) then
          write(*,*) 'Allocating CD and fD'
          allocate(CD(niter*size(dat,2)))
          allocate(fD(niter*size(dat,2)))
          CD=0._rprec;
          fD=0._rprec;
        endif
        nsamples = size(dat,2); ! size(dat,2) is same for all files
        nstart = nsamples*(iter_count-1) !  Keep data in time order
        write(*,*) 'nsamples : ', nsamples
        write(*,*) 'nstart   : ', nstart
        do n=1,nsamples
          CD(nstart + n) = CD(nstart + n) + dat(2,n) ! Sum proc contributions
          fD(nstart + n) = fD(nstart + n) + dat(3,n) ! Sum proc contributions
        enddo
        write(*,*) 'Summed proc', np, ' contribution for gen ', ng
        deallocate(dat)
      endif        
    enddo
  enddo
  write(11,*) 'Gen, Ap, CD, fD : ', ng, Ap, sum(CD)/size(CD,1)/Ap, sum(fD)/size(fD,1)

  CD_tot = CD_tot + sum(CD)/size(CD,1)/Ap
  fD_tot = fD_tot + sum(fD)/size(fD,1)
  Ap_tot = Ap_tot + Ap
 
  deallocate(CD)
  deallocate(fD)

enddo

write(11,*) '# Gen, Ap_tot, CD_tot, fD_tot : ',ngen, Ap_tot, CD_tot, fD_tot
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
