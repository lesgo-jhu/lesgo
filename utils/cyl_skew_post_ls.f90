program cylinder_skew_post
use types, only : rprec
use cylinder_skew_base_ls, only : d, l, scale_fact, ntrunk,ngen,nproc
implicit none

!  in thousands
integer, parameter :: iter_start=100; 
integer, parameter :: iter_step=10;
integer, parameter :: iter_end=170; 
integer, parameter :: niter=(iter_end - iter_start)/iter_step + 1

character(200) :: fdir, fname, temp
integer :: iter, n,ng, np, nsamples, nstart
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
  inquire(file=fname, exist=exst)
  if(.not. exst) then
    write(*,*) 'Directory', fname, ' not found. Please correct iter settings!'
    stop
  endif
enddo

Ap_tot = 0._rprec
CD_tot = 0._rprec
fD_tot = 0._rprec

!  Loop over each generation
do ng=1,ngen
  Ap = (scale_fact**2)**(ng-1)*d*l*(ntrunk**ng) ! Projected area for generation

  do iter=iter_start,iter_end,iter_step

    fdir = 'output.'
    write (temp, '(i0)') iter
    fdir = trim (fdir) // temp
    fdir = trim(fdir) // 'k'

    do np=0,nproc-1 ! Sum over all processors

      fname = trim(fdir) // '/cylinder_skew_CD.dat.g'
      write (temp, '(i0)') ng
      fname = trim (fname) // temp
      fname = trim(fname) // '.c'
      write (temp, '(i0)') np
      fname = trim (fname) // temp
      
      inquire (file=fname, exist=exst)
      if(exst) then ! Generation is associated with proc
        write(*,*) fname
        call load_data()
        if(.not. allocated(CD)) then
          allocate(CD(niter*size(dat,2)))
          allocate(fD(niter*size(dat,2)))
          CD=0._rprec;
          fD=0._rprec;
        endif
        nsamples = size(dat,2); ! size(1,dat) is same for all files
        nstart = nsamples*(iter-1) !  Keep data in time order
        do n=1,nsamples
          CD(nstart + n) = CD(nstart + n) + dat(2,n) ! Sum proc contributions
          fD(nstart + n) = fD(nstart + n) + dat(3,n) ! Sum proc contributions
        enddo
        deallocate(dat)
      endif        
    enddo
  enddo
  write(*,*) 'im here'
  fname ='cylinder_skew_CD.dat.g'
  write (temp, '(i0)') ng
  fname = trim (fname) // temp
  open (unit = 11,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')
  write(11,*) 'Ap, CD, fD : ', Ap, sum(CD)/size(CD,1)/Ap, sum(fD)/size(fD,1)
  close(11)
  write(*,*) 'Gen #, Ap, CD, fD : ',ng, Ap, sum(CD)/size(CD,1)/Ap, sum(fD)/size(fD,1)

  CD_tot = CD_tot + sum(CD)/size(CD,1)/Ap
  fD_tot = fD_tot + sum(fD)/size(fD,1)
  Ap_tot = Ap_tot + Ap
 
  deallocate(CD)
  deallocate(fD)

enddo

write(*,*) '# Gen, Ap_tot, CD_tot, fD_tot : ',ngen, Ap_tot, CD_tot, fD_tot

stop

contains

!**********************************************************************
subroutine load_data()
!**********************************************************************
implicit none

integer :: ios, n, nlines

open (unit = 10,file = fname, status='old',form='formatted', &
  action='read',position='rewind')

nlines = 0
do
  read(10,'(a)',iostat=ios,end=100)
  nlines=nlines+1
enddo

100 close(10)
allocate(dat(4,nlines))

open(unit = 10,file=fname,action='read',status='old',form='formatted', &
  position='rewind')
do n=1,nlines
  read(10,fmt=*,iostat=ios) dat(1,n), dat(2,n), dat(3,n), dat(4,n)
enddo
close(10)

return
end subroutine load_data

end program cylinder_skew_post
