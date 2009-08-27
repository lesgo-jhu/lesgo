program cylinder_skew_post
implicit none
integer, parameter :: DP = kind(0.0D0)

!  Base dimensions for projected area
real(DP), parameter :: d = 28.8_DP*4._DP/185._DP
real(DP), parameter :: l = 50.4_DP*4._DP/185._DP
real(DP), parameter :: scale_fact=0.5_DP

!  in thousands
real(DP), parameter :: iter_start=100; 
real(DP), parameter :: iter_step=10;
real(DP), parameter :: iter_end=170; 

integer, parameter :: ntrunk=3;
integer, parameter :: ngen=3; 
integer, parameter :: nproc=64;

character(200) :: fname, temp
integer :: iter, ng, np, nsamples
real(DP) :: Ap, fD, CD, Uinf, fD_proc, CD_proc, Uinf_proc
real(DP), dimension(:,:), allocatable :: dat
logical :: exst

!  Loop over each generation
do ng=1,ngen
  Ap = (scale_fact**2)**(ng-1)*d*l*(ntrunk**ng)
  CD=0.
  fD=0.
  Uinf=0.
  do np=0,nproc-1 ! Sum over all processors
    CD_proc=0.
    fD_proc=0.
    Uinf_proc=0.
    nsamples=0
    do iter=iter_start,iter_end,iter_step

      fname ='cylinder_skew_CD.dat.g'
      write (temp, '(i0)') ng
      fname = trim (fname) // temp
      fname = trim(fname) // '.c'
      write (temp, '(i0)') np
      fname = trim (fname) // temp
      fname = trim(fname) // '.ngen'
      write (temp, '(i0)') ngen
      fname = trim (fname) // temp
      fname = trim(fname) // '.'
      write (temp, '(i0)') iter
      fname = trim (fname) // temp
      fname = trim(fname) // 'k'
      
      inquire (file=fname, exist=exst)
      if(exst) then
!         write(*,*) fname
        call load_data()
        nsamples= nsamples + size(dat,1);
        CD_proc = CD_proc + sum(dat(:,2))
        fD_proc = fD_proc + sum(dat(:,3))
        deallocate(dat)
      endif        
    enddo
    if(nsamples > 0) then
      CD_proc = CD_proc/nsamples; !  Average over all proc samples
      fD_proc = fD_proc/nsamples;
      Uinf_proc = Uinf_proc/nsamples; !  Average over all proc samples
      CD = CD + CD_proc; !  Added proc component to gen total
      fD = fD + fD_proc;
    endif
  enddo
  write(*,*) 'im here'
  fname ='cylinder_skew_CD.dat.g'
  write (temp, '(i0)') ng
  fname = trim (fname) // temp
  fname = trim(fname) // '.ngen'
  write (temp, '(i0)') ngen
  fname = trim (fname) // temp
  open (unit = 11,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')
  write(11,*) 'Ap, CD, fD : ', Ap, CD, fD
  close(11)
  write(*,*) 'Gen #, Ap, CD, fD : ',ng, Ap, CD, fD
enddo

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
allocate(dat(nlines,4))

open(unit = 10,file=fname,action='read',status='old',form='formatted', &
  position='rewind')
do n=1,nlines
  read(10,fmt=*,iostat=ios) dat(n,1), dat(n,2), dat(n,3), dat(n,4)
enddo
close(10)

return
end subroutine load_data

end program cylinder_skew_post
