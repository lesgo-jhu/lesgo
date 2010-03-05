!**********************************************************************
module rns_ls
!**********************************************************************
use types, only : rprec
use param
use rns_base_ls
use messages
implicit none

save
private

public :: rns_init_ls, rns_u_write_ls

character (*), parameter :: mod_name = 'rns_ls'

contains

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'
character(64), parameter :: fbase= path // 'rns_planes_ls.out'

character(64) :: fname, temp

integer :: nt, np

call mesg ( sub_name, 'reading rns_planes' )

!  Set the number of trees to use
rns_t % ntrees = 1
!  Set plane average write
rns_t % plane_u_calc = .true.

do nt=1, rns_t % ntrees

  write(temp, '(".t",i0)') nt
  fname = trim (fbase) // temp
  
  !open (unit = 2,file = fname, status='old',form='formatted', &
    !action='read',position='rewind')
  open (unit = 2,file = fname, status='old',form='unformatted', &
    action='read',position='rewind')
 
  !  Get the number of planes 
  read(2) rns_t%nplanes
  
  if(.not. allocated(rns_planes_t)) allocate(rns_planes_t(rns_t % nplanes))

  do np=1,rns_t % nplanes
	read(2) rns_planes_t(np)%indx, rns_planes_t(np)%bp
  enddo

  close(2)
enddo
  
return
end subroutine rns_init_ls

!**********************************************************************
subroutine rns_u_write_ls()
!**********************************************************************
use sim_param, only : u
use functions, only : plane_avg_3D
use param, only : jt_total, dt_dim
implicit none

character(*), parameter :: fname= path // 'output/rns_planes_u.out'

character(64) :: fmt

integer :: np


!fmt2 = "'(" // fmt1 // ")'"
!fmt1 = trim(adjustl(fmt1))

!write(*,*) 'fmt1 ', fmt1

!stop

do np = 1, rns_t%nplanes
  rns_planes_t(np)%u = plane_avg_3D(u,rns_planes_t(np)%bp,20,20)
  write(*,*) 'rns_planes_t(np)%bp : ', rns_planes_t(np)%bp
  write(*,*) 'rns_planes_t(np)%u : ', rns_planes_t(np)%u
enddo

if(coord == 0) then
open(unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position='append')
write (fmt, '("(",i0,"f6.3)")') rns_t%nplanes
write(2,fmt) jt_total*dt_dim, rns_planes_t(:)%u

close(2)
endif

return
end subroutine rns_u_write_ls

end module rns_ls


