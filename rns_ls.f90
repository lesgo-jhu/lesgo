!**********************************************************************
module rns_ls
!**********************************************************************
use types, only : rprec
use param
use rns_base_ls
use messages
use strmod
implicit none

save
private

public :: rns_init_ls, rns_u_write_ls

character (*), parameter :: mod_name = 'rns_ls'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'
character(64), parameter :: fbase= path // 'rns_planes_ls.out'

character(64) :: fname

integer :: nt, np

call mesg ( sub_name, 'reading rns_planes' )

!  Set the number of trees to use
rns_t % ntrees = 1
!  Set plane average write
rns_t % plane_u_calc = .true.

do nt=1, rns_t % ntrees
  
  fname=fbase
  call strcat(fname,'.t')
  call strcat(fname,nt)
  
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
use sim_param, only : u, v, w
use functions, only : plane_avg_3D
use param, only : jt_total, dt_dim
use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag
use io, only : write_tecplot_header_xyline, write_real_data_append, interp_to_uv_grid
implicit none

character(*), parameter :: fbase= path // 'output/uvw_rns_planes.dat'

character(64) :: fmt
character(120) :: fname

integer :: np
logical :: exst

fmt=''

call interp_to_uv_grid(w, w_uv, w_uv_tag)

do np = 1, rns_t%nplanes
  
  fname = fbase
  call strcat(fname,'.p')
  call strcat(fname,np)
  
  rns_planes_t(np)%u = plane_avg_3D(u,rns_planes_t(np)%bp,20,20)
  rns_planes_t(np)%v = plane_avg_3D(v,rns_planes_t(np)%bp,20,20)
  rns_planes_t(np)%w = plane_avg_3D(w,rns_planes_t(np)%bp,20,20)

  if(coord == 0) then
  
    inquire (file=fname, exist=exst)
	if(.not. exst) call write_tecplot_header_xyline(fname, 'rewind', &
	    '"t (s)", "u", "v", "w"')
    
	call write_real_data_append(fname, (/ jt_total*dt_dim, &
	  rns_planes_t(np)%u, rns_planes_t(np)%v, rns_planes_t(np)%w /), &
	  4)

  endif
  
enddo

return
end subroutine rns_u_write_ls

end module rns_ls


