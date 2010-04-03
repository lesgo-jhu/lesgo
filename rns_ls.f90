!**********************************************************************
module rns_ls
!**********************************************************************
use types, only : rprec
use param
use rns_base_ls
use messages
use strmod

$if($CYLINDER_SKEW_LS)
use cylinder_skew_base_ls
use cylinder_skew_ls, only : cylinder_skew_fill_tree_array_ls
use cylinder_skew_ls, only : cylinder_skew_fill_ref_plane_array_ls
use cylinder_skew_ls, only : cylinder_skew_get_branch_id_ls
$endif

implicit none

save
private

$if($CYLINDER_SKEW_LS)
logical, parameter :: brindx_based = .false.
logical, parameter :: clindx_based = .true.
$endif

logical :: brindx_initialized = .false.
integer :: brindx(ld, ny, nz)

public :: rns_init_ls !, rns_u_write_ls

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

call mesg ( sub_name, 'setting reference planes' )

$if($CYLINDER_SKEW_LS)

call cylinder_skew_fill_tree_array_ls()
!  Create reference value planes
call cylinder_skew_fill_ref_plane_array_ls()

!  Allocate the branch force to the number of representative branches
if (.not. USE_MPI .or. (USE_MPI .and. coord == 0) then
if(clforce_calc) allocate( clforce_t( tr_t(1)%ncluster ))
if(clforce_calc) allocate( brforce_t( tr_t(1)%nbranch ))
endif

$endif

!  Load the brindx file
call brindx_init()
  
return
end subroutine rns_init_ls

!!**********************************************************************
!subroutine rns_u_write_ls()
!!**********************************************************************
!use sim_param, only : u, v, w
!use functions, only : plane_avg_3D
!use param, only : jt_total, dt_dim, Ny, Nz_tot, L_y, L_z, nproc
!use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag
!use io, only : write_tecplot_header_xyline, write_real_data, interp_to_uv_grid
!implicit none

!character(*), parameter :: fbase= path // 'output/uvw_rns_planes.dat'

!character(64) :: fmt
!character(120) :: fname

!integer :: np, nzeta, neta
!logical :: exst

!fmt=''

!call interp_to_uv_grid(w, w_uv, w_uv_tag)

!do np = 1, rns_t%nplanes
!  
!  fname = fbase
!  call strcat(fname,'.p')
!  call strcat(fname,np)
! 
!!  nzeta = abs(rns_planes_t(np)%bp(1,1)-rns_planes_t(np)%bp(1,2))/L_y*Ny
!!  neta  = abs(rns_planes_t(np)%bp(3,3)-rns_planes_t(np)%bp(3,2))/L_z*nproc*Nz_tot 
!  nzeta=100
!  neta=100

!  rns_planes_t(np)%u = plane_avg_3D(u,rns_planes_t(np)%bp,nzeta,neta)
!  rns_planes_t(np)%v = plane_avg_3D(v,rns_planes_t(np)%bp,nzeta,neta)
!  rns_planes_t(np)%w = plane_avg_3D(w_uv,rns_planes_t(np)%bp,nzeta,neta)

!  if(coord == 0) then
!  
!    inquire (file=fname, exist=exst)
!	if(.not. exst) call write_tecplot_header_xyline(fname, 'rewind', &
!	    '"t (s)", "u", "v", "w"')
!	!if (.not. exst) call write_tecplot_header_ND(fname, 'rewind', 4, (/ Nx /), '"t (s)", "u", "v", "w"', &
!	!	  np, 2)

!	call write_real_data(fname, 'append', 4, (/ jt_total*dt_dim, &
!	  rns_planes_t(np)%u, rns_planes_t(np)%v, rns_planes_t(np)%w /))

!  endif
!  
!enddo

!return
!end subroutine rns_u_write_ls

!**********************************************************************
subroutine rns_CD_ls()
!**********************************************************************
!  This subroutine computes the CD associated with each region dictated
!  by the brindx value; this could be a branch, cluster of branches, etc.
!  Depends only on brindx.
!
use param, only : dx, dy, dz
use sim_param, only : u
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_CD_ls'

integer, pointer :: brindx_p => null(), clindx_p => null()
integer, pointer :: ncluster => null(), nbranch => null()
real(rprec), pointer :: fx_p => null()

ncluster => tr_t(1)%ncluster
nbranch  => tr_t(1)%nbranch

$if ($MPI)
allocate (cl_fD(ncluster))
allocate (br_fD(nbranch))

cl_fD = 0._rprec
bf_fD = 0._rprec

$endif

!  Get reference velocity for all reference planes
do n = 1, size(ref_plane_t)
  ref_plane_t(n) % u = plane_avg_3D( u, ref_plane_t(n) % p1, ref_plane_t(n) % p2, &
    ref_plane_t(n) % p3, ref_plane_t(n) % nzeta, ref_plane_t(n) % neta )
enddo

if(clforce_calc) clforce_t % fx = 0._rprec
if(brforce_calc) brforce_t % fx = 0._rprec


do k=1, nz - 1 ! since nz over laps
  do j = 1, ny
    do i = 1, nx

       if ( brindx(i,j,k) > 0 ) then
       
         fx_p => fx(i,j,k)
       
         $if ($CYLINDER_SKEW_LS)
         !  Map branch id if using clindx)
         if(clforce_calc) then
         
           ! map brindx to clindx
           brindx_p => brindx_to_loc_id(:,brindx(i,j,k))
           clindx_p => tr_t(branch_id(branch_id(1)) % gen_t(branch_id(2)) % cl_t (branch_id(3)) % indx
           
           $if($MPI)
           cl_fD(clindx_p) = cl_fD(clindx_p) - fx_p * dx * dy * dz
           $else
           clforce_t(clindx_p)%fD = clforce_t(clindx_p)%fD - fx_p * dx * dy * dz
           $endif
           
           nullify(brindx_p)

         elseif(brforce_calc) then
         
           brindx_p => brindx(i,j,k)
           
           $if($MPI)
           br_fD(brindx_p) = br_fD(brindx_p) - fx_p * dx * dy * dz
           $else
           brforce_t(brindx_p)%fD = brforce_t(brindx_p)%fD - fx_p * dx * dy * dz
           $endif
           
         endif
         
         $else
         
         if(brforce_calc) then
           brindx_p => brindx(i,j,k)
           
           $if($MPI)
           br_fD(brindx_p) = br_fD(brindx_p) - fx_p * dx * dy * dz
           $else
           brforce_t(brindx_p)%fD = brforce_t(brindx_p)%fD - fx_p * dx * dy * dz
           $endif
           
         endif
         
         $endif
         
         nullify(brindx_p, clindx_p, fx_p)
         
       endif
    enddo
  enddo
enddo



if(clforce_calc) then
  
  $if($MPI)
  !  Need to sync forces
  call mpi_reduce (cl_fD, clforce_t%fD, ncluster, MPI_RPREC, MPI_SUM, rank_of_coord(0), comm, ierr)
  deallocate(cl_fD)
  $endif 

  $if($MPI)
  if(coord == 0) then
  $endif 
  
    do n = 1, ncluster
    
      clforce_t(n)%CD = clforce_t(n)%fD / (0.5_rprec * ref_plane_t(n)%area * (ref_plane_t(n)%u)**2)
      
    enddo
    
  $if($MPI)
  endif
  $endif    
    
endif

if(brforce_calc) then

  $if($MPI)
  call mpi_reduce (br_fD, brforce_t%fD, nbranch, MPI_RPREC, MPI_SUM, rank_of_coord(0), comm, ierr)
  deallocate(br_fD)
  $endif
  
  $if($MPI)
  if(coord == 0) then
  $endif
    
    do n = 1, nbranch
    
      ! map brindx to clindx
      brindx_p => brindx_to_loc_id(:,n)
      clindx_p => tr_t(branch_id(branch_id(1)) % gen_t(branch_id(2)) % cl_t (branch_id(3)) % indx      
    
      brforce_t(n)%CD = brforce_t(n)%fD / (0.5_rprec * ref_plane_t(clindx_p)%area * (ref_plane_t(clindx_p)%u)**2)
      
      nullify(brindx_p, clindx_p)
      
    enddo
  
  $if($MPI)
  endif
  $endif 
  
endif

return
end subroutine rns_CD_ls

!**********************************************************************
subroutine rns_write_CD_ls(jt)
!**********************************************************************
use io, only : write_tecplot_header_xyline, write_real_data

implicit none

integer, intent(in) :: jt

character(*), parameter :: sub_name = mod_name // '.rns_write_CD_ls'
character(*), parameter :: fname = 'rns_CD_ls.dat'


call write_tecplot_header_xyline(point_t%fname(i), 'rewind', var_list)





return
end subroutine rns_write_CD_ls

!**********************************************************************
subroutine brindx_init ()
!**********************************************************************
use param, only : iBOGUS
implicit none

character (*), parameter :: sub_name = mod_name // '.brindx_init'
character (*), parameter :: fbrindx_in = 'brindx.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fbrindx_in_MPI
$endif

integer :: ip

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fbrindx_in_MPI, '(a,a,i0)') fbrindx_in, MPI_suffix, coord
    
  inquire (file=fbrindx_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fbrindx_in_MPI)

  open (1, file=fbrindx_in_MPI, action='read', position='rewind',  &
         form='unformatted')
  read (1) brindx(:, :, 1:nz-1)
  close (1)

  brindx(:, :, nz) = iBOGUS

$else

  inquire (file=fbrindx_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fbrindx_in)

  open (1, file=fbrindx_in, action='read', position='rewind',  &
         form='unformatted')
  read (1) brindx
  close (1)

$endif

brindx_initialized = .true.

end subroutine brindx_init

end module rns_ls


