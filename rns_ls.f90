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
!use cylinder_skew_ls, only : cylinder_skew_fill_cl_ref_plane_array_ls
!use cylinder_skew_ls, only : cylinder_skew_get_branch_id_ls
$endif

implicit none

save
private

!----- Parameters -----
logical, parameter :: brindx_calc = .false.
logical, parameter :: clindx_calc = .true.
!----------------------

logical :: brindx_initialized = .false.
integer :: brindx(ld, ny, nz)

public :: rns_init_ls !, rns_u_write_ls
public :: rns_CD_ls

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
$endif

!  Create cluster reference value plane
call rns_fill_cl_ref_plane_array_ls()

!  Allocate the branch force to the number of representative branches
if (.not. USE_MPI .or. (USE_MPI .and. coord == 0) then
if(clforce_calc) allocate( clforce_t( tr_t(1)%ncluster ))
if(clforce_calc) allocate( brforce_t( tr_t(1)%nbranch ))
endif

!  Load the brindx file
call brindx_init()
  
return
end subroutine rns_init_ls

!**********************************************************************
subroutine rns_fill_cl_ref_plane_array_ls()
!**********************************************************************

use param, only : dy, dz

implicit none

real(rprec), parameter :: alpha=1._rprec

integer :: nt, ng, nc, nb

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)

integer, pointer :: clindx_p, nbranch_p
real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer, dimension(:) :: origin_p

nullify(d_p, l_p, skew_angle_p, clindx_p)

!  allocate ref_array_t
allocate(cl_ref_plane_t( size(clindx_to_loc_id,2) ))

do nt=1, ntree

  do ng=1, tr_t(nt)%ngen
  
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
    
      nbranch_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
      
      clindx_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%indx
           
      do nb = 1, nbranch

        d_p          => tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t(nb)%d
        l_p          => tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t(nb)%l
        skew_angle_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t(nb)%skew_angle
        
        h         = l_p * dcos(skew_angle_p)
        h_m       = h_m + h
        area_proj = area_proj + d_p * h
        
        nullify(d_p, l_p, skew_angle_p)
      
      enddo
      
      !  Mean height of branch cluster  and height of reference area
      h_m = h_m / nbranch_p
      !  width of reference area
      w   = area_proj / h_m     

      cl_ref_plane_t(clindx_p) % area = area_proj
      !  These are defined to be x - planes (no not the NASA experimental planes)
      cl_ref_plane_t(clindx_p) % nzeta = ceiling( w / dy + 1)
      cl_ref_plane_t(clindx_p) % neta  = ceiling( h_m / dz + 1)
      
      origin_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%origin
      
      !  Offset in the upstream x-direction
      zeta_c = origin_p + (/ -alpha * w, 0._rprec, 0._rprec /)
      
      cl_ref_plane_t(clindx_p) % p1    = zeta_c 
      cl_ref_plane_t(clindx_p) % p1(2) = cl_ref_plane_t(clindx_p) % p1(2) + w / 2._rprec
      
      cl_ref_plane_t(clindx_p) % p2    = cl_ref_plane_t(clindx_p) % p1
      cl_ref_plane_t(clindx_p) % p2(2) = cl_ref_plane_t(clindx_p) % p2(2) - w
      
      cl_ref_plane_t(clindx_p) % p3    = cl_ref_plane_t(clindx_p) % p2
      cl_ref_plane_t(clindx_p) % p3(3) = cl_ref_plane_t(clindx_p) % p3(3) + h_m
      
      nullify(clindx_p, origin_p)
      
    enddo
    
  enddo
 
enddo
      
return

end subroutine rns_fill_cl_ref_plane_array_ls

!**********************************************************************
subroutine rns_CD_ls()
!**********************************************************************
!  This subroutine handles all CD calculation within the RNS module; 
!  all CD and force calculations associated with
!
!  tree -> generation -> cluster -> branch
!
!  are handled here
!

implicit none

if(clforce_calc) then
  call rns_cl_CD_ls()
  call rns_write_cl_CD_ls()
endif

return
end subroutine rns_CD_ls


!**********************************************************************
subroutine rns_cl_CD_ls()
!**********************************************************************
!  This subroutine computes the CD of the branch cluster (cl) associated
!  with each region dictated by the brindx value. The cl is mapped from 
!  brindex
!
use param, only : dx, dy, dz
use sim_param, only : u
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_cl_CD_ls'

integer, pointer :: brindx_p => null(), clindx_p => null()
integer :: ncluster_tot
real(rprec), pointer :: fx_p => null()
$if ($MPI)
real(rprec), pointer, dimension(:) :: cl_fD
$endif
!ncluster_p => tr_t(1)%ncluster

!  Get the total number of clusters; needs to match the
!  size of cl_ref_plane_t
ncluster_tot = size( cl_ref_plane_t )

!ncluster_tot = 0
!do n=1, ntree
!  ncluster_tot = ncluster_tot + tr_t(n) % ncluster
!enddo

$if ($MPI)
allocate (cl_fD( ncluster_tot ))

cl_fD = 0._rprec

$endif

!  Get reference velocity for all reference planes
do n = 1, ncluster_tot
  cl_ref_plane_t(n) % u = plane_avg_3D( u(1:nx,:,1:nz), cl_ref_plane_t(n) % p1, cl_ref_plane_t(n) % p2, &
    cl_ref_plane_t(n) % p3, cl_ref_plane_t(n) % nzeta, cl_ref_plane_t(n) % neta )
enddo

clforce_t % fx = 0._rprec

do k=1, nz - 1 ! since nz over laps
  do j = 1, ny
    do i = 1, nx

       if ( brindx(i,j,k) > 0 .and. phi(i,j,k) <= 0._rprec ) then 
       
         fx_p => fx(i,j,k)
     
         ! map brindx to clindx
         brindx_p => brindx_to_loc_id(:,brindx(i,j,k))
         clindx_p => tr_t(branch_id(branch_id(1)) % gen_t(branch_id(2)) % cl_t (branch_id(3)) % indx
           
         $if($MPI)
         cl_fD(clindx_p) = cl_fD(clindx_p) - fx_p * dx * dy * dz
         $else
         clforce_t(clindx_p)%fD = clforce_t(clindx_p)%fD - fx_p * dx * dy * dz
         $endif
         
         nullify(brindx_p, clindx_p, fx_p)
         
       endif
    enddo
  enddo
enddo
  
$if($MPI)
!  Need to sum forces 
call mpi_reduce (cl_fD, clforce_t%fD, ncluster_tot , MPI_RPREC, MPI_SUM, rank_of_coord(0), comm, ierr)
deallocate(cl_fD)
$endif 

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then

  do n = 1, ncluster_tot 
   
    clforce_t(n)%CD = clforce_t(n)%fD / (0.5_rprec * cl_ref_plane_t(n)%area * (cl_ref_plane_t(n)%u)**2)
      
  enddo
    
endif

return
end subroutine rns_cl_CD_ls

!**********************************************************************
subroutine rns_write_cl_CD_ls()
!**********************************************************************
use io, only : write_real_data
use param, only : jt_total, dt

implicit none

integer, intent(in) :: jt

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_CD_ls'
character(*), parameter :: fname = 'rns_cl_CD_ls.dat'

nvars = size( cl_ref_planes_t, 1 ) + 1
call write_real_data(fname, write_pos, nvars, (/ jt_total*dt, clforce_t%CD /))


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


