!**********************************************************************
module rns_ls
!**********************************************************************
use rns_base_ls
$if($CYLINDER_SKEW_LS)
use cylinder_skew_base_ls, only : tr_t, clindx_to_loc_id, brindx_to_loc_id, ntree
use cylinder_skew_ls, only : cylinder_skew_fill_tree_array_ls
!use cylinder_skew_ls, only : cylinder_skew_fill_cl_ref_plane_array_ls
!use cylinder_skew_ls, only : cylinder_skew_get_branch_id_ls
$endif

implicit none

save
private

public :: rns_init_ls !, rns_u_write_ls
public :: rns_CD_ls

character (*), parameter :: mod_name = 'rns_ls'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************
use messages
use param, only : USE_MPI, coord

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'

integer :: nt, np

call mesg ( sub_name, 'setting reference planes' )

$if($CYLINDER_SKEW_LS)
call cylinder_skew_fill_tree_array_ls()
$endif

!  Create cluster reference value plane
call rns_fill_cl_ref_plane_array_ls()

!  Allocate the branch force to the number of representative branches
if (.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
if(clforce_calc) allocate( clforce_t( size(cl_ref_plane_t, 1 ))) !  Assuming they are the same size
if(brforce_calc) allocate( brforce_t( tr_t(1)%nbranch ))
endif

!  Load the brindx file
call brindx_init()
  
return
end subroutine rns_init_ls

!**********************************************************************
subroutine rns_fill_cl_ref_plane_array_ls()
!**********************************************************************
use types, only : rprec
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
           
      do nb = 1, nbranch_p

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
      
      nullify(nbranch_p, clindx_p, origin_p)
      
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
use param, only : jt, USE_MPI, coord
implicit none

if(clforce_calc) then
  if(modulo (jt, clforce_nskip) == 0) then
    call rns_cl_CD_ls()
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then
      call rns_write_cl_CD_ls()
    endif
  endif
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
use types, only : rprec
use param, only : nx, ny, nz, dx, dy, dz
use param, only : USE_MPI, coord
$if($MPI)
use param, only : MPI_RPREC, MPI_SUM, rank_of_coord, comm, ierr
$endif
use sim_param, only : u
use functions, only : plane_avg_3D
use level_set_base, only : phi
use immersedbc, only : fx
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_cl_CD_ls'

integer, pointer, dimension(:) :: br_loc_id_p => null()
integer, pointer :: clindx_p => null()
integer :: i, j, k
integer :: nc, ncluster_tot
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
do nc = 1, ncluster_tot
  cl_ref_plane_t(nc) % u = plane_avg_3D( u(1:nx,:,1:nz), cl_ref_plane_t(nc) % p1, cl_ref_plane_t(nc) % p2, &
    cl_ref_plane_t(nc) % p3, cl_ref_plane_t(nc) % nzeta, cl_ref_plane_t(nc) % neta )
enddo

clforce_t % fD = 0._rprec

do k=1, nz - 1 ! since nz over laps
  do j = 1, ny
    do i = 1, nx

       if ( brindx(i,j,k) > 0 .and. phi(i,j,k) <= 0._rprec ) then 
     
         ! map brindx to clindx
         br_loc_id_p => brindx_to_loc_id(:,brindx(i,j,k))
         clindx_p => tr_t(br_loc_id_p(1)) % gen_t(br_loc_id_p(2)) % cl_t (br_loc_id_p(3)) % indx
           
         $if($MPI)
         cl_fD(clindx_p) = cl_fD(clindx_p) - fx(i,j,k) * dx * dy * dz
         $else
         clforce_t(clindx_p)%fD = clforce_t(clindx_p)%fD - fx(i,j,k) * dx * dy * dz
         $endif
         
         nullify(br_loc_id_p, clindx_p)
         
       endif
    enddo
  enddo
enddo
  
$if($MPI)
!  Need to sum forces over all processors
call mpi_reduce (cl_fD, clforce_t%fD, ncluster_tot , MPI_RPREC, MPI_SUM, rank_of_coord(0), comm, ierr)
deallocate(cl_fD)
$endif 

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) then

  do nc = 1, ncluster_tot 
   
    clforce_t(nc) % CD = clforce_t(nc)%fD / (0.5_rprec * cl_ref_plane_t(nc)%area * (cl_ref_plane_t(nc)%u)**2)
      
  enddo
    
endif

return
end subroutine rns_cl_CD_ls

!**********************************************************************
subroutine rns_write_cl_CD_ls()
!**********************************************************************
use io, only : write_real_data, write_tecplot_header_xyline
use param, only : jt_total, dt, path
use strmod
implicit none

character(*), parameter :: sub_name = mod_name // '.rns_write_cl_CD_ls'
character(*), parameter :: fname = path // 'output/rns_cl_CD_ls.dat'

logical :: exst
character(5000) :: var_list
integer :: nc, nvar, nvar_count
integer, pointer, dimension(:) :: cl_loc_id_p => null()

real(rprec), pointer, dimension(:) :: cl_CD_p

!  Write cluster force (CD) for all trees + time step
nvar = size( clforce_t, 1 ) + 1

if(write_tree_1_only) then

  nvar_count = 0
  
  nv_search : do nc = 1, nvar - 1
  
    cl_loc_id_p => clindx_to_loc_id(:,nc)
    
    if(cl_loc_id_p(1) == 1) then
    
      nvar_count = nvar_count + 1
      exit nv_search 
      
    endif
    
  enddo nv_search
  
  nvar = nvar_count + 1
  
endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  var_list = '"jt"'
  do nc = 1, nvar-1
  
    cl_loc_id_p => clindx_to_loc_id(:,nc)
    !  Create variable list name:
    call strcat(var_list, ',"CD<sub>')
    call strcat(var_list, cl_loc_id_p(1))
    call strcat(var_list, ',')
    call strcat(var_list, cl_loc_id_p(2))
    call strcat(var_list, ',')
    call strcat(var_list, cl_loc_id_p(3))
    call strcat(var_list, '</sub>"')
  enddo
  nullify(cl_loc_id_p)
  call write_tecplot_header_xyline(fname, 'rewind', trim(adjustl(var_list)))
endif

!write(*,*) '----------------'
!write(*,*) '(/ jt_total*dt, clforce_t%CD /) ', (/ jt_total*dt, clforce_t%CD /)
!write(*,*) '----------------'
call write_real_data(fname, 'append', nvar, (/ jt_total*dt, clforce_t(1:nvar)%CD /))



return
end subroutine rns_write_cl_CD_ls

!**********************************************************************
subroutine brindx_init ()
!**********************************************************************
use param, only : iBOGUS, coord
use messages
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


