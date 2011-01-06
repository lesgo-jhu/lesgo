!**********************************************************************
module rns_cyl_skew_ls
!**********************************************************************
!  The purpose of this module is to provide all geometry specific info
!  to the rns module. This module sets and defines all variables used
!  in the rns module and fills all of the element structures used for
!  the rns calculations.
!
!  If additional geometric configurations (foo) are desired to be used
!  with the rns module, then this module can be used as a template
!  with all subsequent names replaced with foo (e.g. rns_foo_ls). Only one
!  rns definitions module may be used at a time.
!
use types, only : rprec
use param, only : ld, ny, nz
use rns_base_ls
use cyl_skew_base_ls

implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

save
private

public rns_init_ls, rns_tree_layout

character (*), parameter :: mod_name = 'rns_cyl_skew_ls'

!  Options: 'default', 'kc-3'
character(*), parameter :: rns_tree_layout = 'kc-3'

integer, pointer, dimension(:) :: cl_to_r_elem_map
integer, pointer, dimension(:) :: cl_to_beta_elem_map
integer, pointer, dimension(:) :: cl_to_b_elem_map

integer, pointer, dimension(:) :: r_elem_to_basecl_map
integer, pointer, dimension(:) :: beta_elem_to_basecl_map
integer, pointer, dimension(:) :: b_elem_to_basecl_map

!integer, pointer, dimension(:) :: rns_tree_map(:) ! This maps the tree number from cyl_skew to the trees considered during rns
integer, pointer, dimension(:) :: rns_to_cyl_skew_tree_map
integer, pointer, dimension(:) :: cyl_skew_to_rns_tree_map

integer :: clindx(ld, ny, $lbz:nz)
logical :: clindx_initialized = .false.

contains

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************
use messages
use param, only : USE_MPI, coord
use cyl_skew_ls, only : fill_tree_array_ls
use rns_ls, only : rns_force_init_ls

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Initializing RNS Data Structure'
  write(*,*) ' '
endif

if ( rns_ntree > ntree ) call error ( sub_name, 'rns_ntree > ntree ')

!----- Load Data -----
! Load clindx
call clindx_init()
! Load filtered indicator function (chi)
call chi_init()
!!----- Load Data -----

!----- Fill CYL_SKEW Data Structures -----
!  Fill the cyl_skew tree array
call fill_tree_array_ls()
!  Fill the mapping for the rns trees
call fill_rns_tree_map()
!!  Get the number of resolved clusters 
!call set_ncluster_reslv()
!!  Get the total number of clusters
!call set_ncluster_tot()
!----- Fill CYL_SKEW Data Structures -----

!----- Fill RNS_CYL_SKEW Data Structures -----
!  Set the number of r elements
call set_nr_elem()
!  Set the number of beta elements
call set_nbeta_elem()
!  Set the number of b elements
call set_nb_elem()

!  Fill the mapping from tree clusters to r regions
call fill_cl_to_r_elem_map()  

!  Fill the mapping from tree clusters to beta regions
call fill_cl_to_beta_elem_map()

!  Fill the mapping from tree clusters to b regions
call fill_cl_to_b_elem_map()
!  Fill the r_elem struct
call fill_r_elem()
!  Fill the beta_elem struct
call fill_beta_elem()
!  Fill the b_elem struct
call fill_b_elem()
!----- Fill RNS_CYL_SKEW Data Structures -----

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) 'ncluster_reslv : ', ncluster_reslv
endif

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'RNS Data Structure Initialized'
  write(*,*) ' '
endif 

!  Now update old force data
call rns_force_init_ls()

return

end subroutine rns_init_ls

!**********************************************************************
subroutine clindx_init ()
!**********************************************************************
use param, only : iBOGUS, coord, nz, USE_MPI
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.clindx_init'
character (*), parameter :: fname_in = 'clindx.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fname_in_MPI
$endif

logical :: opn, exst

!---------------------------------------------------------------------
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> loading clindx'
  write(*,*) 'size(clindx,1) : ', size(clindx,1)
endif



inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fname_in_MPI, '(a,a,i0)') fname_in, MPI_suffix, coord
    
  inquire (file=fname_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fname_in_MPI)

  $if ($READ_BIG_ENDIAN)
  open (1, file=fname_in_MPI, action='read', position='rewind', &
    form='unformatted', convert='big_endian')
  $elseif ($READ_LITTLE_ENDIAN)
  open (1, file=fname_in_MPI, action='read', position='rewind', &
    form='unformatted', convert='little_endian')
  $else
  open (1, file=fname_in_MPI, action='read', position='rewind', form='unformatted')
  $endif                              

  read (1) clindx
  close (1)

  clindx(:, :, nz) = iBOGUS

$else

  inquire (file=fname_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fname_in)

  $if ($READ_BIG_ENDIAN)
  open (1, file=fname_in, action='read', position='rewind', &
    form='unformatted', convert='big_endian')
  $elseif ($READ_LITTLE_ENDIAN)
  open (1, file=fname_in, action='read', position='rewind', &
    form='unformatted', convert='little_endian')
  $else
  open (1, file=fname_in, action='read', position='rewind', form='unformatted')
  $endif   

  read (1) clindx
  close (1)

$endif

clindx_initialized = .true.

end subroutine clindx_init

!**********************************************************************
subroutine chi_init ()
!**********************************************************************
use param, only : iBOGUS, coord, nz
use messages
use types, only : rprec
implicit none

character (*), parameter :: sub_name = mod_name // '.chi_init'
character (*), parameter :: fchi_in = 'chi.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fchi_in_MPI
$endif

real(rprec) :: chi_max

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fchi_in_MPI, '(a,a,i0)') fchi_in, MPI_suffix, coord
    
  inquire (file=fchi_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fchi_in_MPI)

  $if ($READ_BIG_ENDIAN)
  open (1, file=fchi_in_MPI, action='read', position='rewind', &
    form='unformatted', convert='big_endian')
  $elseif ($READ_LITTLE_ENDIAN)
  open (1, file=fchi_in_MPI, action='read', position='rewind', &
    form='unformatted', convert='little_endian')
  $else
  open (1, file=fchi_in_MPI, action='read', position='rewind', form='unformatted')
  $endif                               

  read (1) chi
  close (1)

  chi(:, :, nz) = iBOGUS

$else

  inquire (file=fchi_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fchi_in)
  
  $if ($READ_BIG_ENDIAN)
  open (1, file=fchi_in, action='read', position='rewind', &
    form='unformatted', convert='big_endian')
  $elseif ($READ_LITTLE_ENDIAN)
  open (1, file=fchi_in, action='read', position='rewind', &
    form='unformatted', convert='little_endian')
  $else
  open (1, file=fchi_in, action='read', position='rewind', form='unformatted')
  $endif   
  read (1) chi
  close (1)

$endif

chi_max = maxval(chi)
write(*,*) 'coord, chi_max : ', coord, chi_max

chi_initialized = .true.

end subroutine chi_init

!**********************************************************************
subroutine fill_rns_tree_map ()
!**********************************************************************
use messages
use param, only : USE_MPI, coord
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_rns_tree_map'
integer :: n, nt

!  this is used to map the brindx to correct rns tree
allocate( rns_to_cyl_skew_tree_map( rns_ntree ) )
allocate( cyl_skew_to_rns_tree_map( ntree ) )

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling RNS Tree Mapping'
  write(*,*) ' '
endif

!  Assign rns_tree_map layout
select case (rns_tree_layout)

  case('default')

    do nt = 1, ntree
  
      if( nt <= rns_ntree ) then
          cyl_skew_to_rns_tree_map( nt ) = nt
          rns_to_cyl_skew_tree_map( nt ) = nt
      else
          cyl_skew_to_rns_tree_map(nt) = rns_ntree
      endif
    
    enddo

  case('kc-3')

    cyl_skew_to_rns_tree_map = (/ 1, 2, 3, 3, 4, 4, 4, 4 /)
    rns_to_cyl_skew_tree_map = (/ 1, 2, 3, 5 /)

  case default

    call error(sub_name, 'rns_tree_layout not specified correctly')

end select

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
write(*,*) '----> RNS Tree Mapping : '
  write(*,*) '| ID       | RNS NT   |'
  do n=1, ntree
    write(*,'(i12,i11)') n, cyl_skew_to_rns_tree_map(n)
  enddo
endif

return

end subroutine fill_rns_tree_map

!**********************************************************************
subroutine set_nr_elem()
!**********************************************************************
use param, only : USE_MPI, coord
implicit none

integer :: nt, ng, nc

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(tr_t_p, gen_t_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Setting RNS NR_ELEM'
  write(*,*) ' '
endif

nr_elem = 0
do nt = 1, rns_ntree
  
  !  Point to RNS mapped tree
  tr_t_p => tr_t( rns_to_cyl_skew_tree_map(nt) ) 
  
  do ng=1, tr_t_p % ngen_reslv

    !  Point to the generation of the rns mapped tree  
    gen_t_p => tr_t_p % gen_t( ng )
	
    do nc = 1, gen_t_p % ncluster 
	
      nr_elem = nr_elem + 1
    
    enddo 
    
	  nullify(gen_t_p)    
    
  enddo
  
  nullify( tr_t_p )
  
enddo

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
write(*,*) '----> NR_ELEM : ', nr_elem
endif

return
end subroutine set_nr_elem

!**********************************************************************
subroutine set_nb_elem()
!**********************************************************************
use param, only : USE_MPI, coord
implicit none

integer :: nt, nc

type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(gen_t_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Setting RNS NB_ELEM'
  write(*,*) ' '
endif

nb_elem = 0
do nt = 1, rns_ntree

  !  Point to the resolved generation of the rns mapped tree  
  gen_t_p => tr_t( rns_to_cyl_skew_tree_map( nt ) ) % gen_t( ngen_reslv )
  
  do nc = 1, gen_t_p % ncluster 
  
    nb_elem = nb_elem + 1
    
  enddo    

  nullify(gen_t_p)
  
enddo

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
write(*,*) '----> NB_ELEM : ', nb_elem
endif

return
end subroutine set_nb_elem

!**********************************************************************
subroutine set_nbeta_elem()
!**********************************************************************
use param, only : USE_MPI, coord
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.set_nbeta_elem'

integer :: nt, nc

type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(gen_t_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Setting RNS NBETA_ELEM'
  write(*,*) ' '
endif

!  First check if ngen <= ngen_reslv
if( ngen <= ngen_reslv ) call error(sub_name, 'ngen <= ngen_reslv : RNS cannot be implemented')

nbeta_elem = 0
do nt = 1, rns_ntree

  !  Point to the resolved generation of the rns mapped tree  
  gen_t_p => tr_t( rns_to_cyl_skew_tree_map( nt ) ) % gen_t( ngen_reslv + 1)
	
  do nc = 1, gen_t_p % ncluster 
  
    nbeta_elem = nbeta_elem + 1
    
	enddo    
	
	nullify(gen_t_p)
  
enddo

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
write(*,*) '----> NBETA_ELEM : ', nbeta_elem
endif

return
end subroutine set_nbeta_elem

!!**********************************************************************
!subroutine fill_elem_ref_region(elem, nelem, ref_region_t)
!!**********************************************************************
!use param, only : dy, dz, USE_MPI, coord
!use param, only : nx, ny, nz
!use messages
!use cyl_skew_base_ls, only : tr_t, tree, generation
!implicit none

!character (*), parameter :: sub_name = mod_name // '.fill_elem_ref_region'


!character(*), intent(IN) :: elem
!integer, intent(IN) :: nelem
!type(ref_region), intent(IN), target, dimension(nelem) :: ref_region_t

!integer :: n, nb
!integer :: nxi, neta, nzeta

!real(rprec) :: h, w, xi_c(3)
!real(rprec), dimension(3) :: p1, p2, p3

!integer, pointer, dimension(:) :: cl_loc_id_p

!real(rprec), pointer :: d_p, l_p, skew_angle_p
!real(rprec), pointer :: hbot_p, htop_p

!type(cluster),    pointer :: cl_t_p
!type(branch),     pointer :: br_t_p

!type(ref_region), pointer :: ref_region_t_p

!nullify(cl_loc_id_p)
!nullify(d_p, l_p, skew_angle_p)
!nullify(hbot_p, htop_p)
!nullify(cl_t_p, br_t_p)
!nullify(ref_region_t_p)

!if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
!  write(*,*) ' '
!  write(*,*) '--> Filling R_ELEM Reference Plane Arrays'
!  write(*,*) ' '
!  write(*,*) '----> Reference Region Information : '
!  write(*,'(a116)') '| ID       | NPOINT   | AREA     | P1.X   | P1.Y   | P1.Z   | P2.X   | P2.Y   | P2.Z   | P3.X   | P3.Y   | P3.Z   |'
!endif

!do n=1, nelem

!  select case( elem )
!    case ('r')
!	
!      !  Get the base cluster local id
!      cl_loc_id_p => clindx_to_loc_id(:, r_elem_to_basecl_map(n) )
!  
!      cl_t_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % cl_t( cl_loc_id_p(3) )

!      !  Point to the top and bottom of the plane
!      hbot_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % bplane
!      htop_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % tplane	

!    case ('b')
!	
!	  !  Point to ref_region_t of the resolved element
!      ref_region_t_p => b_elem_t ( n ) % ref_region_t

!      !  Get the base cluster local id
!      cl_loc_id_p => clindx_to_loc_id(:, b_elem_to_basecl_map(n) )
!  
!      cl_t_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % cl_t( cl_loc_id_p(3) )

!      !  Point to the top and bottom of the plane
!      hbot_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % bplane
!      htop_p => tr_t( cl_loc_id_p(1) ) % gen_t( ngen ) % tplane

!    case ('beta')
!	
!      !  Get the base cluster local id
!      cl_loc_id_p => clindx_to_loc_id(:, beta_elem_to_basecl_map(n) )
!  
!      cl_t_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % cl_t( cl_loc_id_p(3) )

!      !  Point to the top and bottom of the plane
!      hbot_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % bplane
!      htop_p => tr_t( cl_loc_id_p(1) ) % gen_t( ngen ) % tplane	

!    case default
!        call error( sub_name, 'Invalid element specification.' )
!  end select  
!  


!  h = htop_p - hbot_p
!  w = alpha_width * h
! 
!  nxi   = ceiling ( 2. * alpha_dist * h / dx ) + 1
!  neta  = ceiling( w / dy ) + 1
!  nzeta = ceiling( h / dz ) + 1  
!      
!  !  Offset in the upstream x-direction
!  xi_c = cl_t_p % origin + (/ -alpha_dist * h, 0._rprec, 0._rprec /)
!      
!  !  Set the ordered corner points of the plane
!  p1    = xi_c 
!  p1(2) = p1(2) + w / 2._rprec
!    
!  p2    = p1
!  p2(2) = p2(2) - w
!      
!  p3    = p2
!  p3(3) = p3(3) + h
!  
!  p4    = p1
!  p4(1) = p4(1) + 2. * alpha_dist * h
!      
!  !  Point to ref_region_t of the resolved element
!  ref_region_t_p => 

!  ref_region_t(n) % area    = h * w
!  ref_region_t(n) % npoint = nxi * neta * nzeta
!    
!  !  Check if the element has been allocated
!  !if( associated( ref_region_t_p % points ) ) then
!  !  call error(sub_name, 'reference region points already allocated.')
!  !else
!  allocate( ref_region_t(n) % points( 3, ref_region_t(n) % npoint ) )
!  !endif
!            
!  !call set_points_in_plane( p1, p2, p3, nxi, neta, ref_region_t_p % points )  
!  call set_points_in_box( p1, p2, p3, p4, nzeta, neta, ref_region_t(n) % points ) 

!  !  Finally initialize velocity reference components
!  ref_region_t_p % u = 0._rprec
!  
!  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
!    write(*,'(i12,i11,f11.6, 9f9.4)') n, r_elem_t(n) % ref_region_t % npoint, r_elem_t(n) % ref_region_t % area, p1(:), p2(:), p3(:)
!  endif  
!  
!      
!  nullify( ref_region_t_p )
!  nullify( htop_p, hbot_p )
!  nullify( cl_t_p )
!  nullify( cl_loc_id_p )
!      
!enddo
!    
!   
!return

!end subroutine fill_elem_ref_region


!!!**********************************************************************
!!subroutine fill_beta_elem_to_cl_map()
!!!**********************************************************************
!!use cyl_skew_base_ls, only : tree, generation, tr_t, ntree
!!implicit none

!!type(tree), pointer :: tr_t_p
!!type(generation), pointer :: gen_t_p

!!integer :: ncount

!!!  Nullify all pointers
!!nullify(tr_t_p, gen_t_p)

!!nbeta = 0
!!do nt = 1, rns_ntree

!!  !  Point to rns mapped tree
!!  tr_t_p => tr_t( rns_tree_map( nt ) )
!!  
!!  do ng = 1, tr_t_p % ngen_reslv + 1
!!  
!!    gen_t_p => tr_t_p % gen_t( ng )
!!	
!!    do nc = 1, gen_t_p % ncluster 
!!	
!!      nbeta = nbeta + 1
!!    
!!	enddo    
!!	
!!	nullify(gen_t_p)
!!	
!!  enddo
!!  
!!  nullify(tr_t_p)
!!  
!!enddo

!!!  Each resolved cluster recieves a unique id
!!allocate( beta_elem_to_beta_map( nbeta ) )
!!r_elem_to_cl_map = -1

!!ncount=0
!!if(coord == 0) write(*,*) 'Setting r_elem_to_cl_map'

!!do nt = 1, rns_ntree

!!  !  Point to rns mapped tree
!!  tr_t_p => tr_t( rns_tree_map( nt ) )
!!  
!!  do ng = 1, tr_t_p % ngen_reslv
!!  
!!    gen_t_p => tr_t_p % gen_t( ng )
!!	
!!    do nc = 1, gen_t_p % ncluster 
!!	
!!      if(coord == 0) write(*,'(a,3i)') 'rns_tree_map( nt ), ng, nc : ', rns_tree_map( nt ), ng, nc
!!    
!!	  ncount = ncount + 1
!!      
!!	  r_elem_to_cl_map( ncount ) = gen_t_p % cl_t (nc) %indx
!!      
!!	  if(coord == 0) then
!!        write(*,'(1a,5i4)') 'rns_tree_map( nt ), ng, nc, gen_t_p % cl_t (nc) %indx, rns_reslv_cl_indx : ', rns_tree_map( nt ), ng, nc, gen_t_p % cl_t (nc) %indx, ncount
!!      endif
!!    
!!	enddo    
!!	
!!	nullify(gen_t_p)
!!	
!!  enddo
!!  
!!  nullify(tr_t_p)
!!  
!!enddo

!!return
!!end subroutine fill_beta_elem_to_cl_map

!**********************************************************************
subroutine fill_cl_to_r_elem_map()
!**********************************************************************
use messages
use param, only : coord, USE_MPI
implicit none

integer :: n, nt, ng, nc
integer :: ncount

integer, pointer, dimension(:) :: cl_loc_id_p

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(cl_loc_id_p)
nullify(tr_t_p, gen_t_p)

!  Each cluster in the RNS trees recieves a unique id
allocate( cl_to_r_elem_map( ncluster_tot ) )
cl_to_r_elem_map = -1
allocate( r_elem_to_basecl_map( nr_elem ) )
r_elem_to_basecl_map=-1


if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling RNS R_ELEM and Cluster Mapping'
  write(*,*) ' '
endif

!  Need to first set for all trees belonging to the RNS trees
ncount=0
do nt = 1, rns_ntree

  !  Point to RNS mapped tree
  tr_t_p => tr_t( rns_to_cyl_skew_tree_map( nt ) )
  
  do ng = 1, tr_t_p % ngen_reslv
  
    gen_t_p => tr_t_p % gen_t( ng )
    
    do nc = 1, gen_t_p % ncluster 
    
      ncount = ncount + 1
      
      cl_to_r_elem_map( gen_t_p % cl_t(nc) % indx ) = ncount
      r_elem_to_basecl_map( ncount ) = gen_t_p % cl_t(nc) % indx
      
    enddo
    
    nullify( gen_t_p )
    
  enddo
  
  nullify( tr_t_p )
  
enddo

!  Output elem indx to cluster mapping
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> RNS Element Base Cluster ID : '
  write(*,*) '| ID       | NT       | NG       | NC       |'
  do n=1, nr_elem
    cl_loc_id_p => clindx_to_loc_id(:,r_elem_to_basecl_map(n))
    write(*,'(i12,3i11)') n, cl_loc_id_p(1), cl_loc_id_p(2), cl_loc_id_p(3)
    nullify(cl_loc_id_p)
  enddo
endif

return
end subroutine fill_cl_to_r_elem_map


!**********************************************************************
subroutine fill_cl_to_beta_elem_map
!**********************************************************************
use param, only : USE_MPI, coord
implicit none

integer :: ncount, nt, ng, nc, n
integer :: iparent, ng_parent

integer, pointer :: clindx_p, trindx_p
integer, pointer, dimension(:) :: cl_loc_id_p

type(tree), pointer :: tr_t_p 
type(generation), pointer :: gen_t_p

! Nullify all pointers
nullify(tr_t_p, gen_t_p)
nullify(clindx_p)
nullify(trindx_p)
nullify(cl_loc_id_p)

allocate( cl_to_beta_elem_map( ncluster_tot ) ) 
cl_to_beta_elem_map = -1
allocate( beta_elem_to_basecl_map( nbeta_elem ) )
beta_elem_to_basecl_map=-1

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling RNS BETA_ELEM and Cluster Mapping'
  write(*,*) ' '
endif

!  Assign cl_to_beta_elem_map; for a given cluster at gen>=g+1 it returns the
!  the unique rns beta id  
ncount=0  
do nt = 1, rns_ntree
    
  tr_t_p => tr_t( rns_to_cyl_skew_tree_map( nt ) ) 
    
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1 )
    
  do nc = 1, gen_t_p % ncluster ! number of clusters in the g+1 generation
    
    ncount = ncount + 1
      
    !  So clindx_p corresponds to ib
    cl_to_beta_elem_map(  gen_t_p % cl_t(nc) % indx ) = ncount
    beta_elem_to_basecl_map( ncount ) = gen_t_p % cl_t(nc) % indx
    
  enddo
    
  nullify(gen_t_p)
    
  !  Now set the beta id for all of the descendant clusters
  do ng = tr_t_p % ngen_reslv + 2, tr_t_p % ngen
    
    iparent = 0
      
    gen_t_p => tr_t_p % gen_t ( ng ) ! point to the current generation
      
    do nc = 1, gen_t_p % ncluster ! loop over all clusters of gen_t_p
      
      clindx_p => gen_t_p % cl_t(nc) % indx
      iparent = gen_t_p % cl_t(nc) % parent
        
      ng_parent = ng - 1
        
      do while ( ng_parent > tr_t_p % ngen_reslv + 1)
        
        cl_loc_id_p => clindx_to_loc_id(:,iparent)
        trindx_p => rns_to_cyl_skew_tree_map( cyl_skew_to_rns_tree_map( cl_loc_id_p(1) ) )
        
        !  Careful - cl_loc_id_p(1) should point to a RNS tree (set above)
        iparent = tr_t( trindx_p ) % gen_t (cl_loc_id_p(2)) % cl_t(cl_loc_id_p(3)) % parent
          
        ng_parent = cl_loc_id_p(2) - 1
          
        nullify(cl_loc_id_p)
        nullify(trindx_p)
          
      enddo
        
      cl_to_beta_elem_map( clindx_p ) = cl_to_beta_elem_map( iparent ) ! cl_to_beta_elem_map(iparent) is set from above

      nullify(clindx_p)
        
    enddo
      
    nullify(gen_t_p)
      
  enddo
    
  nullify(tr_t_p)
    
enddo

!  Output elem indx to cluster mapping
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> RNS Element Base Cluster ID : '
  write(*,*) '| ID       | NT       | NG       | NC       |'
  do n=1, nbeta_elem
    cl_loc_id_p => clindx_to_loc_id(:,beta_elem_to_basecl_map(n))
    write(*,'(i12,3i11)') n, cl_loc_id_p(1), cl_loc_id_p(2), cl_loc_id_p(3)
    nullify(cl_loc_id_p)
  enddo
endif

return

end subroutine fill_cl_to_beta_elem_map

!**********************************************************************
subroutine fill_cl_to_b_elem_map
!**********************************************************************
use cyl_skew_base_ls, only : tr_t, tree, generation
use param, only : USE_MPI, coord
implicit none

integer :: ncount, nt, ng, nc, n
integer :: iparent, ng_parent

integer, pointer :: clindx_p => null(), trindx_p => null()
integer, pointer, dimension(:) :: cl_loc_id_p

type(tree), pointer :: tr_t_p 
type(generation), pointer :: gen_t_p

! Nullify all pointers
nullify(tr_t_p, gen_t_p)
nullify(cl_loc_id_p)

!  Set the number of beta regions
allocate( cl_to_b_elem_map( ncluster_tot ) )
cl_to_b_elem_map = -1
allocate( b_elem_to_basecl_map( nb_elem ) )
b_elem_to_basecl_map=-1

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling B_ELEM and Cluster Mapping'
  write(*,*) ' '
endif
  
!  Assign rns_rbeta_iarray; for a given cluster at gen=g it returns the
!  the unique rns rbeta id
ncount = 0
do nt = 1, rns_ntree
  
  tr_t_p => tr_t( rns_to_cyl_skew_tree_map( nt ) ) 
   
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv )
    
  do nc = 1, gen_t_p % ncluster ! number of clusters in the g generation
    
    ncount = ncount + 1
      
    clindx_p => gen_t_p % cl_t(nc) % indx
      
    !  So clindx_p corresponds to ncount
    cl_to_b_elem_map( clindx_p ) = ncount
    b_elem_to_basecl_map( ncount ) = clindx_p
      
    nullify(clindx_p)
    
  enddo
  
  nullify( gen_t_p )
    
   !  Now set the b id for all of the descendant clusters
  do ng = tr_t_p % ngen_reslv + 1, tr_t_p % ngen
    
    iparent = 0
      
    gen_t_p => tr_t_p % gen_t ( ng ) ! point to the current generation
      
    do nc = 1, gen_t_p % ncluster ! loop over all clusters of gen_t_p
      
      clindx_p => gen_t_p % cl_t(nc) % indx
      iparent = gen_t_p % cl_t(nc) % parent
        
      ng_parent = ng - 1
        
      do while ( ng_parent > tr_t_p % ngen_reslv)
        
        cl_loc_id_p => clindx_to_loc_id(:,iparent)
        trindx_p => rns_to_cyl_skew_tree_map( cyl_skew_to_rns_tree_map( cl_loc_id_p(1) ) )  

        iparent = tr_t( trindx_p ) % gen_t (cl_loc_id_p(2)) % cl_t(cl_loc_id_p(3)) % parent
          
        ng_parent = cl_loc_id_p(2) - 1
          
        nullify(cl_loc_id_p)
        nullify(trindx_p)
          
      enddo
        
      cl_to_b_elem_map( clindx_p ) = cl_to_b_elem_map( iparent ) ! cl_to_beta_elem_map(iparent) is set from above

      nullify(clindx_p)
        
    enddo
      
    nullify(gen_t_p)
      
  enddo
    
  nullify(tr_t_p)
    
enddo

!  Output elem indx to cluster mapping
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> RNS Element Base Cluster ID : '
  write(*,*) '| ID       | NT       | NG       | NC       |'
  do n=1, nb_elem
    cl_loc_id_p => clindx_to_loc_id(:,b_elem_to_basecl_map(n))
    write(*,'(i12,3i11)') n, cl_loc_id_p(1), cl_loc_id_p(2), cl_loc_id_p(3)
    nullify(cl_loc_id_p)
  enddo
endif

return
end subroutine fill_cl_to_b_elem_map

!**********************************************************************
subroutine fill_r_elem()
!**********************************************************************
!  This subroutine allocates and initializes all information needed to 
!  use the resolved elements (r_elem) in the RNS module
!
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_r_elem'
integer :: n

allocate( r_elem_t ( nr_elem ) )

!  Fill reference region
call fill_r_elem_ref_region()
!  Fill index array
call fill_r_elem_indx_array()
!  Initialize force data
do n=1, nr_elem
  r_elem_t( n ) % force_t = force_type_1( fx=0._rprec, fy=0._rprec, CD=0._rprec, kappa=0._rprec )
enddo

return

end subroutine fill_r_elem

!**********************************************************************
subroutine fill_r_elem_ref_region()
!**********************************************************************
use param, only : dx, dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_r_elem_ref_region'

integer :: n, nb
integer :: nxi, neta, nzeta

real(rprec) :: h, w, xi_c(3)
real(rprec), dimension(3) :: p1, p2, p3, p4

integer, pointer, dimension(:) :: cl_loc_id_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p

type(cluster),    pointer :: cl_t_p
type(branch),     pointer :: br_t_p

type(ref_region), pointer :: ref_region_t_p

nullify(cl_loc_id_p)
nullify(d_p, l_p, skew_angle_p)
nullify(hbot_p, htop_p)
nullify(cl_t_p, br_t_p)
nullify(ref_region_t_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling R_ELEM Reference Plane Arrays'
  write(*,*) ' '
  write(*,*) '----> Reference Region Information : '
  write(*,'(a116)') '| ID       | NPOINT   | AREA     | P1.X   | P1.Y   | P1.Z   | P2.X   | P2.Y   | P2.Z   | P3.X   | P3.Y   | P3.Z   |'
endif

do n=1, nr_elem
  
  !  Get the base cluster local id
  cl_loc_id_p => clindx_to_loc_id(:, r_elem_to_basecl_map(n) )
  
  cl_t_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % cl_t( cl_loc_id_p(3) )

  !  Point to the top and bottom of the plane
  hbot_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % bplane
  htop_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % tplane

  h = htop_p - hbot_p
  w = alpha_width * h
 
  nxi   = ceiling ( 2. * alpha_dist * h / dx ) + 1
  neta  = ceiling( w / dy ) + 1
  nzeta = ceiling( h / dz ) + 1  
      
  !  Offset in the upstream x-direction
  xi_c = cl_t_p % origin + (/ -alpha_dist * h, 0._rprec, 0._rprec /)
      
  !  Set the ordered corner points of the box
  p1    = xi_c 
  p1(2) = p1(2) + w / 2._rprec
    
  p2    = p1
  p2(2) = p2(2) - w
      
  p3    = p2
  p3(3) = p3(3) + h
  
  p4    = p2
  p4(1) = p4(1) + 2. * alpha_dist * h
      
  !  Point to ref_region_t of the resolved element
  ref_region_t_p => r_elem_t ( n ) % ref_region_t

  ref_region_t_p % area    = h * w
  ref_region_t_p % npoint = nxi * neta * nzeta
    
  !  Check if the element has been allocated
  !if( associated( ref_region_t_p % points ) ) then
  !  call error(sub_name, 'reference region points already allocated.')
  !else
  allocate( ref_region_t_p % points( 3, ref_region_t_p % npoint ) )
  !endif
        
  !call set_points_in_plane( p1, p2, p3, nxi, neta, ref_region_t_p % points )  
  call set_points_in_box( p1, p2, p3, p4, nxi, neta, nzeta, ref_region_t_p % points ) 

  !  Finally initialize velocity reference components
  ref_region_t_p % u = 0._rprec
  
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,'(i12,i11,f11.6, 9f9.4)') n, r_elem_t(n) % ref_region_t % npoint, r_elem_t(n) % ref_region_t % area, p1(:), p2(:), p3(:)
  endif   
      
  nullify( ref_region_t_p )
  nullify( htop_p, hbot_p )
  nullify( cl_t_p )
  nullify( cl_loc_id_p )
      
enddo
    
   
return

end subroutine fill_r_elem_ref_region

!**********************************************************************
subroutine fill_r_elem_indx_array
!**********************************************************************
!  This subroutine gets the indx_array for r_elem. 
!
use types, only : rprec
use param, only : nx,ny,nz, coord, USE_MPI
use level_set_base, only : phi
use messages
$if($MPI)
use mpi
use param, only : comm, ierr
$endif

implicit none

character (*), parameter :: sub_name = mod_name // '.fill_r_elem_indx_array'

type(indx_array), target, allocatable, dimension(:) :: pre_indx_array_t

integer :: i,j,k, n, np

integer, allocatable, dimension(:) :: npoint_global

integer, pointer :: clindx_p, indx_p, npoint_p, trindx_p
integer, pointer, dimension(:) :: cl_loc_id_p

! ---- Nullify all pointers ----
nullify(clindx_p)
nullify(trindx_p)
nullify(indx_p)
nullify(npoint_p)
nullify(cl_loc_id_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling R_ELEM Index Array'
  write(*,*) ' '
endif

! ---- Allocate all arrays ----
allocate( pre_indx_array_t( nr_elem ) )

do n=1, nr_elem
  allocate(pre_indx_array_t(n) % iarray(3,nx*ny*(nz-1)))
enddo

!  Intialize the number of points assigned to the cluster
do n=1, nr_elem
  pre_indx_array_t(n) % npoint = 0
enddo

if(.not. chi_initialized) call error(sub_name, 'chi not initialized')

do k=1, nz - 1

  do j=1, ny

    do i = 1, nx

      if ( clindx(i,j,k) > 0 ) then
      
        !  Need to map cluster to coorespond cluster in RNS tree
        cl_loc_id_p => clindx_to_loc_id(:, clindx(i,j,k))
        trindx_p => rns_to_cyl_skew_tree_map( cyl_skew_to_rns_tree_map( cl_loc_id_p(1) ) )

        clindx_p => tr_t( trindx_p ) % gen_t( cl_loc_id_p(2) ) % cl_t ( cl_loc_id_p(3) ) % indx
 
        indx_p => cl_to_r_elem_map( clindx_p )

        nullify( trindx_p )
      
        !  Check if cluster belongs to a r_elem
        if( indx_p > 0 ) then

          !  Use only inside points
          if ( phi(i,j,k) <= 0._rprec ) then 


            pre_indx_array_t( indx_p ) % npoint = pre_indx_array_t( indx_p ) % npoint + 1

            npoint_p => pre_indx_array_t( indx_p ) % npoint
        
            pre_indx_array_t( indx_p ) % iarray(1, npoint_p ) = i
            pre_indx_array_t( indx_p ) % iarray(2, npoint_p ) = j
            pre_indx_array_t( indx_p ) % iarray(3, npoint_p ) = k
          
            nullify( npoint_p )

          endif
  
        endif

        nullify( indx_p )
        nullify(clindx_p)
        nullify(cl_loc_id_p)        
       
      endif
      
    enddo
    
  enddo
  
enddo

!  Now set set index array for each r_elem
do n=1, nr_elem

  r_elem_t(n) % indx_array_t % npoint = pre_indx_array_t(n) % npoint
  
  allocate(r_elem_t(n) % indx_array_t % iarray(3, r_elem_t(n) % indx_array_t % npoint ) )
  
  do np=1, r_elem_t(n) % indx_array_t % npoint
    r_elem_t(n) % indx_array_t % iarray(:, np) = pre_indx_array_t( n ) % iarray(:, np)
  enddo
    
enddo

deallocate(pre_indx_array_t)

allocate(npoint_global( nr_elem ))
$if($MPI)
call mpi_allreduce(r_elem_t(:) % indx_array_t % npoint, npoint_global(:), nr_elem, MPI_INTEGER, MPI_SUM, comm, ierr)
$else
npoint_global(:) = r_elem_t(:) % indx_array_t % npoint
$endif

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> Index Array Information : '
  write(*,*) '| ID       | NPOINT   |'
  do n=1, nr_elem
    write(*,'(i12,i11)') n, npoint_global(n)
  enddo
endif
deallocate(npoint_global)

return
end subroutine fill_r_elem_indx_array

!**********************************************************************
subroutine fill_beta_elem()
!**********************************************************************
!  This subroutine allocates and initializes all information needed to 
!  use the BETA elements (beta_elem) in the RNS module
!
implicit none
character (*), parameter :: sub_name = mod_name // '.fill_beta_elem'
integer :: n

allocate( beta_elem_t ( nbeta_elem ) )

!  Fill reference region
call fill_beta_elem_ref_region()
!  Fill index array
call fill_beta_elem_indx_array()
!  Initialize force data
do n=1, nbeta_elem
  beta_elem_t( n ) % force_t = force_type_1( fx=0._rprec, fy=0._rprec, CD=0._rprec, kappa=0._rprec )
enddo

return

end subroutine fill_beta_elem

!**********************************************************************
subroutine fill_beta_elem_ref_region()
!**********************************************************************
use param, only : dx, dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_beta_elem_ref_region'

integer :: n
integer :: nxi, neta, nzeta

real(rprec) :: h, w
real(rprec), dimension(3) :: p1, p2, p3,p4, xi_c

integer, pointer, dimension(:) :: cl_loc_id_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p

type(cluster),    pointer :: cl_t_p

type(ref_region), pointer :: ref_region_t_p

nullify(cl_loc_id_p)
nullify(d_p, l_p, skew_angle_p)
nullify(hbot_p, htop_p)
nullify(cl_t_p)
nullify(ref_region_t_p)
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling BETA_ELEM Reference Plane Arrays'
  write(*,*) ' '
  write(*,*) '----> Reference Region Information : '
  write(*,'(a116)') '| ID       | NPOINT   | AREA     | P1.X   | P1.Y   | P1.Z   | P2.X   | P2.Y   | P2.Z   | P3.X   | P3.Y   | P3.Z   |'    
endif

do n=1, nbeta_elem
  
  !  Get the base cluster local id
  cl_loc_id_p => clindx_to_loc_id(:, beta_elem_to_basecl_map(n) )
  
  cl_t_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % cl_t( cl_loc_id_p(3) )

  !  Point to the top and bottom of the plane
  hbot_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % bplane
  htop_p => tr_t( cl_loc_id_p(1) ) % gen_t( ngen ) % tplane
   
  h = htop_p - hbot_p
  w = alpha_width * h
    
  nxi   = ceiling ( 2. * alpha_dist * h / dx ) + 1
  neta  = ceiling( w / dy ) + 1
  nzeta = ceiling( h / dz ) + 1  
      
  !  Offset in the upstream x-direction
  xi_c = cl_t_p % origin + (/ -alpha_dist * h, 0._rprec, 0._rprec /)
      
  !  Set the ordered corner points of the box
  p1    = xi_c 
  p1(2) = p1(2) + w / 2._rprec
    
  p2    = p1
  p2(2) = p2(2) - w
      
  p3    = p2
  p3(3) = p3(3) + h
  
  p4    = p2
  p4(1) = p4(1) + 2. * alpha_dist * h

  !  Point to ref_region_t of the resolved element
  ref_region_t_p => beta_elem_t ( n ) % ref_region_t

  ref_region_t_p % area = h * w
  ref_region_t_p % npoint = nxi*neta*nzeta
      
    !!  Check if the element has been allocated
    !if( associated( ref_region_t_p % points ) ) then
    !  call error(sub_name, 'reference region points already allocated.')
    !else
  allocate( ref_region_t_p % points( 3, ref_region_t_p % npoint ) )
    !endif
            
  !call set_points_in_plane( p1, p2, p3, nxi, neta, ref_region_t_p % points )  
  call set_points_in_box( p1, p2, p3, p4, nxi, neta, nzeta, ref_region_t_p % points ) 

  !  Finally initialize velocity reference components
  ref_region_t_p % u = 0._rprec
  
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,'(i12,i11,f11.6, 9f9.4)') n, beta_elem_t(n) % ref_region_t % npoint, beta_elem_t(n) % ref_region_t % area, p1(:), p2(:), p3(:)
  endif     
      
  nullify( ref_region_t_p )
  nullify( htop_p, hbot_p)
  nullify( cl_t_p )
  nullify(cl_loc_id_p)
      
enddo
          
return

end subroutine fill_beta_elem_ref_region

!**********************************************************************
subroutine fill_beta_elem_indx_array
!**********************************************************************
!  This subroutine gets the indx_array for beta_elem. 
!
use types, only : rprec
use param, only : nx,ny,nz, coord, USE_MPI
use level_set_base, only : phi
use messages
$if($MPI)
use mpi
use param, only : comm, ierr
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_beta_elem_indx_array'

type(indx_array), target, allocatable, dimension(:) :: pre_indx_array_t

integer :: i,j,k, n, np

integer, allocatable, dimension(:) :: npoint_global

integer, pointer :: clindx_p, indx_p, npoint_p, trindx_p
integer, pointer, dimension(:) :: cl_loc_id_p

! ---- Nullify all pointers ----
nullify(clindx_p)
nullify(trindx_p)
nullify(indx_p)
nullify(npoint_p)
nullify(cl_loc_id_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling BETA_ELEM Index Array'
  write(*,*) ' '
endif

! ---- Allocate all arrays ----
allocate( pre_indx_array_t( nbeta_elem ) )

do n=1, nbeta_elem
  allocate(pre_indx_array_t(n) % iarray(3,nx*ny*(nz-1)))
enddo

!  Intialize the number of points assigned to the cluster
do n=1, nbeta_elem
  pre_indx_array_t(n) % npoint = 0
enddo

if(.not. chi_initialized) call error(sub_name, 'chi not initialized')

do k=1, nz - 1

  do j=1, ny

    do i = 1, nx

           
      if ( clindx(i,j,k) > 0 ) then
      
        !  Need to map cluster to coorespond cluster in RNS tree
        cl_loc_id_p => clindx_to_loc_id(:, clindx(i,j,k))
        trindx_p => rns_to_cyl_skew_tree_map( cyl_skew_to_rns_tree_map( cl_loc_id_p(1) ))

        clindx_p => tr_t( trindx_p ) % gen_t( cl_loc_id_p(2) ) % cl_t ( cl_loc_id_p(3) ) % indx
  
        indx_p => cl_to_beta_elem_map( clindx_p )

        nullify(trindx_p)
      
        !  Check if cluster belongs to beta_elem
        if( indx_p > 0 ) then

          !  Use only points within cutoff
          if ( chi(i,j,k) >= chi_cutoff ) then 

            pre_indx_array_t( indx_p ) % npoint = pre_indx_array_t( indx_p ) % npoint + 1

            npoint_p => pre_indx_array_t( indx_p ) % npoint
        
            pre_indx_array_t( indx_p ) % iarray(1, npoint_p ) = i
            pre_indx_array_t( indx_p ) % iarray(2, npoint_p ) = j
            pre_indx_array_t( indx_p ) % iarray(3, npoint_p ) = k
          
            nullify( npoint_p )

          endif
  
        endif

        nullify( indx_p )
       
      endif
      
      nullify(clindx_p)
	  nullify(cl_loc_id_p)
      
    enddo
    
  enddo
  
enddo

!  Now set set index array for each beta_elem
do n=1, nbeta_elem

  beta_elem_t(n) % indx_array_t % npoint = pre_indx_array_t(n) % npoint
  
  allocate(beta_elem_t(n) % indx_array_t % iarray(3, beta_elem_t(n) % indx_array_t % npoint ) )
  
  do np=1, beta_elem_t(n) % indx_array_t % npoint
    beta_elem_t(n) % indx_array_t % iarray(:, np) = pre_indx_array_t( n ) % iarray(:, np)
  enddo
    
enddo

deallocate(pre_indx_array_t)

allocate(npoint_global( nbeta_elem ))
$if($MPI)
call mpi_allreduce(beta_elem_t(:) % indx_array_t % npoint, npoint_global(:), nbeta_elem, MPI_INTEGER, MPI_SUM, comm, ierr)
$else
npoint_global(:) = beta_elem_t(:) % indx_array_t % npoint
$endif

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> Index Array Information : '
  write(*,*) '| ID       | NPOINT   |'
  do n=1, nbeta_elem
    write(*,'(i12,i11)') n, npoint_global(n)
  enddo
endif

return
end subroutine fill_beta_elem_indx_array

!**********************************************************************
subroutine fill_b_elem()
!**********************************************************************
!  This subroutine allocates and initializes all information needed to 
!  use the B elements (b_elem) in the RNS module
!
implicit none
character (*), parameter :: sub_name = mod_name // '.fill_b_elem'

integer :: n

allocate( b_elem_t ( nb_elem ) )

!  Fill reference region
call fill_b_elem_ref_region()
! Set all child elements belonging to b_elem_t
call set_b_elem_children()
!  Initialize force data
do n=1, nb_elem
  b_elem_t( n ) % force_t = force_type_2( fx=0._rprec, fy=0._rprec, CD=0._rprec, &
    error=0._rprec, LAB=0._rprec, LBB=0._rprec )
enddo

return
end subroutine fill_b_elem

!**********************************************************************
subroutine fill_b_elem_ref_region()
!**********************************************************************
use param, only : dx, dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_b_elem_ref_region'

integer :: n
integer :: nxi, neta, nzeta

real(rprec) :: h, w
real(rprec), dimension(3) :: p1, p2, p3, p4, xi_c

integer, pointer, dimension(:) :: cl_loc_id_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p

type(tree),       pointer :: tr_t_p
type(generation), pointer :: gen_t_p
type(cluster),    pointer :: cl_t_p
type(branch),     pointer :: br_t_p

type(ref_region), pointer :: ref_region_t_p

nullify(cl_loc_id_p)
nullify(d_p, l_p, skew_angle_p)
nullify(hbot_p, htop_p)
nullify(tr_t_p, gen_t_p, cl_t_p, br_t_p)
nullify(ref_region_t_p)
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filling B_ELEM Reference Plane Arrays'
  write(*,*) ' '
  write(*,*) '----> Reference Region Information : '
  write(*,'(a116)') '| ID       | NPOINT   | AREA     | P1.X   | P1.Y   | P1.Z   | P2.X   | P2.Y   | P2.Z   | P3.X   | P3.Y   | P3.Z   |'  
endif

do n=1, nb_elem
  
  !  Point to ref_region_t of the resolved element
  ref_region_t_p => b_elem_t ( n ) % ref_region_t

  !  Get the base cluster local id
  cl_loc_id_p => clindx_to_loc_id(:, b_elem_to_basecl_map(n) )
  
  cl_t_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % cl_t( cl_loc_id_p(3) )

  !  Point to the top and bottom of the plane
  hbot_p => tr_t( cl_loc_id_p(1) ) % gen_t( cl_loc_id_p(2) ) % bplane
  htop_p => tr_t( cl_loc_id_p(1) ) % gen_t( ngen ) % tplane

  h = htop_p - hbot_p
  w = alpha_width * h
    
  nxi   = ceiling ( 2. * alpha_dist * h / dx ) + 1
  neta  = ceiling( w / dy ) + 1
  nzeta = ceiling( h / dz ) + 1  
      
  !  Offset in the upstream x-direction
  xi_c = cl_t_p % origin + (/ -alpha_dist * h, 0._rprec, 0._rprec /)
      
  !  Set the ordered corner points of the box
  p1    = xi_c 
  p1(2) = p1(2) + w / 2._rprec
    
  p2    = p1
  p2(2) = p2(2) - w
      
  p3    = p2
  p3(3) = p3(3) + h
  
  p4    = p2
  p4(1) = p4(1) + 2. * alpha_dist * h


  ref_region_t_p % area = h * w
  ref_region_t_p % npoint = nxi*neta*nzeta
      
  !  Check if the element has been allocated
  !if( associated( ref_region_t_p % points ) ) then
  !  call error(sub_name, 'reference region points already allocated.')
  !else
  allocate( ref_region_t_p % points( 3, ref_region_t_p % npoint ) )
  !endif
            
  !call set_points_in_plane( p1, p2, p3, nxi, neta, ref_region_t_p % points )  
  call set_points_in_box( p1, p2, p3, p4, nxi, neta, nzeta, ref_region_t_p % points )    

    !  Finally initialize velocity reference components
    ref_region_t_p % u = 0._rprec
    
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,'(i12,i11,f11.6, 9f9.4)') n, b_elem_t(n) % ref_region_t % npoint, b_elem_t(n) % ref_region_t % area, p1(:), p2(:), p3(:)
  endif      
      
  nullify( ref_region_t_p )
	nullify( htop_p, hbot_p )
  nullify( cl_t_p )
  nullify( cl_loc_id_p )

enddo
         
return

end subroutine fill_b_elem_ref_region

!**********************************************************************
subroutine set_b_elem_children()
!**********************************************************************
!  This subroutine sets the r_elem and beta_elem (children) which belong to 
!  each B region
use param, only : dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_name = mod_name // '.set_b_elem_children'

integer ::  n

integer, pointer :: clindx_p, b_elem_indx_p

type(child_elem), allocatable, dimension(:) :: pre_r_child_t
type(child_elem), allocatable, dimension(:) :: pre_beta_child_t

nullify(clindx_p, b_elem_indx_p)

allocate(pre_r_child_t( nb_elem ))
allocate(pre_beta_child_t( nb_elem ))

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Setting B_ELEM Children'
  write(*,*) ' '
endif

!  Initialize temporary children
do n=1, nb_elem

  pre_r_child_t(n) % nelem = 0
  allocate(pre_r_child_t(n) % indx( nr_elem ) )
  pre_r_child_t(n) % indx(:) = -1

enddo

!  Loop over all r_elem
do n=1, nr_elem

  !  Get the base cluster id
  clindx_p => r_elem_to_basecl_map(n) 
  
  !  Check if r_elem base cluster belongs to any b_elem
  b_elem_indx_p => cl_to_b_elem_map( clindx_p )
  
  if( b_elem_indx_p  > 0 ) then
    pre_r_child_t( b_elem_indx_p ) % nelem = pre_r_child_t( b_elem_indx_p ) % nelem + 1
    pre_r_child_t( b_elem_indx_p ) % indx( pre_r_child_t( b_elem_indx_p ) % nelem ) = n 
  endif
  
  nullify( b_elem_indx_p )
  nullify( clindx_p )
  
enddo

!  Now allocate and set the actual children
do n=1, nb_elem
  allocate( b_elem_t( n ) % r_child_t % indx( pre_r_child_t( n ) % nelem ))
  b_elem_t( n ) % r_child_t % nelem = pre_r_child_t( n ) % nelem
  b_elem_t( n ) % r_child_t % indx(1:b_elem_t( n ) % r_child_t % nelem) = pre_r_child_t( n ) % indx(1:b_elem_t( n ) % r_child_t % nelem)
enddo

deallocate( pre_r_child_t )

!  Initialize temporary children
do n=1, nb_elem
  
  pre_beta_child_t(n) % nelem = 0
  allocate(pre_beta_child_t(n) % indx( nbeta_elem ) )
  pre_beta_child_t(n) % indx( : ) = -1

enddo

!  Loop over all beta_elem
do n=1, nbeta_elem

  !  Get the base cluster id
  clindx_p => beta_elem_to_basecl_map(n) 
  
  !  Check if r_elem base cluster belongs to any b_elem
  b_elem_indx_p => cl_to_b_elem_map( clindx_p )
  
  if( b_elem_indx_p  > 0 ) then
    pre_beta_child_t( b_elem_indx_p ) % nelem = pre_beta_child_t( b_elem_indx_p ) % nelem + 1
    pre_beta_child_t( b_elem_indx_p ) % indx( pre_beta_child_t( b_elem_indx_p ) % nelem ) = n 
  endif
  nullify( b_elem_indx_p )
  nullify( clindx_p )
  
enddo


!  Now allocate and set the actual children
do n=1, nb_elem
  allocate( b_elem_t( n ) % beta_child_t % indx( pre_beta_child_t( n ) % nelem ))
  b_elem_t( n ) % beta_child_t % nelem = pre_beta_child_t( n ) % nelem
  b_elem_t( n ) % beta_child_t % indx(1:b_elem_t( n ) % beta_child_t % nelem) = pre_beta_child_t( n ) % indx(1:b_elem_t( n ) % beta_child_t % nelem)
enddo

deallocate( pre_beta_child_t )

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '----> Children Element Information : '
  write(*,*) '| ID       | NR_ELEM  | CHILD ID         |'
  do n=1, nb_elem
    write(*,'(i12,i11,i5)') n, b_elem_t(n) % r_child_t % nelem, b_elem_t(n) % r_child_t % indx(:)
  enddo
  write(*,*) '| ID       | NBETA_ELEM | CHILD ID         |'
  do n=1, nb_elem
    write(*,'(i12,i13,3i5)') n, b_elem_t(n) % beta_child_t % nelem, b_elem_t(n) % beta_child_t % indx(:)
  enddo  
endif

return
end subroutine set_b_elem_children

!**********************************************************************
subroutine set_points_in_plane(bp1, bp2, bp3, nzeta, neta, points)
!**********************************************************************
!
!  This subroutine assigns the points in an arbitrary 3D plane
!

use types, only : rprec
use param, only : Nx, Ny, Nz, dx, dy, dz, L_x, L_y
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif
use grid_defs
use messages
implicit none

real(RPREC), intent(IN), dimension(:) :: bp1, bp2, bp3
INTEGER, INTENT(IN) :: nzeta, neta

real(rprec), intent(out), dimension(3,nzeta*neta) :: points

character (*), parameter :: sub_name = mod_name // '.set_points_in_plane'

integer :: i, j, np

REAL(RPREC) :: dzeta, deta, vec_mag

real(RPREC), dimension(3) :: zeta_vec, eta_vec, eta
real(RPREC), dimension(3) :: cell_center

points(:,:) = -huge(1.) ! Initialize to some bogus value

!  vector in zeta direction
zeta_vec = bp1 - bp2
!  vector in eta direction
eta_vec   = bp3 - bp2

!  Normalize to create unit vector
vec_mag = sqrt(zeta_vec(1)*zeta_vec(1) + zeta_vec(2)*zeta_vec(2) + zeta_vec(3)*zeta_vec(3))
dzeta = vec_mag/nzeta
zeta_vec = zeta_vec / vec_mag

vec_mag = sqrt(eta_vec(1)*eta_vec(1) + eta_vec(2)*eta_vec(2) + eta_vec(3)*eta_vec(3))
deta = vec_mag/neta
eta_vec = eta_vec / vec_mag

np=0
!  Compute cell centers
do j=1,neta
  !  Attempt for cache friendliness
  eta = (j - 0.5)*deta*eta_vec
  do i=1,nzeta
  
    np = np + 1
    
    ! Simple vector addition
    cell_center = bp2 + (i - 0.5)*dzeta*zeta_vec + eta
        
    points(:,np) = cell_center

  enddo
enddo

return

end subroutine set_points_in_plane

!**********************************************************************
subroutine set_points_in_box(bp1, bp2, bp3, bp4, nxi, neta, nzeta, points)
!**********************************************************************
!
!  This subroutine assigns the points in an arbitrary 3D box
!

use types, only : rprec
use param, only : Nx, Ny, Nz, dx, dy, dz, L_x, L_y
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif
use grid_defs
use messages
implicit none

real(RPREC), intent(IN), dimension(3) :: bp1, bp2, bp3, bp4
INTEGER, INTENT(IN) :: nxi, neta, nzeta

real(rprec), intent(out), dimension(3,nxi*neta*nzeta) :: points

character (*), parameter :: sub_name = mod_name // '.set_points_in_box'

integer :: i, j, k, np

REAL(RPREC) :: dxi, deta, dzeta, vec_mag

real(RPREC), dimension(3) :: xi_vec, eta_vec, zeta_vec, eta, zeta
real(RPREC), dimension(3) :: cell_center

points(:,:) = -huge(1.) ! Initialize to some bogus value

! bp2 serves as the local origin
!  vector in xi direction
xi_vec = bp4 - bp2
!  vector in eta direction
eta_vec   = bp1 - bp2
!  vector in zeta direction
zeta_vec   = bp3 - bp2

!  Normalize to create unit vectors
vec_mag = sqrt(xi_vec(1)*xi_vec(1) + xi_vec(2)*xi_vec(2) + xi_vec(3)*xi_vec(3))
dxi = vec_mag/nxi
xi_vec = xi_vec / vec_mag

vec_mag = sqrt(eta_vec(1)*eta_vec(1) + eta_vec(2)*eta_vec(2) + eta_vec(3)*eta_vec(3))
deta = vec_mag/neta
eta_vec = eta_vec / vec_mag

vec_mag = sqrt(zeta_vec(1)*zeta_vec(1) + zeta_vec(2)*zeta_vec(2) + zeta_vec(3)*zeta_vec(3))
dzeta = vec_mag/nzeta
zeta_vec = zeta_vec / vec_mag

np=0
do k=1,nzeta
  !  Compute cell centers
  !  Attempt for cache friendliness
  zeta = (k - 0.5)*dzeta*zeta_vec  
  do j=1,neta
    !  Attempt for cache friendliness
    eta = (j - 0.5)*deta*eta_vec
    do i=1,nxi
  
      np = np + 1
    
      ! Simple vector addition
      cell_center = bp2 + (i - 0.5)*dxi*xi_vec + eta + zeta
      ! Autowrap point
      cell_center(1) = modulo(cell_center(1), L_x)
      cell_center(2) = modulo(cell_center(2), L_y) 
      points(:,np) = cell_center

    enddo

  enddo
enddo

return

end subroutine set_points_in_box

end module rns_cyl_skew_ls
