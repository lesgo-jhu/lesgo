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
use rns_base_ls
use cyl_skew_ls, only : fill_tree_array_ls

implicit none

save
private

public rns_init_ls

integer, parameter :: rns_ntree = 2 ! Number of unique trees 
integer, parameter :: rns_tree_layout = 1

integer :: ncluster_reslv ! total number of resolved clusters
integer :: ncluster_tot

integer, pointer, dimension(:) :: cl_to_r_elem_map
integer, pointer, dimension(:) :: cl_to_beta_elem_map
integer, pointer, dimension(:) :: cl_to_b_elem_map

integer, pointer, dimension(:) :: rns_tree_map(:) ! This maps the tree number from cyl_skew to the trees considered during rns

type(indx_array), pointer, dimension(:) :: pre_beta_elem_indx_array_t

contains

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************
use messages
use param, only : USE_MPI, coord
use cyl_skew_base_ls, only : tree, generation, tr_t, ntree, clindx_to_loc_id
use cyl_skew_ls, only : fill_tree_array_ls

implicit none

character (*), parameter :: sub_name = mod_name // '.rns_init_ls'

integer, pointer :: clindx_p
integer, pointer, dimension(:) :: cl_loc_id_p
type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(clindx_p)
nullify(cl_loc_id_p)
nullify(tr_t_p, gen_t_p)

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
!----- Load Data -----

!----- Fill CYL_SKEW Data Structures -----
!  Fill the cyl_skew tree array
call fill_tree_array_ls()
!  Get the number of resolved clusters 
call set_ncluster_reslv()
!  Get the total number of clusters
call set_ncluster_tot()
!----- Fill CYL_SKEW Data Structures -----

!----- Fill RNS_CYL_SKEW Data Structures -----
!  Fill the mapping for the rns trees
call fill_rns_tree_map()
!  Set the number of r elements
call set_nr_elem()
!  Set the number of beta elements
call set_nbeta_elem()
!  Set the number of b elements
call set_nb_elem()
!!  Fill the mapping of resolved rns element ids to cyl_skew cluster ids
!call fill_r_elem_to_cl_map()
!!  Fill the mapping of the beta region from which to base the reference region calculations on
!call fill_beta_elem_to_cl_map()
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

!!  Set the number of beta regions
!allocate( rns_beta_iarray( ncluster_tot ) ) 
!allocate( rns_rbeta_iarray( ncluster_tot ) )
!  
!rns_beta_iarray = -1
!rns_rbeta_iarray = -1

!nbeta = 0
!nrbeta = 0
!do nt = 1, rns_ntree

!  ng_beta = tr_t(rns_tree_iarray(nt)) % ngen_reslv + 1
!  ng_rbeta = ng_beta - 1
!    
!  nbeta = nbeta + tr_t(rns_tree_iarray(nt)) % gen_t(ng_beta) % ncluster
!  nrbeta = nrbeta + tr_t(rns_tree_iarray(nt)) % gen_t(ng_rbeta) % ncluster
!  
!enddo
!  
!!  Assign rns_rbeta_iarray; for a given cluster at gen=g it returns the
!!  the unique rns rbeta id
!irb = 0
!do nt = 1, rns_ntree
!  
!  tr_t_p => tr_t( rns_tree_iarray( nt ) ) 
!   
!  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv )
!    
!  do nc = 1, gen_t_p % ncluster ! number of clusters in the g generation
!    
!    irb = irb + 1
!      
!    clindx_p => gen_t_p % cl_t(nc) % indx
!      
!    !  So clindx_p corresponds to ib
!    rns_rbeta_iarray(clindx_p) = irb
!      
!    nullify(clindx_p)
!    
!  enddo
!    
!  nullify(gen_t_p, tr_t_p)    
!    
!enddo
!  
!!  Assign rns_beta_iarray; for a given cluster at gen=g it returns the
!!  the unique rns beta id  
!ib=0  
!do nt = 1,rns_ntree
!    
!  tr_t_p => tr_t( rns_tree_iarray( nt ) ) 
!    
!  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1 )
!    
!  do nc = 1, gen_t_p % ncluster ! number of clusters in the g+1 generation
!    
!    ib = ib + 1
!      
!    clindx_p => gen_t_p % cl_t(nc) % indx
!      
!    !  So clindx_p corresponds to ib
!    rns_beta_iarray(clindx_p) = ib
!      
!    nullify(clindx_p)
!    
!  enddo
!    
!  nullify(gen_t_p)
!    
!  !  Now set the beta id for all of the descendant clusters
!  do ng = tr_t_p % ngen_reslv + 2, tr_t_p % ngen
!    
!    iparent = 0
!      
!    gen_t_p => tr_t_p % gen_t ( ng ) ! point to the current generation
!      
!    do nc = 1, gen_t_p % ncluster ! loop over all clusters of gen_t_p
!      
!      clindx_p => gen_t_p % cl_t(nc) % indx
!      iparent = gen_t_p % cl_t(nc) % parent
!        
!      ng_parent = ng - 1
!        
!      do while ( ng_parent > tr_t_p % ngen_reslv + 1)
!        
!        cl_loc_id_p => clindx_to_loc_id(:,iparent)
!          
!        iparent = tr_t(rns_tree_iarray(cl_loc_id_p(1))) % gen_t (cl_loc_id_p(2)) % &
!          cl_t(cl_loc_id_p(3)) % parent
!          
!        ng_parent = cl_loc_id_p(2) - 1
!          
!        nullify(cl_loc_id_p)
!          
!      enddo
!        
!      rns_beta_iarray(clindx_p) = rns_beta_iarray(iparent) ! rns_beta_iarray(iparent) is set from above

!      nullify(clindx_p)
!        
!    enddo
!      
!    nullify(gen_t_p)
!      
!  enddo
!    
!  nullify(tr_t_p)
!    
!enddo

!call rns_fill_ref_plane_array_ls()
!  Create cluster index array
!call rns_fill_indx_array_ls()
  
!allocate( clforce_t ( ncluster_reslv ) ) 
!clforce_t = force(parent = 0, CD = 0._rprec, fD = 0._rprec, kappa=0._rprec)
  
!allocate( beta_force_t( nbeta ) )
!beta_force_t = force(parent = 0, CD = 0._rprec, fD = 0._rprec, kappa=0._rprec)
  
!call rns_force_init_ls()
  
!call rns_set_parent_ls()

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'RNS Data Structure Initialized'
  write(*,*) ' '
endif 

return

end subroutine rns_init_ls

!**********************************************************************
subroutine clindx_init ()
!**********************************************************************
use param, only : iBOGUS, coord, ld, ny, nz
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.clindx_init'
character (*), parameter :: fname_in = 'clindx.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fname_in_MPI
$endif

integer :: ip, i,j,k, clindx_max

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)

  write (fname_in_MPI, '(a,a,i0)') fname_in, MPI_suffix, coord
    
  inquire (file=fname_in_MPI, exist=exst)
  if (.not. exst) call error (sub_name,                             &
                              'cannot find file ' // fname_in_MPI)

  open (1, file=fname_in_MPI, action='read', position='rewind',  &
         form='unformatted')
  read (1) clindx
  close (1)

  clindx(:, :, nz) = iBOGUS

$else

  inquire (file=fname_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fname_in)

  open (1, file=fname_in, action='read', position='rewind',  &
         form='unformatted')
  read (1) clindx
  close (1)

$endif

!brindx_max = -100
!do k=0,nz
!  do j=1,ny
!    do i=1,ld
!    
!      if(brindx(i,j,k) > brindx_max) brindx_max = brindx(i,j,k)

!    enddo
!  enddo
!enddo

clindx_initialized = .true.

end subroutine clindx_init

!**********************************************************************
subroutine chi_init ()
!**********************************************************************
use param, only : iBOGUS, coord
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.chi_init'
character (*), parameter :: fchi_in = 'chi.out'
$if ($MPI)
  character (*), parameter :: MPI_suffix = '.c'

  character (128) :: fchi_in_MPI
$endif

integer :: ip
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

  open (1, file=fchi_in_MPI, action='read', position='rewind',  &
         form='unformatted')
  read (1) chi
  close (1)

  chi(:, :, nz) = iBOGUS

$else

  inquire (file=fchi_in, exist=exst)
  if (.not. exst) call error (sub_name, 'cannot find file ' // fchi_in)

  open (1, file=fchi_in, action='read', position='rewind',  &
         form='unformatted')
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
implicit none
!  this is used to map the brindx to correct rns tree
allocate( rns_tree_map( ntree ) )

!  Assign rns_tree_map layout
if(rns_tree_layout == 1) then 

  do nt = 1, ntree
  
    if( nt < rns_ntree ) then
        rns_tree_map( nt ) = nt
    else
        rns_tree_map(nt) = rns_ntree
    endif
    
  enddo
  
else
 
  call error(sub_name, 'rns_tree_layout must be 1 for now')
  
endif

return

end subroutine fill_rns_tree_map

!**********************************************************************
subroutine set_ncluster_relsv ()
!**********************************************************************
implicit none
!  Get the number of resolved clusters 
ncl_reslv = 0
do nt = 1, rns_ntree
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      if(ng <= tr_t(nt) % ngen_reslv) ncl_reslv = ncl_reslv + 1   
    enddo    
  enddo
enddo

return

end subroutine set_ncl_reslv

!**********************************************************************
subroutine set_ncluster_tot ()
!**********************************************************************
implicit none
!  Get the total number of clusters 
ncluster_tot = 0
do nt = 1, ntree
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster 
      ncluster_tot = ncluster_tot + 1   
    enddo    
  enddo
enddo
return

end subroutine set_ncluster_tot

!**********************************************************************
subroutine set_nr_elem()
!**********************************************************************
implicit none
 
nr_elem = ncluster_reslv

return
end subroutine set_nr_elem

!**********************************************************************
subroutine set_nb_elem()
!**********************************************************************
use cyl_skew_base_ls, only : tree, generation, tr_t
implicit none

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(tr_t_p, gen_t_p)

nb_elem = 0
do nt = 1, rns_ntree

  !  Point to rns mapped tree
  tr_t_p => tr_t( rns_tree_map( nt ) )
  
  do ng = 1, tr_t_p % ngen_reslv
  
    gen_t_p => tr_t_p % gen_t( ng )
	
    do nc = 1, gen_t_p % ncluster 
	
      nb_elem = nb_elem + 1
    
	enddo    
	
	nullify(gen_t_p)
	
  enddo
  
  nullify(tr_t_p)
  
enddo

return
end subroutine set_nb_elem

!**********************************************************************
subroutine set_nbeta_elem()
!**********************************************************************
use cyl_skew_base_ls, only : tree, generation, tr_t
implicit none

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

!  Nullify all pointers
nullify(tr_t_p, gen_t_p)

nbeta_elem = 0
do nt = 1, rns_ntree

  !  Point to rns mapped tree
  tr_t_p => tr_t( rns_tree_map( nt ) )
  
  do ng = 1, tr_t_p % ngen_reslv + 1
  
    gen_t_p => tr_t_p % gen_t( ng )
	
    do nc = 1, gen_t_p % ncluster 
	
      nbeta_elem = nbeta_elem + 1
    
	enddo    
	
	nullify(gen_t_p)
	
  enddo
  
  nullify(tr_t_p)
  
enddo

return
end subroutine set_nbeta_elem


!**********************************************************************
subroutine fill_beta_elem_to_cl_map()
!**********************************************************************
use cyl_skew_base_ls, only : tree, generation, tr_t, ntree
implicit none

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

integer :: ncount

!  Nullify all pointers
nullify(tr_t_p, gen_t_p)

nbeta = 0
do nt = 1, rns_ntree

  !  Point to rns mapped tree
  tr_t_p => tr_t( rns_tree_map( nt ) )
  
  do ng = 1, tr_t_p % ngen_reslv + 1
  
    gen_t_p => tr_t_p % gen_t( ng )
	
    do nc = 1, gen_t_p % ncluster 
	
      nbeta = nbeta + 1
    
	enddo    
	
	nullify(gen_t_p)
	
  enddo
  
  nullify(tr_t_p)
  
enddo

!  Each resolved cluster recieves a unique id
allocate( beta_elem_to_beta_map( nbeta ) )
r_elem_to_cl_map = -1

ncount=0
if(coord == 0) write(*,*) 'Setting r_elem_to_cl_map'

do nt = 1, rns_ntree

  !  Point to rns mapped tree
  tr_t_p => tr_t( rns_tree_map( nt ) )
  
  do ng = 1, tr_t_p % ngen_reslv
  
    gen_t_p => tr_t_p % gen_t( ng )
	
    do nc = 1, gen_t_p % ncluster 
	
      if(coord == 0) write(*,'(a,3i)') 'rns_tree_map( nt ), ng, nc : ', rns_tree_map( nt ), ng, nc
    
	  ncount = ncount + 1
      
	  r_elem_to_cl_map( ncount ) = gen_t_p % cl_t (nc) %indx
      
	  if(coord == 0) then
        write(*,'(1a,5i4)') 'rns_tree_map( nt ), ng, nc, gen_t_p % cl_t (nc) %indx, rns_reslv_cl_indx : ', rns_tree_map( nt ), ng, nc, gen_t_p % cl_t (nc) %indx, ncount
      endif
    
	enddo    
	
	nullify(gen_t_p)
	
  enddo
  
  nullify(tr_t_p)
  
enddo

return
end subroutine fill_beta_elem_to_cl_map

!**********************************************************************
subroutine fill_cl_to_r_elem_map()
!**********************************************************************
use cyl_skew_base_ls, only : tree, generation, tr_t, ntree
implicit none

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

integer :: ncount

!  Nullify all pointers
nullify(tr_t_p, gen_t_p)

!  Each resolved cluster recieves a unique id
allocate( cl_to_r_elem_map( ncluster_tot ) )
cl_to_r_elem_map = -1

ncount=0
if(coord == 0) write(*,*) 'Setting cl_to_r_elem_map'

do nt = 1, rns_ntree

  !  Point to rns mapped tree
  tr_t_p => tr_t( rns_tree_map( nt ) )
  
  do ng = 1, tr_t_p % ngen_reslv
  
    gen_t_p => tr_t_p % gen_t( ng )
	
    do nc = 1, gen_t_p % ncluster 
	
      if(coord == 0) write(*,'(a,3i)') 'rns_tree_map( nt ), ng, nc : ', rns_tree_map( nt ), ng, nc
    
	  ncount = ncount + 1
      
	  cl_to_r_elem_map( gen_t_p % cl_t (nc) %indx ) = ncount
      
	  if(coord == 0) then
        write(*,'(1a,5i4)') 'rns_tree_map( nt ), ng, nc, gen_t_p % cl_t (nc) %indx, rns_reslv_cl_indx : ', rns_tree_map( nt ), ng, nc, gen_t_p % cl_t (nc) %indx, ncount
      endif
    
	enddo    
	
	nullify(gen_t_p)
	
  enddo
  
  nullify(tr_t_p)
  
enddo

return
end subroutine fill_cl_to_r_elem_map


!**********************************************************************
subroutine fill_cl_to_beta_elem_map
!**********************************************************************
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

integer :: ib, nt, ng, nc
integer :: iparent, ng_parent

integer, pointer :: clindx_p
integer, pointer, dimension(:) :: cl_loc_id_p

type(tree), pointer :: tr_t_p 
type(generation), pointer :: gen_t_p

! Nullify all pointers
nullify(tr_t_p, gen_t_p)

allocate( cl_to_beta_elem_map( ncluster_tot ) ) 
!  
cl_to_beta_elem_map = -1

!  Assign cl_to_beta_elem_map; for a given cluster at gen>=g+1 it returns the
!  the unique rns beta id  
ib=0  
do nt = 1,rns_ntree
    
  tr_t_p => tr_t( rns_tree_map( nt ) ) 
    
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1 )
    
  do nc = 1, gen_t_p % ncluster ! number of clusters in the g+1 generation
    
    ib = ib + 1
      
    clindx_p => gen_t_p % cl_t(nc) % indx
      
    !  So clindx_p corresponds to ib
    cl_to_beta_elem_map(clindx_p) = ib
      
    nullify(clindx_p)
    
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
          
        iparent = tr_t(rns_tree_map( cl_loc_id_p(1)) ) % gen_t (cl_loc_id_p(2)) % &
          cl_t(cl_loc_id_p(3)) % parent
          
        ng_parent = cl_loc_id_p(2) - 1
          
        nullify(cl_loc_id_p)
          
      enddo
        
      cl_to_beta_elem_map( clindx_p ) = cl_to_beta_elem_map( iparent ) ! cl_to_beta_elem_map(iparent) is set from above

      nullify(clindx_p)
        
    enddo
      
    nullify(gen_t_p)
      
  enddo
    
  nullify(tr_t_p)
    
enddo

end subroutine fill_cl_to_beta_elem_map

!**********************************************************************
subroutine fill_cl_to_b_elem_map
!**********************************************************************
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

integer :: ib, nt, ng, nc
integer :: iparent, ng_parent

integer, pointer :: clindx_p => null()

type(tree), pointer :: tr_t_p 
type(generation), pointer :: gen_t_p

! Nullify all pointers
nullify(tr_t_p, gen_t_p)

!  Set the number of beta regions
allocate( cl_to_b_elem_map( ncluster_tot ) )

cl_to_b_elem_map = -1
  
!  Assign rns_rbeta_iarray; for a given cluster at gen=g it returns the
!  the unique rns rbeta id
ib = 0
do nt = 1, rns_ntree
  
  tr_t_p => tr_t( rns_tree_map( nt ) ) 
   
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv )
    
  do nc = 1, gen_t_p % ncluster ! number of clusters in the g generation
    
    ib = ib + 1
      
    clindx_p => gen_t_p % cl_t(nc) % indx
      
    !  So clindx_p corresponds to ib
    cl_to_b_elem_map( clindx_p ) = ib
      
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
          
        iparent = tr_t(rns_tree_map( cl_loc_id_p(1)) ) % gen_t (cl_loc_id_p(2)) % &
          cl_t(cl_loc_id_p(3)) % parent
          
        ng_parent = cl_loc_id_p(2) - 1
          
        nullify(cl_loc_id_p)
          
      enddo
        
      cl_to_b_elem_map( clindx_p ) = cl_to_b_elem_map( iparent ) ! cl_to_beta_elem_map(iparent) is set from above

      nullify(clindx_p)
        
    enddo
      
    nullify(gen_t_p)
      
  enddo
    
  nullify(tr_t_p)
    
enddo

return
end subroutine fill_cl_to_b_elem_map

!**********************************************************************
subroutine fill_r_elem()
!**********************************************************************
implicit none
!! ---- Secondary structures ----
!type ref_region
!  integer :: npoint
!  real(rprec), pointer, dimension(:,:) :: points !  3 ordered points
!  real(rprec) :: u ! reference values
!  real(rprec) :: area
!end type ref_region

!type force
!  integer :: parent !  parent CD; for resolved branches 
!  real(rprec) :: fD
!  real(rprec) :: CD
!  real(rprec) :: kappa ! Used for unresolved branches
!end type force

!type indx_array
!  integer :: npoint
!  integer, pointer, dimension(:,:) :: iarray
!end type

!! ---- Secondary structures ----

!! ---- Primary structures ----
!type primary_struct_type_1
!  type(ref_region)  :: ref_region_t
!  type(force)       :: force_t
!  type(indx_array)  :: indx_array_t
!end type primary_struct_type_1

type(ref_region), pointer, dimension(:) :: pre_r_elem_ref_region_t
type(indx_array), pointer, dimension(:) :: pre_r_elem_indx_array_t

allocate( r_elem_t ( nr_elem ) )

contains 

!======================================================================
subroutine get_ref_region
!======================================================================
use types, only : rprec, vec2d
use param, only : dy, dz, USE_MPI, coord
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_r_elem.get_ref_region'

real(rprec), parameter :: alpha=1._rprec
real(rprec), parameter :: alpha_beta_width = 2.0_rprec
real(rprec), parameter :: alpha_beta_dist = 1.25_rprec

integer :: nt, ng, nc, nb
integer :: ib, irb

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)

integer, pointer :: clindx_p, nbranch_p
integer, pointer :: r_elem_indx_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p
real(rprec), pointer, dimension(:) :: origin_p

type(vec2d) :: rvec_t

type(tree),       pointer :: tr_t_p
type(generation), pointer :: gen_t_p
type(cluster),    pointer :: cl_t_p
type(branch),     pointer :: br_t_p

nullify(d_p, l_p, skew_angle_p)
nullify(clindx_p, nbranch_p)
nullify(r_elem_indx_p)
nullify(tr_t_p, gen_t_p, cl_t_p, br_t_p)


allocate( pre_r_elem_ref_region_t( nr_elem ) )

!  Define entire unresolved region of tree as a beta region
call mesg(sub_name, 'allocating pre_r_elem_ref_region_t')
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling Reference Plane Arrays'
  write(*,*) ' '
endif

do nt=1, rns_ntree

  tr_t_p => tr_t ( rns_tree_map ( nt ) )

  do ng=1, tr_t_p % ngen_reslv
   
    gen_t_p => tr_t_p % gen_t( ng )
  
    do nc = 1, gen_t_p % ncluster
    
      !  Point to cluster
      cl_t_p => gen_t_p % cl_t(nc)
      !  Point to number of branches
      nbranch_p => cl_t_p % nbranch
      !  Point to the cluster index
      clindx_p => cl_t_p % indx
      !  Point to the corresponding r_elem index
      r_elem_indx_p => cl_to_r_elem_map( clindx_p )
      
      if( r_elem_indx_p == -1 ) call error(sub_name, 'Incorrect cluster referenced')
      
      h_m = 0._rprec
      area_proj = 0._rprec
      
      do nb = 1, nbranch_p
        
        br_t_p => cl_t_p % br_t( nb )

        d_p          => br_t_p % d
        l_p          => br_t_p % l
        skew_angle_p => br_t_p % skew_angle
        
        h         = l_p * dcos(skew_angle_p)
        h_m       = h_m + h
        area_proj = area_proj + d_p * h
        
        nullify(d_p, l_p, skew_angle_p, br_t_p)
      
      enddo
      
      !  Mean height of branch cluster  and height of reference area
      h_m = h_m / nbranch_p
      !  width of reference area
      w   = area_proj / h_m     

      pre_r_elem_ref_region_t( r_elem_indx_p ) % area = area_proj
      !  These are defined to be x - planes (no not the NASA experimental planes)

      nzeta = ceiling( w / dy + 1)
      neta  = ceiling( h_m / dz + 1)
      
      pre_r_elem_ref_region_t( r_elem_indx_p ) % npoints = nzeta*neta

      origin_p => cl_t_p % origin
      
      !  Offset in the upstream x-direction
      zeta_c = origin_p + (/ -alpha * w, 0._rprec, 0._rprec /)
      
      !  Need to set the number of points
      !cl_ref_plane_t(clindx_p) % p1    = zeta_c 
      !cl_ref_plane_t(clindx_p) % p1(2) = cl_ref_plane_t(clindx_p) % p1(2) + w / 2._rprec
      
      !cl_ref_plane_t(clindx_p) % p2    = cl_ref_plane_t(clindx_p) % p1
      !cl_ref_plane_t(clindx_p) % p2(2) = cl_ref_plane_t(clindx_p) % p2(2) - w
      
      !cl_ref_plane_t(clindx_p) % p3    = cl_ref_plane_t(clindx_p) % p2
      !cl_ref_plane_t(clindx_p) % p3(3) = cl_ref_plane_t(clindx_p) % p3(3) + h_m
      
      nullify(nbranch_p, clindx_p, origin_p)
      nullify(cl_t_p)
      
    enddo
    
    nullify( gen_t_p )
    
  enddo
  
  nullify(tr_t_p)
 
enddo
    
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  !write(*,*) '--> Reference Plane Values For All Tree Clusters : '
  !do nt=1, rns_ntree

  !  do ng = 1, tr_t(nt)%ngen_reslv
  !    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster

  !      write(*,*) '-------------------------'
  !      write(*,*) 'nt, ng, nc : ', nt, ng, nc  
  !      write(*,*) 'nzeta, neta : ', cl_ref_plane_t(rns_reslv_cl_iarray( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx )) % nzeta, &
  !         cl_ref_plane_t(rns_reslv_cl_iarray ( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % neta
  !       write(*,*) 'p1 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p1
  !       write(*,*) 'p2 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p2
  !       write(*,*) 'p3 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p3
  !       write(*,*) 'area : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % area
  !      
  !       write(*,*) '-------------------------'
  !    enddo
  !  enddo
  !  enddo

endif
          
return

end subroutine get_ref_region

end subroutine fill_r_elem

!**********************************************************************
subroutine rns_fill_ref_plane_array_ls()
!**********************************************************************
!  This subroutine fills all reference plane arrays
!
use types, only : rprec, vec2d
use param, only : dy, dz, USE_MPI, coord
use messages
$if($CYL_SKEW_LS) 
use cyl_skew_base_ls, only : tr_t, tree, generation
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_fill_cl_ref_plane_array_ls'
real(rprec), parameter :: alpha=1._rprec
real(rprec), parameter :: alpha_beta_width = 2.0_rprec
real(rprec), parameter :: alpha_beta_dist = 1.25_rprec

integer :: nt, ng, nc, nb
integer :: ib, irb

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)

integer, pointer :: clindx_p, nbranch_p
integer, pointer :: rbeta_indx_p, beta_indx_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p
real(rprec), pointer, dimension(:) :: origin_p

type(vec2d) :: rvec_t

type(tree), pointer :: tr_t_p
type(generation), pointer :: gen_t_p

nullify(d_p, l_p, skew_angle_p)
nullify(clindx_p, nbranch_p)
nullify(rbeta_indx_p, beta_indx_p)
nullify(tr_t_p, gen_t_p)

allocate( cl_ref_plane_t( ncluster_reslv ) )

if(use_beta_sub_regions) then

  allocate( rbeta_ref_plane_t ( nrbeta ) )
    !  Define entire unresolved region of tree as a beta region
  allocate( beta_ref_plane_t( nbeta ) )

  call mesg(sub_name, 'allocating rbeta, beta ref plane array')
  
else

  call error(sub_name, 'use_beta_sub_regions must be true')
  
endif

!if(ntree_ref < 1) call error( sub_name, 'ntree_ref not specified correctly')

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling Reference Plane Arrays'
  write(*,*) ' '
endif

do nt=1, rns_ntree

  do ng=1, tr_t(nt) % ngen_reslv
  
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
    
      nbranch_p => tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
      
      clindx_p => rns_reslv_cl_iarray( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx )
      
      h_m = 0._rprec
      area_proj = 0._rprec
      
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

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filled Resolved Cluster Reference Plane Array'
  write(*,*) ' '
endif

!  Need to set rbeta_ref_plane_t

do nt = 1, rns_ntree

  tr_t_p => tr_t(rns_tree_iarray(nt))
  
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv )

  do nc = 1, gen_t_p % ncluster
  
    clindx_p => gen_t_p % cl_t(nc) % indx
    
    rbeta_indx_p => rns_rbeta_iarray( clindx_p )
  
    hbot_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv ) % bplane

    htop_p => tr_t_p % gen_t ( tr_t_p % ngen ) % tplane
   
    origin_p => gen_t_p % cl_t(nc) % origin
   
    h = htop_p - hbot_p
  
   !!  Let w = 2 * (2D distance) from the tree origin to the center of a top most cluster
   !rvec_t % xy = tr_t(nt) % gen_t (2) % cl_t(1) % origin(1:2)
   
   !rvec_t % xy = rvec_t % xy - origin_p(1:2)
   
   !call vector_magnitude_2d( rvec_t % xy, rvec_t % mag )
   
    w = alpha_beta_width * h
   
    rbeta_ref_plane_t( rbeta_indx_p ) % area = h * w 
   
    rbeta_ref_plane_t( rbeta_indx_p ) % nzeta = ceiling( w / dy + 1)
    rbeta_ref_plane_t( rbeta_indx_p ) % neta  = ceiling( h / dz + 1)
      
      !  Offset in the upstream x-direction 
    zeta_c = origin_p + (/ -alpha_beta_dist * h, 0._rprec, 0._rprec /)  
    
    rbeta_ref_plane_t(rbeta_indx_p) % p1    = zeta_c 
    rbeta_ref_plane_t(rbeta_indx_p) % p1(2) = rbeta_ref_plane_t(rbeta_indx_p) % p1(2) + w / 2._rprec
      
    rbeta_ref_plane_t(rbeta_indx_p) % p2    = rbeta_ref_plane_t(rbeta_indx_p) % p1
    rbeta_ref_plane_t(rbeta_indx_p) % p2(2) = rbeta_ref_plane_t(rbeta_indx_p) % p2(2) - w
      
    rbeta_ref_plane_t(rbeta_indx_p) % p3    = rbeta_ref_plane_t(rbeta_indx_p) % p2
    rbeta_ref_plane_t(rbeta_indx_p) % p3(3) = rbeta_ref_plane_t(rbeta_indx_p) % p3(3) + h
    
    nullify(hbot_p, htop_p, origin_p, clindx_p, rbeta_indx_p)
    
  enddo
   
  nullify(tr_t_p, gen_t_p)
  
enddo 

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filled RBETA Reference Plane Array'
  write(*,*) ' '
endif

!  Now need to define beta_ref_plane_t; needs to be consistent with corresponding volume of beta_indx_array
if(use_beta_sub_regions) then

do nt = 1, rns_ntree

  tr_t_p => tr_t( rns_tree_iarray(nt) )
  
  if(tr_t_p % ngen_reslv + 1 > tr_t_p % ngen) call error(sub_name, 'helpppp')
  
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1 )

  do nc = 1, gen_t_p % ncluster
  
    clindx_p => gen_t_p % cl_t(nc) % indx
    
    beta_indx_p => rns_beta_iarray( clindx_p )
  
    hbot_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv) % tplane
    
    htop_p => tr_t_p % gen_t ( tr_t_p % ngen ) % tplane
   
    origin_p => gen_t_p % cl_t(nc) % origin
   
    h = htop_p - hbot_p

   !!  Let w = 2 * (2D distance) from the tree origin to the center of a top most cluster
   !rvec_t % xy = tr_t(nt) % gen_t (2) % cl_t(1) % origin(1:2)
   
   !rvec_t % xy = rvec_t % xy - origin_p(1:2)
   
   !call vector_magnitude_2d( rvec_t % xy, rvec_t % mag )
   
    w = alpha_beta_width * h
   
    beta_ref_plane_t( beta_indx_p ) % area = h * w 
   
    beta_ref_plane_t( beta_indx_p ) % nzeta = ceiling( w / dy + 1)
    beta_ref_plane_t( beta_indx_p ) % neta  = ceiling( h / dz + 1)
      
      !  Offset in the upstream x-direction 
    zeta_c = origin_p + (/ -alpha_beta_dist * h, 0._rprec, 0._rprec /)  
    
    beta_ref_plane_t(beta_indx_p) % p1    = zeta_c 
    beta_ref_plane_t(beta_indx_p) % p1(2) = beta_ref_plane_t(beta_indx_p) % p1(2) + w / 2._rprec
      
    beta_ref_plane_t(beta_indx_p) % p2    = beta_ref_plane_t(beta_indx_p) % p1
    beta_ref_plane_t(beta_indx_p) % p2(2) = beta_ref_plane_t(beta_indx_p) % p2(2) - w
      
    beta_ref_plane_t(beta_indx_p) % p3    = beta_ref_plane_t(beta_indx_p) % p2
    beta_ref_plane_t(beta_indx_p) % p3(3) = beta_ref_plane_t(beta_indx_p) % p3(3) + h
    
    nullify(hbot_p, htop_p, origin_p, clindx_p, beta_indx_p)
    
  enddo
   
  nullify(tr_t_p, gen_t_p)
  
enddo

else

  call error(sub_name, 'use_beta_sub_regions must be true for now.')
  
  !do nt = 1, nbeta

  ! hbot_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen_reslv ) % tplane
  ! htop_p => tr_t(nt) % gen_t ( tr_t(nt) % ngen ) % tplane
  ! 
  ! origin_p => tr_t(nt) % origin
  ! 
  ! h = htop_p - hbot_p

  ! !  Let w = 2 * (2D distance) from the tree origin to the center of a top most cluster
  ! rvec_t % xy = tr_t(nt) % gen_t ( 2 ) % cl_t(1) % origin(1:2)
  ! 
  ! rvec_t % xy = rvec_t % xy - origin_p(1:2)
  ! 
  ! call vector_magnitude_2d( rvec_t % xy, rvec_t % mag )
  ! 
  ! w = 2._rprec * rvec_t % mag
  ! 
  ! beta_ref_plane_t(nt) % area = h * w 
  ! 
  ! beta_ref_plane_t( nt ) % nzeta = ceiling( w / dy + 1)
  ! beta_ref_plane_t( nt ) % neta  = ceiling( h / dz + 1)
  !    
  !    !  Offset in the upstream x-direction and displace in z
  !  zeta_c = origin_p + (/ -alpha_beta * rvec_t % mag, 0._rprec, hbot_p - origin_p(3) /)  
  !  
  !  beta_ref_plane_t(nt) % p1    = zeta_c 
  !  beta_ref_plane_t(nt) % p1(2) = beta_ref_plane_t(nt) % p1(2) + w / 2._rprec
  !    
  !  beta_ref_plane_t(nt) % p2    = beta_ref_plane_t(nt) % p1
  !  beta_ref_plane_t(nt) % p2(2) = beta_ref_plane_t(nt) % p2(2) - w
  !    
  !  beta_ref_plane_t(nt) % p3    = beta_ref_plane_t(nt) % p2
  !  beta_ref_plane_t(nt) % p3(3) = beta_ref_plane_t(nt) % p3(3) + h
  !  
  !  nullify(hbot_p, htop_p, origin_p)
  !  
  !enddo
  
endif

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) '--> Filled BETA Reference Plane Array'
  write(*,*) ' '
endif
    
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) '--> Reference Plane Values For All Tree Clusters : '
  do nt=1, rns_ntree

    do ng = 1, tr_t(nt)%ngen_reslv
      do nc = 1, tr_t(nt)%gen_t(ng)%ncluster

        write(*,*) '-------------------------'
        write(*,*) 'nt, ng, nc : ', nt, ng, nc  
        write(*,*) 'nzeta, neta : ', cl_ref_plane_t(rns_reslv_cl_iarray( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx )) % nzeta, &
           cl_ref_plane_t(rns_reslv_cl_iarray ( tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % neta
         write(*,*) 'p1 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p1
         write(*,*) 'p2 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p2
         write(*,*) 'p3 : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % p3
         write(*,*) 'area : ', cl_ref_plane_t(rns_reslv_cl_iarray(tr_t(nt)%gen_t(ng)%cl_t(nc)%indx)) % area
        
         write(*,*) '-------------------------'
      enddo
    enddo
    enddo
    
  write(*,*) ' --> Reference Plane Values For RBETA Regions : '
  do irb=1, nrbeta

    write(*,*) '-------------------------'
    write(*,*) 'irb : ', irb
    write(*,*) 'nzeta, neta : ', rbeta_ref_plane_t(irb) % nzeta, rbeta_ref_plane_t(irb) % neta
    write(*,*) 'p1 : ', rbeta_ref_plane_t(irb) % p1
    write(*,*) 'p2 : ', rbeta_ref_plane_t(irb) % p2
    write(*,*) 'p3 : ', rbeta_ref_plane_t(irb) % p3
    write(*,*) 'area : ', rbeta_ref_plane_t(irb) % area
    write(*,*) '-------------------------'
  enddo
  
  write(*,*) ' --> Reference Plane Values For BETA Regions : '
  do ib=1, nbeta

    write(*,*) '-------------------------'
    write(*,*) 'ib : ', ib
    write(*,*) 'nzeta, neta : ', beta_ref_plane_t(ib) % nzeta, beta_ref_plane_t(ib) % neta
    write(*,*) 'p1 : ', beta_ref_plane_t(ib) % p1
    write(*,*) 'p2 : ', beta_ref_plane_t(ib) % p2
    write(*,*) 'p3 : ', beta_ref_plane_t(ib) % p3
    write(*,*) 'area : ', beta_ref_plane_t(ib) % area
    write(*,*) '-------------------------'
  enddo  

endif
          
      
return

end subroutine rns_fill_ref_plane_array_ls

!**********************************************************************
subroutine get_indx_array_ls()
!**********************************************************************
!  This subroutine gets the indx_array for r_elem and beta_elem. It assigns
!  all information to temporary index structs until all secondary RNS
!  structs can be allocated
!
use types, only : rprec
use param, only : nx,ny,nz, coord, USE_MPI
$if($CYL_SKEW_LS)
use cyl_skew_base_ls, only : ngen, ngen_reslv, brindx_to_loc_id, tr_t
$endif
use level_set_base, only : phi
use messages
implicit none

character (*), parameter :: sub_name = mod_name // '.rns_fill_indx_array_ls'

integer :: i,j,k, nc, np, nb
integer :: ib
integer, pointer :: clindx_p, brindx_p, beta_indx_p
integer, pointer :: nt_p, ng_p, nc_p, npoint_p

integer, pointer, dimension(:) :: cl_loc_id_p

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling Cluster Index Arrays'
  write(*,*) ' '
endif

! ---- Nullify all pointers ----
nullify(clindx_p, cl_loc_id_p)
nullify(nt_p, ng_p, nc_p, npoint_p)
nullify(pre_r_elem_indx_array_t, pre_beta_elem_indx_array_t)

! ---- Allocate all arrays ----
allocate( pre_r_elem_indx_array_t( nr_elem ) )
allocate( pre_beta_elem_indx_array_t( nbeta_elem ) )

do nr=1, nr_elem
  allocate(pre_r_elem_indx_array_t(nc) % iarray(3,nx*ny*(nz-1)))
enddo

do nb=1, nbeta_elem
  allocate(pre_beta_elem_indx_array_t(nb) % iarray(3,nx*ny*(nz-1)))
enddo

!  Intialize the number of points assigned to the cluster
do nr=1, nr_elem
  pre_r_elem_indx_array_t(nc) % npoint = 0
enddo

do nb=1, nbeta_elem
  pre_beta_elem_indx_array_t(nb) % npoint = 0
enddo

if(.not. chi_initialized) call error(sub_name, 'chi not initialized')

do k=1, nz - 1

  do j=1, ny

    do i = 1, nx

      clindx_p => clindx(i,j,k)
           
      if ( clindx_p > 0 ) then
      
	    !  Check if cluster belongs to r_elem
		if( cl_to_r_elem_map( clindx_p ) > 0 ) then
		
          !  Use only inside points
          if ( phi(i,j,k) <= 0._rprec ) then 
		    
		    r_elem_indx_p => cl_to_r_elem_map( clindx_p )

            pre_r_elem_indx_array_t( r_elem_indx_p ) % npoint = pre_r_elem_indx_array_t( r_elem_indx_p ) % npoint + 1
			
			npoint_p => pre_r_elem_indx_array_t( r_elem_indx_p ) % npoint
        
            pre_r_elem_indx_array_t( r_elem_indx_p ) % iarray(1, npoint_p ) = i
            pre_r_elem_indx_array_t( r_elem_indx_p ) % iarray(2, npoint_p ) = j
            pre_r_elem_indx_array_t( r_elem_indx_p ) % iarray(3, npoint_p ) = k
          
			nullify( r_elem_indx_p, npoint_p )
			
          endif
		  
		endif

        !  Check if cluster belongs to beta_elem
        if ( cl_to_beta_elem_map( clindx_p ) > 0 ) then
 
        ! Point needs to be assigned to beta region
          
          if( chi(i,j,k) > chi_cutoff) then
            
		    beta_elem_indx_p => cl_to_beta_elem_map( clindx_p ) 
              
            beta_pre_indx_array_t( beta_elem_indx_p ) % npoint = beta_pre_indx_array_t( beta_elem_indx_p ) % npoint + 1
			
			npoint_p => beta_pre_indx_array_t( beta_elem_indx_p ) % npoint
          
            beta_pre_indx_array_t( beta_elem_indx_p ) % iarray(1, npoint_p ) = i
            beta_pre_indx_array_t( beta_elem_indx_p ) % iarray(2, npoint_p ) = j
            beta_pre_indx_array_t( beta_elem_indx_p ) % iarray(3, npoint_p ) = k 
              
            nullify(beta_elem_indx_p, npoint_p )
              
          endif
                 
        endif
        
      endif
      
      nullify(clindx_p)
      
    enddo
    
  enddo
  
enddo

!  Allocate true indx_array
allocate(cl_indx_array_t( ncl_reslv ) )


do nc=1, ncl_reslv
  allocate(cl_indx_array_t(nc) % iarray(3, cl_pre_indx_array_t(nc) % npoint))
enddo

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
write(*,*) '------------------------'
write(*,*) 'Setting cluster point array'
endif

do nc=1, ncl_reslv

  cl_indx_array_t(nc) % npoint = cl_pre_indx_array_t(nc) % npoint
  
  do np = 1, cl_indx_array_t(nc) % npoint
  
    cl_indx_array_t(nc) % iarray(:,np) = cl_pre_indx_array_t(nc) % iarray(:,np)
    
  enddo
  
enddo

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
write(*,*) '------------------------'
write(*,*) 'Setting BETA point array'
endif
allocate(beta_indx_array_t ( nbeta ) )

do ib = 1, nbeta
  allocate( beta_indx_array_t(ib) % iarray(3, beta_pre_indx_array_t(ib) % npoint))
enddo

do ib = 1, nbeta
  beta_indx_array_t(ib) % npoint = beta_pre_indx_array_t(ib) % npoint
  
  do np = 1, beta_indx_array_t(ib) % npoint
    beta_indx_array_t(ib) % iarray(:,np) = beta_pre_indx_array_t(ib) % iarray(:,np)
  enddo
  
enddo

!  No longer needed
deallocate(cl_pre_indx_array_t)
deallocate(beta_pre_indx_array_t)

!!  Sort each cl_indx_array_t into column major order on the iarray output
!do nc=1, ncluster_tot

!  call isortcm(cl_indx_array_t(nc) % iarray, 3, cl_indx_array_t(nc) % npoint)
!  
!enddo

!if(coord == 0) call mesg(sub_name, 'Exiting ' // sub_name)

return
end subroutine get_indx_array_ls

!!**********************************************************************
!subroutine rns_set_parent_ls()
!!**********************************************************************
!!  This subroutine sets the parent of each beta region
!!
!use types, only : rprec
!$if($CYL_SKEW_LS)
!use cyl_skew_base_ls, only : tree, generation, tr_t
!$endif
!use messages
!implicit none

!character (*), parameter :: sub_name = mod_name // '.rns_set_parent_ls'

!integer :: nc, nt

!integer, pointer :: beta_indx_p
!type(tree), pointer :: tr_t_p
!type(generation) , pointer :: gen_t_p

!nullify(beta_indx_p)
!nullify(tr_t_p, gen_t_p)

!do nt = 1, rns_ntree

!  tr_t_p => tr_t ( rns_tree_iarray(nt) )
!  
!  gen_t_p => tr_t_p % gen_t( tr_t_p % ngen_reslv + 1 )
!  
!  do nc = 1, gen_t_p % ncluster 
!    
!    beta_indx_p => rns_beta_iarray( gen_t_p % cl_t(nc) % indx )
!  
!    beta_force_t(beta_indx_p) % parent = rns_rbeta_iarray( gen_t_p % cl_t(nc) % parent )
!    
!    nullify(beta_indx_p)
!    
!  enddo
!  
!  nullify(tr_t_p, gen_t_p)
!  
!enddo


!return
!end subroutine rns_set_parent_ls





end module rns_cyl_skew_ls
