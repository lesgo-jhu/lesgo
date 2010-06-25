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

integer, parameter :: rns_tree_layout = 1

!  Parameters for setting reference regions
real(rprec), parameter :: alpha=1._rprec 
real(rprec), parameter :: alpha_width = 2.0_rprec
real(rprec), parameter :: alpha_dist = 1.25_rprec


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
!  This subroutine allocates and initializes all information needed to 
!  use the resolved elements (r_elem) in the RNS module
!
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_r_elem'
integer :: n

allocate( r_elem_t ( nr_elem ) )

!  Fill reference region
call fill_ref_region()
!  Fill index array
call fill_indx_array()
!  Initialize force data
do n=1, nr_elem
  r_elem_t( n ) % force_t = force( fD=0._rprec, CD=0._rprec, kappa=0._rprec )
enddo

contains 

!======================================================================
subroutine fill_ref_region()
!======================================================================
use param, only : dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_sub_name = mod_name // sub_name // '.get_ref_region'

real(rprec), parameter :: alpha=1._rprec
real(rprec), parameter :: alpha_beta_width = 2.0_rprec
real(rprec), parameter :: alpha_beta_dist = 1.25_rprec

integer :: nt, ng, nc, nb

real(rprec) :: h, h_m, w, area_proj, zeta_c(3)
real(rprec), dimension(3) :: p1, p2, p3

integer, pointer :: r_elem_indx_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p

type(tree),       pointer :: tr_t_p
type(generation), pointer :: gen_t_p
type(cluster),    pointer :: cl_t_p
type(branch),     pointer :: br_t_p

nullify(r_elem_indx_p)
nullify(d_p, l_p, skew_angle_p)
nullify(tr_t_p, gen_t_p, cl_t_p, br_t_p)

call mesg(sub_name, 'allocating pre_r_elem_ref_region_t')
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling R_ELEM Reference Plane Arrays'
  write(*,*) ' '
endif

do nt=1, rns_ntree

  tr_t_p => tr_t ( rns_tree_map ( nt ) )

  do ng=1, tr_t_p % ngen_reslv
   
    gen_t_p => tr_t_p % gen_t( ng )
  
    do nc = 1, gen_t_p % ncluster
    
      !  Point to cluster
      cl_t_p => gen_t_p % cl_t(nc)

      !  Point to the corresponding r_elem index
      r_elem_indx_p => cl_to_r_elem_map( cl_t_p % indx )
      if( r_elem_indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')
      
      !  Initialize mean cluster height and projected area
      h_m = 0._rprec
      area_proj = 0._rprec
      
      do nb = 1, cl_t_p % nbranch
        
        br_t_p => cl_t_p % br_t( nb )

        d_p          => br_t_p % d
        l_p          => br_t_p % l
        skew_angle_p => br_t_p % skew_angle
        
        h         = l_p * dcos(skew_angle_p)
        h_m       = h_m + h
        area_proj = area_proj + d_p * h
        
        nullify( d_p, l_p, skew_angle_p )
        nullify( br_t_p )
      
      enddo
      
      !  Mean height of branch cluster  and height of reference area
      h_m = h_m / nbranch_p
      !  width of reference area
      w   = area_proj / h_m     
      
      nzeta = ceiling( w / dy + 1)
      neta  = ceiling( h_m / dz + 1)
      
      !  Offset in the upstream x-direction
      zeta_c = cl_t_p % origin + (/ -alpha * w, 0._rprec, 0._rprec /)
      
      !  Set the ordered corner points of the plane
      p1    = zeta_c 
      p1(2) = p1(2) + w / 2._rprec
      
      p2    = p1
      p2(2) = p2(2) - w
      
      p3    = p2
      p3(3) = p3(3) + h_m      
      
      !  Point to ref_region_t of the resolved element
      ref_region_t_p => r_elem_t ( r_elem_indx_p ) % ref_region_t

      ref_region_t_p % area    = area_proj
      ref_region_t_p % npoints = nzeta*neta
    
      !  Check if the element has been allocated
      if( allocated( ref_region_t_p % points ) ) then
        call error(sub_name, 'reference region points already allocated.')
      else
        allocate( ref_region_t_p % points( 3, ref_region_t_p % npoints )
      endif
            
      call set_points_in_plane( p1, p2, p3, nzeta, neta, ref_region_t_p % points )    

      !  Finally initialize velocity reference components
      ref_region_t_p % u = 0._rprec
      
      nullify( ref_region_t_p )
      nullify( r_elem_indx_p )
      nullify( cl_t_p )
      
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

end subroutine fill_ref_region

!======================================================================
subroutine fill_indx_array
!======================================================================
!  This subroutine gets the indx_array for r_elem. 
!
use types, only : rprec
use param, only : nx,ny,nz, coord, USE_MPI
use level_set_base, only : phi
use messages
implicit none

character (*), parameter :: sub_sub_name = mod_name // sub_name // '.get_indx_array'

type(indx_array), allocatable, dimension(:) :: pre_indx_array_t

integer :: i,j,k, n, np

integer, pointer :: clindx_p, indx_p, npoint_p

! ---- Nullify all pointers ----
nullify(clindx_p)
nullify(indx_p)
nullify(npoint_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling R_ELEM Index Array'
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

      clindx_p => clindx(i,j,k)
           
      if ( clindx_p > 0 ) then
  
        indx_p => cl_to_r_elem_map( clindx_p )
      
	    !  Check if cluster belongs to r_elem
        if( indx_p > 0 ) then
		
          !  Use only inside points
          if ( phi(i,j,k) <= 0._rprec ) then 


            pre_indx_array_t( indx_p ) % npoint = pre_indx_array_t( indx_p ) % npoint + 1
			
	        npoint_p => pre_indx_array_t( indx_p ) % npoint
        
            pre_indx_array_t( indx_p ) % iarray(1, npoint_p ) = i
            pre_indx_array_t( indx_p ) % iarray(2, npoint_p ) = j
            pre_indx_array_t(_indx_p ) % iarray(3, npoint_p ) = k
          
	        nullify( npoint_p )
			
          endif
		  
        endif

        nullify( indx_p )
       
      endif
      
      nullify(clindx_p)
      
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

return
end subroutine fill_indx_array

end subroutine fill_r_elem

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
call fill_ref_region()
!  Fill index array
call fill_indx_array()
!  Initialize force data
do n=1, nbeta_elem
  beta_elem_t( n ) % force_t = force( fD=0._rprec, CD=0._rprec, kappa=0._rprec )
enddo

return

contains 

!======================================================================
subroutine fill_ref_region()
!======================================================================
use param, only : dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_sub_name = mod_name // sub_name // '.fill_ref_region'

integer :: nt, ng, nc, nb

real(rprec) :: h, h_m, w, area_proj
real(rprec), dimension(3) :: p1, p2, p3, zeta_c

integer, pointer :: indx_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p

type(vec2d) :: rvec_t

type(tree),       pointer :: tr_t_p
type(generation), pointer :: gen_t_p
type(cluster),    pointer :: cl_t_p
type(branch),     pointer :: br_t_p

nullify(indx_p)
nullify(d_p, l_p, skew_angle_p)
nullify(tr_t_p, gen_t_p, cl_t_p, br_t_p)
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling BETA_ELEM Reference Plane Arrays'
  write(*,*) ' '
endif

do nt=1, rns_ntree

  tr_t_p => tr_t ( rns_tree_map ( nt ) )
  
  if(tr_t_p % ngen_reslv + 1 > tr_t_p % ngen) call error(sub_sub_name, 'ngen_reslv + 1 > ngen')
  
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1 )

  do nc = 1, gen_t_p % ncluster
  
    !  Point to cluster
    cl_t_p => gen_t_p % cl_t(nc)
    
    !  Point to the corresponding r_elem index
    indx_p => cl_to_b_elem_map( cl_t_p % indx )    
    if( indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')    
  
    !  Point to the top and bottom of the plane
    hbot_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1) % bplane
    htop_p => tr_t_p % gen_t ( tr_t_p % ngen ) % tplane
   
    origin_p => cl_t_p(nc) % origin
   
    h = htop_p - hbot_p
    w = alpha_width * h
    
    nzeta = ceiling( w / dy + 1)
    neta  = ceiling( h / dz + 1)     
    
    !  Offset in the upstream x-direction 
    zeta_c = cl_t_p % origin + (/ -alpha_dist * h, 0._rprec, 0._rprec /)  

    !  Need to set the number of points
    p1    = zeta_c 
    p1(2) = p1(2) + w / 2._rprec
      
    p2    = p1
    p2(2) = p2(2) - w
      
    p3    = p2
    p3(3) = p3(3) + h

    !  Point to ref_region_t of the resolved element
    ref_region_t_p => beta_elem_t ( indx_p ) % ref_region_t

    ref_region_t_p % area = h * w
    ref_region_t_p % npoints = nzeta*neta
      
    !  Check if the element has been allocated
    if( allocated( ref_region_t_p % points ) ) then
      call error(sub_sub_name, 'reference region points already allocated.')
    else
      allocate( ref_region_t_p % points( 3, ref_region_t_p % npoints )
    endif
            
    call set_points_in_plane( p1, p2, p3, nzeta, neta, ref_region_t_p % points )    

    !  Finally initialize velocity reference components
    ref_region_t_p % u = 0._rprec
      
    nullify( ref_region_t_p )
    nullify( indx_p )
    nullify( cl_t_p )
      
  enddo
    
  nullify( gen_t_p )
  nullify( tr_t_p )
    
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

end subroutine fill_ref_region


!======================================================================
subroutine fill_indx_array
!======================================================================
!  This subroutine gets the indx_array for r_elem. 
!
use types, only : rprec
use param, only : nx,ny,nz, coord, USE_MPI
use level_set_base, only : phi
use messages
implicit none

character (*), parameter :: sub_sub_name = mod_name // sub_name // '.get_indx_array'

type(indx_array), allocatable, dimension(:) :: pre_indx_array_t

integer :: i,j,k, n, np

integer, pointer :: clindx_p, indx_p, npoint_p

! ---- Nullify all pointers ----
nullify(clindx_p)
nullify(indx_p)
nullify(npoint_p)

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling BETA_ELEM Index Array'
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

      clindx_p => clindx(i,j,k)
           
      if ( clindx_p > 0 ) then
  
        indx_p => cl_to_beta_elem_map( clindx_p )
      
	    !  Check if cluster belongs to r_elem
        if( indx_p > 0 ) then
		
          !  Use only points within cutoff
          if ( chi(i,j,k) >= chi_cutoff ) then 


            pre_indx_array_t( indx_p ) % npoint = pre_indx_array_t( indx_p ) % npoint + 1
			
	        npoint_p => pre_indx_array_t( indx_p ) % npoint
        
            pre_indx_array_t( indx_p ) % iarray(1, npoint_p ) = i
            pre_indx_array_t( indx_p ) % iarray(2, npoint_p ) = j
            pre_indx_array_t(_indx_p ) % iarray(3, npoint_p ) = k
          
	        nullify( npoint_p )
			
          endif
		  
        endif

        nullify( indx_p )
       
      endif
      
      nullify(clindx_p)
      
    enddo
    
  enddo
  
enddo

!  Now set set index array for each r_elem
do n=1, nbeta_elem

  beta_elem_t(n) % indx_array_t % npoint = pre_indx_array_t(n) % npoint
  
  allocate(beta_elem_t(n) % indx_array_t % iarray(3, beta_elem_t(n) % indx_array_t % npoint ) )
  
  do np=1, beta_elem_t(n) % indx_array_t % npoint
    beta_elem_t(n) % indx_array_t % iarray(:, np) = pre_indx_array_t( n ) % iarray(:, np)
  enddo
    
enddo

deallocate(pre_indx_array_t)

return
end subroutine fill_indx_array

end subroutine fill_beta_elem

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
call fill_ref_region()
! Set all child elements belonging to b_elem_t
call set_children()
!  Initialize force data
do n=1, nb_elem
  b_elem_t( n ) % force_t = force( fD=0._rprec, CD=0._rprec, kappa=0._rprec )
enddo

return

contains 

!======================================================================
subroutine fill_ref_region()
!======================================================================
use param, only : dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_sub_name = mod_name // sub_name // '.fill_ref_region'

integer :: nt, ng, nc, nb

real(rprec) :: h, h_m, w, area_proj
real(rprec), dimension(3) :: p1, p2, p3, zeta_c

integer, pointer :: indx_p

real(rprec), pointer :: d_p, l_p, skew_angle_p
real(rprec), pointer :: hbot_p, htop_p

type(vec2d) :: rvec_t

type(tree),       pointer :: tr_t_p
type(generation), pointer :: gen_t_p
type(cluster),    pointer :: cl_t_p
type(branch),     pointer :: br_t_p

nullify(indx_p)
nullify(d_p, l_p, skew_angle_p)
nullify(tr_t_p, gen_t_p, cl_t_p, br_t_p)
  
if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Filling B_ELEM Reference Plane Arrays'
  write(*,*) ' '
endif

do nt=1, rns_ntree

  tr_t_p => tr_t ( rns_tree_map ( nt ) )
  
  if(tr_t_p % ngen_reslv > tr_t_p % ngen) call error(sub_sub_name, 'ngen_reslv > ngen')
  
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv )

  do nc = 1, gen_t_p % ncluster
  
    !  Point to cluster
    cl_t_p => gen_t_p % cl_t(nc)
    
    !  Point to the corresponding r_elem index
    indx_p => cl_to_b_elem_map( cl_t_p % indx )    
    if( indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')    
  
    !  Point to the top and bottom of the plane
    hbot_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv ) % bplane
    htop_p => tr_t_p % gen_t ( tr_t_p % ngen ) % tplane
   
    origin_p => cl_t_p(nc) % origin
   
    h = htop_p - hbot_p
    w = alpha_width * h
    
    nzeta = ceiling( w / dy + 1)
    neta  = ceiling( h / dz + 1)     
    
    !  Offset in the upstream x-direction 
    zeta_c = cl_t_p % origin + (/ -alpha_dist * h, 0._rprec, 0._rprec /)  

    !  Need to set the number of points
    p1    = zeta_c 
    p1(2) = p1(2) + w / 2._rprec
      
    p2    = p1
    p2(2) = p2(2) - w
      
    p3    = p2
    p3(3) = p3(3) + h

    !  Point to ref_region_t of the resolved element
    ref_region_t_p => b_elem_t ( indx_p ) % ref_region_t

    ref_region_t_p % area = h * w
    ref_region_t_p % npoints = nzeta*neta
      
    !  Check if the element has been allocated
    if( allocated( ref_region_t_p % points ) ) then
      call error(sub_sub_name, 'reference region points already allocated.')
    else
      allocate( ref_region_t_p % points( 3, ref_region_t_p % npoints )
    endif
            
    call set_points_in_plane( p1, p2, p3, nzeta, neta, ref_region_t_p % points )    

    !  Finally initialize velocity reference components
    ref_region_t_p % u = 0._rprec
      
    nullify( ref_region_t_p )
    nullify( indx_p )
    nullify( cl_t_p )
      
  enddo
    
  nullify( gen_t_p )
  nullify( tr_t_p )
    
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

end subroutine fill_ref_region

!======================================================================
subroutine set_children()
!======================================================================
!  This subroutine sets the r_elem and beta_elem (children) which belong to 
!  each B region
use param, only : dy, dz, USE_MPI, coord
use param, only : nx, ny, nz
use messages
use cyl_skew_base_ls, only : tr_t, tree, generation
implicit none

character (*), parameter :: sub_sub_name = mod_name // sub_name // '.set_children'

integer :: nt, ng, nc, nb

integer, pointer :: r_elem_indx_p
integer, pointer :: beta_elem_indx_p
integer, pointer :: b_elem_indx_p

type(tree),       pointer :: tr_t_p
type(generation), pointer :: gen_t_p
type(cluster),    pointer :: cl_t_p

type(child_elem), allocatable, dimension(:) :: pre_r_child_t
type(child_elem), allocatable, dimension(:) :: pre_beta_child_t

nullify(r_elem_indx_p, beta_elem_indx_p, b_elem_indx_p)
nullify(tr_t_p, gen_t_p, cl_t_p)

allocate(pre_r_child_t( nb_elem ))
allocate(pre_beta_child_t( nb_elem ))

if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
  write(*,*) ' '
  write(*,*) 'Setting B_ELEM Children'
  write(*,*) ' '
endif

!  Initialize temporary children
do n=1, nb_elem

  pre_r_child_t(n) % nelem = 0
  allocate(pre_r_child_t(n) % indx( nr_elem ) )
  pre_r_child_t(n) % indx(:) = -1

enddo

do nt=1, rns_ntree

  tr_t_p => tr_t ( rns_tree_map ( nt ) )
  
  if(tr_t_p % ngen_reslv > tr_t_p % ngen) call error(sub_sub_name, 'ngen_reslv > ngen')
  
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv )

  do nc = 1, gen_t_p % ncluster
  
    !  Point to cluster
    cl_t_p => gen_t_p % cl_t(nc)
    
    !  Point to the corresponding r_elem index
    r_elem_indx_p => cl_to_r_elem_map( cl_t_p % indx )   
    if( r_elem_indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')  
	
    b_elem_indx_p => cl_to_beta_elem_map(  cl_t_p % indx ) 
    if( b_elem_indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')  
	
	pre_r_child_t( b_elem_indx_p ) % nelem = pre_r_child_t( b_elem_indx_p ) % nelem + 1
	pre_r_child_t( b_elem_indx_p ) % indx( pre_r_child_t( b_elem_indx_p ) % nelem ) = r_elem_indx_p
	
	nullify( b_elem_indx_p )
	nullify( r_elem_indx_p ) 
	
	nullify( cl_t_p )
	
  enddo
  
  nullify( gen_t_p )
  nullify( tr_t_p )
  
enddo

!  Now allocate and set the actual children
do n=1, nb_elem
  allocate( b_elem_t( nb ) % r_child_t % indx( pre_r_child_t( n ) % nelem ))
  b_elem_t( nb ) % r_child_t % nelem = pre_r_child_t( n ) % nelem
  b_elem_t( nb ) % r_child_t % indx(:) = pre_r_child_t( n ) % indx(:)
enddo

deallocate( pre_r_child_t )

!  Initialize temporary children
do n=1, nb_elem
  
  pre_beta_child_t(n) % nelem = 0
  allocate(pre_beta_child_t(n) % indx( nbeta_elem ) )
  pre_beta_child_t(n) % indx( : ) = -1

enddo

do nt=1, rns_ntree

  tr_t_p => tr_t ( rns_tree_map ( nt ) )
  
  if(tr_t_p % ngen_reslv + 1 > tr_t_p % ngen) call error(sub_sub_name, 'ngen_reslv + 1 > ngen')
  
  gen_t_p => tr_t_p % gen_t ( tr_t_p % ngen_reslv + 1)

  do nc = 1, gen_t_p % ncluster
  
    !  Point to cluster
    cl_t_p => gen_t_p % cl_t(nc)
    
    !  Point to the corresponding r_elem index
    beta_elem_indx_p => cl_to_beta_elem_map( cl_t_p % indx )   
    if( beta_elem_indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')  
	
    b_elem_indx_p => cl_to_beta_elem_map( cl_t_p % indx ) 
    if( b_elem_indx_p == -1 ) call error(sub_sub_name, 'Incorrect cluster referenced')  
	
	pre_beta_child_t( b_elem_indx_p ) % nelem = pre_beta_child_t( b_elem_indx_p ) % nelem + 1
	pre_beta_child_t( b_elem_indx_p ) % indx( pre_beta_child_t( b_elem_indx_p ) % nelem ) = beta_elem_indx_p
	
	nullify( b_elem_indx_p )
	nullify( beta_elem_indx_p ) 
	
	nullify( cl_t_p )
	
  enddo
  
  nullify( gen_t_p )
  nullify( tr_t_p )
  
enddo

!  Now allocate and set the actual children
do n=1, nb_elem
  allocate( b_elem_t( n ) % beta_child_t % indx( pre_beta_child_t( n ) % nelem ))
  b_elem_t( n ) % beta_child_t % nelem = pre_beta_child_t( n ) % nelem
  b_elem_t( n ) % beta_child_t % indx(:) = pre_beta_child_t( n ) % indx(:)
enddo

deallocate( pre_beta_child_t )

return
end subroutine set_children

end subroutine fill_b_elem

!!**********************************************************************
!subroutine get_indx_array_ls()
!!**********************************************************************
!!  This subroutine gets the indx_array for r_elem and beta_elem. It assigns
!!  all information to temporary index structs until all secondary RNS
!!  structs can be allocated
!!
!use types, only : rprec
!use param, only : nx,ny,nz, coord, USE_MPI
!$if($CYL_SKEW_LS)
!use cyl_skew_base_ls, only : ngen, ngen_reslv, brindx_to_loc_id, tr_t
!$endif
!use level_set_base, only : phi
!use messages
!implicit none

!character (*), parameter :: sub_name = mod_name // '.rns_fill_indx_array_ls'

!integer :: i,j,k, nc, np, nb
!integer :: ib
!integer, pointer :: clindx_p, brindx_p, beta_indx_p
!integer, pointer :: nt_p, ng_p, nc_p, npoint_p

!integer, pointer, dimension(:) :: cl_loc_id_p

!if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
!  write(*,*) ' '
!  write(*,*) 'Filling Cluster Index Arrays'
!  write(*,*) ' '
!endif

!! ---- Nullify all pointers ----
!nullify(clindx_p, cl_loc_id_p)
!nullify(nt_p, ng_p, nc_p, npoint_p)
!nullify(pre_r_elem_indx_array_t, pre_beta_elem_indx_array_t)

!! ---- Allocate all arrays ----
!allocate( pre_r_elem_indx_array_t( nr_elem ) )
!allocate( pre_beta_elem_indx_array_t( nbeta_elem ) )

!do nr=1, nr_elem
!  allocate(pre_r_elem_indx_array_t(nc) % iarray(3,nx*ny*(nz-1)))
!enddo

!do nb=1, nbeta_elem
!  allocate(pre_beta_elem_indx_array_t(nb) % iarray(3,nx*ny*(nz-1)))
!enddo

!!  Intialize the number of points assigned to the cluster
!do nr=1, nr_elem
!  pre_r_elem_indx_array_t(nc) % npoint = 0
!enddo

!do nb=1, nbeta_elem
!  pre_beta_elem_indx_array_t(nb) % npoint = 0
!enddo

!if(.not. chi_initialized) call error(sub_name, 'chi not initialized')

!do k=1, nz - 1

!  do j=1, ny

!    do i = 1, nx

!      clindx_p => clindx(i,j,k)
!           
!      if ( clindx_p > 0 ) then
!      
!	    !  Check if cluster belongs to r_elem
!		if( cl_to_r_elem_map( clindx_p ) > 0 ) then
!		
!          !  Use only inside points
!          if ( phi(i,j,k) <= 0._rprec ) then 
!		    
!		    r_elem_indx_p => cl_to_r_elem_map( clindx_p )

!            pre_r_elem_indx_array_t( r_elem_indx_p ) % npoint = pre_r_elem_indx_array_t( r_elem_indx_p ) % npoint + 1
!			
!			npoint_p => pre_r_elem_indx_array_t( r_elem_indx_p ) % npoint
!        
!            pre_r_elem_indx_array_t( r_elem_indx_p ) % iarray(1, npoint_p ) = i
!            pre_r_elem_indx_array_t( r_elem_indx_p ) % iarray(2, npoint_p ) = j
!            pre_r_elem_indx_array_t( r_elem_indx_p ) % iarray(3, npoint_p ) = k
!          
!			nullify( r_elem_indx_p, npoint_p )
!			
!          endif
!		  
!		endif

!        !  Check if cluster belongs to beta_elem
!        if ( cl_to_beta_elem_map( clindx_p ) > 0 ) then
! 
!        ! Point needs to be assigned to beta region
!          
!          if( chi(i,j,k) > chi_cutoff) then
!            
!		    beta_elem_indx_p => cl_to_beta_elem_map( clindx_p ) 
!              
!            beta_pre_indx_array_t( beta_elem_indx_p ) % npoint = beta_pre_indx_array_t( beta_elem_indx_p ) % npoint + 1
!			
!			npoint_p => beta_pre_indx_array_t( beta_elem_indx_p ) % npoint
!          
!            beta_pre_indx_array_t( beta_elem_indx_p ) % iarray(1, npoint_p ) = i
!            beta_pre_indx_array_t( beta_elem_indx_p ) % iarray(2, npoint_p ) = j
!            beta_pre_indx_array_t( beta_elem_indx_p ) % iarray(3, npoint_p ) = k 
!              
!            nullify(beta_elem_indx_p, npoint_p )
!              
!          endif
!                 
!        endif
!        
!      endif
!      
!      nullify(clindx_p)
!      
!    enddo
!    
!  enddo
!  
!enddo

!!  Allocate true indx_array
!allocate(cl_indx_array_t( ncl_reslv ) )


!do nc=1, ncl_reslv
!  allocate(cl_indx_array_t(nc) % iarray(3, cl_pre_indx_array_t(nc) % npoint))
!enddo

!if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
!write(*,*) '------------------------'
!write(*,*) 'Setting cluster point array'
!endif

!do nc=1, ncl_reslv

!  cl_indx_array_t(nc) % npoint = cl_pre_indx_array_t(nc) % npoint
!  
!  do np = 1, cl_indx_array_t(nc) % npoint
!  
!    cl_indx_array_t(nc) % iarray(:,np) = cl_pre_indx_array_t(nc) % iarray(:,np)
!    
!  enddo
!  
!enddo

!if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
!write(*,*) '------------------------'
!write(*,*) 'Setting BETA point array'
!endif
!allocate(beta_indx_array_t ( nbeta ) )

!do ib = 1, nbeta
!  allocate( beta_indx_array_t(ib) % iarray(3, beta_pre_indx_array_t(ib) % npoint))
!enddo

!do ib = 1, nbeta
!  beta_indx_array_t(ib) % npoint = beta_pre_indx_array_t(ib) % npoint
!  
!  do np = 1, beta_indx_array_t(ib) % npoint
!    beta_indx_array_t(ib) % iarray(:,np) = beta_pre_indx_array_t(ib) % iarray(:,np)
!  enddo
!  
!enddo

!!  No longer needed
!deallocate(cl_pre_indx_array_t)
!deallocate(beta_pre_indx_array_t)

!!!  Sort each cl_indx_array_t into column major order on the iarray output
!!do nc=1, ncluster_tot

!!  call isortcm(cl_indx_array_t(nc) % iarray, 3, cl_indx_array_t(nc) % npoint)
!!  
!!enddo

!!if(coord == 0) call mesg(sub_name, 'Exiting ' // sub_name)

!return
!end subroutine get_indx_array_ls

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

!**********************************************************************
subroutine get_points_in_plane(bp1, bp2, bp3, nzeta, neta, points)
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

character (*), parameter :: func_name = mod_name // '.get_points_in_plane'

integer :: i, j, np

REAL(RPREC) :: dzeta, deta, Lzeta, Leta, vec_mag, zmin, zmax

real(RPREC), dimension(3) :: zeta_vec, eta_vec, eta
real(RPREC), dimension(3) :: bp4, cell_center

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

    endif

  enddo
enddo

return

end function get_points_in_plane



end module rns_cyl_skew_ls
