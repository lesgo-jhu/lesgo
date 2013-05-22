!!
!!  Copyright 2009,2010,2011,2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

!**********************************************************************
module cyl_skew_ls
!**********************************************************************
use types, only : rprec
use param
use cyl_skew_base_ls
implicit none

save
private

public :: fill_tree_array_ls

character (*), parameter :: mod_name = 'cyl_skew_base_ls'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine fill_tree_array_ls()
!**********************************************************************
!
!  This subroutine sets all values for the tree struct - tr_t; it also
!  defines the key arrays clindex_to_loc_id and brindx_to_loc_id 
!  which are used to get the local id of a global index. This subroutine
!  should be called any time the tree struct and its settings as defined
!  in cyl_skew_base_ls are needed (only once though)
!
use param, only : coord
implicit none

integer :: nt, ng, nc, nb, nc_g1, nb_cl
integer :: clindx, brindx
integer :: nbcount, nbcount_tot, nccount, nccount_tot
integer :: nreslv_ccount_tot, reslv_clindx
integer :: nunreslv_ccount_tot, unreslv_clindx

!integer, pointer, dimension(:) :: cl_loc_id

real(rprec) :: angle
real(rprec) :: gen_scale_fact

integer, pointer :: pclindx_p
type(cluster), pointer :: cl_t_p, pcl_t_p, ccl_t_p

!  Nullify all pointers
nullify(pclindx_p, cl_t_p)
nullify(pcl_t_p, ccl_t_p)

allocate(tr_t(ntree))

!  Set the number of generations in the tree; for now 
!  they are all the same
do nt=1, ntree
  tr_t(nt) % ngen = ngen
  tr_t(nt) % ngen_reslv = ngen_reslv
enddo

!  Allocate the number of clusters in the generation
do nt=1, ntree
  
  allocate(tr_t(nt)%gen_t( tr_t(nt)%ngen ))

  do ng=1, tr_t(nt)%ngen
    
    !  Set the number of clusters for the generation
    tr_t(nt)%gen_t(ng)%ncluster = nbranch**(ng - 1)
    
    !  Allocate space for the clusters    
    allocate( tr_t(nt)%gen_t(ng)%cl_t( tr_t(nt)%gen_t(ng)%ncluster ))
    
    do nc=1, tr_t(nt)%gen_t(ng)%ncluster
    
      !  Set the number of branches for the cluster - here they are all the same
      tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch = nbranch
      !  Allocate space for the branches
      allocate( tr_t(nt)%gen_t(ng)%cl_t(nc)%br_t( tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch ))
     
    enddo
    
  enddo
 
enddo

!if(coord == 0) then
!do nt=1, ntree
!  write(*,*) 'tr_t(nt)%ngen : ', tr_t(nt)%ngen
!  do ng=1, tr_t(nt)%ngen
!    write(*,*) 'tr_t(nt)%gen_t(ng)%ncluster : ', tr_t(nt)%gen_t(ng)%ncluster
!    do nc=1, tr_t(nt)%gen_t(ng)%ncluster
!      write(*,*) 'tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch : ', tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
!    enddo
!  enddo
!enddo
!endif
!  ----- Start setting values for the tree struc -----

!  Must allocate all arrays before setting data; may be that
!  memory addresses move as needed during the allocation of 
!  sub types
do nt=1, ntree
  tr_t(nt) % origin = tree_location(nt) % xyz
enddo

do nt = 1, ntree
    
  do ng = 1, tr_t(nt)%ngen
    
    gen_scale_fact = scale_fact**(ng-1)
          
    do nc = 1, tr_t(nt)%gen_t(ng)%ncluster

      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % offset = gen_scale_fact*offset
           
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % d = gen_scale_fact * d
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % l = gen_scale_fact * l
            
      ! Ellipse minor axis
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % b = &
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % d / 2._rprec 
            
      ! Ellipse major axis
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % a = &
      tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t % b/cos(skew_angle) 
            
      do nb = 1, tr_t(nt) % gen_t(ng) % cl_t(nc) % nbranch

        angle =  zrot_angle + &
          2.*pi*(nb-1)/(tr_t(nt) % gen_t(ng) % cl_t(nc) % nbranch) + &
          (ng - 1)*pi ! Rotate 180 degrees for each generation
                    
        tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % angle = angle
                   
                
        tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % skew_axis = &
          (/ cos(angle +pi/2.), sin(angle + pi/2.), 0._rprec/)
                        
        tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % skew_angle = skew_angle         
                  
                 
                 
      enddo
        
    enddo
        
  enddo
    
enddo

!  Get the number of resolved clusters on all trees
ncluster_reslv = 0
do nt = 1, ntree
  do ng = 1, tr_t( nt ) % ngen
    do nc = 1, tr_t( nt ) % gen_t(ng) % ncluster 
      if(ng <= tr_t( nt ) % ngen_reslv) ncluster_reslv = ncluster_reslv + 1   
    enddo    
  enddo
enddo

!  Get the total number of clusters on all trees
ncluster_tot = 0
do nt = 1, ntree
  do ng = 1, tr_t( nt ) % ngen
    do nc = 1, tr_t( nt ) % gen_t(ng) % ncluster 
      ncluster_tot = ncluster_tot + 1   
    enddo    
  enddo
enddo

!  Get a global count of clusters and branches (all may not be the same)
nccount_tot = 0
nbcount_tot = 0
nreslv_ccount_tot = 0
nunreslv_ccount_tot = 0

do nt = 1, ntree
  nccount = 0
  nbcount = 0
  do ng = 1, tr_t(nt) % ngen
    do nc = 1, tr_t(nt) % gen_t(ng) % ncluster
    
      if(ng <= tr_t(nt) % ngen_reslv) then
        nreslv_ccount_tot = nreslv_ccount_tot + 1
      elseif(ng <= tr_t(nt) % ngen) then
        nunreslv_ccount_tot = nunreslv_ccount_tot + 1
      endif      
      
      nccount_tot = nccount_tot + 1
      nccount = nccount + 1

      do nb = 1, tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
        nbcount_tot = nbcount_tot + 1
        nbcount = nbcount + 1
      enddo
    enddo
  enddo
  !  Set the total number of clusters and branches of the trees
  tr_t(nt) % ncluster = nccount
  tr_t(nt) % nbranch  = nbcount

enddo

allocate(clindx_to_loc_id(3,nccount_tot))
allocate(brindx_to_loc_id(4,nbcount_tot))
allocate(reslv_clindx_to_loc_id(3,nreslv_ccount_tot))
allocate(unreslv_clindx_to_loc_id(3,nunreslv_ccount_tot))


!  Initialize global indexes for clusters and branches
clindx = 0
brindx = 0
reslv_clindx = 0
unreslv_clindx = 0

if(coord == 0) write(*,*) 'Filling Tree Array'

do nt = 1, ntree

    do ng=1, tr_t(nt)%ngen
 
        !  Set cluster id for ng+1 generation
        nc_g1 = 0
        
        do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
        
            !  Update global cluster index
            clindx = clindx + 1  
            tr_t(nt)%gen_t(ng)%cl_t(nc)%indx = clindx
            clindx_to_loc_id(:,clindx) = (/ nt, ng, nc /)
            
            if(ng <= tr_t(nt) % ngen_reslv) then
              reslv_clindx = reslv_clindx + 1
              reslv_clindx_to_loc_id(:, reslv_clindx) = (/ nt, ng, nc /)
            elseif(ng <= tr_t(nt) % ngen) then
              unreslv_clindx = unreslv_clindx + 1
              unreslv_clindx_to_loc_id(:, unreslv_clindx) = (/ nt, ng, nc /)
            endif
                       
            !  Set cluster origin to tree origin
            if( ng == 1 ) tr_t(nt) % gen_t(ng) % cl_t(nc) % origin = tr_t(nt) % origin
            
            do nb = 1, tr_t(nt)%gen_t(ng)%cl_t(nc)%nbranch
                !  Update global branch index
                brindx = brindx + 1
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % indx = brindx
                brindx_to_loc_id(:,brindx) = (/ nt, ng, nc, nb /)
                
                !  Update cluster id for ng+1 generation
                nc_g1 = nc_g1 + 1
                    
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % origin
                
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot(1) = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % origin(1) + &          
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % offset * &
                    cos(tr_t(nt) % gen_t(ng) %cl_t(nc) % br_t(nb) % angle)
                
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot(2) = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % origin(2) + &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % offset * &
                    sin(tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % angle)
                   
                call rotation_axis_vector_3d( tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % skew_axis, &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % skew_angle, &
                    (/ 0._rprec, 0._rprec, tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % l /), &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % top )
                
                tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % top = & 
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % top + &
                    tr_t(nt)%gen_t(ng)%cl_t(nc)% br_t(nb) % bot
                
            !  Now set the cluster origin of the ng+1 cluster (with nc = nb)
                if ( ng < tr_t(nt)%ngen ) then
                    tr_t(nt) % gen_t(ng+1) % cl_t(nc_g1) % origin = &
                    tr_t(nt) % gen_t(ng) % cl_t(nc) % br_t(nb) % top
                endif

 
            enddo
                   
        enddo
        
        !  Set the top and bottom plane of the generation - this assumes that all 
        !  branches are the same height!      
        tr_t(nt) % gen_t(ng) % bplane = tr_t(nt) % gen_t(ng) % cl_t(1) % br_t(1) % bot(3)
        tr_t(nt) % gen_t(ng) % tplane = tr_t(nt) % gen_t(ng) % cl_t(1) % br_t(1) % top(3)
    
    enddo
    
enddo

!do nb = 1, brindx
!  write(*,'(a,5i)') 'coord, brindx_to_loc_id(:,nb) : ', coord, brindx_to_loc_id(:,nb)
!enddo

!  Now compute the parent cluster of each cluster
!  setting all generation 1 to have no parent
do nt = 1, ntree

  do nc = 1, tr_t(nt) % gen_t(1) % ncluster
  
    tr_t(nt) % gen_t(1) % cl_t(nc) % parent = 0
    
  enddo
  
enddo

!  Loop over all trees and set parent for generations > 1
if(ngen > 1) then

  
  do nt = 1, ntree

    if(coord == 0) write(*,*) 'tr_t(nt) % ncluster : ', tr_t(nt) % ncluster

    do ng=1, tr_t(nt)%ngen - 1
    
      nb_cl = 0 ! number of branches in g = number of clusters in g + 1
      
      do nc = 1, tr_t(nt)%gen_t(ng)%ncluster
      
        !  Parent cluster
        pcl_t_p => tr_t(nt) % gen_t(ng) % cl_t(nc)
      
        do nb = 1, pcl_t_p % nbranch
        
          nb_cl = nb_cl + 1
      
          !  Child cluster of the parent
          ccl_t_p => tr_t(nt) % gen_t(ng+1) % cl_t(nb_cl)
        
          ccl_t_p % parent = pcl_t_p % indx
        
          nullify(ccl_t_p)
        
        enddo
      
        nullify(pcl_t_p)
      
      enddo
    
    enddo

  enddo
  
endif
        
return
end subroutine fill_tree_array_ls

end module cyl_skew_ls
