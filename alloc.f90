!**********************************************************************
subroutine alloc()
!**********************************************************************
!  This subroutine is used to allocate global arrays. It
!  does not include all of the allocatable arrays in the code
!  however.

use param2,only:nx,ny,nz
use bottombc, only : zo,patch,num_patch,trash
implicit none
character(len=3)::trash

!  Allocate arrays for patches() in bottombc.f90
allocate(zo(nx,ny))
if(.not. use_default_patch) then
  allocate(patch(nx,ny))
  allocate(T_s(nx,ny))
  allocate(q_s(nxy,ny))
  
!  Load patch data
  open(1,file='patch.dat',position='rewind')
  read(1,*)trash
  read(1,*)trash
  read(1,*)trash
  read(1,*)num_patch  
  
!  Allocate appropriate arrays  
  allocate(zot(num_patch,5))
  allocate(patchnum(num_patch))  
  
endif

return
end subroutine alloc
