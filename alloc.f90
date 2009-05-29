!**********************************************************************
subroutine alloc()
!**********************************************************************
!  This subroutine is used to allocate global arrays. It
!  does not include all of the allocatable arrays in the code
!  however.

!$if (LVLSET)
use Sij_defs
!$endif

use convec_defs
use test_filtermodule_defs
use param2,only:nx,ny,nz,ld,lh,ld_big,ny2
!use bottombc, only : zo,patch,num_patch,use_default_patch
use bottombc
use fft, only : kx, ky, k2
use immersedbc, only : alloc_immersedbc

implicit none

character(len=3)::trash

!  Check if running in parallel; allocate lower bounds
!  of z index accordingly (by virtue of fpx3)
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

!  Allocate arrays for patches() in bottombc.f90
allocate(zo(nx,ny))
if(.not. use_default_patch) then
  allocate(patch(nx,ny))
  allocate(T_s(nx,ny))
  allocate(q_s(nx,ny))
  
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

!  Check if using Sij during level set calculations
!$if (LVLSET)
allocate(S11(ld, ny, nz), S12(ld, ny, nz))
allocate(S22(ld, ny, nz), S33(ld, ny, nz))
allocate(S13(ld, ny, nz), S23(ld, ny, nz))
!$endif

!  Allocate arrays for fft module
allocate(kx(lh,ny),ky(lh,ny),k2(lh,ny))

!  Allocate arrays for convec_defs
allocate(cc_big(ld_big,ny2,nz))
allocate(u1_big(ld_big, ny2, $lbz:nz))
allocate(u2_big(ld_big, ny2, $lbz:nz))
allocate(u3_big(ld_big, ny2, $lbz:nz))
allocate(vort1_big(ld_big, ny2, nz))
allocate(vort2_big(ld_big, ny2, nz))
allocate(vort3_big(ld_big, ny2, nz))

!  Allocate arrays in test_filtermodule_defs
allocate(G_test(lh,ny),G_test_test(lh,ny))


!------ Check above for consistency with below ------
!  allocate arrays for immersedbc
call alloc_immersedbc()

!  Allocate arrays for sgsmodule
call alloc_sgsmodule()

return
end subroutine alloc
