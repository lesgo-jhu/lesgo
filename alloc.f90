!**********************************************************************
subroutine alloc()
!**********************************************************************
!  This subroutine is used to allocate global arrays. It
!  does not include all of the allocatable arrays in the code
!  however.

use convec_defs
use test_filtermodule_defs
use param2,only:nx,ny,nz,ld,lh,ld_big,ny2,ubc
!use bottombc, only : zo,patch,num_patch,use_default_patch
use bottombc
use fft, only : kx, ky, k2
use immersedbc, only : alloc_immersedbc
use io, only : alloc_io

$if ($LVLSET)
use Sij_defs
use level_set_base, only : alloc_level_set_base
use level_set, only : alloc_level_set
$endif

$if ($TREES_LS)
use trees_global_fmask_ls, only : alloc_trees_global_fmask_ls
use trees_ls, only : alloc_trees_ls
$endif

implicit none

!  Check if running in parallel; allocate lower bounds
!  of z index accordingly (by virtue of fpx3)
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

!  allocate arrays for bottombc
call alloc_bottombc()

!  Allocate arrays for topbc
if(ubc==1)  call alloc_topbc()

!  Allocate arrays for io
call alloc_io()

!  Check if using Sij during level set calculations
$if ($LVLSET)
allocate(S11(ld, ny, nz), S12(ld, ny, nz))
allocate(S22(ld, ny, nz), S33(ld, ny, nz))
allocate(S13(ld, ny, nz), S23(ld, ny, nz))
$endif

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

!  Allocate arrays for scalars_module
call alloc_scalars_module()

!  Allocate arrays for level_set_base and level_set
$if ($LVLSET)
call alloc_level_set_base()
call alloc_level_set()
$endif


$if ($TREES_LS)
!  Allocate arrays for trees_global_fmask_ls
  call alloc_trees_global_fmask_ls ()
!  Allocate arrays for trees_ls
  call alloc_trees_ls()
$endif

!return
end subroutine alloc
