!**********************************************************************
module bottombc
!**********************************************************************
use types,only:rprec
use param,only:S_FLAG
use param2,only:nx,ny,ld

implicit none

!private
!public zo_avg,num_patch,zot,zo,phi_m,psi_m,phi_h,psi_h,T_s,q_s,q_mix

logical, parameter :: use_default_patch = .true.

integer::num_patch
integer,allocatable::patchnum(:)
integer, allocatable, dimension(:,:) :: patch

real (rprec), parameter :: zo_default = 0.0001_rprec  !--nondimensional

! num_patch= numbr of patches, zo?=surface roughness for the patches types
! ptypes=number of pacthes types to be used, usually we use 2
real(kind=rprec)::zo_avg,q_mix
real(kind=rprec), allocatable, dimension(:,:) :: zo,T_s,q_s
!TS add for non-neutral case
real(kind=rprec),allocatable, dimension(:,:) :: phi_m,psi_m,phi_h,psi_h
!VK The obukhov similarity functions are computed using obukhov(scalars_module.f90) 
!VK for non-neutral scenario
real(kind=rprec),allocatable::zot(:,:)

end module bottombc

