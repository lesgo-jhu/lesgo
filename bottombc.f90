module bottombc
use types,only:rprec
use param,only:nx,ny,ld,S_FLAG,use_default_patch,zo_default
implicit none
private
public zo_avg,num_patch,zot,zo,phi_m,psi_m,phi_h,psi_h,&
     patches,avgpatch,T_s,q_s,q_mix

! num_patch= numbr of patches, zo?=surface roughness for the patches types
! ptypes=number of pacthes types to be used, usually we use 2
real(kind=rprec)::zo_avg,q_mix
real(kind=rprec),dimension(nx,ny)::zo,T_s,q_s
!TS add for non-neutral case
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m,phi_h,psi_h
!VK The obukhov similarity functions are computed using obukhov(scalars_module.f90) 
!VK for non-neutral scenario
integer,dimension(nx,ny)::patch
integer::num_patch
integer,allocatable::patchnum(:)
real(kind=rprec),allocatable::zot(:,:)

contains

subroutine patches()
!VK This assigns momentum roughness, temperature and wetness 
!VK for the different patches
!VK and fills the lookup tables patch and patchnum
!VK This is called from routines patch_or_remote.f90 depending
!VK whether to use remotely-sensed data or patches
!TS NOTES ZOT(:,2:4) should locate in the same places for zo, T_s, q_s
!TS OTHERWISE INCONSISTENCY
!use types,only:rprec
use param,only:z_i,t_scale
implicit none
integer::i,j,jp
character(len=3)::trash

!--inserted quick option for default: single patch with specified roughness
if (use_default_patch) then
  num_patch = 1
  zo = zo_default
  zo_avg = zo_default
else
  open(1,file='patch.dat')
  read(1,*)trash
  read(1,*)trash
  read(1,*)trash
  read(1,*)num_patch
  allocate(zot(num_patch,5))
  allocate(patchnum(num_patch))
  patchnum=0
  zo_avg=0._rprec
  do i=1,num_patch
     read(1,*)zot(i,1:5)
  enddo
  do jp=1,num_patch
  do j=int(zot(jp,4)),int(zot(jp,5))
  do i=int(zot(jp,2)),int(zot(jp,3))
     zo(i,j)=zot(jp,1)/z_i
     patch(i,j)=jp
     patchnum(jp)=patchnum(jp)+1
  enddo
  enddo
  zo_avg=zo_avg+zot(jp,1)
  enddo
  zo_avg=zo_avg/real(num_patch,kind=rprec)
  patchnum(1)=patchnum(1)-sum(patchnum(2:num_patch))

  if(S_FLAG)then
   read(1,*)trash
   do i=1,num_patch
      read(1,*)zot(i,1:5)
   enddo
   do jp=1,num_patch
   do i=int(zot(jp,2)),int(zot(jp,3))
   do j=int(zot(jp,4)),int(zot(jp,5))
      T_s(i,j)=zot(jp,1)/t_scale
   enddo
   enddo
   enddo
   read(1,*)trash
   do i=1,num_patch
      read(1,*)zot(i,1:5)
   enddo
   q_mix=0._rprec
   do jp=1,num_patch
   do i=int(zot(jp,2)),int(zot(jp,3))
   do j=int(zot(jp,4)),int(zot(jp,5))
      q_s(i,j)=zot(jp,1)
   enddo
   enddo
   q_mix=q_mix+zot(jp,1)
   enddo
   q_mix=q_mix/real(num_patch,kind=rprec)
  endif
  close(1)
end if

end subroutine patches

subroutine avgpatch(u,u_avg)
! computes the averaged value of a variable (at the wall) over a patch
! and assigns it to an nx X ny array

!use types,only:rprec
implicit none
real(kind=rprec),dimension(nx,ny),intent(in)::u
real(kind=rprec),dimension(nx,ny),intent(out)::u_avg
integer::jp,is,ie,js,je
real(kind=rprec)::temp

!TS calculate the average for the background values
!TS NEED SPECIAL TREATMENT
jp=1;is=int(zot(jp,2));ie=int(zot(jp,3));js=int(zot(jp,4));je=int(zot(jp,5))
temp=sum(u(is:ie,js:je))
do jp=2,num_patch
is=int(zot(jp,2));ie=int(zot(jp,3));js=int(zot(jp,4));je=int(zot(jp,5))
temp=temp-sum(u(is:ie,js:je))
enddo
jp=1;is=int(zot(jp,2));ie=int(zot(jp,3));js=int(zot(jp,4));je=int(zot(jp,5))
u_avg(is:ie,js:je)=temp/real(patchnum(1),kind=rprec)
do jp=2,num_patch
is=int(zot(jp,2));ie=int(zot(jp,3));js=int(zot(jp,4));je=int(zot(jp,5))
u_avg(is:ie,js:je)=sum(u(is:ie,js:je))/real(patchnum(jp),kind=rprec)
enddo
end subroutine avgpatch
end module bottombc
