subroutine dns_stress(txx,txy,txz,tyy,tyz,tzz)
! using the 'sgs' sign convention for stress, so there is a - sign
use types,only:rprec
use param,only:ld,ld_big,nx,ny,nz,z_i,u_star,nu_molec,  &
               coord, nproc, BOGUS
use sim_param,only:dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
implicit none
real(kind=rprec),dimension(ld,ny,nz),intent(out)::txx,txy,txz,tyy, tyz,tzz
real(kind=rprec)::S11,S12,S13,S22,S23,S33
real(kind=rprec)::nu
integer::jx,jy,jz
integer :: jz_min

! non-dimensional molecular viscosity
nu=nu_molec/(z_i*u_star)

! uvp-nodes
do jz=1, nz-1
do jy=1,ny
do jx=1,nx
   S11=dudx(jx,jy,jz)
   S12=0.5_rprec*(dudy(jx,jy,jz)+dvdx(jx,jy,jz))
   S22=dvdy(jx,jy,jz)
   S33=dwdz(jx,jy,jz)
   txx(jx,jy,jz)=-2._rprec*nu*S11
   txy(jx,jy,jz)=-2._rprec*nu*S12
   tyy(jx,jy,jz)=-2._rprec*nu*S22
   tzz(jx,jy,jz)=-2._rprec*nu*S33
end do
end do
end do

!--if values not needed, set to bogus value (easier to catch errors)
$if ($MPI)
if (coord == nproc-1) then
$endif
  ! top values of txx, txy, tyy, tzz not needed for stress free bc's
  !--if values not needed, set to bogus value (easier to catch errors)
  txx(:,:,nz) = BOGUS
  txy(:,:,nz) = BOGUS
  tyy(:,:,nz) = BOGUS
  tzz(:,:,nz) = BOGUS
$if ($MPI)
endif
$endif

! w-nodes
if (coord == 0) then
  ! leave the wall level alone: taken care of with wall stress
  !--assume here that wall stress has already been set (for MPI)
  jz_min = 2
else
  jz_min = 1
end if

do jz = jz_min, nz-1
  do jy = 1, ny
    do jx = 1, nx
      S13 = 0.5_rprec*(dudz(jx,jy,jz)+dwdx(jx,jy,jz))
      S23 = 0.5_rprec*(dvdz(jx,jy,jz)+dwdy(jx,jy,jz))
      txz(jx,jy,jz) = -2._rprec*nu*S13
      tyz(jx,jy,jz) = -2._rprec*nu*S23
    end do
  end do
end do

$if ($MPI) 
  if (coord == nproc-1) then
    ! stress-free lid: this should check which bc options we use
    txz(:,:,nz)=0._rprec
    tyz(:,:,nz)=0._rprec
  else
    !--nz here saves communication in MPI version: can only do this since
    !  dudz, dwdx, dvdz, dwdy are available at nz (not true w/ all components)
    do jy = 1, ny
      do jx = 1, nx
        S13 = 0.5_rprec*(dudz(jx,jy,nz)+dwdx(jx,jy,nz))
        S23 = 0.5_rprec*(dvdz(jx,jy,nz)+dwdy(jx,jy,nz))
        txz(jx,jy,nz) = -2._rprec*nu*S13
        tyz(jx,jy,nz) = -2._rprec*nu*S23
      enddo
    enddo
  endif
$else
  ! stress-free lid: this should check which bc options we use
  txz(:,:,nz)=0._rprec
  tyz(:,:,nz)=0._rprec
$endif

end subroutine dns_stress
