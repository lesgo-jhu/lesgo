!///////////////////////////////////////////////////////////////////////////////
module immersedbc
!///////////////////////////////////////////////////////////////////////////////
!
! This module provides definitions used for applying body forces including both
! applied and level set, immersed boundary forces.
!
! In the future this module should probably be removed or renamed more
! appropriately.
!
use types,only:rprec
use param,only:ld,ny,nz
implicit none
private
public fx,fy,fz,fxa,fya,fza

real(kind=rprec),dimension(ld,ny,nz) :: fx = 0._rprec, &
     fy = 0._rprec, &
     fz = 0._rprec
real(kind=rprec),dimension(ld,ny,nz) :: fxa = 0._rprec, &
     fya = 0._rprec, &
     fza = 0._rprec

end module immersedbc
