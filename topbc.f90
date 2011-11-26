module topbc
use types,only:rprec
use param,only:nz
implicit none
real(kind=rprec),dimension(nz)::sponge

contains
subroutine setsponge()
use param
implicit none
real(kind=rprec)::factor
integer::k
$if ($MPI)
  integer, parameter :: nz_global = nz * nproc
  integer :: k_global
  real (rprec) :: sponge_top
$endif
! sets relaxation term to vertical momentum equation in top quarter
! of domain, relaxation timescale on top is 50s with a factor of 5 if we
! had Nz=40 for the layers 40...31, Nieuwstadt et al. 1991, turbulent shear 
! flows

$if ($MPI)
  !--the non-MPI recursive form in inconvenient, so instead replace with
  !  analytically evaluated version (should be this way anyway...)
  sponge=0._rprec
  factor=9._rprec/(nz_global - 3*nz_global/4 + 1)
  sponge_top = z_i / (50._rprec * u_star)
  do k = 1, nz
    k_global = k + coord * nz
    if (k_global > 3*nz_global/4 + 1) then
      sponge(k) = sponge_top * 5._rprec**((k_global-nz_global) * factor)
    end if
  end do
$else
  sponge=0._rprec
  factor=9._rprec/(nz-3*nz/4+1)
  sponge(nz)=z_i/(50._rprec*u_star)
  do k=nz-1,3*nz/4+1,-1
     sponge(k)=sponge(k+1)/5._rprec**factor
  end do
$endif

end subroutine setsponge
end module topbc
