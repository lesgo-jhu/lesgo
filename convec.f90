!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
subroutine convec
!*******************************************************************************
!
! Computes the rotation convective term in physical space
!       c = - (u X vort)
! Uses 3/2-rule for dealiasing for more info see Canuto 1991 Spectral Methods,
! chapter 7
!
use types, only : rprec
use param
use sim_param, only : u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy
use sim_param, only : fxa, fya, fza, RHSx, RHSy, RHSz
use fft

implicit none

integer :: jz
integer :: jz_min
integer :: jzLo, jzHi, jz_max  ! added for full channel capabilities

real(rprec), save, allocatable, dimension(:,:,:) :: cc_big,                    &
    u_big, v_big, w_big, vort1_big, vort2_big, vort3_big
logical, save :: arrays_allocated = .false.

real(rprec) :: const

if (sgs) then
    jzLo = 2        !! necessary for LES or else blows up ....?
    jzHi = nz-1     !! can remove after testing
else
    jzLo = 1        !! for DNS
    jzHi = nz-1     !! can remove after testing
endif

if( .not. arrays_allocated ) then
   allocate( cc_big( ld_big,ny2,nz ) )
   allocate( u_big(ld_big, ny2, 0:nz) )
   allocate( v_big(ld_big, ny2, 0:nz) )
   allocate( w_big(ld_big, ny2, 0:nz) )
   allocate( vort1_big( ld_big,ny2,nz ) )
   allocate( vort2_big( ld_big,ny2,nz ) )
   allocate( vort3_big( ld_big,ny2,nz ) )
   arrays_allocated = .true.
endif

! Recall dudz, and dvdz are on UVP node for k=1 only
! So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! Loop through horizontal slices
! MPI: u_big, v_big needed at jz = 0, w_big not needed though
! MPI: could get u{1,2}_big
const = 1._rprec/(nx*ny)
do jz = 0, nz
    ! use fxa, fya, fza for temp storage
    fxa(1:nx,1:ny,jz)=const*u(1:nx,1:ny,jz)
    fya(1:nx,1:ny,jz)=const*v(1:nx,1:ny,jz)
    fza(1:nx,1:ny,jz)=const*w(1:nx,1:ny,jz)

    ! do forward fft on normal-size arrays
    call dfftw_execute_dft_r2c(forw, fxa(:,:,jz), fxa(:,:,jz))
    call dfftw_execute_dft_r2c(forw, fya(:,:,jz), fya(:,:,jz))
    call dfftw_execute_dft_r2c(forw, fza(:,:,jz), fza(:,:,jz))

    ! zero pad: padd takes care of the oddballs
    call padd(u_big(:,:,jz), fxa(:,:,jz))
    call padd(v_big(:,:,jz), fya(:,:,jz))
    call padd(w_big(:,:,jz), fza(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back_big, u_big(:,:,jz), u_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, v_big(:,:,jz), v_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, w_big(:,:,jz), w_big(:,:,jz))
end do


! Do the same thing with the vorticity
do jz = 1, nz
    ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
    ! definition of the vorticity in that case which also interpolates
    ! dwdx, dwdy to the u-node at jz=1
    if ( (coord == 0) .and. (jz == 1) ) then

        select case (lbc_mom)
        ! Stress free
        case (0)
            fxa(:, :, 1) = 0._rprec
            fya(:, :, 1) = 0._rprec

        ! Wall (all cases >= 1)
        case (1:)
            ! dwdy(jz=1) should be 0, so we can use this
            fxa(1:nx,1:ny,1) = const * ( 0.5_rprec * (dwdy(1:nx,1:ny,1) +             &
                dwdy(1:nx,1:ny,2))  - dvdz(1:nx,1:ny,1) )
            ! dwdx(jz=1) should be 0, so we can use this
            fya(1:nx,1:ny,1) = const * ( dudz(1:nx,1:ny,1) -                          &
                0.5_rprec * (dwdx(1:nx,1:ny,1) + dwdx(1:nx,1:ny,2)) )

        end select
  endif

  if ( (coord == nproc-1) .and. (jz == nz) ) then

     select case (ubc_mom)

     ! Stress free
     case (0)

         fxa(:, :, nz) = 0._rprec
         fya(:, :, nz) = 0._rprec

      ! No-slip and wall model
      case (1:)

         ! dwdy(jz=1) should be 0, so we could use this
         ! this fxa = vort1 is actually uvp nz-1 but stored as w nz
         fxa(1:nx,1:ny,nz) = const * ( 0.5_rprec * (dwdy(1:nx,1:ny,nz-1) +            &
            dwdy(1:nx,1:ny,nz)) - dvdz(1:nx,1:ny,nz-1) )
         ! dwdx(jz=1) should be 0, so we could use this
         ! this fya = vort2 is actually uvp nz-1 but stored as w nz
         fya(1:nx,1:ny,nz) = const * ( dudz(1:nx,1:ny,nz-1) -                         &
            0.5_rprec * (dwdx(1:nx,1:ny,nz-1) + dwdx(1:nx,1:ny,nz)) )

      end select
   endif

    ! very kludgy -- fix later      !! channel
    if (.not.(coord==0 .and. jz==1) .and. .not. (ubc_mom>0 .and.               &
        coord==nproc-1 .and. jz==nz)  ) then
        fxa(1:nx,1:ny,jz)=const*(dwdy(1:nx,1:ny,jz)-dvdz(1:nx,1:ny,jz))
        fya(1:nx,1:ny,jz)=const*(dudz(1:nx,1:ny,jz)-dwdx(1:nx,1:ny,jz))
    end if

    fza(1:nx,1:ny,jz)=const*(dvdx(1:nx,1:ny,jz)-dudy(1:nx,1:ny,jz))

    ! do forward fft on normal-size arrays
    call dfftw_execute_dft_r2c(forw, fxa(:,:,jz), fxa(:,:,jz))
    call dfftw_execute_dft_r2c(forw, fya(:,:,jz), fya(:,:,jz))
    call dfftw_execute_dft_r2c(forw, fza(:,:,jz), fza(:,:,jz))
    call padd(vort1_big(:,:,jz), fxa(:,:,jz))
    call padd(vort2_big(:,:,jz), fya(:,:,jz))
    call padd(vort3_big(:,:,jz), fza(:,:,jz))

    ! Back to physical space
    ! the normalization should be ok...
    call dfftw_execute_dft_c2r(back_big, vort1_big(:,:,jz), vort1_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, vort2_big(:,:,jz), vort2_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, vort3_big(:,:,jz), vort3_big(:,:,jz))
end do

! fxa
! redefinition of const
const=1._rprec/(nx2*ny2)

if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,1)=const*(v_big(:,:,1)*(-vort3_big(:,:,1))&
       +0.5_rprec*w_big(:,:,2)*(vort2_big(:,:,jzLo)))   ! (default index was 2)
    !--vort2(jz=1) is located on uvp-node        ^  try with 1 (experimental)
    !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    !  above the wall (could arguably be 0.25 * w(:,:,2))
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1 ) then  ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,nz-1)=const*(v_big(:,:,nz-1)*(-vort3_big(:,:,nz-1))&
        +0.5_rprec*w_big(:,:,nz-1)*(vort2_big(:,:,jzHi)))   ! channel
    !--vort2(jz=1) is located on uvp-node        ^  try with nz-1 (experimental)
    !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    !  below the wall (could arguably be 0.25 * w(:,:,2))

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max    !nz-1   ! channel
    cc_big(:,:,jz)=const*(v_big(:,:,jz)*(-vort3_big(:,:,jz))&
        +0.5_rprec*(w_big(:,:,jz+1)*(vort2_big(:,:,jz+1))&
        +w_big(:,:,jz)*(vort2_big(:,:,jz))))
end do

! Loop through horizontal slices
do jz=1,nz-1
    call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz),cc_big(:,:,jz))
    ! un-zero pad
    ! note: cc_big is going into fxa
    call unpadd(fxa(:,:,jz),cc_big(:,:,jz))
    ! Back to physical space
    call dfftw_execute_dft_c2r(back, fxa(:,:,jz), fxa(:,:,jz))
end do

! fya
! const should be 1./(nx2*ny2) here
if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,1)=const*(u_big(:,:,1)*(vort3_big(:,:,1))&
        +0.5_rprec*w_big(:,:,2)*(-vort1_big(:,:,jzLo)))   ! channel
    !--vort1(jz=1) is uvp-node                    ^ try with 1 (experimental)
    !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    !  above the wall (could arguably be 0.25 * w(:,:,2))
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then   ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,nz-1)=const*(u_big(:,:,nz-1)*(vort3_big(:,:,nz-1))&
        +0.5_rprec*w_big(:,:,nz-1)*(-vort1_big(:,:,jzHi)))    ! channel
    !--vort1(jz=1) is uvp-node                    ^ try with nz-1 (experimental)
    !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    !  below the wall

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max  !nz - 1   ! channel
   cc_big(:,:,jz)=const*(u_big(:,:,jz)*(vort3_big(:,:,jz))&
        +0.5_rprec*(w_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
        +w_big(:,:,jz)*(-vort1_big(:,:,jz))))
end do

do jz=1,nz-1
    call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz), cc_big(:,:,jz))
    ! un-zero pad
    ! note: cc_big is going into fya
    call unpadd(fya(:,:,jz), cc_big(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back, fya(:,:,jz), fya(:,:,jz))
end do

! fza

if (coord == 0) then
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
    !--earlier comment is also wrong, it is true that fza = 0 at both walls and
    ! slip BC
    cc_big(:,:,1)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then     ! channel
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--but now we do solve z-eqn (w-eqn) at top wall --pj
    !--earlier comment is also wrong, it is true that fza = 0 at both walls and
    ! slip BC
    cc_big(:,:,nz)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_max = nz-1
else
    jz_max = nz-1   !! or nz ?       ! channel
end if

!#ifdef PPMPI
!  if (coord == nproc-1) then
!    cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!    jz_max = nz - 1
!  else
!    jz_max = nz
!  endif
!#else
!  cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!  jz_max = nz - 1
!#endif

! channel
do jz = jz_min, jz_max    !nz - 1
    cc_big(:,:,jz) = const*0.5_rprec*(                                         &
        (u_big(:,:,jz)+u_big(:,:,jz-1))*(-vort2_big(:,:,jz))                   &
        +(v_big(:,:,jz)+v_big(:,:,jz-1))*(vort1_big(:,:,jz)))
end do

! Loop through horizontal slices
do jz=1,nz !nz - 1
    call dfftw_execute_dft_r2c(forw_big,cc_big(:,:,jz),cc_big(:,:,jz))

    ! un-zero pad
    ! note: cc_big is going into fza!!!!
    call unpadd(fza(:,:,jz),cc_big(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back,fza(:,:,jz),   fza(:,:,jz))
end do

#ifdef PPMPI
#ifdef PPSAFETYMODE
fxa(:, :, 0) = BOGUS
fya(:, :, 0) = BOGUS
fza(: ,:, 0) = BOGUS
#endif
#endif

!--top level is not valid
#ifdef PPSAFETYMODE
fxa(:, :, nz) = BOGUS
fya(:, :, nz) = BOGUS
if(coord<nproc-1) fza(:, :, nz) = BOGUS
#endif

RHSx = fxa(1:nx,1:ny,0:nz)
RHSy = fya(1:nx,1:ny,0:nz)
RHSz = fza(1:nx,1:ny,0:nz)

end subroutine convec
