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
use sim_param, only : RHSx, RHSy, RHSz
use sim_param, only : RHSx_var, RHSy_var, RHSz_var
use sim_param, only : u_big, v_big, w_big, u_big_var, v_big_var, w_big_var
use sim_param, only : u_var, v_var, w_var, vortx_big_var, vorty_big_var, vortz_big_var
use sim_param, only : vortx_big, vorty_big, vortz_big, RHSx_var, RHSy_var, RHSz_var
use sim_param, only : cc_big, cc_big_var
use fft

implicit none

integer :: jz
integer :: jz_min
integer :: jzLo, jzHi, jz_max  ! added for full channel capabilities

! real(rprec), save, allocatable, dimension(:,:,:) :: cc_big!,                    &
     ! vortx_big, vorty_big, vortz_big
! logical, save :: arrays_allocated = .false.

real(rprec) :: const

if (sgs) then
    jzLo = 2        !! necessary for LES or else blows up ....?
    jzHi = nz-1     !! can remove after testing
else
    jzLo = 1        !! for DNS
    jzHi = nz-1     !! can remove after testing
endif

! if( .not. arrays_allocated ) then
!    allocate( cc_big( ld_big,ny2,nz ) )
!    ! allocate( vortx_big( ld_big,ny2,nz ) )
!    ! allocate( vorty_big( ld_big,ny2,nz ) )
!    ! allocate( vortz_big( ld_big,ny2,nz ) )
!    arrays_allocated = .true.
! endif

! Recall dudz, and dvdz are on UVP node for k=1 only
! So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! Loop through horizontal slices
! MPI: u_big, v_big needed at jz = 0, w_big not needed though
! MPI: could get u{1,2}_big
! const = 1._rprec/(nx*ny)
! do jz = 0, nz
!     use RHSx, RHSy, RHSz for temp storage
!     RHSx(1:nx,1:ny,jz)=const*u(1:nx,1:ny,jz)
!     RHSy(1:nx,1:ny,jz)=const*v(1:nx,1:ny,jz)
!     RHSz(1:nx,1:ny,jz)=const*w(1:nx,1:ny,jz)
!
!     do forward fft on normal-size arrays
!     call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
!     call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
!     call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
!
!     zero pad: padd takes care of the oddballs
!     call padd(u_big(:,:,jz), RHSx(:,:,jz))
!     call padd(v_big(:,:,jz), RHSy(:,:,jz))
!     call padd(w_big(:,:,jz), RHSz(:,:,jz))
!
!     Back to physical space
!     call dfftw_execute_dft_c2r(back_big, u_big(:,:,jz), u_big(:,:,jz))
!     call dfftw_execute_dft_c2r(back_big, v_big(:,:,jz), v_big(:,:,jz))
!     call dfftw_execute_dft_c2r(back_big, w_big(:,:,jz), w_big(:,:,jz))
! end do

! Calculate vorticity, using RHSx, RHSy, and RHSz as temporary storage locations
RHSx(:,:,1:nz) = dwdy(:,:,1:nz) - dvdz(:,:,1:nz)
RHSy(:,:,1:nz) = dudz(:,:,1:nz) - dwdx(:,:,1:nz)
RHSz(:,:,1:nz) = dvdx(:,:,1:nz) - dudy(:,:,1:nz)

! Bottom boundary condition
if (coord == 0) then
    if (lbc_mom == 0) then
        ! Stress free
        RHSx(:,:,1) = 0._rprec
        RHSy(:,:,1) = 0._rprec
    else
        ! Wall
        ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
        ! definition of the vorticity in that case which also interpolates
        ! dwdx, dwdy to the u-node at jz=1
        ! dwdy(1) and dwdx(1) should be 0
        RHSx(:,:,1) = 0.5_rprec*(dwdy(:,:,1) + dwdy(:,:,2)) - dvdz(:,:,1)
        RHSy(:,:,1) = dudz(:,:,1) - 0.5_rprec*(dwdx(:,:,1) + dwdx(:,:,2))
    endif
end if

! Top boundary condition
if (coord == nproc-1) then
    if (ubc_mom == 0) then
        ! Stress free
        RHSx(:,:,nz) = 0._rprec
        RHSy(:,:,nz) = 0._rprec
    else
        ! Wall
        ! if dudz, dvdz are on u-nodes for jz=nz, then we need a special
        ! definition of the vorticity in that case which also interpolates
        ! dwdx, dwdy to the u-node at jz=nz
        ! dwdy(1) and dwdx(1) should be 0
        RHSx(:,:,nz) = 0.5_rprec*(dwdy(:,:,nz-1) + dwdy(:,:,nz)) - dvdz(:,:,nz-1)
        RHSy(:,:,nz) = dudz(:,:,nz-1) - 0.5_rprec*(dwdx(:,:,nz-1) + dwdx(:,:,nz))
    endif
end if

! Transfer onto larger grid
call u_var%padd(u_big_var)
call v_var%padd(v_big_var)
call w_var%padd(w_big_var)
call RHSx_var%padd(vortx_big_var)
call RHSy_var%padd(vorty_big_var)
call RHSz_var%padd(vortz_big_var)

! ! Do the same thing with the vorticity
! do jz = 1, nz
!     ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
!     ! definition of the vorticity in that case which also interpolates
!     ! dwdx, dwdy to the u-node at jz=1
!     if ( (coord == 0) .and. (jz == 1) ) then
!
!         select case (lbc_mom)
!         ! Stress free
!         case (0)
!             RHSx(:, :, 1) = 0._rprec
!             RHSy(:, :, 1) = 0._rprec
!
!         ! Wall (all cases >= 1)
!         case (1:)
!             ! dwdy(jz=1) should be 0, so we can use this
!             RHSx(1:nx,1:ny,1) = const * ( 0.5_rprec * (dwdy(1:nx,1:ny,1) +             &
!                 dwdy(1:nx,1:ny,2))  - dvdz(1:nx,1:ny,1) )
!             ! dwdx(jz=1) should be 0, so we can use this
!             RHSy(1:nx,1:ny,1) = const * ( dudz(1:nx,1:ny,1) -                          &
!                 0.5_rprec * (dwdx(1:nx,1:ny,1) + dwdx(1:nx,1:ny,2)) )
!
!         end select
!   endif
!
!   if ( (coord == nproc-1) .and. (jz == nz) ) then
!
!      select case (ubc_mom)
!
!      ! Stress free
!      case (0)
!
!          RHSx(:, :, nz) = 0._rprec
!          RHSy(:, :, nz) = 0._rprec
!
!       ! No-slip and wall model
!       case (1:)
!
!          ! dwdy(jz=1) should be 0, so we could use this
!          ! this RHSx = vort1 is actually uvp nz-1 but stored as w nz
!          RHSx(1:nx,1:ny,nz) = const * ( 0.5_rprec * (dwdy(1:nx,1:ny,nz-1) +            &
!             dwdy(1:nx,1:ny,nz)) - dvdz(1:nx,1:ny,nz-1) )
!          ! dwdx(jz=1) should be 0, so we could use this
!          ! this RHSy = vort2 is actually uvp nz-1 but stored as w nz
!          RHSy(1:nx,1:ny,nz) = const * ( dudz(1:nx,1:ny,nz-1) -                         &
!             0.5_rprec * (dwdx(1:nx,1:ny,nz-1) + dwdx(1:nx,1:ny,nz)) )
!
!       end select
!    endif
!
!     ! very kludgy -- fix later      !! channel
!     if (.not.(coord==0 .and. jz==1) .and. .not. (ubc_mom>0 .and.               &
!         coord==nproc-1 .and. jz==nz)  ) then
!         RHSx(1:nx,1:ny,jz)=const*(dwdy(1:nx,1:ny,jz)-dvdz(1:nx,1:ny,jz))
!         RHSy(1:nx,1:ny,jz)=const*(dudz(1:nx,1:ny,jz)-dwdx(1:nx,1:ny,jz))
!     end if
!
!     RHSz(1:nx,1:ny,jz)=const*(dvdx(1:nx,1:ny,jz)-dudy(1:nx,1:ny,jz))

    ! ! do forward fft on normal-size arrays
    ! call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
    ! call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
    ! call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    ! call padd(vortx_big(:,:,jz), RHSx(:,:,jz))
    ! call padd(vorty_big(:,:,jz), RHSy(:,:,jz))
    ! call padd(vortz_big(:,:,jz), RHSz(:,:,jz))
    !
    ! ! Back to physical space
    ! ! the normalization should be ok...
    ! call dfftw_execute_dft_c2r(back_big, vortx_big(:,:,jz), vortx_big(:,:,jz))
    ! call dfftw_execute_dft_c2r(back_big, vorty_big(:,:,jz), vorty_big(:,:,jz))
    ! call dfftw_execute_dft_c2r(back_big, vortz_big(:,:,jz), vortz_big(:,:,jz))
! end do

! RHSx
! redefinition of const
! const=1._rprec/(nx2*ny2)
!
! if (coord == 0) then
!     ! the cc's contain the normalization factor for the upcoming fft's
!     cc_big(1:nx2,1:ny2,1)=const*(v_big(1:nx2,1:ny2,1)*(-vortz_big(1:nx2,1:ny2,1))&
!        +0.5_rprec*w_big(1:nx2,1:ny2,2)*(vorty_big(1:nx2,1:ny2,jzLo)))   ! (default index was 2)
!     !--vort2(jz=1) is located on uvp-node        ^  try with 1 (experimental)
!     !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
!     !  above the wall (could arguably be 0.25 * w(:,:,2))
!     jz_min = 2
! else
!     jz_min = 1
! end if
!
! if (coord == nproc-1 ) then  ! channel
!     ! the cc's contain the normalization factor for the upcoming fft's
!     cc_big(1:nx2,1:ny2,nz-1)=const*(v_big(1:nx2,1:ny2,nz-1)*(-vortz_big(1:nx2,1:ny2,nz-1))&
!         +0.5_rprec*w_big(1:nx2,1:ny2,nz-1)*(vorty_big(1:nx2,1:ny2,jzHi)))   ! channel
!     !--vort2(jz=1) is located on uvp-node        ^  try with nz-1 (experimental)
!     !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
!     !  below the wall (could arguably be 0.25 * w(:,:,2))
!
!     jz_max = nz-2
! else
!     jz_max = nz-1
! end if
!
! do jz = jz_min, jz_max    !nz-1   ! channel
!     cc_big(1:nx2,1:ny2,jz)=const*(v_big(1:nx2,1:ny2,jz)*(-vortz_big(1:nx2,1:ny2,jz))&
!         +0.5_rprec*(w_big(1:nx2,1:ny2,jz+1)*(vorty_big(1:nx2,1:ny2,jz+1))&
!         +w_big(1:nx2,1:ny2,jz)*(vorty_big(1:nx2,1:ny2,jz))))
! end do
!
! ! Loop through horizontal slices
! do jz=1,nz-1
!     call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz),cc_big(:,:,jz))
!     ! un-zero pad
!     ! note: cc_big is going into RHSx
!     call unpadd(RHSx(:,:,jz),cc_big(:,:,jz))
!     ! Back to physical space
!     call dfftw_execute_dft_c2r(back, RHSx(:,:,jz), RHSx(:,:,jz))
! end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute (u x omega)_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc_big(:,:,0:nz-1) = -v_big(:,:,0:nz-1) * vortz_big(:,:,0:nz-1)                &
    + 0.5_rprec * (w_big(:,:,1:nz)*(vorty_big(:,:,1:nz))                       &
    + w_big(:,:,0:nz-1) * (vorty_big(:,:,0:nz-1)))

! Boundary conditions
if (coord == 0) then
    ! vort(1) is located on uvp-node
    ! the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    ! above the wall (could arguably be 0.25 * w(:,:,2))
    cc_big(:,:,1) = -v_big(:,:,1) * vortz_big(:,:,1)                           &
       + 0.5_rprec * w_big(:,:,2) * vorty_big(:,:,jzLo)
end if
if (coord == nproc-1) then
    ! vort(nz-1) is located on uvp-node
    ! the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    cc_big(:,:,nz-1) = -v_big(:,:,nz-1) * vortz_big(:,:,nz-1)                 &
       + 0.5_rprec * w_big(:,:,nz-1) * vorty_big(:,:,jzHi)
end if

! Dealias
call cc_big_var%unpadd(RHSx_var)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute (u x omega)_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc_big(:,:,0:nz-1) = u_big(:,:,0:nz-1) * vortz_big(:,:,0:nz-1)                &
    - 0.5_rprec * (w_big(:,:,1:nz)*(vortx_big(:,:,1:nz))                       &
    - w_big(:,:,0:nz-1) * (vortx_big(:,:,0:nz-1)))

! Boundary conditions
if (coord == 0) then
    ! vort(1) is located on uvp-node
    ! the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    ! above the wall (could arguably be 0.25 * w(:,:,2))
    cc_big(:,:,1) = u_big(:,:,1) * vortz_big(:,:,1)                           &
       - 0.5_rprec * w_big(:,:,2) * vortx_big(:,:,jzLo)
end if
if (coord == nproc-1) then
    ! vort(nz-1) is located on uvp-node
    ! the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    cc_big(:,:,nz-1) = u_big(:,:,nz-1) * vortz_big(:,:,nz-1)                 &
       - 0.5_rprec * w_big(:,:,nz-1) * vortx_big(:,:,jzHi)
end if

! Dealias
call cc_big_var%unpadd(RHSy_var)

! ! RHSy
! ! const should be 1./(nx2*ny2) here
! if (coord == 0) then
!     ! the cc's contain the normalization factor for the upcoming fft's
!     cc_big(1:nx2,1:ny2,1)=const*(u_big(1:nx2,1:ny2,1)*(vortz_big(1:nx2,1:ny2,1))&
!         +0.5_rprec*w_big(1:nx2,1:ny2,2)*(-vortx_big(1:nx2,1:ny2,jzLo)))   ! channel
!     !--vort1(jz=1) is uvp-node                    ^ try with 1 (experimental)
!     !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
!     !  above the wall (could arguably be 0.25 * w(:,:,2))
!     jz_min = 2
! else
!     jz_min = 1
! end if
!
! if (coord == nproc-1) then   ! channel
!     ! the cc's contain the normalization factor for the upcoming fft's
!     cc_big(1:nx2,1:ny2,nz-1)=const*(u_big(1:nx2,1:ny2,nz-1)*(vortz_big(1:nx2,1:ny2,nz-1))&
!         +0.5_rprec*w_big(1:nx2,1:ny2,nz-1)*(-vortx_big(1:nx2,1:ny2,jzHi)))    ! channel
!     !--vort1(jz=1) is uvp-node                    ^ try with nz-1 (experimental)
!     !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
!     !  below the wall
!
!     jz_max = nz-2
! else
!     jz_max = nz-1
! end if
!
! do jz = jz_min, jz_max  !nz - 1   ! channel
!    cc_big(1:nx2,1:ny2,jz)=const*(u_big(1:nx2,1:ny2,jz)*(vortz_big(1:nx2,1:ny2,jz))&
!         +0.5_rprec*(w_big(1:nx2,1:ny2,jz+1)*(-vortx_big(1:nx2,1:ny2,jz+1))&
!         +w_big(1:nx2,1:ny2,jz)*(-vortx_big(1:nx2,1:ny2,jz))))
! end do
!
! do jz=1,nz-1
!     call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz), cc_big(:,:,jz))
!     ! un-zero pad
!     ! note: cc_big is going into RHSy
!     call unpadd(RHSy(:,:,jz), cc_big(:,:,jz))
!
!     ! Back to physical space
!     call dfftw_execute_dft_c2r(back, RHSy(:,:,jz), RHSy(:,:,jz))
! end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute (u x omega)_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc_big(:,:,0:nz-1) = 0.5_rprec * (                                             &
    - (u_big(:,:,0:nz-1) + u_big(:,:,1:nz)) * vorty_big(:,:,0:nz-1)            &
    + (v_big(:,:,0:nz-1) + v_big(:,:,1:nz)) * vortx_big(:,:,0:nz-1) )

! Boundary conditions
! Some previous musings on this:
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
    !--earlier comment is also wrong, it is true that RHS = 0 at both walls and
    ! slip BC
if (coord == 0) then
    cc_big(:,:,1) = 0._rprec
end if
if (coord == nproc-1) then
    cc_big(:,:,nz-1) = 0._rprec
end if

! Dealias
call cc_big_var%unpadd(RHSz_var)
!
! ! RHSz
!
! if (coord == 0) then
!     ! There is no convective acceleration of w at wall or at top.
!     !--not really true at wall, so this is an approximation?
!     !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
!     !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
!     !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
!     ! slip BC
!     cc_big(:,:,1)=0._rprec
!     !! ^must change for Couette flow ... ?
!     jz_min = 2
! else
!     jz_min = 1
! end if
!
! if (coord == nproc-1) then     ! channel
!     ! There is no convective acceleration of w at wall or at top.
!     !--not really true at wall, so this is an approximation?
!     !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
!     !--but now we do solve z-eqn (w-eqn) at top wall --pj
!     !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
!     ! slip BC
!     cc_big(:,:,nz)=0._rprec
!     !! ^must change for Couette flow ... ?
!     jz_max = nz-1
! else
!     jz_max = nz-1   !! or nz ?       ! channel
! end if

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
!
! ! channel
! do jz = jz_min, jz_max    !nz - 1
!     cc_big(1:nx2,1:ny2,jz) = const*0.5_rprec*(                                         &
!         (u_big(1:nx2,1:ny2,jz)+u_big(1:nx2,1:ny2,jz-1))*(-vorty_big(1:nx2,1:ny2,jz))                   &
!         +(v_big(1:nx2,1:ny2,jz)+v_big(1:nx2,1:ny2,jz-1))*(vortx_big(1:nx2,1:ny2,jz)))
! end do
!
! ! Loop through horizontal slices
! do jz=1,nz !nz - 1
!     call dfftw_execute_dft_r2c(forw_big,cc_big(:,:,jz),cc_big(:,:,jz))
!
!     ! un-zero pad
!     ! note: cc_big is going into RHSz!!!!
!     call unpadd(RHSz(:,:,jz),cc_big(:,:,jz))
!
!     ! Back to physical space
!     call dfftw_execute_dft_c2r(back,RHSz(:,:,jz),   RHSz(:,:,jz))
! end do

#ifdef PPMPI
#ifdef PPSAFETYMODE
RHSx(:, :, 0) = BOGUS
RHSy(:, :, 0) = BOGUS
RHSz(: ,:, 0) = BOGUS
#endif
#endif

!--top level is not valid
#ifdef PPSAFETYMODE
RHSx(:, :, nz) = BOGUS
RHSy(:, :, nz) = BOGUS
if(coord<nproc-1) RHSz(:, :, nz) = BOGUS
#endif

end subroutine convec
