!***********************************************************************
subroutine convec ()
!***********************************************************************
!
! c = - (u X vort)
!...Compute rotational convective term in physical space  (X-mom eq.)
!...For staggered grid
!...u2, du2, and du4 are on W nodes, rest on UVP nodes
!...(Fixed k=Nz plane on 1/24)
!...Corrected near wall on 4/14/96 {w(DZ/2)=0.5*w(DZ)}
! sc: better to have 0.25*w(dz), since really parabolic.
!--provides cx, cy, cz at 1:nz-1
!
! uses 3/2-rule for dealiasing 
!-- for more info see Canuto 1991 Spectral Methods (0387522050), chapter 7
!
use types,only:rprec
use param
use sim_param, only : u1=>u, u2=>v, u3=>w, du1d2=>dudy, du1d3=>dudz,   &
                      du2d1=>dvdx, du2d3=>dvdz, du3d1=>dwdx, du3d2=>dwdy
use sim_param, only : cx => RHSx, cy => RHSy, cz => RHSz
use fft

$if ($DEBUG)
use debug_mod
$endif

implicit none

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer::jz
integer :: jz_min

! !--save forces heap storage
! real(kind=rprec), save, dimension(ld_big,ny2,nz)::cc_big
! !--save forces heap storage
! real (rprec), save, dimension (ld_big, ny2, lbz:nz) :: u1_big, u2_big, u3_big
! !--MPI: only u1_big(0:nz-1), u2_big(0:nz-1), u3_big(1:nz) are used
! !--save forces heap storage 
! real (rprec), save, dimension (ld_big, ny2, nz) :: vort1_big, vort2_big,  &
!                                                    vort3_big

real(rprec), save, allocatable, dimension(:,:,:)::cc_big
real (rprec), save, allocatable, dimension (:,:,:) :: u1_big, u2_big, u3_big
!--MPI: only u1_big(0:nz-1), u2_big(0:nz-1), u3_big(1:nz) are used
real (rprec), save, allocatable, dimension (:, :, :) :: vort1_big, vort2_big, vort3_big
logical, save :: arrays_allocated = .false. 

real(kind=rprec) :: const

$if ($VERBOSE)
write (*, *) 'started convec'
$endif

if( .not. arrays_allocated ) then

   allocate( cc_big( ld_big,ny2,nz ) )
   allocate( u1_big(ld_big, ny2, lbz:nz), &
        u2_big(ld_big, ny2, lbz:nz), &
        u3_big(ld_big, ny2, lbz:nz) )
   allocate( vort1_big( ld_big,ny2,nz ), &
        vort2_big( ld_big,ny2,nz ), &
        vort3_big( ld_big,ny2,nz ) )

   arrays_allocated = .true. 

endif

!...Recall dudz, and dvdz are on UVP node for k=1 only
!...So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! sc: it would be nice to NOT have to loop through slices here
! Loop through horizontal slices

const=1._rprec/(nx*ny)
!$omp parallel do default(shared) private(jz)		
do jz = lbz, nz
  !--MPI: u1_big, u2_big needed at jz = 0, u3_big not needed though
  !--MPI: could get u{1,2}_big
! use cx,cy,cz for temp storage here!      
  $if ($FFTW3)
  call dfftw_execute_dft_r2c(plan_forward,const*u1(1:nx,1:ny,jz),cx(1:nx+2,1:ny,jz))
  call dfftw_execute_dft_r2c(plan_forward,const*u2(1:nx,1:ny,jz),cy(1:nx+2,1:ny,jz))
  call dfftw_execute_dft_r2c(plan_forward,const*u3(1:nx,1:ny,jz),cz(1:nx+2,1:ny,jz))
  $else
! do forward fft on normal-size arrays
   cx(:,:,jz)=const*u1(:,:,jz)
   cy(:,:,jz)=const*u2(:,:,jz)
   cz(:,:,jz)=const*u3(:,:,jz)
   call rfftwnd_f77_one_real_to_complex(forw,cx(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_real_to_complex(forw,cy(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_real_to_complex(forw,cz(:,:,jz),fftwNull_p)
  $endif 
! zero pad: padd takes care of the oddballs
   call padd(u1_big(:,:,jz),cx(:,:,jz))
   call padd(u2_big(:,:,jz),cy(:,:,jz))
   call padd(u3_big(:,:,jz),cz(:,:,jz))
! Back to physical space
! the normalization should be ok...
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(plan_backward_big,u1_big(1:nx2+2,1:ny2,jz),   u1_big(1:nx2,1:ny2,jz))
  call dfftw_execute_dft_c2r(plan_backward_big,u2_big(1:nx2+2,1:ny2,jz),   u2_big(1:nx2,1:ny2,jz))
  call dfftw_execute_dft_c2r(plan_backward_big,u3_big(1:nx2+2,1:ny2,jz),   u3_big(1:nx2,1:ny2,jz))    
  $else
  call rfftwnd_f77_one_complex_to_real(back_big,u1_big(:,:,jz),fftwNull_p)
  call rfftwnd_f77_one_complex_to_real(back_big,u2_big(:,:,jz),fftwNull_p)
  call rfftwnd_f77_one_complex_to_real(back_big,u3_big(:,:,jz),fftwNull_p)
  $endif
end do

do jz = 1, nz
! now do the same, but with the vorticity!
!--if du1d3, du2d3 are on u-nodes for jz=1, then we need a special
!  definition of the vorticity in that case which also interpolates
!  du3d1, du3d2 to the u-node at jz=1
   if ( (coord == 0) .and. (jz == 1) ) then
        
     select case (lbc_mom)

     ! Stress free
     case (0)

         cx(:, :, 1) = 0._rprec
         cy(:, :, 1) = 0._rprec

      ! Wall
      case (1)

         !--du3d2(jz=1) should be 0, so we could use this
         cx(:, :, 1) = const * ( 0.5_rprec * (du3d2(:, :, 1) +  &
                                              du3d2(:, :, 2))   &
                                 - du2d3(:, :, 1) )
         !--du3d1(jz=1) should be 0, so we could use this
         cy(:, :, 1) = const * ( du1d3(:, :, 1) -               &
                                 0.5_rprec * (du3d1(:, :, 1) +  &
                                              du3d1(:, :, 2)) )

     end select

   else
   
     cx(:,:,jz)=const*(du3d2(:,:,jz)-du2d3(:,:,jz))
     cy(:,:,jz)=const*(du1d3(:,:,jz)-du3d1(:,:,jz))

   end if

   cz(:,:,jz)=const*(du2d1(:,:,jz)-du1d2(:,:,jz))

  $if ($FFTW3)
  in2(1:nx,1:ny)=cx(1:nx,1:ny,jz)
  call dfftw_execute_dft_r2c(plan_forward,in2(1:nx,1:ny),cx(1:nx+2,1:ny,jz))
  in2(1:nx,1:ny)=cy(1:nx,1:ny,jz)
  call dfftw_execute_dft_r2c(plan_forward,in2(1:nx,1:ny),cy(1:nx+2,1:ny,jz))
  in2(1:nx,1:ny)=cz(1:nx,1:ny,jz)
  call dfftw_execute_dft_r2c(plan_forward,in2(1:nx,1:ny),cz(1:nx+2,1:ny,jz))
  $else
! do forward fft on normal-size arrays
   call rfftwnd_f77_one_real_to_complex(forw,cx(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_real_to_complex(forw,cy(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_real_to_complex(forw,cz(:,:,jz),fftwNull_p)
  $endif 
   call padd(vort1_big(:,:,jz),cx(:,:,jz))
   call padd(vort2_big(:,:,jz),cy(:,:,jz))
   call padd(vort3_big(:,:,jz),cz(:,:,jz))

! Back to physical space
! the normalization should be ok...
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(plan_backward_big,vort1_big(1:nx2+2,1:ny2,jz),   vort1_big(1:nx2,1:ny2,jz))
  call dfftw_execute_dft_c2r(plan_backward_big,vort2_big(1:nx2+2,1:ny2,jz),   vort2_big(1:nx2,1:ny2,jz))
  call dfftw_execute_dft_c2r(plan_backward_big,vort3_big(1:nx2+2,1:ny2,jz),   vort3_big(1:nx2,1:ny2,jz))
  $else
  call rfftwnd_f77_one_complex_to_real(back_big,vort1_big(:,:,jz),fftwNull_p)
  call rfftwnd_f77_one_complex_to_real(back_big,vort2_big(:,:,jz),fftwNull_p)
  call rfftwnd_f77_one_complex_to_real(back_big,vort3_big(:,:,jz),fftwNull_p)
  $endif
end do
!$omp end parallel do

! CX
! redefinition of const
const=1._rprec/(nx2*ny2)

if (coord == 0) then
  ! the cc's contain the normalization factor for the upcoming fft's
  cc_big(:,:,1)=const*(u2_big(:,:,1)*(-vort3_big(:,:,1))&
       +0.5_rprec*u3_big(:,:,2)*(vort2_big(:,:,2)))
  !--vort2(jz=1) is located on uvp-node        ^  try with 1 (experimental)
  !--the 0.5 * u3(:,:,2) is the interpolation of u3 to the first uvp node
  !  above the wall (could arguably be 0.25 * u3(:,:,2))

  jz_min = 2
else
  jz_min = 1
end if

!$omp parallel do default(shared) private(jz)
do jz=jz_min,nz-1
   cc_big(:,:,jz)=const*(u2_big(:,:,jz)*(-vort3_big(:,:,jz))&
        +0.5_rprec*(u3_big(:,:,jz+1)*(vort2_big(:,:,jz+1))&
        +u3_big(:,:,jz)*(vort2_big(:,:,jz))))
end do
!$omp end parallel do

! Loop through horizontal slices
!$omp parallel do default(shared) private(jz)	
do jz=1,nz-1
  $if ($FFTW3)
  inbig2(1:nx2,1:ny2)=cc_big(1:nx2,1:ny2,jz)
  call dfftw_execute_dft_r2c(plan_forward_big,inbig2(1:nx2,1:ny2),cc_big(1:nx2+2,1:ny2,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw_big,cc_big(:,:,jz),fftwNull_p)
  $endif   
! un-zero pad
! note: cc_big is going into cx!!!!
   call unpadd(cx(:,:,jz),cc_big(:,:,jz))
! Back to physical space
   $if ($FFTW3)
   call dfftw_execute_dft_c2r(plan_backward,cx(1:nx+2,1:ny,jz),   cx(1:nx,1:ny,jz))   
   $else
   call rfftwnd_f77_one_complex_to_real(back,cx(:,:,jz),fftwNull_p)
   $endif
end do
!$omp end parallel do

! CY
! const should be 1./(nx2*ny2) here

if (coord == 0) then
  ! the cc's contain the normalization factor for the upcoming fft's
  cc_big(:,:,1)=const*(u1_big(:,:,1)*(vort3_big(:,:,1))&
       +0.5_rprec*u3_big(:,:,2)*(-vort1_big(:,:,2)))
  !--vort1(jz=1) is uvp-node                    ^ try with 1 (experimental)
  !--the 0.5 * u3(:,:,2) is the interpolation of u3 to the first uvp node
  !  above the wall (could arguably be 0.25 * u3(:,:,2))

  jz_min = 2
else
  jz_min = 1
end if

!$omp parallel do default(shared) private(jz)
do jz = jz_min, nz - 1
   cc_big(:,:,jz)=const*(u1_big(:,:,jz)*(vort3_big(:,:,jz))&
        +0.5_rprec*(u3_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
        +u3_big(:,:,jz)*(-vort1_big(:,:,jz))))
end do
!$omp end parallel do

!$omp parallel do default(shared) private(jz)		
do jz=1,nz-1
  $if ($FFTW3)
  inbig2(1:nx2,1:ny2)=cc_big(1:nx2,1:ny2,jz)
  call dfftw_execute_dft_r2c(plan_forward_big,inbig2(1:nx2,1:ny2),cc_big(1:nx2+2,1:ny2,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw_big,cc_big(:,:,jz),fftwNull_p)
  $endif
 ! un-zero pad
! note: cc_big is going into cy!!!!
   call unpadd(cy(:,:,jz),cc_big(:,:,jz))

! Back to physical space
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(plan_backward,cy(1:nx+2,1:ny,jz),   cy(1:nx,1:ny,jz))     
  $else
  call rfftwnd_f77_one_complex_to_real(back,cy(:,:,jz),fftwNull_p)
  $endif
end do
!$omp end parallel do

! CZ

if (coord == 0) then
  ! There is no convective acceleration of w at wall or at top.
  !--not really true at wall, so this is an approximation?
  !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
  cc_big(:,:,1)=0._rprec

  jz_min = 2
else
  jz_min = 1
end if

!$if ($MPI)
!  if (coord == nproc-1) then
!    cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!    jz_max = nz - 1
!  else
!    jz_max = nz
!  endif
!$else
!  cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!  jz_max = nz - 1
!$endif

!$omp parallel do default(shared) private(jz)
do jz=jz_min, nz - 1
   cc_big(:,:,jz)=const*0.5_rprec*(&
        (u1_big(:,:,jz)+u1_big(:,:,jz-1))*(-vort2_big(:,:,jz))&
        +(u2_big(:,:,jz)+u2_big(:,:,jz-1))*(vort1_big(:,:,jz))&
         )
end do
!$omp end parallel do

! Loop through horizontal slices
!$omp parallel do default(shared) private(jz)		
do jz=1,nz - 1
  $if ($FFTW3)
  inbig2(1:nx2,1:ny2)=cc_big(1:nx2,1:ny2,jz)
  call dfftw_execute_dft_r2c(plan_forward_big,inbig2(1:nx2,1:ny2),cc_big(1:nx2+2,1:ny2,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw_big,cc_big(:,:,jz),fftwNull_p)
  $endif

! un-zero pad
! note: cc_big is going into cz!!!!
   call unpadd(cz(:,:,jz),cc_big(:,:,jz))

! Back to physical space
   $if ($FFTW3)
   call dfftw_execute_dft_c2r(plan_backward,cz(1:nx+2,1:ny,jz),   cz(1:nx,1:ny,jz))
   $else
   call rfftwnd_f77_one_complex_to_real(back,cz(:,:,jz),fftwNull_p)
   $endif
end do
!$omp end parallel do

$if ($MPI)
$if ($SAFETYMODE)
  cx(:, :, 0) = BOGUS
  cy(:, :, 0) = BOGUS
  cz(: ,:, 0) = BOGUS
$endif  
$endif

!--top level is not valid
$if ($SAFETYMODE)
cx(:, :, nz) = BOGUS
cy(:, :, nz) = BOGUS
cz(:, :, nz) = BOGUS
$endif

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (cx(:, :, 1:nz), 'convec.z.cx')
  call DEBUG_write (cy(:, :, 1:nz), 'convec.z.cy')
  call DEBUG_write (cz(:, :, 1:nz), 'convec.z.cz')
endif
$endif

$if ($VERBOSE)
write (*, *) 'finished convec'
$endif

end subroutine convec

