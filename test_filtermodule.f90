module test_filtermodule
use types,only:rprec
use param,only:lh,ny

private lh,ny
!TS Truely grid refinement test needs to keep the filter_size
!TS the same as that in coarse grid (Double grid size: filter_size=2. etc)
integer,parameter::filter_size=1
real(kind=rprec),dimension(lh,ny)::G_test,G_test_test
end module test_filtermodule

!--THIS IS NOT IN THE MODULE
subroutine test_filter(f,G_test)
! note: this filters in-place, so input is ruined
use types,only:rprec
use param,only:lh,ny
use fft
implicit none
! note we're treating as complex here
complex(kind=rprec), dimension(lh,ny),intent(inout)::f
real(kind=rprec), dimension(lh,ny),intent(in) :: G_test

real (rprec) :: ignore_me
!f_c=cmplx(f,0.,kind=rprec)
call rfftwnd_f77_one_real_to_complex(forw,f,ignore_me)
! the normalization is "in" G_test
f = G_test*f  ! Nyquist taken care of with G_test
call rfftwnd_f77_one_complex_to_real(back,f,ignore_me)
!f=real(f_c,kind=rprec)
end subroutine test_filter

subroutine test_filter_init(alpha,G_test)
! spectral cutoff filter at width alpha*delta
! note the normalization for FFT's is already in G! (see 1/(nx*ny))
use types,only:rprec
use param,only:lh,nx,ny,dx,dy,pi,ifilter,model
use fft
implicit none
real(kind=rprec):: alpha, delta, kc2
real(kind=rprec),dimension(lh,ny) :: G_test
G_test=1._rprec/(nx*ny)  ! normalization for the forward FFT
delta = alpha*sqrt(dx*dy)  ! "2d-delta", not full 3d one
if(ifilter==1) then ! spectral cutoff filter
   if (model==6.OR.model==7) then
      print *, 'Use Gaussian or Top-hat filter for mixed models'
      stop
   endif
! alpha is ratio of test filter to grid filter widths
   kc2 = (pi/(delta))**2
   where (real(k2, rprec) >= kc2) G_test = 0._rprec
else if(ifilter==2) then ! Gaussian filter
   G_test=exp(-(delta)**2*k2/(4._rprec*6._rprec))*G_test       
else if(ifilter==3) then ! Top-hat (Box) filter
   G_test= (sin(kx*delta/2._rprec)*sin(ky*delta/2._rprec)+1E-8)/&
        (kx*delta/2._rprec*ky*delta/2._rprec+1E-8)*G_test
endif
! since our k2 has zero at Nyquist, we have to do this by hand
   G_test(lh,:) = 0._rprec
   G_test(:,ny/2+1) = 0._rprec
end subroutine test_filter_init
