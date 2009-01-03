subroutine tridag(a,b,c,r,u,n)
! sc: modified for complex!
use types,only:rprec
implicit none
integer,intent(in)::n
real(kind=rprec),dimension(n),intent(in)::a,b,c
complex(kind=rprec),dimension(n),intent(in)::r
complex(kind=rprec),dimension(n),intent(out)::u
integer::j
real(kind=rprec)::bet
real(kind=rprec),dimension(n)::gam

if (b(1).eq.0._rprec) print*,'tridag: rewrite equations'
bet=b(1)
u(1)=r(1)/bet
do j=2,n
   gam(j)=c(j-1)/bet
   bet=b(j)-a(j)*gam(j)
   if(bet.eq.0._rprec) then
      print *, 'tridag failed at k=',j  
      print *,'a, b, c, gam, and bet=',a(j),b(j),c(j),gam(j),bet
   end if
   u(j)=(r(j)-a(j)*u(j-1))/bet
enddo
do j=n-1,1,-1
   u(j)=u(j)-gam(j+1)*u(j+1)
enddo
END subroutine tridag
!  (C) Copr. 1986-92 Numerical Recipes Software ]2#"0>Ya%.
