!!
!!  Copyright 2009,2010,2011 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

subroutine tridag(a,b,c,r,u,n)
! sc: modified for complex!
use types,only:rprec
implicit none
integer,intent(in)::n
real(kind=rprec),dimension(n),intent(in)::a,b,c
!complex(kind=rprec),dimension(n),intent(in)::r
!complex(kind=rprec),dimension(n),intent(out)::u
real(rprec), dimension(2*n), intent(in) :: r
real(rprec), dimension(2*n), intent(out) :: u

integer::j, ir, ii

real(kind=rprec)::bet
real(kind=rprec),dimension(n)::gam

if (b(1).eq.0._rprec) print*,'tridag: rewrite equations'
bet=b(1)
!u(1)=r(1)/bet 
u(1:2) = r(1:2)/bet 
do j=2,n

  ii = 2*j
  ir = ii - 1

  gam(j)=c(j-1)/bet
  bet=b(j)-a(j)*gam(j)
  if(bet.eq.0._rprec) then
    print *, 'tridag failed at k=',j  
    print *,'a, b, c, gam, and bet=',a(j),b(j),c(j),gam(j),bet
  end if
  
  !u(j)=(r(j)-a(j)*u(j-1))/bet
  u(ir:ii) = (r(ir:ii) - a(j)*u(ir:ii))/bet

enddo

do j=n-1,1,-1
  ii= 2*j
  ir = ii - 1
  !u(j)=u(j)-gam(j+1)*u(j+1)
  u(ir:ii) = u(ir:ii) - gam(j+1)*u(ir:ii)
enddo
END subroutine tridag
!  (C) Copr. 1986-92 Numerical Recipes Software ]2#"0>Ya%.
