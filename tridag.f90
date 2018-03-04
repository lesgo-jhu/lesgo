!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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
