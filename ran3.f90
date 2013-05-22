!!
!!  Copyright 2009,2012 Johns Hopkins University
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

!*******************************************************************************
FUNCTION ran3(idum)
!*******************************************************************************
INTEGER*4 idum
INTEGER*4 MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
REAL*8 ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER*4 i,iff,ii,inext,inextp,k
INTEGER*4 mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
if(idum.lt.0.or.iff.eq.0)then
  iff=1
  mj=MSEED-iabs(idum)
  mj=mod(mj,MBIG)
  ma(55)=mj
  mk=1
  do 11 i=1,54
    ii=mod(21*i,55)
    ma(ii)=mk
    mk=mj-mk
    if(mk.lt.MZ)mk=mk+MBIG
    mj=ma(ii)
11      continue
  do 13 k=1,4
    do 12 i=1,55
      ma(i)=ma(i)-ma(1+mod(i+30,55))
      if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
  inext=0
  inextp=31
  idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
return
END
