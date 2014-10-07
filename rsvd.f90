subroutine rsvd(vegin1,vegin2) 

! This subroutine reads SSiB vegetation data in the format of the original 
! SSiB3 code given to me by Xue.  The outputs include vegin1 and vegin2 arrays
! which simply contain the variables stored sequentially.

implicit none

integer j,i
real,intent(inout) :: vegin1(58),vegin2(12,12)

! CODE TO READ VEGIN1 DATA

! here is the old way that the vegetation data was read...
!      include 'comsib'
!      READ(1,*)
!      READ(1,*) (TRAN(1,IW,1), IW=1,3), (TRAN(1,IW,2), IW=1,3)
!      READ(1,*) (TRAN(2,IW,1), IW=1,3), (TRAN(2,IW,2), IW=1,3)
!      READ(1,*) (REF (1,IW,1), IW=1,3), (REF (1,IW,2), IW=1,3)
!      READ(1,*) (REF (2,IW,1), IW=1,3), (REF (2,IW,2), IW=1,3)
!      READ (1,*) (RSTPAR(1,IWAVE), IWAVE=1,3),&
!        (RSTPAR(2,IWAVE), IWAVE=1,3)
!      READ (1,*) (SOREF(IWAVE), IWAVE=1,3),(CHIL(IV), IV=1,2)
!      READ (1,*) (TOPT(IV), IV=1,2),
!        (TLL(IV),  IV=1,2), (TU(IV),   IV=1,2)
!      READ (1,*) (DEFAC(IV),IV=1,2), (PH1(IV),  IV=1,2),
!        (PH2(IV),  IV=1,2)
!      READ (1,*) (ROOTD(IV),IV=1,2),BEE,PHSAT
!      READ (1,*) SATCO, POROS, SLOPE
!      READ (1,*) (ZDEPTH(IDEP), IDEP=1,3)
!      READ(1,*) ZWIND


!open(unit=1,file='sahel11a.data',status='old')
open(unit=1,file='veg_const.in',status='old')
read(1,*)
read(1,*) (vegin1(i),i=1,6)     ! tran(1,iw,il)
read(1,*) (vegin1(i),i=7,12)    ! tran(2,iw,il)
read(1,*) (vegin1(i),i=13,18)   ! ref(1,iw,il)
read(1,*) (vegin1(i),i=19,24)   ! ref(2,iw,il)
read(1,*) (vegin1(i),i=25,30)   ! rstpar(iv,iw)
read(1,*) (vegin1(i),i=31,35)   ! soref(iw),chil(iv) 
read(1,*) (vegin1(i),i=36,41)   ! topt(iv),tll(iv),tu(iv)
read(1,*) (vegin1(i),i=42,47)   ! defac(iv),ph1(iv),ph2(iv)
read(1,*) (vegin1(i),i=48,51)   ! rootd(iv),bee,phsat
read(1,*) (vegin1(i),i=52,54)   ! satco,poros,slope
read(1,*) (vegin1(i),i=55,57)   ! zdepth(idep)
read(1,*) vegin1(58)            ! zwind
close(1)

! CODE TO READ VEGIN2 DATA

! This is the old code...
!      READ(2,*)
!      READ(2,*) (ZLT(IV), IV=1,2),(GREEN(IV), IV=1,2),Z2,Z1
!      READ(2,*) (VCOVER(IV), IV=1,2),Z0,D
!      READ(2,*) RBC,  RDC

!open(unit=2,file='sahel11.datm',status='old')
open(unit=2,file='veg_dynamic.in',status='old')

do j=1,12
  read(2,*)
  read(2,*) (vegin2(j,i),i=1,6)   !zlt(iv),green(iv),z2,z1
  read(2,*) (vegin2(j,i),i=7,10)  !vcover(iv),z0,d
  read(2,*) (vegin2(j,i),i=11,12) !rbc,rdc
end do
close(2)

end subroutine rsvd


