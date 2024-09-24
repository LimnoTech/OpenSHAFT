subroutine ReadRatingCurves
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Reads rating curves (optional)	
!%=====================================================================
    
use GlobalVariables

implicit none

integer(i4) :: kt,kdummy,m

!> scan rating curve file for number of curve tables and longest table
MaxNRatingCurveLines=0
NRatingCurves=0
open (unit=14,file='RatingCurves.inp',status='old',err=903)
do iskip=1,6
  read (14,*)
end do
do while(not(eof(14)))
  read (14,*) kt,idummy
  if (idummy > MaxNRatingCurveLines) MaxNRatingCurveLines=idummy
  NRatingCurves=NRatingCurves+1
  if (kt < 0) read (14,*)
  do iread=1,idummy
    read (14,*) FloatDummy,FloatDummy
  end do
end do
!write (*,*) 'MaxNRatingCurveLines =',MaxNRatingCurveLines
!write (*,*)
close (unit=14)

!> allocate rating curve arrays then read data
allocate(RatingCurveID(NRatingCurves),NRatingCurveLines(NRatingCurves),NTailWaterH(NRatingCurves))
allocate(TailWaterH(10,NRatingCurves),RatingCurveH(MaxNRatingCurveLines,NRatingCurves),RatingCurveQ(MaxNRatingCurveLines,NRatingCurves),RatingCurveQ2(MaxNRatingCurveLines,NRatingCurves))
allocate(TWRatingCurveQ(MaxNRatingCurveLines,10,NRatingCurves))
NTailwaterH=0
TailwaterH=0.0_wp
RatingCurveH=0.0_wp
RatingCurveQ=0.0_wp
RatingCurveQ2=0.0_wp
TWRatingCurveQ=0.0_wp

open (unit=14,file='RatingCurves.inp',status='old',err=903)
do iskip=1,6
  read (14,*)
end do
do kt=1,NRatingCurves
  read (14,*) RatingCurveID(kt),NRatingCurveLines(kt)
  if (RatingCurveID(kt) < 0) then
    RatingCurveID(kt)=abs(RatingCurveID(kt))
    if (NRatingCurveLines(kt) == 0) then
      read (14,*) Junc(RatingCurveID(kt))%JuncPartner,Junc(RatingCurveID(kt))%RatingCurveType,Junc(RatingCurveID(kt))%WeirLength
    else
      read (14,*) Junc(RatingCurveID(kt))%JuncPartner,Junc(RatingCurveID(kt))%RatingCurveType
    end if
  else
    Junc(RatingCurveID(kt))%RatingCurveType=1
    Junc(RatingCurveID(kt))%JuncPartner=0
  end if
  if (Junc(RatingCurveID(kt))%RatingCurveType == 4) then
    read (14,*) NTailWaterH(kt),(TailWaterH(m,kt),m=1,NTailWaterH(kt))
  end if
  do iread=1,NRatingCurveLines(kt)
    if ((Junc(RatingCurveID(kt))%RatingCurveType == 2) .or. (Junc(RatingCurveID(kt))%RatingCurveType == 3)) then
      read (14,*) RatingCurveH(iread,kt),RatingCurveQ(iread,kt),RatingCurveQ2(iread,kt)
    else if (Junc(RatingCurveID(kt))%RatingCurveType == 4) then
      read (14,*) RatingCurveH(iread,kt),(TWRatingCurveQ(iread,m,kt),m=1,NTailWaterH(kt))
    else
      read (14,*) RatingCurveH(iread,kt),RatingCurveQ(iread,kt)
    end if
  end do
!  write (*,612) NRatingCurveLines(kt),RatingCurveID(kt)
end do
close (unit=14)

!> assign rating curves to proper junctions

do kt=1,NRatingCurves
  do kdummy=1,NJuncs
    if (RatingCurveID(kt) == kdummy) then
	  Junc(kdummy)%RatingCurveMap=kt
      exit
	end if
  end do
end do
    
!> check to see if there are rating curves for each junction that expects one

do kdummy=1,NJuncs
  if ((Junc(kdummy)%IsRatingCurve) .and. (Junc(kdummy)%RatingCurveMap == 0)) then
    write (*,*) 'rating curve for junction',kdummy,' is missing'
	IOError=.true.
  end if
end do

go to 999

903 write(*,*) 'RatingCurves.inp not found'
IOError=.true.

999 continue
    
return

612 format ('+read ',i6,' entries for junction ',i5)

end subroutine