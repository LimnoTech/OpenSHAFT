subroutine ReadStorageCurves
    
use GlobalVariables

implicit none

integer(i4) :: kt,kdummy,kk

!> scan storage curve file for number of curve tables and longest table
MaxNStorageCurveLines=0
NStorageCurves=0
open (unit=14,file='StorageCurves.inp',status='old',err=903)
do iskip=1,6
  read (14,*)
end do
do while(not(eof(14)))
  read (14,*) kt,idummy
  if (idummy > MaxNStorageCurveLines) MaxNStorageCurveLines=idummy
  NStorageCurves=NStorageCurves+1
  do iread=1,idummy
    read (14,*) FloatDummy,FloatDummy
  end do
end do
!write (*,*) 'MaxNStorageCurveLines =',MaxNStorageCurveLines
!write (*,*)
close (unit=14)

!> allocate storage curve arrays then read data
allocate(StorageCurveID(NStorageCurves),NStorageCurveLines(NStorageCurves))
allocate(StorageCurveH(MaxNStorageCurveLines,NStorageCurves),StorageCurveA(MaxNStorageCurveLines,NStorageCurves))
allocate(ChamberVol(MaxNStorageCurveLines,NStorageCurves),ChamberVolH(MaxNStorageCurveLines,NStorageCurves),ChamberVolA(MaxNStorageCurveLines,NStorageCurves))
open (unit=15,file='StorageCurves.inp',status='old',err=903)
do iskip=1,6
  read (15,*)
end do
do kt=1,NStorageCurves
  read (15,*) StorageCurveID(kt),NStorageCurveLines(kt)
  do iread=1,NStorageCurveLines(kt)
    read (15,*) StorageCurveH(iread,kt),StorageCurveA(iread,kt)
  end do
!  write (*,612) NStorageCurveLines(kt),StorageCurveID(kt)
end do
close (unit=15)

!> assign storage curves to proper junctions

do kt=1,NStorageCurves
  do kdummy=1,Njuncs
    if (StorageCurveID(kt) == kdummy) then
	  Junc(kdummy)%StorageMap=kt
      exit
	end if
  end do
end do
    
!> check to see if there are storage curves for each junction that expects one

do kdummy=1,Njuncs
  if ((Junc(kdummy)%IsStorageCurve) .and. (Junc(kdummy)%StorageMap == 0)) then
    write (*,*) 'shaft storage curve for junction',kdummy,' is missing'
	IOError=.true.
  end if
end do

!> create chamber volume curves if needed

do kdummy=1,Njuncs
  if ((Junc(kdummy)%IsStorageCurve) .and. (Junc(kdummy)%Option == -8)) then
    kk=Junc(kdummy)%StorageMap
    do kt=1,NStorageCurveLines(kk)
       ChamberVolH(NStorageCurveLines(kk)-kt+1,kk)=max(Junc(kdummy)%OverflowElev-Junc(kdummy)%Elev-StorageCurveH(kt,kk),0.0_wp)
       ChamberVolA(NStorageCurveLines(kk)-kt+1,kk)=StorageCurveA(kt,kk)
    end do
    do kt=2,NStorageCurveLines(kk)
      ChamberVol(kt,kk)=(ChamberVolH(kt,kk)-ChamberVolH(kt-1,kk))/3.0_wp*(ChamberVolA(kt,kk)+ChamberVolA(kt-1,kk)+sqrt(ChamberVolA(kt,kk)*ChamberVolA(kt-1,kk)))+ChamberVol(kt-1,kk)
    end do
  end if
end do

go to 999

903 write(*,*) 'StorageCurves.inp not found'
IOError=.true.

999 continue
    
return

612 format ('+read ',i6,' entries for junction ',i5)

end subroutine