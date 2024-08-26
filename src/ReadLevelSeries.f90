subroutine ReadLevelSeries
    
use GlobalVariables

implicit none

integer(i4) :: kk,kt

allocate(LevelSeriesID(NLevelSer),NLevelTimes(NLevelSer),LevelTSIndex(NLevelSer))

! open level series file and scan to find longest time series
MaxNLevels=0
open (unit=13,file='LevelTimeSeries.inp',status='old',err=54)
do iskip=1,6
  read (13,*)
end do
do kk=1,NLevelSer
  read (13,*) LevelSeriesID(kk),NLevelTimes(kk)
  if (NLevelTimes(kk) > MaxNLevels) MaxNLevels=NLevelTimes(kk)
  do i=1,NLevelTimes(kk)
    read (13,*) FloatDummy,FloatDummy
  end do
end do
!write (*,*) 'MaxNLevels =',MaxNLevels
!write (*,*)
close (unit=13)

! now we know what to allocate for level series
allocate(LevelTime(MaxNLevels,NLevelSer),Level(MaxNLevels,NLevelSer))

LevelTSIndex=1
LevelSeriesID=0
NLevelTimes=0
LevelTime=0.0D+0
Level=0.0D+0

! open the file again and actually read the series this time
open (unit=13,file='LevelTimeSeries.inp',status='old',err=54)
do iskip=1,6
  read (13,*)
end do
do kk=1,NLevelSer
  write (*,*)
  read (13,*) LevelSeriesID(kk),NLevelTimes(kk)
  do i=1,NLevelTimes(kk)
    read (13,*) LevelTime(i,kk),Level(i,kk)
  end do
!  write (*,612) NLevelTimes(kk),LevelSeriesID(kk)
end do
close (unit=13)
write (*,*)

! assign level time series to proper junctions

do kk=1,NLevelSer
  do k=1,Njuncs
    if (LevelSeriesID(kk) == k) then
	  Junc(k)%LevelTSMap=kk
	end if
  end do
end do
    
! check to see if level time series are long enough
do k=1,Njuncs
  if (Junc(k)%LevelTSMap /= 0) then
    if (LevelTime(NLevelTimes(Junc(k)%LevelTSMap),Junc(k)%LevelTSMap) < Tsim) then
      write (*,*) 'level time series for junction',k,' ends before simulation'
	  IOError=.true.
    end if
  end if
end do

! check to see if there are level time series for each junction that expects one

do k=1,Njuncs
  if ((Junc(k)%Option == -6) .and. (Junc(k)%LevelTSMap == 0)) then
    write (*,*) 'level time series for junction',k,' is missing'
	IOError=.true.
  end if
end do

go to 999

54 write(*,*) 'LevelTimeSeries.inp not found'
IOError=.true.

999 continue
    
return

end subroutine ReadLevelSeries