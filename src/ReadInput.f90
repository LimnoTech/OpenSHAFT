subroutine ReadInput
    
use GlobalVariables

implicit none

integer :: nc_temp,kt,jj,jdummy,StartPos
real(wp) :: fragment

IOError=.false.

open (10,file='RunControl.inp',status='old',err=900)

do iskip=1,4
    read (10,*)
end do
read (10,*) RunID
read (10,*) Tsim
read (10,*) Twrite
read (10,*) StartWrite
read (10,*) StopWrite
read (10,*) Crf
read (10,*) HFnf
read (10,*) Ratio_dx_dpipe
read (10,*) initHfrac
read (10,*) amax
read (10,*) AirExhaustSwitch
read (10,*) AirPocketSwitch
read (10,*) StartCondition
read (10,*) ShowJunc(1)
read (10,*) ShowJunc(2)
read (10,*) ShowJunc(3)
read (10,*) alpha1
read (10,*) alpha2
read (10,*) alpha4
read (10,*) ThresholdVol
read (10,*) PocketBuffer
read (10,*) NThreadsUser
read (10,*) HotStartID
close (10)

write (sj1,'(i2)') ShowJunc(1)
write (sj2,'(i2)') ShowJunc(2)
write (sj3,'(i2)') ShowJunc(3)

if (AirExhaustSwitch == 1) then
  CalcAirExhaust=.true.
  write (*,*) 'air exhaust calculation enabled'
  write (*,*)
else if (AirExhaustSwitch == 0) then
  CalcAirExhaust=.false.
  write (*,*) 'air exhaust calculation disabled'
  write (*,*)
else
  write (*,*) 'air exhaust switch set incorrectly'
  IOError=.true.
  go to 999
end if

if (AirPocketSwitch == 1) then
  TrapAirPockets=.true.
  write (*,*) 'air pocket calculation enabled'
  write (*,*)
else if (AirPocketSwitch == 0) then
  TrapAirPockets=.false.
  write (*,*) 'air pocket calculation disabled'
  write (*,*)
!else if (AirPocketSwitch == 2) then
!  DoAirPockets=.true.
!  write (*,*) 'moving air pocket calculation enabled - nice!'
!  write (*,*)
else
  write (*,*) 'air pocket switch set incorrectly'
  IOError=.true.
  go to 999
end if

if (StartCondition == 0) then
  IsHGLInit=.false.
  IsHotStart=.false.
else if (StartCondition == 1) then
  IsHGLInit=.true.
  IsHotStart=.false.
else if (StartCondition == 2) then
  IsHGLInit=.false.
  IsHotStart=.true.
else
  write (*,*) 'StartCondition not specified correctly'
  IOError=.true.
  go to 999
end if

! scan files to get number of reaches and junctions

open (12,file='Reaches.inp',status='old',err=901)
Nreaches=0
do iskip=1,4
    read (12,*)
end do
do while(not(eof(12)))
    read (12,*)
    Nreaches=Nreaches+1
end do
close (12)

open (12,file='Junctions.inp',status='old',err=902)
Njuncs=0
do iskip=1,4
    read (12,*)
end do
do while(not(eof(12)))
    read (12,*)
    Njuncs=Njuncs+1
end do
close (12)

allocate(Reach(Nreaches))
allocate(Junc(Njuncs))

open (12,file='Reaches.inp',status='old',err=901)
do iskip=1,4
    read (12,*)
end do
NShapeTables=0
NNN=0
MaxNCells=0
do iread=1,Nreaches
    read (12,*) idummy,Reach(iread)%UpJunc,Reach(iread)%UpElev,Reach(iread)%DownJunc,Reach(iread)%DownElev,Reach(iread)%L,Reach(iread)%D,Reach(iread)%n,Reach(iread)%Xsec,Reach(iread)%Width,Reach(iread)%Kup,Reach(iread)%Kdown
    nc_temp=aint(Reach(iread)%L/(Reach(iread)%D*Ratio_dx_Dpipe))
    fragment=Reach(iread)%L/(Reach(iread)%D*Ratio_dx_Dpipe)-float(nc_temp)
   if (fragment < 0.5) then
     Reach(iread)%Ncells=nc_temp
   else
     Reach(iread)%Ncells=nc_temp+1
   end if
   Reach(iread)%Ncells=max(Reach(iread)%Ncells,MinCellCt)
   NNN=NNN+Reach(iread)%Ncells
   Reach(iread)%dX=Reach(iread)%L/Reach(iread)%Ncells
   if (Reach(iread)%Ncells > MaxNCells) MaxNCells=Reach(iread)%Ncells
   if (Reach(iread)%XSec == 3) NShapeTables=NShapeTables+1
   Reach(iread)%So=(Reach(iread)%UpElev-Reach(iread)%DownElev)/Reach(iread)%L
   select case(Reach(iread)%Xsec)
     case(1)
       Reach(iread)%Apipe=pi/4.0_wp*Reach(iread)%D**2
     case(2)
       Reach(iread)%Apipe=Reach(iread)%D*Reach(iread)%Width
   end select
end do
close (12)

! write cell count information for plotting output
open (13,file='ReachCellCount.csv',status='replace')
write (13,800)
StartPos=0
do jj=1,Nreaches
    StartPos=StartPos+Reach(jj)%NCells
    write (13,801) jj,Reach(jj)%NCells,StartPos
end do
close (13)
800 format ('ReachNo,NCells,StartPos')
801 format (i4,',',i8,',',i8)
    
if (NShapeTables > 0) call ReadCustomShapeTables

open (12,file='Junctions.inp',status='old',err=902)
do iskip=1,4
    read (12,*)
end do
SomeRatingCurves=.false.
SomeStorageCurves=.false.
SomeVariableInflows=.false.
SomeInflowControls=.false.
SomeSlideGates=.false.
SomeLevelSeries=.false.
Junc(:)%IsRatingCurve=.false.
Junc(:)%IsStorageCurve=.false.
Junc(:)%IsVariableInflow=.false.
Junc(:)%IsVented=.true.
Junc(:)%InflowMap=0
Junc(:)%StorageMap=0
Junc(:)%RatingCurveMap=0
Junc(:)%InflowControlMap=0
Junc(:)%RatingCurveType=0
Junc(:)%LevelTSMap=0
NGates=0
NLevelSer=0

do iread=1,Njuncs
    read (12,*) idummy,idummy,Junc(iread)%Elev,Junc(iread)%ShaftArea,Junc(iread)%OverflowElev,Junc(iread)%WeirLength,Junc(iread)%InflowInit,Junc(iread)%Option
    if (Junc(iread)%WeirLength < 0.0_wp) then
        SomeRatingCurves=.true.
        Junc(iread)%IsRatingCurve=.true.
    end if
    if (Junc(iread)%InflowInit < 0.0_wp) then
        SomeVariableInflows=.true.
        Junc(iread)%IsVariableInflow=.true.
    end if
    if (Junc(iread)%ShaftArea < 0.0_wp) then
        SomeStorageCurves=.true.
        Junc(iread)%IsStorageCurve=.true.
    end if
    if (Junc(iread)%Option == -1) then
        SomeSlideGates=.true.
        NGates=NGates+1
    end if
    if ((Junc(iread)%Option == -4) .or. (Junc(iread)%Option == -5)) SomeInflowControls=.true.
    if (Junc(iread)%Option == -6) then
        SomeLevelSeries=.true.
        NLevelSer=NLevelSer+1
    end if
end do

!> determine connectivity
do kt=1,Njuncs
  do jj=1,4
      Junc(kt)%ReachNo(jj)=0
      Junc(kt)%ReachCell(jj)=0
      Junc(kt)%ReachElev(jj)=0.0_wp
  end do
  jj=1
  do jdummy=1,Nreaches
    if (Reach(jdummy)%DownJunc == kt) then
	  Junc(kt)%ReachNo(jj)=jdummy
	  Junc(kt)%ReachCell(jj)=Reach(jdummy)%Ncells
	  Junc(kt)%ReachElev(jj)=Reach(jdummy)%DownElev
	  jj=jj+1
	end if
  end do
  do jdummy=1,Nreaches
    if (Reach(jdummy)%UpJunc == kt) then
	  Junc(kt)%ReachNo(jj)=jdummy
	  Junc(kt)%ReachCell(jj)=1
	  Junc(kt)%ReachElev(jj)=Reach(jdummy)%UpElev
	  jj=jj+1
	end if
  end do
  Junc(kt)%ReachCount=jj-1
end do

if (SomeVariableInflows) call ReadVariableInflows
if (SomeRatingCurves) call ReadRatingCurves
if (SomeStorageCurves) call ReadStorageCurves
if (SomeInflowControls) call ReadInflowControls
if (SomeSlideGates) call ReadSlideGates
if (SomeLevelSeries) call ReadLevelSeries

DiagOutUnits=0
QDiagOutUnits=0
open (unit=18,file='JunctionDiagnostics.inp',status='old',err=999)
kt=0
read (18,*)
do while(not(eof(18)))
  kt=kt+1
  read (18,*) JuncDiagIndex(kt),JuncDiagTime(kt)
  DiagOutUnits(kt)=800+kt
  QDiagOutUnits(kt)=1800+kt
  write (JuncD,'(i2)') JuncDiagIndex(kt)
  open (unit=DiagOutUnits(kt),form='formatted',file='Junction'//trim(JuncD)//'.dia',status='replace')
  open (unit=QDiagOutUnits(kt),form='formatted',file='JunctionQ'//trim(JuncD)//'.dia',status='replace')
end do
close(18)

WriteJuncDiagnostics=.true.
NJuncDiag=kt

go to 999

900 write(*,*) 'RunControl.inp not found'
go to 999
901 write(*,*) 'Reaches.inp not found'
go to 999
902 write(*,*) 'Junctions.inp not found'
go to 999

999 continue
    
end subroutine    