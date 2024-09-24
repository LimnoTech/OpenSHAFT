subroutine WriteOutput
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Interpolates output to specific write intervals and then writes said output	
!%=====================================================================

use GlobalVariables

implicit none

real(wp) :: deltaT,TempVelSum,MaxVel,MinVel,AvgVel,TempDepthSum,AvgDepth,Goof,AvgAirFlux
integer :: kk,m,ldump
! perform interpolation on all variables to be written
deltaT=(WriteTimes(IWrite)-(T-dT))/dT

do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    QWrite(i,j)=Qsave(i,j)*(1.0_wp-deltaT) + Q(i,j)*deltaT
    VWrite(i,j)=Vold(i,j)*(1.0_wp-deltaT) + V(i,j)*deltaT
    AWrite(i,j)=Aold(i,j)*(1.0_wp-deltaT) + A(i,j)*deltaT
    yWrite(i,j)=yold(i,j)*(1.0_wp-deltaT) + y(i,j)*deltaT
    hsWrite(i,j)=hsold(i,j)*(1.0_wp-deltaT) + hs(i,j)*deltaT
  end do
end do

JuncHeadWrite=Junc(:)%HeadOld*(1.0_wp-deltaT) + Junc(:)%Head*deltaT
JuncInflowWrite=JuncInflowOld*(1.0_wp-deltaT) + Junc(:)%Inflow*deltaT
JuncOutflowWrite=JuncOutflowOld*(1.0_wp-deltaT) + Junc(:)%Outflow*deltaT
!OfflineStoHeadWrite=OfflineStoHeadOld*(1.0_wp-deltaT) + OfflineStoHead*deltaT
!JuncStoFlowWrite=JuncStoFlowOld*(1.0_wp-deltaT) + JuncStoFlow*deltaT
!JuncStoVolWrite=JuncStoVolOld*(1.0_wp-deltaT) + JuncStoVol*deltaT
!HydrographQWrite=HydrographQOld*(1.0_wp-deltaT) +HydrographQ*deltaT
!RejectedInflowWrite=RejectedInflowOld*(1.0_wp-deltaT) + RejectedInflow*deltaT
ReachVolTotWrite=ReachVolTotOld*(1.0_wp-deltaT) + ReachVolTot*deltaT
JuncVolWrite=Junc(:)%VolOld*(1.0_wp-deltaT) + Junc(:)%Volume*deltaT
JuncInfWrite=JuncCumInfOld*(1.0_wp-deltaT) + Junc(:)%CumInf*deltaT
JuncOutWrite=JuncCumOutOld*(1.0_wp-deltaT) + Junc(:)%CumOut*deltaT
JuncClosureWrite=JuncClosureOld*(1.0_wp-deltaT) + Junc(:)%Closure*deltaT
AvgFluxSumWrite=AvgFluxSumOld*(1.0_wp-deltaT) + Junc(:)%AvgFluxSum*deltaT
!ReachVolWrite=ReachVolSave*(1.0_wp-deltaT) + ReachVol*deltaT

HaWrite=HaOld*(1.0_wp-deltaT) + Reach(:)%Ha*deltaT
VolPocketWrite=VolPocketOld*(1.0_wp-deltaT) + Reach(:)%VolPocket*deltaT
QColumnWrite=QColumnOld*(1.0_wp-deltaT) + Reach(:)%QColumn*deltaT
YUnderWrite=YUnderOld*(1.0_wp-deltaT) + Reach(:)%YUnder*deltaT
QColDiffWrite=QColDiffOld*(1.0_wp-deltaT) + Reach(:)%QColDiff*deltaT
LpocketWrite=LpocketOld*(1.0_wp-deltaT) + Reach(:)%Lpocket*deltaT
QPocketDownWrite=QPocketDownOld*(1.0_wp-deltaT) + Reach(:)%QPocketDown*deltaT
QPocketUpWrite=QPocketUpOld*(1.0_wp-deltaT) + Reach(:)%QPocketUp*deltaT

HchWrite=HchOld*(1.0_wp-deltaT) + Junc(:)%Hch*deltaT
LchWrite=LchOld*(1.0_wp-deltaT) + Junc(:)%Lch*deltaT
VolChamberWrite=VolChamberOld*(1.0_wp-deltaT) + Junc(:)%VolChamber*deltaT
QairWrite=QairOld*(1.0_wp-deltaT) + Junc(:)%Qair*deltaT

!if (NOrif > 0) GatePositionWrite=GatePositionOld*(1.0_wp-deltaT) + GatePosition*deltaT

! write junction output

do k=1,Njuncs
  write (54,504) WriteTimes(IWrite),k,Junc(k)%BCType,JuncHeadWrite(k),JuncInflowWrite(k),QWrite(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1)),QWrite(Junc(k)%ReachCell(2),Junc(k)%ReachNo(2)), &
                 QWrite(Junc(k)%ReachCell(3),Junc(k)%ReachNo(3)),QWrite(Junc(k)%ReachCell(4),Junc(k)%ReachNo(4)),JuncOutflowWrite(k),Junc(k)%CurrentArea
!  if (IsOffStoCurve(k)) write (59,506) WriteTimes(IWrite),k,JuncHeadWrite(k),OfflineStoHeadWrite(OffStoCurveMap(k)),JuncStoFlowWrite(k),JuncStoVolWrite(k)
!  if (IsInflowLimit(k)) write (61,506) WriteTimes(IWrite),k,JuncHeadWrite(k),HydrographQWrite(k),RejectedInflowWrite(k),JuncInflowWrite(k)
  if (Junc(k)%BCType == 4) write (720,702) WriteTimes(IWrite),k,(ThreeWayCases(k,m),m=1,20)
  if (Junc(k)%Option == -8) write (53,506) WriteTimes(IWrite),k,HchWrite(k),LchWrite(k),VolChamberWrite(k),QairWrite(k)
!  if (JuncOption(k) == -11) write (953,804) WriteTimes(IWrite),k,(PumpCurrentTau(kk),kk=PumpCurveMap(k,1),PumpCurveMap(k,NPumps(k)))
end do

! calculate average air flux and write
if (CalcAirExhaust) then
  do k=1,Njuncs
    if (IWrite == 1) then
      AvgAirFlux=AvgFluxSumWrite(k)/StartWrite
    else
      AvgAirFlux=AvgFluxSumWrite(k)/Twrite
    end if
    if (Junc(k)%IsVented) write (60,610) WriteTimes(IWrite),k,AvgAirFlux
    Junc(k)%AvgFluxSum=Junc(k)%AvgFluxSum-AvgFluxSumWrite(k)
  end do
end if

! calculate reach-average max and min velocity
do j=1,Nreaches
!  TempVelSum=0.0_wp
!  TempDepthSum=0.0_wp
!  MaxVel=-9999.0_wp
!  MinVel=9999.0_wp
!  do i=1,Ncell(j)
!    TempVelSum=TempVelSum+VWrite(i,j)
!    TempDepthSum=TempDepthSum+yWrite(i,j)
!    if (VWrite(i,j) > MaxVel) then
!      MaxVel=VWrite(i,j)
!    end if
!    if (VWrite(i,j) < MinVel) then
!      MinVel=VWrite(i,j)
!    end if
!  end do
!  AvgVel=TempVelSum/Ncell(j)
!  AvgDepth=TempDepthSum/Ncell(j)
!  write (65,605) WriteTimes(IWrite),j,VWrite(MiddleCell(j),j),AvgVel,MinVel,MaxVel,AvgDepth
  if (Reach(j)%HasPocket) then
    if (AirPocketSwitch == 2) then
      write (85,806) WriteTimes(IWrite),j,HaWrite(j),QColumnWrite(j),VolPocketWrite(j),YUnderWrite(j),QColDiffWrite(j),LpocketWrite(j),QPocketDownWrite(j),QPocketUpWrite(j),Reach(j)%iFRT1,Reach(j)%iFRT2
    else
      write (85,805) WriteTimes(IWrite),j,HaWrite(j),QColumnWrite(j),VolPocketWrite(j),YUnderWrite(j),QColDiffWrite(j)
    end if
  end if
end do

! write gate position if there are any
!if (NOrif > 0) then
!    do k=1,NOrif
!        write (88,609) WriteTimes(IWrite),JuncOrifID(k),GatePositionWrite(k)
!    end do
!end if

! binary dumps
ldump=1
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    TempDumpBin(ldump)=yWrite(i,j)+z(i,j)
	ldump=ldump+1
  end do 
end do
write (72,rec=IWrite) TempDumpBin

ldump=1
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    TempDumpBin(ldump)=VWrite(i,j)
	ldump=ldump+1
  end do 
end do
write (75,rec=IWrite) TempDumpBin

ldump=1
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    Goof=AWrite(i,j)
    TempDumpBin(ldump)=AWrite(i,j)
	ldump=ldump+1
  end do 
end do
write (74,rec=IWrite) TempDumpBin

ldump=1
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    TempDumpBin(ldump)=QWrite(i,j)
	ldump=ldump+1
  end do 
end do
write (73,rec=IWrite) TempDumpBin

ldump=1
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    TempDumpBin(ldump)=hsWrite(i,j)
	ldump=ldump+1
  end do 
end do
write (76,rec=IWrite) TempDumpBin

! write continuity checking variables
do k=1,Njuncs
  write (52,609) WriteTimes(IWrite),k,ReachVolTotWrite,JuncVolWrite(k),JuncInfWrite(k),JuncOutWrite(k),JuncClosureWrite(k)
end do

! write reach volumes for air calcs
!do j=1,Nreaches
!  write (152,609) WriteTimes(IWrite),j,ReachVolWrite(j)
!end do

! if last write, write restart file
if (IWrite == NWrites) then
  open (unit=444,file='hotstart.inp',status='replace')

  do j=1,Nreaches
    do i=1,Reach(j)%Ncells
      write (444,*) z(i,j),yWrite(i,j),VWrite(i,j),AWrite(i,j),QWrite(i,j),MaxVacuum(i,j),MaxHGL(i,j)
    end do
  end do
  do k=1,Njuncs
    write (444,*) JuncHeadWrite(k)
  end do
  if (TrapAirPockets) then
    do j=1,Nreaches
      if (Reach(j)%HasPocket) then
        write (444,*) j,Reach(j)%iFRT1,Reach(j)%iFRT2,Reach(j)%LPocket,Reach(j)%LColumn,HaWrite(j),QColumnWrite(j),VolPocketWrite(j)
      end if
    end do
  end if
  close (444)
end if

return

504 format (es13.6,',',i3,',',i3,12(',',es12.5))
506 format (es13.6,',',i3,4(',',es12.5))
609 format (es13.6,',',i3,7(',',es16.8))
610 format (es13.6,',',i3,3(',',es16.8))
605 format (es13.6,',',i3,5(',',es16.8))
805 format (es13.6,',',i3,6(',',es16.8))
806 format (es13.6,',',i3,8(',',es16.8),2(',',i5))
702 format (f14.4,i3,20i8)
804 format (es13.6,',',i3,20(',',es16.8))

end subroutine WriteOutput