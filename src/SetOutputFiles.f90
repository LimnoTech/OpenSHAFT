subroutine SetOutputFiles
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Opens output files	
!%=====================================================================

use GlobalVariables

implicit none

integer :: m

!open (unit=71,form='unformatted',file=trim(RunID)//'.dump.bin',status='replace')
open (unit=72,file=trim(RunID)//'.HDump.bin',form='unformatted',access='direct',recl=2*NNN,status='replace')
open (unit=73,file=trim(RunID)//'.QDump.bin',form='unformatted',access='direct',recl=2*NNN,status='replace')
open (unit=74,file=trim(RunID)//'.ADump.bin',form='unformatted',access='direct',recl=2*NNN,status='replace')
open (unit=75,file=trim(RunID)//'.VDump.bin',form='unformatted',access='direct',recl=2*NNN,status='replace')
open (unit=76,file=trim(RunID)//'.hsDump.bin',form='unformatted',access='direct',recl=2*NNN,status='replace')
allocate (TempDumpBin(NNN))
TempDumpBin=0.0_wp

open (unit=52,form='formatted',file=trim(RunID)//'.Continuity.csv',status='replace')
open (unit=54,form='formatted',file=trim(RunID)//'.JunctionVariables.csv',status='replace')
!open (unit=65,form='formatted',file=trim(RunID)//'.ReachVelocity.csv',status='replace')
!open (unit=152,form='formatted',file=trim(RunID)//'.ReachVolume.csv',status='replace')
open (unit=720,form='formatted',file='ThreeWayCases.dia',status='replace')

if (CalcAirExhaust) open (unit=60,form='formatted',file=trim(RunID)//'.AirExhaust.csv',status='replace')

if (TrapAirPockets) open (unit=85,form='formatted',file=trim(RunID)//'.PocketOutput.csv',status='replace')
if (TrapAirPockets) open (unit=222,form='formatted',file='PocketFormation.dia',status='replace')
if (TrapAirPockets) open (unit=223,form='formatted',file='PocketTrapping.dia',status='replace')
if (TrapAirPockets) open (unit=227,form='formatted',file='PocketFlows.dia',status='replace')

write (52,512)
write (54,514)
!write (152,1512)
!write (65,525)
if (CalcAirExhaust) write (60,520)
if (AirPocketSwitch == 1) write (85,585)
if (AirPocketSwitch == 2) write (85,586)

WriteFilesOpen=.true.

!> create vector of output writing times
NWrites=aint((StopWrite-StartWrite)/Twrite)+1
allocate (WriteTimes(NWrites))
do IWrite=1,NWrites
  WriteTimes(IWrite)=StartWrite+float(IWrite-1)*Twrite
end do
IWrite=1

!> create vector of control evaluation times
NControlEvalTimes=aint(Tsim/Twrite)
allocate (ControlEvalTimes(NControlEvalTimes))
do m=1,NControlEvalTimes
    ControlEvalTimes(m)=float(m)*Twrite
end do
IControlEval=1

return

512 format ('T(s),JuncID,ReachVolTot,JuncVol,JuncInf,JuncOut,FlowClosure')
513 format ('T(s),JuncID,Hch(m),Lch(m),Vol(m3),Qair(m3/s)')
514 format ('T(s),JuncID,JuncType,H(m),Qjun(m3/s),Q_1(m3/s),Q_2(m3/s),Q_3(m3/s),Q_4(m3/s),Qoverflow(m3/s),ShaftArea(m2)')
519 format ('T(s),Junction No,H(m),StoH(m),StoFlow(m3/s),StoVol(m3)')
520 format ('T(s),JuncID,AvgAirFlux(m3/s)')
525 format ('T(s),ReachID,MidV(m/s),AvgV(m/s),MinV(m/s),Max(m/s),Y(m)')
585 format ('T(s),ReachID,Ha(m),QColumn(m3/s),VolPocket(m3),YUnder(m),QColDiff(m3/s)')
586 format ('T(s),ReachID,Ha(m),QColumn(m3/s),VolPocket(m3),YUnder(m),QColDiff(m3/s),Lpocket(m),QPocketDown(m3/s),QPocketUp(m3/s),iFRT1,iFRT2')
688 format ('T(s),JuncID,GatePos(--)')
788 format ('T(s),JuncID,PumpTau(--)')
1512 format ('T(s),ReachNo,ReachVol(m3)')

end subroutine SetOutputFiles