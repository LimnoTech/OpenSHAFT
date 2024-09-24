subroutine SaveOutputVariables
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% saves previous values of all output variables	
!%=====================================================================
use GlobalVariables
implicit none
! PRK 1/28/2015 this routine saves old values of all variables that are written to output, necessary for time interpolation at exact write times

do j=1,Nreaches
  do i=1,Reach(j)%NCells
    Qsave(i,j)=Q(i,j)
    Vold(i,j)=V(i,j)
    Aold(i,j)=A(i,j)
    yold(i,j)=y(i,j)
    hsold(i,j)=hs(i,j)
  end do
end do
Junc(:)%HeadOld=Junc(:)%Head
JuncInflowOld=Junc(:)%Inflow
JuncOutflowOld=Junc(:)%Outflow
!OfflineStoHeadOld=OfflineStoHead
!JuncStoFlowOld=JuncStoFlow
!JuncStoVolOld=JuncStoVol
!HydrographQOld=HydrographQ
!RejectedInflowOld=RejectedInflow
ReachVolTotOld=ReachVolTot
Junc(:)%VolOld=Junc(:)%Volume
JuncCumInfOld=Junc(:)%CumInf
JuncCumOutOld=Junc(:)%CumOut
JuncClosureOld=Junc(:)%Closure
AvgFluxSumOld=Junc(:)%AvgFluxSum
HchOld=Junc(:)%Hch
LchOld=Junc(:)%Lch
VolChamberOld=Junc(:)%VolChamber
QairOld=Junc(:)%Qair
!ReachVolSave=ReachVol

!> save output variables associated with reach objects
HaOld=Reach(:)%Ha
VolPocketOld=Reach(:)%VolPocket
QColumnOld=Reach(:)%QColumn
YUnderOld=Reach(:)%YUnder
QColDiffOld=Reach(:)%QColDiff
LpocketOld=Reach(:)%Lpocket
QPocketDownOld=Reach(:)%QPocketDown
QPocketUpOld=Reach(:)%QPocketUp

!if (NOrif > 0) GatePositionOld=GatePosition

return

end subroutine SaveOutputVariables