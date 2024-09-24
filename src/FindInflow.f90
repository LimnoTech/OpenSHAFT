subroutine FindInflow
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Looks up/interpolates time series inflows; adjusts to simulate gate closure (optional)	
!%=====================================================================

use GlobalVariables

implicit none

real(wp) :: JuncHeadTemp,frac

if (dT == 0.0_wp) then
  go to 99
end if

if ((Junc(k)%Option == -4) .or. (Junc(k)%Option == -5)) call SetInflowControl

if (not(Junc(k)%IsVariableInflow)) then
    Junc(k)%Inflow=Junc(k)%InflowInit
else
  if (InflowTime(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap) == T) then
    Junc(k)%Inflow=Inflow(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap)
  else if (InflowTime(TimeSeriesIndex(Junc(k)%InflowMap)+1,Junc(k)%InflowMap) <= T) then
    TimeSeriesIndex(Junc(k)%InflowMap)=TimeSeriesIndex(Junc(k)%InflowMap)+1 
    frac=(T-InflowTime(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap))/ &
	(InflowTime(TimeSeriesIndex(Junc(k)%InflowMap)+1,Junc(k)%InflowMap)-InflowTime(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap))
	Junc(k)%Inflow=Inflow(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap) + &
	frac*(Inflow(TimeSeriesIndex(Junc(k)%InflowMap)+1,Junc(k)%InflowMap)-Inflow(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap))
  else if (InflowTime(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap) < T) then
    frac=(T-InflowTime(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap))/ &
	(InflowTime(TimeSeriesIndex(Junc(k)%InflowMap)+1,Junc(k)%InflowMap)-InflowTime(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap))
	Junc(k)%Inflow=Inflow(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap) + &
	frac*(Inflow(TimeSeriesIndex(Junc(k)%InflowMap)+1,Junc(k)%InflowMap)-Inflow(TimeSeriesIndex(Junc(k)%InflowMap),Junc(k)%InflowMap))
  else
    write (*,*) 'inflow interpolation problem at junction',k
	interpoproblem=.true.
  end if
end if    

!> adjust inflow according to inflow control settings if necessary
if ((Junc(k)%Option == -4) .or. (Junc(k)%Option == -5)) then
  if (Control(Junc(k)%InflowControlMap)%IsOpen .and. (not(Control(Junc(k)%InflowControlMap)%IsClosing))) then
     Junc(k)%Inflow=Control(Junc(k)%InflowControlMap)%InflowFrac*Junc(k)%Inflow
  else
     Junc(k)%Inflow=min(max(Control(Junc(k)%InflowControlMap)%InflowFrac*Junc(k)%Inflow,Control(Junc(k)%InflowControlMap)%MinInflow),Junc(k)%Inflow)
  end if
end if

99 return

end subroutine FindInflow