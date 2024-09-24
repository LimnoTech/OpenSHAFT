subroutine CalcAirFlux
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Calculates air displacement from reaches to junctions; called from ExecuteTimeLoop
!%=====================================================================

use GlobalVariables

implicit none

real(kind=8) :: OpenSum,VolSum,OldVolSum,TempFluxSum,TempFlux,TempA
integer :: m,ti,tj,mm,mk,kk,CurrentReach,NextReach,CurrentJunc,NextJunc,mr,mj,ms,mrt,mst,msb

if (dT == 0.0_wp) go to 800

IsSurcharged=.false.
TempFlux=0.0_wp

do m=1,NAirUnits
    OpenSum=0.0_wp
    do mm=1,AirUnitNEndJunc(m)
      ti=AirUnitIJ(m,mm,1)
	  tj=AirUnitIJ(m,mm,2)
      if (ti == 0) then
        OpenSum=OpenSum+VentArea(m,mm)
      else
        TempA=A(ti,tj)
        OpenSum=OpenSum+max((Reach(tj)%Apipe-TempA),0.0_wp)
      end if
    end do
    if (OpenSum <= 0.001_wp) then
      IsSurcharged(m)=.true.
	  cycle
    end if
    VolSum=0.0_wp
    OldVolSum=0.0_wp
    do mm=1,AirUnitNreaches(m)
      VolSum=VolSum+Reach(AirUnitReaches(m,mm))%Volume
      OldVolSum=OldVolSum+Reach(AirUnitReaches(m,mm))%VolOld
    end do
    do mm=1,AirUnitNEndJunc(m)
      ti=AirUnitIJ(m,mm,1)
	  tj=AirUnitIJ(m,mm,2)
      if (ti == 0) then
        AirUnitOpeningFlux(m,mm)=(VolSum-OldVolSum)*VentArea(m,mm)/OpenSum
      else
        TempA=A(ti,tj)
	    AirUnitOpeningFlux(m,mm)=(VolSum-OldVolSum)*max((Reach(tj)%Apipe-TempA),0.0_wp)/OpenSum
      end if
    end do
end do

do m=1,NAirUnits
  if (IsSurcharged(m)) then
    do mm=1,AirUnitNEndJunc(m)
	  AirUnitOpeningFlux(m,mm)=0.0_wp
	end do
  end if
end do

do k=1,Njuncs
  TempFluxSum=0.0_wp
  do kk=1,Junc(k)%AirFluxCount
    mk=Junc(k)%AirFluxMap(kk)
	TempFluxSum=TempFluxSum+AirUnitOpeningFlux(OpeningMap(mk,1),OpeningMap(mk,2))
  end do
!  JuncAirFlux(k)=TempFluxSum/dT
  Junc(k)%AvgFluxSum=Junc(k)%AvgFluxSum+TempFluxSum
end do

800 return

end subroutine