subroutine CalcContinuity

use GlobalVariables

implicit none

real(kind=8) :: TempSumR
integer :: i1,i2,i3,j1,j2,j3

if (dT == 0.0_wp) go to 601

Reach(:)%VolOld=Reach(:)%Volume
Junc(:)%VolOld=Junc(:)%Volume
ReachVolTot=0.0_wp
JuncVolTot=0.0_wp
JuncInfTot=0.0_wp
JuncOutTot=0.0_wp

do j=1,Nreaches
  TempSumR=0.0_wp
  do i=1,Reach(j)%Ncells
    TempSumR=TempSumR+Reach(j)%dX*A(i,j)
  end do
  Reach(j)%Volume=TempSumR
  ReachVolTot=ReachVolTot+Reach(j)%Volume
end do

do k=1,Njuncs
  Junc(k)%Volume=Junc(k)%VolOld+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)
  Junc(k)%CumInf=Junc(k)%CumInf+dT*Junc(k)%Inflow
  Junc(k)%CumOut=Junc(k)%CumOut+dT*Junc(k)%Outflow
  JuncVolTot=JuncVolTot+Junc(k)%Volume
  JuncInfTot=JuncInfTot+Junc(k)%CumInf
  JuncOutTot=JuncOutTot+Junc(k)%CumOut
end do

PercentVolError=100.0_wp*(ReachVolTot+JuncVolTot+JuncOutTot-(InitialSystemVol+JuncInfTot))/(InitialSystemVol+JuncInfTot)

if (nts > 1) then
  do k=1,Njuncs
    select case(Junc(k)%BCType)
	  case(1)
	    Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)-Q(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1)))
	  case(2)
	    Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1)))
      case(3)
	    Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1))-Q(Junc(k)%ReachCell(2),Junc(k)%ReachNo(2)))
      case(4)
	    i1=Junc(k)%ReachCell(1)
        i3=Junc(k)%ReachCell(2)
        i2=Junc(k)%ReachCell(3)
        j1=Junc(k)%ReachNo(1)
        j3=Junc(k)%ReachNo(2)
        j2=Junc(k)%ReachNo(3)
		Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(i1,j1)+Q(i3,j3)-Q(i2,j2))
      case(5)
	    i1=Junc(k)%ReachCell(1)
        i3=Junc(k)%ReachCell(2)
        i2=Junc(k)%ReachCell(3)
        j1=Junc(k)%ReachNo(1)
        j3=Junc(k)%ReachNo(2)
        j2=Junc(k)%ReachNo(3)
		Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(i1,j1)-Q(i3,j3)-Q(i2,j2))
      case(6)
	    Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1))-Q(Junc(k)%ReachCell(2),Junc(k)%ReachNo(2)))
	  case(7)
	    Junc(k)%Closure=Junc(k)%Closure+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1)))
      case default
        write (*,*) 'unsupported boundary condition type',Junc(k)%BCType
        stop
	end select
  end do
end if

601 return

end subroutine CalcContinuity