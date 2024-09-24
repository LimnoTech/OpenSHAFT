subroutine InflowControlTrigger
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Determine whether to initiate inflow control based on junction depth	
!%=====================================================================

use GlobalVariables

implicit none

real(wp) :: HSCNodeSum
integer(i4) :: kk

do k=1,NJuncs
  if ((Junc(k)%Option == -4) .or. (Junc(k)%Option == -5)) then
    kk=Junc(k)%InflowControlMap
    if (IControlEval < Control(kk)%NControlSteps) then
      Control(kk)%HSCNode(IControlEval)=Junc(Control(kk)%TriggerNode)%Head
    else
      Control(kk)%HSCNode(Control(kk)%NControlSteps)=Junc(Control(kk)%TriggerNode)%Head
      HSCNodeSum=0.0_wp
      do i=1,Control(kk)%NControlSteps
        HSCNodeSum=HSCNodeSum+Control(kk)%HSCNode(i)
      end do
      Control(kk)%HSCNodeAvg=HSCNodeSum/float(Control(kk)%NControlSteps)
      do i=1,Control(kk)%NControlSteps-1
        Control(kk)%HSCNode(i)=Control(kk)%HSCNode(i+1)
      end do
    end if
  end if
end do

do k=1,Njuncs
  if ((Junc(k)%Option == -4) .or. (Junc(k)%Option == -5)) then
    kk=Junc(k)%InflowControlMap
    if (Control(kk)%Fails) cycle
    if (Control(kk)%IsOpen) then
      if ((Control(kk)%HSCNodeAvg > Control(kk)%TriggerElev) .and. (not(Control(kk)%IsClosing))) then
        Control(kk)%IsClosing=.true.
        Control(kk)%StartCloseTime=T
        Control(kk)%EndCloseTime=T+Control(kk)%ClosingTime
        write (609,*) 'control at junction',k,'began closing at',T
      end if
    end if
  end if
end do

return

end subroutine InflowControlTrigger