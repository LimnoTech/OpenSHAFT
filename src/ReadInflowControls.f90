subroutine ReadInflowControls
        
use GlobalVariables

implicit none

integer(i4) :: kt,kdummy,NGroups
logical :: HasControlGroups

!> scan inflow controls file for number of control locations
NInflowControl=0
HasControlGroups=.false.
open (unit=17,file='InflowControl.inp',status='old',err=903)
do iskip=1,4
  read (17,*)
end do
do while(not(eof(17)))
    read (17,*)
    NInflowControl=NInflowControl+1
end do
close (unit=17)

!> allocate arrays for inflow control parameters
allocate(Control(NInflowControl))
Control(:)%IsOpen=.true.
Control(:)%IsClosing=.false.
IsCentralControl=.false.
Control(:)%StartCloseTime=0.0_wp
Control(:)%EndCloseTime=0.0_wp
Control(:)%ClosingTime=0.0_wp
Control(:)%InflowFrac=1.0_wp
Control(:)%TriggerElev=0.0_wp
Control(:)%TriggerNode=0
Control(:)%MinInflow=0.0_wp
Control(:)%MinOpenFrac=0.0_wp
Control(:)%Fails=.false.
NTriggerSteps=0
TriggerStepCount=0
Control(:)%NControlSteps=0
Control(:)%TriggerGroup=0
MaxControlSteps=0
do kt=1,NInflowControl
  do kdummy=1,60
    Control(kt)%HSCNode(kdummy)=0.0_wp
  end do
end do
Control(:)%HSCNodeAvg=0.0_wp

open (unit=609,file='ControlActionsLog.dia',status='replace')

open (unit=17,file='InflowControl.inp',status='old',err=903)
do iskip=1,4
  read (17,*)
end do
do kt=1,NInflowControl
  read (17,*) Control(kt)%InflowLoc,Control(kt)%TriggerNode,Control(kt)%TriggerElev,Control(kt)%ClosingTime,Control(kt)%MinInflow,Control(kt)%FailFlag,Control(kt)%NControlSteps,Control(kt)%MinOpenFrac,Control(kt)%TriggerGroup
  if (Control(kt)%FailFlag == 1) Control(kt)%Fails=.true.
  if (Control(kt)%TriggerGroup > 0) HasControlGroups=.true.
  MaxControlSteps=max(MaxControlSteps,Control(kt)%NControlSteps)
  do kdummy=1,Njuncs
    if (Control(kt)%InflowLoc == kdummy) then
      Junc(kdummy)%InflowControlMap=kt
      exit
    end if
  end do
end do
close(unit=17)

!> check to see if inflow control data exists for each junction that expects it

do kdummy=1,NJuncs
  if ((Junc(kdummy)%Option == -4) .or. (Junc(kdummy)%Option == -5)) then
      if (Junc(kdummy)%InflowControlMap == 0) then
      write (*,*) 'inflow control data for junction',kdummy,' is missing'
	  IOError=.true.
      go to 999
      end if
  end if
end do

! optionally read control group information if needed
if (HasControlGroups) then
    open (unit=17,file='ControlGroupData.inp',status='old',err=904)
    do iskip=1,4
      read (17,*)
    end do
    read (17,*) NGroups
    allocate (CCGroup(NGroups))
    do kdummy=1,NGroups
        read (17,*) CCGroup(kdummy)%NJuncs
        read (17,*) (CCGroup(kdummy)%JuncID(kt),kt=1,CCGroup(kdummy)%NJuncs)
        read (17,*) (CCGroup(kdummy)%TriggerH(kt),kt=1,CCGroup(kdummy)%NJuncs)
    end do
end if

go to 999

903 write(*,*) 'InflowControls.inp not found'
IOError=.true.

go to 999

904 write(*,*) 'ControlGroupData.inp not found'
IOError=.true.

999 continue
    
return

end subroutine
