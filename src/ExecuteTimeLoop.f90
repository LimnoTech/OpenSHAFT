subroutine ExecuteTimeLoop
    
use GlobalVariables
use GlobalFunctions
use omp_lib

implicit none

integer(i4) :: kdummy,jt,tid

T=0.0_wp
dT=MaxdT

do while (T < Tsim)
    
    call SaveOutputVariables

!$omp parallel private(maxcroechange,ceil,jt,tid)
    tid=omp_get_thread_num()+1    
    do jt=1,WorkByThread(tid)%jcount
      j=WorkByThread(tid)%jlist(jt)
!    do j=1,Nreaches
      call CalcFluxSources
      call CalcCellFluxes
      if (TrapAirPockets) call ScanForPockets
      call UpdateAQ
      call TPACalc
    end do
!$omp end parallel 

    if (TrapAirPockets) then
      do j=1,Nreaches
        if (Reach(j)%PocketFound) then
          write (223,'(f11.4,3i5,f9.4)') T,j,Reach(j)%iFRT1,Reach(j)%iFRT2,Reach(j)%VolPocket
        end if
        if ((Reach(j)%HasPocket) .and. (Reach(j)%InitializePocket)) then
          Reach(j)%InitializePocket=.false.
          do i=Reach(j)%iFRT1,Reach(j)%iFRT2
            write (227,'(i3,i6,3f10.4)') j,i,hs(i,j),Reach(j)%Apipe-A(i,j),Q(i,j)
          end do
          write (222,'(f11.4,3i5,f9.4)') Reach(j)%PocketStart,j,Reach(j)%iFRT1,Reach(j)%iFRT2,Reach(j)%VolPocket
          write (85,805) T,j,Reach(j)%Ha,Reach(j)%QColumn,Reach(j)%VolPocket,Reach(j)%YUnder,Reach(j)%QColDiff
        end if
      end do
    end if
    
    do k=1,NJuncs
        call FindCurrentShaftArea
        call FindOutflow
        call FindInflow
        if (Junc(k)%Option == -4) then
            if (not(Control(Junc(k)%InflowControlMap)%IsOpen)) Junc(k)%Outflow=0.0_wp
        end if
    end do
    
    do k=1,NJuncs
      do kdummy=1,Njuncs
        if (kdummy == k) cycle
        if (Junc(kdummy)%JuncPartner == k) Junc(k)%Inflow=Junc(k)%Inflow+Junc(kdummy)%Outflow
      end do
    end do
!=======================================================
!
!  the junction loop calls appropriate boundary condition routines
!  based on the value of Junc(k)%BCType.  currently supported options are:
!
!    1   single upstream dropshaft
!    2   single downstream dropshaft
!    3   two-way dropshaft
!    4   three-way dropshaft
!    5   three-way split
!    6   slope change
!    7   level time series
!
!=======================================================
    do k=1,Njuncs
      select case (Junc(k)%BCType)
	    case (1)
          call SingleDropshaftUp 
        case (2)
	      call SingleDropshaftDown
	    case (3)
	      call TwoWayDropshaft
	    case (4)
    	  call ThreeWayDropshaft
        case (5)
          call ThreeWaySplit
        case (6)
          call InnerSlopeChange
        case (7)
          call LevelSeriesShaft
        case default
          write (*,*) 'unsupported boundary condition type',Junc(k)%BCType,' for junction',k
          stop
      end select
    
    end do
    
    call CalcContinuity
    
    if (CalcAirExhaust) call CalcAirFlux
    
    if (WriteFilesOpen) then
      if (T >= WriteTimes(IWrite)) then
        call WriteOutput
        IWrite=IWrite+1
        LastWriteNSteps=nts
        if (IWrite > NWrites) WriteFilesOpen=.false.
      end if
    end if
    
    if (SomeInflowControls) then
      if (T >= ControlEvalTimes(IControlEval)) then
        call InflowControlTrigger
        IControlEval=IControlEval+1
        if (IControlEval > NControlEvalTimes) SomeInflowControls=.false.
      end if
    end if
    
    do j=1,Nreaches
      do i=1,Reach(j)%Ncells
        if (not(PocketZone(i,j))) then
	      if (hs(i,j) < MaxVacuum(i,j)) MaxVacuum(i,j)=hs(i,j)
        else
          if (T > Reach(j)%PocketStart+5.0_wp) then
            if (hs(i,j) < MaxVacuum(i,j)) MaxVacuum(i,j)=hs(i,j)
          end if
        end if
	    if (y(i,j)+z(i,j) > MaxHGL(i,j)) MaxHGL(i,j)=y(i,j)+z(i,j)
      end do
    end do

    call UpdateTime

    write (*,334) T,dT,PercentVolError,Junc(ShowJunc(1))%Head,Junc(ShowJunc(2))%Head,Junc(ShowJunc(3))%Head

    T=T+dT
    nts=nts+1

end do    
    
return

334 format ('+',f10.3,es13.4,f10.3,3f9.1)
805 format (es13.6,',',i3,5(',',es16.8))

end subroutine ExecuteTimeLoop