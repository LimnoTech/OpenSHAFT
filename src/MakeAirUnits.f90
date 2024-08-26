subroutine MakeAirUnits
    
use GlobalVariables

implicit none

integer(i4) :: m,mm,mk,OpeningCount
integer(i4), allocatable, dimension(:) :: TempReachID
real(wp), allocatable, dimension(:) :: EndCode

  open (22,file='AirFluxInput.inp',status='old')
  do iskip=1,2
    read (22,*)
  end do
  read (22,*) NAirUnits
  read (22,*) OpeningCount
  allocate(AirUnitNreaches(NAirUnits),AirUnitReaches(NAirUnits,NReaches),AirUnitNEndJunc(NAirUnits),AirUnitEndJunc(NAirUnits,OpeningCount),AirUnitOpeningFlux(NAirUnits,OpeningCount))
  allocate(OpeningID(NAirUnits,OpeningCount),OpeningMap(NReaches*OpeningCount,2),AirUnitIJ(NAirUnits,OpeningCount,2))
  allocate(TempReachID(OpeningCount),EndCode(OpeningCount),IsSurcharged(NAirUnits))
  allocate(VentArea(NAirUnits,OpeningCount))
  AirUnitOpeningFlux=0.0_wp
  EndCode=0.0_wp
  VentArea=0.0_wp
  do iskip=1,2
    read (22,*)
  end do
  do m=1,NAirUnits
    read (22,*) AirUnitNreaches(m),(AirUnitReaches(m,mm),mm=1,AirUnitNreaches(m))
  end do
  do iskip=1,2
    read (22,*)
  end do
  do m=1,NAirUnits
    read (22,*) AirUnitNEndJunc(m),(AirUnitEndJunc(m,mm),mm=1,AirUnitNEndJunc(m))
  end do
  do iskip=1,2
    read (22,*)
  end do
  do m=1,NAirUnits
    read (22,*) (TempReachID(mm),mm=1,AirUnitNEndJunc(m))
    read (22,*) (EndCode(mm),mm=1,AirUnitNEndJunc(m))
    do mm=1,AirUnitNEndJunc(m)
      if (EndCode(mm) == 1.0_wp) then
        AirUnitIJ(m,mm,1)=1
        AirUnitIJ(m,mm,2)=TempReachID(mm)
      end if
      if (EndCode(mm) == 0.0_wp) then
        AirUnitIJ(m,mm,1)=Reach(TempReachID(mm))%Ncells
        AirUnitIJ(m,mm,2)=TempReachID(mm)
      end if
      if (EndCode(mm) < 0.0_wp) then
        AirUnitIJ(m,mm,1)=0
        AirUnitIJ(m,mm,2)=0
        VentArea(m,mm)=abs(EndCode(mm))
      end if
    end do
  end do
  mk=1
  do m=1,NAirUnits
    do mm=1,AirUnitNEndJunc(m)
      OpeningID(m,mm)=mk
      OpeningMap(mk,1)=m
      OpeningMap(mk,2)=mm
      mk=mk+1
    end do
  end do

Junc(:)%AirFluxCount=0
do k=1,Njuncs
  do mk=1,4
    Junc(k)%AirFluxMap(mk)=0
  end do
  do m=1,NAirUnits
    do mm=1,AirUnitNEndJunc(m)
	  if (AirUnitEndJunc(m,mm) == k) then
	    Junc(k)%AirFluxCount=Junc(k)%AirFluxCount+1
		Junc(k)%AirFluxMap(Junc(k)%AirFluxCount)=OpeningID(m,mm)
        Junc(k)%IsVented=.true.
	  end if
	end do
  end do
end do

return

end subroutine