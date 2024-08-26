subroutine ScanForPockets

use GlobalVariables
use GlobalFunctions

implicit none

real(wp), save :: TempPocketVol,APocket
!$omp threadprivate (TempPocketVol,APocket)

Reach(j)%PocketFound=.false.

!do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    Aold(i,j)=A(i,j)
	Qold(i,j)=Q(i,j)
	PocketZone(i,j)=.false.
  end do
!end do

! scan reach for air pockets
!do j=1,Nreaches
  if ((TrapAirPockets) .and. (not(Reach(j)%HasPocket)) .and. (T > Reach(j)%PocketReleaseT+30.0_wp)) then
    Reach(j)%iFRT1=0
    Reach(j)%iFRT2=0
    do i=2*PocketBuffer+1,Reach(j)%Ncells-2*PocketBuffer-1
      if ((hs(i,j) /= 0) .and. (hs(i+1,j) == 0) .and. (Reach(j)%iFRT1 == 0)) Reach(j)%iFRT1=i !PRK 1/2/2013 added third clause to prevent reassignment of iFRT1 further downstream
      if ((hs(i,j) == 0) .and. (hs(i+1,j) /= 0) .and. (Reach(j)%iFRT1 > 0)) Reach(j)%iFRT2=i  !PRK 1/2/2013 third clause says do not close a pocket unless one has been identified
    end do
    if (Reach(j)%iFRT1*Reach(j)%iFRT2 > 0) then
      TempPocketVol=0.0_wp
      do i=Reach(j)%iFRT1,Reach(j)%iFRT2
        TempPocketVol=TempPocketVol+Reach(j)%dX*(Reach(j)%Apipe-A(i,j))
      end do
      Reach(j)%PocketFound=.true.
      Reach(j)%VolPocket=TempPocketVol
      if (Reach(j)%VolPocket > Reach(j)%MinPocketVol) then
        Reach(j)%HasPocket=.true.
        Reach(j)%InitializePocket=.true.
        Reach(j)%PocketStart=T
        Reach(j)%VolPocket=TempPocketVol
        Reach(j)%VolPocketInit=Reach(j)%VolPocket
        Reach(j)%Ha=0.0_wp
        Reach(j)%QColumn=0.0_wp
        do i=Reach(j)%iFRT1,Reach(j)%iFRT2
          Reach(j)%Ha=Reach(j)%Ha+hs(i,j)
        end do
        Reach(j)%Ha=Reach(j)%Ha/(Reach(j)%iFRT2-Reach(j)%iFRT1+1)
        Reach(j)%Ha=Reach(j)%Ha+atm
        Reach(j)%LPocket=Reach(j)%dX*abs(Reach(j)%iFRT2-Reach(j)%iFRT1+1)
        Reach(j)%LColumn=Reach(j)%dX*abs(Reach(j)%iFRT2-Reach(j)%iFRT1+2*PocketBuffer)
        do i=Reach(j)%iFRT1-PocketBuffer,Reach(j)%iFRT2+PocketBuffer
          Reach(j)%QColumn=Reach(j)%QColumn+Q(i,j)
        end do
        Reach(j)%QColumn=Reach(j)%QColumn/(Reach(j)%iFRT2-Reach(j)%iFRT1+2*PocketBuffer+1)
        Reach(j)%QColDiff=Q(Reach(j)%iFRT1-PocketBuffer-1,j)-Q(Reach(j)%iFRT2+PocketBuffer+1,j)
        Reach(j)%xDown=Reach(j)%iFRT2*Reach(j)%dX
        Reach(j)%xUp=(Reach(j)%iFRT1-1)*Reach(j)%dx
        APocket=Reach(j)%VolPocket/Reach(j)%LPocket
        select case(Reach(j)%XSec)
          case(1)
            tta=FindTheta(Reach(j)%Apipe-APocket,Reach(j)%D)
            Reach(j)%YUnder=0.5*Reach(j)%D*(1-cos(0.5*tta))
          case(2)
            Reach(j)%YUnder=(Reach(j)%Apipe-APocket)/Reach(j)%Width
          case(3)
            call GetLookupValue(Reach(j)%Apipe-APocket,Reach(j)%YUnder,Reach(j)%ReachShapeMap,2,1)
        end select
      end if
    end if
  end if
!end do
! end scan for air pockets

! PRK 4/22/2015 adding check to see if pocket can be released, based on either end returning to free surface flow
!do j=1,Nreaches
  if ((Reach(j)%HasPocket) .and. (T > Reach(j)%PocketStart+Reach(j)%ReleaseTime)) then
    if ((hs(Reach(j)%iFRT1-PocketBuffer-1,j) == 0.0_wp) .or. (hs(Reach(j)%iFRT2+PocketBuffer+1,j) == 0.0_wp)) then
      Reach(j)%HasPocket=.false.
      Reach(j)%PocketReleaseT=T
    end if
  end if
!end do
! PRK 4/22/2015 end

! optionally update pocket length unless this is the first call
!if ((Reach(j)%HasPocket) .and. (not(Reach(j)%InitializePocket)) .and. (AirPocketSwitch == 2)) call AirPocketDispCalc


return


end subroutine ScanForPockets