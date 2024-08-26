subroutine SlideGate

use GlobalVariables
use GlobalFunctions

implicit none

integer :: LeftCase,RightCase,iL,jL,iR,jR,kk
logical :: IsDownstreamQ,CaseIsHandled
real(wp) :: GateHead,VolCorr,Hdrop,ygate,Hcl
JunctionCase=0
LeftCase=0
RightCase=0
Hdrop=0.0_wp

if (dT == 0.0_wp) then
  go to 900
end if

iL=Junc(k)%ReachCell(1)
jL=Junc(k)%ReachNo(1)
iR=Junc(k)%ReachCell(2)
jR=Junc(k)%ReachNo(2)
kk=Junc(k)%GateMap
!Kgate=0.8_wp
Junc(k)%HeadOld=Junc(k)%Head
HcL=FindCriticalDepth(iL-1,jL)

! set gate position - this may be replaced with a subroutine call, for now time series is only option

!call FindCurrentGatePosition




if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,3i4,23f9.3)') T,JunctionCase,LeftCase,RightCase,Junc(k)%Head,y(iL-1,jL),y(iL,jL),y(iR,jR),y(iR+1,jR),Q(iL-1,jL), &
	  Q(iL,jL),Q(iR,jR),Q(iR+1,jR),V(iL-1,jL),V(iL,jL),V(iR,jR),V(iR+1,jR),c(iL-1,jL),c(iL,jL),c(iR,jR),c(iR+1,jR),hs(iL-1,jL),hs(iL,jL),hs(iR,jR),hs(iR+1,jR),ygate,HcL
	end if
  end do
end if

if (not(CaseIsHandled)) then
    write (*,*) 'case',JunctionCase, 'is not handled'
    stop
end if

900 continue

return

end subroutine