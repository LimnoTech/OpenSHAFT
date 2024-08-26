subroutine TPACalcBC(idum,jdum,JuncCase)

use GlobalVariables
use GlobalFunctions

implicit none

integer, intent(in) :: idum,jdum,JuncCase

if (A(idum,jdum) >= Reach(jdum)%Apipe) then
  hs(idum,jdum)=amax*amax*(A(idum,jdum)-Reach(jdum)%Apipe)/(g*Reach(jdum)%Apipe)
  y(idum,jdum)=Reach(jdum)%D+hs(idum,jdum)
  hc(idum,jdum)=FindCentroid(idum,jdum)
  Imom(idum,jdum)=A(idum,jdum)*hc(idum,jdum)+A(idum,jdum)*hs(idum,jdum)
  Tfs(idum,jdum)=FindTfs(idum,jdum)
  Rh(idum,jdum)=FindRh(idum,jdum)
else if (((Aold(idum,jdum) >= A(idum,jdum)) .and. (((Aold(idum+1,jdum) > 0.9999_wp*A(idum,jdum)) .or. &
  (vacuum(idum+1,jdum))) .and. &
  ((Aold(idum-1,jdum) > 0.9999_wp*A(idum,jdum)) .or. (vacuum(idum-1,jdum))))) .or. &
  ((vacuum(idum,jdum)) .and. (((Aold(idum+1,jdum) > 0.9999_wp*A(idum,jdum)) .or. (vacuum(idum+1,jdum))) .and. &
  ((Aold(idum-1,jdum) > 0.9999_wp*A(idum,jdum)) .or. (vacuum(idum-1,jdum)))))) then

  hs(idum,jdum)=amax*amax*(A(idum,jdum)-Reach(jdum)%Apipe)/(g*Reach(jdum)%Apipe)
  y(idum,jdum)=Reach(jdum)%D+hs(idum,jdum)
  hc(idum,jdum)=FindCentroid(idum,jdum)
  Imom(idum,jdum)=A(idum,jdum)*hc(idum,jdum)+A(idum,jdum)*hs(idum,jdum)
  Tfs(idum,jdum)=FindTfs(idum,jdum)
  Rh(idum,jdum)=FindRh(idum,jdum)
else
  if (A(idum,jdum) <= 0.0_dp) write (333,*) T,k,y(idum,jdum),A(idum,jdum),V(idum,jdum),idum,jdum,JuncCase
  hs(idum,jdum)=0.0_wp
  hc(idum,jdum)=FindCentroid(idum,jdum)
  if (isnan(hc(idum,jdum))) then
      write (6161,*) y(idum,jdum),hs(idum,jdum),vacuum(idum,jdum),idum,jdum,JuncCase,A(idum,jdum),Reach(jdum)%Apipe
      stop
  end if
  Rh(idum,jdum)=FindRh(idum,jdum)  
  Imom(idum,jdum)=A(idum,jdum)*hc(idum,jdum)
  Tfs(idum,jdum)=FindTfs(idum,jdum)
end if

if ((A(idum,jdum) > 0.0_wp) .and. (Tfs(idum,jdum) > 0.0_wp)) then
  c(idum,jdum)=min(amax,sqrt(g*A(idum,jdum)/Tfs(idum,jdum)))
end if
V(idum,jdum)=Q(idum,jdum)/A(idum,jdum)

if (hs(idum,jdum) < 0.0_wp) then
  vacuum(idum,jdum)=.true.
else
  vacuum(idum,jdum)=.false.
end if
	  
99 return
end subroutine TPACalcBC