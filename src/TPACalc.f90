subroutine TPACalc
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Performs two-component pressure approach calculations for interior reach cells	
!%=====================================================================

use GlobalVariables
use GlobalFunctions

implicit none

!do j=1,Nreaches
  do i=2,Reach(j)%Ncells-1
    if (not(PocketZone(i,j))) then
      if (A(i,j) >= Reach(j)%Apipe) then
        hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
        y(i,j)=Reach(j)%D+hs(i,j)
        hc(i,j)=FindCentroid(i,j)
        Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
        Tfs(i,j)=FindTfs(i,j)
        Rh(i,j)=FindRh(i,j)
      else if (((Aold(i,j) >= Reach(j)%Apipe) .and. (((Aold(i+1,j) > 0.9999*Reach(j)%Apipe) .or. (vacuum(i+1,j))) .and. &
        ((Aold(i-1,j) > 0.9999*Reach(j)%Apipe) .or. (vacuum(i-1,j))))) .or. &
            ((vacuum(i,j)) .and. (((Aold(i+1,j) > 0.9999*Reach(j)%Apipe) .or. (vacuum(i+1,j))) .and. &
        ((Aold(i-1,j) > 0.9999*Reach(j)%Apipe) .or. (vacuum(i-1,j)))))) then

        hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
        y(i,j)=Reach(j)%D+hs(i,j)
        hc(i,j)=FindCentroid(i,j)
        Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
        Tfs(i,j)=FindTfs(i,j)
        Rh(i,j)=FindRh(i,j)
	  else
        hs(i,j)=0.0_wp
        y(i,j)=FindDepth(i,j)
        hc(i,j)=FindCentroid(i,j)
        Imom(i,j)=A(i,j)*hc(i,j)
        Rh(i,j)=FindRh(i,j)
        Tfs(i,j)=FindTfs(i,j)
	  end if
	  if (A(i,j) > 0.0_dp) then
	    c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
	  else
	    c(i,j)=0.
	  end if
	  V(i,j)=Q(i,j)/A(i,j)
    end if
  end do
!end do

!do j=1,Nreaches
  do i=1,Reach(j)%Ncells-1
    if (hs(i,j) < 0.0_dp) then
	  vacuum(i,j)=.true.
	else
	  vacuum(i,j)=.false.
	end if
  end do
!end do

return

end subroutine TPACalc