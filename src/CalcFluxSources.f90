subroutine CalcFluxSources
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Calculates source terms in St. Venant eqs	
!%=====================================================================
    
use GlobalVariables

implicit none

!do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    FA(i,j)=Q(i,j)
	FQ(i,j)=(Q(i,j)**2)/A(i,j) + g*Imom(i,j)
	Rey=abs(V(i,j))*y(i,j)/Nua + 0.001_wp
	if ((y(i,j) > Reach(j)%D) .and. (Rey < 2000.0_wp)) then
	  ndum=sqrt(8.0_wp/g*(Reach(j)%D/4.0_wp)**0.333_wp/Rey)
	  nn=max(Reach(j)%n,ndum)
	else
	  nn=Reach(j)%n
	end if
	if (Rh(i,j) > 0.0001_wp) then
	  Sf(i,j)=nn**2*V(i,j)*abs(V(i,j))/(Rh(i,j)**1.333_wp)
	else
	  Sf(i,j)=0.0_wp
	end if
	SA(i,j)=0._wp
	SQ(i,j)=g*A(i,j)*(Reach(j)%So-Sf(i,j))
  end do
!end do

return

end subroutine CalcFluxSources