subroutine CalcCellFluxes
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Updates flux terms in St. Venant eqs	
!%=====================================================================

use GlobalVariables

implicit none

!integer :: tid,nthreads
real(wp) :: tempA

!do j=1,NReaches
  do i=0,Reach(j)%Ncells
    f_c(i,j)=0.0_wp
  end do
!end do

!do j=1,Nreaches
  do i=1,Reach(j)%Ncells-1
    Ahat(i,j)=sqrt(A(i,j)*A(i+1,j))
!--> diagnostic error dump
    if (A(i,j) <= 0.0_dp) then
	  write (334,*) T,i,j,A(i,j)
	else if (A(i+1,j) <= 0.0_wp) then
	  write (334,*) T,i+1,j,A(i+1,j)
	end if
    if (isnan(Ahat(i,j))) then
      tempA=Ahat(i,j)
    end if
!--> end diagnostic dump
	Qhat(i,j)=(sqrt(A(i,j))*Q(i+1,j)+sqrt(A(i+1,j))*Q(i,j))/(sqrt(A(i,j))+sqrt(A(i+1,j)))
	if (abs(A(i,j)-A(i+1,j))/A(i,j) < 0.0001_wp) then
	  chat(i,j)=sqrt(g*(A(i,j)+A(i+1,j))/(Tfs(i,j)+Tfs(i+1,j)))
	else if ((Imom(i+1,j)-Imom(i,j))/(A(i+1,j)-A(i,j)) > 0) then
	  chat(i,j)=sqrt(g*(Imom(i+1,j)-Imom(i,j))/(A(i+1,j)-A(i,j)))
	else
	  chat(i,j)=sqrt(g*(A(i,j)+A(i+1,j))/(Tfs(i,j)+Tfs(i+1,j)))
	end if
	Lambda1(i,j)=Qhat(i,j)/Ahat(i,j) + chat(i,j)
	Lambda2(i,j)=Qhat(i,j)/Ahat(i,j) - chat(i,j)

    dw1(i,j)=+( (Q(i+1,j)-Q(i,j)) + (-Qhat(i,j)/Ahat(i,j) + chat(i,j))*(A(i+1,j)-A(i,j)) )
	dw2(i,j)=-( (Q(i+1,j)-Q(i,j)) + (-Qhat(i,j)/Ahat(i,j) - chat(i,j))*(A(i+1,j)-A(i,j)) )

  end do
!end do

!do j=1,Nreaches
  Reach(j)%HFnflocal=1.0_wp
  maxcroechange=0.0_wp
  do i=1,Reach(j)%Ncells-1
    if (abs(chat(i+1,j)-chat(i,j)) > maxcroechange) then
	  maxcroechange=abs(chat(i+1,j)-chat(i,j))
	  if (chat(i,j) > 0.5_wp*amax) then
	    Reach(j)%HFnflocal=HFnf
	  end if
	end if
  end do
  if (maxcroechange > 0.0_wp) then
	do i=1,Reach(j)%Ncells-1
	  f_c(i,j)=abs(chat(i+1,j)-chat(i,j))/maxcroechange
	end do
!  else
!    do i=1,Reach(j)%Ncells-1
!	  f_c(i,j)=0.0_wp
!	end do
  end if
!  f_c(0,j)=0.0_wp
!  f_c(Reach(j)%Ncells,j)=0.0_wp
!end do

if (dT == 0.0_wp) go to 99

!do j=1,Nreaches
    if (Reach(j)%HFnflocal == 0.0_wp) then
      ceil=0.0_wp
    else
      ceil=1.0_wp
    end if
    do i=1,Reach(j)%Ncells-1
      FintA(i,j)=0.5_wp*((FA(i,j)+FA(i+1,j)) - &
  	  (min(Reach(j)%dX/dT,(ceil*abs(Lambda1(i,j)) + (1-Reach(j)%HFnflocal)*Reach(j)%dX/dT*(f_c(i,j))**Reach(j)%HFnflocal))*dw1(i,j)*0.5_wp/chat(i,j)) - &
	  (min(Reach(j)%dX/dT,(ceil*abs(Lambda2(i,j)) + (1-Reach(j)%HFnflocal)*Reach(j)%dX/dT*(f_c(i,j))**Reach(j)%HFnflocal))*dw2(i,j)*0.5_wp/chat(i,j)))
      FintQ(i,j)=0.5_wp*((FQ(i,j)+FQ(i+1,j)) - &
	  (min(Reach(j)%dX/dT,(ceil*abs(Lambda1(i,j)) + (1-Reach(j)%HFnflocal)*Reach(j)%dX/dT*(f_c(i,j))**Reach(j)%HFnflocal))*dw1(i,j)*0.5_wp* &
	  Lambda1(i,j)/chat(i,j)) - &
	  (min(Reach(j)%dX/dT,(ceil*abs(Lambda2(i,j)) + (1-Reach(j)%HFnflocal)*Reach(j)%dX/dT*(f_c(i,j))**Reach(j)%HFnflocal))*dw2(i,j)*0.5_wp* &
	  Lambda2(i,j)/chat(i,j)))
    end do
!end do

    
99 continue


    
return
    
end subroutine CalcCellFluxes
