subroutine InnerSlopeChange
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Boundary condition to accomodate change in slope for otherwise constant cross sections	
!%=====================================================================

use GlobalVariables
use GlobalFunctions

implicit none

integer :: LeftCase,RightCase,iL,iR,jL,jR,kk
real(wp) :: r

LeftCase=0
RightCase=0
iL=Junc(k)%ReachCell(1)
jL=Junc(k)%ReachNo(1)
iR=Junc(k)%ReachCell(2)
jR=Junc(k)%ReachNo(2)

i=iL
j=jL

Ahat(i,j)=sqrt(A(i,j)*A(iR,jR))
Qhat(i,j)=(sqrt(A(i,j))*Q(iR,jR)+sqrt(A(iR,jR))*Q(i,j))/(sqrt(A(i,j))+sqrt(A(iR,jR)))
if (abs(A(i,j)-A(iR,jR))/A(i,j) < 0.0001_wp) then
  chat(i,j)=sqrt(g*(A(i,j)+A(iR,jR))/(Tfs(i,j)+Tfs(iR,jR)))
else if ((Imom(iR,jR)-Imom(i,j))/(A(iR,jR)-A(i,j)) > 0) then
  chat(i,j)=sqrt(g*(Imom(iR,jR)-Imom(i,j))/(A(iR,jR)-A(i,j)))
else
  chat(i,j)=sqrt(g*(A(i,j)+A(iR,jR))/(Tfs(i,j)+Tfs(iR,jR)))
end if
Lambda1(i,j)=Qhat(i,j)/Ahat(i,j) + chat(i,j)
Lambda2(i,j)=Qhat(i,j)/Ahat(i,j) - chat(i,j)

dw1(i,j)=+( (Q(iR,jR)-Q(i,j)) + (-Qhat(i,j)/Ahat(i,j) + chat(i,j))*(A(iR,jR)-A(i,j)) )
dw2(i,j)=-( (Q(iR,jR)-Q(i,j)) + (-Qhat(i,j)/Ahat(i,j) - chat(i,j))*(A(iR,jR)-A(i,j)) )

! update last cell of left reach

j=jL
i=iL
Aold(i,j)=A(i,j)
Qold(i,j)=Q(i,j)
r=dT/Reach(j)%dX

A(i,j)=A(i,j)-0.5_wp*r*(((FA(i,j)+FA(iR,jR))-(abs(Lambda1(i,j))*dw1(i,j)*0.5_wp/chat(i,j))-(abs(Lambda2(i,j))*dw2(i,j)*0.5_wp/chat(i,j))) &
       -((FA(i,j)+FA(i-1,j))-(abs(Lambda1(i-1,j))*dw1(i-1,j)*0.5_wp/chat(i-1,j))-(abs(Lambda2(i-1,j))*dw2(i-1,j)*0.5_wp/chat(i-1,j)))) &
	   +dT*SA(i,j)
Q(i,j)=Q(i,j)-0.5_wp*r*(((FQ(i,j)+FQ(iR,jR))-(abs(Lambda1(i,j))*dw1(i,j)*0.5_wp*Lambda1(i,j)/chat(i,j))- &
       (abs(Lambda2(i,j))*dw2(i,j)*0.5_wp*Lambda2(i,j)/chat(i,j))) &
       -((FQ(i,j)+FQ(i-1,j))-(abs(Lambda1(i-1,j))*dw1(i-1,j)*0.5_wp*Lambda1(i-1,j)/chat(i-1,j))- &
	   (abs(Lambda2(i-1,j))*dw2(i-1,j)*0.5_wp*Lambda2(i-1,j)/chat(i-1,j)))) &
       +dT*g*A(i,j)*Reach(j)%So

V(i,j)=Q(i,j)/A(i,j)
if (A(i,j) >= Reach(j)%Apipe) then 
  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
  y(i,j)=Reach(j)%D+hs(i,j)
else
  y(i,j)=FindDepth(i,j)
end if

if (A(i,j) >= Reach(j)%Apipe) then
      LeftCase=1
      hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
      y(i,j)=Reach(j)%D+hs(i,j)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
else if (((Aold(i,j) >= Reach(j)%Apipe) .and. (((Aold(iR,jR) > 0.9999_wp*Reach(jR)%Apipe) .or. (vacuum(iR,jR))) .and. &
        ((Aold(i-1,j) > 0.9999_wp*Reach(j)%Apipe) .or. (vacuum(i-1,j))))) .or. &
            ((vacuum(i,j)) .and. (((Aold(iR,jR) > 0.9999_wp*Reach(jR)%Apipe) .or. (vacuum(iR,jR))) .and. &
        ((Aold(i-1,j) > 0.9999_wp*Reach(j)%Apipe) .or. (vacuum(i-1,j)))))) then

      LeftCase=2
      hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
      y(i,j)=Reach(j)%D+hs(i,j)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
else
      LeftCase=3
      hs(i,j)=0.0_wp
      y(i,j)=FindDepth(i,j)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)
      Rh(i,j)=FindRh(i,j)
      Tfs(i,j)=FindTfs(i,j)
end if
if ((A(i,j) > 0.0_wp) .and. (Tfs(i,j) > 0.0_dp)) then
	  c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
end if
if (A(i,j) > 0.0_wp) then
      V(i,j)=Q(i,j)/A(i,j)
end if
if (hs(i,j) < 0.0_wp) then
	  vacuum(i,j)=.true.
else
	  vacuum(i,j)=.false.
end if

! update first cell of right reach

j=jR
i=iR

Aold(i,j)=A(i,j)
Qold(i,j)=Q(i,j)
r=dT/Reach(j)%dX

A(i,j)=A(i,j)-0.5_wp*r*(((FA(i,j)+FA(i+1,j))-(abs(Lambda1(i,j))*dw1(i,j)*0.5_wp/chat(i,j))-(abs(Lambda2(i,j))*dw2(i,j)*0.5_wp/chat(i,j))) &
       -((FA(i,j)+FA(iL,jL))-(abs(Lambda1(iL,jL))*dw1(iL,jL)*0.5_wp/chat(iL,jL))-(abs(Lambda2(iL,jL))*dw2(iL,jL)*0.5_wp/chat(iL,jL)))) &
	   +dT*SA(i,j)
Q(i,j)=Q(i,j)-0.5_wp*r*(((FQ(i,j)+FQ(i+1,j))-(abs(Lambda1(i,j))*dw1(i,j)*0.5_wp*Lambda1(i,j)/chat(i,j))- &
       (abs(Lambda2(i,j))*dw2(i,j)*0.5_wp*Lambda2(i,j)/chat(i,j))) &
       -((FQ(i,j)+FQ(iL,jL))-(abs(Lambda1(iL,jL))*dw1(iL,jL)*0.5_wp*Lambda1(iL,jL)/chat(iL,jL))- &
	   (abs(Lambda2(iL,jL))*dw2(iL,jL)*0.5_wp*Lambda2(iL,jL)/chat(iL,jL)))) &
       +dT*g*A(i,j)*Reach(j)%So

V(i,j)=Q(i,j)/A(i,j)
if (A(i,j) >= Reach(j)%Apipe) then 
  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
  y(i,j)=Reach(j)%D+hs(i,j)
else
  y(i,j)=FindDepth(i,j)
end if

if (A(i,j) >= Reach(j)%Apipe) then
      RightCase=4
      hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
      y(i,j)=Reach(j)%D+hs(i,j)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
else if (((Aold(i,j) >= Reach(j)%Apipe) .and. (((Aold(i+1,j) > 0.9999_wp*Reach(j)%Apipe) .or. (vacuum(i+1,j))) .and. &
        ((Aold(iL,jL) > 0.9999_wp*Reach(jL)%Apipe) .or. (vacuum(iL,jL))))) .or. &
            ((vacuum(i,j)) .and. (((Aold(i+1,j) > 0.9999_wp*Reach(j)%Apipe) .or. (vacuum(i+1,j))) .and. &
        ((Aold(iL,jL) > 0.9999_wp*Reach(jL)%Apipe) .or. (vacuum(iL,jL)))))) then

      RightCase=5
      hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
      y(i,j)=Reach(j)%D+hs(i,j)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
else
      RightCase=6
      hs(i,j)=0.0_wp
      y(i,j)=FindDepth(i,j)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)
      Rh(i,j)=FindRh(i,j)
      Tfs(i,j)=FindTfs(i,j)
end if
if ((A(i,j) > 0.0_wp) .and. (Tfs(i,j) > 0.0_wp)) then
	  c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
end if
if (A(i,j) > 0.0_wp) then
     V(i,j)=Q(i,j)/A(i,j)
end if
if (hs(i,j) < 0.0_wp) then
	  vacuum(i,j)=.true.
else
	  vacuum(i,j)=.false.
end if

Junc(k)%HeadOld=Junc(k)%Head
Junc(k)%Head=0.5_wp*(y(iL,jL)+y(iR,jR))

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,2i4,21f11.5)') T,LeftCase,RightCase,Junc(k)%Head,y(iL-1,jL),y(iL,jL),y(iR,jR),y(iR+1,jR),Q(iL-1,jL), &
	  Q(iL,jL),Q(iR,jR),Q(iR+1,jR),A(iL-1,jL),A(iL,jL),A(iR,jR),A(iR+1,jR),c(iL-1,jL),c(iL,jL),c(iR,jR),c(iR+1,jR),Tfs(iL-1,jL),Tfs(iL,jL),Tfs(iR,jR),Tfs(iR+1,jR)
	end if
  end do
end if

return

end subroutine