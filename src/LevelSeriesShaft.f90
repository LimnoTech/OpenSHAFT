subroutine LevelSeriesShaft
    
use GlobalVariables
use GlobalFunctions

implicit none

real(wp), save :: JuncHeadTemp,TempDiff,D3,Hc3,H3O,Z3,JElev,r3,Vr3,Cr3,Sfr3,yr3,kr3,dV31,dV32,dV33,dV34
integer(i4) :: i3,j3,kk
logical :: IsSurch3,IsDrowned3,IsWeir3,CaseProblem

if (dT == 0.0_wp) then
  go to 900
end if

i3=Junc(k)%ReachCell(1)
j3=Junc(k)%ReachNo(1)
if (Junc(k)%ReachCount /= 0) D3=Reach(j3)%D

JElev=Junc(k)%Elev
Junc(k)%HeadOld=Junc(k)%Head
CaseProblem=.false.

Z3=Junc(k)%ReachElev(1)-JElev

call FindCurrentJuncHead

if (abs(Q(i3-1,j3)) < 4.0D-05*sqrt(g*D3**5)) then
  go to 900
end if

! Find critical depth for upstream reach

Hc3=FindCriticalDepth(i3-1,j3)

H3O=min(Hc3,y(i3-1,j3))

! set values for characteristic equations

r3=dT/Reach(j3)%dX
Vr3=(V(i3,j3) + r3*(-V(i3,j3)*c(i3-1,j3)+c(i3,j3)*V(i3-1,j3))) / (1 +r3*(V(i3,j3)-V(i3-1,j3)+c(i3,j3)-c(i3-1,j3)))
cr3=(c(i3,j3) + r3*Vr3*(c(i3-1,j3) - c(i3,j3))) / (1 + r3*(c(i3,j3) - c(i3-1,j3)))
Sfr3=Reach(j3)%n**2*Vr3*abs(Vr3)/Rh(i3,j3)**1.333
yr3=y(i3-1,j3) + r3*(Vr3+cr3)*(y(i3,j3) - y(i3-1,j3))
Kr3=Vr3 + g*yr3/cr3 - g*(Sfr3-Reach(j3)%So)*dT

! establish criteria for case selection
if (Junc(k)%HeadOld < Hc3+Z3) then
  IsWeir3=.true.
else
  IsWeir3=.false.
end if

if (V(i3-1,j3) > c(i3-1,j3)) then
  if (Junc(k)%HeadOld < alpha4*D3+Z3) then
    IsSurch3=.false.
  else
    IsSurch3=.true.
  end if
else if (y(i3,j3) < 0.999_wp*D3) then
  IsSurch3=.false.
else
  IsSurch3=.true.
end if

!if (JuncHeadOld(k) < Z3+Hc3) then
if (Junc(k)%HeadOld < Z3+D3) then
  IsDrowned3=.false.
else
  IsDrowned3=.true.
end if

JunctionCase=0

if ((IsWeir3) .and. (not(IsSurch3))) JunctionCase=1
if ((not(IsWeir3)) .and. (not(IsSurch3))) JunctionCase=2
if (IsSurch3) JunctionCase=3

if (JunctionCase == 0) then
  CaseProblem=.true.
end if

if (CaseProblem) then
  write (538,*) k,JunctionCase,IsSurch3,IsWeir3
  go to 900
end if

select case(JunctionCase)
    
  case(1)
    
    if  (Q(i3-1,j3) > 0.6D+0*(2.0D+0/3.0D+0)*sqrt(2.0D+0*g)*0.5D+0*Tfs(i3,j3)*(max(0.,y(i3,j3))**1.5D+0)) then
      Q(i3,j3)=Q(i3-1,j3)
      A(i3,j3)=A(i3-1,j3)
      V(i3,j3)=V(i3-1,j3)
      if (V(i3-1,j3) > c(i3-1,j3)) then
        y(i3,j3)=y(i3-1,j3)
      else
        y(i3,j3)=Hc3
      end if
      c(i3,j3)=FindCel(i3,j3)
    else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
      Q(i3,j3)=0.6D+0*(2.0D+0/3.0D+0)*sqrt(2.0D+0*g)*0.5D+0*Tfs(i3-1,j3)*(max(0.,y(i3-1,j3))**1.5D+0)
      V(i3,j3)=c(i3-1,j3)
      A(i3,j3)=Q(i3,j3)/V(i3,j3)
      c(i3,j3)=c(i3-1,j3)
      y(i3,j3)=FindDepth(i3,j3)
    end if

  case(2)

    y(i3,j3)=Junc(k)%Head-abs(Z3) !assume that the previous value of the junction head will be the final depth for reach 1
!    y(i3,j3)=JuncHead(k)-Z3+V(i3-1,j3)*sqrt(g/y(i3-1,j3))
    hs(i3,j3)=0.0D+0
    if (y(i3,j3) < D3) then
	  A(i3,j3)=FindArea(i3,j3)
    else
      hs(i3,j3)=y(i3,j3)-D3
      A(i3,j3)=Reach(j3)%Apipe + hs(i3,j3)*g*Reach(j3)%Apipe/amax**2
    end if

    V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
    c(i3,j3)=FindCel(i3,j3)
    if ((V(i3,j3) > 0.0D+0) .and. (V(i3,j3) > c(i3,j3))) V(i3,j3)=c(i3,j3)
    if ((V(i3,j3) < 0.0D+0) .and. (abs(V(i3,j3)) > c(i3,j3))) V(i3,j3)=-c(i3,j3)
    Q(i3,j3)=V(i3,j3)*A(i3,j3)

! last correction for bore invading junction
    if ((y(i3,j3)+y(i3-1,j3) > 100*initHfrac*D3) .and. (y(i3,j3) < 0.7_wp*y(i3-1,j3))) then
      Q(i3,j3)=Q(i3-1,j3)
      A(i3,j3)=A(i3-1,j3)
      V(i3,j3)=V(i3-1,j3)
      y(i3,j3)=y(i3-1,j3)
    end if

  case(3)
    
    if (IsDrowned3) then
      dV31=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      V(i3,j3)=V(i3,j3) + (dV31+2*dV32+2*dV33+dV34)/6.0D+0
	  y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*D3)
      A(i3,j3)=FindArea(i3,j3)
      c(i3,j3)=FindCel(i3,j3)
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    else
!      dV31=dT*g/dX(j3)*(-Hc3 + y(i3-1,j3) - Kdow.n(j3)*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
!      dV32=dT*g/dX(j3)*(-Hc3 + y(i3-1,j3) - Kdown(j3)*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
!      dV33=dT*g/dX(j3)*(-Hc3 + y(i3-1,j3) - Kdown(j3)*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
!      dV34=dT*g/dX(j3)*(-Hc3 + y(i3-1,j3) - Kdown(j3)*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
!      V(i3,j3)=V(i3,j3) + (dV31+2*dV32+2*dV33+dV34)/6.0D+0
!	  y(i3,j3)=max(Hc3+Kdown(j3)*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*D3)
      y(i3,j3)=y(i3-1,j3)-Reach(j3)%dX*Reach(j3)%n**2*V(i3-1,j3)*abs(V(i3-1,j3))/Rh(i3-1,j3)**1.333
      Q(i3,j3)=Q(i3-1,j3)
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    end if


  end select

Junc(k)%Inflow=0.0_wp
Junc(k)%Outflow=Q(i3,j3)

call TPACalcBC(i3,j3,JunctionCase)

98 continue

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,i5,26f9.3)') T,JunctionCase,Junc(k)%Head,Junc(k)%Inflow,y(i3,j3),y(i3-1,j3),V(i3,j3),V(i3-1,j3),c(i3,j3),c(i3-1,j3),Q(i3,j3),Q(i3-1,j3),Hc3
	end if
  end do
end if

900 continue

return

end subroutine LevelSeriesShaft