subroutine ThreeWaySplit
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Boundary condition for three reaches: two at downstream end, two at upstream end	
!%=====================================================================

use GlobalVariables
use GlobalFunctionsone

implicit none

real(wp) ::    JElev,Kloss,Hup,Hdown,Hc1,A1,V1,C1,Tfs1,Hc3,A3,V3,C3,Tfs3,Qlim1,Qlim2,Qlim3,Crown1,Crown2,Crown3,Hconj1
real(wp) ::    Z1,Z2,Z3,D1,D2,D3,Vr1,cr1,Sfr1,yr1,Kr1,r1,Vs2,cs2,Sfs2,ys2,Ks2,r2,Vs3,cs3,Sfs3,ys3,Ks3,r3
real(wp) ::    dV11,dV12,dV13,dV14,dV21,dV22,dV23,dV24,dV31,dV32,dV33,dV34,dY1,dY2,dY3,dY4,V1old,V2old,V3old
real(wp) ::    B1,B2,B3,Rf1,Rf2,Rf3,f1,f2,f3,CP1,CP3,CM2,BP1,BP3,BM2,SB,SC
integer ::    niter,i1,i2,i3,j1,j2,j3,kk
logical ::    converged,IsWeir1,IsWeir3,IsSurch1,IsSurch3,IsSurch2,IsNoShaft,CaseProblem,IsDrowned1

if (dT == 0.0_dp) then
  go to 99
end if

i1=Junc(k)%ReachCell(1)    !adjacent cell from upstream branch
i2=Junc(k)%ReachCell(2)    !adjacent cell from first downstream branch
i3=Junc(k)%ReachCell(3)    !adjacent cell from second downstream branch

j1=Junc(k)%ReachNo(1)      !reach ID of upstream branch
j2=Junc(k)%ReachNo(2)      !reach ID of first downstream branch
j3=Junc(k)%ReachNo(3)      !reach ID of second downstream branch

D1=Reach(j1)%D
D2=Reach(j2)%D
D3=Reach(j3)%D

JElev=Junc(k)%Elev
Junc(k)%HeadOld=Junc(k)%Head
CaseProblem=.false.
niter=0

Z1=Junc(k)%ReachElev(1)-JElev
Z2=Junc(k)%ReachElev(2)-JElev
Z3=Junc(k)%ReachElev(3)-JElev

Crown1=Z1+D1
Crown2=Z2+D2
Crown3=Z3+D3

V1old=V(i1,j1)
V2old=V(i2,j2)
V3old=V(i3,j3)

if (IsHGLInit) then
  if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0.0_wp) .and. (abs(Q(i1-1,j1)) < 4e-4_wp*sqrt(g*Reach(j1)%D**5)) .and. (abs(Q(i2+1,j2)) < 4e-4_wp*sqrt(g*Reach(j2)%D**5)) &
      .and. (abs(Q(i3+1,j3)) < 4e-4_wp*sqrt(g*Reach(j3)%D**5))) then
    go to 99
  end if
end if

if (not(IsHGLInit)) then
  if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0._wp) .and. (abs(Q(i1-1,j1)) < 4e-5_wp*sqrt(g*Reach(j1)%D**5)) .and. (abs(Q(i2+1,j2)) < 4e-5_wp*sqrt(g*Reach(j2)%D**5)) &
      .and. (abs(Q(i3+1,j3)) < 4e-5_wp*sqrt(g*Reach(j3)%D**5))) then
    go to 99
  end if
end if

!find critical depth for upstream reach

Hc1=FindCriticalDepth(i1-1,j1)

if (V(i1-1,j1) > c(i1-1,j1)) then
  Hconj1=y(i1-1,j1)*0.5_wp*(sqrt(1.0_wp +8.0_wp*V(i1-1,j1)*V(i1-1,j1)/(g*y(i1-1,j1)))-1.0_wp)
else
  Hconj1=0.0_wp
end if
! find limiting flow for convergence criterion

Qlim1=4e-4*sqrt(g*Reach(j1)%D**5)
Qlim2=4e-4*sqrt(g*Reach(j2)%D**5)
Qlim3=4e-4*sqrt(g*Reach(j3)%D**5)
Qlim=min(max(Qlim1,abs(Q(i1-1,j1))),max(Qlim2,abs(Q(i2+1,j2))),max(Qlim3,abs(Q(i3+1,j3))))

! set values for characteristic equations

r1=dT/Reach(j1)%dX
Vr1=(V(i1,j1) + r1*(-V(i1,j1)*c(i1-1,j1)+c(i1,j1)*V(i1-1,j1))) / (1.0_wp +r1*(V(i1,j1)-V(i1-1,j1)+c(i1,j1)-c(i1-1,j1)))
cr1=(c(i1,j1) + r1*Vr1*(c(i1-1,j1) - c(i1,j1))) / (1.0_wp + r1*(c(i1,j1) - c(i1-1,j1)))
Sfr1=Reach(j1)%n**2*Vr1*abs(Vr1)/Rh(i1,j1)**1.333_wp
yr1=y(i1-1,j1) + r1*(Vr1+cr1)*(y(i1,j1) - y(i1-1,j1))
Kr1=Vr1 + g*yr1/cr1 - g*(Sfr1-Reach(j1)%So)*dT

r2=dT/Reach(j2)%dX
Vs2=(V(i2,j2) + r2*(c(i2,j2)*V(i2+1,j2)-c(i2+1,j2)*V(i2,j2))) / (1.0_wp + r2*(-V(i2,j2)+V(i2+1,j2)+c(i2,j2)-c(i2+1,j2)))
cs2=(c(i2,j2) + r2*Vs2*(c(i2,j2) - c(i2+1,j2))) / (1.0_wp + r2*(c(i2,j2) - c(i2+1,j2)))
Sfs2=Reach(j2)%n**2*Vs2*abs(Vs2)/Rh(i2,j2)**1.333_wp
ys2=y(i2,j2) - r2*(Vs2-cs2)*(y(i2+1,j2) - y(i2,j2))
Ks2=Vs2 - g*ys2/cs2 - g*(Sfs2-Reach(j2)%So)*dT

r3=dT/Reach(j3)%dX
Vs3=(V(i3,j3) + r3*(c(i3,j3)*V(i3+1,j3)-c(i3+1,j3)*V(i3,j3))) / (1.0_wp + r3*(-V(i3,j3)+V(i3+1,j3)+c(i3,j3)-c(i3+1,j3)))
cs3=(c(i3,j3) + r3*Vs3*(c(i3,j3) - c(i3+1,j3))) / (1.0_wp + r3*(c(i3,j3) - c(i3+1,j3)))
Sfs3=Reach(j3)%n**2*Vs3*abs(Vs3)/Rh(i3,j3)**1.333_wp
ys3=y(i3,j3) - r3*(Vs3-cs3)*(y(i3+1,j3) - y(i3,j3))
Ks3=Vs3 - g*ys3/cs3 - g*(Sfs3-Reach(j3)%So)*dT

! establish criteria for case selection


if (V(i1-1,j1) > c(i1-1,j1)) then
  if ((Junc(k)%HeadOld < Hconj1+Z1) .and. (Q(i1-1,j1) > 0.0_wp)) then
    IsWeir1=.true.
  else
    IsWeir1=.false.
  end if
else if ((Junc(k)%HeadOld < Hc1+Z1) .and. (Q(i1-1,j1) > 0.0_wp)) then
  IsWeir1=.true.
else
  IsWeir1=.false.
end if

! establish surcharge on branch 1
if (V(i1-1,j1) > c(i1-1,j1)) then
  if (Junc(k)%HeadOld < alpha4*Reach(j1)%D+Z1) then
    IsSurch1=.false.
  else
    IsSurch1=.true.
  end if
else if (y(i1,j1) < 0.999*Reach(j1)%D) then
  IsSurch1=.false.
else
  IsSurch1=.true.
end if

! establish surcharge on branch 2
if (Junc(k)%HeadOld < 0.999*Reach(j2)%D+Z2) then
  IsSurch2=.false.
else
  IsSurch2=.true.
end if

! establish surcharge on branch 3
if  (Crown3 == Crown2) then
  if (IsSurch2) then
    IsSurch3=.true.
  else
    IsSurch3=.false.
  end if
else if (Junc(k)%HeadOld < 0.999*Reach(j3)%D+Z3) then
  IsSurch3=.false.
else
  IsSurch3=.true.
end if

! check for drowned condition on branch 1

if (Z1 == 0.0) then
  if (IsSurch1) then
    IsDrowned1=.true.
  else
    IsDrowned1=.false.
  end if
else if (Junc(k)%HeadOld < Z1+Hc1) then
  IsDrowned1=.false.
else
  IsDrowned1=.true.
end if
if (Q(i1,j1) < 0.0) IsDrowned1=.true.

if (Junc(k)%CurrentArea == 0.0_dp) then
  IsNoShaft=.true.
  Junc(k)%CurrentArea=0.25*pi*D2**2
else
  IsNoShaft=.false.
end if

JunctionCase=0

if ((IsWeir1) .and. (not(IsSurch2)) .and. (not(IsSurch3))) JunctionCase=1
if ((IsWeir1) .and. (IsSurch2) .and. (not(IsSurch3))) JunctionCase=2
if ((IsWeir1) .and. (not(IsSurch2)) .and. (IsSurch3)) JunctionCase=3
if ((IsWeir1) .and. (IsSurch2) .and. (IsSurch3)) JunctionCase=4
if ((not(IsWeir1)) .and. (not(IsSurch2)) .and. (not(IsSurch3))) JunctionCase=5
if ((not(IsWeir1)) .and. (IsSurch2) .and. (not(IsSurch3))) JunctionCase=6
if ((not(IsWeir1)) .and. (not(IsSurch2)) .and. (IsSurch3)) JunctionCase=7
if ((not(IsWeir1)) .and. (IsSurch2) .and. (IsSurch3)) JunctionCase=8
if ((IsSurch1) .and. (not(IsSurch2)) .and. (not(IsSurch3)) .and. (not(IsNoShaft))) JunctionCase=9
if ((IsSurch1) .and. (IsSurch2) .and. (not(IsSurch3)) .and. (not(IsNoShaft))) JunctionCase=10
if ((IsSurch1) .and. (not(IsSurch2)) .and. (IsSurch3) .and. (not(IsNoShaft))) JunctionCase=11
if ((IsSurch1) .and. (IsSurch2) .and. (IsSurch3) .and. (not(IsNoShaft))) JunctionCase=12

if (JunctionCase == 0) then
  CaseProblem=.true.
end if

if (CaseProblem) then
  write (538,*) k,JunctionCase,IsWeir1,IsWeir3,IsSurch1,IsSurch3,IsSurch2
  write (538,'(4f11.5)')Junc(k)%HeadOld,y(i1,j1),y(i2,j2),y(i3,j3)
  go to 99
end if

select case(JunctionCase)

  case(1)

  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    y(i1,j1)=Hc1
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

!    bisection to find pressure head at junction
  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
   Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Z2-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)
    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < 1.05*y(i2,j2)) .or. (y(i2+2,j2) < 1.10*y(i2,j2))) then
        if (V(i2,j2) > 0) then
          V(i2,j2)=c(i2,j2)
          Q(i2,j2)=A(i2,j2)*V(i2,j2)
        else
          V(i2,j2)=-c(i2,j2)
          Q(i2,j2)=A(i2,j2)*V(i2,j2)
        end if
      else
        A(i2,j2)=A(i2+1,j2)
        V(i2,j2)=V(i2+1,j2)
        Q(i2,j2)=Q(i2+1,j2)
      end if
    else
      Q(i2,j2)=A(i2,j2)*V(i2,j2)
    end if

    y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3-Reach(j3)%dX/2*Reach(j3)%So)
    c(i3,j3)=FindCel(i3,j3)
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Ks3 + y(i3,j3)*g/cs3
    Q(i3,j3)=V(i3,j3)*A(i3,j3)
    if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3+1,j3) < 1.05*y(i3,j3)) .or. (y(i3+2,j3) < 1.10*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3+1,j3)
        V(i3,j3)=V(i3+1,j3)
        Q(i3,j3)=Q(i3+1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)+ Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
      Hup=Junc(k)%Head
    else
      Hdown=Junc(k)%Head
    end if
    if (Hup < Hdown) then
     Junc(k)%Head=Hup
      Hup=Hdown
      Hdown=Junc(k)%Head
    end if
    if (((abs((Hup-Hdown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

  end do
  
  case(2)
      
  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    y(i1,j1)=Hc1
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

!    bisection to find pressure head at junction
  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
   Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3-Reach(j3)%dX/2*Reach(j3)%So)
    c(i3,j3)=FindCel(i3,j3)
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Ks3 + y(i3,j3)*g/cs3
    Q(i3,j3)=V(i3,j3)*A(i3,j3)
    if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3+1,j3) < 1.05*y(i3,j3)) .or. (y(i3+2,j3) < 1.10*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3+1,j3)
        V(i3,j3)=V(i3+1,j3)
        Q(i3,j3)=Q(i3+1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if

    dV21= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dV22= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dV23= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dV24= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    V(i2,j2)=V2old + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)

    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=findcel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)
    
    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)+ Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
      Hup=Junc(k)%Head
    else
      Hdown=Junc(k)%Head
    end if
    if (Hup < Hdown) then
     Junc(k)%Head=Hup
      Hup=Hdown
      Hdown=Junc(k)%Head
    end if
    if (((abs((Hup-Hdown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

  end do
  
  case(3)

  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    y(i1,j1)=Hc1
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

!    bisection to find pressure head at junction
  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
   Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Z2-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)
    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < 1.05*y(i2,j2)) .or. (y(i2+2,j2) < 1.10*y(i2,j2))) then
        if (V(i2,j2) > 0) then
          V(i2,j2)=c(i2,j2)
          Q(i2,j2)=A(i2,j2)*V(i2,j2)
        else
          V(i2,j2)=-c(i2,j2)
          Q(i2,j2)=A(i2,j2)*V(i2,j2)
        end if
      else
        A(i2,j2)=A(i2+1,j2)
        V(i2,j2)=V(i2+1,j2)
        Q(i2,j2)=Q(i2+1,j2)
      end if
    else
      Q(i2,j2)=A(i2,j2)*V(i2,j2)
    end if

    dV31= dT*g/Reach(j3)%dX * (Junc(k)%Head-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
    dV32= dT*g/Reach(j3)%dX * ( (Junc(k)%Head+dY1/2-Z3) - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
    dV33= dT*g/Reach(j3)%dX * ( (Junc(k)%Head+dY2/2-Z3) - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
    dV34= dT*g/Reach(j3)%dX * ( (Junc(k)%Head+dY3-Z3) - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
    V(i3,j3)=V3old + (dV31+2*(dV32+dV33)+dV34)/6.0_wp
    y(i3,j3)=max(Junc(k)%Head-Z3-Reach(j3)%Kup*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=findcel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)+ Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
      Hup=Junc(k)%Head
    else
      Hdown=Junc(k)%Head
    end if
    if (Hup < Hdown) then
     Junc(k)%Head=Hup
      Hup=Hdown
      Hdown=Junc(k)%Head
    end if
    if (((abs((Hup-Hdown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

  end do
  
  case(4)
      
  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    y(i1,j1)=Hc1
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

  dV21= dT*g/Reach(j2)%dX * (Junc(k)%Head-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
  dV31= dT*g/Reach(j3)%dX * (Junc(k)%Head-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
  dY1=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

  dV22= dT*g/Reach(j2)%dX * (Junc(k)%Head+dY1/2-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
  dV32= dT*g/Reach(j3)%dX * (Junc(k)%Head+dY1/2-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
  dY2=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV31/2 ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

  dV23= dT*g/Reach(j2)%dX * (Junc(k)%Head+dY2/2-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
  dV33= dT*g/Reach(j3)%dX * (Junc(k)%Head+dY2/2-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
  dY3=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3+dV32/2) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

  dV24= dT*g/Reach(j2)%dX * (Junc(k)%Head+dY3-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
  dV34= dT*g/Reach(j3)%dX * (Junc(k)%Head+dY3-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
  dY4=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV33 ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

  V(i2,j2)=V2old + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
  V(i3,j3)=V3old + (dV31+2*(dV32+dV33)+dV34)/6.0_wp
  Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

  y(i2,j2)=Junc(k)%Head - Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g)
  A(i2,j2)=FindArea(i2,j2)
  c(i2,j2)=FindCel(i2,j2)
  Q(i2,j2)=A(i2,j2)*V(i2,j2)
	
  y(i3,j3)=max(Junc(k)%Head-Z3-Reach(j3)%Kup*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
  A(i3,j3)=FindArea(i3,j3)
  c(i3,j3)=findcel(i3,j3)
  Q(i3,j3)=A(i3,j3)*V(i3,j3)

  case(5)
  
 Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
   Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i1,j1)=max(initHfrac*Reach(j1)%D,Junc(k)%Head-Z1+Reach(j1)%dX/2*Reach(j1)%So)
    c(i1,j1)=FindCel(i1,j1)
    A(i1,j1)=FindArea(i1,j1)
    V(i1,j1)=Kr1 - y(i1,j1)*g/cr1
    Q(i1,j1)=V(i1,j1)*A(i1,j1)
    if (V(i1-1,j1) > c(i1-1,j1)) then
      Q(i1,j1)=Q(i1-1,j1)
      V(i1,j1)=Q(i1,j1)/A(i1,j1)
    else if (abs(V(i1,j1)) > c(i1,j1)) then
      if ((y(i1-1,j1) < 1.05*y(i1,j1)) .or. (y(i1-2,j1)<1.10*y(i1,j1))) then
        if (V(i1,j1) > 0) then
          V(i1,j1)=c(i1,j1)
          Q(i1,j1)=A(i1,j1)*V(i1,j1)
        else
          V(i1,j1)=-c(i1,j1)
          Q(i1,j1)=A(i1,j1)*V(i1,j1)
        end if
      else
        A(i1,j1)=A(i1-1,j1)
        y(i1,j1)=y(i1-1,j1)
        V(i1,j1)=V(i1-1,j1)
        Q(i1,j1)=Q(i1-1,j1)
      end if
    else
      Q(i1,j1)=A(i1,j1)*V(i1,j1)
    end if

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Z2-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)

    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < 1.05*y(i2,j2)) .or. (y(i2+2,j2) < 1.10*y(i2,j2))) then
        if (V(i2,j2) > 0) then
          V(i2,j2)=c(i2,j2)
          Q(i2,j2)=A(i2,j2)*V(i2,j2)
        else
          V(i2,j2)=-c(i2,j2)
          Q(i2,j2)=A(i2,j2)*V(i2,j2)
        end if
      else
        A(i2,j2)=A(i2+1,j2)
        V(i2,j2)=V(i2+1,j2)
        Q(i2,j2)=Q(i2+1,j2)
      end if
    else
      Q(i2,j2)=A(i2,j2)*V(i2,j2)
    end if

    y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3-Reach(j3)%dX/2*Reach(j3)%So)
    c(i3,j3)=FindCel(i3,j3)
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Ks3 + y(i3,j3)*g/cs3
    Q(i3,j3)=V(i3,j3)*A(i3,j3)
    if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3+1,j3) < 1.05*y(i3,j3)) .or. (y(i3+2,j3) < 1.10*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3+1,j3)
        V(i3,j3)=V(i3+1,j3)
        Q(i3,j3)=Q(i3+1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)+ Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
      Hup=Junc(k)%Head
    else
      Hdown=Junc(k)%Head
    end if
    if (Hup < Hdown) then
     Junc(k)%Head=Hup
      Hup=Hdown
      Hdown=Junc(k)%Head
    end if
    if (((abs((Hup-Hdown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

  end do

  case(12)
  
    dV11= dT*g/Reach(j1)%dX * (-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
    dV21= dT*g/Reach(j2)%dX * (Junc(k)%Head-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dV31= dT*g/Reach(j3)%dX * (Junc(k)%Head-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
    dY1=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

    dV12= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY1/2-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
    dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2-Z2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dV32= dT*g/Reach(j3)%dX * ( (Junc(k)%Head+dY1/2-Z3) - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
    dY2=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV31/2 ) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

    dV13= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY2/2-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
    dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2-Z2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dV33= dT*g/Reach(j3)%dX * ( (Junc(k)%Head+dY2/2-Z3) - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
    dY3=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV32/2 ) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

    dV14= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY3-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
    dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3-Z2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    dV34= dT*g/Reach(j3)%dX * ( (Junc(k)%Head+dY3-Z3) - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
    dY4=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV33 ) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

    V(i1,j1)=V1old + (dV11+2*(dV12+dV13)+dV14)/6.0_wp
    V(i2,j2)=V2old + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    V(i3,j3)=V3old + (dV31+2*(dV32+dV33)+dV34)/6.0_wp
   Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

	y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=findcel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=findcel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

	y(i3,j3)=max(Junc(k)%Head-Z3-Reach(j3)%Kup*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=findcel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

  
  case(9)
      
   Junc(k)%Head=0.001
    converged=.false.

    Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
    Hdown=0.001

    do while(not(converged))
     Junc(k)%Head=0.5*(Hup+Hdown)
      niter=niter+1

      y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Z2-Reach(j2)%dX/2*Reach(j2)%So)
      c(i2,j2)=FindCel(i2,j2)
      A(i2,j2)=FindArea(i2,j2)
      V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
      Q(i2,j2)=V(i2,j2)*A(i2,j2)
      if (abs(V(i2,j2)) > c(i2,j2)) then
        if ((y(i2+1,j2) < 1.05*y(i2,j2)) .or. (y(i2+2,j2) < 1.10*y(i2,j2))) then
          if (V(i2,j2) > 0) then
            V(i2,j2)=c(i2,j2)
            Q(i2,j2)=A(i2,j2)*V(i2,j2)
          else
            V(i2,j2)=-c(i2,j2)
            Q(i2,j2)=A(i2,j2)*V(i2,j2)
          end if
        else
          A(i2,j2)=A(i2+1,j2)
          V(i2,j2)=V(i2+1,j2)
          Q(i2,j2)=Q(i2+1,j2)
        end if
      else
        Q(i2,j2)=A(i2,j2)*V(i2,j2)
      end if

      y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3-Reach(j3)%dX/2*Reach(j3)%So)
      c(i3,j3)=FindCel(i3,j3)
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Ks3 + y(i3,j3)*g/cs3
      Q(i3,j3)=V(i3,j3)*A(i3,j3)
      if (abs(V(i3,j3)) > c(i3,j3)) then
        if ((y(i3+1,j3) < 1.05*y(i3,j3)) .or. (y(i3+2,j3) < 1.10*y(i3,j3))) then
          if (V(i3,j3) > 0) then
            V(i3,j3)=c(i3,j3)
            Q(i3,j3)=A(i3,j3)*V(i3,j3)
         else
            V(i3,j3)=-c(i3,j3)
            Q(i3,j3)=A(i3,j3)*V(i3,j3)
          end if
        else
          A(i3,j3)=A(i3+1,j3)
          V(i3,j3)=V(i3+1,j3)
          Q(i3,j3)=Q(i3+1,j3)
        end if
      else
        Q(i3,j3)=A(i3,j3)*V(i3,j3)
      end if

      if (IsDrowned1) then
        dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
        dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
        dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
        dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
        V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
	    y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
      else
        dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
        dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
        dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
        dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
        V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
	    y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
      end if    

      A(i1,j1)=FindArea(i1,j1)
      c(i1,j1)=findcel(i1,j1)
      Q(i1,j1)=A(i1,j1)*V(i1,j1)

      if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)+ Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
        Hup=Junc(k)%Head
      else
        Hdown=Junc(k)%Head
      end if
      if (Hup < Hdown) then
       Junc(k)%Head=Hup
        Hup=Hdown
        Hdown=Junc(k)%Head
      end if
      if (((abs((Hup-Hdown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

  end do

  case(10)
      
   Junc(k)%Head=0.001
    converged=.false.

    Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
    Hdown=0.001

    do while(not(converged))
     Junc(k)%Head=0.5*(Hup+Hdown)
      niter=niter+1

      y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3-Reach(j3)%dX/2*Reach(j3)%So)
      c(i3,j3)=FindCel(i3,j3)
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Ks3 + y(i3,j3)*g/cs3
      Q(i3,j3)=V(i3,j3)*A(i3,j3)
      if (abs(V(i3,j3)) > c(i3,j3)) then
        if ((y(i3+1,j3) < 1.05*y(i3,j3)) .or. (y(i3+2,j3) < 1.10*y(i3,j3))) then
          if (V(i3,j3) > 0) then
            V(i3,j3)=c(i3,j3)
            Q(i3,j3)=A(i3,j3)*V(i3,j3)
         else
            V(i3,j3)=-c(i3,j3)
            Q(i3,j3)=A(i3,j3)*V(i3,j3)
          end if
        else
          A(i3,j3)=A(i3+1,j3)
          V(i3,j3)=V(i3+1,j3)
          Q(i3,j3)=Q(i3+1,j3)
        end if
      else
        Q(i3,j3)=A(i3,j3)*V(i3,j3)
      end if

      if (IsDrowned1) then
        dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
        dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
        dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
        dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
        V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
	    y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
      else
        dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
        dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
        dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
        dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
        V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
	    y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
      end if    

      dV21= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dV22= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dV23= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dV24= dT*g/Reach(j2)%dX * (Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      V(i2,j2)=V2old + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
      y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)

      A(i1,j1)=FindArea(i1,j1)
      c(i1,j1)=findcel(i1,j1)
      Q(i1,j1)=A(i1,j1)*V(i1,j1)

      A(i2,j2)=FindArea(i2,j2)
      c(i2,j2)=findcel(i2,j2)
      Q(i2,j2)=A(i2,j2)*V(i2,j2)

      if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)+ Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
        Hup=Junc(k)%Head
      else
        Hdown=Junc(k)%Head
      end if
      if (Hup < Hdown) then
       Junc(k)%Head=Hup
        Hup=Hdown
        Hdown=Junc(k)%Head
      end if
      if (((abs((Hup-Hdown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

    end do
    
  case(8)
    
    dV21= dT*g/Reach(j2)%dX * (Junc(k)%Head-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dV31= dT*g/Reach(j3)%dX * (Junc(k)%Head-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
    dY1=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

    dV22= dT*g/Reach(j2)%dX * (Junc(k)%Head+dY1/2-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dV32= dT*g/Reach(j3)%dX * (Junc(k)%Head+dY1/2-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
    dY2=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV31/2 ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

    dV23= dT*g/Reach(j2)%dX * (Junc(k)%Head+dY2/2-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dV33= dT*g/Reach(j3)%dX * (Junc(k)%Head+dY2/2-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
    dY3=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3+dV32/2) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

    dV24= dT*g/Reach(j2)%dX * (Junc(k)%Head+dY3-Z2 - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    dV34= dT*g/Reach(j3)%dX * (Junc(k)%Head+dY3-Z3 - y(i3+1,j3) - Reach(j3)%Kup* ( V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
    dY4=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) - A(i3,j3)*( V(i3,j3)+dV33 ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

    V(i2,j2)=V2old + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    V(i3,j3)=V3old + (dV31+2*(dV32+dV33)+dV34)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    y(i2,j2)=Junc(k)%Head - Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)
	
    y(i3,j3)=max(Junc(k)%Head-Z3-Reach(j3)%Kup*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=findcel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

    y(i1,j1)=Junc(k)%Head-abs(Z1) !assume that the previous value of the junction head will be the final depth for reach 1
    hs(i1,j1)=0.0_wp
    if (y(i1,j1) < Reach(j1)%D) then
	  A(i1,j1)=FindArea(i1,j1)
    else
      hs(i1,j1)=y(i1,j1)-Reach(j1)%D
      A(i1,j1)=Reach(j1)%Apipe + hs(i1,j1)*g*Reach(j1)%Apipe/amax**2
    end if

    V(i1,j1)=Kr1 - y(i1,j1)*g/cr1
    c(i1,j1)=FindCel(i1,j1)
    if ((V(i1,j1) > 0.0_wp) .and. (V(i1,j1) > c(i1,j1))) V(i1,j1)=c(i1,j1)
    if ((V(i1,j1) < 0.0_wp) .and. (abs(V(i1,j1)) > c(i1,j1))) V(i1,j1)=-c(i1,j1)
    Q(i1,j1)=V(i1,j1)*A(i1,j1)

! last correction for bore invading junction
    if ((y(i1,j1)+y(i1-1,j1) > 100*initHfrac*Reach(j1)%D) .and. (y(i1,j1) < 0.7*y(i1-1,j1))) then
      Q(i1,j1)=Q(i1-1,j1)
      A(i1,j1)=A(i1-1,j1)
      V(i1,j1)=V(i1-1,j1)
      y(i1,j1)=y(i1-1,j1)
    end if

end select

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,i5,26f9.3,i5,f9.3)') T,JunctionCase,Junc(k)%Head,Junc(k)%Inflow,y(i1,j1),y(i2,j2),y(i3,j3),y(i1-1,j1),y(i2+1,j2),y(i3+1,j3), &
	  V(i1,j1),V(i2,j2),V(i3,j3),V(i1-1,j1),V(i2+1,j2),V(i3+1,j3),c(i1,j1),c(i2,j2),c(i3,j3),c(i1-1,j1),c(i2+1,j2),c(i3+1,j3),Q(i1,j1),Q(i2,j2),Q(i3,j3), &
	  Q(i1-1,j1),Q(i2+1,j2),Q(i3+1,j3),niter,Hconj1
	end if
  end do
end if


call TPACalcBC(i1,j1,JunctionCase)
call TPACalcBC(i2,j2,JunctionCase)
call TPACalcBC(i3,j3,JunctionCase)

99 continue

return

end subroutine ThreeWaySplit
