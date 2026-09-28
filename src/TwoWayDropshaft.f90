subroutine TwoWayDropshaft
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Boundary condition for two reaches, one downstream and one upstream	
!%=====================================================================

use GlobalVariables
use GlobalFunctions

implicit none

real(wp) :: Kmid,Hup,Hdown,dVn1,dVn2,dVn3,dVn4,dVm1,dVm2,dVm3,dVm4,dYm1,dYm2,dYm3,dYm4
real(wp) :: ZL,ZR,Zcorr,ReachZL,ReachZR,AL,CL,VL,HcL,TfsL,Kloss,JElev,Vlold,QlimL,QlimR,Hconj1,Rf1,Rf2,B1,B2,f1,f2,CP1,CM2,BP1,BM2,SB,SC
real(wp) :: Vs,cs,Sfs,ys,Ks,Vr,cr,Sfr,yr,Kr,r,Htotal

integer :: niter,iL,jL,iR,jR,kk

logical :: converged,IsSurchLow,IsSurchHigh,IsWeirFlow,IsDrownedUp,IsCrownMatch

real(wp) :: OldVolChamber,HchO,dVolChamber,Pdrop,VolExhaust,dHch1,dHch2,dHch3,dHch4
integer :: LeftCase,RightCase

LeftCase=0
RightCase=0

if (dT == 0.0_wp) then
  go to 99
end if

iL=Junc(k)%ReachCell(1)
jL=Junc(k)%ReachNo(1)
iR=Junc(k)%ReachCell(2)
jR=Junc(k)%ReachNo(2)
JElev=Junc(k)%Elev
ReachZL=Junc(k)%ReachElev(1)
ReachZR=Junc(k)%ReachElev(2)
niter=0

! Zcorr is the amount that we will remove from the junction head to make the ODE
! calculations at each node that arrives at the surge shaft

ZL=ReachZL
ZR=ReachZR
Zcorr=ZL-Jelev

!call FindCurrentShaftArea

Junc(k)%HeadOld=Junc(k)%Head
OldVolChamber=Junc(k)%VolChamber
HchO=Junc(k)%Hch

r=dT/Reach(jL)%dX
Vr=(V(iL,jL)+r*(-V(iL,jL)*c(iL-1,jL)+c(iL,jL)*V(iL-1,jL)))/(1.0_wp +r*(V(iL,jL)-V(iL-1,jL)+c(iL,jL)-c(iL-1,jL)))
cr=(c(iL,jL)+r*Vr*(c(iL-1,jL)-c(iL,jL)))/(1.0_wp +r*(c(iL,jL)-c(iL-1,jL)))
Sfr=Reach(jL)%n**2*Vr*abs(Vr)/Rh(iL,jL)**1.333
yr=y(iL,jL) + r*(Vr+cr)*(y(iL-1,jL) - y(iL,jL))
Kr=Vr + g*yr/cr - g*(Sfr-Reach(jL)%So)*dT

r=dT/Reach(jR)%dX
Vs=(V(iR,jR)+r*(-V(iR,jR)*c(iR+1,jR)+c(iR,jR)*V(iR+1,jR)))/(1.0_wp +r*(-V(iR,jR)+V(iR+1,jR)+c(iR,jR)-c(iR+1,jR)))
cs=(c(iR,jR)+r*Vs*(c(iR,jR)-c(iR+1,jR)))/(1.0_wp +r*((c(iR,jR)-c(iR+1,jR))))
Sfs=Reach(jR)%n**2*Vs*abs(Vs)/Rh(iR,jR)**1.333
ys=y(iR,jR)-r*(Vs-cs)*(y(iR+1,jR)-y(iR,jR))
Ks=Vs-g*ys/cs-g*(Sfs-Reach(jR)%So)*dT

!call FindOverflowValue

if (IsHGLInit) then
  if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0.0_dp) .and. (abs(Q(iR+1,jR)) < 1.0e-4_wp*sqrt(g*Reach(jR)%D**5)) .and. (abs(Q(iL-1,jL)) < 1.0e-4_wp*sqrt(g*Reach(jL)%D**5))) then
    go to 99
  end if
end if

QlimL=4e-4*sqrt(g*Reach(jL)%D**5)
QlimR=4e-4*sqrt(g*Reach(jR)%D**5)
!Qlim=min(max(QlimL,abs(Q(iL-1,jL))),max(QlimR,abs(Q(iR+1,jR))))
Qlim=min(QlimL,QlimR)

HcL=FindCriticalDepth(iL-1,jL)

if (V(iL-1,jL) > c(iL-1,jL)) then
  Hconj1=y(iL-1,jL)*0.5_wp*(sqrt(1.0_wp +8.0_wp*V(iL-1,jL)*V(iL-1,jL)/(g*y(iL-1,jL)))-1.0_wp)
else
  Hconj1=0.0_wp
end if

! establish case criteria

if (Junc(k)%HeadOld < 0.999*Reach(jR)%D) then
  IsSurchLow=.false.
else
  IsSurchLow=.true.
end if

if (V(iL-1,jL) > c(iL-1,jL)) then
  if (Junc(k)%HeadOld < alpha4*Reach(jL)%D+Zcorr) then
    IsSurchHigh=.false.
  else
    IsSurchHigh=.true.
  end if
else if (y(iL,jL) < 0.999_wp*Reach(jL)%D) then
  IsSurchHigh=.false.
else
  IsSurchHigh=.true.
end if

! PRK 6/17/2014 weir flow for left-hand side is not relevant in case of reverse flow so check this first
!if (Q(iL,jL) < 0.0) then
!  IsWeirFlow=.false.
!else if ((Junc(k)%HeadOld < Zcorr+HcL) .and. (V(iL-1,jL) < c(iL-1,jL))) then
!  IsWeirFlow=.true.
!    else
!  IsWeirFlow=.false.
!end if

if (V(iL-1,jL) > c(iL-1,jL)) then
  if ((Junc(k)%HeadOld < Hconj1+Zcorr) .and. (Q(iL-1,jL) > 0.0)) then
    IsWeirFlow=.true.
  else
    IsWeirFlow=.false.
  end if
else if ((Junc(k)%HeadOld < HcL+Zcorr) .and. (Q(iL-1,jL) > 0.0)) then
  IsWeirFlow=.true.
else
  IsWeirFlow=.false.
end if


if ((Reach(jL)%D+ReachZL) == (Reach(jR)%D+ReachZR)) then
    IsCrownMatch=.true.
else
    IsCrownMatch=.false.
end if

! PRK 6/13/2014 check if spur is actually drowned or just surcharged at downstream end (same as in ThreeWayDropshaft)
if ((IsCrownMatch) .or. (not(IsWeirFlow))) then
  if (IsSurchHigh) then
    IsDrownedUp=.true.
  else
    IsDrownedUp=.false.
  end if
else if (Junc(k)%HeadOld < 0.999_wp*(Reach(jL)%D+Zcorr)) then
!else if (JuncHeadOld(k) < (Reach(jL)%D+Zcorr)) then
  IsDrownedUp=.false.
else
  IsDrownedUp=.true.
end if
! PRK 6/13/2014 end

! choose appropriate case

JunctionCase=0

if ((not(IsSurchLow)) .and. (not(IsSurchHigh)) .and. (IsWeirFlow)) JunctionCase=1
if ((not(IsSurchLow)) .and. (not(IsSurchHigh)) .and. (not(IsWeirFlow))) JunctionCase=2
if ((IsSurchLow) .and. (not(IsSurchHigh)) .and. (IsWeirFlow)) JunctionCase=3
if ((IsSurchLow) .and. (not(IsSurchHigh)) .and. (not(IsWeirFlow))) JunctionCase=4
if ((IsSurchLow) .and. (IsSurchHigh) .and. (Junc(k)%CurrentArea > 0.0_dp)) JunctionCase=5
if ((not(IsSurchLow)) .and. (IsSurchHigh)) JunctionCase=6
if ((IsSurchLow) .and. (IsSurchHigh) .and. (Junc(k)%CurrentArea == 0.0_dp)) JunctionCase=7
if ((IsSurchLow) .and. (IsSurchHigh) .and. (Junc(k)%Option == -8)) JunctionCase=8

! PRK 4/20/2015 depressurize chamber if necessary
if ((not(Junc(k)%InitializeHch)) .and. (JunctionCase /= 8)) then
  Junc(k)%Hch=atm
  Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
  call FindCurrentChamberVolume
  Junc(k)%InitializeHch=.true.
  Junc(k)%Qair=0.0_wp
end if
! PRK 4/20/2015 end

! PRK 8/16/2016 bail on chamber pressure calc if WSE is too close to top of chamber
if ((JunctionCase == 8) .and. (Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head < 1.0D-02)) then
    JunctionCase=5
    Junc(k)%InitializeHch=.true.
end if
! PRK 8/16/2016 end

if (Junc(k)%CurrentArea == 0.0_dp) Junc(k)%CurrentArea=0.25*pi*Reach(jR)%D*Reach(jR)%D

! PRK 2/17/2017 test corrective action for no-shaft junctions to eas transition to surcharge mode
!if ((JunctionCase == 7) .and. (abs(Q(iL,jL)-Q(iR,jR))/sqrt(g*D(jL)**5) < 0.005_wp)) JunctionCase=5
! PRK 2/17/2016 end

select case(JunctionCase)

  case(1)

    Junc(k)%Head=0.001
    converged=.false.

    Hup=1.00*max(Reach(jL)%D,Reach(jR)%D)
    Hdown=0.001

    do while(not(converged))
      Junc(k)%Head=0.5*(Hup+Hdown)
      niter=niter+1

      y(iR,jR)=max(initHfrac*Reach(jR)%D,Junc(k)%Head-0.5*Reach(jR)%dX*Reach(jR)%So)
      c(iR,jR)=FindCel(iR,jR)
      A(iR,jR)=FindArea(iR,jR)
      V(iR,jR)=Ks + y(iR,jR)*g/cs
      Q(iR,jR)=V(iR,jR)*A(iR,jR)
      if (abs(V(iR,jR)) > c(iR,jR)) then
        if ((y(iR+1,jR) < 1.05*y(iR,jR)) .or. (y(iR+2,jR)<1.10*y(iR,jR))) then
          if (V(iR,jR) > 0) then
            V(iR,jR)=c(iR,jR)
            Q(iR,jR)=A(iR,jR)*V(iR,jR)
          else
            V(iR,jR)=-c(iR,jR)
            Q(iR,jR)=A(iR,jR)*V(iR,jR)
          end if
        else
          A(iR,jR)=A(iR+1,jR)
          V(iR,jR)=V(iR+1,jR)
          Q(iR,jR)=Q(iR+1,jR)
        end if
      else
        Q(iR,jR)=A(iR,jR)*V(iR,jR)
      end if

      if  (Q(iL-1,jL) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(iL,jL)*(max(0.0_wp,y(iL,jL)))**1.5) then
        Q(iL,jL)=Q(iL-1,jL)
        V(iL,jL)=V(iL-1,jL)
        if (V(iL-1,jL) > c(iL-1,jL)) then
          y(iL,jL)=y(iL-1,jL)
        else
          y(iL,jL)=HcL
        end if
        A(iL,jL)=FindArea(iL,jL)
        c(iL,jL)=FindCel(iL,jL)
      else  !if flow from upstream tunnel is small, force behaviour like a weir discharge
        Q(iL,jL)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*(0.5_wp*Tfs(iL,jL))*(max(0.0_wp,y(iL,jL)))**1.5
        A(iL,jL)=A(iL-1,jL)
        V(iL,jL)=Q(iL,jL)/A(iL,jL)
        c(iL,jL)=FindCel(iL,jL)
      end if

      if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(iR,jR)-Q(iL,jL)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)

  case(2)

    Junc(k)%Head=0.001
    converged=.false.

    Hup=1.00*max(Reach(jL)%D,Reach(jR)%D)
    Hdown=0.001

    do while(not(converged))
      Junc(k)%Head=0.5*(Hup+Hdown)
      niter=niter+1

      y(iR,jR)=max(initHfrac*Reach(jR)%D,Junc(k)%Head-0.5*Reach(jR)%dX*Reach(jR)%So)
      c(iR,jR)=FindCel(iR,jR)
      A(iR,jR)=FindArea(iR,jR)
      V(iR,jR)=Ks + y(iR,jR)*g/cs
      Q(iR,jR)=V(iR,jR)*A(iR,jR)
      if (abs(V(iR,jR)) > c(iR,jR)) then
        if ((y(iR+1,jR) < alpha1*y(iR,jR)) .or. (y(iR+2,jR) < alpha2*y(iR,jR))) then
          if (V(iR,jR) > 0) then
            V(iR,jR)=c(iR,jR)
            Q(iR,jR)=A(iR,jR)*V(iR,jR)
          else
            V(iR,jR)=-c(iR,jR)
            Q(iR,jR)=A(iR,jR)*V(iR,jR)
          end if
        else
          A(iR,jR)=A(iR+1,jR)
          V(iR,jR)=V(iR+1,jR)
          Q(iR,jR)=Q(iR+1,jR)
        end if
      else
        Q(iR,jR)=A(iR,jR)*V(iR,jR)
      end if

! PRK 3/10/2015 insert new handling of supercritical flow from LHS -- just pass values from upstream cell
!      if (V(iL-1,jL) > c(iL-1,jL)) then
!          A(iL,jL)=A(iL-1,jL)
!          y(iL,jL)=y(iL-1,jL)
!          V(iL,jL)=V(iL-1,jL)
!          Q(iL,jL)=Q(iL-1,jL)
!          c(iL,jL)=FindCel(iL,jL)
!      else
! PRK 3/10/2015
          y(iL,jL)=max(initHfrac*Reach(jL)%D,Junc(k)%Head-Zcorr+0.5*Reach(jL)%dX*Reach(jL)%So)
          c(iL,jL)=FindCel(iL,jL)
          A(iL,jL)=FindArea(iL,jL)
          V(iL,jL)=Kr - y(iL,jL)*g/cr
          Q(iL,jL)=V(iL,jL)*A(iL,jL)
          if (V(iL-1,jL) > c(iL-1,jL)) then
            Q(iL,jL)=Q(iL-1,jL)
            V(iL,jL)=Q(iL,jL)/A(iL,jL)
          else if (abs(V(iL,jL)) > c(iL,jL)) then
            if ((y(iL-1,jL) < alpha1*y(iL,jL)) .or. (y(iL-2,jL) < alpha2*y(iL,jL))) then
              if (V(iL,jL) > 0) then
                V(iL,jL)=c(iL,jL)
                Q(iL,jL)=A(iL,jL)*V(iL,jL)
              else
                V(iL,jL)=-c(iL,jL)
                Q(iL,jL)=A(iL,jL)*V(iL,jL)
              end if
            else
              A(iL,jL)=A(iL-1,jL)
              y(iL,jL)=y(iL-1,jL)
              V(iL,jL)=V(iL-1,jL)
              Q(iL,jL)=Q(iL-1,jL)
            end if
          else
            Q(iL,jL)=A(iL,jL)*V(iL,jL)
          end if
!      end if

      if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(iR,jR)-Q(iL,jL)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)

  case(3)
    
    if (Q(iL-1,jL) > 0.005*sqrt(g*(Reach(jL)%D**5))) then !either critical or supercritical
      Q(iL,jL)=Q(iL-1,jL)
      A(iL,jL)=A(iL-1,jL)
      V(iL,jL)=V(iL-1,jL)
      y(iL,jL)=HcL
!      y(iL,jL)=max(initHfrac*D(jL),y(iL-1,jL))
      c(iL,jL)=FindCel(iL,jL)
    else if (Q(iL-1,jL) > 0) then !if flow from upstream spur tunnel is small but positive, force behaviour like a weir discharge
      Q(iL,jL)=(2/3.)*sqrt(2*g)*(0.5*Tfs(iL,jL))*(max(0.0_wp,y(iL,jL)))**1.5
      A(iL,jL)=A(iL-1,jL)
      V(iL,jL)=Q(iL,jL)/A(iL,jL)
      c(iL,jL)=FindCel(iL,jL)
    else
      y(iL,jL)=max(initHfrac*Reach(jL)%D,Junc(k)%Head-Zcorr+0.5*Reach(jL)%dX*Reach(jL)%So)
      c(iL,jL)=FindCel(iL,jL)
      A(iL,jL)=FindArea(iL,jL)
      V(iL,jL)=Kr - y(iL,jL)*g/cr
      Q(iL,jL)=V(iL,jL)*A(iL,jL)
      if (V(iL-1,jL) > c(iL-1,jL)) then
        Q(iL,jL)=Q(iL-1,jL)
        V(iL,jL)=Q(iL,jL)/A(iL,jL)
      else if (abs(V(iL,jL)) > c(iL,jL)) then
        if ((y(iL-1,jL) < 1.05*y(iL,jL)) .or. (y(iL-2,jL) < 1.10*y(iL,jL))) then
          if (V(iL,jL) > 0) then
            V(iL,jL)=c(iL,jL)
            Q(iL,jL)=A(iL,jL)*V(iL,jL)
          else
            V(iL,jL)=-c(iL,jL)
            Q(iL,jL)=A(iL,jL)*V(iL,jL)
          end if
        else
          A(iL,jL)=A(iL-1,jL)
          y(iL,jL)=y(iL-1,jL)
          V(iL,jL)=V(iL-1,jL)
          Q(iL,jL)=Q(iL-1,jL)
        end if
      else
        Q(iL,jL)=A(iL,jL)*V(iL,jL)
      end if
    end if

    dVn1= dT*g/Reach(jR)%dX * ( Junc(k)%Head - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) )*abs( V(iR,jR) )/(2*g) )
    dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR) ) )

    dVn2= dT*g/Reach(jR)%dX * ( (Junc(k)%Head+dYm1/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn1/2  )*abs( V(iR,jR) + dVn1/2 )/(2*g) )
    dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR)+dVn1/2 ) )

    dVn3= dT*g/Reach(jR)%dX * ( (Junc(k)%Head+dYm2/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn2/2  )*abs( V(iR,jR) + dVn2/2 )/(2*g) )
    dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR)+dVn2/2 ) )

    dVn4= dT*g/Reach(jR)%dX * ( (Junc(k)%Head+dYm3) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn3  )*abs( V(iR,jR) + dVn3 )/(2*g) )
    dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR)+dVn3 ) )

    V(iR,jR)=V(iR,jR) + (dVn1+2*(dVn2+dVn3)+dVn4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dYm1+2*(dYm2+dYm3)+dYm4)/6.0_wp

    y(iR,jR)=Junc(k)%Head - Reach(jR)%Kup*V(iR,jR)*abs(V(iR,jR))/(2*g)
    A(iR,jR)=FindArea(iR,jR)
    c(iR,jR)=findcel(iR,jR)
    Q(iR,jR)=A(iR,jR)*V(iR,jR)

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)

  case(4)

    y(iL,jL)=max(initHfrac*Reach(jL)%D,Junc(k)%Head-Zcorr+0.5*Reach(jL)%dX*Reach(jL)%So)
    c(iL,jL)=FindCel(iL,jL)
    A(iL,jL)=FindArea(iL,jL)
    V(iL,jL)=Kr - y(iL,jL)*g/cr
    Q(iL,jL)=V(iL,jL)*A(iL,jL)
    if (V(iL-1,jL) > c(iL-1,jL)) then
      Q(iL,jL)=Q(iL-1,jL)
      V(iL,jL)=Q(iL,jL)/A(iL,jL)
    else if (abs(V(iL,jL)) > c(iL,jL)) then
      if ((y(iL-1,jL) < 1.05*y(iL,jL)) .or. (y(iL-2,jL) < 1.10*y(iL,jL))) then
        if (V(iL,jL) > 0) then
          V(iL,jL)=c(iL,jL)
          Q(iL,jL)=A(iL,jL)*V(iL,jL)
        else
          V(iL,jL)=-c(iL,jL)
          Q(iL,jL)=A(iL,jL)*V(iL,jL)
        end if
      else
        A(iL,jL)=A(iL-1,jL)
        y(iL,jL)=y(iL-1,jL)
        V(iL,jL)=V(iL-1,jL)
        Q(iL,jL)=Q(iL-1,jL)
      end if
    else
      Q(iL,jL)=A(iL,jL)*V(iL,jL)
    end if

    dVn1= dT*g/Reach(jR)%dX * ( Junc(k)%Head - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) )*abs( V(iR,jR) )/(2*g) )
    dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR) ) )

    dVn2= dT*g/Reach(jR)%dX * ( (Junc(k)%Head+dYm1/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn1/2  )*abs( V(iR,jR) + dVn1/2 )/(2*g) )
    dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR)+dVn1/2 ) )

    dVn3= dT*g/Reach(jR)%dX * ( (Junc(k)%Head+dYm2/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn2/2  )*abs( V(iR,jR) + dVn2/2 )/(2*g) )
    dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR)+dVn2/2 ) )

    dVn4= dT*g/Reach(jR)%dX * ( (Junc(k)%Head+dYm3) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn3  )*abs( V(iR,jR) + dVn3 )/(2*g) )
    dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*V(iL,jL) - A(iR,jR)*( V(iR,jR)+dVn3 ) )

    V(iR,jR)=V(iR,jR) + (dVn1+2*(dVn2+dVn3)+dVn4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dYm1+2*(dYm2+dYm3)+dYm4)/6.0_wp

    y(iR,jR)=Junc(k)%Head - Reach(jR)%Kup*V(iR,jR)*abs(V(iR,jR))/(2*g)
    A(iR,jR)=FindArea(iR,jR)
    c(iR,jR)=findcel(iR,jR)
    Q(iR,jR)=A(iR,jR)*V(iR,jR)

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)

  case(5)

    if (IsDrownedUp) then
      dVm1= dT*g/Reach(jL)%dX * (-(Junc(k)%Head-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) )*abs( V(iL,jL) )/(2*g) )
      dVn1= dT*g/Reach(jR)%dX * (Junc(k)%Head - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) )*abs( V(iR,jR) )/(2*g) )
      dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL) ) - A(iR,jR)*( V(iR,jR) ) )

      dVm2= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm1/2-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm1/2  )*abs( V(iL,jL) + dVm1/2 )/(2*g) )
      dVn2= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm1/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn1/2  )*abs( V(iR,jR) + dVn1/2 )/(2*g) )
      dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm1/2 ) - A(iR,jR)*( V(iR,jR)+dVn1/2 ) )

      dVm3= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm2/2-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm2/2  )*abs( V(iL,jL) + dVm2/2 )/(2*g) )
      dVn3= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm2/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn2/2  )*abs( V(iR,jR) + dVn2/2 )/(2*g) )
      dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm2/2 ) - A(iR,jR)*( V(iR,jR)+dVn2/2 ) )

      dVm4= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm3-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm3  )*abs( V(iL,jL) + dVm3 )/(2*g) )
      dVn4= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm3) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn3  )*abs( V(iR,jR) + dVn3 )/(2*g) )
      dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm3 ) - A(iR,jR)*( V(iR,jR)+dVn3 ) )
    
    else

      dVm1= dT*g/Reach(jL)%dX * (-HcL + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) )*abs( V(iL,jL) )/(2*g) )
      dVn1= dT*g/Reach(jR)%dX * (Junc(k)%Head - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) )*abs( V(iR,jR) )/(2*g) )
      dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL) ) - A(iR,jR)*( V(iR,jR) ) )

      dVm2= dT*g/Reach(jL)%dX * (-(HcL+dYm1/2) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm1/2  )*abs( V(iL,jL) + dVm1/2 )/(2*g) )
      dVn2= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm1/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn1/2  )*abs( V(iR,jR) + dVn1/2 )/(2*g) )
      dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm1/2 ) - A(iR,jR)*( V(iR,jR)+dVn1/2 ) )

      dVm3= dT*g/Reach(jL)%dX * (-(HcL+dYm2/2) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm2/2  )*abs( V(iL,jL) + dVm2/2 )/(2*g) )
      dVn3= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm2/2) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn2/2  )*abs( V(iR,jR) + dVn2/2 )/(2*g) )
      dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm2/2 ) - A(iR,jR)*( V(iR,jR)+dVn2/2 ) )

      dVm4= dT*g/Reach(jL)%dX * (-(HcL+dYm3) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm3  )*abs( V(iL,jL) + dVm3 )/(2*g) )
      dVn4= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm3) - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn3  )*abs( V(iR,jR) + dVn3 )/(2*g) )
      dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm3 ) - A(iR,jR)*( V(iR,jR)+dVn3 ) )

    end if
    
    V(iL,jL)=V(iL,jL) + (dVm1+2*(dVm2+dVm3)+dVm4)/6.0_wp
    V(iR,jR)=V(iR,jR) + (dVn1+2*(dVn2+dVn3)+dVn4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dYm1+2*(dYm2+dYm3)+dYm4)/6.0_wp

    y(iR,jR)=Junc(k)%Head - Reach(jR)%Kup*V(iR,jR)*abs(V(iR,jR))/(2*g)
    A(iR,jR)=FindArea(iR,jR)
    c(iR,jR)=findcel(iR,jR)
    Q(iR,jR)=A(iR,jR)*V(iR,jR)

    if (IsDrownedUp) then
      y(iL,jL)=Junc(k)%Head - Zcorr + Reach(jL)%Kdown*V(iL,jL)*abs(V(iL,jL))/(2*g)
    else
      y(iL,jL)=HcL + Reach(jL)%Kdown*V(iL,jL)*abs(V(iL,jL))/(2*g)
    end if
    
    A(iL,jL)=FindArea(iL,jL)
    c(iL,jL)=findcel(iL,jL)
    Q(iL,jL)=A(iL,jL)*V(iL,jL)

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)

  case(6)
    
    y(iR,jR)=max(initHfrac*Reach(jR)%D,Junc(k)%Head-0.5*Reach(jR)%dX*Reach(jR)%So)
    c(iR,jR)=FindCel(iR,jR)
    A(iR,jR)=FindArea(iR,jR)
    V(iR,jR)=Ks + y(iR,jR)*g/cs
    Q(iR,jR)=V(iR,jR)*A(iR,jR)
    if (abs(V(iR,jR)) > c(iR,jR)) then
      if ((y(iR+1,jR) < 1.05*y(iR,jR)) .or. (y(iR+2,jR)<1.10*y(iR,jR))) then
        if (V(iR,jR) > 0) then
          V(iR,jR)=c(iR,jR)
          Q(iR,jR)=A(iR,jR)*V(iR,jR)
        else
          V(iR,jR)=-c(iR,jR)
          Q(iR,jR)=A(iR,jR)*V(iR,jR)
        end if
      else
        A(iR,jR)=A(iR+1,jR)
        V(iR,jR)=V(iR+1,jR)
        Q(iR,jR)=Q(iR+1,jR)
      end if
    else
      Q(iR,jR)=A(iR,jR)*V(iR,jR)
    end if
! PRK (6/13/2014) alternate R-K to distinguish between drowned and freely discharging upstream pipe
    if (IsDrownedUp) then
      dVm1= dT*g/Reach(jL)%dX * (-(Junc(k)%Head-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) )*abs( V(iL,jL) )/(2*g) )
      dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL) ) - A(iR,jR)*( V(iR,jR) ) )

      dVm2= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm1/2-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm1/2  )*abs( V(iL,jL) + dVm1/2 )/(2*g) )
      dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm1/2 ) - A(iR,jR)*V(iR,jR))

      dVm3= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm2/2-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm2/2  )*abs( V(iL,jL) + dVm2/2 )/(2*g) )
      dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm2/2 ) - A(iR,jR)*V(iR,jR))

      dVm4= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm3-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm3  )*abs( V(iL,jL) + dVm3 )/(2*g) )
      dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm3 ) - A(iR,jR)*V(iR,jR))

      V(iL,jL)=V(iL,jL) + (dVm1+2*(dVm2+dVm3)+dVm4)/6.0_wp
      Junc(k)%Head=Junc(k)%Head + (dYm1+2*(dYm2+dYm3)+dYm4)/6.0_wp
! PRK 5/13/2020 experimental handling of loss term to behave better under reverse flow conditions
      if (V(iL,jL) > 0.0_wp) then
          y(iL,jL)=Junc(k)%Head - Zcorr + Reach(jL)%Kdown*V(iL,jL)*abs(V(iL,jL))/(2*g)
      else
          y(iL,jL)=max(Junc(k)%Head - Zcorr + Reach(jL)%Kdown*V(iL,jL)*abs(V(iL,jL))/(2*g),0.5_wp*(y(iL-1,jL)+Junc(k)%Head-Zcorr))
      end if
! PRK 5/13/2020 end
    else
      dVm1= dT*g/Reach(jL)%dX * (-HcL + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) )*abs( V(iL,jL) )/(2*g) )
      dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL) ) - A(iR,jR)*( V(iR,jR) ) )

      dVm2= dT*g/Reach(jL)%dX * (-(HcL+dYm1/2) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm1/2  )*abs( V(iL,jL) + dVm1/2 )/(2*g) )
      dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm1/2 ) - A(iR,jR)*V(iR,jR))

      dVm3= dT*g/Reach(jL)%dX * (-(HcL+dYm2/2) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm2/2  )*abs( V(iL,jL) + dVm2/2 )/(2*g) )
      dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm2/2 ) - A(iR,jR)*V(iR,jR))

      dVm4= dT*g/Reach(jL)%dX * (-(HcL+dYm3) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm3  )*abs( V(iL,jL) + dVm3 )/(2*g) )
      dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm3 ) - A(iR,jR)*V(iR,jR))

      V(iL,jL)=V(iL,jL) + (dVm1+2*(dVm2+dVm3)+dVm4)/6.0_wp
      Junc(k)%Head=Junc(k)%Head + (dYm1+2*(dYm2+dYm3)+dYm4)/6.0_wp

      y(iL,jL)=HcL + Reach(jL)%Kdown*V(iL,jL)*abs(V(iL,jL))/(2*g)
    end if
! PRK (6/13/2014) end
      
    A(iL,jL)=FindArea(iL,jL)
    c(iL,jL)=findcel(iL,jL)
    Q(iL,jL)=A(iL,jL)*V(iL,jL)

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)
    
  case(7)
      
! we follow here the calculation sequence outlined in Wylie and Streeter
    B1=amax/(g*Reach(jL)%Apipe)
    B2=amax/(g*Reach(jR)%Apipe)

! conversion to f (Darcy) from n (Manning)

    f1=12.69920844327177*g*Reach(jL)%n**2/(Reach(jL)%D)**0.3333
    f2=12.69920844327177*g*Reach(jR)%n**2/(Reach(jR)%D)**0.3333

    Rf1=f1*Reach(jL)%dX/(2*g*Reach(jL)%D*Reach(jL)%Apipe**2)
    Rf2=f2*Reach(jR)%dX/(2*g*Reach(jR)%D*Reach(jR)%Apipe**2)

! assume one C+ characteristic line arriving into the
! junction, and one C- characteristic line leaving

    CP1=y(iL-1,jL) + z(iL-1,jL) + B1*Q(iL-1,jL)
    CM2=y(iR+1,jR) + z(iR+1,jR) - B2*Q(iR+1,jR)

    BP1=B1 + Rf1*abs(Q(iL-1,jL))
    BM2=B2 - Rf2*abs(Q(iR+1,jR))

    SC=CP1/BP1 + CM2/BM2
    SB=1/BP1 + 1/BM2

! overall pressure is uniform for all cells and junction, the value
! calculated by the expression. I am assuming that Qn=0
  
    Htotal=SC/SB + (1/SB)*0
    Junc(k)%Head=Htotal - Junc(k)%Elev

    Q(iL,jL)= (CP1 - Htotal)/BP1
    Q(iR,jR)= (Htotal - CM2)/BM2 

    A(iL,jL)=FindArea(iL,jL)
    V(iL,jL)=Q(iL,jL)/A(iL,jL)

    A(iR,jR)=FindArea(iR,jR)
    V(iR,jR)=Q(iR,jR)/A(iR,jR)
    
    y(iL,jL)=Junc(k)%Head-Zcorr
    y(iR,jR)=Junc(k)%Head
    
! PRK 10/4/2012 insert corrected handling of boundary cells instead of call to TPACalcBC (which does not handle negative depth although it is legal in this case)    
    i=iL
    j=jL
    if (A(i,j) >= Reach(j)%Apipe) then
      LeftCase=1
   	  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
!    else if (((Aold(i,j) >= A(i,j)) .and. ((Aold(i-1,j) > A(i,j)*0.9999) .or. (vacuum(i-1,j)))) .or. &
!             ((vacuum(i,j)) .and. ((Aold(i-1,j) > A(i,j)*0.9999) .or. (vacuum(i-1,j))))) then
    else if (((Aold(i,j) >= Reach(j)%Apipe) .and. ((Aold(i-1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i-1,j)))) .or. &
             ((vacuum(i,j)) .and. ((Aold(i-1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i-1,j))))) then
      LeftCase=2
      hs(i,j)=y(i,j)-Reach(j)%D
      A(i,j)=Reach(j)%Apipe*(1.0_wp + g*hs(i,j)/amax**2)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
    else
      LeftCase=3
      hs(i,j)=0.0_wp
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
	end if
	if ((A(i,j) > 0.0_dp) .and. (Tfs(i,j) > 0.0_dp)) then
	  c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
    end if
	if (A(i,j) > 0.0_dp) then
      V(i,j)=Q(i,j)/A(i,j)
	end if
	if (hs(i,j) < 0.0_dp) then
	  vacuum(i,j)=.true.
	else
	  vacuum(i,j)=.false.
    end if

    i=iR
    j=jR

    if (A(i,j) >= Reach(j)%Apipe) then
      RightCase=4
   	  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
!    else if (((Aold(i,j) >= A(i,j)) .and. ((Aold(i+1,j) > A(i,j)*0.9999) .or. (vacuum(i+1,j)))) .or. &
!             ((vacuum(i,j)) .and. ((Aold(i+1,j) > A(i,j)*0.9999) .or. (vacuum(i+1,j))))) then
    else if (((Aold(i,j) >= Reach(j)%Apipe) .and. ((Aold(i+1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i+1,j)))) .or. &
             ((vacuum(i,j)) .and. ((Aold(i+1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i+1,j))))) then
      RightCase=5
      hs(i,j)=y(i,j)-Reach(j)%D
      A(i,j)=Reach(j)%Apipe*(1.0_wp + g*hs(i,j)/amax**2)
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
    else
      RightCase=6
      hs(i,j)=0.0_wp
      hc(i,j)=FindCentroid(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Rh(i,j)=FindRh(i,j)
	end if
	if ((A(i,j) > 0.0_dp) .and. (Tfs(i,j) > 0.0_dp)) then
	  c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
    end if
	if (A(i,j) > 0.0_dp) then
      V(i,j)=Q(i,j)/A(i,j)
	end if
	if (hs(i,j) < 0.0_dp) then
	  vacuum(i,j)=.true.
	else
	  vacuum(i,j)=.false.
    end if
! PRK 10/4/2012 end 

  case(8)
!if this is the first time that this routine is used, we need to initialize the trapped air chamber
    if (Junc(k)%InitializeHch) then
      Junc(k)%Hch=atm
      Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
      call FindCurrentChamberVolume
      Junc(k)%InitializeHch=.false.
    end if

    dVm1= dT*g/Reach(jL)%dX * (-(Junc(k)%Head-Zcorr+ Junc(k)%Hch-atm) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) )*abs( V(iL,jL) )/(2*g) )
    dVn1= dT*g/Reach(jR)%dX * (Junc(k)%Head + Junc(k)%Hch-atm- y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) )*abs( V(iR,jR) )/(2*g) )
    dYm1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL) ) - A(iR,jR)*( V(iR,jR) ) )
    Pdrop=(Junc(k)%Hch-atm)/(Junc(k)%Hch)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    dHch1=-npoly*(Junc(k)%Hch)*(-dYm1*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber)

    dVm2= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm1/2-Zcorr + Junc(k)%Hch-atm+dHch1/2) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm1/2  )*abs( V(iL,jL) + dVm1/2 )/(2*g) )
    dVn2= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm1/2) + Junc(k)%Hch-atm+dHch1/2 - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn1/2  )*abs( V(iR,jR) + dVn1/2 )/(2*g) )
    dYm2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm1/2 ) - A(iR,jR)*( V(iR,jR)+dVn1/2 ) )
    Pdrop=(Junc(k)%Hch+dHch1/2-atm)/(Junc(k)%Hch+dHch1/2)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)/2
    dHch2=-npoly*(Junc(k)%Hch+dHch1/2)*(-(dYm1/2)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dYm1/2*Junc(k)%CurrentArea)

    dVm3= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm2/2 + Junc(k)%Hch-atm+dHch2/2-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm2/2  )*abs( V(iL,jL) + dVm2/2 )/(2*g) )
    dVn3= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm2/2) + Junc(k)%Hch-atm+dHch2/2 - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn2/2  )*abs( V(iR,jR) + dVn2/2 )/(2*g) )
    dYm3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm2/2 ) - A(iR,jR)*( V(iR,jR)+dVn2/2 ) )
    Pdrop=(Junc(k)%Hch+dHch2/2-atm)/(Junc(k)%Hch+dHch2/2)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)/2
    dHch3=-npoly*(Junc(k)%Hch+dHch2/2)*(-(dYm2/2)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dYm2/2*Junc(k)%CurrentArea)

    dVm4= dT*g/Reach(jL)%dX * (-(Junc(k)%Head+dYm3 + Junc(k)%Hch-atm+dHch3-Zcorr) + y(iL-1,jL) - Reach(jL)%Kdown* ( V(iL,jL) + dVm3  )*abs( V(iL,jL) + dVm3 )/(2*g) )
    dVn4= dT*g/Reach(jR)%dX * ((Junc(k)%Head+dYm3) + Junc(k)%Hch-atm+dHch3 - y(iR+1,jR) - Reach(jR)%Kup* ( V(iR,jR) + dVn3  )*abs( V(iR,jR) + dVn3 )/(2*g) )
    dYm4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(iL,jL)*( V(iL,jL)+dVm3 ) - A(iR,jR)*( V(iR,jR)+dVn3 ) )
    Pdrop=(Junc(k)%Hch+dHch3-atm)/(Junc(k)%Hch+dHch3)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    dHch4=-npoly*(Junc(k)%Hch+dHch3)*(-(dYm3)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dYm3*Junc(k)%CurrentArea)

    V(iL,jL)=V(iL,jL) + (dVm1+2*(dVm2+dVm3)+dVm4)/6.0_wp
    V(iR,jR)=V(iR,jR) + (dVn1+2*(dVn2+dVn3)+dVn4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dYm1+2*(dYm2+dYm3)+dYm4)/6.0_wp
    Junc(k)%Hch=Junc(k)%Hch + (dHch1+2*(dHch2+dHch3)+dHch4)/6.0_wp
    Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
    call FindCurrentChamberVolume
    dVolChamber=Junc(k)%VolChamber-OldVolChamber
    Pdrop=(Junc(k)%Hch-atm)/(Junc(k)%Hch)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    Junc(k)%Qair=VolExhaust/dT

    y(iR,jR)=Junc(k)%Head - Reach(jR)%Kup*V(iR,jR)*abs(V(iR,jR))/(2*g)+Junc(k)%Hch-atm
    A(iR,jR)=FindArea(iR,jR)
    c(iR,jR)=findcel(iR,jR)
    Q(iR,jR)=A(iR,jR)*V(iR,jR)

    y(iL,jL)=Junc(k)%Head - Zcorr + Reach(jL)%Kdown*V(iL,jL)*abs(V(iL,jL))/(2*g)+Junc(k)%Hch-atm
    A(iL,jL)=FindArea(iL,jL)
    c(iL,jL)=findcel(iL,jL)
    Q(iL,jL)=A(iL,jL)*V(iL,jL)

    call TPACalcBC(iR,jR,JunctionCase)
    call TPACalcBC(iL,jL,JunctionCase)
    
  case default

     write (*,*) 'indeterminate case in two-way dropshaft routine for junction',k
     stop
     
end select

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,i4,19es16.8,i8,es16.8)') T,JunctionCase,Junc(k)%Head,Junc(k)%Inflow,Junc(k)%Outflow,y(iL-1,jL),y(iL,jL),y(iR,jR),y(iR+1,jR),Q(iL-1,jL), &
	  Q(iL,jL),Q(iR,jR),Q(iR+1,jR),V(iL-1,jL),V(iL,jL),V(iR,jR),V(iR+1,jR),c(iL-1,jL),c(iL,jL),c(iR,jR),c(iR+1,jR),niter,HcL
      if (Junc(k)%Option == -8) write (6444,'(f11.4,i4,14f12.6)') T,JunctionCase,Junc(k)%Head,Junc(k)%Hch,dYm1,dHch1,dYm2,dHch2,dYm3,dHch3,dYm4,dHch4,Junc(k)%Lch,Junc(k)%VolChamber,Junc(k)%Qair,VolExhaust
      write (6445,*) T,HcL,IsWeirFlow,IsDrownedUp
	end if
  end do
end if

99 continue

return

end subroutine TwoWayDropshaft

