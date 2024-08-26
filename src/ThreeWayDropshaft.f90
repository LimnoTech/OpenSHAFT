subroutine ThreeWayDropshaft

use GlobalVariables
use GlobalFunctions

implicit none

integer :: LeftCase,RightCase,i1,i2,i3,j1,j2,j3,niter,kk
real(wp) :: D1,D2,D3,JElev,Z1,Z2,Z3,Crown1,Crown2,Crown3,V1old,V3old,Qlim1,Qlim2,Qlim3,Qlim4,Pdrop,VolExhaust,Htotal
real(wp) :: Qtemp,Vtemp,Atemp,ytemp,ytemp1,ytemp2,OldVolChamber,HchO,Hc1,H1O,Hc3,H3O,dV11,dV12,dV13,dV14,dV21,dV22,dV23,dV24,dV31,dV32,dV33,dV34,dY1,dY2,dY3,dY4,dHch1,dHch2,dHch3,dHch4
real(wp) :: r1,Vr1,cr1,Sfr1,yr1,kr1,r2,Vs2,cs2,Sfs2,ys2,ks2,r3,Vr3,cr3,Sfr3,yr3,kr3,Hconj1,Hup,Hdown,Rf1,Rf2,Rf3,B1,B2,B3,f1,f2,f3,CP1,CM2,CP3,BP1,BM2,BP3,SB,SC
logical :: IsSurch1,IsSurchMain,IsSurch3,IsWeir1,IsWeir3,IsNoShaft,IsDrowned1,IsDrowned3,converged

if (dT == 0.0_wp) then
  go to 99
end if

i1=Junc(k)%ReachCell(1)    !first reach cell connected downstream witn this junction - main tunnel
i3=Junc(k)%ReachCell(2)    !second reach cell connected downstream witn this junction - spur tunnel
i2=Junc(k)%ReachCell(3)    !first reach cell connected upstream witn this junction - main tunnel

j1=Junc(k)%ReachNo(1)      !first reach cell connected downstream witn this junction - main tunnel
j3=Junc(k)%ReachNo(2)      !second reach cell connected downstream witn this junction - spur tunnel
j2=Junc(k)%ReachNo(3)      !first reach cell connected upstream witn this junction - main tunnel

D1=Reach(j1)%D         !this is upstream tunnel 1
D2=Reach(j2)%D         !this is MAIN TUNNEL DOWNSTREAM
D3=Reach(j3)%D         !this is upstream tunnel 3

Aold(i1,j1)=A(i1,j1)
Aold(i2,j2)=A(i2,j2)
Aold(i3,j3)=A(i3,j3)

JElev=Junc(k)%Elev
Junc(k)%HeadOld=Junc(k)%Head
OldVolChamber=Junc(k)%VolChamber
HchO=Junc(k)%Hch
niter=0
LeftCase=0
RightCase=0

Z1=Junc(k)%ReachElev(1)-JElev
Z2=Junc(k)%ReachElev(3)-JElev
Z3=Junc(k)%ReachElev(2)-JElev

Crown1=Z1+D1
Crown2=Z2+D2
Crown3=Z3+D3

V1old=V(i1,j1)
V3old=V(i3,j3)

if (IsHGLInit) then
  if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0.0_dp) .and. (abs(Q(i1-1,j1)) < 4e-4*sqrt(g*Reach(j1)%D**5)) .and. (abs(Q(i2+1,j2)) < 4e-4*sqrt(g*Reach(j2)%D**5)) &
      .and. (abs(Q(i3-1,j3)) < 4e-4*sqrt(g*Reach(j3)%D**5))) then
    go to 99
  end if
end if

if (not(IsHGLInit)) then
  if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0.0_dp) .and. (abs(Q(i1-1,j1)) < 4e-5*sqrt(g*Reach(j1)%D**5)) .and. (abs(Q(i2+1,j2)) < 4e-5*sqrt(g*Reach(j2)%D**5)) &
      .and. (abs(Q(i3-1,j3)) < 4e-5*sqrt(g*Reach(j3)%D**5))) then
    go to 99
  end if
end if

Hc1=FindCriticalDepth(i1-1,j1)
H1O=min(Hc1,y(i1-1,j1))

if (V(i1-1,j1) > c(i1-1,j1)) then
  Hconj1=y(i1-1,j1)*0.5_wp*(sqrt(1.0_wp +8.0_wp*V(i1-1,j1)*V(i1-1,j1)/(g*y(i1-1,j1)))-1.0_wp)
else
  Hconj1=0.0_wp
end if

Hc3=FindCriticalDepth(i3-1,j3)
H3O=min(Hc3,y(i3-1,j3))

! find limiting flow for convergence criterion

Qlim1=4e-4*sqrt(g*Reach(j1)%D**5)
Qlim2=4e-4*sqrt(g*Reach(j2)%D**5)
Qlim3=4e-4*sqrt(g*Reach(j3)%D**5)
Qlim=min(max(Qlim1,abs(Q(i1-1,j1))),max(Qlim2,abs(Q(i2+1,j2))),max(Qlim3,abs(Q(i3-1,j3))))

! set values for characteristic equations

r1=dT/Reach(j1)%dX
Vr1=(V(i1,j1) + r1*(-V(i1,j1)*c(i1-1,j1)+c(i1,j1)*V(i1-1,j1))) / (1.0_wp +r1*(V(i1,j1)-V(i1-1,j1)+c(i1,j1)-c(i1-1,j1)))
cr1=(c(i1,j1) + r1*Vr1*(c(i1-1,j1) - c(i1,j1))) / (1.0_wp + r1*(c(i1,j1) - c(i1-1,j1)))
Sfr1=Reach(j1)%n**2*Vr1*abs(Vr1)/Rh(i1,j1)**1.333
yr1=y(i1-1,j1) + r1*(Vr1+cr1)*(y(i1,j1) - y(i1-1,j1))
Kr1=Vr1 + g*yr1/cr1 - g*(Sfr1-Reach(j1)%So)*dT

r2=dT/Reach(j2)%dX
Vs2=(V(i2,j2) + r2*(c(i2,j2)*V(i2+1,j2)-c(i2+1,j2)*V(i2,j2))) / (1.0_wp + r2*(-V(i2,j2)+V(i2+1,j2)+c(i2,j2)-c(i2+1,j2)))
cs2=(c(i2,j2) + r2*Vs2*(c(i2,j2) - c(i2+1,j2))) / (1.0_wp + r2*(c(i2,j2) - c(i2+1,j2)))
Sfs2=Reach(j2)%n**2*Vs2*abs(Vs2)/Rh(i2,j2)**1.333
ys2=y(i2,j2) - r2*(Vs2-cs2)*(y(i2+1,j2) - y(i2,j2))
Ks2=Vs2 - g*ys2/cs2 - g*(Sfs2-Reach(j2)%So)*dT

r3=dT/Reach(j3)%dX
Vr3=(V(i3,j3) + r3*(-V(i3,j3)*c(i3-1,j3)+c(i3,j3)*V(i3-1,j3))) / (1.0_wp +r3*(V(i3,j3)-V(i3-1,j3)+c(i3,j3)-c(i3-1,j3)))
cr3=(c(i3,j3) + r3*Vr3*(c(i3-1,j3) - c(i3,j3))) / (1.0_wp + r3*(c(i3,j3) - c(i3-1,j3)))
Sfr3=Reach(j3)%n**2*Vr3*abs(Vr3)/Rh(i3,j3)**1.333
yr3=y(i3-1,j3) + r3*(Vr3+cr3)*(y(i3,j3) - y(i3-1,j3))
Kr3=Vr3 + g*yr3/cr3 - g*(Sfr3-Reach(j3)%So)*dT

! establish criteria for case selection

! PRK 11/6/2014 changing criteria for IsWeir1 to exclude cases of reverse flow
!if (JuncHeadOld(k) < Hc1+Z1) then
if (V(i1-1,j1) > c(i1-1,j1)) then
  if ((Junc(k)%HeadOld < Hconj1+Z1) .and. (Q(i1-1,j1) > 0.0)) then
    IsWeir1=.true.
  else
    IsWeir1=.false.
  end if
else if ((Junc(k)%HeadOld < Hc1+Z1) .and. (Q(i1-1,j1) > 0.0)) then
  IsWeir1=.true.
else
  IsWeir1=.false.
end if

if (Junc(k)%HeadOld < Hc3+Z3) then
  IsWeir3=.true.
else
  IsWeir3=.false.
end if

!if (Junc(k)HeadOld < 0.999*D(j1)+Z1) then
! PRK 11/7/2014 enforce surcharge case for matching main tunnel crowns
if (Crown1 == Crown2) then
!  if ((Junc(k)HeadOld < 0.999*D(j1)+Z1) .and. (hs(i1,j1) == 0.0_wp) .and. (hs(i2,j2) == 0.0_wp)) then
  if ((Junc(k)%HeadOld < 0.999*Reach(j1)%D+Z1) .and. ((hs(i1-1,j1) == 0.0_wp) .or. (hs(i2+1,j2) == 0.0_wp))) then
    IsSurch1=.false.
    IsSurchMain=.false.
  else
    IsSurch1=.true.
    IsSurchMain=.true.
  end if
else if (V(i1-1,j1) > c(i1-1,j1)) then
  if (Junc(k)%HeadOld < alpha4*Reach(j1)%D+Z1) then
    IsSurch1=.false.
  else
    IsSurch1=.true.
  end if
else if (y(i1,j1) < 0.999*Reach(j1)%D) then
!if (y(i1,j1) < 0.999*D(j1)) then
  IsSurch1=.false.
else
  IsSurch1=.true.
end if

!if (JuncHeadOld(k) < 0.999*D(j3)+Z3) then
!PRK 11/7/2014 add special logic for all matching crown elevations
if ((Crown1 == Crown2) .and. (Crown3 == Crown2)) then
  if ((Junc(k)%HeadOld < 0.999*Reach(j3)%D+Z3) .and. (hs(i3,j3) == 0.0_wp)) then
    IsSurch3=.false.
  else
    IsSurch3=.true.
  end if
! PRK 11/7/2014 now continue with the other stuff
else if (V(i3-1,j3) > c(i3-1,j3)) then
  if (Junc(k)%HeadOld < alpha4*Reach(j3)%D+Z3) then
    IsSurch3=.false.
  else
    IsSurch3=.true.
  end if
else if (y(i3,j3) < 0.999*Reach(j3)%D) then
!if (y(i3,j3) < 0.999*D(j3)) then
  IsSurch3=.false.
else
  IsSurch3=.true.
end if

! PRK 7/31/2013 check if spur is actually drowned or just surcharged at downstream end
! PRK 3/18/2014 modify for case when spur is not offset
if (Z3 == 0.0) then
  if (IsSurch3) then
    IsDrowned3=.true.
  else
    IsDrowned3=.false.
  end if
!else if (JuncHeadOld(k) < 0.999*(D(j3)+Z3)) then
else if (Junc(k)%HeadOld < Z3+Hc3) then
  IsDrowned3=.false.
else
  IsDrowned3=.true.
end if
! PRK 7/31/2013 end
! PRK 12/31/2014 enforce drowned spur whenever flow reverses
if (Q(i3,j3) < 0.0) IsDrowned3=.true.
! PRK 12/31/2014 end

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

if (Crown2 /= Crown1) then
  if (Junc(k)%HeadOld < 0.999*Reach(j2)%D+Z2) then
  !if (y(i2,j2) < 0.999*D(j2)) then
    IsSurchMain=.false.
  else
    IsSurchMain=.true.
  end if
end if

if (Junc(k)%CurrentArea == 0.0_wp) then
  IsNoShaft=.true.
  Junc(k)%CurrentArea=0.25*pi*D2**2
else
  IsNoShaft=.false.
end if

! PRK 11/7/2014 test for subatmospheric pressure in blind junctions
if (IsNoShaft) then
  if (vacuum(i1-1,j1)) IsSurch1=.true.
  if (vacuum(i2+1,j2)) IsSurchMain=.true.
  if (vacuum(i3-1,j3)) IsSurch3=.true.
end if
! PRK 11/7/2014

if ((IsWeir1) .and. (IsWeir3) .and. (not(IsSurchMain))) JunctionCase=1
if ((IsWeir1) .and. (not(IsWeir3)) .and. (not(IsSurchMain))) JunctionCase=2
if ((IsWeir1) .and. (IsSurch3) .and. (not(IsSurchMain))) JunctionCase=3
if ((not(IsWeir1)) .and. (IsWeir3) .and. (not(IsSurchMain))) JunctionCase=4
if ((not(IsWeir1)) .and. (not(IsWeir3)) .and. (not(IsSurchMain))) JunctionCase=5
if ((not(IsWeir1)) .and. (IsSurch3) .and. (not(IsSurchMain))) JunctionCase=6
if ((IsSurch1) .and. (IsWeir3) .and. (not(IsSurchMain))) JunctionCase=7
if ((IsSurch1) .and. (not(IsWeir3)) .and. (not(IsSurchMain))) JunctionCase=8
if ((IsSurch1) .and. (IsSurch3) .and. (not(IsSurchMain))) JunctionCase=9
if ((IsWeir1) .and. (IsWeir3) .and. (IsSurchMain)) JunctionCase=10
if ((IsWeir1) .and. (not(IsWeir3)) .and. (IsSurchMain)) JunctionCase=11
if ((IsWeir1) .and. (IsSurch3) .and. (IsSurchMain)) JunctionCase=12
if ((not(IsWeir1)) .and. (IsWeir3) .and. (IsSurchMain)) JunctionCase=13
if ((not(IsWeir1)) .and. (not(IsWeir3)) .and. (IsSurchMain)) JunctionCase=14
if ((not(IsWeir1)) .and. (IsSurch3) .and. (IsSurchMain)) JunctionCase=15
if ((IsSurch1) .and. (IsWeir3) .and. (IsSurchMain)) JunctionCase=16
if ((IsSurch1) .and. (not(IsWeir3)) .and. (IsSurchMain)) JunctionCase=17
if ((IsSurch1) .and. (IsSurch3) .and. (IsSurchMain) .and. (not(IsNoShaft))) JunctionCase=18
if ((IsSurch1) .and. (IsSurch3) .and. (IsSurchMain) .and. (IsNoShaft)) JunctionCase=19
if ((IsSurch1) .and. (IsSurch3) .and. (IsSurchMain) .and. (Junc(k)%Option == -8)) JunctionCase=20

! PRK 4/20/2015 depressurize chamber if necessary
if ((not(Junc(k)%InitializeHch)) .and. (JunctionCase /= 20)) then
  Junc(k)%Hch=atm
  Junc(k)%Lch=0.0_wp
  Junc(k)%VolChamber=0.0_wp
  Junc(k)%InitializeHch=.true.
  Junc(k)%Qair=0.0_wp
end if
! PRK 4/20/2015 end

select case(JunctionCase)
  
  case(1)
  ThreeWayCases(k,1)=ThreeWayCases(k,1)+1

  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    if (V(i1-1,j1) > c(i1-1,j1)) then
      y(i1,j1)=y(i1-1,j1)
    else
      y(i1,j1)=Hc1
    end if
    if (vacuum(i1-1,j1)) then
      A(i1,j1)=FindArea(i1,j1)
      V(i1,j1)=Q(i1,j1)/A(i1,j1)
    end if
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

  if  (Q(i3-1,j3) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)) then
    Q(i3,j3)=Q(i3-1,j3)
    A(i3,j3)=A(i3-1,j3)
    V(i3,j3)=V(i3-1,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      y(i3,j3)=y(i3-1,j3)
    else
      y(i3,j3)=Hc3
    end if
    if (vacuum(i3-1,j3)) then
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    end if
    c(i3,j3)=FindCel(i3,j3)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i3,j3)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)
    V(i3,j3)=c(i3-1,j3)
    A(i3,j3)=Q(i3,j3)/V(i3,j3)
    c(i3,j3)=c(i3-1,j3)
    y(i3,j3)=FindDepth(i3,j3)
  end if

!    bisection to find pressure head at junction
  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
    Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)
!         // do not admit the characteristic eq. based velocity if flow is supercritical. Make it critical then
!         // attempt to account for the returning bore - let flow be influenced by downstream cell via characteristic equation
    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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
  ThreeWayCases(k,2)=ThreeWayCases(k,2)+1

  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    if (V(i1-1,j1) > c(i1-1,j1)) then
      y(i1,j1)=y(i1-1,j1)
    else
      y(i1,j1)=Hc1
    end if
    if (vacuum(i1-1,j1)) then
      A(i1,j1)=FindArea(i1,j1)
      V(i1,j1)=Q(i1,j1)/A(i1,j1)
    end if
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

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)
!         // do not admit the characteristic eq. based velocity if flow is supercritical. Make it critical then
!         // attempt to account for the returning bore - let flow be influenced by downstream cell via characteristic equation
    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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

    y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3+Reach(j3)%dX/2*Reach(j3)%So)
    c(i3,j3)=FindCel(i3,j3)
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
    Q(i3,j3)=V(i3,j3)*A(i3,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      Q(i3,j3)=Q(i3-1,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    else if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3-1,j3) < 1.05*y(i3,j3)) .or. (y(i3-2,j3) < 1.10*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3-1,j3)
        V(i3,j3)=V(i3-1,j3)
        Q(i3,j3)=Q(i3-1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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
  ThreeWayCases(k,3)=ThreeWayCases(k,3)+1

  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    if (V(i1-1,j1) > c(i1-1,j1)) then
      y(i1,j1)=y(i1-1,j1)
    else
      y(i1,j1)=Hc1
    end if
    if (vacuum(i1-1,j1)) then
      A(i1,j1)=FindArea(i1,j1)
      V(i1,j1)=Q(i1,j1)/A(i1,j1)
    end if
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
    Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)
!         // do not admit the characteristic eq. based velocity if flow is supercritical. Make it critical then
!         // attempt to account for the returning bore - let flow be influenced by downstream cell via characteristic equation
    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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

    if (IsDrowned3) then
      dV31=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      V(i3,j3)=V3old + (dV31+2*dV32+2*dV33+dV34)/6.0_wp
	  y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
    else
      dV31=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      V(i3,j3)=V3old + (dV31+2*dV32+2*dV33+dV34)/6.0_wp
	  y(i3,j3)=max(Hc3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
    end if
    
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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
  ThreeWayCases(k,4)=ThreeWayCases(k,4)+1

  if (Q(i3-1,j3) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)) then
    Q(i3,j3)=Q(i3-1,j3)
    A(i3,j3)=A(i3-1,j3)
    V(i3,j3)=V(i3-1,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      y(i3,j3)=y(i3-1,j3)
    else
      y(i3,j3)=Hc3
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    end if
    if (vacuum(i3-1,j3)) then
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    end if
    c(i3,j3)=FindCel(i3,j3)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i3,j3)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5)
    V(i3,j3)=c(i3-1,j3)
    A(i3,j3)=Q(i3,j3)/V(i3,j3)
    c(i3,j3)=c(i3-1,j3)
    y(i3,j3)=FindDepth(i3,j3)
  end if

  Qtemp=Q(i3,j3)
  Vtemp=V(i3,j3)
  Atemp=A(i3,j3)
  ytemp=y(i3,j3)
  
  
!    Find junction head using bisection
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

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)

    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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

  case(5)
  ThreeWayCases(k,5)=ThreeWayCases(k,5)+1

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
      if ((y(i1-1,j1) < alpha1*y(i1,j1)) .or. (y(i1-2,j1) < alpha2*y(i1,j1))) then
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

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)

    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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

    y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3+Reach(j3)%dX/2*Reach(j3)%So)
    c(i3,j3)=FindCel(i3,j3)
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
    Q(i3,j3)=V(i3,j3)*A(i3,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      Q(i3,j3)=Q(i3-1,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    else if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3-1,j3) < 1.05*y(i3,j3)) .or. (y(i3-2,j3) < 1.10*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3-1,j3)
        V(i3,j3)=V(i3-1,j3)
        Q(i3,j3)=Q(i3-1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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

  case(6)
  ThreeWayCases(k,6)=ThreeWayCases(k,6)+1

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

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)

    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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
! PRK (7/31/2013) inserting alternate R-K formulation to distinguish between surcharged and drowned spur
    if (IsDrowned3) then
      dV31=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )

   	  y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    else
      dV31=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )

      y(i3,j3)=max(Hc3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
      
    end if

    V(i3,j3)=V3old + (dV31+2*dV32+2*dV33+dV34)/6.0_wp
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)
    
    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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
  
  case(7)
  ThreeWayCases(k,7)=ThreeWayCases(k,7)+1

  if  (Q(i3-1,j3) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)) then
    Q(i3,j3)=Q(i3-1,j3)
    A(i3,j3)=A(i3-1,j3)
    V(i3,j3)=V(i3-1,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      y(i3,j3)=y(i3-1,j3)
    else
      y(i3,j3)=Hc3
    end if
    if (vacuum(i3-1,j3)) then
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    end if
    c(i3,j3)=FindCel(i3,j3)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i3,j3)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)
    V(i3,j3)=c(i3-1,j3)
    A(i3,j3)=Q(i3,j3)/V(i3,j3)
    c(i3,j3)=c(i3-1,j3)
    y(i3,j3)=FindDepth(i3,j3)
  end if

  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
    Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
    c(i2,j2)=FindCel(i2,j2)
    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Ks2 + y(i2,j2)*g/cs2
    Q(i2,j2)=V(i2,j2)*A(i2,j2)
!         // do not admit the characteristic eq. based velocity if flow is supercritical. Make it critical then
!         // attempt to account for the returning bore - let flow be influenced by downstream cell via characteristic equation
    if (abs(V(i2,j2)) > c(i2,j2)) then
      if ((y(i2+1,j2) < alpha1*y(i2,j2)) .or. (y(i2+2,j2) < alpha2*y(i2,j2))) then
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

    if (IsDrowned1) then
      dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
! PRK (11/11/2021) fix minor loss expression for reverse flow
      if (V(i1,j1) > 0.0_wp) then
        y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
      else
        y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),0.5_wp*(y(i1-1,j1)+Junc(k)%Head-Z1))
      end if
! PRK (11/11/2021) end      
    else
      dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
      y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    end if
    
    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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
  ThreeWayCases(k,8)=ThreeWayCases(k,8)+1

  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(Reach(j1)%D+Z1,Reach(j3)%D+Z3)
  Hdown=0.001

  do while(not(converged))
    Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3+Reach(j3)%dX/2*Reach(j3)%So)
    c(i3,j3)=FindCel(i3,j3)
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
    Q(i3,j3)=V(i3,j3)*A(i3,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      Q(i3,j3)=Q(i3-1,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    else if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3-1,j3) < 1.05*y(i3,j3)) .or. (y(i3-2,j3) < 1.10*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3-1,j3)
        V(i3,j3)=V(i3-1,j3)
        Q(i3,j3)=Q(i3-1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
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

    if (IsDrowned1) then
      dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
! PRK (11/11/2021) fix minor loss expression for reverse flow
      if (V(i1,j1) > 0.0_wp) then
        y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
      else
        y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),0.5_wp*(y(i1-1,j1)+Junc(k)%Head-Z1))
      end if
! PRK (11/11/2021) end      
    else
      dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
      y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    end if

    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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

  case(9)
  ThreeWayCases(k,9)=ThreeWayCases(k,9)+1

  Junc(k)%Head=0.001
  converged=.false.

  Hup=max(D1+Z1,D2,D3+Z3)
  Hdown=0.001

  do while(not(converged))
    Junc(k)%Head=0.5*(Hup+Hdown)
    niter=niter+1

    y(i2,j2)=max(initHfrac*Reach(j2)%D,Junc(k)%Head-Reach(j2)%dX/2*Reach(j2)%So)
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

    if (IsDrowned1) then
      dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )

	  y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    
    else
      dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )

	  y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    
    end if    
    
    V(i1,j1)=V1old + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    if (IsDrowned3) then
      dV31=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )

   	  y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    else
      dV31=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )

      y(i3,j3)=max(Hc3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
      
    end if
   
    V(i3,j3)=V3old + (dV31+2*dV32+2*dV33+dV34)/6.0_wp
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

    if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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
  ThreeWayCases(k,10)=ThreeWayCases(k,10)+1

  if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
    Q(i1,j1)=Q(i1-1,j1)
    A(i1,j1)=A(i1-1,j1)
    V(i1,j1)=V(i1-1,j1)
    if (V(i1-1,j1) > c(i1-1,j1)) then
      y(i1,j1)=y(i1-1,j1)
    else
      y(i1,j1)=Hc1
    end if
    if (vacuum(i1-1,j1)) then
      A(i1,j1)=FindArea(i1,j1)
      V(i1,j1)=Q(i1,j1)/A(i1,j1)
    end if
    c(i1,j1)=FindCel(i1,j1)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
    V(i1,j1)=c(i1-1,j1)
    A(i1,j1)=Q(i1,j1)/V(i1,j1)
    c(i1,j1)=c(i1-1,j1)
    y(i1,j1)=FindDepth(i1,j1)
  end if

  if  (Q(i3-1,j3) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)) then
    Q(i3,j3)=Q(i3-1,j3)
    A(i3,j3)=A(i3-1,j3)
    V(i3,j3)=V(i3-1,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      y(i3,j3)=y(i3-1,j3)
    else
      y(i3,j3)=Hc3
    end if
    if (vacuum(i3-1,j3)) then
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    end if
    c(i3,j3)=FindCel(i3,j3)
  else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
    Q(i3,j3)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)
    V(i3,j3)=c(i3-1,j3)
    A(i3,j3)=Q(i3,j3)/V(i3,j3)
    c(i3,j3)=c(i3-1,j3)
    y(i3,j3)=FindDepth(i3,j3)
  end if

  dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
  dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

  dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
  dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

  dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
  dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

  dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
  dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

  V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
  Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

  y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
  A(i2,j2)=FindArea(i2,j2)
  c(i2,j2)=FindCel(i2,j2)
  Q(i2,j2)=A(i2,j2)*V(i2,j2)

  case(11)
    
    ThreeWayCases(k,11)=ThreeWayCases(k,11)+1

    y(i3,j3)=Junc(k)%Head-abs(Z3) !assume that the previous value of the junction head will be the final depth for reach 1
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
    c(i3,j3)=FindCel(i3,j3)
    if ((V(i3,j3) > 0) .and. (V(i3,j3) > c(i3,j3))) V(i3,j3)=c(i3,j3)
    if ((V(i3,j3) < 0) .and. (abs(V(i3,j3)) > c(i3,j3))) V(i3,j3)=-c(i3,j3)
    Q(i3,j3)=V(i3,j3)*A(i3,j3)

    if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
      Q(i1,j1)=Q(i1-1,j1)
      A(i1,j1)=A(i1-1,j1)
      V(i1,j1)=V(i1-1,j1)
      if (V(i1-1,j1) > c(i1-1,j1)) then
        y(i1,j1)=y(i1-1,j1)
      else
        y(i1,j1)=Hc1
      end if
      c(i1,j1)=FindCel(i1,j1)
      if (vacuum(i1-1,j1)) then
        A(i1,j1)=FindArea(i1,j1)
        V(i1,j1)=Q(i1,j1)/A(i1,j1)
      end if
    else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
      Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
      V(i1,j1)=c(i1-1,j1)
      A(i1,j1)=Q(i1,j1)/V(i1,j1)
      c(i1,j1)=c(i1-1,j1)
      y(i1,j1)=FindDepth(i1,j1)
    end if

    dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

    dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

    dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

    dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

  case(12)
  
    ThreeWayCases(k,12)=ThreeWayCases(k,12)+1

    if  (Q(i1-1,j1) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1)))**1.5_wp) then
      Q(i1,j1)=Q(i1-1,j1)
      A(i1,j1)=A(i1-1,j1)
      V(i1,j1)=V(i1-1,j1)
      if (V(i1-1,j1) > c(i1-1,j1)) then
        y(i1,j1)=y(i1-1,j1)
      else
        y(i1,j1)=Hc1
      end if
      if (vacuum(i1-1,j1)) then
        A(i1,j1)=FindArea(i1,j1)
        V(i1,j1)=Q(i1,j1)/A(i1,j1)
      end if
      c(i1,j1)=FindCel(i1,j1)
    else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
      Q(i1,j1)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i1-1,j1)*(max(0.0_wp,y(i1-1,j1))**1.5_wp)
      V(i1,j1)=c(i1-1,j1)
      A(i1,j1)=Q(i1,j1)/V(i1,j1)
      c(i1,j1)=c(i1-1,j1)
      y(i1,j1)=FindDepth(i1,j1)
    end if

    if (IsDrowned3) then
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dV31=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV31/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV32/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV33) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
      
      y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    else
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dV31=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV31/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV32/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV33) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
      
  	  y(i3,j3)=max(Hc3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    end if
    
    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    V(i3,j3)=V(i3,j3) + (dV31+2*dV32+2*dV33+dV34)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

!	y(i3,j3)=max(JuncHead(k)-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*D(j3))
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

  case(13)
   
   ThreeWayCases(k,13)=ThreeWayCases(k,13)+1
    
	y(i1,j1)=Junc(k)%Head-abs(Z1)+0.5*Reach(j1)%dX*Reach(j1)%So !assume that the previous value of the junction head will be the final depth for reach 1
! PRK 6/14/2022 trapping a strange case that has come up in ALCOSURGE 
    if (y(i1,j1) < 0.0_wp) y(i1,j1)=y(i1-1,j1)
! PRK 6/14/2022 end
    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=findcel(i1,j1)
    V(i1,j1)=Kr1 - y(i1,j1)*g/cr1
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
    if ((V(i1,j1) > 0) .and. (V(i1,j1) > c(i1,j1))) V(i1,j1)=c(i1,j1)
    if ((V(i1,j1) < 0) .and. (abs(V(i1,j1)) > c(i1,j1))) V(i1,j1)=-c(i1,j1)
    Q(i1,j1)=V(i1,j1)*A(i1,j1)

    if  (Q(i3-1,j3) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)) then
      Q(i3,j3)=Q(i3-1,j3)
      A(i3,j3)=A(i3-1,j3)
      V(i3,j3)=V(i3-1,j3)
      if (V(i3-1,j3) > c(i3-1,j3)) then
        y(i3,j3)=y(i3-1,j3)
      else
        y(i3,j3)=Hc3
      end if
      if (vacuum(i3-1,j3)) then
        A(i3,j3)=FindArea(i3,j3)
        V(i3,j3)=Q(i3,j3)/A(i3,j3)
      end if
      c(i3,j3)=FindCel(i3,j3)
    else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
      Q(i3,j3)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)
      V(i3,j3)=c(i3-1,j3)
      A(i3,j3)=Q(i3,j3)/V(i3,j3)
      c(i3,j3)=c(i3-1,j3)
      y(i3,j3)=FindDepth(i3,j3)
    end if

    dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

    dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

    dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

    dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

  case(14)
   
   ThreeWayCases(k,14)=ThreeWayCases(k,14)+1
    
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
      if (V(i1-1,j1) > c(i1-1,j1)) then
        Q(i1,j1)=Q(i1-1,j1)
        V(i1,j1)=Q(i1,j1)/A(i1,j1)
      else if (abs(V(i1,j1)) > c(i1,j1)) then
        if ((y(i1-1,j1) < alpha1*y(i1,j1)) .or. (y(i1-2,j1) < alpha2*y(i1,j1))) then
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

      y(i3,j3)=max(initHfrac*Reach(j3)%D,Junc(k)%Head-Z3+Reach(j3)%dX/2*Reach(j3)%So)
      c(i3,j3)=FindCel(i3,j3)
      A(i3,j3)=FindArea(i3,j3)
      V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
      Q(i3,j3)=V(i3,j3)*A(i3,j3)
      if (V(i3-1,j3) > c(i3-1,j3)) then
        Q(i3,j3)=Q(i3-1,j3)
        V(i3,j3)=Q(i3,j3)/A(i3,j3)
      else if (abs(V(i3,j3)) > c(i3,j3)) then
        if ((y(i3-1,j3) < alpha1*y(i3,j3)) .or. (y(i3-2,j3) < alpha2*y(i3,j3))) then
          if (V(i3,j3) > 0) then
            V(i3,j3)=c(i3,j3)
            Q(i3,j3)=A(i3,j3)*V(i3,j3)
          else
            V(i3,j3)=-c(i3,j3)
            Q(i3,j3)=A(i3,j3)*V(i3,j3)
          end if
        else
          A(i3,j3)=A(i3-1,j3)
          y(i3,j3)=y(i3-1,j3)
          V(i3,j3)=V(i3-1,j3)
          Q(i3,j3)=Q(i3-1,j3)
        end if
      else
        Q(i3,j3)=A(i3,j3)*V(i3,j3)
      end if

      if ((Junc(k)%Inflow-Junc(k)%Outflow) < (Q(i2,j2) - Q(i1,j1)- Q(i3,j3)+Junc(k)%CurrentArea*(Junc(k)%Head-Junc(k)%HeadOld)/dT)) then
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

    dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

    dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

    dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

    dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

  case(15)
   
    ThreeWayCases(k,15)=ThreeWayCases(k,15)+1

	y(i1,j1)=Junc(k)%Head-abs(Z1)+0.5*Reach(j1)%dX*Reach(j1)%So !assume that the previous value of the junction head will be the final depth for reach 1
    A(i1,j1)=FindArea(i1,j1)
    V(i1,j1)=Kr1 - y(i1,j1)*g/cr1
    c(i1,j1)=FindCel(i1,j1)
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
    if ((V(i1,j1) > 0) .and. (V(i1,j1) > c(i1,j1))) V(i1,j1)=c(i1,j1)
    if ((V(i1,j1) < 0) .and. (abs(V(i1,j1)) > c(i1,j1))) V(i1,j1)=-c(i1,j1)
    Q(i1,j1)=V(i1,j1)*A(i1,j1)

    if (IsDrowned3) then
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dV31=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV31/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV32/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV33) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
      
      y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    else
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dV31=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dV32=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV31/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dV33=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV32/2) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dV34=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*(V(i3,j3)+dV33) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
      
  	  y(i3,j3)=max(Hc3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)

    end if
    
    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    V(i3,j3)=V(i3,j3) + (dV31+2*dV32+2*dV33+dV34)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

!	y(i3,j3)=max(JuncHead(k)-Z3+Kdown(j3)*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*D(j3))
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)
    
  case(16)
   
   ThreeWayCases(k,16)=ThreeWayCases(k,16)+1

    if  (Q(i3-1,j3) > 0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)) then
      Q(i3,j3)=Q(i3-1,j3)
      A(i3,j3)=A(i3-1,j3)
      V(i3,j3)=V(i3-1,j3)
      if (V(i3-1,j3) > c(i3-1,j3)) then
        y(i3,j3)=y(i3-1,j3)
      else
        y(i3,j3)=Hc3
      end if
      if (vacuum(i3-1,j3)) then
        A(i3,j3)=FindArea(i3,j3)
        V(i3,j3)=Q(i3,j3)/A(i3,j3)
      end if
      c(i3,j3)=FindCel(i3,j3)
    else  !if flow from upstream spur tunnel is small, force behaviour like a weir discharge
      Q(i3,j3)=0.6_wp*(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.5_wp*Tfs(i3-1,j3)*(max(0.0_wp,y(i3-1,j3))**1.5_wp)
      V(i3,j3)=c(i3-1,j3)
      A(i3,j3)=Q(i3,j3)/V(i3,j3)
      c(i3,j3)=c(i3-1,j3)
      y(i3,j3)=FindDepth(i3,j3)
    end if

    if (IsDrowned1) then
      dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

 	  y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    
    else
        
      dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

 	  y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    
    end if
    
    V(i1,j1)=V(i1,j1) + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

  case(17)
   
   ThreeWayCases(k,17)=ThreeWayCases(k,17)+1
   
    y(i3,j3)=Junc(k)%Head-abs(Z3)+0.5*Reach(j3)%dX*Reach(j3)%So !assume that the previous value of the junction head will be the final depth for reach 1
    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Kr3 - y(i3,j3)*g/cr3
    c(i3,j3)=FindCel(i3,j3)
    if (V(i3-1,j3) > c(i3-1,j3)) then
      Q(i3,j3)=Q(i3-1,j3)
      V(i3,j3)=Q(i3,j3)/A(i3,j3)
    else if (abs(V(i3,j3)) > c(i3,j3)) then
      if ((y(i3-1,j3) < alpha1*y(i3,j3)) .or. (y(i3-2,j3) < alpha2*y(i3,j3))) then
        if (V(i3,j3) > 0) then
          V(i3,j3)=c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        else
          V(i3,j3)=-c(i3,j3)
          Q(i3,j3)=A(i3,j3)*V(i3,j3)
        end if
      else
        A(i3,j3)=A(i3-1,j3)
        y(i3,j3)=y(i3-1,j3)
        V(i3,j3)=V(i3-1,j3)
        Q(i3,j3)=Q(i3-1,j3)
      end if
    else
      Q(i3,j3)=A(i3,j3)*V(i3,j3)
    end if
    if ((V(i3,j3) > 0) .and. (V(i3,j3) > c(i3,j3))) V(i3,j3)=c(i3,j3)
    if ((V(i3,j3) < 0) .and. (abs(V(i3,j3)) > c(i3,j3))) V(i3,j3)=-c(i3,j3)
    Q(i3,j3)=V(i3,j3)*A(i3,j3)

	if (IsDrowned1) then
      dV11=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV12=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV13=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV14=dT*g/Reach(j1)%dX*(-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
	  
      y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    
    else
      dV11=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
      dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
      dY1= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

      dV12=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
      dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
      dY2= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

      dV13=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
      dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
      dY3= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

      dV14=dT*g/Reach(j1)%dX*(-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown*(V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
      dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
      dY4= dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow) + Q(i3,j3) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
      
      y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)

   end if        
    V(i1,j1)=V(i1,j1) + (dV11+2*dV12+2*dV13+dV14)/6.0_wp
    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

  case(18)
   
   ThreeWayCases(k,18)=ThreeWayCases(k,18)+1

    if (IsDrowned1) then
      dV11= dT*g/Reach(j1)%dX * (-(Junc(k)%Head-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
    else
      dV11= dT*g/Reach(j1)%dX * (-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
    end if        
    dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    if (IsDrowned3) then
      dV31= dT*g/Reach(j3)%dX * (-(Junc(k)%Head-Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
    else
      dV31=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
    end if
    dY1=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )

    if (IsDrowned1) then
      dV12= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY1/2-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
    else
      dV12= dT*g/Reach(j1)%dX * (-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
    end if
    dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    if (IsDrowned3) then
      dV32= dT*g/Reach(j3)%dX * (-(Junc(k)%Head+dY1/2 -Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
    else
      dV32=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
    end if
    dY2=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3)+dV31/2 ) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )

    if (IsDrowned1) then
      dV13= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY2/2-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
    else
      dV13= dT*g/Reach(j1)%dX * (-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
    end if
    dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    if (IsDrowned3) then
      dV33= dT*g/Reach(j3)%dX * (-(Junc(k)%Head+dY2/2 -Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
    else
      dV33=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
    end if
    dY3=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3)+dV32/2 ) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )

    if (IsDrowned1) then
      dV14= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY3-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
    else
      dV14= dT*g/Reach(j1)%dX * (-Hc1 + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
    end if
    dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    if (IsDrowned3) then
      dV34= dT*g/Reach(j3)%dX * (-(Junc(k)%Head+dY3 -Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
    else
      dV34=dT*g/Reach(j3)%dX*(-Hc3 + y(i3-1,j3) - Reach(j3)%Kdown*(V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
    end if
    dY4=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3)+dV33 ) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )

    V(i1,j1)=V(i1,j1) + (dV11+2*(dV12+dV13)+dV14)/6.0_wp
    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    V(i3,j3)=V(i3,j3) + (dV31+2*(dV32+dV33)+dV34)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp

	if (IsDrowned1) then
      y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    else
      y(i1,j1)=max(Hc1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)
    end if
    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

    if (IsDrowned3) then
	  y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
    else
      y(i3,j3)=max(Hc3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)
    end if    
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

  case(19)
   
   ThreeWayCases(k,19)=ThreeWayCases(k,19)+1
    
	Junc(k)%CurrentArea=0.0_wp

! we follow here the calculation sequence outlined in Wylie and Streeter
    B1=amax/(g*Reach(j1)%Apipe)
    B2=amax/(g*Reach(j2)%Apipe)
    B3=amax/(g*Reach(j3)%Apipe)

! conversion to f (Darcy) from n (Manning)

    f1=12.69920844327177*g*Reach(j1)%n**2/(Reach(j1)%D)**0.3333
    f2=12.69920844327177*g*Reach(j2)%n**2/(Reach(j2)%D)**0.3333
    f3=12.69920844327177*g*Reach(j3)%n**2/(Reach(j3)%D)**0.3333

    Rf1=f1*Reach(j1)%dX/(2*g*Reach(j1)%D*Reach(j1)%Apipe**2)
    Rf2=f2*Reach(j2)%dX/(2*g*Reach(j2)%D*Reach(j2)%Apipe**2)
    Rf3=f3*Reach(j3)%dX/(2*g*Reach(j3)%D*Reach(j3)%Apipe**2)

! assuming there are two C+ characteristic lines arriving into the
! junction, and one C- characteristic line

    CP1=y(i1-1,j1) + z(i1-1,j1) + B1*Q(i1-1,j1)
    CM2=y(i2+1,j2) + z(i2+1,j2) - B2*Q(i2+1,j2)
    CP3=y(i3-1,j3) + z(i3-1,j3) + B3*Q(i3-1,j3)

    BP1=B1 + Rf1*abs(Q(i1-1,j1))
    BM2=B2 + Rf2*abs(Q(i2+1,j2))
    BP3=B3 + Rf3*abs(Q(i3-1,j3))

    SC=CP1/BP1 + CP3/BP3 + CM2/BM2
    SB=1/BP1 + 1/BP3 + 1/BM2

! overall pressure is uniform for all cells and junction, the value
! calculated by the expression. I am assuming that Qn=0
  
    Htotal=SC/SB + (1/SB)*0.0_dp
    Junc(k)%Head=Htotal - Junc(k)%Elev
    Q(i1,j1)= (CP1 - Htotal)/BP1
    Q(i2,j2)= (Htotal - CM2)/BM2
    Q(i3,j3)= (CP3 - Htotal)/BP3

!	Hdown=0.
!	Hup=Z3
!	converged=.false.
!	niter=0
!	do while(not(converged))
!	  Hcorr=0.5*(Hdown+Hup)
!	  niter=niter+1

!      Q(i1,j1)= CP1/BP1 - (Junc(k)%Head+Hcorr)/BP1
!      Q(i2,j2)=-CM2/BM2 + (Junc(k)%Head+Hcorr)/BM2
!      Q(i3,j3)= CP3/BP3 - (Junc(k)%Head-Z3+Hcorr)/BP3

!	  if (Q(i2,j2) > Q(i1,j1)+Q(i3,j3)) then
!	    Hup=Hcorr
!	  else
!	    Hdown=Hcorr
!	  end if
!     if (Hup < Hdown) then
!        Hcorr=Hup
!        Hup=Hdown
!        Hdown=Hcorr
!      end if
!      if (((abs((Hup-Hdown)) < 0.0001) .and. (niter > 5)) .or. (niter > 100)) converged=.true.

!    end do

!	Junc(k)%Head=Junc(k)%Head+Hcorr

    y(i1,j1)=Junc(k)%Head
    y(i2,j2)=Junc(k)%Head
    y(i3,j3)=Junc(k)%Head-Z3

    A(i1,j1)=FindArea(i1,j1)
    V(i1,j1)=Q(i1,j1)/A(i1,j1)

    A(i2,j2)=FindArea(i2,j2)
    V(i2,j2)=Q(i2,j2)/A(i2,j2)

    A(i3,j3)=FindArea(i3,j3)
    V(i3,j3)=Q(i3,j3)/A(i3,j3)

! PRK 11/2/2015 important update regarding TPA calcs for blind junctions, so that negative depths are handled correctly
    i=i1
    j=j1
    if (A(i,j) >= Reach(j)%Apipe) then
   	  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
    else if (((Aold(i,j) >= Reach(j)%Apipe) .and. ((Aold(i-1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i-1,j)))) .or. &
             ((vacuum(i,j)) .and. ((Aold(i-1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i-1,j))))) then
      hs(i,j)=y(i,j)-Reach(j)%D
      A(i,j)=Reach(j)%Apipe*(1.0_wp + g*hs(i,j)/amax**2)
    else
      hs(i,j)=0.0_wp
	end if
    hc(i,j)=FindCentroid(i,j)
    Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
    Tfs(i,j)=FindTfs(i,j)
    Rh(i,j)=FindRh(i,j)
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

    i=i2
    j=j2

    if (A(i,j) >= Reach(j)%Apipe) then
   	  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
    else if (((Aold(i,j) >= Reach(j)%Apipe) .and. ((Aold(i+1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i+1,j)))) .or. &
             ((vacuum(i,j)) .and. ((Aold(i+1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i+1,j))))) then
      hs(i,j)=y(i,j)-Reach(j)%D
      A(i,j)=Reach(j)%Apipe*(1.0_wp + g*hs(i,j)/amax**2)
    else
      hs(i,j)=0.0_wp
	end if
    hc(i,j)=FindCentroid(i,j)
    Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
    Tfs(i,j)=FindTfs(i,j)
    Rh(i,j)=FindRh(i,j)
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

    i=i3
    j=j3
    if (A(i,j) >= Reach(j)%Apipe) then
   	  hs(i,j)=amax**2*(A(i,j)-Reach(j)%Apipe)/(g*Reach(j)%Apipe)
    else if (((Aold(i,j) >= Reach(j)%Apipe) .and. ((Aold(i-1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i-1,j)))) .or. &
             ((vacuum(i,j)) .and. ((Aold(i-1,j) > Reach(j)%Apipe*0.9999) .or. (vacuum(i-1,j))))) then
      hs(i,j)=y(i,j)-Reach(j)%D
      A(i,j)=Reach(j)%Apipe*(1.0_wp + g*hs(i,j)/amax**2)
    else
      hs(i,j)=0.0_wp
	end if
    hc(i,j)=FindCentroid(i,j)
    Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
    Tfs(i,j)=FindTfs(i,j)
    Rh(i,j)=FindRh(i,j)
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


  case(20)
   
    ThreeWayCases(k,20)=ThreeWayCases(k,20)+1
!if this is the first time that this routine is used, we need to initialize the trapped air chamber
    if (Junc(k)%InitializeHch) then
      Junc(k)%Hch=atm
      Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
      call FindCurrentChamberVolume
      Junc(k)%InitializeHch=.false.
    end if

    dV11= dT*g/Reach(j1)%dX * (-(Junc(k)%Head + Junc(k)%Hch-atm -Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1) )*abs( V(i1,j1) )/(2*g) )
    dV21= dT*g/Reach(j2)%dX * ( Junc(k)%Head  + Junc(k)%Hch-atm - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2) )*abs( V(i2,j2) )/(2*g) )
    dV31= dT*g/Reach(j3)%dX * (-(Junc(k)%Head + Junc(k)%Hch-atm -Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3) )*abs( V(i3,j3) )/(2*g) )
    dY1=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3) ) + A(i1,j1)*( V(i1,j1) ) - A(i2,j2)*( V(i2,j2) ) )
    Pdrop=(Junc(k)%Hch-atm)/Junc(k)%Hch
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    dHch1=-npoly*(Junc(k)%Hch)*(-dY1*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber)

    dV12= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY1/2 + Junc(k)%Hch-atm+dHch1/2-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV11/2 )*abs( V(i1,j1)+dV11/2 )/(2*g) )
    dV22= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY1/2 + Junc(k)%Hch-atm+dHch1/2)- y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV21/2 )*abs( V(i2,j2)+dV21/2 )/(2*g) )
    dV32= dT*g/Reach(j3)%dX * (-(Junc(k)%Head+dY1/2 + Junc(k)%Hch-atm+dHch1/2-Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3)+dV31/2 )*abs( V(i3,j3)+dV31/2 )/(2*g) )
    dY2=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3)+dV31/2 ) + A(i1,j1)*( V(i1,j1)+dV11/2 ) - A(i2,j2)*( V(i2,j2)+dV21/2 ) )
    Pdrop=(Junc(k)%Hch+dHch1/2-atm)/(Junc(k)%Hch+dHch1/2)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)/2
    dHch2=-npoly*(Junc(k)%Hch+dHch1/2)*(-(dY1/2)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY1/2*Junc(k)%CurrentArea)

    dV13= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY2/2 + Junc(k)%Hch-atm+dHch2/2-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV12/2 )*abs( V(i1,j1)+dV12/2 )/(2*g) )
    dV23= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY2/2 + Junc(k)%Hch-atm+dHch2/2) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV22/2 )*abs( V(i2,j2)+dV22/2 )/(2*g) )
    dV33= dT*g/Reach(j3)%dX * (-(Junc(k)%Head+dY2/2 + Junc(k)%Hch-atm+dHch2/2-Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3)+dV32/2 )*abs( V(i3,j3)+dV32/2 )/(2*g) )
    dY3=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3)+dV32/2 ) + A(i1,j1)*( V(i1,j1)+dV12/2 ) - A(i2,j2)*( V(i2,j2)+dV22/2 ) )
    Pdrop=(Junc(k)%Hch+dHch2/2-atm)/(Junc(k)%Hch+dHch2/2)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)/2
    dHch3=-npoly*(Junc(k)%Hch+dHch2/2)*(-(dY2/2)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY2/2*Junc(k)%CurrentArea)

    dV14= dT*g/Reach(j1)%dX * (-(Junc(k)%Head+dY3 + Junc(k)%Hch-atm+dHch3-Z1) + y(i1-1,j1) - Reach(j1)%Kdown* ( V(i1,j1)+dV13 )*abs( V(i1,j1)+dV13 )/(2*g) )
    dV24= dT*g/Reach(j2)%dX * ( (Junc(k)%Head+dY3 + Junc(k)%Hch-atm+dHch3) - y(i2+1,j2) - Reach(j2)%Kup* ( V(i2,j2)+dV23 )*abs( V(i2,j2)+dV23 )/(2*g) )
    dV34= dT*g/Reach(j3)%dX * (-(Junc(k)%Head+dY3 + Junc(k)%Hch-atm+dHch3-Z3) + y(i3-1,j3) - Reach(j3)%Kdown* ( V(i3,j3)+dV33 )*abs( V(i3,j3)+dV33 )/(2*g) )
    dY4=  dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow) + A(i3,j3)*( V(i3,j3)+dV33 ) + A(i1,j1)*( V(i1,j1)+dV13 ) - A(i2,j2)*( V(i2,j2)+dV23 ) )
    Pdrop=(Junc(k)%Hch+dHch3-atm)/(Junc(k)%Hch+dHch3)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    dHch4=-npoly*(Junc(k)%Hch+dHch3)*(-(dY3)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY3*Junc(k)%CurrentArea)

    V(i1,j1)=V(i1,j1) + (dV11+2*(dV12+dV13)+dV14)/6.0_wp
    V(i2,j2)=V(i2,j2) + (dV21+2*(dV22+dV23)+dV24)/6.0_wp
    V(i3,j3)=V(i3,j3) + (dV31+2*(dV32+dV33)+dV34)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1+2*(dY2+dY3)+dY4)/6.0_wp
    Junc(k)%Hch=Junc(k)%Hch + (dHch1+2*(dHch2+dHch3)+dHch4)/6.0_wp
    Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
    call FindCurrentChamberVolume
    Pdrop=(Junc(k)%Hch-atm)/(Junc(k)%Hch)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    Junc(k)%Qair=VolExhaust/dT

	y(i1,j1)=max(Junc(k)%Head-Z1+Reach(j1)%Kdown*V(i1,j1)*abs(V(i1,j1))/(2*g),initHfrac*Reach(j1)%D)+Junc(k)%Hch-atm
    A(i1,j1)=FindArea(i1,j1)
    c(i1,j1)=FindCel(i1,j1)
    Q(i1,j1)=A(i1,j1)*V(i1,j1)

    y(i2,j2)=max(Junc(k)%Head-Z2-Reach(j2)%Kup*V(i2,j2)*abs(V(i2,j2))/(2*g),initHfrac*Reach(j2)%D)+Junc(k)%Hch-atm
    A(i2,j2)=FindArea(i2,j2)
    c(i2,j2)=FindCel(i2,j2)
    Q(i2,j2)=A(i2,j2)*V(i2,j2)

	y(i3,j3)=max(Junc(k)%Head-Z3+Reach(j3)%Kdown*V(i3,j3)*abs(V(i3,j3))/(2*g),initHfrac*Reach(j3)%D)+Junc(k)%Hch-atm
    A(i3,j3)=FindArea(i3,j3)
    c(i3,j3)=FindCel(i3,j3)
    Q(i3,j3)=A(i3,j3)*V(i3,j3)

  case default
    
    write(*,*) 'undefined case in three-way dropshaft routine for junction',k
    write (538,*) k,JunctionCase,IsWeir1,IsWeir3,IsSurch1,IsSurch3,IsSurchMain
    stop
    
end select

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,i5,33es16.8,i5,es16.8)') T,JunctionCase,Junc(k)%Head,Junc(k)%Inflow,Junc(k)%Outflow,y(i1,j1),y(i2,j2),y(i3,j3),y(i1-1,j1),y(i2+1,j2),y(i3-1,j3), &
	  V(i1,j1),V(i2,j2),V(i3,j3),V(i1-1,j1),V(i2+1,j2),V(i3-1,j3),c(i1,j1),c(i2,j2),c(i3,j3),c(i1-1,j1),c(i2+1,j2),c(i3-1,j3),Q(i1,j1),Q(i2,j2),Q(i3,j3), &
	  Q(i1-1,j1),Q(i2+1,j2),Q(i3-1,j3),hs(i1,j1),hs(i2,j2),hs(i3,j3),hs(i1-1,j1),hs(i2+1,j2),hs(i3-1,j3),niter,Hconj1
      write (1234,*) T,Hc1,Hc3,IsDrowned1,IsDrowned3
	end if
  end do
end if

if (JunctionCase /= 19) then
  call TPACalcBC(i1,j1,JunctionCase)
  call TPACalcBC(i2,j2,JunctionCase)
  call TPACalcBC(i3,j3,JunctionCase)
end if

99 continue

return

end subroutine