subroutine SingleDropshaftDown
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Boundary condition for downstream end of a single reach	
!%=====================================================================

use GlobalVariables
use GlobalFunctions

implicit none

integer :: niter,kk
real(kind=8) :: dV1,dV2,dV3,dV4,dY1,dY2,dY3,dY4,Yup,Ydown,Vr,cr,Sfr,yr,Kr,r,Zcorr,AL,CL,VL,HcL,TfsL,Pdrop,VolExhaust,dHch1,dHch2,dHch3,dHch4
logical :: converged,IsWeirFlow,IsSurch,CaseProblem

converged=.false.
niter=0

if (dT == 0.0_wp) then
  go to 99
end if

j=Junc(k)%ReachNo(1)
i=Junc(k)%ReachCell(1)
Zcorr=Junc(k)%ReachElev(1)-Junc(k)%Elev
r=dT/Reach(j)%dX
Junc(k)%HeadOld=Junc(k)%Head

! PRK 1/27/2015 modify bailout check so that calcs always performed for surcharged junction
if (Junc(k)%Head < 0.999*Reach(j)%D) then
  if (IsHGLInit) then
    if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0.0_dp) .and. (abs(Q(i-1,j)) < 1.0e-4_wp*sqrt(g*Reach(j)%D**5))) then
      go to 99
    end if
  else
    if (((Junc(k)%Inflow-Junc(k)%Outflow) == 0.0_dp) .and. (abs(Q(i-1,j)) < 1.0e-4_wp*sqrt(g*Reach(j)%D**5))) then
      go to 99
    end if
  end if
end if
! PRK 1/27/2015 end

Qlim=4.0D-4*sqrt(g*Reach(j)%D**5)

HcL=FindCriticalDepth(i-1,j)
if (Zcorr > 0.0_dp) then
  if (Junc(k)%HeadOld < Zcorr+HcL) then
    IsWeirFlow=.true.
  else
    IsWeirFlow=.false.
  end if
else 
  IsWeirFlow=.false.
end if

if (y(i-1,j) < Reach(j)%D) then
  IsSurch=.false.
else
  IsSurch=.true.
end if

if ((IsWeirFlow) .and. (not(IsSurch))) JunctionCase=1
if ((not(IsWeirFlow)) .and. (not(IsSurch))) JunctionCase=2
if (IsSurch) JunctionCase=3
if ((IsSurch) .and. (Junc(k)%Option == -8)) JunctionCase=4

if ((not(Junc(k)%InitializeHch)) .and. (JunctionCase /= 4)) then
  Junc(k)%Hch=atm
  Junc(k)%Lch=0.0_wp
  Junc(k)%VolChamber=0.0_wp
  Junc(k)%InitializeHch=.true.
end if

select case(JunctionCase)

  case(1)

    if (Q(i-1,j) > 0.6*(2.0_wp/3.)*sqrt(2.0_wp*g)*(0.5*Tfs(i,j))*(max(0.0_wp,y(i,j)))**1.5) then
      Q(i,j)=Q(i-1,j)
      V(i,j)=V(i-1,j)
      if (V(i-1,j) > c(i-1,j)) then
        y(i,j)=y(i-1,j)
      else
        y(i,j)=HcL
      end if
      A(i,j)=FindArea(i,j)
      c(i,j)=FindCel(i,j)
    else  !if flow from upstream tunnel is small, force behaviour like a weir discharge
      Q(i,j)=0.6*(2.0_wp/3.)*sqrt(2.0_wp*g)*(0.5*Tfs(i,j))*(max(0.0_wp,y(i,j)))**1.5
      A(i,j)=A(i-1,j)
      V(i,j)=Q(i,j)/A(i,j)
      c(i,j)=FindCel(i,j)
      y(i,j)=FindDepth(i,j)
    end if

    Junc(k)%Head=Junc(k)%HeadOld+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(i,j))/Junc(k)%CurrentArea

  case(2)

    Vr=(V(i,j) + r*(-V(i,j)*c(i-1,j)+c(i,j)*V(i-1,j))) / (1.0_wp +r*(V(i,j)-V(i-1,j)+c(i,j)-c(i-1,j)))
    cr=(c(i,j) + r*Vr*(c(i-1,j) - c(i,j))) / (1.0_wp + r*(c(i,j) - c(i-1,j)))
    Sfr=Reach(j)%n**2*Vr*abs(Vr)/Rh(i,j)**1.333
    yr=y(i-1,j) + r*(Vr+cr)*(y(i,j) - y(i-1,j))
    Kr=Vr+g*yr/cr-g*(Sfr-Reach(j)%So)*dT

    Yup=1.00*max(y(i,j),Reach(j)%D)
    Ydown=initHfrac*Reach(j)%D
!    y(i,j)=0.001
    niter=0
    converged=.false.

    do while (not(converged))
      if ((V(i-1,j) > c(i-1,j)) .and. (y(i-1,j)+V(i-1,j)**2/(2*g) > Junc(k)%HeadOld-Zcorr)) then
        y(i,j)=y(i-1,j)
        V(i,j)=V(i-1,j)
        c(i,j)=c(i-1,j)
        A(i,j)=A(i-1,j)
        Q(i,j)=Q(i-1,j)
        Junc(k)%Head=Junc(k)%HeadOld + dT*((Junc(k)%Inflow-Junc(k)%Outflow) + Q(i,j))/Junc(k)%CurrentArea
        converged=.true.
      else
	    y(i,j)=0.5*(Yup+Ydown)
!        y(i,j)=max(initHfrac*Reach(j)%D,JuncHead(k)-Zcorr-0.5*dX(j)*So(j))
        niter=niter+1
        A(i,j)=FindArea(i,j)
        c(i,j)=FindCel(i,j)
        V(i,j)=Kr - g*y(i,j)/cr ! applicable if flow is subcritical
        if (abs(V(i,j)) > c(i,j)) then
          if ((y(i-1,j) < 1.15*y(i,j)) .or. (y(i-2,j) < 1.25*y(i,j))) then
            if (V(i,j) > 0.0_dp) then
              V(i,j)=c(i,j)
            else
              V(i,j)=-c(i,j)
            end if
	      end if
        end if
        Q(i,j)=A(i,j)*V(i,j)

        Junc(k)%Head=Junc(k)%HeadOld+dT*((Junc(k)%Inflow-Junc(k)%Outflow)+Q(i,j))/Junc(k)%CurrentArea
        if (Junc(k)%Head > 0.0_wp) then
          if (y(i,j)+0.5_wp*Reach(j)%dX*Reach(j)%So > Junc(k)%Head-Zcorr) then
            Yup=Junc(k)%Head-0.5_wp*Reach(j)%dX*Reach(j)%So-Zcorr
          else
            Ydown=Junc(k)%Head-0.5_wp*Reach(j)%dX*Reach(j)%So-Zcorr
		  end if
        else
          Yup=0.5_wp*Yup ! if junctionhead is negative, decrease the Yup abruptly
        end if
        if (Yup < Ydown) then
          y(i,j)=Yup
          Yup=Ydown
          Ydown=y(i,j)
        end if
!        if ((abs((Yup-Ydown)/Ydown) < epsilon) .and. (niter > 5)) converged=.true.
        if ((abs(Yup-Ydown) < epsilon) .and. (niter > 5)) converged=.true.
        if ((abs((Yup-Ydown)/Ydown) < 0.01) .and. (niter > 99)) converged=.true.
	    if (niter > 10000) then
	      divergence=.true.
	      go to 99
	    end if
      end if
    end do

    if (y(i,j) < initHfrac*Reach(j)%D) y(i,j)=initHfrac*Reach(j)%D

! last correction for bore invading junction
    if ((y(i,j)+y(i-1,j) > 100*initHfrac*Reach(j)%D) .and. (y(i,j) < 0.7*y(i-1,j))) then
      Q(i,j)=Q(i-1,j)
      A(i,j)=A(i-1,j)
      V(i,j)=V(i-1,j)
      y(i,j)=y(i-1,j)
      Junc(k)%Head=Junc(k)%HeadOld + dT*((Junc(k)%Inflow-Junc(k)%Outflow) + Q(i,j))/Junc(k)%CurrentArea
    end if

  case(3)

    dV1=dT*g/Reach(j)%dX * (-(Junc(k)%Head-Zcorr) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)) * abs(V(i,j)) )
    dY1=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)) )

    dV2=dT*g/Reach(j)%dX * (-(Junc(k)%Head+dY1/2.0_wp-Zcorr) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)+dV1/2.0_wp) * abs(V(i,j)+dV1/2.0_wp) )
    dY2=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)+dV1/2.0_wp) )

    dV3=dT*g/Reach(j)%dX * (-(Junc(k)%Head+dY2/2.0_wp-Zcorr) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)+dV2/2.0_wp) * abs(V(i,j)+dV2/2.0_wp) )
    dY3=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)+dV2/2.0_wp) )

    dV4=dT*g/Reach(j)%dX * (-(Junc(k)%Head+dY3-Zcorr) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)+dV3) * abs(V(i,j)+dV3) )
    dY4=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)+dV3) )

    V(i,j)=V(i,j) + (dV1 + 2.0_wp*(dV2 + dV3)+dV4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1 + 2.0_wp*(dY2 + dY3) +dY4)/6.0_wp
    y(i,j)=Junc(k)%Head - Zcorr + Reach(j)%Kdown*V(i,j)*abs(V(i,j))/(2.0_wp*g)
    A(i,j)=FindArea(i,j)
    Q(i,j)=A(i,j)*V(i,j)
    c(i,j)=FindCel(i,j)

  case(4)

    if (Junc(k)%InitializeHch) then
      Junc(k)%Hch=atm
      Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
      call FindCurrentChamberVolume
      Junc(k)%InitializeHch=.false.
    end if

    dV1=dT*g/Reach(j)%dX * (-(Junc(k)%Head-Zcorr+ Junc(k)%Hch-atm) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)) * abs(V(i,j)) )
    dY1=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)) )
    Pdrop=max((Junc(k)%Hch-atm)/(Junc(k)%Hch),0.0_wp)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2.0_wp*RT*abs(Pdrop)),Pdrop)
    dHch1=-npoly*(Junc(k)%Hch)*(-dY1*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber)

    dV2=dT*g/Reach(j)%dX * (-(Junc(k)%Head+dY1/2.0_wp-Zcorr+ Junc(k)%Hch-atm+dHch1/2.0_wp) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)+dV1/2.0_wp) * abs(V(i,j)+dV1/2.0_wp) )
    dY2=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)+dV1/2.0_wp) )
    Pdrop=max((Junc(k)%Hch+dHch1/2.0_wp-atm)/(Junc(k)%Hch+dHch1/2.0_wp),0.0_wp)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2.0_wp*RT*abs(Pdrop)),Pdrop)/2.0_wp
    dHch2=-npoly*(Junc(k)%Hch+dHch1/2.0_wp)*(-(dY1/2.0_wp)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY1/2.0_wp*Junc(k)%CurrentArea)

    dV3=dT*g/Reach(j)%dX * (-(Junc(k)%Head+dY2/2.0_wp-Zcorr+ Junc(k)%Hch-atm+dHch2/2.0_wp) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)+dV2/2.0_wp) * abs(V(i,j)+dV2/2.0_wp) )
    dY3=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)+dV2/2.0_wp) )
    Pdrop=max((Junc(k)%Hch+dHch2/2.0_wp-atm)/(Junc(k)%Hch+dHch2/2.0_wp),0.0_wp)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2.0_wp*RT*abs(Pdrop)),Pdrop)/2.0_wp
    dHch3=-npoly*(Junc(k)%Hch+dHch2/2.0_wp)*(-(dY2/2.0_wp)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY2/2.0_wp*Junc(k)%CurrentArea)

    dV4=dT*g/Reach(j)%dX * (-(Junc(k)%Head+dY3-Zcorr+ Junc(k)%Hch-atm+dHch3) + y(i-1,j) - Reach(j)%Kdown/(2.0_wp*g) * (V(i,j)+dV3) * abs(V(i,j)+dV3) )
    dY4=dT*(1/Junc(k)%CurrentArea) * ( (Junc(k)%Inflow-Junc(k)%Outflow)  + A(i,j)*(V(i,j)+dV3) )
    Pdrop=max((Junc(k)%Hch+dHch3-atm)/(Junc(k)%Hch+dHch3),0.0_wp)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2.0_wp*RT*abs(Pdrop)),Pdrop)
    dHch4=-npoly*(Junc(k)%Hch+dHch3)*(-(dY3)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY3*Junc(k)%CurrentArea)

    V(i,j)=V(i,j) + (dV1 + 2.0_wp*(dV2 + dV3)+dV4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1 + 2.0_wp*(dY2 + dY3) +dY4)/6.0_wp
    Junc(k)%Hch=max(Junc(k)%Hch + (dHch1+2.0_wp*(dHch2+dHch3)+dHch4)/6.0_wp,atm)
    Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
    call FindCurrentChamberVolume
    Pdrop=max((Junc(k)%Hch-atm)/(Junc(k)%Hch),0.0_wp)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2.0_wp*RT*abs(Pdrop)),Pdrop)
    Junc(k)%Qair=VolExhaust/dT

    y(i,j)=Junc(k)%Head - Zcorr + Reach(j)%Kdown*V(i,j)*abs(V(i,j))/(2.0_wp*g)
    A(i,j)=FindArea(i,j)
    Q(i,j)=A(i,j)*V(i,j)

  case default

     write (*,*) 'indeterminate case in single downstream dropshaft routine for junction',k
     write (538,*) k,JunctionCase,IsSurch,IsWeirFlow
     stop
     
end select

if ((y(i,j) >= Reach(j)%D) .and. (A(i,j) < Reach(j)%Apipe)) then
  write (812,*) T,JunctionCase
end if
call TPACalcBC(i,j,JunctionCase)

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,i4,11es16.8,i8)') T,JunctionCase,Junc(k)%Head,Junc(k)%Inflow,Junc(k)%Outflow,y(i-1,j),y(i,j),Q(i-1,j),Q(i,j),V(i-1,j),V(i,j),c(i-1,j),c(i,j),niter
	end if
  end do
end if


99 continue
return
end subroutine
