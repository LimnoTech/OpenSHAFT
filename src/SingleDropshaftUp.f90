subroutine SingleDropshaftUp

use GlobalVariables
use GlobalFunctions

implicit none

integer :: niter,kk
real(wp) :: dV1,dV2,dV3,dV4,dY1,dY2,dY3,dY4,Yup,Ydown,Vs,cs,Sfs,ys,Ks,r,Pdrop,VolExhaust,dHch1,dHch2,dHch3,dHch4
logical :: converged

converged=.false.

if (dT == 0.0_dp) then
  go to 99
end if

j=Junc(k)%ReachNo(1)
i=Junc(k)%ReachCell(1)
r=dT/Reach(j)%dX
Junc(k)%HeadOld=Junc(k)%Head
niter=0
JunctionCase=0

! PRK 1/27/2015 modify bailout check so that calcs always performed for surcharged junction
if (Junc(k)%Head < 0.999_wp*Reach(j)%D) then
  if (IsHGLInit) then
    if ((Junc(k)%Inflow == 0.0_wp) .and. (abs(Q(i+1,j)) < 1.0e-4_wp*sqrt(g*Reach(j)%D**5))) then
      go to 99
    end if
  else
! PRK 12/6/2013 testing higher flow threshold
    if ((Junc(k)%Inflow == 0.0_wp) .and. (abs(Q(i+1,j)) < 1.0e-4_wp*sqrt(g*Reach(j)%D**5))) then
!  if ((JuncInflow(k) == 0.0_dp) .and. (abs(Q(i+1,j)) < 1e-5*sqrt(g*D(j)**5))) then
! PRK 12/6/2013 end
      go to 99
    end if
  end if
end if
! PRK 1/27/2015 end

Qlim=max(4e-4_wp*sqrt(g*Reach(j)%D**5),abs(Q(i+1,j)))

if (y(i+1,j) < 0.999_wp*Reach(j)%D) then
!if (JuncHead(k) < 0.999_wp*D(j)) then
    JunctionCase=1
else if (Junc(k)%Option == -8) then
    JunctionCase=3
else
    JunctionCase=2
end if

if ((not(Junc(k)%InitializeHch)) .and. (JunctionCase /= 3)) then
  Junc(k)%Hch=atm
  Junc(k)%Lch=0.0_wp
  Junc(k)%VolChamber=0.0_wp
  Junc(k)%Qair=0.0_wp
  Junc(k)%InitializeHch=.true.
end if

select case(JunctionCase)

  case(1)

    Vs=(V(i,j) + r*(c(i,j)*V(i+1,j)-c(i+1,j)*V(i,j)))/(1.0_wp + r*(-V(i,j)+V(i+1,j)+c(i,j)-c(i+1,j)))
    cs=(c(i,j) + r*Vs*(c(i,j) - c(i+1,j)))/(1.0_wp + r*(c(i,j)-c(i+1,j)))
    Sfs=Reach(j)%n**2*Vs*abs(Vs)/Rh(i,j)**1.333
    ys=y(i,j) + r*abs(Vs-cs)*(y(i+1,j)-y(i,j))
    Ks=Vs-g*ys/cs-g*(Sfs-Reach(j)%So)*dT

    Yup=max(y(i,j),Reach(j)%D)
!    Ydown=initHfrac*Reach(j)%D
    Ydown=0.0_wp

    niter=0
    converged=.false.

    do while(not(converged))
      y(i,j)=0.5_wp*(Yup+Ydown)
!      Junc(k)%Head=y(i,j)+Reach(j)%dX/2*Reach(j)%So
!      JuncHead(k)=0.5*(Yup+Ydown)
!      y(i,j)=max(initHfrac*D(j),JuncHead(k)-dX(j)/2*So(j))
      niter=niter+1
      A(i,j)=FindArea(i,j)
      c(i,j)=FindCel(i,j)
      V(i,j)=Ks+g*y(i,j)/cs
      if (abs(V(i,j)) > c(i,j)) then
        if ((y(i+1,j) < 1.05_wp*y(i,j)) .or. (y(i+2,j) < 1.10_wp*y(i,j))) then
          if (V(i,j) > 0.0_wp) then
            V(i,j)=c(i,j)
          else
            V(i,j)=-c(i,j)
          end if
        else
          A(i,j)=A(i+1,j)
          V(i,j)=V(i+1,j)
          Q(i,j)=Q(i+1,j)
	      end if
	    end if
	    Q(i,j)=A(i,j)*V(i,j)


	    Junc(k)%Head=Junc(k)%HeadOld+dT*((Junc(k)%Inflow-Junc(k)%Outflow)-Q(i,j))/Junc(k)%CurrentArea
	    if (Junc(k)%Head > 0.0_wp) then
	      if (y(i,j) - 0.5_wp*Reach(j)%dX*Reach(j)%So > Junc(k)%Head) then
	        Yup=Junc(k)%Head + 0.5_wp*Reach(j)%dX*Reach(j)%So
	      else
	        Ydown=Junc(k)%Head + 0.5_wp*Reach(j)%dX*Reach(j)%So
	      end if
        else
!        write (2330,'(f11.4,i5,15es16.8)') T,k,Junc(k)%Head,Yup,Ydown
	      Yup=0.5_wp*Yup
 	    end if
      if (Yup < Ydown) then
	      y(i,j)=Yup
	      Yup=Ydown
	      Ydown=y(i,j)
	    end if

!      if ((abs((Yup-Ydown)/Ydown) < epsilon) .and. (niter > 5)) converged=.true.
!      if ((abs((Yup-Ydown)/Ydown) < 0.01) .and. (niter  >99)) converged=.true.
!	  if ((niter > 50) .and. (Yup > 0.999*JuncHead(k)+dX(j)/2*So(j)) .and. ( Yup < 1.001*JuncHead(k)+dX(j)/2*So(j))) converged=.true.
      if (((abs((Yup-Ydown)) < Qlim*dT/Junc(k)%CurrentArea*tol) .and. (niter > 5)) .or. (niter > 100)) converged=.true.
!      if ((abs(Yup-Ydown) < epsilon) .and. (niter > 5)) converged=.true.
	    if (niter > 10000) then
	      divergence=.true.
	      go to 99
	    end if
    end do
    if (y(i,j) < initHfrac*Reach(j)%D) y(i,j)=initHfrac*Reach(j)%D

! last correction for bore invading junction
    if ((y(i,j)+y(i+1,j) > 100*initHfrac*Reach(j)%D) .and. (y(i,j) < 0.7*y(i+1,j)) .and. (V(i+1,j) < 0.0_dp)) then
      Q(i,j)=Q(i+1,j)
	    A(i,j)=A(i+1,j)
	    V(i,j)=V(i+1,j)
	    y(i,j)=y(i+1,j)
	    Junc(k)%Head=Junc(k)%HeadOld + dT*((Junc(k)%Inflow-Junc(k)%Outflow) - Q(i,j))/Junc(k)%CurrentArea
    end if
! PRK 1/4/2024 For some zero inflow cases with steep slopes, JuncHead may (harmlessly) end up negative. Here we correct this for purely cosmetic considerations
    Junc(k)%Head=max(Junc(k)%Head,0.0_wp)

  case(2)   !this is for pressurized flow
    
    dV1=dT*g/Reach(j)%dX*(Junc(k)%Head - y(i+1,j) - Reach(j)%Kup/(2*g)*V(i,j)*abs(V(i,j)))
    dY1=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*V(i,j))

    dV2=dT*g/Reach(j)%dX*((Junc(k)%Head+dY1/2.0_wp) - y(i+1,j) - Reach(j)%Kup/(2*g)*(V(i,j)+dV1/2.0_wp)*abs((V(i,j)+dV1/2.0_wp)))
    dY2=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*(V(i,j)+dV1/2.0_wp))

    dV3=dT*g/Reach(j)%dX*((Junc(k)%Head+dY2/2.0_wp) - y(i+1,j) - Reach(j)%Kup/(2*g)*(V(i,j)+dV2/2.0_wp)*abs((V(i,j)+dV2/2.0_wp)))
    dY3=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*(V(i,j)+dV2/2.0_wp))

    dV4=dT*g/Reach(j)%dX*((Junc(k)%Head+dY3) - y(i+1,j) - Reach(j)%Kup/(2*g)*(V(i,j)+dV3)*abs((V(i,j)+dV3)))
    dY4=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*(V(i,j)+dV3))

    V(i,j)=V(i,j) + (dV1 + 2*(dV2 + dV3) + dV4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1 + 2*(dY2 + dY3) + dY4)/6.0_wp
    y(i,j)=Junc(k)%Head - Reach(j)%Kup*V(i,j)*abs(V(i,j))/(2*g)
    A(i,j)=FindArea(i,j)
    c(i,j)=FindCel(i,j)
    Q(i,j)=A(i,j)*V(i,j)

  case(3)   !this is for pressurized flow and pressurized shaft air space
    
    if (Junc(k)%InitializeHch) then
      Junc(k)%Hch=atm
      Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
      call FindCurrentChamberVolume
      Junc(k)%InitializeHch=.false.
    end if

    dV1=dT*g/Reach(j)%dX*(Junc(k)%Head + Junc(k)%Hch - atm - y(i+1,j) - Reach(j)%Kup/(2*g)*V(i,j)*abs(V(i,j)))
    dY1=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*V(i,j))
    Pdrop=(Junc(k)%Hch-atm)/(Junc(k)%Hch)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    if (dY1 < 0.0_wp) then
      VolExhaust=max(dY1*Junc(k)%CurrentArea,VolExhaust)
    else
      VolExhaust=min(dY1*Junc(k)%CurrentArea,VolExhaust)
    end if
    dHch1=-npoly*(Junc(k)%Hch)*(-dY1*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber)

    dV2=dT*g/Reach(j)%dX*((Junc(k)%Head+dY1/2.0_wp) + Junc(k)%Hch + dHch1/2 - atm - y(i+1,j) - Reach(j)%Kup/(2*g)*(V(i,j)+dV1/2.0_wp)*abs((V(i,j)+dV1/2.0_wp)))
    dY2=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*(V(i,j)+dV1/2.0_wp))
    Pdrop=(Junc(k)%Hch+dHch1/2-atm)/(Junc(k)%Hch+dHch1/2)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)/2
    if (dY1 < 0.0_wp) then
      VolExhaust=max((dY1/2)*Junc(k)%CurrentArea,VolExhaust)
    else
      VolExhaust=min((dY1/2)*Junc(k)%CurrentArea,VolExhaust)
    end if
    dHch2=-npoly*(Junc(k)%Hch+dHch1/2)*(-(dY1/2)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY1/2*Junc(k)%CurrentArea)

    dV3=dT*g/Reach(j)%dX*((Junc(k)%Head+dY2/2.0_wp) + Junc(k)%Hch + dHch2/2 - atm - y(i+1,j) - Reach(j)%Kup/(2*g)*(V(i,j)+dV2/2.0_wp)*abs((V(i,j)+dV2/2.0_wp)))
    dY3=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*(V(i,j)+dV2/2.0_wp))
    Pdrop=(Junc(k)%Hch+dHch2/2-atm)/(Junc(k)%Hch+dHch2/2)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)/2
    if (dY2 < 0.0_wp) then
      VolExhaust=max((dY2/2)*Junc(k)%CurrentArea,VolExhaust)
    else
      VolExhaust=min((dY2/2)*Junc(k)%CurrentArea,VolExhaust)
    end if
    dHch3=-npoly*(Junc(k)%Hch+dHch2/2)*(-(dY2/2)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY2/2*Junc(k)%CurrentArea)

    dV4=dT*g/Reach(j)%dX*((Junc(k)%Head+dY3) + Junc(k)%Hch + dHch3 - atm - y(i+1,j) - Reach(j)%Kup/(2*g)*(V(i,j)+dV3)*abs((V(i,j)+dV3)))
    dY4=dT*(1/Junc(k)%CurrentArea)*((Junc(k)%Inflow-Junc(k)%Outflow)-A(i,j)*(V(i,j)+dV3))
    Pdrop=(Junc(k)%Hch+dHch3-atm)/(Junc(k)%Hch+dHch3)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    if (dY3 < 0.0_wp) then
      VolExhaust=max(dY3*Junc(k)%CurrentArea,VolExhaust)
    else
      VolExhaust=min(dY3*Junc(k)%CurrentArea,VolExhaust)
    end if
    dHch4=-npoly*(Junc(k)%Hch+dHch3)*(-(dY3)*Junc(k)%CurrentArea+VolExhaust)/(Junc(k)%VolChamber-dY3*Junc(k)%CurrentArea)

    V(i,j)=V(i,j) + (dV1 + 2*(dV2 + dV3) + dV4)/6.0_wp
    Junc(k)%Head=Junc(k)%Head + (dY1 + 2*(dY2 + dY3) + dY4)/6.0_wp
    Junc(k)%Hch=Junc(k)%Hch + (dHch1+2*(dHch2+dHch3)+dHch4)/6.0_wp
    Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head
    call FindCurrentChamberVolume
    Pdrop=(Junc(k)%Hch-atm)/(Junc(k)%Hch)
    VolExhaust=dT*Cd*Junc(k)%Aopen*sign(sqrt(2*RT*abs(Pdrop)),Pdrop)
    if (Junc(k)%Head-Junc(k)%HeadOld < 0.0_wp) then
      VolExhaust=max((Junc(k)%Head-Junc(k)%HeadOld)*Junc(k)%CurrentArea,VolExhaust)
    else
      VolExhaust=min((Junc(k)%Head-Junc(k)%HeadOld)*Junc(k)%CurrentArea,VolExhaust)
    end if
    Junc(k)%Qair=VolExhaust/dT

    y(i,j)=Junc(k)%Head - Reach(j)%Kup*V(i,j)*abs(V(i,j))/(2*g)
    A(i,j)=FindArea(i,j)
    c(i,j)=FindCel(i,j)
    Q(i,j)=A(i,j)*V(i,j)

  case default

     write (*,*) 'indeterminate case in single upstream dropshaft routine for junction',k
     stop
     
end select

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (DiagOutUnits(kk),'(f11.4,i5,15es16.8,i6)') T,JunctionCase,Junc(k)%Head,Junc(k)%Inflow,Junc(k)%Outflow,y(i,j),y(i+1,j),y(i+2,j),Q(i,j), &
	  Q(i+1,j),Q(i+2,j),V(i,j),V(i+1,j),V(i+2,j),c(i,j),c(i+1,j),c(i+2,j),niter
	end if
  end do
end if


call TPACalcBC(i,j,JunctionCase)

99 continue
return
end subroutine SingleDropshaftUp