subroutine UpdateAQPocket
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Updates air pocket region - expect some changes to this routine at some point	
!%=====================================================================

use GlobalVariables
use GlobalFunctions

implicit none

real(wp), save :: Qin,Qout,Qvent,APocket,AColumn,VColumn
integer(i4), save :: iStart,iEnd
!$omp threadprivate (Qin,Qout,Qvent,APocket,AColumn,VColumn,iStart,iEnd)

iStart=Reach(j)%iFRT1-PocketBuffer
iEnd=Reach(j)%iFRT2+PocketBuffer
APocket=Reach(j)%VolPocket/Reach(j)%LPocket
AColumn=(Reach(j)%LPocket*(Reach(j)%Apipe-APocket)+2*PocketBuffer*Reach(j)%dX*Reach(j)%Apipe)/Reach(j)%LColumn
VColumn=Reach(j)%QColumn/AColumn
Rey=VColumn*Reach(j)%D/Nua+1.0_wp
if (Rey > 2000.0_wp) then
  f=0.25_wp/(dlog10(kss/(3.7_wp*Reach(j)%D)+5.74_wp/Rey**0.9_wp))**2
else
  f=64.0_wp/Rey
end if

Qin=Q(iStart-1,j)
Qout=Q(iEnd+1,j)
Reach(j)%QColDiff=Qin-Qout

!dQ1=dT*(g*AColumn/LColumn(j))*(y(iStart-1,j)-y(iEnd,j)+(1-f*LColumn(j)/D(j))*QColumn(j)*abs(QColumn(j))/(2*g*AColumn*AColumn))
dQ1=dT*(g*AColumn/Reach(j)%LColumn)*(y(iStart-1,j)-y(iEnd+1,j)+(iEnd-iStart+2)*Reach(j)%So*Reach(j)%dX-f*Reach(j)%LColumn/Reach(j)%D*Reach(j)%QColumn*abs(Reach(j)%QColumn)/(2.0_wp*g*AColumn*AColumn))
dHa1=dT*(npoly*Reach(j)%Ha/Reach(j)%VolPocket)*(Qin-Qout-Qvent)
dVa1=-dT*(Qin-Qout)

!dQ2=dT*(g*AColumn/LColumn(j))*(y(iStart-1,j)-y(iEnd,j)+(1-f*LColumn(j)/D(j))*(QColumn(j)+dQ1/2.0_wp)*abs(QColumn(j)+dQ1/2.0_wp)/(2*g*AColumn*AColumn))
dQ2=dT*(g*AColumn/Reach(j)%LColumn)*(y(iStart-1,j)-y(iEnd+1,j)+(iEnd-iStart+2)*Reach(j)%So*Reach(j)%dX-f*Reach(j)%LColumn/Reach(j)%D*(Reach(j)%QColumn+dQ1/2.0_wp)*abs(Reach(j)%QColumn+dQ1/2.0_wp)/(2.0_wp*g*AColumn*AColumn))
dHa2=dT*(npoly*(Reach(j)%Ha+dHa1/2.0_wp)/(Reach(j)%VolPocket+dVa1/2.0_wp))*(Qin-Qout-Qvent)
dVa2=-dT*(Qin-Qout)

!dQ3=dT*(g*AColumn/LColumn(j))*(y(iStart-1,j)-y(iEnd,j)+(1-f*LColumn(j)/D(j))*(QColumn(j)+dQ2/2.0_wp)*abs(QColumn(j)+dQ2/2.0_wp)/(2*g*AColumn*AColumn))
dQ3=dT*(g*AColumn/Reach(j)%LColumn)*(y(iStart-1,j)-y(iEnd+1,j)+(iEnd-iStart+2)*Reach(j)%So*Reach(j)%dX-f*Reach(j)%LColumn/Reach(j)%D*(Reach(j)%QColumn+dQ2/2.0_wp)*abs(Reach(j)%QColumn+dQ2/2.0_wp)/(2*g*AColumn*AColumn))
dHa3=dT*(npoly*(Reach(j)%Ha+dHa2/2.0_wp)/(Reach(j)%VolPocket+dVa2/2.0_wp))*(Qin-Qout-Qvent)
dVa3=-dT*(Qin-Qout)

!dQ4=dT*(g*AColumn/LColumn(j))*(y(iStart-1,j)-y(iEnd,j)+(1-f*LColumn(j)/D(j))*(QColumn(j)+dQ3)*abs(QColumn(j)+dQ3)/(2*g*AColumn*AColumn))
dQ4=dT*(g*AColumn/Reach(j)%LColumn)*(y(iStart-1,j)-y(iEnd+1,j)+(iEnd-iStart+2)*Reach(j)%So*Reach(j)%dX  -f*Reach(j)%LColumn/Reach(j)%D*(Reach(j)%QColumn+dQ3)*abs(Reach(j)%QColumn+dQ3)/(2*g*AColumn*AColumn))
dHa4=dT*(npoly*(Reach(j)%Ha+dHa3)/(Reach(j)%VolPocket+dVa3))*(Qin-Qout-Qvent)
dVa4=-dT*(Qin-Qout)

Reach(j)%QColumn=Reach(j)%QColumn+(dQ1+2*(dQ2+dQ3)+dQ4)/6.
Reach(j)%Ha=Reach(j)%Ha+(dHa1+2*(dHa2+dHa3)+dHa4)/6.
Reach(j)%VolPocket=Reach(j)%VolPocket+(dVa1+2*(dVa2+dVa3)+dVa4)/6.

! check if pocket has gotten too small to continue tracking (because of venting)
if ((Reach(j)%Ha*Reach(j)%VolPocket)/(atm*Reach(j)%VolPocketInit) < 1.0D-03) then
    Reach(j)%HasPocket=.false.
    Reach(j)%PocketReleaseT=T
    Reach(j)%VolPocket=0.0_wp
    dYp=y(Reach(j)%iFRT1-PocketBuffer-1,j)-y(Reach(j)%iFRT2+Pocketbuffer+1,j)
    do i=Reach(j)%iFRT1-PocketBuffer,Reach(j)%iFRT2+Pocketbuffer
        y(i,j)=y(Reach(j)%iFRT1-PocketBuffer-1,j)-dYp*float(i)/float(Reach(j)%iFRT2-Reach(j)%iFRT1+1)
        Q(i,j)=Reach(j)%QColumn
        A(i,j)=FindArea(i,j)
    end do
    go to 800
end if
! now do TPA-like calculations for these cells

APocket=Reach(j)%VolPocket/Reach(j)%LPocket
select case(Reach(j)%XSec)
  case(1)
    tta=FindTheta(Reach(j)%Apipe-APocket,Reach(j)%D)
    Reach(j)%YUnder=0.5*Reach(j)%D*(1-cos(0.5*tta))
  case(2)
    Reach(j)%YUnder=(Reach(j)%Apipe-APocket)/Reach(j)%Width
  case(3)
    call GetLookupValue(Reach(j)%Apipe-APocket,Reach(j)%YUnder,Reach(j)%ReachShapeMap,2,1)
end select

do i=iStart,Reach(j)%iFRT1-1
  Q(i,j)=Reach(j)%QColumn
!  y(i,j)=D(j)+Ha(j)-atm
!  hs(i,j)=Ha(j)-atm
!  y(i,j)=D(j)+Ha(j)-atm+(i-iStart)*dX(j)*So(j)-f*((i-iStart)*dX(j)/D(j))*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
!  hs(i,j)=Ha(j)-atm+i*dX(j)*So(j)-f*((i-iStart)*dX(j)/D(j))*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
  if (i == iStart) then
    y(i,j)=Reach(j)%D/2+Reach(j)%Ha-atm+(Reach(j)%D-Reach(j)%YUnder)+(i-iStart)*Reach(j)%dX*Reach(j)%So-f*(i-iStart)*Reach(j)%dX/Reach(j)%D*(Reach(j)%QColumn+dQ2/2)*abs(Reach(j)%QColumn+dQ2/2)/(2*g*Acolumn*Acolumn)
!    y(i,j)=Ha(j)-atm+YUnder(j)+dX(j)*So(j)-f*dX(j)/D(j)*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
  else
    y(i,j)=y(i-1,j)-f*Reach(j)%dX/Reach(j)%D*Reach(j)%QColumn*abs(Reach(j)%QColumn)/(2*g*Acolumn*Acolumn)+Reach(j)%dX*Reach(j)%So
  end if
  hs(i,j)=y(i,j)-Reach(j)%D
!  hc(i,j)=D(j)/2.
  hc(i,j)=FindCentroid(i,j)
  A(i,j)=Reach(j)%Apipe*(1.0_wp +hs(i,j)*g/amax**2)
  V(i,j)=Q(i,j)/A(i,j)
  Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
  Tfs(i,j)=g*A(i,j)/amax**2
!  Rh(i,j)=0.25*D(j)
  Rh(i,j)=FindRh(i,j)
  c(i,j)=amax
  if (isnan(hc(i,j))) write (6262,*) T,i,j,y(i,j),hs(i,j),A(i,j),'start buffer'
end do

do i=Reach(j)%iFRT1,Reach(j)%iFRT2
  Q(i,j)=Reach(j)%QColumn
!  y(i,j)=D(j)+Ha(j)-atm
!  hs(i,j)=Ha(j)
!  y(i,j)=D(j)+Ha(j)-atm + (i-iStart)*dX(j)*So(j)-f*((i-iStart)*dX(j)/D(j))*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
!  hs(i,j)=Ha(j)+(i-iStart)*dX(j)*So(j)-f*((i-iStart)*dX(j)/D(j))*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
  y(i,j)=y(i-1,j)-f*Reach(j)%dX/Reach(j)%D*Reach(j)%QColumn*abs(Reach(j)%QColumn)/(2*g*Acolumn*Acolumn)+Reach(j)%dX*Reach(j)%So
  hs(i,j)=y(i,j)-Reach(j)%D
  APocket=Reach(j)%VolPocket/Reach(j)%LPocket
  select case(Reach(j)%XSec)
    case(1)
      tta=FindTheta(Reach(j)%Apipe-APocket,Reach(j)%D)
      Reach(j)%YUnder=0.5*Reach(j)%D*(1-cos(0.5*tta))
      hc(i,j)=Reach(j)%YUnder-Reach(j)%D*(0.5_wp-(2.0_wp/3.0_wp)*((sin(tta/2.0_wp))**3)/(tta-sin(tta)))
      A(i,j)=0.125*Reach(j)%D*Reach(j)%D*(tta-sin(tta))
      V(i,j)=Q(i,j)/A(i,j)
      Imom(i,j)=A(i,j)*(hc(i,j)+hs(i,j))
      Tfs(i,j)=2.*sqrt(Reach(j)%YUnder*(Reach(j)%D-Reach(j)%YUnder))
      Rh(i,j)=0.25*Reach(j)%D*(1-sin(tta)/tta)
      c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
    case(2)
      Reach(j)%YUnder=(Reach(j)%Apipe-APocket)/Reach(j)%Width
      hc(i,j)=0.5*Reach(j)%YUnder
      A(i,j)=Reach(j)%Width*Reach(j)%YUnder
      V(i,j)=Q(i,j)/A(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)
      Tfs(i,j)=Reach(j)%Width
	  Rh(i,j)=A(i,j)/(Reach(j)%Width+2.*y(i,j))
      c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
    case(3)
      call GetLookupValue(Reach(j)%Apipe-APocket,Reach(j)%YUnder,Reach(j)%ReachShapeMap,2,1)
      call GetLookupValue(Reach(j)%YUnder,hc(i,j),Reach(j)%ReachShapeMap,1,5)
      call GetLookupValue(Reach(j)%YUnder,A(i,j),Reach(j)%ReachShapeMap,1,2)
      V(i,j)=Q(i,j)/A(i,j)
      Imom(i,j)=A(i,j)*hc(i,j)
      call GetLookupValue(Reach(j)%YUnder,Tfs(i,j),Reach(j)%ReachShapeMap,1,4)
      call GetLookupValue(Reach(j)%YUnder,Rh(i,j),Reach(j)%ReachShapeMap,1,3)
      c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
  end select
  if (isnan(hc(i,j))) then
      write (6262,*) T,i,j,y(i,j),hs(i,j),A(i,j),'inside pocket'
  end if
end do

do i=Reach(j)%iFRT2+1,iEnd
  Q(i,j)=Reach(j)%QColumn
!  y(i,j)=D(j)+Ha(j)-atm
!  hs(i,j)=Ha(j)-atm
!  y(i,j)=D(j)+Ha(j)-atm + (i-iStart)*dX(j)*So(j)-f*((i-iStart)*dX(j)/D(j))*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
!  hs(i,j)=Ha(j)-atm+(i-iStart)*dX(j)*So(j)-f*((i-iStart)*dX(j)/D(j))*(Qcolumn(j)+dQ2/2)*abs(Qcolumn(j)+dQ2/2)/(2*g*Acolumn*Acolumn)
  y(i,j)=y(i-1,j)-f*Reach(j)%dX/Reach(j)%D*Reach(j)%QColumn*abs(Reach(j)%QColumn)/(2*g*Acolumn*Acolumn)+Reach(j)%dX*Reach(j)%So
  hs(i,j)=y(i,j)-Reach(j)%D
  hc(i,j)=FindCentroid(i,j)
  A(i,j)=Reach(j)%Apipe*(1.0_wp +hs(i,j)*g/amax**2)
  V(i,j)=Q(i,j)/A(i,j)
  Imom(i,j)=A(i,j)*hc(i,j)+A(i,j)*hs(i,j)
  Tfs(i,j)=g*A(i,j)/amax**2
  Rh(i,j)=FindRh(i,j)
  c(i,j)=amax
  if (isnan(hc(i,j))) write (6262,*) T,i,j,y(i,j),hs(i,j),A(i,j),'end buffer'
end do

800 continue
    
return

end subroutine UpdateAQPocket