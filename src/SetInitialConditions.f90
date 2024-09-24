subroutine SetInitialConditions
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Initializes the depthj in reach cells and junctions, or reads a hot start file	
!%=====================================================================
        
use GlobalVariables
use GlobalFunctions

implicit none

! check for user-specified initial HGL elevations
if (StartCondition == 1) then
  open (unit=77,file='InitialHGL.inp',status='old',err=903)
  read (77,*)
  read (77,*)
  do j=1,Nreaches
      read (77,*) Reach(j)%InitHGL
  end do
  read (77,*)
  read (77,*)
  do k=1,NJuncs
      read (77,*) Junc(k)%InitHGL
  end do
end if

!> set reach values
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    if (i == 1) then
	  z(i,j)=Reach(j)%UpElev
	else
      z(i,j)=z(i-1,j)-Reach(j)%So*Reach(j)%dX
	end if
	y(i,j)=Reach(j)%D*initHfrac
	V(i,j)=0.0_wp
	hs(i,j)=0.0_wp
    hc(i,j)=FindCentroid(i,j)
    A(i,j)=FindArea(i,j)
    Rh(i,j)=FindRh(i,j)
    Tfs(i,j)=FindTfs(i,j)
    Q(i,j)=A(i,j)*V(i,j)
    Imom(i,j)=hc(i,j)*A(i,j)
	c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
  end do
end do

if (IsHotStart) then
  open (unit=444,file=trim(HotStartID)//'.hotstart.inp',status='old',err=60)
  do j=1,Nreaches
    do i=1,Reach(j)%Ncells
	  read (444,*) z(i,j),y(i,j),V(i,j),A(i,j),Q(i,j),MaxVacuum(i,j),MaxHGL(i,j)
      if (y(i,j) >= Reach(j)%D) then
	    select case(Reach(j)%XSec)
	      case(1)
            hs(i,j)=y(i,j)-Reach(j)%D
	        hc(i,j)=0.5_wp*Reach(j)%D
	        Imom(i,j)=A(i,j)*(hc(i,j)+hs(i,j))
	        Tfs(i,j)=g*A(i,j)/amax**2
	        Rh(i,j)=0.25_wp*Reach(j)%D
	        c(i,j)=amax
		  case(2)
            hs(i,j)=y(i,j)-Reach(j)%D
	        hc(i,j)=0.5_wp*Reach(j)%D
	        Imom(i,j)=A(i,j)*(hc(i,j)+hs(i,j))
	        Tfs(i,j)=g*A(i,j)/amax**2
	        Rh(i,j)=Reach(j)%Width*Reach(j)%D/(2.0_wp*(Reach(j)%Width+Reach(j)%D))
	        c(i,j)=amax
		end select
      else
	    select case(Reach(j)%XSec)
	      case(1)
            hs(i,j)=0.0_wp
            hc(i,j)=FindCentroid(i,j)
            Imom(i,j)=hc(i,j)*A(i,j)
            Tfs(i,j)=20_wp*sqrt(y(i,j)*(Reach(j)%D-y(i,j)))
            tta=2.0_wp*acos(1.0_wp-2.0_wp*y(i,j)/Reach(j)%D)
            Rh(i,j)=(1-sin(tta)/tta)*0.25_wp*Reach(j)%D
            c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
		  case(2)
		    hs(i,j)=0.0_wp
			hc(i,j)=0.5_wp*y(i,j)
			Imom(i,j)=hc(i,j)*A(i,j)
			Tfs(i,j)=Reach(j)%Width
			Rh(i,j)=Reach(j)%Width+2.0_wp*y(i,j)
			c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
		end select
      end if
	end do
  end do
  do k=1,Njuncs
    read (444,*) Junc(k)%Head
  end do
  MaxHGL=0.0_wp
  MaxVacuum=0.0_wp
!  if (DoAirPockets) then
!    do while(not(eof(444)))
!      read (444,*) j,iFRT1(j),iFRT2(j),LPocket(j),LColumn(j),Ha(j),QColumn(j),VolPocket(j)
!      HasPocket(j)=.true.
!      InitializePocket(j)=.false.
!    end do
!  end if
  close (444)
  go to 47
end if

!> set initial HGL if desired
if (IsHGLInit) then
  do j=1,Nreaches
    do i=1,Reach(j)%Ncells
	  y(i,j)=max(Reach(j)%D*initHfrac,Reach(j)%InitHGL-z(i,j))
	  V(i,j)=0.0_wp
      A(i,j)=FindArea(i,j)
      hc(i,j)=FindCentroid(i,j)
      Rh(i,j)=FindRh(i,j)
      Tfs(i,j)=FindTfs(i,j)
      Q(i,j)=A(i,j)*V(i,j)
      Imom(i,j)=hc(i,j)*A(i,j)
	  c(i,j)=min(amax,sqrt(g*A(i,j)/Tfs(i,j)))
      Q(i,j)=V(i,j)*A(i,j)
	end do
  end do
end if

!> set initial junction heads

do k=1,NJuncs
  if (IsHGLInit) then
    if (Junc(k)%ReachCount == 0) then
      Junc(k)%Head=initHfrac
    else
      Junc(k)%Head=max(Junc(k)%InitHGL-Junc(k)%Elev,y(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1)))
    end if
  else if (Junc(k)%ReachCount == 0) then
    Junc(k)%Head=initHfrac
  else
    Junc(k)%Head=y(Junc(k)%ReachCell(1),Junc(k)%ReachNo(1))
  end if
!  if (Junc(k)%Option == -8) call FindCurrentChamberVolume
end do


47 continue
   
!> calculate initial volume
Reach(:)%Volume=0.0_wp
ReachVolTot=0.0_wp
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    Reach(j)%Volume=Reach(j)%Volume+Reach(j)%dX*A(i,j)
  end do
  ReachVolTot=ReachVolTot+Reach(j)%Volume
end do

Junc(:)%Volume=0.0_wp
JuncVolTot=0.0_wp
do k=1,NJuncs
  call FindCurrentShaftArea
  Junc(k)%Volume=Junc(k)%CurrentArea*Junc(k)%Head
  JuncVolTot=JuncVolTot+Junc(k)%Volume
end do

InitialSystemVol=ReachVolTot+JuncVolTot

do k=1,NJuncs
  write (52,1609) T,k,ReachVolTot,Junc(k)%Volume,Junc(k)%Inflow,Junc(k)%Outflow,Junc(k)%Closure
end do
1609 format (es13.6,',',i3,7(',',es16.8))

go to 61

903 write(*,*) 'InitHGL.inp not found'
IOError=.true.
go to 61
60 write (*,*) 'restart file missing'
IOError=.true.

61 continue

return

end subroutine