module GlobalFunctions
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Utility functions called by core calculations (from Open MP parallel region)	
!%=====================================================================

use GlobalVariables

implicit none

contains
  real(wp) recursive function FindArea(idum,jdum)
  implicit none
  integer :: idum,jdum
  real(wp), save :: xval,yval
!$omp threadprivate(xval,yval)

  if (y(idum,jdum) < 0.0_wp) then
    FindArea=A(idum,jdum)
  else if (y(idum,jdum) > Reach(jdum)%D) then
    hs(idum,jdum)=y(idum,jdum)-Reach(jdum)%D
    FindArea=Reach(jdum)%Apipe+hs(idum,jdum)*g*Reach(jdum)%Apipe/amax**2
  else
    hs(idum,jdum)=0.0_wp
    select case(Reach(jdum)%XSec)
	  case(1)
        tta=2.0_wp*acos(1.0_wp-2.0_wp*y(idum,jdum)/Reach(jdum)%D)
        FindArea=0.125_wp*Reach(jdum)%D**2*(tta-sin(tta))
	  case(2)
	    FindArea=y(idum,jdum)*Reach(jdum)%Width
      case(3)
        xval=y(idum,jdum)
        yval=0.0_wp
        call GetLookupValue(xval,yval,Reach(jdum)%ReachShapeMap,1,2)
        FindArea=yval
	end select
  end if

  end function

  real(wp) recursive function FindCel(idum,jdum)
  implicit none
  integer :: idum,jdum
  real(wp), save :: Adum2,Tfsdum,ydum
!$omp threadprivate(Adum2,Tfsdum,ydum)
  ydum=y(idum,jdum)

  if (y(idum,jdum) == 0.0_dp) then
    FindCel=c(idum,jdum)
  else if (ydum < Reach(jdum)%D) then
    Tfsdum=FindTfs(idum,jdum)
    Adum2=FindArea(idum,jdum)
    FindCel=min(amax,sqrt(g*Adum2/Tfsdum))
  else
!    Tfsdum=g*Adum/amax**2
    FindCel=amax
  end if

  end function

  real(wp) recursive function FindTfs(idum,jdum)
  implicit none
  integer :: idum,jdum
  real(wp), save :: ydum,xval,yval
!$omp threadprivate(ydum,xval,yval)
  ydum=y(idum,jdum)

  if (hs(idum,jdum) /= 0.0_wp) then
    FindTfs=g*A(idum,jdum)/amax**2
  else
    select case(Reach(jdum)%XSec)
	  case(1)
        FindTfs=2.0_wp*sqrt(ydum*(Reach(jdum)%D-ydum))
	  case(2)
	    FindTfs=Reach(jdum)%Width
      case(3)
        xval=y(idum,jdum)
        yval=0.0_wp
        call GetLookupValue(xval,yval,Reach(jdum)%ReachShapeMap,1,4)
        FindTfs=yval
	end select
  end if
  
  end function

  real(wp) recursive function FindCentroid(idum,jdum)
  implicit none
  integer :: idum,jdum
  real(wp), save :: ydum,tta,xval,yval
!$omp threadprivate(ydum,tta,xval,yval)
  ydum=y(idum,jdum)

  if (hs(idum,jdum) /= 0.0_wp) then
    select case(Reach(jdum)%XSec)
	  case(1)
        FindCentroid=Reach(jdum)%D/2.0_wp
	  case(2)
	    FindCentroid=Reach(jdum)%D/2.0_wp
      case(3)
        FindCentroid=CustomShapeTables(TableRowCount(Reach(jdum)%ReachShapeMap),5,Reach(jdum)%ReachShapeMap)
      end select
  else if (ydum == 0.0_wp) then
      FindCentroid=0.0_wp
  else
    select case(Reach(jdum)%XSec)
	  case(1)
        tta=2.0_wp*acos(1.0_wp-2.0_wp*ydum/Reach(jdum)%D)
        FindCentroid=ydum-Reach(jdum)%D*(0.5_wp-(2.0_wp/3.0_wp)*((sin(tta/2.0_wp))**3)/(tta-sin(tta)))
	  case(2)
	    FindCentroid=ydum/2.0_wp
      case(3)
        xval=y(idum,jdum)
        yval=0.0_wp
        call GetLookupValue(xval,yval,Reach(jdum)%ReachShapeMap,1,5)
        FindCentroid=yval
	end select
  end if

  if (isnan(FindCentroid)) then
        FindCentroid=Reach(jdum)%D/2.0_wp
  end if
  
  end function

  real recursive function findTheta(Adum,Ddum)
  implicit none
  real(wp) :: Adum,Ddum
  real(wp), save :: Aunit,ThetaDum,Func,FuncPrime,emax,correction
!$omp threadprivate(Aunit,ThetaDum,Func,FuncPrime,emax,correction)
  Aunit=Adum/(0.25_wp*pi*Ddum**2)
  emax=0.00000001_wp
  correction=1.0_wp
  ThetaDum=pi

  if (Aunit <= 0.0_wp) then
    findTheta=0.00001_wp
  else if (Aunit <= 0.9999_wp) then
    do while (abs(correction) > emax)
      Func=Ddum*Ddum/8.0_wp*(ThetaDum-sin(ThetaDum))-Adum
      FuncPrime=Ddum*Ddum/8.0_wp*(1.0_wp-cos(ThetaDum))
      correction=Func/FuncPrime
	  ThetaDum=ThetaDum-correction
    end do
    findTheta=ThetaDum
  else
    findTheta=2.0_wp*pi+10000.0_wp*(Aunit-1.0_wp)*(2.0_wp*pi-6.12748610682723_wp)  ! used in version 3.03
  end if

  end function

  real(wp) recursive function FindRh(idum,jdum)
  implicit none
  integer :: idum,jdum
  real(wp), save :: ttadum,xval,yval,Ddum
!$omp threadprivate(ttadum,xval,yval)
  
  dDum=Reach(jdum)%D
  select case(Reach(jdum)%XSec)
    case(1)
      if ((y(idum,jdum) < Ddum) .and. (vacuum(idum,jdum) == .false.)) then
        ttadum=findTheta(A(idum,jdum),Reach(jdum)%D)
        FindRh=0.25_wp*Reach(jdum)%D*(1-sin(ttadum)/ttadum)
      else
        FindRh=0.25_wp*Reach(jdum)%D
      end if
    case(2)
      if ((y(idum,jdum) < Ddum) .and. (vacuum(idum,jdum) == .false.)) then
        FindRh=A(idum,jdum)/(Reach(jdum)%Width+2.0_wp*y(idum,jdum))
      else
        FindRh=Reach(jdum)%Width*Reach(jdum)%D/(2.0_wp*(Reach(jdum)%Width+Reach(jdum)%D))
      end if
    case(3)
      if ((y(idum,jdum) < Reach(jdum)%D) .and. (vacuum(idum,jdum) == .false.)) then
        xval=y(idum,jdum)
        yval=0.0_wp
        call GetLookupValue(xval,yval,Reach(jdum)%ReachShapeMap,1,3)
        FindRh=yval
      else
        FindRh=CustomShapeTables(TableRowCount(Reach(jdum)%ReachShapeMap),3,Reach(jdum)%ReachShapeMap)
      end if
  end select
  end function
  
  real(wp) recursive function FindDepth(idum,jdum)
  implicit none
  integer :: idum,jdum
  real(wp), save :: ttadum,xval,yval
!$omp threadprivate(ttadum,xval,yval)
  
  select case(Reach(jdum)%XSec)
    case(1)
      ttadum=findTheta(A(idum,jdum),Reach(jdum)%D)
      FindDepth=0.5_wp*Reach(jdum)%D*(1.0_wp-cos(0.5_wp*ttadum))
    case(2)
      FindDepth=A(idum,jdum)/Reach(jdum)%Width
    case(3)
      xval=A(idum,jdum)
      yval=0.0_wp
      call GetLookupValue(xval,yval,Reach(jdum)%ReachShapeMap,2,1)
      FindDepth=yval
  end select
  end function

  real(wp) recursive function FindCriticalDepth(idum,jdum)
  implicit none
  integer :: idum,jdum,nit
  real(wp), save :: Yupper,Ylower,Ytest,Aupper,Tupper,Fupper,Atest,Ttest,Ftest,LittleQ,xval,yval,Alower,Tlower,Flower,Qmin
  logical :: converge
  
    select case(Reach(jdum)%XSec)
    case(1)
      Yupper=0.999_wp*Reach(jdum)%D
      Ylower=0.001_wp
      converge=.false.
      nit=0
      Qmin=4.0D-4*sqrt(g*Reach(jdum)%D**5)
      tta=2.0_wp*acos(1.0_wp-2.0_wp*Ylower/Reach(jdum)%D)
      Alower=0.125_wp*Reach(jdum)%D**2*(tta-sin(tta))
      Tlower=2.0_wp*sqrt(Ylower*(Reach(jdum)%D-Ylower))
      tta=2.0_wp*acos(1.0_wp-2.0_wp*Yupper/Reach(jdum)%D)
      Aupper=0.125_wp*Reach(jdum)%D**2*(tta-sin(tta))
      Tupper=2.0_wp*sqrt(Yupper*(Reach(jdum)%D-Yupper))
      Fupper=Q(idum,jdum)/Aupper-sqrt(g*Aupper/Tupper)
      Flower=Q(idum,jdum)/Alower-sqrt(g*Alower/Tlower)
      if ((Q(idum,jdum) > Qmin) .and. (Fupper*Flower > 0.0_wp)) then
        Ytest=1.01_wp*Reach(jdum)%D
        converge=.true.
      end if
      do while(not(converge))
        tta=2.0_wp*acos(1.0_wp-2.0_wp*Yupper/Reach(jdum)%D)
        Aupper=0.125_wp*Reach(jdum)%D**2*(tta-sin(tta))
        Tupper=2.0_wp*sqrt(Yupper*(Reach(jdum)%D-Yupper))
        Ytest=(Yupper+Ylower)/2
        tta=2.0_wp*acos(1.0_wp-2.0_wp*Ytest/Reach(jdum)%D)
        Atest=0.125_wp*Reach(jdum)%D**2*(tta-sin(tta))
        Ttest=2.0_wp*sqrt(Ytest*(Reach(jdum)%D-Ytest))
        Fupper=Q(idum,jdum)/Aupper-sqrt(g*Aupper/Tupper)
        Ftest=Q(idum,jdum)/Atest-sqrt(g*Atest/Ttest)
        if (Fupper*Ftest > 0.0_wp) then
          Yupper=Ytest
        else
          Ylower=Ytest
        end if
        nit=nit+1
        if (((Yupper-Ylower) < epsilon) .or. (nit > 20)) converge =.true.
      end do
      FindCriticalDepth=Ytest
      FindCriticalDepth=min(FindCriticalDepth,Reach(jdum)%D)
    case(2)
      LittleQ=Q(idum,jdum)/Reach(jdum)%Width
      FindCriticalDepth=(LittleQ**2/g)**(1.0_wp/3.0_wp)
    case(3)
      xval=Q(idum,jdum)
      yval=0.0_wp
      call GetLookupValue(xval,yval,Reach(jdum)%ReachShapeMap,6,1)
      FindCriticalDepth=yval
    end select
    
    
  
  end function
  function wtime ( )

!*****************************************************************************80
!
!! WTIME returns a reading of the wall clock time.
!
!  Discussion:
!
!    To get the elapsed wall clock time, call WTIME before and after a given
!    operation, and subtract the first reading from the second.
!
!    This function is meant to suggest the similar routines:
!
!      "omp_get_wtime ( )" in OpenMP,
!      "MPI_Wtime ( )" in MPI,
!      and "tic" and "toc" in MATLAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) WTIME, the wall clock reading, in seconds.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer clock_max
  integer clock_rate
  integer clock_reading
  real ( kind = rk ) wtime

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real ( clock_reading, kind = rk ) &
      / real ( clock_rate, kind = rk )

  return
  end function 
end module GlobalFunctions