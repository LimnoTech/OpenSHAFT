subroutine FindOutflow
    
use GlobalVariables

implicit none

real(wp) :: Hoverflow,FlowDirection,OpeningDepth,OpeningTopElev,slope,Qt1,Qt2,Qf,H2
integer(i4) :: index,TWIndex,k1,id,k2,kk

!> set Hoverflow based on rating curve type
select case(Junc(k)%RatingCurveType)
  case(0)
      Hoverflow=max(Junc(k)%Head-(Junc(k)%OverflowElev-Junc(k)%Elev),0.0_wp)
  case(1)
      Hoverflow=max(Junc(k)%Head-(Junc(k)%OverflowElev-Junc(k)%Elev),0.0_wp)
  case(2)
      if (Junc(k)%Head+Junc(k)%Elev < Junc(k)%OverflowElev) then
        Hoverflow=0.0_wp
      else if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then  ! free overflow
        Hoverflow=max(Junc(k)%Head-(Junc(k)%OverflowElev-Junc(k)%Elev),0.0_wp)
      else
        Hoverflow=Junc(k)%Head+Junc(k)%Elev-(Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev)
        FlowDirection=sign(1.0_wp,Hoverflow)
        Hoverflow=abs(Hoverflow)
      end if      
  case(3)
      if (Junc(k)%Head+Junc(k)%Elev < Junc(k)%OverflowElev) then
        Hoverflow=0.0_wp
      else if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then  ! free overflow
        Hoverflow=max(Junc(k)%Head-(Junc(k)%OverflowElev-Junc(k)%Elev),0.0_wp)
      else
        Hoverflow=max(Junc(k)%Head+Junc(k)%Elev-(Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev),0.0_wp)
      end if      
  case(4)
      OpeningTopElev=Junc(k)%OverflowElev+RatingCurveH(NRatingCurveLines(Junc(k)%RatingCurveMap),Junc(k)%RatingCurveMap)
      if (Junc(k)%Head+Junc(k)%Elev < Junc(k)%OverflowElev) then
        Hoverflow=0.0_wp
      else if (Junc(k)%Head+Junc(k)%Elev > OpeningTopElev) then
        OpeningDepth=min(OpeningTopElev-(Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev),OpeningTopElev-Junc(k)%OverflowElev)
        if (OpeningDepth > 0.0_wp) OpeningDepth=OpeningDepth/2.0_wp
        Hoverflow=Junc(k)%Head+Junc(k)%Elev+OpeningDepth
        FlowDirection=sign(1.0_wp,Hoverflow)
        Hoverflow=abs(Hoverflow)
      else if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then  ! free overflow
        Hoverflow=max(Junc(k)%Head-(Junc(k)%OverflowElev-Junc(k)%Elev),0.0_wp)
      else
        Hoverflow=max(Junc(k)%Head+Junc(k)%Elev-(Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev),0.0_wp)
        FlowDirection=sign(1.0_wp,Hoverflow)
        Hoverflow=abs(Hoverflow)
      end if 
  case(5)
      Hoverflow=max(Junc(k)%Head-(Junc(k)%OverflowElev-Junc(k)%Elev),0.0_wp)
end select

!> if Hoverflow is zero, set outflow to zero and skip the rest of the routine
if (Hoverflow == 0.0_wp) then
    Junc(k)%Outflow=0.0_wp
    go to 99
end if

!>find position in rating curve table
index=0
TWIndex=0
if ((Junc(k)%IsRatingCurve) .and. (Junc(k)%RatingCurveType /= 5)) then
  k1=Junc(k)%RatingCurveMap
  do id=1,NRatingCurveLines(k1)
    if (Hoverflow <= RatingCurveH(id,k1)) then
	  index=id
	  exit
	end if
  end do
  if (index == 0) then
    if (Junc(k)%JuncPartner /= 0) then
      write (757,'(es13.6,i3,3es16.8)') T,k,Hoverflow,Junc(k)%Head,Junc(Junc(k)%JuncPartner)%Head
    else
      write (757,'(es13.6,i3,3es16.8)') T,k,Hoverflow,Junc(k)%Head
    end if  
    index=NRatingCurveLines(k1)
  end if
  if (Junc(k)%RatingCurveType == 4) then
    do id=1,NTailWaterH(k1)
      if (Junc(Junc(k)%JuncPartner)%Head <= TailWaterH(id,k1)) then
        TWIndex=id
        exit
      end if
    end do
    if (TWIndex == 0) TWIndex=NTailWaterH(k1)
  end if
end if
  
!>calculate overflow based on rating curve type  
select case(Junc(k)%RatingCurveType)
  case(0)
    Junc(k)%Outflow=(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.611_wp*Junc(k)%WeirLength*Hoverflow**1.5_wp
  case(1)
    slope=(RatingCurveQ(index,k1)-RatingCurveQ(index-1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	Junc(k)%Outflow=RatingCurveQ(index-1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
  case(2)
    if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then
      slope=(RatingCurveQ2(index,k1)-RatingCurveQ2(index-1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=RatingCurveQ2(index-1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
    else
      slope=(RatingCurveQ(index,k1)-RatingCurveQ(index-1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=RatingCurveQ(index-1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
      Junc(k)%Outflow=FlowDirection*Junc(k)%Outflow
    end if
  case(3)
    if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then
      slope=(RatingCurveQ2(index,k1)-RatingCurveQ2(index-1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=RatingCurveQ2(index-1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
    else
      slope=(RatingCurveQ(index,k1)-RatingCurveQ(index-1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=RatingCurveQ(index-1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
    end if
  case(4)
    if (Junc(k)%Head+Junc(k)%Elev > OpeningTopElev) then
      slope=(TWRatingCurveQ(index,TWIndex,k1)-TWRatingCurveQ(index-1,TWIndex,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=TWRatingCurveQ(index-1,TWIndex,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
    else if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then
      slope=(TWRatingCurveQ(index,1,k1)-TWRatingCurveQ(index-1,1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=TWRatingCurveQ(index-1,1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
    else if (Junc(Junc(k)%JuncPartner)%Head > TailWaterH(NTailWaterH(k),k1)) then
      slope=(TWRatingCurveQ(index,TWIndex,k1)-TWRatingCurveQ(index-1,TWIndex,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Junc(k)%Outflow=TWRatingCurveQ(index-1,TWIndex,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
    else
      slope=(TWRatingCurveQ(index,TWIndex-1,k1)-TWRatingCurveQ(index-1,TWIndex-1,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Qt1=TWRatingCurveQ(index-1,TWIndex-1,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
      slope=(TWRatingCurveQ(index,TWIndex,k1)-TWRatingCurveQ(index-1,TWIndex,k1))/(RatingCurveH(index,k1)-RatingCurveH(index-1,k1))
	  Qt2=TWRatingCurveQ(index-1,TWIndex,k1)+slope*(Hoverflow-RatingCurveH(index-1,k1))
      slope=(Qt2-Qt1)/(TailWaterH(TWIndex,k1)-TailWaterH(TWIndex-1,k1))
      Junc(k)%Outflow=Qt1+slope*(Junc(Junc(k)%JuncPartner)%Head-TailWaterH(TWIndex-1,k1))
    end if
  case(5)
    Qf=(2.0_wp/3.0_wp)*sqrt(2.0_wp*g)*0.611_wp*Junc(k)%WeirLength*Hoverflow**1.5_wp
    if (Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev < Junc(k)%OverflowElev) then
	  Junc(k)%Outflow=Qf
    else
      H2=Junc(Junc(k)%JuncPartner)%Head+Junc(Junc(k)%JuncPartner)%Elev-Junc(k)%OverflowElev
	  Junc(k)%Outflow=Qf*(1.0_wp-(H2/Hoverflow)**1.5_wp)**0.385_wp
    end if
  end select

if (WriteJuncDiagnostics) then
  do kk=1,NJuncDiag
    if ((T > JuncDiagTime(kk)) .and. (JuncDiagIndex(kk) == k)) then
      write (QDiagOutUnits(kk),'(f11.4,i5,15es16.8,i6)') T,Junc(k)%RatingCurveType,Hoverflow,H2,Qf,Junc(k)%Outflow
	end if
  end do
end if
  
99 return

end subroutine FindOutflow
    
    