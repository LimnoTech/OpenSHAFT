subroutine UpdateAQ
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Performs final step of St. Venant solution for interior reach cells	
!%=====================================================================
    
use GlobalVariables

implicit none

!do j=1,Nreaches
  if (not(Reach(j)%HasPocket)) then
    do i=2,Reach(j)%Ncells-1
      A(i,j)=A(i,j)-dT/Reach(j)%dX*(FintA(i,j)-FintA(i-1,j))+dT*SA(i,j)
      Q(i,j)=Q(i,j)-dT/Reach(j)%dX*(FintQ(i,j)-FintQ(i-1,j))+dT*SQ(i,j)
    end do
  else
    do i=2,Reach(j)%iFRT1-PocketBuffer-1
      A(i,j)=A(i,j)-dT/Reach(j)%dX*(FintA(i,j)-FintA(i-1,j))+dT*SA(i,j)
      Q(i,j)=Q(i,j)-dT/Reach(j)%dX*(FintQ(i,j)-FintQ(i-1,j))+dT*SQ(i,j)
    end do
    do i=Reach(j)%iFRT2+PocketBuffer+1,Reach(j)%Ncells-1
      A(i,j)=A(i,j)-dT/Reach(j)%dX*(FintA(i,j)-FintA(i-1,j))+dT*SA(i,j)
      Q(i,j)=Q(i,j)-dT/Reach(j)%dX*(FintQ(i,j)-FintQ(i-1,j))+dT*SQ(i,j)
    end do
    do i=Reach(j)%iFRT1-PocketBuffer,Reach(j)%iFRT2+PocketBuffer
      PocketZone(i,j)=.true.
    end do
    call UpdateAQPocket
  end if
!end do



return

end subroutine UpdateAQ