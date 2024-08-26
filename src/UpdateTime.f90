subroutine UpdateTime

use GlobalVariables
implicit none

!integer :: i,j
real(wp) :: tempthing,dTTemp

if (abs(V(1,1)) + c(1,1) == 0.0_wp) then
  dTTemp=0.5_wp
else
  dTTemp=CrF*Reach(1)%dX/(abs(V(1,1)) + c(1,1))
end if

do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    if (abs(V(1,1)) + c(1,1) == 0.0_wp) cycle
    if (dTTemp > CrF*Reach(j)%dX/(abs(V(i,j)) + c(i,j))) then
      dTTemp=CrF*Reach(j)%dX/(abs(V(i,j)) + c(i,j))
    end if
    if (isnan(A(i,j))) then
      tempthing=A(i,j)
    end if
  end do
end do

dT=dTTemp

if (dT > MaxDT) dT=MaxdT

return

end subroutine UpdateTime

