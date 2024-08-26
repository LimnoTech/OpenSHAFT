subroutine FindCurrentShaftArea

use GlobalVariables

implicit none

integer(i4) :: k1,id,index
real(wp) :: slope
logical :: exceed

exceed=.false.

if (Junc(k)%IsStorageCurve) then
  k1=Junc(k)%StorageMap
  do id=1,NStorageCurveLines(k1)
    if (Junc(k)%Head <= StorageCurveH(id,k1)) then
	  index=id
	  exit
	end if
  end do
  if (Junc(k)%Head > StorageCurveH(NStorageCurveLines(k1),k1)) exceed=.true.
end if

if (Junc(k)%IsStorageCurve) then
  if (Junc(k)%Head <= 0.0_wp) then
    Junc(k)%CurrentArea=StorageCurveA(1,k1)
  else if (exceed) then
    Junc(k)%CurrentArea=StorageCurveA(NStorageCurveLines(k1),k1)
  else
    slope=(StorageCurveA(index,k1)-StorageCurveA(index-1,k1))/(StorageCurveH(index,k1)-StorageCurveH(index-1,k1))
    Junc(k)%CurrentArea=StorageCurveA(index-1,k1)+slope*(Junc(k)%Head-StorageCurveH(index-1,k1))
  end if
else
  Junc(k)%CurrentArea=Junc(k)%ShaftArea
end if

return

end subroutine

