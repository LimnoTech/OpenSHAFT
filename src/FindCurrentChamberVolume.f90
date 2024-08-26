subroutine FindCurrentChamberVolume

use GlobalVariables

implicit none

real(wp) :: slope
integer(i4) :: k1,index,id
logical :: exceed

exceed=.false.

Junc(k)%Lch=Junc(k)%OverflowElev-Junc(k)%Elev-Junc(k)%Head

if (Junc(k)%IsStorageCurve) then
  k1=Junc(k)%StorageMap
  do id=1,NStorageCurveLines(k1)
    if (Junc(k)%Lch <= ChamberVolH(id,k1)) then
	  index=id
	  exit
	end if
  end do
  if (Junc(k)%Lch > ChamberVolH(NStorageCurveLines(k1),k1)) exceed=.true.
end if
    
if (Junc(k)%IsStorageCurve) then
  if (Junc(k)%Lch <= 0.0_wp) then
    Junc(k)%VolChamber=ChamberVol(1,k1)
  else if (exceed) then
    Junc(k)%VolChamber=ChamberVol(NStorageCurveLines(k1),k1)
  else
    slope=(ChamberVol(index,k1)-ChamberVol(index-1,k1))/(ChamberVolH(index,k1)-ChamberVolH(index-1,k1))
    Junc(k)%VolChamber=ChamberVol(index-1,k1)+slope*(Junc(k)%Lch-ChamberVolH(index-1,k1))
  end if
else
  Junc(k)%VolChamber=Junc(k)%Lch*Junc(k)%CurrentArea
end if

if (Junc(k)%VolChamber < 0.0_wp) Junc(k)%VolChamber=0.0_wp

return
    
end subroutine FindCurrentChamberVolume