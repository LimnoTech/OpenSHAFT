subroutine FindCurrentJuncHead

use GlobalVariables

implicit none

integer(i4) :: TSIndex,TSId
real(wp) :: frac

TSId=Junc(k)%LevelTSMap
TSIndex=LevelTSIndex(TSId)

if (LevelTime(TSIndex,TSId) == T) then
    Junc(k)%Head=Level(TSIndex,TSId)
else if (LevelTime(TSIndex+1,TSId) <= T) then
    LevelTSIndex(TSId)=LevelTSIndex(TSId)+1
    TSIndex=TSIndex+1 
    frac=(T-LevelTime(TSIndex,TSId))/(LevelTime(TSIndex+1,TSId)-LevelTime(TSIndex,TSId))
	Junc(k)%Head=Level(TSIndex,TSId) + frac*(Level(TSIndex+1,TSId)-Level(TSIndex,TSId))
else if (LevelTime(TSIndex,TSId) < T) then
    frac=(T-LevelTime(TSIndex,TSId))/(LevelTime(TSIndex+1,TSId)-LevelTime(TSIndex,TSId))
	Junc(k)%Head=Level(TSIndex,TSId) + frac*(Level(TSIndex+1,TSId)-Level(TSIndex,TSId))
else
    write (*,*) 'level interpolation problem at junction',k
	interpoproblem=.true.
end if

99 return

end subroutine
