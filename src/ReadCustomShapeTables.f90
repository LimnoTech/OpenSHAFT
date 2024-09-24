subroutine ReadCustomShapeTables
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Sets up custom shape tables	
!%=====================================================================

use GlobalVariables

implicit none

integer :: MaxRowCount,Ntemp,jtemp
MaxRowCount=0
allocate(TableRowCount(NShapeTables))
Reach(:)%ReachShapeMap=0
TableRowCount=0

open (unit=12,file='CustomShapeTables.inp',status='old',err=2000)
do iskip=1,4
    read (12,*)
end do
read (12,*) Ntemp
if (Ntemp /= NShapeTables) then
    write (*,*) 'number of shape tables is different from number of reaches requiring custom shapes'
    write (*,*) 'one table for each shape, please'
    stop
end if

do i=1,NShapeTables
    read (12,*) jtemp,TableRowCount(i),Reach(jtemp)%Apipe
    if (TableRowCount(i) > MaxRowCount) MaxRowCount=TableRowCount(i)
    do iskip=1,TableRowCount(i)
        read (12,*)
    end do
end do
write (*,*) 'MaxRowCount=',MaxRowCount
write (*,*)
close(12)

allocate(CustomShapeTables(MaxRowCount,6,NShapeTables))
CustomShapeTables=0.0_wp
open (unit=12,file='CustomShapeTables.inp',status='old',err=2000)
do iskip=1,4
    read (12,*)
end do
read (12,*) Ntemp
do iread=1,NShapeTables
    read (12,*) jtemp,TableRowCount(iread),Reach(jtemp)%Apipe
    Reach(jtemp)%ReachShapeMap=iread
    do i=1,TableRowCount(iread)
        read (12,*) (CustomShapeTables(i,j,iread),j=1,6)
    end do
end do

go to 3000

2000 write (*,*) 'CustomShapes.inp missing'
     stop
     
3000 continue

return

end subroutine
