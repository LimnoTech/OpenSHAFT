recursive subroutine GetLookupValue(LookupIn,LookupOut,TableID,XCol,YCol)
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Generalized table lookup for custom cross sections	
!%=====================================================================
use GlobalVariables

implicit none

integer, intent(in) :: TableID,XCol,YCol
integer(i4), save :: RowIndex,m1
real(wp), intent(in) :: LookupIn
real(wp), intent(inout) :: LookupOut
real(wp), save :: LocalSlope
!$omp threadprivate (RowIndex,m1,LocalSlope)

RowIndex=0
do m1=1,TableRowCount(TableID)
    if (LookupIn < CustomShapeTables(m1,XCol,TableID)) then
        RowIndex=m1
        exit
    end if
end do

if (RowIndex == 1) then
  LookupOut=CustomShapeTables(1,YCol,TableID)
else if (RowIndex == 0) then
  LookupOut=CustomShapeTables(TableRowCount(TableID),YCol,TableID)
else
  LocalSlope=(CustomShapeTables(RowIndex,YCol,TableID)-CustomShapeTables(RowIndex-1,YCol,TableID))/(CustomShapeTables(RowIndex,XCol,TableID)-CustomShapeTables(RowIndex-1,XCol,TableID))
  LookupOut=CustomShapeTables(RowIndex-1,YCol,TableID)+LocalSlope*(LookupIn-CustomShapeTables(RowIndex-1,XCol,TableID))
end if

return

end subroutine
