subroutine ReadSlideGates
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Reads slide gate data
!%
!% currently under development	
!%=====================================================================
    
use GlobalVariables

implicit none

integer(i4) :: kk

open (unit=17,file='SlideGates.inp',status='old',err=51)
do iskip=1,4
  read (17,*)
end do
do kk=1,NGates

    
end do
51 write (*,*) 'SlideGates.inp not found'
IOError=.true.
go to 999
   
999 continue   

end subroutine ReadSlideGates