subroutine ReadSlideGates
    
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