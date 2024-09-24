subroutine AbnormalStop
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!%=====================================================================

implicit none

write (*,*) 'simulation not performed'
write (*,*) 'there were input errors'
write (*,*) 'please see messages above, which are meant to be helpful'

stop

end subroutine AbnormalStop