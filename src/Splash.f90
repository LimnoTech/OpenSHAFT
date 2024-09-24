subroutine Splash
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Description:
!% Writes a header to the console to keep the modeler interested
!%=====================================================================

use GlobalVariables, only: sj1,sj2,sj3

implicit none

write (*,*)
write (*,*)  '             _____  __  __ ___     ______ ______'
write (*,*)  '            / ___/ / / / //   |   / ____//_  __/'
write (*,*)  '            \__ \ / /_/ // /| |  / /_     / /   '
write (*,*)  '           ___/ // __  // ___ | / __/    / /    '
write (*,*)  '          /____//_/ /_//_/  |_|/_/      /_/     '
write (*,*)  '                                                '
write (*,*)  '=============================================================='
write (*,*)  '|           |          |          |     junction depths      |'
write (*,*)  '|  T (sec)  |   dT     | % error  |--------------------------|'
write (*,*)  '|           |          |          |   '//sj1//'   |   '//sj2//'   |   '//sj3//'   |'
write (*,*)  '=============================================================='
write (*,*)
write (*,*)

return

end subroutine