subroutine Splash

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