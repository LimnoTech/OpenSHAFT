subroutine AbnormalStop

implicit none

write (*,*) 'simulation not performed'
write (*,*) 'there were input errors'
write (*,*) 'please see messages above, which are meant to be helpful'

stop

end subroutine AbnormalStop