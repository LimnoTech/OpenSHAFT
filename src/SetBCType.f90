subroutine SetBCType
!=======================================================
!
!  This subroutine sets the boundary condition type. Currently supported options are:
!
!    1   single upstream dropshaft
!    2   single downstream dropshaft
!    3   two-way dropshaft
!    4   three-way dropshaft
!    5   three-way split
!    6   slope change
!    7   level time series
!
!=======================================================

use GlobalVariables

implicit none

integer(i4) :: kt,kdummy,jj

NGates=0
do kt=1,Njuncs
  if (Junc(kt)%ReachCount == 1) then
    if (Junc(kt)%ReachCell(1) == 1) then
	  Junc(kt)%BCType=1
    else if (Junc(kt)%ReachCell(1) == Reach(Junc(kt)%ReachNo(1))%Ncells) then
      Junc(kt)%BCType=2
	else
	  Junc(kt)%BCType=0
    end if
  else if (Junc(kt)%ReachCount == 2) then
    if (Junc(kt)%ShaftArea == 0.0_wp) then
      if ((Reach(Junc(kt)%ReachNo(1))%D /= Reach(Junc(kt)%ReachNo(2))%D) .or. (Reach(Junc(kt)%ReachNo(1))%XSec /= Reach(Junc(kt)%ReachNo(2))%XSec) .or. (Junc(kt)%ReachElev(1) /= Junc(kt)%ReachElev(2))) then
        Junc(kt)%BCType=3
      else
        Junc(kt)%BCType=6
      end if
    else if (Junc(kt)%Option == -1) then
      Junc(kt)%BCType=7
      NGates=NGates+1
	else
	  Junc(kt)%BCType=3
	end if
  else if (Junc(kt)%ReachCount == 3) then
! PRK 12/3/2010 add check for splitter-type 3-way junction
    if ((Junc(kt)%ReachCell(2) == 1) .and. (Junc(kt)%ReachCell(3) == 1)) then
      Junc(kt)%BCType=5
    else
      Junc(kt)%BCType=4
    end if
! PRK 12/3/2010 end
!  else if (ReachCount(k) == 4) then
!    if ((JuncReachCell(k,3) == 1) .and. (JuncReachCell(k,4) == 1)) then
!      BCType(k)=18
!    else
!      BCType(k)=5
!    end if
!  else if (ReachCount(k) == 0) then
!    if (JuncOption(k) == -6) then
!      BCType(k)=15
!    else
!      BCType(k)=13
!    end if
  else
    Junc(kt)%BCType=0
  end if
  if ((Junc(kt)%BCType == 2) .and. (Junc(kt)%Option == -6)) Junc(kt)%BCType=7
end do

!> check to see that types have been assigned to each junction

do kt=1,Njuncs
  if (Junc(kt)%BCType == 0) then
    write (*,*) 'boundary condition not assigned for junction',kt
	IOError=.true.
  end if
end do

!> write system connectivity diagnostic file
open (unit=101,file='system_map.dia',status='replace')

write (101,*) 'JuncReachNo'
do kt=1,Njuncs
  write (101,'(i4,8i6)') kt,(Junc(kt)%ReachNo(jj),jj=1,Junc(kt)%ReachCount)
end do

write (101,*)
write (101,*) 'JuncReachCell'
do kt=1,Njuncs
  write (101,'(i4,8i6)') kt,(Junc(kt)%ReachCell(jj),jj=1,Junc(kt)%ReachCount)
end do

write (101,*)
write (101,*) 'JuncReachElev'
do kt=1,Njuncs
  write (101,'(i4,8f10.3)') kt,(Junc(kt)%ReachElev(jj),jj=1,Junc(kt)%ReachCount)
end do

write (101,*)
do jj=1,Nreaches
  write (101,'(i5,f8.4,2f14.4)') jj,Reach(jj)%dX
end do

write (101,*)
write (101,*) 'BCType'
do kt=1,Njuncs
  write (101,'(2i4)') kt,Junc(kt)%BCType
end do

close (unit=101)



return

end subroutine