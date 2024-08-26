subroutine SetInflowControl

use GlobalVariables

implicit none

real(wp), save :: Tstar,MinInflowFrac
integer(i4), save :: InflowControlCase,kk
!$omp threadprivate (Tstar,MinInflowFrac,InflowControlCase,kk)


kk=Junc(k)%InflowControlMap
InflowControlCase=0

MinInflowFrac=((-0.6765_wp*(1.0_wp-Control(kk)%MinOpenFrac)-0.035_wp)*(1.0_wp-Control(kk)%MinOpenFrac)-0.1509_wp)*(1.0_wp-Control(kk)%MinOpenFrac)+0.8686_wp

if ((Control(kk)%IsOpen) .and. (not(Control(kk)%IsClosing))) InflowControlCase=1
if ((Control(kk)%IsOpen) .and. (Control(kk)%IsClosing)) InflowControlCase=2
if (not(Control(kk)%IsOpen)) InflowControlCase=3

select case (InflowControlCase)
  
  case(1)
    Control(kk)%InflowFrac=1.0_wp
  case(2)
      Tstar=(T-Control(kk)%StartCloseTime)/Control(kk)%ClosingTime
      Control(kk)%InflowFrac=((-0.6765_wp*Tstar-0.035_wp)*Tstar-0.1509_wp)*Tstar+0.8686_wp
!      JuncInflowFrac(kk)=(ControlEndCloseTime(kk)-T)/ControlClosingTime(kk)
    if (Control(kk)%InflowFrac < max(MinInflowFrac,0.001_wp)) then
      Control(kk)%InflowFrac=MinInflowFrac
      Control(kk)%IsOpen=.false.
      Control(kk)%IsClosing=.false.
      write (609,*) 'control at junction',k,'closed completely at',T
    end if
  case(3)
    Control(kk)%InflowFrac=max(0.0_wp,MinInflowFrac)
  case default
    write (*,*) 'undefined inflow control case...uh oh'
    stop
  
end select
    
return

end subroutine