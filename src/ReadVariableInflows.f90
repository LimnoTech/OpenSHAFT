subroutine ReadVariableInflows

use GlobalVariables

implicit none

integer(i4) :: kt,kdummy

!> scan inflow time series file for number of series and longest series
MaxNFlows=0
NInflowSeries=0
open (unit=13,file='JuncInflows.inp',status='old',err=903)
do iskip=1,6
  read (13,*)
end do
do while(not(eof(13)))
  read (13,*) idummy,idummy,FloatDummy
  if (idummy > MaxNFlows) MaxNFlows=idummy
  NInflowSeries=NInflowSeries+1
  do iread=1,idummy
    read (13,*) FloatDummy,FloatDummy
  end do
end do
!write (*,*) 'MaxNFlows =',MaxNFlows
!write (*,*)
close (unit=13)

!> allocate inflow time series arrays then read data
allocate(InflowTime(MaxNFlows,NInflowSeries),Inflow(MaxNFlows,NInflowSeries),Ntimes(NInflowSeries),QMult(NInflowSeries),JuncSeriesID(NInflowSeries))
open (unit=13,file='JuncInflows.inp',status='old',err=903)
do iskip=1,6
  read (13,*)
end do
do iread=1,NInflowSeries
!  write (*,*)
  read (13,*) JuncSeriesID(iread),Ntimes(iread),QMult(iread)
  do kt=1,Ntimes(iread)
    read (13,*) InflowTime(kt,iread),Inflow(kt,iread)
	Inflow(kt,iread)=Inflow(kt,iread)*QMult(iread)
  end do
!  write (*,612) Ntimes(iread),JuncSeriesID(iread)
end do
close (unit=13)

!> initialize time series index for interpolation
allocate (TimeSeriesIndex(NInflowSeries))
TimeSeriesIndex=1

!> assign time series to proper junctions

do kt=1,NInflowSeries
  do kdummy=1,Njuncs
    if (JuncSeriesID(kt) == kdummy) then
	  Junc(kdummy)%InflowMap=kt
      exit
	end if
  end do
end do
    
!> check to see if inflow time series are long enough
do kdummy=1,Njuncs
  if (Junc(kdummy)%InflowMap /= 0) then
    if (InflowTime(Ntimes(Junc(kdummy)%InflowMap),Junc(kdummy)%InflowMap) < Tsim) then
      write (*,*) 'not enough inflow data for junction',kdummy
	  IOError=.true.
    end if
  end if
end do

!> check to see if there are time series for each junction that expects one

do kdummy=1,Njuncs
  if ((Junc(kdummy)%IsVariableInflow) .and. (Junc(kdummy)%InflowMap == 0)) then
    write (*,*) 'inflow time series for junction',kdummy,' is missing'
	IOError=.true.
  end if
end do

go to 999

903 write(*,*) 'JuncInflows.inp not found'
IOError=.true.

999 continue
    
return

612 format ('+read ',i6,' entries for junction ',i5)

end subroutine