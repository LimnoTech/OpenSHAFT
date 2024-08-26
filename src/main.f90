program OpenSHAFT

use GlobalVariables
use GlobalFunctions, only: wtime
use omp_lib

implicit none

integer :: seconds,ndummy,NThreads

seconds = wtime()

call ReadInput

if (IOError) call AbnormalStop

call SetBCType

if (IOError) call AbnormalStop

call SetOutputFiles

call InitializeVariables

if (CalcAirExhaust) call MakeAirUnits

call SetInitialConditions

call Splash

!ndummy=min(NNN/MaxNCells,NNN/MinCellsPerThread)
NThreads=min(omp_get_num_procs(),NThreadsUser)
call omp_set_num_threads(NThreads)
allocate (WorkByThread(NThreads))

do i=1,NThreads
  WorkByThread(i)%jcount=0
  WorkByThread(i)%jlist(:)=0
end do

call WorkBalance(NThreads)

if (IOError) call AbnormalStop

call ExecuteTimeLoop

open (unit=115,form='formatted',file=trim(RunID)//'.MaxVacuum.csv',status='replace')
write (115,'(A)') 'cell,reach,invert,crown,Minhs,MaxHGL'
do j=1,Nreaches
  do i=1,Reach(j)%Ncells
    write (115,239) i,j,z(i,j),z(i,j)+Reach(j)%D,MaxVacuum(i,j),MaxHGL(i,j)
  end do
end do
239 format (i4,',',i6,4(',',es12.5))

seconds=wtime()-seconds

write (*,*)
write (*,*) 'Elapsed time = ',seconds,' seconds'
write (*,*)

end program    