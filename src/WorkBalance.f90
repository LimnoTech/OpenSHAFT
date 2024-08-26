subroutine WorkBalance(NThreads)

use GlobalVariables, only: WorkByThread,Reach,NNN,NReaches,NThreadsUser,IOError

implicit none
integer, intent(in) :: NThreads
integer :: TotalCells,chunk,tempsum,FirstReach,NextReach,dif,mindif,savemin,jj,MinReach,it,jt,jchk,ReachesLeft
integer, allocatable, dimension(:) :: CellCount
logical, allocatable, dimension(:) :: ReachNotUsed

TotalCells=0
chunk=0
dif=0
mindif=0
savemin=0
FirstReach=0
NextReach=0
MinReach=0

allocate (CellCount(NThreads))
allocate (ReachNotUsed(Nreaches))

CellCount=0
ReachNotUsed=.true.
chunk=NNN/NThreads

do it=1,NThreads-1
  tempsum=0
  do jt=1,Nreaches
    if (ReachNotUsed(jt)) then
      FirstReach=jt
      ReachNotUsed(jt)=.false.
      exit
    end if
  end do
  jj=1
  WorkByThread(it)%jlist(jj)=FirstReach
  WorkByThread(it)%jcount=jj
  tempsum=Reach(FirstReach)%Ncells
  mindif=abs(Reach(FirstReach)%Ncells-chunk)
  savemin=mindif
  do
    jj=jj+1
    do jt=1,Nreaches
      if (ReachNotUsed(jt)) then
        NextReach=jt
        exit
      end if
    end do
    do jt=NextReach,Nreaches
      if (ReachNotUsed(jt)) then
        dif=abs(tempsum+Reach(jt)%Ncells-chunk)
        if (dif < mindif) then
          mindif=dif
          MinReach=jt
        end if
      else
        cycle
      end if
    end do
    if (mindif < savemin) then
      ReachNotUsed(MinReach)=.false.
      WorkByThread(it)%jlist(jj)=MinReach
      WorkByThread(it)%jcount=jj
      savemin=mindif
      tempsum=tempsum+Reach(MinReach)%Ncells
    else
      exit
    end if
  end do
! before going on to next thread, check if all reaches have been assigned; if so, exit because the user has chosen too many threads for a balanced work load
  ReachesLeft=0
  do jchk=1,NReaches
      if (ReachNotUsed(jchk)) ReachesLeft=ReachesLeft+1
  end do
  if (ReachesLeft == 0) then
      write (*,*) 'You selected',NThreadsUser,' threads, but that is too many for optimal work balancing.'
      write (*,*) 'Please select a smaller number of threads.'
      IOError=.true.
  end if
end do

! pack last thread with remaining reaches

jj=0
do jt=1,Nreaches
  if (ReachNotUsed(jt)) then
    ReachNotUsed(jt)=.false.
    jj=jj+1
    WorkByThread(Nthreads)%jlist(jj)=jt
    WorkByThread(Nthreads)%jcount=jj
  end if
end do

open (12,file='WorkBalance.out',status='replace')
write (12,*) 'there will be',NThreads,' threads executed in the parallel region'

do it=1,Nthreads
  do jt=1,WorkByThread(it)%jcount
    CellCount(it)=CellCount(it)+Reach(WorkByThread(it)%jlist(jt))%Ncells
  end do
  write (12,*) CellCount(it),(WorkByThread(it)%jlist(jt),jt=1,WorkByThread(it)%jcount)
end do

close (12)

return

end subroutine WorkBalance
