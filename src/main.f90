program OpenSHAFT
!%=====================================================================
!% OpenSHAFT release, version 1.0.0
!% 20240924
!% A one-dimensional finite volume solver for mixed-flow hydraulic systems
!% September 24, 2024
!%
!% Primary contact: Peter Klaver, email: pklaver@limno.com
!% 
!% Code authors
!% Jose Vasconselos -- original design and core calculations (in Delphi Pascal)
!% Peter Klaver -- conversion to Fortran, additional boundary conditions    
!%   
!% This code is made available under the GNU General Public License, version 3
!%     
!% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!% IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
!% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
!% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!% OTHER DEALINGS IN THE SOFTWARE.
!%
!%=====================================================================
  
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