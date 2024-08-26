module GlobalVariables

use kind_parameter

implicit none
save

integer(i4), parameter :: wp = dp

type :: ReachObject
    real(wp) :: L
    real(wp) :: D
    real(wp) :: Apipe
    real(wp) :: So
    real(wp) :: n
    real(wp) :: dX
    real(wp) :: UpElev
    real(wp) :: DownElev
    real(wp) :: Width
    real(wp) :: Kup
    real(wp) :: Kdown
    real(wp) :: InitHGL
    real(wp) :: HFnflocal
    integer(i4) :: UpJunc
    integer(i4) :: DownJunc
    integer(i4) :: NCells
    integer(i4) :: XSec
    integer(i4) :: ReachShapeMap
    real(wp) :: Volume
    real(wp) :: VolOld
    logical :: HasPocket
    logical :: InitializePocket
    logical :: PocketFound
    integer(i4) :: iFRT1
    integer(i4) :: iFRT2
    real(wp) :: PocketStart
    real(wp) :: PocketReleaseT
    real(wp) :: VolPocket
    real(wp) :: VolPocketInit
    real(wp) :: MinPocketVol
    real(wp) :: ReleaseTime
    real(wp) :: Ha
    real(wp) :: LPocket
    real(wp) :: LColumn
    real(wp) :: QColumn
    real(wp) :: xDown
    real(wp) :: xUp
    real(wp) :: QColDiff
    real(wp) :: YUnder
    real(wp) :: QPocketDown
    real(wp) :: QPocketUp
end type

type :: JunctionObject
    real(wp) :: Elev
    real(wp) :: ShaftArea
    real(wp) :: OverflowElev
    real(wp) :: WeirLength
    real(wp) :: InflowInit
    real(wp) :: Inflow
    real(wp) :: Outflow
    real(wp) :: CurrentArea
    real(wp) :: Head
    real(wp) :: HeadOld
    real(wp) :: ReachElev(4)
    real(wp) :: InitHGL
    real(wp) :: Volume
    real(wp) :: VolOld
    real(wp) :: Closure
    real(wp) :: CumInf
    real(wp) :: CumOut
    real(wp) :: AvgFluxSum
    integer(i4) :: Option
    integer(i4) :: InflowMap
    integer(i4) :: StorageMap
    integer(i4) :: RatingCurveMap
    integer(i4) :: InflowControlMap
    integer(i4) :: GateMap
    integer(i4) :: ReachNo(4)
    integer(i4) :: ReachCount
    integer(i4) :: ReachCell(4)
    integer(i4) :: BCType
    integer(i4) :: JuncPartner
    integer(i4) :: RatingCurveType
    integer(i4) :: AirFluxCount
    integer(i4) :: AirFluxMap(4)
    integer(i4) :: LevelTSMap
    logical :: IsVented
    logical :: IsRatingCurve
    logical :: IsStorageCurve
    logical :: IsVariableInflow
    logical :: InitializeHch
    real(wp) :: Hch
    real(wp) :: Lch
    real(wp) :: VolChamber
    real(wp) :: Qair
    real(wp) :: Aopen
end type

type ControlObject
    integer(i4) :: InflowLoc
    integer(i4) :: TriggerNode
    real(wp) :: TriggerElev
    real(wp) :: ClosingTime
    real(wp):: MinInflow
    real(wp) :: MinOpenFrac
    integer(i4) :: FailFlag
    integer(i4) :: NControlSteps
    integer(i4) :: TriggerGroup
    logical :: Fails
    logical :: IsOpen
    logical :: IsClosing
    real(wp) :: StartCloseTime
    real(wp) :: EndCloseTime
    real(wp) :: InflowFrac
    real(wp) :: HSCNode(60)
    real(wp) :: HSCNodeAvg
end type

type ControlGroupObject
    integer(i4) :: NJuncs
    integer(i4) :: JuncID(5)
    real(wp) :: TriggerH(5)
end type

type(ReachObject), allocatable, dimension(:) :: Reach
type(JunctionObject), allocatable, dimension(:) :: Junc
type(ControlObject), allocatable, dimension(:) :: Control
type(ControlGroupObject), allocatable, dimension(:) :: CCGroup

real(wp) :: Crf,HFnf,Ratio_dx_Dpipe,initHfrac,Tsim,Twrite,StartWrite,StopWrite,amax,FloatDummy,tta,ReachVolTot,ReachVolTotOld,ReachVolTotWrite,JuncVolTot,InitialSystemVol,T,Rey,dT,Qlim,alpha1,alpha2,alpha4,ThresholdVol
real(wp), allocatable, dimension(:) :: QMult,WriteTimes,ControlEvalTimes
real(wp), allocatable, dimension(:,:) :: Inflow,InflowTime,TailwaterH,RatingCurveH,RatingCurveQ,RatingCurveQ2,StorageCurveH,StorageCurveA,AirUnitOpeningFlux,VentArea,Level,LevelTime
real(wp), allocatable, dimension(:,:,:) :: TWRatingCurveQ,CustomShapeTables
real(wp), parameter :: pi=3.1415926_wp
real(wp), parameter :: g=9.806_wp
real(wp), parameter :: MaxdT=0.05_wp
real(wp), parameter :: Nua=1.0e-6_wp
real(wp), parameter :: kss=1.0e-4_wp
real(wp), parameter :: tol=1.0e-3_wp
real(wp), parameter :: npoly=1.4_wp
real(wp), parameter :: atm=10.33_wp
real(wp), parameter :: RT=81264.5_wp
real(wp), parameter :: Cd=0.60_wp
real(wp), parameter :: epsilon=5.0e-3_wp

real(wp) :: maxcroechange,ceil,nn,ndum
character(len=2) :: sj1,sj2,sj3,sj1t,sj2t,sj3t,JuncD

real(wp), allocatable, dimension(:,:) :: z,y,V,hs,hc,A,Rh,Tfs,Q,Imom,c,Sf,FA,FQ,SA,SQ,FintA,FintQ,Aold,Qold,Ahat,Qhat,chat,f_c,Lambda1,Lambda2,dw1,dw2,MaxVacuum,MaxHGL,Qsave
real(wp), allocatable, dimension(:,:) :: hsold,hsWrite

integer(i4) :: i,j,k,nts,NThreadsUser
integer(i4) :: Njuncs,Nreaches,iskip,iread,iwrite,idummy,NShapeTables,NNN,MaxNCells,MaxNFlows,NInflowSeries,NRatingCurves,MaxNRatingCurveLines,NStorageCurves,MaxNStorageCurveLines,NAirUnits,NWrites,MaxNLevels
integer(i4) :: StartCondition,AirExhaustSwitch,AirPocketSwitch,PocketBuffer
integer(i4) :: TriggerStepCount,NInflowControl,NTriggerSteps,MaxControlSteps,NGates,JunctionCase,NControlEvalTimes,LastWriteNSteps,IControlEval,NLevelSer
integer(i4), allocatable, dimension(:) :: JuncSeriesID,Ntimes,TimeSeriesIndex,RatingCurveID,NTailwaterH,NRatingCurveLines,StorageCurveID,AirUnitNreaches,AirUnitNEndJunc
integer(i4), allocatable, dimension(:) :: NStorageCurveLines,TableRowCount,LevelSeriesID,NLevelTimes,LevelTSIndex
integer(i4), allocatable, dimension(:,:) :: AirUnitReaches,AirUnitEndJunc,OpeningMap,OpeningID,ThreeWayCases
integer(i4), allocatable, dimension(:,:,:) :: AirUnitIJ
integer(i4), parameter :: MinCellCt=11
integer(i4), parameter :: MinCellsPerThread=1500
integer, dimension(3) :: ShowJunc

!$omp threadprivate (i,j,k,nn,Rey,ndum,tta)

real(wp) :: dQ1,dQ2,dQ3,dQ4,dHa1,dHa2,dHa3,dHa4,dVa1,dVa2,dVa3,dVa4,dYp,f
!$omp threadprivate (dQ1,dQ2,dQ3,dQ4,dHa1,dHa2,dHa3,dHa4,dVa1,dVa2,dVa3,dVa4,dYp,f)

character(len=40) :: RunID,HotStartID


logical :: SomeRatingCurves,SomeStorageCurves,SomeVariableInflows,SomeInflowControls,IOError,IsCentralControl,CalcAirExhaust,TrapAirPockets,IsHGLInit,IsHotStart,WriteFilesOpen,interpoproblem,divergence,SomeSlideGates,SomeLevelSeries
logical, allocatable, dimension(:,:) :: PocketZone,vacuum
logical, allocatable, dimension(:) :: IsSurcharged          !defines whether an air unit is still able to exhaust air to a connected shaft; allocated to length NAirUnits

!> variables needed for junction diagnostics
real(wp), dimension(10) :: JuncDiagTime
integer(i4) :: NJuncDiag
integer(i4), dimension(10) :: JuncDiagIndex,DiagOutUnits,QDiagOutUnits
logical :: WriteJuncDiagnostics

!> variables for chamber volume storage curve calculations
real(wp), allocatable, dimension(:,:) :: ChamberVol
real(wp), allocatable, dimension(:,:) :: ChamberVolA
real(wp), allocatable, dimension(:,:) :: ChamberVolH

!> variables for continuity error calculation
real(wp) :: PercentVolError
real(wp) :: JuncInfTot
real(wp) :: JuncOutTot

!> variables for writing time-interpolated output
real(wp), allocatable, dimension(:) :: JuncHeadWrite,JuncInflowWrite,JuncOutflowWrite,JuncCumInfOld,JuncCumOutOld,JuncClosureOld,JuncInflowOld,JuncOutflowOld,HchOld,LchOld,VolChamberOld,QairOld,TempDumpBin
real(wp), allocatable, dimension(:) :: HchWrite,LchWrite,VolChamberWrite,QairWrite,AvgFluxSumWrite,JuncVolWrite,JuncInfWrite,JuncOutWrite,JuncClosureWrite,AvgFluxSumOld
real(wp), allocatable, dimension(:,:) :: Vold,yold,Awrite,Qwrite,Vwrite,ywrite
real(wp), allocatable, dimension(:) :: HaOld,VolPocketOld,QColumnOld,YUnderOld,QColDiffOld,LpocketOld,QPocketDownOld,QPocketUpOld
real(wp), allocatable, dimension(:) :: HaWrite,VolPocketWrite,QColumnWrite,YUnderWrite,QColDiffWrite,LpocketWrite,QPocketDownWrite,QPocketUpWrite

type :: WorkPacket
  integer :: jcount
  integer, dimension(100) :: jlist
end type

type(WorkPacket), allocatable, dimension(:) :: WorkByThread

end module GlobalVariables
