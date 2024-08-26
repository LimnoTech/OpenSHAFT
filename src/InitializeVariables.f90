subroutine InitializeVariables
    
use GlobalVariables

implicit none

allocate (z(MaxNCells,Nreaches),y(MaxNCells,Nreaches),V(MaxNCells,Nreaches),hs(MaxNCells,Nreaches))
allocate (hc(MaxNCells,Nreaches),A(0:MaxNCells+1,Nreaches),Q(0:MaxNCells+1,0:Nreaches+1),Imom(MaxNCells,Nreaches))
allocate (Tfs(MaxNCells,Nreaches),Rh(MaxNCells,Nreaches),c(MaxNCells,Nreaches),Sf(MaxNCells,Nreaches))
allocate (FA(MaxNCells,Nreaches),FQ(MaxNCells,Nreaches),SA(MaxNCells,Nreaches),SQ(MaxNCells,Nreaches))
allocate (FintA(MaxNCells,Nreaches),FintQ(MaxNCells,Nreaches),Aold(0:MaxNCells+1,Nreaches),Qold(MaxNCells,0:Nreaches+1),Qsave(MaxNCells,0:Nreaches+1))
allocate (Ahat(MaxNCells,Nreaches),Qhat(MaxNCells,Nreaches),chat(MaxNCells,Nreaches),f_c(0:MaxNCells+1,Nreaches))
allocate (Lambda1(MaxNCells,Nreaches),Lambda2(MaxNCells,Nreaches),dw1(MaxNCells,Nreaches),dw2(MaxNCells,Nreaches))
allocate (hsold(MaxNCells,Nreaches),hsWrite(MaxNCells,Nreaches))
allocate (vacuum(0:MaxNCells+1,Nreaches),MaxVacuum(MaxNCells,Nreaches),MaxHGL(MaxNCells,Nreaches),PocketZone(MaxNCells,Nreaches))
allocate (ThreeWayCases(NJuncs,20))

allocate (JuncHeadWrite(NJuncs),JuncInflowWrite(NJuncs),JuncOutflowWrite(NJuncs),Vold(MaxNCells,Nreaches),yold(MaxNCells,Nreaches))
allocate (JuncCumInfOld(NJuncs),JuncCumOutOld(NJuncs),JuncClosureOld(NJuncs),JuncInflowOld(NJuncs),JuncOutflowOld(NJuncs),HchOld(NJuncs),LchOld(NJuncs),VolChamberOld(NJuncs),QairOld(NJuncs))
allocate (HaOld(Nreaches),VolPocketOld(Nreaches),QColumnOld(Nreaches),YUnderOld(Nreaches),QColDiffOld(Nreaches),LpocketOld(Nreaches),QPocketDownOld(Nreaches),QPocketUpOld(Nreaches))
allocate (HchWrite(NJuncs),LchWrite(NJuncs),VolChamberWrite(NJuncs),QairWrite(NJuncs),AvgFluxSumWrite(NJuncs))
allocate (Awrite(0:MaxNCells+1,Nreaches),Qwrite(0:MaxNCells+1,0:Nreaches+1),Vwrite(MaxNCells,Nreaches),ywrite(MaxNCells,Nreaches))
allocate (HaWrite(Nreaches),VolPocketWrite(Nreaches),QColumnWrite(Nreaches),YUnderWrite(Nreaches),QColDiffWrite(Nreaches),LpocketWrite(Nreaches),QPocketDownWrite(Nreaches),QPocketUpWrite(Nreaches))
allocate (JuncVolWrite(NJuncs),JuncInfWrite(NJuncs),JuncOutWrite(NJuncs),JuncClosureWrite(NJuncs),AvgFluxSumOld(NJuncs))

PocketZone=.false.
interpoproblem=.false.
vacuum=.false.
divergence=.false.

Reach(:)%VolOld=0.0_wp
Reach(:)%HasPocket=.false.
Reach(:)%InitializePocket=.false.
Reach(:)%MinPocketVol=ThresholdVol
Reach(:)%VolPocket=0.0_wp
Reach(:)%ReleaseTime=10.0_wp
Reach(:)%PocketReleaseT=0.0_wp

Junc(:)%Inflow=0.0_wp
Junc(:)%Outflow=0.0_wp
Junc(:)%CurrentArea=0.0_wp
Junc(:)%Head=0.0_wp
Junc(:)%HeadOld=0.0_wp
Junc(:)%InitHGL=0.0_wp
Junc(:)%Volume=0.0_wp
Junc(:)%VolOld=0.0_wp
Junc(:)%AvgFluxSum=0.0_wp
Junc(:)%CumInf=0.0_wp
Junc(:)%CumOut=0.0_wp
Junc(:)%Closure=0.0_wp
!Junc(:)%JuncPartner=0
Junc(:)%Hch=0.0_wp
Junc(:)%Lch=0.0_wp
Junc(:)%VolChamber=0.0_wp
Junc(:)%Qair=0.0_wp
Junc(:)%Aopen=0.1_wp
Junc(:)%InitializeHch=.true.

Q=0.0_wp
Qold=0.0_wp
Qsave=0.0_wp
Qwrite=0.0_wp
hsold=0.0_wp
hsWrite=0.0_wp

ThreeWayCases=0

return

end subroutine