!******************************************************************************
!****m* /decisionScore
! NAME
! module decisionScore
! PURPOSE
! feature-022: emit the log-probability of individual random Monte-Carlo
! decisions taken during a run (decay hazard/branching, cascade collision
! acceptance, 2-body phase-space weight, neutrino-vertex channel choice),
! gated by inputGeneral/emitDecisionLogp (default .false. => no-op, no
! change to output or performance). Consumed by feature-022's PYTHIA/DUNE
! decision-count-audit pipeline; see docs/reports/2026-08-07-feature022-a4-decision-count-audit.md.
!******************************************************************************
module decisionScore

  implicit none
  private
  public :: writeDecisionLogp

contains

  !****************************************************************************
  !****s* decisionScore/writeDecisionLogp
  ! NAME
  ! subroutine writeDecisionLogp(itemTag, firstEvent, logp, outcome,
  ! sigmaPreClamp, kinPrefactor, catVec)
  ! PURPOSE
  ! Append one record to decisionLogProb.dat: which instrumented decision
  ! class (itemTag), which cascade (firstEvent, matches collisionList.txt/
  ! ROOT-ntuple first_event), the log-probability of the outcome actually
  ! sampled (logp), which branch was taken (outcome, issue-073), and
  ! optionally the sufficient statistic needed to recompute the decision
  ! at a different theta post-hoc (sigmaPreClamp/kinPrefactor, issue-073) or
  ! the full categorical distribution over alternatives that was sampled
  ! from (catVec).
  ! INPUTS
  ! * character(len=*)      :: itemTag     -- short tag for the decision class
  ! * integer                :: firstEvent -- cascade id (perturbative particle)
  ! * real                   :: logp       -- log-probability of the outcome taken
  ! * integer, optional      :: outcome    -- branch discriminator, meaning is
  !   itemTag-specific (0 where not applicable, kept uniform across tags):
  !   collAcc: 0 reject, 1 accept (genuine draw), 2 forced accept (p>1,
  !   logp=0 exactly); decayHaz: 0 no decay, 1 decay.
  ! * real, optional         :: sigmaPreClamp -- collAcc only: the pre-clamp
  !   cross section sigmaTot*sigmaBoost [mb], i.e. before min(350.,.), so
  !   p(theta') = min(350., sigmaPreClamp * s'/s) * kinPrefactor is exact and
  !   the 350 mb saturation rate is measurable from this column alone.
  ! * real, optional         :: kinPrefactor  -- collAcc only: the kinematic
  !   prefactor vRel*deltaT*weightLocal/(10*numEnsembles*d^3x) multiplying
  !   min(350.,sigmaPreClamp) to form the accept probability.
  ! * real, dimension(:), optional :: catVec -- full categorical distribution
  ! NOTES
  ! No-op unless inputGeneral/emitDecisionLogp is set; the file is opened
  ! lazily on first use so a default run never creates decisionLogProb.dat.
  ! Column layout is uniform across itemTag so the file can be parsed with a
  ! single fixed-width reader: itemTag firstEvent logp outcome sigmaPreClamp
  ! kinPrefactor [catVec...]. sigmaPreClamp/kinPrefactor read 0. for every
  ! tag except collAcc (0. is unambiguous: both quantities are strictly
  ! positive whenever they apply).
  !****************************************************************************
  subroutine writeDecisionLogp(itemTag, firstEvent, logp, outcome, &
       sigmaPreClamp, kinPrefactor, catVec)
    use inputGeneral, only: emitDecisionLogp

    character(len=*), intent(in) :: itemTag
    integer, intent(in) :: firstEvent
    real, intent(in) :: logp
    integer, intent(in), optional :: outcome
    real, intent(in), optional :: sigmaPreClamp
    real, intent(in), optional :: kinPrefactor
    real, dimension(:), intent(in), optional :: catVec

    integer, parameter :: iFile = 309
    logical, save :: isOpen = .false.
    integer :: i, outcomeOut
    real :: sigmaOut, prefOut

    if (.not. emitDecisionLogp) return

    if (.not. isOpen) then
       open(file="decisionLogProb.dat", UNIT=iFile, Status='unknown', &
            Position='Append', Action='Write')
       isOpen = .true.
    end if

    outcomeOut = 0
    if (present(outcome)) outcomeOut = outcome
    sigmaOut = 0.
    if (present(sigmaPreClamp)) sigmaOut = sigmaPreClamp
    prefOut = 0.
    if (present(kinPrefactor)) prefOut = kinPrefactor

    if (present(catVec)) then
       write(iFile,'(A8,1X,I9,1X,ES14.6,1X,I2,1X,ES14.6,1X,ES14.6)', advance='no') &
            itemTag, firstEvent, logp, outcomeOut, sigmaOut, prefOut
       do i=1,size(catVec)
          write(iFile,'(1X,ES14.6)', advance='no') catVec(i)
       end do
       write(iFile,*)
    else
       write(iFile,'(A8,1X,I9,1X,ES14.6,1X,I2,1X,ES14.6,1X,ES14.6)') &
            itemTag, firstEvent, logp, outcomeOut, sigmaOut, prefOut
    end if

  end subroutine writeDecisionLogp

end module decisionScore
