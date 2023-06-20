!******************************************************************************
!****m* /neutrino2piBack
! NAME
! module neutrino2piBack
!
! PURPOSE
! provide some routines for generating a 2 pion background for neutrino
! induced events
!******************************************************************************
module neutrino2piBack
  use particleDefinition
  use eN_eventDefinition
  use Neutrinoparms

  implicit none

  private

  public :: DoNu2piBack
  public :: BGStruct
  public :: BGDualXS

 
	
 contains   

  !****************************************************************************
  !****s* neutrino2piBack/Donu2PiBack
  ! NAME
  ! subroutine DoNu2piBack(eN,outPart,XS,new2piBG,indBG)
  !
  ! PURPOSE
  ! This subroutine produces the 2pi background cross section in nanobarn
  ! in the mass range between 1.3 GeV and 2.0 GeV
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- the Lepton-Nucleon Event
  ! * logical  :: new2piBG - switch for using new 2piBG calculation (if T)
  ! * integer :: indBG - Index to choose BG-like background parametrization
  !
  ! OUTPUT
  ! * real  :: XS -- calculated cross section dsigma/(dOmega dE')
  !            in nb/GeV for electrons, in 10^(-38) cm^2)/GeV for neutrinos
  ! *       :: Outpart -- outgoing articles: 1 nucleon, 2 pions
  !****************************************************************************
  subroutine DoNu2piBack(eN,outPart,XS,new2piBG,indBG)

    use constants, only: pi, mN, mPi, alphaQED, GF
    use IdTable, only: pion,nucleon
    use eventGenerator_eN_lowEnergy, only: init_2Pi_getBG
    use eN_event, only: eNeV_GetKinV
    use ParamEP, only: CalcParamEP
    use nBodyPhaseSpace, only: momenta_in_3BodyPS
    use energyCalc, only: energyCorrection
    use mediumDefinition
    use mediumModule, only: mediumAt
    use monteCarlo, only: monteCarloChoose
    use offShellPotential, only: setOffShellParameter
    use Neutrinoparms

    type(electronNucleon_event), intent(in) :: eN
    type(particle),dimension(:), intent(out) :: OutPart
    real, intent(out) :: XS
    logical, intent(in) :: new2piBG
    integer, intent(in) :: indBG

    real :: nu,Q2,W,Wfree,eps,fT,Wrec,x
    integer :: i,qnuk, iProc
    real, dimension(0:3) :: sig2Pi
    real :: dum1,dum2
    real, dimension(1:3,1:3) :: p3       ! momenta of three particles
    logical :: flagOK
    real, dimension(0:3) :: ptot
    real, dimension(1:3) :: betaCMToLab = 0.0
    type(medium)         :: mediumAtPosition
    real :: sumWeight
    real :: Q2fak
    real :: sig2piBG

    real, dimension(0:3)        :: mom_in,mom_out
    integer :: j
	
	real, dimension(0:1, -3:3, 3), save:: RelFak2pi ! (qnuk,eN%idProcess,iProc)
    integer, dimension(0:1, -3:3, 3, 3), save:: Proc2piCharge


    ! Relfak2pi determines the probability for a given 2pi final state,
    ! here taken for simplicity to be the same (1/3). The first index in
    ! Relfak2pi describes the charge of the hit nucleon (n,p), the second
    ! index the processID (EM, CC, or NC)

    ! In Proc2piCharge on the lhs the first index gives the charge
    ! of the hit nucleon, the second index describes the process
    ! (1 = em, 2 = CC, 3 = NC, negative signs for antiparticles)
    ! the third index numbers the process
    ! On the rhs the first index is the charge of the outgoing nucleon
    ! the second and third indices describe the charges of the two
    ! produced pions

    ! EM:

    RelFak2pi(0, 1, 1:3) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(0, 1, 1, 1:3) = (/0, 0, 0/)
    Proc2piCharge(0, 1, 2, 1:3) = (/0, 1,-1/)
    Proc2piCharge(0, 1, 3, 1:3) = (/1, 0,-1/)

    RelFak2pi(1, 1, 1:3) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(1, 1, 1, 1:3) = (/0, 1, 0/)
    Proc2piCharge(1, 1, 2, 1:3) = (/1, 1,-1/)
    Proc2piCharge(1, 1, 3, 1:3) = (/1, 0, 0/)

    ! antiEM
    RelFak2pi(0, -1, 1:3) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(0, -1, 1, 1:3) = (/0, 1,-1/)
    Proc2piCharge(0, -1, 2, 1:3) = (/0, 0, 0/)
    Proc2piCharge(0, -1, 3, 1:3) = (/1, 0,-1/)

    RelFak2pi(1, -1, 1:3) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(1, -1, 1, 1:3) = (/0, 1, 0/)
    Proc2piCharge(1, -1, 2, 1:3) = (/1, 1,-1/)
    Proc2piCharge(1, -1, 3, 1:3) = (/1, 0, 0/)

    ! CC:

    RelFak2pi( 0, 2, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(0, 2, 1, 1:3) = (/0, 1, 0/)
    Proc2piCharge(0, 2, 2, 1:3) = (/1, 0, 0/)
    Proc2piCharge(0, 2, 3, 1:3) = (/1, 1,-1/)

    RelFak2pi( 1, 2, 1:3 ) = (/ 1./2., 1./2., 0. /)
    Proc2piCharge(1, 2, 1, 1:3) = (/0, 1, 1/)
    Proc2piCharge(1, 2, 2, 1:3) = (/1, 1, 0/)
    Proc2piCharge(1, 2, 3, 1:3) = (/9, 9, 9/) ! does not exist


    ! antiCC:

    RelFak2pi( 0,-2, 1:3 ) = (/ 1./2., 1./2., 0. /)
    Proc2piCharge(0,-2, 1, 1:3) = (/0,-1, 0/)
    Proc2piCharge(0,-2, 2, 1:3) = (/1,-1,-1/)
    Proc2piCharge(0,-2, 3, 1:3) = (/9, 9, 9/) ! does not exist

    RelFak2pi( 1,-2, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(1,-2, 1, 1:3) = (/0, 0, 0/)
    Proc2piCharge(1,-2, 2, 1:3) = (/0, 1,-1/)
    Proc2piCharge(1,-2, 3, 1:3) = (/1,-1, 0/)

    ! NC
    RelFak2pi( 0, 3, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(0, 3 , 1, 1:3) = (/0, 0, 0/)
    Proc2piCharge(0, 3 , 2, 1:3) = (/0, 1,-1/)
    Proc2piCharge(0, 3 , 3, 1:3) = (/1,-1, 0/)
	
	RelFak2pi( 1, 3, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(1, 3 , 1, 1:3) = (/1, 0, 0/)
    Proc2piCharge(1, 3 , 2, 1:3) = (/1, 1,-1/)
    Proc2piCharge(1, 3 , 3, 1:3) = (/0, 1, 0/)

    ! antiNC
    RelFak2pi( 0,-3, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(0,-3, 1, 1:3) = (/0, 0, 0/)
    Proc2piCharge(0,-3, 2, 1:3) = (/0, 1,-1/)
    Proc2piCharge(0,-3, 3, 1:3) = (/1, 0, 0/)
	
	RelFak2pi( 1,-3, 1:3 ) = (/ 1./3., 1./3., 1./3. /)
    Proc2piCharge(1,-3, 1, 1:3) = (/1, 0, 0/)
    Proc2piCharge(1,-3, 2, 1:3) = (/1, 1,-1/)
    Proc2piCharge(1,-3, 3, 1:3) = (/0, 1, 0/)

    flagOK = .false.
    XS = 0.0
    OutPart%ID=0
    OutPart%perweight=0.0

    call eNeV_GetKinV(eN, nu,Q2,W,Wfree,eps,fT)
    ! corrected 2piBG determined only for 1.3 <W < 2.2 GeV region,
    ! turns on roughly at 2pi threshold
    x = Q2/(2*mN*nu)
    Wrec = sqrt(mN**2 + 2*mN*nu - Q2)
    if (W < 1.30 .or. W > 2.2)  return

    if (new2piBG) then

       call BGDualXS(eN,sig2piBG,indBG)	! cross section sig2piBG in nanobarn
       if (sig2PiBG.lt.10e-20) return
       ! Now attenuation of background X-section to account for the fact
       ! that 1pi BG is separately taken into account,
       ! Fit function Att is available only for Christy-Bosted (indBG=3)
       if (indBG==3) then
          XS = Att(Wrec,x) * sig2piBG
       else
          XS = sig2PiBG
       end if

    else

       if (eN%QSquared >= 5.0) return
       call init_2Pi_getBG(eN%nucleon_free, Wfree, sig2pi)
       if (sig2Pi(0) <= 10e-20) return

       !===== Scale the event acording Q2:

       XS = sig2pi(0)

       call CalcParamEP(Wfree,0.0,0.0, dum1)
       call CalcParamEP(Wfree,Q2, eps, dum2)
       !    call CalcParamEP(Wfree,Q2, 0.0, dum2)
       XS = XS * dum2/dum1

       !=======================================================================
       XS = XS*fT/ ( 1e3* pi/(eN%lepton_out%momentum(0)*eN%lepton_in%momentum(0)))

       XS = XS * (GF/(2*pi*alphaQED))**2 ! scaling according DIS cross sections
       XS = XS*Norm2pi

       Q2fak = Q2/(1+Q2/(80.399**2)) ! mW=80.399 GeV

       !    XS = XS*Q2 ! DUMMY !!!
       XS = XS*(pi**2/2)*Q2fak**2  ! DUMMY !!!

    end if

    !===== 2: Select Process
    call setToDefault(OutPart)
    qnuk = eN%nucleon_free%charge
    iProc = monteCarloChoose(RelFak2pi(qnuk,eN%idProcess,1:3), sumWeight)
    if (iProc.eq.0) then
       write(*,*) 'neutrino 2piBack.f90:  Failure!'
       write(*,*) 'qnuk=',qnuk,'Proc=',eN%idProcess,'RelFak2pi=',RelFak2pi(qnuk,eN%idProcess,1:3)
       stop
    end if
    OutPart(1:3)%charge = Proc2piCharge(qnuk,eN%idProcess,iProc,1:3)

    XS = XS*sumWeight

    !===== 3: generate an event (Q2=0):

    OutPart(1:3)%ID=(/nucleon, pion, pion /)
    OutPart(1)%mass=mN
    OutPart(2:3)%mass=mPi

    OutPart(1:3)%antiparticle=.false.
    OutPart(1:3)%scaleCS=1.

    p3 = momenta_in_3BodyPS (Wfree, OutPart(1:3)%mass)
    !   momenta are returned in 3-body cm frame

    do i=1,3
       OutPart(i)%momentum(1:3)=p3(:,i)
       OutPart(i)%momentum(0)=FreeEnergy(OutPart(i))
    end do

    !===== 4: boost the event to the gamma* N system:

    mediumAtPosition=mediumAt(eN%nucleon%position)

    ptot = eN%nucleon%momentum+eN%boson%momentum

    betaCMToLab = ptot(1:3)/ptot(0)
!   energy correction brings final state particles back to calculational frame
    call energyCorrection(W, (/0.,0.,0./), betaCMToLab, mediumAtPosition,OutPart,flagOK)
!   the following line adds a new ID to OutPart
!   that flags problems with energy correction

    if (.not.flagOK) then
       OutPart(1:3)%ID = 0
       return
    end if

    !===== 6: Take care of offshellness
    call setOffShellParameter(OutPart,flagOK)
    if (.not.flagOK) XS = 0.0

    OutPart(1:3)%perWeight=XS


    !    write(*,*) ' ----> XS = ',XS


  end subroutine DoNu2piBack


  !****************************************************************************
  !****s* neutrino2piBack/BGDualXS
  ! NAME
  ! subroutine BGDualXS(eN,sig2PiBG,indBG)
  !
  ! PURPOSE
  ! This subroutine produces the duality cross section in nanobarn
  ! in the mass range between 2 GeV and 3 GeV, based on parametrizations
  ! of structure function F2
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- the Lepton-Nucleon Event
  ! * integer :: indBG  -- index to choose Bloom-Gilman like parametriation for background
  !
  ! OUTPUT
  ! * real :: sig2piBG -- cross section d2sigma/(dE dOmega) in mb/GeV
  !****************************************************************************
  subroutine BGDualXS(eN,sig2piBG,indBG)

    use particleDefinition
    use leptonicID
    use constants, only: alphaQED,mN,GF,coscab,pi,hbarc,g_A
    use FF_QE_nucleonScattering, only: MA_in,MV2,mup,mun
!	use lepton2p2h, only: ME_ODW
    use neutrinoparms, only: ME_ODW
    

    type(electronNucleon_event), intent(in)  :: eN

    integer, intent(in)  :: indBG    ! parameter for choosing BG parametrization
    real, intent(out) :: sig2piBG

    real :: Ein,Eout,sinthetahalfsq,costhetahalfsq,nu,Q2,F1,F2,x,Wrec,qvec2
    real :: couplprop,GA0,GMV0,GA,GMV,GA2,GMV2,MA,F1nu,F2nu,F3nu,qnuk,RF2nF2p
    integer :: nuswitch,ID
	real :: GMp,GMn,GMp2,GMn2,GM2
	real :: kinfact

	sig2piBG = 0.0
	
    Ein = eN%lepton_in%momentum(0)
    Eout = eN%lepton_out%momentum(0)
    nu = Ein - Eout
    Q2 = eN%Qsquared
    Wrec = sqrt(mN**2 + 2*mN*nu - Q2)
    x = Q2/(2.*mN*nu)
    qvec2 = Q2 + nu**2
	qnuk = eN%nucleon_free%charge
	
	 ! now get structure functions F1 and F2
    ! from Bloom-Gilman Duality parametrizations
    call BGstruct(Q2,nu,F1,F2,indBG)  

!write(*,*) 'Q2=',Q2,'nu=',nu,'F1=',F1,'F2=',F2
	
	!   vector and axial coupling formfactors
    ! vector and axial coupling constants and cutoff masses
    GMV0 = mup - mun           ! vector component of magnetic coupling GM
    GA0 = - g_A
    MA = MA_in
    GMV = GMV0 * 1./(1 + Q2/MV2)**2
    GA = GA0 * 1./(1 + Q2/MA**2)**2
    GMV2 = GMV**2
    GA2 = GA**2
	
	
	if(qnuk == 1) then ! for protons:
		GMp = (1.0 + 1.793)/(1 + Q2/MV2)**2
		GMp2 = GMp**2
		GM2 = GMp2
	else !for neutrons:
		GMn = -1.913/(1 + Q2/MV2)**2
		GMn2 = GMn**2
		GM2 = GMn2
		! Now get structure functions for neutron (from Baillie et al, PRL 108, 2012)
		RF2nF2p = 0.89 - 0.65*x     ! ratio of neutron structure function to proton	
		F1 = F1 * RF2nF2p           ! structure functions for neutrons
        F2 = F2 * RF2nF2p
	end if

    select case (ME_ODW)
      case (1)       ! Martini et al
         kinfact = nu**2/qvec2
      case (2)        ! O'Connellet al.
         kinfact = qvec2/(2*mN)**2
      case (3)        ! O'Connellet al., relativistic
         kinfact = qvec2/(mN + sqrt(mN**2 + qvec2))**2
      case default
         write(*,*) 'kinfact error'
      end select
    
    ! Now get structure functions for neutrons from those for electrons
    F1nu = ((GMV2 + 1./kinfact * GA2)/GMV2)* F1
    F2nu = ((GMV2 + 1./kinfact * GA2)/GMV2)* F2
    F3nu =  4.0 * 1./kinfact * GA/sqrt(GMV2)* F1
	
!	write(*,*) 'qnuk=',qnuk
!	write(*,*) 'F1=',F1,'F2=',F2,'F1nu=',F1nu,'F2nu=',F2nu,'F3nu=',F3nu

	! now calculate cross sections
    ID = abs(eN%idProcess)     ! 1= electron, 2 = CC, 3 = NC, negative for antiparticles
    sinthetahalfsq = Q2/(4*Ein * Eout)
    costhetahalfsq = 1 - sinthetahalfsq
	
    select case (ID)

    case (1)      ! Electrons
       couplProp = (4*alphaQED**2)/Q2**2 * Eout**2
       sig2piBG = couplprop * (1./nu*F2*costhetahalfsq + 2./mN*F1*sinthetahalfsq)
       ! cross section in millibarn/GeV:
       sig2piBG = 10 * hbarc**2 * sig2piBG

    case (2)      ! Neutrinos CC
       couplProp = (GF*coscab)**2/(2*pi**2) * Eout**2
       nuswitch = sign(1,eN%IdProcess)

       sig2piBG = couplprop * (costhetahalfsq * 1./nu * F2nu &
            &  + sinthetahalfsq * 2 * 1./mN * F1nu  &
            &  - nuswitch*(Ein + Eout)/mN  &
            &  * F3nu/mN * sinthetahalfsq)
       ! T0 = target isospin
       sig2piBG = sig2piBG * 2. * (T0 + 1)

       if (sig2piBG < 0)   then
          write(*,*) 'enforce positivity constraint in 2pi BG, BG = 0'
          sig2piBG = 0.
       end if
       !X-section in mb/GeV
       sig2piBG = 10 * hbarc**2 * sig2piBG

    case (3)         ! Neutrinos NC
       couplProp = GF**2/(2*pi**2) * Eout**2 *1./2.
       ! Coupling for NC according to (C.12) in Leitner thesis
       nuswitch = sign(1,eN%IdProcess)

       sig2piBG = couplprop * (costhetahalfsq * 1./nu * F2nu &
            &  + sinthetahalfsq * 2 * 1./mN * F1nu  &
            &  - nuswitch*(Ein + Eout)/mN  &
            &  * F3nu/mN * sinthetahalfsq)
       ! T0 target isospin
       sig2piBG = sig2piBG * 2. * (T0 + 1)
	   if (sig2piBG < 0)   then
          write(*,*) 'enforce positivity constraint in 2pi BG, BG = 0'
!write(*,*) 'Q2=',Q2,'nu=',nu,'costhalf2=',costhetahalfsq,'sinthalf2=',sinthetahalfsq		  
!write(*,*) 'Wrec=',Wrec,'F1nu=',F1nu,'F2nu=',F2nu,'F3nu=',F3nu
          sig2piBG = 0.
       end if
       !X-section in mb/GeV
       sig2piBG = 10 * hbarc**2 * sig2piBG
	   
    case default
       write(*,*) 'reaction type must be EM, CC or NC'
       stop

    end select
			  

  end subroutine BGDualXS



  !****************************************************************************
  !****s* neutrino2piBack/BGStruct
  ! NAME
  ! subroutine BGStruct(Q2,nu,F1,F2,indBG)
  !
  ! PURPOSE
  ! This subroutine provides the structure functions F1 and F2
  ! obtained from various parametrizations:
  ! * 1. Bloom, Gilman, PRD 4, 2901 (1971)
  ! * 2. Niculescu et al, PRL 85, 1186 (2000)
  ! * 3. Christy, Bosted, PRC 81, 055213 (2010) (non-resonant)
  !
  ! INPUTS
  ! * real, intent(in)   :: Q2 = Q^2 of the incoming boson
  ! * real, intent(in)   :: nu = energy transfer
  ! * integer, intent(in):: indBG = index to choose structure function
  !   parametrizations according to numbers above
  ! OUTPUT
  ! * real, intent(out)  :: F1, F2 structure functions
  !****************************************************************************
  subroutine BGStruct(Q2,nu,F1,F2,indBG)

    use constants, only: mN,mPi,alphaQED,Pi,hbarc
    use ParamEP, only: ParamEP_BostedNonRes,ParamEP_R_BY

    integer, intent(in) :: indBG
    real, intent(in) :: Q2,nu
    real, intent(out) :: F1,F2
    real :: xi,xi1,nupr,Wrec,sigmaT,sigmaL,fak,RLT,K
    integer :: i


    F1 = 0.
    F2 = 0.
    Wrec = sqrt(mN**2 + 2*mN*nu - Q2)

    select case (indBG)

    case (1)         ! === original BG parametrization PR D4, 2901 (1971)
       if (Q2 == 0) return
       if(Wrec < 2.0) return ! limit of Bloom-Gilman fit

       nupr = nu + mN**2/Q2
       F2 = 0.557*(1. - 1./nupr)**3 + 2.1978*(1. - 1./nupr)**4 &
            - 2.5954*(1. - 1./nupr)**5
       F2 = F2 * 0.19732**2 * 1.e7   ! factor to get X-sections in nanobarn

       RLT = ParamEP_R_BY(Wrec,Q2)
       F1 = F2/(1. + RLT) * (Q2 + nu**2)/Q2 * mN/nu

    case (2)        ! === param. from Niculescu et al, PRL 85, 1186 (2000)
       if(Q2 == 0) return
       if(Wrec > 2.0 .or. Q2 < 0.3 .or. Q2 > 8.0) return ! limits of Niculescu fits

       xi = 2*Q2/(2*mN*nu)/(1. + sqrt(1. + 4*mN**2 * (Q2/(2*mN*nu))**2 /Q2**2))
       xi1 = 1 - xi
       F2 = xi**0.87*xi1**0.006*(0.005 - 0.058*xi1 - 0.017*xi1**2 &
            + 2.469*xi1**3  - 0.24*xi1**4)
       F2 = F2 * 0.19732**2 * 1.e7   ! factor to get X-sections in nanobarn

       RLT = ParamEP_R_BY(Wrec,Q2)
       F1 = F2/(1. + RLT) * (Q2 + nu**2)/Q2 * mN/nu

    case (3)        ! === param. of Christy and Bosted PR C81, 055213-1 (2010)

       call ParamEP_BostedNonRes(Wrec,Q2, sigmaL,sigmaT)

       K = (Wrec**2 - mN**2)/(2*mN)
       fak = K * 1/(4.*Pi**2 * alphaQED) * 1.e-4/hbarc**2
       F1 = fak * mN * sigmaT
       F2 = fak * Q2/(Q2+nu**2) * nu * (sigmaL + sigmaT)
       ! numerical factor 1.e-4 /hbarc**2 to make F1, F2 dimensionless
    end select

  end subroutine BGStruct



  !****************************************************************************
  !****if* neutrino2piBack/Att
  ! NAME
  ! real function Att(W,x)
  !
  ! PURPOSE
  ! Attenuation factor that multiplies BG-like background to give cross section,
  ! obtained from fitting sigma_incl(exp) - sigma_incl(1pi_theor)
  !****************************************************************************
  real function Att(W,x)

    use constants, only: mN
    use distributions, only: sigmoid
	use Neutrinoparms

    real, intent(in) :: W,x
    real :: dW,dx,Q2

    !  This fit is obtained for the Christy-Bosted non-resonant background
    !  fitted is difference between data - N* - 1pi
    !  in the mass region 1.3 < W < 2.0 GeV

    Q2 = x*(mN**2 - W**2)/(x-1.)

    dW = W - W00
    dx = x - x00

    Att = 1. + sigmoid(Q2,Q2cut,-0.1) * (A00 + A10*dW + A01*dx + A20*dW**2 &
         + A11*dW*dx + A02*dx**2 + A30*dW**3 + A21*dW**2*dx + A12*dW*dx**2 &
         + A03*dx**3 )

    if (Att > 1.) Att = 1.0
    if (Att < 0.) Att = 0.0

  end function Att

end module neutrino2piBack
