!******************************************************************************
!****m*  /XsectionRatios
! NAME
! module XsectionRatios
! NOTES
! Computes and stores the ratios of the in-medium
! and vacuum BB -> BB + mesons cross sections
! (see M. Wagner et al., PRC 71, 034910 (2005).
!******************************************************************************
module XsectionRatios

  implicit none
  private

  public :: accept_event, getShift0, getSigmaTotal, getSigmaScreened

  ! Global variables:

  integer, parameter :: nMshift = 20         ! number of the mass shift bins
  real, parameter :: MshiftBin = 0.04        ! mass shift bin (GeV)
  integer, parameter :: nSrtsStar = 100      ! number of the SrtsStar bins
  real, parameter :: SrtsStarBin = 0.1       ! SrtsStar bin (GeV)

  real, save, allocatable, dimension(:,:) :: SigmaTotal ! total in-medium pp cross section (mbarn)


  !****************************************************************************
  !****g* XsectionRatios/flagScreen
  ! SOURCE
  !
  logical, save :: flagScreen = .false.
  ! PURPOSE
  ! * If .true. -- in-medium screening is applied to the input cross section.
  ! * If .false. -- no cross section modification.
  !****************************************************************************


  !****************************************************************************
  !****g* XsectionRatios/ScreenMode
  ! SOURCE
  !
  integer, save :: ScreenMode = 1
  ! PURPOSE
  ! possible values:
  ! * 1: in-medium screening of NN total cross section according to
  !   Li and Machleidt, PRC 48 (1993) 1702 and PRC 49 (1994) 566
  ! * 2: in-medium screening of BB total cross section according to
  !   P. Daniewlewicz, NPA 673, 375 (2000); Acta. Phys. Pol. B 33, 45 (2002)
  ! NOTES
  ! relevant when  flagScreen = .true.
  !****************************************************************************


  !****************************************************************************
  !****g* XsectionRatios/flagInMedium
  ! SOURCE
  !
  logical, save :: flagInMedium = .false.
  ! PURPOSE
  ! * If .true. -- In-medium ratios are used to decide whether an event is
  !   accepted or not.
  ! * If .false. -- The event is always accepted
  !****************************************************************************


  !****************************************************************************
  !****g* XsectionRatios/InMediumMode
  ! SOURCE
  !
  integer, save :: InMediumMode = 2
  !
  ! PURPOSE
  ! possible values:
  ! * 1: all events of the type BB -> BB (+ mesons) are subject to
  !   in-medium reduction following Eqs.(194),(195) of GiBUU review paper
  !   [currently works in RMF mode only]
  ! * 2: BB -> BB events (except NN -> NN elastic scattering) are subject to
  !   in-medium reduction according to Eq. (33) from
  !   T. Song, C.M. Ko, PRC 91, 014901 (2015)
  !   [works in all modes (Skyrme, RMF, cascade)]
  ! NOTES
  ! relevant when flaginMedium = .true.
  !****************************************************************************


  !****************************************************************************
  !****g* XsectionRatios/alpha
  ! SOURCE
  !
  real, save :: alpha = 1.2
  ! PURPOSE
  ! Parameter which controls the density dependence of the NN <-> N Delta
  ! cross section.
  !
  ! for the density dependence from: Song/Ko, arXiv:1403.7363 (InMediumMode=2)
  !****************************************************************************

  !****************************************************************************
  !****g* XsectionRatios/shift0
  ! SOURCE
  !
  real, save :: shift0 = 0.
  ! PURPOSE
  ! Mass shift m-m^* (GeV) for using in elementary particle collision mode.
  !****************************************************************************


  logical, save :: iniFlag = .true.

contains

  !****************************************************************************
  !****f* XsectionRatios/getSigmaScreened
  ! NAME
  ! real function getSigmaScreened(pair,sigma,srts,mediumAtColl)
  !
  ! PURPOSE
  ! Returns in-medium screened/reduced cross section (mbarn).
  ! See:
  ! * Li and Machleidt, PRC 48 (1993) 1702 and PRC 49 (1994) 566,
  ! * P. Daniewlewicz, NPA 673, 375 (2000); Acta. Phys. Pol. B 33, 45 (2002)
  !
  ! INPUTS
  ! * type(particle), dimension(1:2) :: pair     -- incoming particles
  ! * real, intent(in) :: sigma                  -- not screened cross section
  ! * real, intent(in) :: srts                   -- c.m. energy
  ! * type(medium), intent(in) :: mediumAtColl   -- medium info
  !
  ! NOTES
  ! srts and mediumAtColl only needed for ScreenMode=1
  !****************************************************************************
  real function getSigmaScreened(pair,sigma,srts,mediumAtColl)

    use particleDefinition
    use IdTable, only: isMeson,nucleon
    use densitymodule, only: densityAt
    use dichtedefinition
    use mediumDefinition
    use lorentzTrafo
    use constants, only: mN

    type(particle), dimension(1:2), intent(in) :: pair
    real, intent(in) :: sigma
    real, intent(in) :: srts
    type(medium), intent(in) :: mediumAtColl

    real, parameter :: y = 1.2    ! coefficient at rho**(-2/3) in the maximal cross section
    real, dimension(1:3)  :: position, beta
    real, dimension(0:3)  :: j_baryon, momentumStar
    type(dichte) :: density
    real :: Elab,densfact,sigma_0

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    getSigmaScreened = sigma

    if ( .not.flagScreen ) return

    ! Currently only baryon-baryon collisions are included:
    if (isMeson(pair(1)%ID) .or. isMeson(pair(2)%ID)) return

    ! Exclude antibaryons:
    if ( pair(1)%antiparticle .or. pair(2)%antiparticle ) return

    select case(ScreenMode)

    case(1)  ! Li and Machleidt

         if(pair(1)%ID==nucleon .and. pair(2)%ID==nucleon) then
            Elab =((srts**2 - 2*mN**2)/(2*mN) - mN)*1000 ! kinetic lab energy in MeV

            if (pair(1)%charge == pair(2)%charge) then ! nn or pp
               densfact = (1. + 0.1667*Elab**1.05 * mediumATColl%density**3  &
                    /(1+exp((Elab - 350)/20.)))/          &
                    (1. + 9.704*mediumATColl%density**1.2/(1+exp((Elab - 350)/20.)))
            else  ! pn
               densfact = (1. + 0.0034*Elab**1.51 * mediumATColl%density**2  &
                    /(1+exp((Elab - 350)/20.)))/                       &
                    (1. + 21.55*mediumATColl%density**1.34/(1+exp((Elab - 350)/20.)))
            end if

            getSigmaScreened = sigma * densfact

         end if

    case(2)  ! Daniewlewicz

       position = (pair(1)%position+pair(2)%position)/2.
       density = densityAt(position)
       j_baryon = density%baryon
       momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
       beta = lorentzCalcBeta (momentumStar, 'XsectionRatios/getSigmaScreened')
       call lorentz(beta, j_baryon,'XsectionRatios/getSigmaScreened')

       if ( j_baryon(0) > 1.e-03 ) then
          sigma_0 = 10. * y / j_baryon(0)**0.666667
          getSigmaScreened = sigma_0 * tanh( sigma / sigma_0 )
       else
          getSigmaScreened = sigma
       end if

    end select


  end function getSigmaScreened


  !****************************************************************************
  !****f* XsectionRatios/getSigmaTotal
  ! NAME
  ! real function getSigmaTotal(pair)
  !
  ! PURPOSE
  ! Returns the total in-medium cross section (mbarn).
  !
  ! INPUTS
  ! * type(particle), dimension(1:2) :: pair        -- incoming particles
  !
  ! NOTES
  ! The total in-medium pp cross section is used for all types of
  ! colliding baryons/antibaryons.
  !
  ! Returns zero if there is at least one incoming meson.
  !****************************************************************************
  real function getSigmaTotal(pair)

    use particleDefinition
    use IdTable, only: isMeson

    type(particle), dimension(1:2), intent(in) :: pair

    real :: sqrtsStar, threshold, Q, shift
    real, dimension(0:3) :: momentumStar
    real, dimension(1:2) :: mstar
    integer :: i, iQ, ishift

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    if (isMeson(pair(1)%ID) .or. isMeson(pair(2)%ID)) then
       ! Use asympotic high energy value of the pi^+ p total cross section:
       getSigmaTotal = 20.
       return
    end if

    if ( .not.flagInMedium ) then
       ! Use asympotic high energy value of the pp total cross section:
       getSigmaTotal = 40.
       return
    end if

    ! Determine the sqrtsStar bin number:

    momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
    sqrtsStar = momentumStar(0)**2 - Dot_Product(momentumStar(1:3),momentumStar(1:3))
    sqrtsStar = sqrt( max(0.,sqrtsStar) )
    do i=1,2
       mstar(i) = pair(i)%momentum(0)**2 - dot_product( pair(i)%momentum(1:3),&
            & pair(i)%momentum(1:3) )
       mstar(i) = sqrt( max(0.,mstar(i)) )
    end do
    threshold = mstar(1) + mstar(2)
    Q = sqrtsStar - threshold
    iQ = nint(Q/SrtsStarBin)
    iQ = max(1,min(iQ,nSrtsStar))

    ! Determine the mass shift bin number:

    shift = pair(1)%mass - mstar(1)
    ishift = nint(shift/MshiftBin)
    ishift = max(0,min(ishift,nMshift))

    ! Use total in-medium pp cross section:

    getSigmaTotal = SigmaTotal(ishift,iQ)

  end function getSigmaTotal

  !****************************************************************************
  !****************************************************************************
  real function getShift0()
    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if
    getShift0 = shift0
  end function getShift0


  !****************************************************************************
  !****f* XsectionRatios/accept_event
  ! NAME
  ! logical function accept_event(pair,finalState)
  !
  ! PURPOSE
  ! Accepts or rejects the event randomly, using the in-medium cross section
  ! ratios:
  ! * returns .true.  --- if the event is accepted
  ! * returns .false. --- if the event is rejected
  !
  ! INPUTS
  ! * type(particle), dimension(1:2) :: pair        -- incoming particles
  ! * type(particle), dimension(:)   :: finalState  -- outgoing particles
  !
  ! NOTES
  ! the procedure is as follows:
  ! first, a particle collision is chosen with the default cross sections.
  ! Then, in a second step, the ratio of the in-medium and the vacuum cross
  ! section is calculated. If smaller than 1, then the event generated in
  ! step 1 is rejected by a Monte Carlo decision.
  !****************************************************************************
  logical function accept_event(pair,finalState)

    use IdTable
    use particleDefinition
    use particleProperties, only: hadron
    use random
    use inputGeneral, only: eventtype
    use eventtypes, only: elementary
    use dichtedefinition
    use densitymodule, only: densityAt
    use minkowski, only: abs4
    use constants, only: rhoNull

    type(particle), dimension(1:2), intent(in) :: pair
    type(particle), dimension(:), intent(in) :: finalState

    type(dichte) :: density

    real :: sqrtsStar, factor , densityNuc
    real, dimension(0:3) :: momentumStar
    real, dimension(1:2) :: mstar
    real, dimension(1:3)  :: position

    integer :: nFinal, nBaryon, i

    integer, dimension(1:2) :: idBaryon    ! Id's of outgoing baryons
    integer, dimension(1:2) :: iBaryon     ! position numbers of outgoing baryons in the finalState

    integer, parameter :: nmax=100
    real, dimension(1:nmax) :: mstar_final

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    if ( .not.flagInMedium ) then
       accept_event = .true.
       return
    end if

    if ( isMeson(pair(1)%Id) .or. isMeson(pair(2)%Id) ) then

       ! No medium modifications implemented for meson-baryon and meson-meson collisions:
       accept_event = .true.
       return

    else if ( pair(1)%antiparticle .and. pair(2)%antiparticle ) then

       ! No medium modifications implemented for antibaryon-antibaryon collisions:
       accept_event = .true.
       return

    else if (      hadron(pair(1)%Id)%strangeness.ne.0 &
         &  .or. hadron(pair(1)%Id)%charm.ne.0 &
         &  .or. hadron(pair(2)%Id)%strangeness.ne.0 &
         &  .or. hadron(pair(2)%Id)%charm.ne.0  ) then

       !  No medium modifications implemented for collisions involving
       !  a strange or a charmed particle:
       accept_event = .true.
       return

    end if

    ! Incoming baryons are nonstrange and noncharmed

    ! Analyse the outgoing particles:

    nFinal = 0
    nBaryon = 0
    idBaryon = 0
    do i=1,size(finalState,dim=1)

       if (finalState(i)%id <= 0) exit

       nFinal = nFinal + 1

       if(nFinal.gt.nmax) then
          accept_event = .true.
          return
       end if

       if ( isBaryon(finalState(i)%id) ) then

          if ( finalState(i)%antiparticle ) then
             !  No medium modifications implemented for baryon-antibaryon pair production:
             accept_event = .true.
             return
          end if

          nBaryon = nBaryon + 1
          if ( nBaryon <= 2 ) then
             idBaryon(nBaryon) = finalState(i)%Id
             iBaryon(nBaryon) = i
          else
             !  No medium modifications implemented for more than 2 outgoing baryons
             !  (this also must be baryon-antibaryon pair production):
             accept_event = .true.
             return
          end if

       end if

    end do


    if ((pair(1)%antiparticle.neqv.pair(2)%antiparticle) .and. nBaryon.eq.0) then

       ! Baryon-antibaryon annihilation:

       if ( rn() .lt.  Ratio_BaB(pair,finalState) ) then
          accept_event = .true.
       else
          accept_event = .false.
       end if

       return

    end if


    if (        hadron(idBaryon(1))%charm .ne. 0 &
         & .or.   hadron(idBaryon(2))%charm .ne. 0 &
         & .or.   hadron(idBaryon(1))%strangeness &
         & + hadron(idBaryon(2))%strangeness.le.-2  ) then
         ! No medium modifications for outgoing charmed baryons or
         ! outgoing pair of strange baryons, or a baryon with S=-2
         accept_event = .true.
         return
    end if

    ! BB -> BB (+mesons) event:

    select case(InMediumMode)

       case(1) ! modify all events of the type BB -> BB (+ mesons)

           momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
           sqrtsStar = momentumStar(0)**2 - Dot_Product(momentumStar(1:3),momentumStar(1:3))
           sqrtsStar = sqrt( max(0.,sqrtsStar) )

           if ( .not. eventtype==elementary ) then

              do i = 1,2
                 mstar(i) = pair(i)%momentum(0)**2 - dot_product( pair(i)%momentum(1:3),&
                                                                & pair(i)%momentum(1:3) )
                 mstar(i) = sqrt( max(0.,mstar(i)) )
              end do

              do i = 1,nFinal
                 mstar_final(i)=finalState(i)%momentum(0)**2 &
                              & - dot_product( finalState(i)%momentum(1:3),&
                              &                finalState(i)%momentum(1:3) )
                 mstar_final(i) = sqrt( max(0.,mstar_final(i)) )
              end do

           else

              mstar(1:2) = pair(1:2)%mass - shift0

              do i = 1,nFinal
                 if ( isBaryon(finalState(i)%id) ) then
                    mstar_final(i) = finalState(i)%mass - shift0
                 else
                    mstar_final(i) = finalState(i)%mass
                 end if
              end do

              sqrtsStar = sqrtsStar - 2.*shift0

           end if


           if (2.le.nFinal .and. nFinal.le.6) then
               factor=Ratio(sqrtsStar,(/mstar(1:2),mstar_final(1:nFinal)/),&
                           &(/pair(1:2)%mass,finalState(1:nFinal)%mass/))
           else
               factor=Ratio(sqrtsStar,(/mstar(1:2),mstar_final(iBaryon(1)),mstar_final(iBaryon(2))/),&
                          &(/pair(1:2)%mass,finalState(iBaryon(1))%mass,finalState(iBaryon(2))%mass/))
           end if

       case(2)  ! modify only BB -> BB events (except NN -> NN elastic scattering)

           if(pair(1)%Id+pair(2)%Id+idBaryon(1)+idBaryon(2).ge.5 .and.  nFinal.eq.2) then
               position=(pair(1)%position+pair(2)%position)/2.
               density=densityAt(position)
               densityNuc=abs4(density%proton)+abs4(density%neutron)
               factor=exp(-alpha*densityNuc/rhoNull)
           else
               accept_event = .true.
               return
           end if

    end select

    if ( factor > 1 ) then
       write(*,*) ' Warning: in accept_event factor=', factor
       write(*,*) pair(:)%Id, finalState(1:nFinal)%Id
    end if

    if ( rn() .le. factor ) then

       accept_event = .true.

    else

       accept_event = .false.

    end if

  end function accept_event


  !****************************************************************************
  !****s* XsectionRatios/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in input switches from namelist XsectionRatios_input.
  !****************************************************************************
  subroutine init

    use output
    use inputGeneral, only: path_To_Input

    integer :: ios, i, j, k
    real :: shift, Q
    character(8) :: dummy1, dummy2

    !**************************************************************************
    !****n* XsectionRatios/XsectionRatios_input
    ! NAME
    ! NAMELIST /XsectionRatios_input/
    ! PURPOSE
    ! Includes the switches:
    ! * flagScreen
    ! * ScreenMode
    ! * flagInMedium
    ! * InMediumMode
    ! * alpha
    ! * shift0
    !**************************************************************************
    NAMELIST /XsectionRatios_input/ flagScreen, ScreenMode, flagInMedium, &
         InMediumMode, alpha, shift0

    call Write_ReadingInput('XsectionRatios_input',0)
    rewind(5)
    read(5,nml=XsectionRatios_input,iostat=ios)
    call Write_ReadingInput('XsectionRatios_input',0,ios)
    write(*,*) ' Set flagScreen to', flagScreen , '.'
    write(*,*) ' Set ScreenMode to', ScreenMode , '.'
    write(*,*) ' Set flagInMedium to', flagInMedium ,'.'
    write(*,*) ' Set InMediumMode to', InMediumMode , '.'
    write(*,*) ' Set alpha to', alpha , '.'
    write(*,*) ' Set shift0 to', shift0 ,'.'
    call Write_ReadingInput('XsectionRatios_input',1)

    if ( flagScreen .and. flagInMedium ) then
       write(*,*) ' In-medium screening and in-medium event acceptance '
       write(*,*) ' can not be done simultaneously !!!'
       stop
    end if

    if (flagInMedium) then

       allocate( SigmaTotal(0:nMshift,1:nSrtsStar) )

       ! Read-in the total in-medium pp cross section:
       write(*,*) ' Total in-medium pp cross section read-in starts ...'
       ios=0
       open(1,file=trim(path_to_input)//'/XsectionTotal.dat',status='old',iostat=ios)
       do i = 0,nMshift
          read(1,*) dummy1, dummy2, shift
          do j = 1,nSrtsStar
             read(1,*) Q, SigmaTotal(i,j)
          end do
       end do
       close(1)
       write(*,*) ' Total in-medium pp cross section read-in is successfully finished'

!                open(1,file='XsectionTotal_chk.dat',status='unknown')
!                do i = 0,nMshift
!                  shift = MshiftBin * float(i)
!                  write(1,*)'# shift:', shift
!                  do j = 1,nSrtsStar
!                    Q = SrtsStarBin * float(j)
!                    write(1,'(2(e13.6,1x))') Q, SigmaTotal(i,j)
!                  end do
!                end do
!                close(1)

    end if

  end subroutine init



  !****************************************************************************
  !****f* XsectionRatios/Ratio
  ! NAME
  ! real function Ratio(srtsStar,mStar,mass)
  !
  ! PURPOSE
  ! Compute the ratio of the in-medium and vacuum
  ! B_1 + B_2 -> B_3 + B_4 + M_5 + M_6 + M_7 + M_8  cross sections
  !
  ! INPUTS
  ! * real :: srtsStar           -- c.m. energy (GeV),
  ! * real dimension(:) :: mStar -- Dirac masses of all particles (GeV)
  ! * real dimension(:) :: mass  -- vacuum masses of all particles (GeV)
  !
  ! OUTPUT
  ! * real :: Ratio  -- ratio of the in-medium and vacuum cross sections
  !
  ! NOTES
  ! The particles 1 and 2 are incoming baryons, the particles 3 and 4 are
  ! outgoing baryons, the particles M_5, ..., M_8 are outgoing mesons.
  ! The routine will compute one of the ratios:
  ! * B_1 + B_2 -> B_3 + B_4,
  ! * B_1 + B_2 -> B_3 + B_4 + M_5,
  ! * ...,
  ! * B_1 + B_2 -> B_3 + B_4 + M_5 + ...+ M_8
  ! depending on the the size n of the arrays mStar and mass, which must be in
  ! the interval [4:8].
  !****************************************************************************
  real function Ratio(srtsStar,mStar,mass)

    use nBodyPhaseSpace
    use twoBodyTools, only: pCM

    real, intent(in) :: srtsStar
    real, intent(in), dimension(:) :: mStar, mass

    integer :: i,n
    real :: srts, Ifac, IfacStar, phaseSpace, phaseSpaceStar

    n = size(mStar,dim=1)

    if ( n < 4 .or. n > 8 ) then
       write(*,*) ' In XsectionRatios/Ratio: wrong size of input array, n= ', n
       stop
    end if

    if ( n .ne. size(mass,dim=1) ) then
       write(*,*) ' In XsectionRatios/Ratio: different sizes of input arrays ', &
            & n, size(mass,dim=1)
       stop
    end if

    if ( srtsStar <= sum(mStar(1:2)) ) then
       write(*,*) ' In XsectionRatios/Ratio: srtsStar, mstar(1:2) : ', srtsStar, mstar(1:2)
       write(*,*) ' Sum of incoming Dirac masses', sum(mStar(1:2))
       stop
    end if

    if ( srtsStar <= sum(mStar(3:n)) + 0.001 ) then
       Ratio = 0.
       return
    end if

    ! c.m. energy in vacuum:
    srts = srtsStar - sum(mStar(3:n)) + sum(mass(3:n))

    ! Vacuum and in-medium flux factors:
    Ifac     = pcm(srts,mass(1),mass(2)) * srts
    IfacStar = pcm(srtsStar,mStar(1),mStar(2)) * srtsStar

    ! Vacuum and in-medium phase spaces:
    phaseSpace     = integrate_nBodyPS(srts,mass(3:n))
    phaseSpaceStar = integrate_nBodyPS(srtsStar,mStar(3:n))

    ! Cross section ratio = sigma_med / sigma_vac:
    if (IfacStar.gt.0. .and. phaseSpace.gt.0.) then
        Ratio = mStar(1)*mStar(2)/(mass(1)*mass(2)) &
              & *Ifac/IfacStar * phaseSpaceStar/phaseSpace
        do i = 3,n
           if (mass(i).gt.0.) Ratio = Ratio * mStar(i)/mass(i)
        end do
    end if

  end function Ratio


  !****************************************************************************
  !****f* XsectionRatios/Ratio_BaB
  ! NAME
  ! real function Ratio_BaB(pair,finalState)
  ! PURPOSE
  ! Compute the ratio of the in-medium and vacuum cross sections
  ! baryon + antibaryon -> mesons
  ! INPUTS
  ! * type(particle), dimension(1:2) :: pair        -- incoming particles
  ! * type(particle), dimension(:) :: finalState    -- outgoing particles
  ! OUTPUT
  ! * real :: Ratio_BaB  -- ratio of the in-medium and vacuum cross sections
  ! NOTES
  ! Maximum 6 mesons can be treated (otherwice returns Ratio_BaB=1.).
  !****************************************************************************
  real function Ratio_BaB(pair,finalState)

    use particleDefinition, only: particle, sqrtS
    use nBodyPhaseSpace
    use RMF, only: getRMF_flag
    use densitymodule, only: Particle4Momentum_RMF
    use twoBodyTools, only: sqrtS_Free,pCM

    type(particle), dimension(1:2), intent(in) :: pair
    type(particle), dimension(:), intent(in) :: finalState

    real :: srtS_star, srtS_vacuum, srtS    ! c.m. energies (GeV),
    real :: mstar1, mstar2       ! Dirac masses of incoming particles (GeV)
    real, dimension(0:3) :: momentum1, momentum2, momentum_tot

    integer :: n, i
    real :: Ifac_vacuum, Ifac, phaseSpace_vacuum, phaseSpace

    logical, parameter :: debugFlag = .false.

    ! Determine how many particles are in final state:
    n=0
    do i=1,size(finalState,dim=1)
       if (finalState(i)%Id.eq.0) then
          cycle
       else if (finalState(i)%Id.lt.0) then
          exit
       else
          n=n+1
       end if
    end do

    if ( n < 2 .or. n > 6 ) then
       Ratio_BaB=1.
       return
    end if

    srtS_star = sqrtS(pair,"Ratio_BaB, srtS_star")

    mstar1 = sqrtS(pair(1),'Ratio_BaB, mstar(1)')
    mstar2 = sqrtS(pair(2),'Ratio_BaB, mstar(2)')

    if ( .not.getRMF_flag() ) then
       srtS_vacuum=sqrtS_free(pair)
       srtS=srtS_star
    else
       srtS_vacuum  = srtS_star - mstar1 - mstar2 + pair(1)%mass + pair(2)%mass
       call Particle4Momentum_RMF(pair(1),momentum1)
       call Particle4Momentum_RMF(pair(2),momentum2)
       ! Compute srtS with canonical momenta:
       momentum_tot(0:3) = momentum1(0:3) + momentum2(0:3)
       srtS = momentum_tot(0)**2 - dot_product(momentum_tot(1:3),momentum_tot(1:3))
       srtS = sqrt(max(0.,srtS))
    end if

    if (debugFlag) then
       write(*,*) ' In Ratio_BaB: incoming particles: ', pair%Id,pair%antiparticle,pair%charge
       write(*,*) ' In Ratio_BaB: srtS_star, srtS, srtS_vacuum: ', srtS_star, srtS, srtS_vacuum
       write(*,*) ' In Ratio_BaB: outgoing particles: ', finalState(1:n)%Id
       write(*,*) ' In Ratio_BaB: sum of outgoing masses: ', sum(finalState(1:n)%mass)
    end if


    if ( srtS <= sum(finalState(1:n)%mass) + 0.001 ) then
       Ratio_BaB = 0.
       return
    end if

    ! Vacuum and in-medium flux factors:
    Ifac_vacuum = pcm(srtS_vacuum,pair(1)%mass,pair(2)%mass) * srtS_vacuum
    Ifac        = pcm(srtS_star,mstar1,mstar2) * srtS_star

    ! Vacuum and in-medium phase spaces:
    phaseSpace_vacuum = integrate_nBodyPS (srtS_vacuum, finalState(1:n)%mass)
    phaseSpace        = integrate_nBodyPS (srtS, finalState(1:n)%mass)

    ! Cross section ratio = sigma_med / sigma_vac:
    Ratio_BaB = mstar1*mstar2/(pair(1)%mass*pair(2)%mass) &
                * Ifac_vacuum/Ifac * phaseSpace/phaseSpace_vacuum

    if (debugFlag) then
       write(*,*) ' In Ratio_BaB: mstar1, mstar2: ', mstar1, mstar2
       write(*,*) ' In Ratio_BaB: Ifac_vacuum, Ifac: ',Ifac_vacuum, Ifac
       write(*,*) ' In Ratio_BaB: phaseSpace_vacuum, phaseSpace: ', phaseSpace_vacuum, phaseSpace
       write(*,*) ' In Ratio_BaB: Ratio_BaB', Ratio_BaB
    end if

  end function Ratio_BaB


end module XsectionRatios
