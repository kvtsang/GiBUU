!******************************************************************************
!****m* /EventOutput
! NAME
! module EventOutput
!
! PURPOSE
! This module provides classes for writing event output to disk
! in different formats.
!
! Currently the following formats are supported:
! * "Les Houches" event files
! * "OSCAR 2013" event files
! * "Shanghai 2014"
! * "Root"
!
! For a description of the Les Houches format, please refer to:
! * http://arxiv.org/abs/hep-ph/0609017
! * https://gibuu.hepforge.org/trac/wiki/LesHouches
!
! For a description of the OSCAR 2013 format, see:
! * http://phy.duke.edu/~jeb65/oscar2013
!
! INPUTS
! (none)
!******************************************************************************
module EventOutput

  implicit none
  private

  !****************************************************************************
  !****t* EventOutput/EventOutputFile
  ! NAME
  ! type EventOutputFile
  ! PURPOSE
  ! This is an abstract base type to represent a file for event output.
  ! It is used as a common interface for LHOutputFile and OscarOutputFile.
  !
  ! SOURCE
  !
  type, abstract, public :: EventOutputFile
     integer, private :: iFile = 0  ! private file handle
  contains
    ! deferred type-bound procedures that need to be implemented in the derived classes
    procedure(open_ifc),       deferred :: open
    procedure(close_ifc),      deferred :: close
    procedure(write_EH_ifc),   deferred :: write_event_header
    procedure(write_EF_ifc),   deferred :: write_event_footer
    procedure(write_part_ifc), deferred :: write_particle
    procedure :: write_additionalInfo
    ! type-bound procedures that are implemented in the base class
    procedure :: write_real
    procedure :: write_pert
  end type
  !****************************************************************************


  ! interfaces for the deferred methods
  abstract interface
    subroutine open_ifc(this, pert, nCall, nTimeStep)
      import :: EventOutputFile
      class(EventOutputFile) :: this
      logical, intent(in) :: pert
      integer, intent(in) :: nCall, nTimeStep
    end subroutine
    subroutine close_ifc(this)
      import :: EventOutputFile
      class(EventOutputFile), intent(in) :: this
    end subroutine
    subroutine write_EH_ifc(this, nParts, nEvent, wgt, iFE)
      import :: EventOutputFile
      class(EventOutputFile) :: this
      integer, intent(in) :: nParts
      integer, intent(in) :: nEvent            ! number of current event
      real, intent(in), optional :: wgt
      integer, intent(in), optional :: iFE
    end subroutine
    subroutine write_EF_ifc(this)
      import :: EventOutputFile
      class(EventOutputFile), intent(in) :: this
    end subroutine
    subroutine write_part_ifc(this, part)
      use particleDefinition
      import :: EventOutputFile
      class(EventOutputFile), intent(in) :: this
      type(particle), intent(in) :: part
    end subroutine
  end interface


  !****************************************************************************
  !****t* EventOutput/LHOutputFile
  ! NAME
  ! type LHOutputFile
  ! PURPOSE
  ! This is an extended type for event output in LesHouches format.
  ! It is derived from the base type EventOutputFile and implements its
  ! interfaces.
  !
  ! For a description of the Les Houches format, please refer to:
  ! * http://arxiv.org/abs/hep-ph/0609017
  ! * https://gibuu.hepforge.org/trac/wiki/LesHouches
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: LHOutputFile
  contains
    procedure :: open                 => LH_open
    procedure :: close                => LH_close
    procedure :: write_event_header   => LH_write_event_header
    procedure :: write_event_footer   => LH_write_event_footer
    procedure :: write_particle       => LH_write_particle
    procedure :: write_additionalInfo => LH_write_additionalInfo
  end type
  !****************************************************************************


  !****************************************************************************
  !****t* EventOutput/OscarOutputFile
  ! NAME
  ! type OscarOutputFile
  ! PURPOSE
  ! This is an extended type for event output in OSCAR 2013 format.
  ! It is derived from the base type EventOutputFile and implements its
  ! interfaces.
  !
  ! For a description of the OSCAR 2013 format, see:
  ! * http://phy.duke.edu/~jeb65/oscar2013
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: OscarOutputFile
     integer, private :: nEvent = 0   ! store the event number
     integer, private :: iFE = 0      ! store the iFE value
  contains
    procedure :: open                 => Oscar_open
    procedure :: close                => Oscar_close
    procedure :: write_event_header   => Oscar_write_event_header
    procedure :: write_event_footer   => Oscar_write_event_footer
    procedure :: write_particle       => Oscar_write_particle
  end type
  !****************************************************************************

  !****************************************************************************
  !****t* EventOutput/OscarExtOutputFile
  ! NAME
  ! type OscarExtOutputFile
  ! PURPOSE
  ! This is an extended type for event output in OSCAR 2013 format.
  ! It is derived from the base type EventOutputFile and implements its
  ! interfaces.
  !
  ! For a description of the OSCAR 2013 format, see:
  ! * http://phy.duke.edu/~jeb65/oscar2013
  !
  ! This is an extended version whith much more output columns
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: OscarExtOutputFile
     integer, private :: nEvent = 0   ! store the event number
     integer, private :: iFE = 0      ! store the iFE value
   contains
    procedure :: open                 => OscarExt_open
    procedure :: close                => OscarExt_close
    procedure :: write_event_header   => OscarExt_write_event_header
    procedure :: write_event_footer   => OscarExt_write_event_footer
    procedure :: write_particle       => OscarExt_write_particle
  end type
  !****************************************************************************


  !****************************************************************************
  !****t* EventOutput/ShanghaiOutputFile
  ! NAME
  ! type ShanghaiOutputFile
  ! PURPOSE
  ! This is an extended type for event output in Shanghai2014 format.
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: ShanghaiOutputFile
  contains
    procedure :: open                 => Shanghai_open
    procedure :: close                => Shanghai_close
    procedure :: write_event_header   => Shanghai_write_event_header
    procedure :: write_event_footer   => Shanghai_write_event_footer
    procedure :: write_particle       => Shanghai_write_particle
  end type
  !****************************************************************************

  !****************************************************************************
  !****t* EventOutput/RootOutputFile
  ! NAME
  ! type RootOutputFile
  ! PURPOSE
  ! This is an extended type for event output in Root format.
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: RootOutputFile
    real, private :: weight
  contains
    procedure :: open                 => Root_open
    procedure :: close                => Root_close
    procedure :: write_event_header   => Root_write_event_header
    procedure :: write_event_footer   => Root_write_event_footer
    procedure :: write_particle       => Root_write_particle
    procedure :: write_additionalInfo => Root_write_additionalInfo
  end type RootOutputFile
  !****************************************************************************


contains

!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/write_additionalInfo
  ! NAME
  ! subroutine write_additionalInfo(this, iFE, pNode)
  ! PURPOSE
  ! Write additional info about the event, depending on eventtype.
  !****************************************************************************
  subroutine write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition, only: tParticleListNode
    class(EventOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode
    ! do nothing here by default
  end subroutine write_additionalInfo

!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/LH_open
  ! NAME
  ! subroutine LH_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output according to the
  ! "Les Houches Event Files" standard.
  !****************************************************************************
  subroutine LH_open(this, pert, nCall, nTimeStep)
    class(LHOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (pert) then
      buf1 = '.Pert.'
      this%iFile = 721
    else
      buf1 = '.Real.'
      this%iFile = 722
    end if
    write(buf2,'(I8.8)') nCall
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.lhe'

    ! open file
    open(this%iFile, file=fName, status='unknown')
    rewind(this%iFile)

    write(this%iFile,'(A)') '<LesHouchesEvents version="1.0">'
    write(this%iFile,'(A)') '<!-- File generated by GiBUU. For documentation see ' //  &
                            'https://gibuu.hepforge.org/trac/wiki/LesHouches -->'

    write(this%iFile,'(A)') '<header>'
    write(this%iFile,'(A)') '     <!-- individual XML tags may follow -->'
    write(this%iFile,'(A)') '</header>'

    write(this%iFile,'(A)') '<init>'
    write(this%iFile,'(1P,2I8,2E14.6,6I6)') 0,0, 0.,0., 0,0,0,0,0,0
    write(this%iFile,'(A)') '</init>'

  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_close
  ! NAME
  ! subroutine LH_close(this)
  ! PURPOSE
  ! Close a file after outputting Les-Houches event information.
  !****************************************************************************
  subroutine LH_close(this)
    class(LHOutputFile), intent(in) :: this

    write(this%iFile, '(A)') '</LesHouchesEvents>'
    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_event_header
  ! NAME
  ! subroutine LH_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for a Les-Houches event, including the number of particles
  ! and the event weight.
  !****************************************************************************
  subroutine LH_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(LHOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(len=15), parameter :: f1 = '(1P,2I6,4E14.6)'
    integer, parameter :: IDPRUP = 0
    real, parameter :: SCALUP = 0.0, AQEDUP = 0.0, AQCDUP = 0.0
    real :: weight

    if (present(wgt)) then
      weight = wgt
    else
      weight = 1.
    end if

    write(this%iFile,'(A)') '<event>'
    write(this%iFile,f1) nParts, IDPRUP, weight, SCALUP, AQEDUP, AQCDUP
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_event_footer
  ! NAME
  ! subroutine LH_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a Les-Houches event.
  !****************************************************************************
  subroutine LH_write_event_footer(this)
    class(LHOutputFile), intent(in) :: this
    write(this%iFile,'(A)') '</event>'
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_particle
  ! NAME
  ! subroutine LH_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in Les-Houches format.
  !****************************************************************************
  subroutine LH_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU

    class(LHOutputFile), intent(in) :: this
    type(particle), intent(in) :: part

    character(len=22), parameter :: f2 = '(1P,I8,5I5,5E18.10,A6)'
    integer :: KF

    KF = KFfromBUU(part)
    write(this%iFile,f2) KF, 0, 0,0, 0,0, &
                         part%momentum(1:3), part%momentum(0), &
                         sqrts(part), '0. 9.'
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_additionalInfo
  ! NAME
  ! subroutine LH_write_additionalInfo(this, iFE)
  ! PURPOSE
  ! Write additional info about the event, depending on eventtype.
  !
  ! This routine tries to find additional information about the event.
  ! It tries routines for different event types, which only return
  ! some information, if it was really stored.
  !
  ! The following cases are handled:
  ! * For eventtype "HiLep", the following line is added:
  !     # 14 nu Q2 eps phiLepton Eventtype
  !   (14 is the magic number of "HiLepton")
  ! * For eventtype "neutrino", the following line is added:
  !     # 5 Eventtype Weight momLepIn(0:3) momLepOut(0:3) momNuc(0:3)
  !   (5 is the magic number for neutrino events)
  ! * For eventtype "heavyIon", the following line is added:
  !     # 1 b
  !   (1 is the magic number of "heavyIon", b is the impact parameter in fm)
  ! * For eventtype "hadron", the following line is added:
  !     # 300 b
  !   (300 is the magic number of "hadron", b is the impact parameter in fm)
  !****************************************************************************
  subroutine LH_write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use inputGeneral, only: eventType
    use eventtypes, only: hiLepton, neutrino, heavyIon, hadron
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use FreezeoutAnalysis, only: getFreezeoutAnalysis_Pert
    use PIL_freezeout, only: PIL_freezeout_GET

    class(LHOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode

    real :: weight,nu,Q2,eps,phiL
    integer :: evtType, chrg_nuc
    real,dimension(0:3) :: momLepIn, momLepOut, momBos, momNuc
    type(particle), pointer :: pPart
    real, dimension(0:3) :: pos
    integer :: history
    logical :: escaped

    select case (eventType)
    case (heavyIon)
      write(this%iFile,'(A,ES13.4)') '# 1 ', b_HI
    case (hadron)
      write(this%iFile,'(A,ES13.4)') '# 300 ', b_had
    case (neutrino)
      if (.not. present(iFE)) return
      if (NeutrinoProdInfo_Get(iFE,evtType,Weight,momLepIn,momLepOut,momBos,momNuc,chrg_nuc)) &
        write(this%iFile,'(A,I5,1P,e18.10,1P,3(" ",4e18.10),0P,A)') &
           '# 5 ', evtType, Weight, momLepIn, momLepOut, momNuc
    case (hiLepton)
      if (.not. present(iFE)) return
      if (EventInfo_HiLep_Get(0,iFE,Weight,nu,Q2,eps,evtType,phi_Lepton=phiL)) &
        write(this%iFile,'(A,1P,4e13.4,0P,I8)') '# 14 ', nu,Q2,eps,phiL,evtType
    end select

    if (getFreezeoutAnalysis_Pert() .and. present(pNode)) then
      do
        if (.not. associated(pNode)) exit
        pPart => pNode%V
        if (PIL_freezeout_GET(pPart%number, pos, history, escaped)) then
          write(this%iFile,'(A,1P,4e13.4,0P,I12,L3)') '# 1001 ', pos, history, escaped
        else
          write(this%iFile,'(A,1P,4e13.4,0P,I12,L3)') '# 1002 ', pos, history, escaped
        end if

        pNode => pNode%next
      end do
    end if

  end subroutine LH_write_additionalInfo


!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/Oscar_open
  ! NAME
  ! subroutine Oscar_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output according to the
  ! "OSCAR 2013" standard.
  !****************************************************************************
  subroutine Oscar_open(this, pert, nCall, nTimeStep)
    class(OscarOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=4)  :: buf

    if (pert) then
      buf = 'Pert'
      this%iFile = 723
    else
      buf = 'Real'
      this%iFile = 724
    end if
    fName = 'EventOutput.' // trim(buf) // '.oscar'

    if (nCall == 1) then
      ! open file for the first time
      open(this%iFile, file=fName, status='unknown')
      ! write header
      write(this%iFile,'(A)') '#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID'
      write(this%iFile,'(A)') '# Units: fm fm fm fm GeV GeV GeV GeV GeV none none'
      write(this%iFile,'(A)') '# File generated by GiBUU (https://gibuu.hepforge.org)'
    else
      ! append to exiting file
      open(this%iFile, file=fName, status='old', position='append')
    end if

  end subroutine


  !****************************************************************************
  !****s* EventOutput/Oscar_close
  ! NAME
  ! subroutine Oscar_close(this)
  ! PURPOSE
  ! Close a file after outputting OSCAR 2013 event information.
  !****************************************************************************
  subroutine Oscar_close(this)
    class(OscarOutputFile), intent(in) :: this

    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Oscar_write_event_header
  ! NAME
  ! subroutine Oscar_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an OSCAR 2013 event, including the number of particles
  ! and the event weight.
  !****************************************************************************
  subroutine Oscar_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(OscarOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(*), parameter :: f1 = '("# event ",I9," out",I4," weight",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," out",I4)'

    this%nEvent = nEvent
    this%iFE = -1
    if (present(iFE)) this%iFE = iFE

    if (present(wgt)) then
       write(this%iFile,f1) this%nEvent, nParts, wgt
    else
       write(this%iFile,f0) this%nEvent, nParts
    end if

  end subroutine

  !****************************************************************************
  !****s* EventOutput/Oscar_write_event_footer
  ! NAME
  ! subroutine Oscar_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes an OSCAR 2013 event.
  !****************************************************************************
  subroutine Oscar_write_event_footer(this)
    use inputGeneral, only: eventType
    use eventtypes, only: heavyIon, hadron, LoPion, HiPion
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use initPion, only: getImpact_Lo => getImpact
    use initHiPion, only: getImpact_Hi => getImpact

    class(OscarOutputFile), intent(in) :: this

    character(*), parameter :: f1 = '("# event ",I9," end 0 impact",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," end 0")'

    select case (eventType)
    case (heavyIon)
       write(this%iFile,f1) this%nEvent, b_Hi
    case (hadron)
       write(this%iFile,f1) this%nEvent, b_had
    case (LoPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Lo(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case (HiPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Hi(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case default
       write(this%iFile,f0) this%nEvent
    end select

  end subroutine

  !****************************************************************************
  !****s* EventOutput/Oscar_write_particle
  ! NAME
  ! subroutine Oscar_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in OSCAR 2013 format.
  !****************************************************************************
  subroutine Oscar_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU

    class(OscarOutputFile), intent(in) :: this
    type(particle), intent(in) :: part

    character(len=22), parameter :: f = '(9ES17.9,2I9)'
    integer :: PDGcode

    PDGcode = KFfromBUU(part)
    write(this%iFile,f) part%lastCollisionTime, part%position(1:3), &
                        sqrts(part), part%momentum(0:3), PDGcode, part%number

  end subroutine


!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/OscarExt_open
  ! NAME
  ! subroutine OscarExt_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output according to the
  ! "OSCAR 2013" standard.
  !****************************************************************************
  subroutine OscarExt_open(this, pert, nCall, nTimeStep)
    class(OscarExtOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=4)  :: buf

    if (pert) then
      buf = 'Pert'
      this%iFile = 723
    else
      buf = 'Real'
      this%iFile = 724
    end if
    fName = 'EventOutput.' // trim(buf) // '.oscar'

    if (nCall == 1) then
      ! open file for the first time
      open(this%iFile, file=fName, status='unknown')
      ! write header
      write(this%iFile,'(A)') '#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID mass0 lastCollisionTime mother1 mother2 mother3 generation'
      write(this%iFile,'(A)') '# Units: fm fm fm fm GeV GeV GeV GeV GeV none none GeV fm none none none none'
      write(this%iFile,'(A)') '# File generated by GiBUU (https://gibuu.hepforge.org)'
    else
      ! append to exiting file
      open(this%iFile, file=fName, status='old', position='append')
    end if

  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_close
  ! NAME
  ! subroutine OscarExt_close(this)
  ! PURPOSE
  ! Close a file after outputting OSCAR 2013 event information.
  !****************************************************************************
  subroutine OscarExt_close(this)
    class(OscarExtOutputFile), intent(in) :: this

    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_write_event_header
  ! NAME
  ! subroutine OscarExt_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an OSCAR 2013 event, including the number of particles
  ! and the event weight.
  !****************************************************************************
  subroutine OscarExt_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(OscarExtOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(*), parameter :: f1 = '("# event ",I9," out",I4," weight",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," out",I4)'

    this%nEvent = nEvent
    this%iFE = -1
    if (present(iFE)) this%iFE = iFE

    if (present(wgt)) then
       write(this%iFile,f1) this%nEvent, nParts, wgt
    else
       write(this%iFile,f0) this%nEvent, nParts
    end if
  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_write_event_footer
  ! NAME
  ! subroutine OscarExt_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes an OSCAR 2013 event.
  !****************************************************************************
  subroutine OscarExt_write_event_footer(this)
    use inputGeneral, only: eventType
    use eventtypes, only: heavyIon, hadron, LoPion, HiPion
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use initPion, only: getImpact_Lo => getImpact
    use initHiPion, only: getImpact_Hi => getImpact

    class(OscarExtOutputFile), intent(in) :: this

    character(*), parameter :: f1 = '("# event ",I9," end 0 impact",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," end 0")'

    select case (eventType)
    case (heavyIon)
       write(this%iFile,f1) this%nEvent, b_Hi
    case (hadron)
       write(this%iFile,f1) this%nEvent, b_had
    case (LoPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Lo(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case (HiPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Hi(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case default
       write(this%iFile,f0) this%nEvent
    end select

  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_write_particle
  ! NAME
  ! subroutine OscarExt_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in OSCAR 2013 extended format.
  !****************************************************************************
  subroutine OscarExt_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU
    use history, only: history_getParents,history_getGeneration

    class(OscarExtOutputFile), intent(in) :: this
    type(particle), intent(in) :: part

    character(len=30), parameter :: f = '(9ES17.9,2I9,2ES17.9,4I9)'
    integer :: PDGcode, generation, k
    integer, dimension(1:3) :: parents

    PDGcode = KFfromBUU(part)
    parents = history_getParents(part%history)
    generation=history_getGeneration(part%history)
    do k=1,2
       if (parents(k)>200) parents(k)=200-parents(k)
    end do
    write(this%iFile,f) -99.9, part%position(1:3), &
         sqrts(part), part%momentum(0:3), PDGcode, part%number, &
         part%mass, part%lastCollisionTime, &
         parents(1:3), generation

  end subroutine


!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/Shanghai_open
  ! NAME
  ! subroutine Shanghai_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_open(this, pert, nCall, nTimeStep)
    class(ShanghaiOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (pert) then
      buf1 = '.Pert.'
      this%iFile = 725
    else
      buf1 = '.Real.'
      this%iFile = 726
    end if
    write(buf2,'(I8.8)') nTimeStep
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.dat'

    ! open file
    open(this%iFile, file=fName, status='unknown', position='append')

  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_close
  ! NAME
  ! subroutine Shanghai_close(this)
  ! PURPOSE
  ! Close a file after outputting event information in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_close(this)
    class(ShanghaiOutputFile), intent(in) :: this

    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_write_event_header
  ! NAME
  ! subroutine Shanghai_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an event in Shanghai2014 format,
  ! including the number of particles and the event weight.
  !****************************************************************************
  subroutine Shanghai_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(ShanghaiOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(len=20) :: f

    if (present(wgt)) then
      f = '(A,I4,A,I4,A,E14.6)'
      write(this%iFile,f) '# event ', nEvent, ' out ', nParts, ' weight', wgt
    else
      f = '(A,I4,A,I4)'
      write(this%iFile,f) '# event ', nEvent, ' out ', nParts
    end if
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_write_event_footer
  ! NAME
  ! subroutine Shanghai_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a event in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_write_event_footer(this)
    class(ShanghaiOutputFile), intent(in) :: this
    ! empty, no footer!
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_write_particle
  ! NAME
  ! subroutine Shanghai_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_write_particle(this, part)
    use particleDefinition
    use minkowski, only: abs4

    class(ShanghaiOutputFile), intent(in) :: this
    type(particle), intent(in) :: part

    character(len=22), parameter :: f = '(2I9,7ES17.9)'
    integer :: fact

    if (part%antiparticle) then
      fact = -1
    else
      fact = 1
    end if

    write(this%iFile,f) fact*part%ID, part%charge, abs4(part%momentum), &
                        part%position(1:3), part%momentum(1:3)

  end subroutine

!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/Root_open
  ! NAME
  ! subroutine Root_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output in Root format.
  !****************************************************************************
  subroutine Root_open(this, pert, nCall, nTimeStep)
    class(RootOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (pert) then
       buf1 = '.Pert.'
    else
       buf1 = '.Real.'
    end if
    write(buf2,'(I8.8)') nCall
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.root'

    ! open file
    call rootinit(fName)

  end subroutine Root_open


  !****************************************************************************
  !****s* EventOutput/Root_close
  ! NAME
  ! subroutine Root_close(this)
  ! PURPOSE
  ! Close a file after outputting event information in Root format.
  !****************************************************************************
  subroutine Root_close(this)
    class(RootOutputFile), intent(in) :: this

    call rootclose()

  end subroutine Root_close


  !****************************************************************************
  !****s* EventOutput/Root_write_event_header
  ! NAME
  ! subroutine Root_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an event in Root format,
  ! including the number of particles and the event weight.
  !****************************************************************************
  subroutine Root_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(RootOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    this%weight = 1.0
    if (present(wgt)) this%weight=wgt

  end subroutine Root_write_event_header


  !****************************************************************************
  !****s* EventOutput/Root_write_event_footer
  ! NAME
  ! subroutine Root_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a event in Root format.
  !****************************************************************************
  subroutine Root_write_event_footer(this)
    class(RootOutputFile), intent(in) :: this

    call rootaddevent(this%weight)

  end subroutine Root_write_event_footer


  !****************************************************************************
  !****s* EventOutput/Root_write_particle
  ! NAME
  ! subroutine Root_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in Root format.
  !****************************************************************************
  subroutine Root_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU
    use history, only: history_getGeneration, history_getParents
    
    class(RootOutputFile), intent(in) :: this
    type(particle), intent(in) :: part

    real, dimension(1:3) :: pos

    integer :: KF
    integer :: parents(1:3)
    integer :: gen
    
    KF = KFfromBUU(part)
    gen = history_getGeneration(part%history)
    parents = history_getParents(part%history)
    
    if (useProductionPos) then
       pos = getProductionPos(part)
    else
       pos = part%position
    end if

    call rootaddparticle(KF, part%ID, part%charge, &
      part%number, part%history, &
      part%momentum(1),part%momentum(2),part%momentum(3), part%momentum(0), &
      pos(1),pos(2),pos(3), &
      part%event(1),part%event(2),part%firstEvent)

  end subroutine Root_write_particle

  !****************************************************************************
  !****s* EventOutput/Root_write_additionalInfo
  ! NAME
  ! subroutine Root_write_additionalInfo(this, iFE)
  ! PURPOSE
  ! add additional info about the event, depending on eventtype.
  !
  ! This routine tries to find additional information about the event.
  ! It tries routines for different event types, which only return
  ! some information, if it was really stored.
  !
  ! The following cases are handled:
  ! * For eventtype "HiLep", the following line is added:
  !     # 14 nu Q2 eps phiLepton Eventtype
  !   (14 is the magic number of "HiLepton")
  ! * For eventtype "neutrino", the following line is added:
  !     # 5 Eventtype Weight momLepIn(0:3) momLepOut(0:3) momNuc(0:3)
  !   (5 is the magic number for neutrino events)
  ! * For eventtype "heavyIon", the following line is added:
  !     # 1 b
  !   (1 is the magic number of "heavyIon", b is the impact parameter in fm)
  ! * For eventtype "hadron", the following line is added:
  !     # 300 b
  !   (300 is the magic number of "hadron", b is the impact parameter in fm)
  !****************************************************************************
  subroutine Root_write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use inputGeneral, only: eventType, numEnsembles, num_runs_SameEnergy
    use initNeutrino, only: process_ID, flavor_ID
    use nucleusDefinition
    use nucleus, only: getTarget
    use eventtypes, only: hiLepton, neutrino, heavyIon, hadron
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use FreezeoutAnalysis, only: getFreezeoutAnalysis_Pert
    use PIL_freezeout, only: PIL_freezeout_GET

    class(RootOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode
    type(tNucleus), pointer :: targetNuc

    ! variable have to have attribute 'save' due to rootadd... routine:
    real, save :: weight,nu,Q2,eps,phiL
    integer, save :: evtType, chrg_nuc
    real,dimension(0:3),save :: momLepIn, momLepOut, momBos, momNuc

!    integer :: history
!    logical :: escaped
    real,save :: tmp

    targetNuc => getTarget()

    select case (eventType)
    case (heavyIon)
       tmp = b_HI
       call rootadddouble(tmp,"b")
    case (hadron)
       tmp = b_had
       call rootadddouble(tmp,"b")
    case (neutrino)
      if (.not. present(iFE)) return
      if (NeutrinoProdInfo_Get(iFE,evtType,Weight,momLepIn,momLepOut,momBos,momNuc,chrg_nuc)) then
         call rootaddint(numEnsembles, "numEnsembles")
         call rootaddint(num_runs_SameEnergy, "numRuns")
         call rootaddint(targetNuc%mass, "nucleus_A")
         call rootaddint(targetNuc%charge, "nucleus_Z")
         call rootaddint(flavor_ID, "flavor_ID")
         call rootaddint(process_ID, "process_ID")
         call rootaddint(evtType, "evType")
         call rootadddouble(momLepIn(0), "lepIn_E")
         call rootadddouble(momLepIn(1), "lepIn_Px")
         call rootadddouble(momLepIn(2), "lepIn_Py")
         call rootadddouble(momLepIn(3), "lepIn_Pz")
         call rootadddouble(momLepOut(0), "lepOut_E")
         call rootadddouble(momLepOut(1), "lepOut_Px")
         call rootadddouble(momLepOut(2), "lepOut_Py")
         call rootadddouble(momLepOut(3), "lepOut_Pz")
         call rootadddouble(momNuc(0), "nuc_E")
         call rootadddouble(momNuc(1), "nuc_Px")
         call rootadddouble(momNuc(2), "nuc_Py")
         call rootadddouble(momNuc(3), "nuc_Pz")
         call rootaddint(chrg_nuc, "nuc_charge")
         
      end if
    case (hiLepton)
       if (.not. present(iFE)) return
       if (EventInfo_HiLep_Get(0,iFE,Weight,nu,Q2,eps,evtType,phi_Lepton=phiL)) then
          call rootaddint(evtType, "evType")
          call rootadddouble(nu, "nu")
          call rootadddouble(Q2, "Q2")
          call rootadddouble(eps, "eps")
          call rootadddouble(phiL, "phiL")

       end if
    end select

  end subroutine Root_write_additionalInfo



!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/write_real
  ! NAME
  ! subroutine write_real(this, parts)
  ! PURPOSE
  ! Do the actual printout for real particles.
  ! NOTES
  ! For the case of real particles, one event simply corresponds to one
  ! ensemble.
  !****************************************************************************
  subroutine write_real(this, parts)
    use particleDefinition
    use IdTable, only: EOV, NOP
    use inputGeneral, only: current_run_number

    class(EventOutputFile), intent(in) :: this
    type(particle), intent(in), dimension(:,:), target :: Parts

    integer :: iEns, iPart, NUP, nEnsembles, iEvent
    integer, dimension(:), allocatable :: nParts

    allocate(nParts(1:size(Parts,dim=1)))
    nParts=0

    ! count particles per ensemble
    do iEns = 1,size(Parts,dim=1)
      do iPart = 1,size(Parts,dim=2)
        if (Parts(iEns,iPart)%ID==EOV) exit
        if (Parts(iEns,iPart)%ID==NOP) cycle
        nParts(iEns) = nParts(iEns) + 1
      end do
    end do

    ! Loop over all events and write them to file:
    nEnsembles = size(Parts,dim=1)
    do iEns = 1, nEnsembles
       NUP = nParts(iEns) ! number of particles
       if (NUP == 0) cycle

       iEvent = (current_run_number-1)*nEnsembles + iEns
       call this%write_event_header(NUP, iEvent)  ! no weight here

       do iPart = 1,size(Parts,dim=2)
         if (Parts(iEns,iPart)%ID==EOV) exit
         if (Parts(iEns,iPart)%ID==NOP) cycle
         call this%write_particle(Parts(iEns,iPart))
       end do

       call this%write_additionalInfo()

       call this%write_event_footer()
    end do

  end subroutine write_real


  !****************************************************************************
  !****s* EventOutput/ValueListAllocate
  ! NAME
  ! subroutine ValueListAllocate(n1)
  ! PURPOSE
  ! Do the allocation stuff for the Particle Info List.
  !****************************************************************************
  subroutine ValueListAllocate(ValueList, n1)
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_INIT

    type(tParticleList), allocatable :: ValueList(:)
    integer, intent(in) :: n1

    integer :: n0,i
    type(tParticleList),allocatable :: L0(:)

    if (.not.allocated(ValueList)) then
        allocate(ValueList(n1))
        do i=1,n1
          call ParticleList_INIT(ValueList(i))
        end do
        return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    do i=1,n0
        L0(i)%first => ValueList(i)%first
        L0(i)%last  => ValueList(i)%last
        L0(i)%nEntries = ValueList(i)%nEntries
    end do
    deallocate(ValueList)
    allocate(ValueList(n1))
    do i=1,n0
        ValueList(i)%first => L0(i)%first
        ValueList(i)%last  => L0(i)%last
        ValueList(i)%nEntries = L0(i)%nEntries
    end do
    do i=n0+1,n1
        call ParticleList_INIT(ValueList(i))
    end do
    deallocate(L0)

  end subroutine ValueListAllocate


  !****************************************************************************
  !****s* EventOutput/write_pert
  ! NAME
  ! subroutine write_pert(this, parts)
  ! PURPOSE
  ! Do the actual printout for perturbative particles.
  ! NOTES
  ! We have to sort the particles according their "firstevent" field.
  ! Therefore we allocate an array of "tParticleList". Unfortunately we can
  ! not use the "firstevent" entry directly as array index, since this
  ! is not starting with 1 and continously increasing for all kind
  ! of eventtypes. Therefore we (ab)use the module "PILIndex", which
  ! implements methods of "indexing". (We do not use the possibility of
  ! reallocating as provided by the module "PILIndex".)
  !****************************************************************************
  subroutine write_pert(this, parts)
    use particleDefinition
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_APPEND, ParticleList_CLEAR
    use PILIndex, only: tIndexList, PILIndex_PUT
    use inputGeneral, only: current_run_number

    class(EventOutputFile), intent(in) :: this
    type(particle), intent(in), dimension(:,:), target :: Parts

    type(tIndexList), save :: IndexList
    type(tParticleList), allocatable, save :: ValueList(:)

    integer :: i,iEns,iPart, iFE,iiFE, NUP, iEvent

    type(particle), pointer :: pPart
    type(tParticleListNode),Pointer  :: pNode

    ! Clean up the arrays:
    if (allocated(ValueList)) then
       do i=1,size(ValueList)
          call ParticleList_CLEAR(ValueList(i))
       end do
    end if
    IndexList%nEntry = 0

    ! Loop over all particles and group them according their first event:
    do iEns = 1,size(Parts,dim=1)
       PartLoop:do iPart = 1,size(Parts,dim=2)
          pPart => Parts(iEns,iPart)
          if (pPart%ID <= 0) cycle PartLoop

          iFE = pPart%firstEvent
          if (iFE.eq.0) cycle PartLoop ! particle did not interact !
          iiFE = PILIndex_PUT(IndexList, iFE, 'LesHouches')
          if (iiFE>0) then
             call ParticleList_APPEND(ValueList(iiFE), pPart)
          else
             call ValueListAllocate(ValueList, size(IndexList%PartNumber))
             call ParticleList_APPEND(ValueList(-iiFE), pPart)
          end if

       end do PartLoop
    end do

    ! Loop over all events and write them to file:
    if (.not.allocated(ValueList)) return

    do iiFE=1,size(ValueList)
       NUP = ValueList(iiFE)%nEntries ! number of particles
       if (NUP == 0) cycle

       pNode => ValueList(iiFE)%first
       iFE = pNode%V%firstEvent
       iEvent = (current_run_number-1)*size(ValueList) + iiFE
       call this%write_event_header(NUP, iEvent, pNode%V%perweight, iFE)

       do
          if (.not. associated(pNode)) exit
          pPart => pNode%V
          call this%write_particle(pPart)
          pNode => pNode%next
       end do

       call this%write_additionalInfo(iFE, ValueList(iiFE)%first)

       call this%write_event_footer()
    end do

  end subroutine write_pert


end module EventOutput
