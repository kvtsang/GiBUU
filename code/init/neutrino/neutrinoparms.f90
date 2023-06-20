module neutrinoParms

use CALLSTACK, only: TRACEBACK

implicit none

  !****************************************************************************
  !****g* neutrinoParms/Wtrans
  ! SOURCE
  real, save:: Wtrans = 3.0
  ! PURPOSE
  ! W for transition from Bloom-Gilman Paramerization to PYTHIA DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/Norm2pi
  ! SOURCE
  real, save:: Norm2pi = 1.2
  ! PURPOSE
  ! overall normalization factor for 2pi BG and Bloom-Gilman X-section
  !****************************************************************************

!****************************************************************************
  !****g* neutrinoParms/normRES
  ! SOURCE
  real, save:: normRES = 1.6
  ! PURPOSE
  ! overall normalization factor for neutrino-induced resonance contributions
  ! beyond the Delta
  !****************************************************************************
  
  !****************************************************************************
  !****g* neutrinoParms/normCB
  ! SOURCE
  real, save:: normCB = 2.0
  ! PURPOSE
  ! overall normalization factor for neutrino-induced Christy-Bosted contributions
  ! between 2 GeV and DIS onset
  !****************************************************************************
  
  !****************************************************************************
  !****g* neutrinoParms/T0
  ! SOURCE
  real, save:: T0 = 0.
  ! PURPOSE
  ! isospin or strength factor for 2pi BG
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/Q2cut
  ! SOURCE
  real, save::  Q2cut = 4.5
  ! PURPOSE
  ! threshold value
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/A00
  ! SOURCE
  real, save::  A00 = -0.43
  real, save::  A10 = -0.97
  real, save::  A01 = 0.87
  real, save::  A20 = 8.01
  real, save::  A11 = 4.59
  real, save::  A02 = 3.18
  real, save::  A30 = -10.42
  real, save::  A21 = 2.04
  real, save::  A12 = -9.14
  real, save::  A03 = -15.15

  real, save::  W00 = 1.55
  real, save::  x00 = 0.35
  ! PURPOSE
  ! coefficients in Taylor expansion of attenuation function Att(W,x)
  !****************************************************************************
   
  ! Now parameters used pimarily in 2p2h
 
  !****************************************************************************
  !****g* lepton2p2h/ME_Version
  ! PURPOSE
  ! indicate the type of matrix element parametrisation
  !
  ! SOURCE
  integer, save :: ME_Version = 4
  ! NOTES
  ! possible values:
  ! * 1: const ME_Norm_XX  ! const for CC  fitted to MiniBooNE is 1.8e-6
  ! * 2: constant transverse and decreasing with Enu
  ! * 3: "Dipole transverse" transverse,  fall with Q2 as 4-th power
  ! * 4: MEC from E. Christy (8/2015), with parametrization for longitudinal
  ! * 5: MEC from Bosted arXiV:1203.2262, with parametrization for longitudinal
  ! * 6: MEC additional parametrization, with parametrization for longitudinal
  ! *    not yet implemented
  !****************************************************************************


  ! The following are all tunable strength parameters for 2p2h hadronic
  ! structure functions. Default is no tuning, i.e. all parameters = 1.
  ! except for ME_Long, for which default is =0 (no longitudinal component)



  !****************************************************************************
  !****g* lepton2p2h/ME_Norm_QE
  ! PURPOSE
  ! Overall strength of 2p2h matrix element with 2N out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Norm_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 gives the coded strength
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Norm_Delta
  ! PURPOSE
  ! Overall strength of 2p2h matrix element with NDelta out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Norm_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Mass_QE
  ! PURPOSE
  ! Cutoff-mass in some parametrizations of 2p2h matrix element for NN out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Mass_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Mass_Delta
  ! PURPOSE
  ! Cutoff-mass in some parametrizations of matrix element for NDelta out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Mass_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_Transversity
  ! PURPOSE
  ! Parametrisation of structure functions
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_Transversity = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value = 1 chooses structure function W2 so that 2p2h is pure transverse
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_LONG
  ! PURPOSE
  ! Parametrization of structure functions
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_LONG = (/0.0, 0.0, 0.0/)
  ! NOTES
  ! The value = 0 turns any additional longitudinal contribution
  ! to structure funct. W2 off
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_W3
  ! PURPOSE
  ! Overall strength factor for structure function W3
  !
  ! only for (CC,NC)
  ! SOURCE
  real,dimension(1:3), save :: ME_W3 = (/0.,1.0, 1.0/)
  ! NOTES
  ! overall strength parameter for structure function W3
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/ME_ODW
  ! PURPOSE
  ! switch for choosing the connection between structure functions
  ! W1(electron) and W1(neutrino) and W3(neutrino):
  ! * 1: for expressions from Martini et al
  ! * 2: for expressions from O'Connell et al
  ! * 3: for expression with relativistic version of O'Connell
  !
  ! only for (CC,NC)
  ! SOURCE
  integer, save :: ME_ODW = 2
  ! NOTES
  ! * O'Connell et al: PR C6 (1972) 719
  ! * Martini et al: PR C80 (2009) 065501
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/inmedW
  ! PURPOSE
  ! Controls which inv mass W is used in parametrization of 2p2h W1
  !
  ! SOURCE
  integer, save :: inmedW = 1
  ! NOTES
  ! * 1: W = static inv. mass in 2p2h parametrization of W1
  ! * 2: W = inv mass for Fermi moving nucleons in potential
  ! * 3: W = inv mass for Fermi moving nucleons without potential
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/T
  ! PURPOSE
  ! target isospin, affects only neutrino 2p2h structure function
  !
  ! SOURCE
  real, save ::  T = 99
  ! NOTES
  ! * T = 0, 1 , ...
  ! * T = 99 gives T = (N-Z)/2
  !****************************************************************************

  !****************************************************************************
  !****g* lepton2p2h/Adep
  ! PURPOSE
  ! Switch for A-dependence of 2p2h structure function
  !
  ! SOURCE
  integer, save ::  Adep = 2
  ! NOTES
  ! * 1: A-dependence for zero-range force (Mosel, Gallmeister, 2016)
  ! * 2: linear A-dependence, normalized to C12
  !****************************************************************************



contains

  !****************************************************************************
  !****is* Neutrinoparms/neutrino_readInput
  ! NAME
  ! subroutine neutrino_readInput
  !****************************************************************************
  subroutine neutrino_readInput
    use output

    integer :: ios
    logical, save:: initflag = .true.
    !************************************************************************
    !****n* neutrinoParms/nl_Neutrino2piBack
    ! NAME
    ! NAMELIST /nl_Neutrino2piBack/
    ! PURPOSE
	! Parameters for transition from RES to DIS:
	! * Wtrans
    ! * Norm2pi
	! * normRES
	! * normCB
    ! Includes parameters of neutrino 2 pion background:
    ! * Q2cut
    ! * A00
    ! * A10
    ! * A01
    ! * A02
    ! * A11
    ! * A20
    ! * A03
    ! * A12
    ! * A21
    ! * A30
    ! * W00
    ! * x00
	! * T0
    !************************************************************************
    NAMELIST /nl_neutrino2piBack/ Wtrans, &
         & Norm2pi,normRES,normCB,Q2cut,A00,A10,A01,A02,A11,A20,A03,A12,A21,A30,W00,x00,T0
		 
	if (.not.initFlag) return	 
		 
    call Write_ReadingInput('neutrinoParms',0)
    rewind(5)
    read(5,nml=nl_neutrino2piBack,IOSTAT=ios)
    call Write_ReadingInput("neutrinoParms",0,ios)
    write(*,*) 'Wtrans=',Wtrans,'Norm2pi= ',Norm2pi,'normRES=',normRES,'normCB=',normCB
    write(*,fmt='(12(F7.2))') Q2cut,A00,A10,A01,A20,A11,A02,A30,A21,A12,A03,T0
    call Write_ReadingInput('neutrinoParms',1)

    initFlag=.false.
  end subroutine neutrino_readInput
  
   !****************************************************************************
  !****is* lepton2p2h/readInput
  ! NAME
  ! subroutine readInput
  !****************************************************************************
  subroutine twoptwoh_readInput

    use output
    use AZN

    integer :: ios, i
	
	logical, save:: initflag = .true.
    


    !**************************************************************************
    !****n* lepton2p2h/Lepton2p2h
    ! NAME
    ! NAMELIST /Lepton2p2h/
    ! PURPOSE
    ! Includes parameters for 2p2h events:
    ! * ME_Version
    ! * ME_Norm_QE
    ! * ME_Norm_Delta
    ! * ME_Mass_QE
    ! * ME_Mass_Delta
    ! * ME_Transversity
    ! * ME_LONG
    ! * ME_W3
    ! * ME_ODW
    ! * inmedW
    ! * T
    ! * Adep
    !**************************************************************************
    NAMELIST /lepton2p2h/ ME_Version, &
         ME_Norm_QE, ME_Norm_Delta, &
         ME_Mass_QE, ME_Mass_Delta,ME_Transversity,ME_LONG, &
         ME_W3,ME_ODW,inmedW,T,Adep

    if (.not.initFlag) return

    call Write_ReadingInput('lepton2p2h',0)
    rewind(5)
    read(5,nml=lepton2p2h,IOSTAT=ios)
    call Write_ReadingInput("lepton2p2h",0,ios)

    call AZNsub(Atarget,Ntarget,Ztarget)

    select case (ME_Version)

    case (1)
       write(*,'(A)') 'ME1  const =  4.0e-6 * ME_Norm_XX'
       !  case 1 is model-I in "Lalakulich Gallmeister Mosel PRC86(2012)014614"

    case (2)
       write(*,'(A)') 'ME2  4.8e4 * 0.635 / Enu^2 * ME_Norm_XX  &
            &  in transverse part only, decreasing with Enu'
       !  case 2 is model-II from "Lalakulich Gallmeister Mosel PRC86(2012)014614"

    case (3)
       write(*,'(A)') 'ME3  dipole parameterization of the "form factor" &
            & ME=8e4*(1+Q2/MA2)^{-4},transverse part only'
       ! case 3 gives a good description of MiniBooNE data with MA ~ 1.5 GeV

    case (4)
       write(*,'(A)') 'ME 4, parametrization of structure functions W1,W2,W3,&
            & W1 from E. Christy'
       !       write(*,*) 'ME_W3   =', ME_W3

       if(T == 99) T = (Ntarget - Ztarget)/2.0

       write(*,*) 'inmedW=',inmedW,'   T=',T
       write(*,*) 'ME_ODW=',ME_ODW
       write(*,*) 'Adep=',Adep

    case (5)
       write(*,'(A)') 'parametrization of structure functions W1,W2,W3,&
            & W1 from Bosted'

    case (6)
       write(*,'(A)') 'parametrization of structure functions W1,W2,W3,&
            & W1 new parametrization'

    case default
       write(*,*) 'ME_Version = ',ME_Version
       call TRACEBACK('wrong value for ME_Version')
    end select


    write(*,'(A)') 'parameters ( N N final state ):  [i=EM,CC,NC]'
    do i=1,3
       write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
            & ME_Norm_QE(i),ME_Mass_QE(i),ME_Transversity(i),ME_LONG(i)
    end do
    write(*,'(A)') 'parameters ( N Delta final state ):'
    do i=1,3
       write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
            & ME_Norm_Delta(i),ME_Mass_Delta(i)
    end do

    call Write_ReadingInput('lepton2p2h',1)

    initFlag=.false.

  end subroutine twoptwoh_readInput
  
  
  end module neutrinoParms