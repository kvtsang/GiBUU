!******************************************************************************
!****m* /formfactors_A_main
! NAME
! module formfactors_A_main
!
! PURPOSE
! This module administrates the formfactors for gamma* N -> N pi.
!
! INPUTS
! None
!******************************************************************************
module formfactors_A_main

 use CallStack, only: TraceBack


  implicit none
  private

  public :: getA
  public :: cleanUp
  public :: getMinW

  !***************************************************************************
  !****t* formfactors_A_main/tFormFacArr
  ! NAME
  ! type tFormFacArr
  ! PURPOSE
  ! store informations from MAID form factor files
  !
  ! SOURCE
  !
  type, public :: tFormFacArr
     private

     complex, dimension(:,:), allocatable :: A
     real, dimension(:), allocatable :: W, theta, Q2
     real :: dW, dTheta, dQ2 = 0.
     real :: maxW, maxTheta, maxQ2 = 0.
     integer :: nW, nTheta, nQ2 = 0
     character(len= 4) :: sShort
     character(len=16) :: sLong
     logical :: initFlag = .true.

   contains

     private
     procedure, public :: getA => FormFacArr_getA
     procedure, public :: cleanup => FormFacArr_cleanup
     procedure :: init
     procedure :: read
     procedure :: get_Filename
     procedure :: getMinW => FormFacArr_getMinW

  end type tFormFacArr

  interface tFormFacArr
     module procedure constructor
  end interface tFormFacArr
  !****************************************************************************


  !****************************************************************************
  !****g* formfactors_A_main/which_MaidVersion
  ! SOURCE
  !
  integer, save :: which_MaidVersion=2
  !
  ! PURPOSE
  ! choice of MAID version:  1=2003, 2=2007
  !****************************************************************************


  integer, parameter :: MAID2003=1
  integer, parameter :: MAID2007=2

  !****************************************************************************
  !****g* formfactors_A_main/fineGridQ2
  ! SOURCE
  !
  logical, save :: fineGridQ2 = .false.
  !
  ! PURPOSE
  ! if .true. use MAID tables with a finer grid in Q^2 (dQ2 = 0.05GeV^2)
  ! instead of default dQ2 = 0.1GeV^2
  !****************************************************************************


  logical, save :: initFlag = .true.
  type(tFormFacArr) :: FormFacArr_pi0p
  type(tFormFacArr) :: FormFacArr_pi0n
  type(tFormFacArr) :: FormFacArr_pipn
  type(tFormFacArr) :: FormFacArr_pimp


contains
  !****************************************************************************
  !****if* formfactors_A_main/constructor
  ! NAME
  ! function constructor(sShort, sLong) result(FormFacArr)
  !
  ! PURPOSE
  ! Initialies a FormFacArr, but only set the strings.
  ! Array are not allocated, this will be done later while reading in from the
  ! file
  !
  ! INPUTS
  ! * character(len=*), intent(in) :: sShort -- the 4 letters used in the
  !   filename
  ! * character(len=*), intent(in) :: sLong -- a string for nice printout
  !****************************************************************************
  function constructor(sShort, sLong) result(FormFacArr)
    type(tFormFacArr) :: FormFacArr

    character(len=*), intent(in) :: sShort
    character(len=*), intent(in) :: sLong

    FormFacArr%sShort = sShort(1:4)
    FormFacArr%sLong = sLong(1:16)
    FormFacArr%initFlag = .true.

    ! arrays stay unallocated at this point

  end function constructor

  !****************************************************************************
  !****f* formfactors_A_main/FormFacArr_getA
  ! NAME
  ! function FormFacArr_getA(this, theta, W2, Q2)
  !
  ! PURPOSE
  ! Initialies a FormFacArr, but only set the strings.
  ! Array are not allocated, this will be done later while reading in from the
  ! file
  !
  ! INPUTS
  ! * class(tFormFacArr), intent(inOut) :: this -- the instance to be used
  ! * real, intent(in) :: theta,W2,Q2 -- the parameters in the array
  !
  ! OUTPUT
  ! complex, dimension(1:6) --  The complex amplitudes A_1 to A_6 in units of:
  ! * [GeV**-2] for A(1)
  ! * [GeV**-4] for A(2) and A(5)
  ! * [GeV**-3] for A(3), A(4) and A(5)
  !****************************************************************************
  function FormFacArr_getA(this, theta, W2, Q2) result(getA)
    class(tFormFacArr), intent(inOut) :: this
    real, intent(in) :: theta,W2,Q2
    complex, dimension(1:6) :: getA

    logical, parameter:: linearInterpolation=.true.
    real :: W
    real :: W_,theta_,Q2_ ! input values resticted to valid range
    integer :: ind, i
    integer, dimension(8) :: iCube
    real, dimension(8) :: cubeW, cubeTheta, cubeQ2
    real :: fac, facSum

    if (this%initFlag) call this%init()

    W = sqrt(W2)

    if (linearInterpolation) then
       ! restrict to tabulated region:
       W_ = min(max(this%W(1),W),this%maxW)
       theta_ = min(max(this%theta(1),theta),this%maxTheta)
       Q2_ = min(max(this%Q2(1),Q2),this%maxQ2)

       ! get corner of embedding cube:
       ind = get_i_lower(theta_, W_, Q2_, iCube)

       ! copy values at the corners:
       cubeW = this%W(iCube)
       cubeTheta = this%Theta(iCube)
       cubeQ2 = this%Q2(iCube)

       facSum = 0
       getA = 0
       do i=1,8
          fac =  ( 1. - abs(cubeTheta(i)-theta_)/this%dTheta ) &
               * ( 1. - abs(cubeW(i)-W_)/this%dW ) &
               * ( 1. - abs(cubeQ2(i)-Q2_)/this%dQ2 )
          facSum = facSum + fac
          getA = getA + this%A(:,iCube(i))*fac
       end do
       if (abs(facSum - 1.) > 0.001) then
          write(*,*) 'facSum=',facSum
          write(*,*) W_,theta_,Q2_
          do i=1,8
             write(*,*) cubeW(i),cubeTheta(i),cubeQ2(i), &
                  ( 1. - abs(cubeTheta(i)-theta_)/this%dTheta ) &
                  * ( 1. - abs(cubeW(i)-W_)/this%dW ) &
                  * ( 1. - abs(cubeQ2(i)-Q2_)/this%dQ2 )
          end do
          call Traceback('sum of factors not 1!')
       end if

    else
       ! all restriction is done by the index calculation routine
       ind = get_i_closest(theta, W, Q2)
       getA = this%A(:,ind)
    end if

  contains
    !*************************************************************************
    !****if* FormFacArr_getA/get_i_closest
    ! NAME
    ! integer function get_i_closest(theta, W, Q2)
    ! PURPOSE
    ! return the bin, which is closest to the given parameters
    !
    ! lower and upper bounds are respected, so there is no extrapolation, but
    ! the range is truncated at the borders.
    !*************************************************************************
    integer function get_i_closest(theta, W, Q2)
      real, intent(in) :: theta, W, Q2

      integer :: iTheta, iW, iQ2, ii

      iTheta = nint((theta-this%theta(1))/this%dTheta)+1
      iW =     nint((W-this%W(1))/this%dW)
      iQ2 =    nint((Q2-this%Q2(1))/this%dQ2)

      ! we allow nTheta,nW,nQ2 different values.
      ! iW and iQ2 start at 0, therefore "-1"
      iTheta = min(max(1,iTheta),this%nTheta)
      iW =     min(max(0,iW),this%nW-1)
      iQ2 =    min(max(0,iQ2),this%nQ2-1)

      ii = (iQ2*this%nW + iW)*this%nTheta + iTheta
      if (ii>size(this%W)) &
           call Traceback('calculated index too large!')
      get_i_closest = ii

    end function get_i_closest

    !*************************************************************************
    !****if* FormFacArr_getA/get_i_lower
    ! NAME
    ! integer function get_i_lower(theta, W, Q2, iCube)
    ! PURPOSE
    ! return the bin, which is smaller in all directions
    !
    ! lower and upper bounds are respected, so there is no extrapolation, but
    ! the range is truncated at the borders.
    !
    ! INPUTS
    ! * real :: theta, W, Q2 -- parameters
    ! * integer, dimension(8), OPTIONAL :: iCube -- indices of corners of cube
    !*************************************************************************
    integer function get_i_lower(theta, W, Q2, iCube)
      real, intent(in) :: theta, W, Q2
      integer, dimension(8), intent(inOut), OPTIONAL :: iCube

      integer :: iTheta, iW, iQ2, ii
      integer :: i,j,k

      iTheta = int((theta-this%theta(1))/this%dTheta)+1
      iW =     int((W-this%W(1))/this%dW)
      iQ2 =    int((Q2-this%Q2(1))/this%dQ2)

      ! we allow nTheta-1,nW-1,nQ2-1 different values.
      ! iW and iQ2 start at 0, therefore "-2"
      iTheta = min(max(1,iTheta),this%nTheta-1)
      iW = min(max(0,iW),this%nW-2)
      iQ2 = min(max(0,iQ2),this%nQ2-2)

      ii = (iQ2*this%nW + iW)*this%nTheta + iTheta
      if (ii>size(this%W)) &
           call Traceback('calculated index too large!')
      get_i_lower = ii

      if (present(iCube)) then
         do i=0,1
            do j=0,1
               do k=0,1
                  iCube(4*i+2*j+k+1) = &
                       ((iQ2+i)*this%nW + (iW+j))*this%nTheta + (iTheta+k)
               end do
            end do
         end do
         if (iCube(8)>size(this%W)) then
            write(*,*) 'iCube=',iCube
            call Traceback('calculated index too large!')
         end if
      end if

    end function get_i_lower

  end function FormFacArr_getA

  !****************************************************************************
  !****is* formfactors_A_main/init
  ! NAME
  ! subroutine init(this)
  !
  ! PURPOSE
  ! Initialies a FormFacArr by also calling the read routine.
  !****************************************************************************
  subroutine init(this)
     use output, only: Write_InitStatus

    class(tFormFacArr), intent(inOut) :: this

    if (.not.this%initFlag) return
    if (initFlag) call readInput ! global init
    this%initFlag = .false.

    call Write_InitStatus('MAID Amplitudes for '//this%sLong, 0)
    call this%read()
    write(*,*)
    write(*,'(A,3I9)')   'nSteps in W,Theta,Q2 :', &
         this%nW,this%nTheta,this%nQ2
    write(*,'(A,3F9.3)') 'min W, Theta, Q2:     ', &
         this%W(1),this%theta(1),this%Q2(1)
    write(*,'(A,3F9.3)') 'max W, Theta, Q2:     ', &
         this%maxW, this%maxTheta, this%maxQ2
    write(*,'(A,3F9.3)') 'dW, dTheta, dQ2:      ', &
         this%dW,this%dTheta,this%dQ2

    call Write_InitStatus('MAID Amplitudes for '//this%sLong, 1)

  end subroutine init

  !****************************************************************************
  !****is* formfactors_A_main/read
  ! NAME
  ! subroutine read(this)
  !
  ! PURPOSE
  ! Do the real reading of the file
  !****************************************************************************
  subroutine read(this)
    use bzip

    class(tFormFacArr), intent(inOut) :: this

    type(bzFile) :: f
    character(len=300) :: line
    integer :: ios, i, ll, nLines

    real, dimension(:,:), allocatable :: spalte
    real, dimension(:), allocatable :: rW, rTheta, rQ2
    real :: rDummy1, rDummy2
    logical, parameter :: convertToGeV=.true.

    integer :: maxGridPoints=100000

    if (fineGridQ2) maxGridPoints = 200000
    allocate(spalte(1:12,1:maxGridPoints))
    allocate(rW(1:maxGridPoints))
    allocate(rTheta(1:maxGridPoints))
    allocate(rQ2(1:maxGridPoints))

    !===== Step 1: read the input file

    f = bzOpenR(trim(this%get_Filename()))

    this%nW = 1
    this%nTheta = 1
    this%nQ2 = 1
    i = 1
    do
       ! for speedup reasons, we assume, that all lines in the file have
       ! the same length ll (this is an output, but also an input parameter)
       ! if we found a line containing data, i is increased and ll is constant
       if (i==1) ll = 0
       !ll = 0 ! for temporary use for non constant line length
       call bzReadLine(f,line,ll)
       if (ll==0) cycle

       read(line(1:ll),*,IOSTAT=ios)  rW(i), rTheta(i), rQ2(i), &
            rDummy1, rDummy2, spalte(1:12,i)

       if (ios.eq.0) then
          ! determine grid spacings:
          if (i.eq.1) then
             this%maxW = rW(i)
             this%maxTheta = rTheta(i)
             this%maxQ2 =  rQ2(i)
          else if (i.gt.1) then
             if (rW(i)-this%maxW.gt.0.000001) then
                this%dW = rW(i)-this%maxW
                this%maxW = rW(i)
                this%nW = this%nW+1
             end if

             if (rTheta(i)-this%maxTheta.gt.0.000001) then
                this%dTheta = rTheta(i)-this%maxTheta
                this%maxTheta = rTheta(i)
                this%nTheta = this%nTheta+1
             end if

             if (rQ2(i)-this%maxQ2.gt.0.000001) then
                this%dQ2 = rQ2(i)-this%maxQ2
                this%maxQ2 = rQ2(i)
                this%nQ2 = this%nQ2+1
             end if
          end if
          i = i+1
          if (i > maxGridPoints) &
               call Traceback('Too many lines in file!')
       else
          write(*,'(i3,i3,A)') i, ios, line(1:ll)
       end if
       if (f%EOF) exit
    end do
    nLines = i-1

    call bzCloseR(f)

    !===== Step 2: copy and convert the values

    allocate(this%W(1:nLines))
    this%W = rW(1:nLines) / 1000.
    this%maxW = this%maxW / 1000.
    this%dW = this%dW / 1000.

    allocate(this%Theta(1:nLines))
    this%Theta = rTheta(1:nLines)

    allocate(this%Q2(1:nLines))
    this%Q2 = rQ2(1:nLines)

    allocate(this%A(1:6,1:nLines))
    this%A(1,1:nLines) = CMPLX(Spalte( 1,1:nLines),Spalte( 2,1:nLines))
    this%A(2,1:nLines) = CMPLX(Spalte( 3,1:nLines),Spalte( 4,1:nLines))
    this%A(3,1:nLines) = CMPLX(Spalte( 5,1:nLines),Spalte( 6,1:nLines))
    this%A(4,1:nLines) = CMPLX(Spalte( 7,1:nLines),Spalte( 8,1:nLines))
    this%A(5,1:nLines) = CMPLX(Spalte( 9,1:nLines),Spalte(10,1:nLines))
    this%A(6,1:nLines) = CMPLX(Spalte(11,1:nLines),Spalte(12,1:nLines))

    if (convertToGeV) then
       ! Convert everything from fm to GeV^-1
       ! 197 MeV fm =1 => fm=1/(0.197 GeV)
       this%A(1,:) = this%A(1,:)/(0.197**2)
       this%A(2,:) = this%A(2,:)/(0.197**4)
       this%A(3,:) = this%A(3,:)/(0.197**3)
       this%A(4,:) = this%A(4,:)/(0.197**3)
       this%A(5,:) = this%A(5,:)/(0.197**4)
       this%A(6,:) = this%A(6,:)/(0.197**3)
    end if

    if (allocated(spalte)) deallocate(spalte)
    if (allocated(rW)) deallocate(rW)
    if (allocated(rTheta)) deallocate(rTheta)
    if (allocated(rQ2)) deallocate(rQ2)

  end subroutine read

  !****************************************************************************
  !****is* formfactors_A_main/FormFacArr_cleanup
  ! NAME
  ! subroutine FormFacArr_cleanup(this)
  !
  ! PURPOSE
  ! Deallocate memory
  !****************************************************************************
  subroutine FormFacArr_cleanup(this)
    class(tFormFacArr), intent(inOut) :: this

    if (allocated(this%A)) deallocate(this%A)
    if (allocated(this%W)) deallocate(this%W)
    if (allocated(this%theta)) deallocate(this%theta)
    if (allocated(this%Q2)) deallocate(this%Q2)
    this%initFlag = .true.

  end subroutine FormFacArr_cleanup

  !****************************************************************************
  !****if* formfactors_A_main/get_Filename
  ! NAME
  ! function get_Filename(this)
  !
  ! PURPOSE
  ! helper routine to construct the file name to be read in
  !****************************************************************************
  function get_Filename(this)

    use inputGeneral, only: path_to_input

    class(tFormFacArr), intent(inOut) :: this
    character(len=1000) :: get_Filename

    if (this%initFlag) call this%init()

    select case (which_MAIDVersion)
    case (MAID2003)
      get_Filename = trim(path_to_Input) // '/electronNucleon/' // trim(this%sShort) // '_03_21_91_19.out.bz2'
   case (MAID2007)
      if (.not.fineGridQ2) then
         get_Filename = trim(path_to_Input) // '/electronNucleon/' // trim(this%sShort) // '_07_51_93_19.out.bz2'
      else
         get_Filename = trim(path_to_Input) // '/electronNucleon/' // trim(this%sShort) // '_07_101_93_19.out.bz2'
      end if
    case default
       call Traceback()
    end select

  end function get_Filename

  !****************************************************************************
  !****f* formfactors_A_main/FormFacArr_getMinW
  ! NAME
  ! real function FormFacArr_getMinW(this)
  !
  ! PURPOSE
  ! return the minimal W value used for the tabulation
  !
  ! It is the lowest data point - dW/2, even if interpolation is used.
  ! If arrays are not allocated, some large value is returned
  !****************************************************************************
  pure real function FormFacArr_getMinW(this)
    class(tFormFacArr), intent(in) :: this

    if (allocated(this%W)) then
       FormFacArr_getMinW = this%W(1) - this%dW/2
    else
       FormFacArr_getMinW = 1e20 ! some large number
    end if
  end function FormFacArr_getMinW

  !****************************************************************************
  !****is* formfactors_A_main/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! read in the corresponding namelist.
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* formfactors_A_main/formfactors_pion
    ! NAME
    ! NAMELIST formfactors_pion
    ! PURPOSE
    ! This namelist includes the following switches:
    ! * which_MaidVersion
    ! * fineGridQ2
    !**************************************************************************
    NAMELIST /formfactors_pion/ which_MaidVersion, fineGridQ2

    integer :: ios

    if (.not.initFlag) return

    call Write_ReadingInput("formfactors_pion",0)
    rewind(5)
    read(5,nml=formfactors_pion,IOSTAT=ios)
    call Write_ReadingInput("formfactors_pion",0,ios)
    select case (which_MaidVersion)
    case (MAID2003)
       write(*,*) 'MAID 2003 is being used for gamma* N -> N pi'
    case (MAID2007)
       write(*,*) 'MAID 2007 is being used for gamma* N -> N pi'
    case default
       write(*,*) 'Wrong MAID switch!! which_MAIDVersion', which_MaidVersion
       write(*,*) 'Must be :'
       write(*,*) MAID2003, 'for MAID 2003 or ', MAID2007,' for MAID2007'
       call Traceback()
    end select
    write(*,*) 'fineGridQ2 = ', fineGridQ2
    call Write_ReadingInput("formfactors_pion",1)

    FormFacArr_pi0p = tFormFacArr('pi0p', 'gamma p -> p pi0')
    FormFacArr_pi0n = tFormFacArr('pi0n', 'gamma n -> n pi0')
    FormFacArr_pipn = tFormFacArr('pipn', 'gamma p -> n pi+')
    FormFacArr_pimp = tFormFacArr('pimp', 'gamma n -> p pi-')

    initFlag = .false.

  end subroutine readInput

  !============================================================================
  !============================================================================
  !============================================================================


  !****************************************************************************
  !****f* formfactors_A_main/getA
  ! NAME
  ! function getA(pionCharge_out,nucCharge_out,thetaIn,sIn,Q2In)
  !
  ! PURPOSE
  ! This function returns the invariant Amplitudes A_1, ... A_6 of the MAID
  ! analysis.
  !
  ! INPUTS
  ! * real :: thetaIn  -- Theta scattering angle of the pion relative to q
  !   in the CM system of the hadronic vertex
  ! * real :: Q2In -- Q^2=-q^mu q_mu  of the gamma
  ! * real :: sIn -- Mandelstam s of the hadronic vertex: gamma* N-> pi N
  ! * integer :: pionCharge_out,nucCharge_out
  !   -- charges of outgoing pion and nucleon
  !
  ! OUTPUT
  ! complex, dimension(1:6) -- The complex amplitudes A_1 to A_6 in units of:
  ! * [GeV**-2] for A(1)
  ! * [GeV**-4] for A(2) and A(5)
  ! * [GeV**-3] for A(3), A(4) and A(5)
  !
  !****************************************************************************
  function getA(pionCharge_out,nucCharge_out,thetaIn,sIn,Q2In)

    complex, dimension(1:6) :: getA
    real, intent(in) :: thetaIn
    real, intent(in) :: Q2In
    real, intent(in) :: sIn
    integer, intent(in) :: pionCharge_out,nucCharge_out

    real :: Q2
    getA = 0.

    if (Q2In.le.0) then
       if (Q2In.lt.-1E-4) then
          write(*,*) "Warning Q^2 less than zero in getA!", Q2In
       end if
       Q2 = 1E-12
    else
       Q2 = Q2In
    end if

    select case (pionCharge_out)
    case (1)
       select case (nucCharge_out)
       case (0)
          getA = FormFacArr_pipn%getA(thetaIn,sIn,Q2)
       case default
          call error()
       end select
    case (-1)
       select case (nucCharge_out)
       case (1)
          getA = FormFacArr_pimp%getA(thetaIn,sIn,Q2)
       case default
          call error()
       end select
    case (0)
       select case (nucCharge_out)
       case (0)
          getA = FormFacArr_pi0n%getA(thetaIn,sIn,Q2)
       case (1)
          getA = FormFacArr_pi0p%getA(thetaIn,sIn,Q2)
       case default
          call error()
       end select
    case default
       call error()
    end select

  contains

    subroutine error()
      write(*,'(A)') 'Invalid charges in formfactors_A_main/getA'
      write(*,'(A,I4)') 'Outgoing pion charge   =', pionCharge_out
      write(*,'(A,I4)') 'Outgoing nucleon charge=', nucCharge_out
      call Traceback('Severe problem: Stopping')
    end subroutine error

  end function getA

  !****************************************************************************
  !****s* formfactors_A_main/cleanUp
  ! NAME
  ! subroutine cleanUp
  !
  ! PURPOSE
  ! dealloctae all possible formfactor arrays
  !****************************************************************************
  subroutine cleanUp
    call FormFacArr_pi0p%cleanUp
    call FormFacArr_pi0n%cleanUp
    call FormFacArr_pipn%cleanUp
    call FormFacArr_pimp%cleanUp
  end subroutine

  !****************************************************************************
  !****f* formfactors_A_main/getMinW
  ! NAME
  ! real function getMinW
  !
  ! PURPOSE
  ! return the minimal W value used in the tables
  !****************************************************************************
  real function getMinW()
    real :: minW
    minW = 1e20
    minW = min(minW,FormFacArr_pi0p%getMinW())
    minW = min(minW,FormFacArr_pi0n%getMinW())
    minW = min(minW,FormFacArr_pipn%getMinW())
    minW = min(minW,FormFacArr_pimp%getMinW())

    getMinW = minW

  end function getMinW



end module formfactors_A_main
