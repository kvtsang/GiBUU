program test
  use formfactors_A_main
  use particleProperties, only: initParticleProperties
  use inputGeneral, only: readinputGeneral
  implicit none

  real, parameter :: dW = 0.01, dTheta = 10.
  integer, parameter :: nW = 2, nTheta = 19
  real, dimension(nW*nTheta) :: arrW,arrTheta

  call readinputGeneral
  call initParticleProperties

  call test0
!  call test1

contains

  subroutine test0
    complex, dimension(1:6) :: A

    !  A = FormFacArr_pi0p%getAA( 90., 1.2, 0.005 )
    A = getA( 1, 0, 90., 1.44, 0.005 )
    write(*,*) A
    A = getA( 0, 1, 90., 1.44, 0.005 )
    write(*,*) A
    A = getA( 1, 1, 90., 1.44, 0.005 )
    write(*,*) A
  end subroutine test0

  subroutine test1
    real :: W,Theta
    integer :: iW, iTheta, ii

    ii = 1
    do iW=1,nW
       W = 1.1 + (iW-1)*dW
       do iTheta=1,nTheta
          Theta = 0. + (iTheta-1)*dTheta

          arrW(ii) = W
          arrTheta(ii) = Theta
          ii = ii+1
       end do
    end do

    do ii=1,nW*nTheta
       write(*,*) ii,arrW(ii),arrTheta(ii)
    end do

    ii = get_i(1.1, 2.)
    ii = get_i(1.1, 175.)
    ii = get_i(1.1, 180.)
    ii = get_i(1.1, 185.)

    W = 1.1
    do iTheta =-12,202,2
       Theta = iTheta*1.0
       ii = get_i(W,Theta)
       write(112,*) W,theta,ii,arrW(ii),arrTheta(ii)
    end do

    Theta = 5.
    do iW = 1070,1135,2
       W = iW/1000.
       ii = get_i(W,Theta)
       write(113,*) W,theta,ii,arrW(ii),arrTheta(ii)
    end do

  end subroutine test1

  integer function get_i(W,Theta)
    real, intent(in) :: W,Theta
    integer :: iW,iTheta,ii

    ! int --> lower

!!$    iTheta = int( (Theta-arrTheta(1))/dTheta )+1
!!$    iW = int( (W-arrW(1))/dW )
!!$    iTheta = min(max(1,iTheta),nTheta-1) ! attention: -1 !!!
!!$    iW = min(max(0,iW),nW-2) ! attention: -2 !!!

    ! nint --> closest

    iTheta = nint( (Theta-arrTheta(1))/dTheta )+1
    iW = nint( (W-arrW(1))/dW )
    iTheta = min(max(1,iTheta),nTheta)
    iW = min(max(0,iW),nW-1)


    ii = iW*nTheta + iTheta
    write(*,*) 'W,Theta=',W,Theta
    write(*,*) 'i::',iTheta,iW,'-->',ii
    write(*,*) '-->', arrW(ii),arrTheta(ii)
    get_i = ii

  end function get_i



end program test
