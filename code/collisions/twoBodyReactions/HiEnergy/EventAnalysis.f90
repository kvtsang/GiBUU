!******************************************************************************
!****m* /EventAnalysis
! PURPOSE
! This module contains routines for determination of the jet structure
! of events generated by PYTHIA
!******************************************************************************
module EventAnalysis

  IMPLICIT NONE
  private

  public :: TwoJets

contains
  
  !****************************************************************************
  !****s* EventAnalysis/TwoJets
  ! logical function TwoJets()
  !
  ! PURPOSE
  ! * Determine whether an event contains two jets with large transverse momenta
  ! OUTPUT
  ! *  true   -- two jets are found
  ! *  false  -- not found
  !****************************************************************************
  logical function TwoJets()

      use constants, only: pi
    
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/
      
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      SAVE /PYDAT1/

      integer :: NJET
      real, dimension(1:2) :: kt1,kt2,qt,kt
      real :: phi1,phi2,dphi,x
      
      TwoJets=.false.
      
      MSTU(46)=4     ! scaled JADE distance y_ij,
                     ! can be also set in the namelist pythia (file DoCollTools.f90)
      
      call pyclus(NJET)

      if(NJET.ne.2) return

      kt1=(/p(N+1,1),p(N+1,2)/)
      kt2=(/p(N+2,1),p(N+2,2)/)
      
      phi1=atan2(kt1(2),kt1(1))
      phi2=atan2(kt2(2),kt2(1))      
      dphi=abs(phi2-phi1)
      
      if(abs(dphi-pi).gt.pi/9.) return
      
      qt=kt1+kt2

      x=(p(N+1,4)+p(N+1,3))/(p(N+1,4)+p(N+2,4)+p(N+1,3)+p(N+2,3))   ! LC (+) momentum fraction of the 1st jet

      kt=kt1-qt*x    ! Galilean transformation to the c.m.s of a two-jet system
                     ! (valid with high accuracy for transverse momenta)

      if(sqrt(kt(1)**2+kt(2)**2).lt.1.5) return



      
      TwoJets=.true.

  end function TwoJets

end module EventAnalysis