!******************************************************************************
!****m* /propagation_RMF
! NAME
! module propagation_RMF
! PURPOSE
! Module which includes the propagation of the test-particles
! in relativistic mean field. It substitutes the module propagation
! in case if RMF is used.
!******************************************************************************
module propagation_RMF

  implicit none

  private
  !****************************************************************************
  !****n* propagation_RMF/propagation_RMF_input
  ! NAME
  ! NAMELIST propagation_RMF_input
  ! PURPOSE
  ! Namelist which includes the input switches:
  ! * predictorCorrector
  !****************************************************************************

  !****************************************************************************
  !****g* propagation_RMF/predictorCorrector
  ! PURPOSE
  ! Switch for predictor-corrector method in the propagation.
  ! If .false. then simple Euler method is used
  ! (i.e. only predictor step is done)
  ! SOURCE
  !
  logical, save :: predictorCorrector=.true.
  !****************************************************************************

  logical, save :: startFlag=.true.
  logical, parameter :: debugFlag=.false.

  public :: propagate_RMF ! subroutine to propagate particles

contains


  !****************************************************************************
  subroutine init
    use output

    integer :: ios
    NAMELIST /propagation_RMF_input/ predictorCorrector

    call Write_ReadingInput('propagation_RMF_input',0)
    rewind(5)
    read(5,nml=propagation_RMF_input,iostat=ios)
    call Write_ReadingInput('propagation_RMF_input',0,ios)

    write(*,*) ' Set Predictor-Corrector mode to',predictorCorrector,'.'

    call Write_ReadingInput('propagation_RMF_input',1)
  end subroutine init

  !****************************************************************************
  !****s* propagation_RMF/propagate_RMF
  ! NAME
  ! subroutine propagate_RMF (realTeilchen, pertTeilchen, delta_T, TimeStep)
  ! PURPOSE
  ! This routine propagates particles in a relativistic mean field.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realTeilchen
  ! * type(particle), dimension(:,:) :: pertTeilchen
  !   -- real and perturbative particle arrays which should be propagated
  ! * real :: delta_T -- size of time step (fm/c)
  ! * integer :: TimeStep -- number of time step
  ! NOTES
  ! Particle coordinates and kinetic momentum are propagated.
  ! Uses predictor corrector scheme or simple Euler time stepping.
  !****************************************************************************
  subroutine propagate_RMF (realTeilchen, pertTeilchen, delta_T, TimeStep)

    use particleDefinition
    use inputGeneral, only: freezeRealParticles
    use densitymodule, only: UpdateRMF, energyDeterminationRMF, storeFields
    use coulomb, only: UpdateCoulomb
    use offShellPotential, only: treatParticleOffShell,SetOffShellEnergy
    use output, only: WriteParticle
    use propagation, only: checkVelos
    
    type(particle),intent(inOut),dimension(:,:) :: realTeilchen
    type(particle),intent(inOut),dimension(:,:) :: pertTeilchen
    real, intent(in) :: delta_T
    integer, intent(in) :: TimeStep

    ! Working variables:
    real, save, Allocatable, dimension(:,:,:) :: rk_re, pk_re, rk_pe, pk_pe    ! arrays to store the coordinates and momenta
                                                                               ! at the current time
    real, save, Allocatable, dimension(:,:,:) :: pdot_old_re, pdot_old_pe      ! arrays to store dp^*/dt calculated
                                                                               ! at the current time
    real, save, Allocatable, dimension(:,:,:) :: rdot_old_re, rdot_old_pe    ! arrays to store the particle velocities
                                                                                     ! at the current time

    integer, save :: n_ensembles, n_particles_re, n_particles_pe
    integer :: i, j
    real, dimension(1:3) :: pdot, rdot
    real, save :: delta_T_old  ! Previous time step size

    if (startFlag) then !Initialisierung im ersten Aufruf
       call init
       n_ensembles = Size(realTeilchen,dim=1)
       n_particles_re = Size(realTeilchen,dim=2)
       n_particles_pe = Size(pertTeilchen,dim=2)
       allocate(rk_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(rk_pe(1:n_ensembles,1:n_particles_pe,1:3))
       allocate(pk_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(pk_pe(1:n_ensembles,1:n_particles_pe,1:3))
       allocate(pdot_old_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(pdot_old_pe(1:n_ensembles,1:n_particles_pe,1:3))
       allocate(rdot_old_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(rdot_old_pe(1:n_ensembles,1:n_particles_pe,1:3))
       startFlag=.false.
    end if

    if ( TimeStep == 1 ) delta_T_old = delta_T

    ! Predictor Step:

    Loop_over_ensembles_1 : do i=1,n_ensembles

      if (.not.freezeRealParticles) then
          do j=1,n_particles_re

              call SetOffShellEnergy(realTeilchen(i,j))
             
              if ( realTeilchen(i,j)%id == 0 ) then
                 cycle
              else if ( realTeilchen(i,j)%id < 0 ) then
                 exit
              end if

              call rhs(realTeilchen(i,j),delta_T_old,pdot,rdot)

              ! if(realTeilchen(i,j)%id == 101) write(*,*)' pion velo: ', rdot
              
              if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
                   write(*,*)' superluminous real particle in RMF-predictor step'
                   call WriteParticle(6,i,j,realTeilchen(i,j))
                   write(*,*)' offshellparameter :', realTeilchen(i,j)%offshellparameter                  
                   rdot(1:3)=realTeilchen(i,j)%momentum(1:3)/realTeilchen(i,j)%momentum(0)
                   pdot(1:3)=0.
              end if                

              rk_re(i,j,1:3) = realTeilchen(i,j)%position(1:3)
              pk_re(i,j,1:3) = realTeilchen(i,j)%momentum(1:3)
              
              realTeilchen(i,j)%position(1:3) = realTeilchen(i,j)%position(1:3)+delta_T*rdot(1:3)

              realTeilchen(i,j)%momentum(1:3) = realTeilchen(i,j)%momentum(1:3)+delta_T*pdot(1:3)

              call SetOffShellEnergy(realTeilchen(i,j))
              realTeilchen(i,j)%velocity(1:3)=rdot(1:3)
             
              pdot_old_re(i,j,1:3) = pdot(1:3)
              rdot_old_re(i,j,1:3) = rdot(1:3)

          end do
      end if

      do j=1,n_particles_pe

          call SetOffShellEnergy(pertTeilchen(i,j))
         
          if ( pertTeilchen(i,j)%id == 0 ) then
             cycle
          else if ( pertTeilchen(i,j)%id < 0 ) then
             exit
          end if

          call rhs(pertTeilchen(i,j),delta_T_old,pdot,rdot)

          if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
               write(*,*)' superluminous perturbative particle in RMF-predictor step'
               call WriteParticle(6,i,j,pertTeilchen(i,j))
               write(*,*)' offshellparameter :', pertTeilchen(i,j)%offshellparameter                  
               rdot(1:3)=pertTeilchen(i,j)%momentum(1:3)/pertTeilchen(i,j)%momentum(0)
               pdot(1:3)=0.
          end if

          rk_pe(i,j,1:3) = pertTeilchen(i,j)%position(1:3)
          pk_pe(i,j,1:3) = pertTeilchen(i,j)%momentum(1:3)
            
          pertTeilchen(i,j)%position(1:3) = pertTeilchen(i,j)%position(1:3)+delta_T*rdot(1:3)

          pertTeilchen(i,j)%momentum(1:3) = pertTeilchen(i,j)%momentum(1:3)+delta_T*pdot(1:3)

          call SetOffShellEnergy(pertTeilchen(i,j))
          pertTeilchen(i,j)%velocity(1:3)=rdot(1:3) 
          
          pdot_old_pe(i,j,1:3) = pdot(1:3)
          rdot_old_pe(i,j,1:3) = rdot(1:3)
          
      end do

    end do Loop_over_ensembles_1

    
    ! Store the space components of the omega- and rho-field:
    if (.not.freezeRealParticles) call storeFields

    ! Corrector step (optional):

    if (predictorCorrector) then

       if (.not.freezeRealParticles) then
           call updateRMF(realTeilchen)
           call updateCoulomb
       end if

       Loop_over_ensembles_2 : do i=1,n_ensembles

          if (.not.freezeRealParticles) then
              do j=1,n_particles_re

                  if ( realTeilchen(i,j)%id == 0 ) then
                     cycle
                  else if ( realTeilchen(i,j)%id < 0 ) then
                     exit
                  end if

                  call rhs(realTeilchen(i,j),delta_T,pdot,rdot)

                  rdot(1:3)=0.5*(rdot_old_re(i,j,1:3)+rdot(1:3))
                  pdot(1:3)=0.5*(pdot_old_re(i,j,1:3)+pdot(1:3))
                  
                  if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
                     write(*,*)' superluminous real particle in RMF-corrector step'
                     call WriteParticle(6,i,j,realTeilchen(i,j))
                     write(*,*)' offshellparameter :', realTeilchen(i,j)%offshellparameter
                     rdot(1:3)=realTeilchen(i,j)%momentum(1:3)/realTeilchen(i,j)%momentum(0)
                     pdot(1:3)=0.
                  end if
                  
                  realTeilchen(i,j)%position(1:3) = rk_re(i,j,1:3)+delta_T*rdot(1:3)

                  realTeilchen(i,j)%momentum(1:3) = pk_re(i,j,1:3)+delta_T*pdot(1:3)

                  call SetOffShellEnergy(realTeilchen(i,j))
                  realTeilchen(i,j)%velocity(1:3)=rdot(1:3)
                  
              end do
          end if

          do j=1,n_particles_pe

              if ( pertTeilchen(i,j)%id == 0 ) then
                 cycle
              else if ( pertTeilchen(i,j)%id < 0 ) then
                 exit
              end if

              if(.not.treatParticleOffShell(pertTeilchen(i,j)%id,pertTeilchen(i,j)%OffShellParameter)) &
              & call energyDeterminationRMF( pertTeilchen(i,j) )
              
              call rhs(pertTeilchen(i,j),delta_T,pdot,rdot)

              rdot(1:3)=0.5*(rdot_old_pe(i,j,1:3)+rdot(1:3))
              pdot(1:3)=0.5*(pdot_old_pe(i,j,1:3)+pdot(1:3))

              if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
                  write(*,*)' superluminous perturbative particle in RMF-corrector step'
                  call WriteParticle(6,i,j,pertTeilchen(i,j))
                  write(*,*)' offshellparameter :', pertTeilchen(i,j)%offshellparameter
                  rdot(1:3)=pertTeilchen(i,j)%momentum(1:3)/pertTeilchen(i,j)%momentum(0)
                  pdot(1:3)=0.
              end if             

              pertTeilchen(i,j)%position(1:3) = rk_pe(i,j,1:3)+delta_T*rdot(1:3)

              pertTeilchen(i,j)%momentum(1:3) = pk_pe(i,j,1:3)+delta_T*pdot(1:3)

              call SetOffShellEnergy(pertTeilchen(i,j))
              pertTeilchen(i,j)%velocity(1:3)=rdot(1:3)
              
          end do

       end do Loop_over_ensembles_2

    end if

    if (.not.freezeRealParticles) then
        call updateRMF(realTeilchen)
        call updateCoulomb
    end if

    Loop_over_ensembles_3 : do i=1,n_ensembles
        do j=1,n_particles_pe

           if ( pertTeilchen(i,j)%id == 0 ) then
               cycle
           else if ( pertTeilchen(i,j)%id < 0 ) then
               exit
           end if

           if(treatParticleOffShell(pertTeilchen(i,j)%id,pertTeilchen(i,j)%OffShellParameter)) cycle
           
           call energyDeterminationRMF( pertTeilchen(i,j) )
           pertTeilchen(i,j)%velocity(1:3) = pertTeilchen(i,j)%momentum(1:3) &
                                        & / pertTeilchen(i,j)%momentum(0)
        end do
    end do Loop_over_ensembles_3

    delta_T_old = delta_T

    call checkVelos(realTeilchen)
    call checkVelos(pertTeilchen)

  end subroutine propagate_RMF



  subroutine rhs(teilchen,dt,pdot,rdot)

    use particleDefinition
    use IdTable
    use densitymodule
    use RMF
    use coulomb, only: emfoca
    use propagation, only: gradients
    use eventtypes, only: HeavyIon
    use inputGeneral, only: eventtype
    
    real,                 intent(in) :: dt
    type(particle),       intent(in) :: teilchen
    real, dimension(1:3), intent(out) :: pdot, rdot

    real, dimension(0:3) :: V             ! vector component of the particle self-energy (vector field) [GeV]

    real, dimension(0:3,1:3) :: dV_dr     ! space derivatives of the vector field V (GeV/fm),
                                          ! 1-st index -> field component,
                                          ! 2-nd index -> space index of the derivative

    real,  dimension(1:3) :: dsigma_dr ! space derivatives of the scalar field (GeV/fm)

    real, dimension(1:3) :: Grad_P     ! auxiliary, needed for off-shell particles
    real, dimension(0:3) :: Grad_R     ! auxiliary, needed for off-shell particles
    
    integer :: i1, k, k_max,id,charge
    integer :: index1, index2, index3
    real, dimension(1:3) :: pdot_space, pdot_time, place, impuls, emForce
    real :: cpot, factor, mass
    real :: V1, V2, S1, S2
    logical :: anti, AddCoulomb


    if (lorentz_flag) then
      k_max=3    ! With space components of the omega-field
    else
      k_max=0    ! W/o space components of the omega-field
    end if

    id = teilchen%id
    charge = teilchen%charge
    anti = teilchen%antiparticle
    mass = teilchen%mass

    ! Modification factor for the coupling constants:
    factor = ModificationFactor(Id,anti)

    pdot_space(1:3) = 0.
    pdot_time(1:3) = 0.
    rdot = teilchen%velocity
    place(1:3) = teilchen%position(1:3)
    impuls(1:3) = teilchen%momentum(1:3)

    AddCoulomb = .true.
    
    ! nuclear part (only for baryons, kaons and antikaons):
    if ( (isBaryon(id)  .or. id.eq.Kaon .or. id.eq.kaonBar) .and. factor.gt.0. ) then

        index1 = nint( place(1) / gridSpacing(1) )
        index2 = nint( place(2) / gridSpacing(2) )
        index3 = nint( place(3) / gridSpacing(3) )

        if ( abs(index1) <= gridPoints(1)-1 .and. &
           &abs(index2) <= gridPoints(2)-1 .and. &
           &abs(index3) <= gridPoints(3)-1 ) then     ! Inside grid:

           do k = 0,k_max

              V(k) = Selfenergy_vector(Index1,Index2,Index3,k,id,charge,anti)

              V1 = Selfenergy_vector(Index1+1,Index2,Index3,k,id,charge,anti)
              V2 = Selfenergy_vector(Index1-1,Index2,Index3,k,id,charge,anti)

              dV_dr(k,1) = ( V1 - V2) / (2.*gridSpacing(1))

              V1 = Selfenergy_vector(Index1,Index2+1,Index3,k,id,charge,anti)
              V2 = Selfenergy_vector(Index1,Index2-1,Index3,k,id,charge,anti)

              dV_dr(k,2) = ( V1 - V2) / (2.*gridSpacing(2))

              V1 = Selfenergy_vector(Index1,Index2,Index3+1,k,id,charge,anti)
              V2 = Selfenergy_vector(Index1,Index2,Index3-1,k,id,charge,anti)

              dV_dr(k,3) = ( V1 - V2) / (2.*gridSpacing(3))

           end do

           !isoscalar, scalar contributions (dm*/dr_i):

           S1 = Selfenergy_scalar(Index1+1,Index2,Index3,id,anti)
           S2 = Selfenergy_scalar(Index1-1,Index2,Index3,id,anti)

           dsigma_dr(1) = ( S1 - S2 ) / (2.*gridSpacing(1))

           S1 = Selfenergy_scalar(Index1,Index2+1,Index3,id,anti)
           S2 = Selfenergy_scalar(Index1,Index2-1,Index3,id,anti)

           dsigma_dr(2) = ( S1 - S2 ) / (2.*gridSpacing(2))

           S1 = Selfenergy_scalar(Index1,Index2,Index3+1,id,anti)
           S2 = Selfenergy_scalar(Index1,Index2,Index3-1,id,anti)

           dsigma_dr(3) = ( S1 - S2 ) / (2.*gridSpacing(3))

           if ( kaonpot_flag .and. (teilchen%id == kaon .or. teilchen%id == kaonBar) ) then
              do i1=1,3
                 dsigma_dr(i1) = -dsigma_dr(i1) &
                      &          + ( dot_product( V(0:k_max),dV_dr(0:k_max,i1) ) ) / (2.*gridSpacing(i1))
              end do
              dsigma_dr(1:3) = dsigma_dr(1:3) / DiracMass(Index1,Index2,Index3,mass,id,charge,anti)
           end if


           do i1 = 1,3

             !**** Space derivative contributions
             !**** to the rhs of equation dp^*/dt=pdot
             !**** pdot_space(i1) = -domega_dr(0,i1) -drho_dr(0,i1)

             !**** vector field contributions:
             pdot_space(i1) = -dV_dr(0,i1)

             if (lorentz_flag) then
               pdot_space(i1) = pdot_space(i1) &
                    & + dot_product( teilchen%velocity(1:3), (dV_dr(1:3,i1)-dV_dr(i1,1:3)) )
             end if

             ! For antibaryons or antikaons the vector field changes sign:
             if ( anti .or. id.eq.kaonBar ) pdot_space(i1) = -pdot_space(i1)

             !**** Add-up scalar field contributions:
             pdot_space(i1) = pdot_space(i1) -(DiracMass(Index1,Index2,Index3,mass,id,charge,anti)) &
                                           &  / teilchen%momentum(0) *  dsigma_dr(i1)

             !**** Time derivative contribution
             !**** to the rhs of equation dp^*/dt=pdot

             if (lorentz_flag) then
                pdot_time(i1) = ( &
                     &   SelfEnergy_vector_old(Index1,Index2,Index3,i1,id,charge,anti) &
                     & - SelfEnergy_vector(Index1,Index2,Index3,i1,id,charge,anti)  &
                     &           ) / dt
                ! For antibaryons or antikaons the vector field changes sign:
                if ( anti .or. id.eq.kaonBar ) pdot_time(i1) = -pdot_time(i1)
             else
                pdot_time(i1) = 0.
             end if

           end do

        end if

        rdot = teilchen%velocity
        
        AddCoulomb = .true.
        
    else if(eventtype==HeavyIon .and. (id.eq.rho .or. id.eq.omegaMeson .or. id.eq.phi)) then

        call gradients(teilchen,Grad_P,Grad_R)

        rdot = Grad_P        
        pdot_space(1:3) = -Grad_R(1:3)

        AddCoulomb = .false.
       
    end if

    if(AddCoulomb) then
       !**** electromagnetic part
       cpot = emfoca(place,impuls,teilchen%charge,teilchen%ID,emForce)
       pdot_space(1:3) = pdot_space(1:3) + emForce(1:3)
    end if
    
    !**** Total derivative dp^*/dt
    pdot = pdot_space + pdot_time
    
  end subroutine rhs


end module propagation_RMF
