program test
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  use version, only: PrintVersion
  use inMediumWidth
  use inMediumWidth_rho
  
  call PrintVersion

  call readinputGeneral
  call initParticleProperties

  !call tabulate_inMediumWidth_mesons
  !call tabulate_inMediumWidth_baryons

  call tabulate_GammaColl_rho
  
end program test
