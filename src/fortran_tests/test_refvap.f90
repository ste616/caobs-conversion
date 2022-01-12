program test_refvap
  implicit none

  ! Set up a test case for an execution of refvap.
  COMPLEX refvap ! Function definition
  COMPLEX nvap ! Value coming from refvap

  ! Observing frequency in Hz, as a real
  REAL nu
  ! Atmospheric temperature in K, as a real
  REAL t
  ! Partial pressure of dry components in Pa, as a real
  REAL pdry
  ! Partial pressure of water vapour in Pa, as a real
  REAL pvap

  nu = 22000000000.0 ! 22 GHz
  t = 293.0 ! Roughly 20 Celcius
  pvap = 2317.0
  pdry = 101325.0 - pvap
  
  nvap = refvap(nu, t, pdry, pvap)
  print *, nu
  print *, t
  print *, pdry
  print *, pvap
  print *, nvap
  
end program test_refvap

