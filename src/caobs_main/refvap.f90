
COMPLEX FUNCTION refvap(nu,t,pdry,pvap)

IMPLICIT NONE

REAL, INTENT(IN)                         :: nu
REAL, INTENT(IN)                         :: t
REAL, INTENT(IN)                         :: pdry
REAL, INTENT(IN)                         :: pvap


!  Determine the complex refractivity of the water vapour monomers.

!  Inputs:
!    nu  Observing frequency (Hz).
!    T  Temperature (Kelvin).
!    Pdry Partial pressure of dry components (Pascals).
!    Pvap Partial pressure of water vapour (Pascals).

!  Reference:
!    Liebe, An updated model for millimeter wave propogation in moist air,
!    Radio Science, 20, 1069-1089 (1985).
!------------------------------------------------------------------------
INTEGER              :: i
INTEGER, PARAMETER   :: nl = 30

REAL                 :: theta, p, e, f, b1(nl),b2(nl),b3(nl)
REAL (KIND=8)        :: nr, ni, x, y, z, gamma, s, nu0(nl)

!  Table of the microwave water lines.
DATA (nu0(i),b1(i),b2(i),b3(i),i=1,18)/  &
    22.235080D00,  0.1090, 2.143, 27.84E-3,  &
    67.813960D00,  0.0011, 8.730, 27.60E-3,  &
    119.995940D00,  0.0007, 8.347, 27.00E-3,  &
    183.310117D00,  2.3000, 0.653, 28.35E-3,  &
    321.225644D00,  0.0464, 6.156, 21.40E-3,  &
    325.152919D00,  1.5400, 1.515, 27.00E-3,  &
    336.187000D00,  0.0010, 9.802, 26.50E-3,  &
    380.197372D00, 11.9000, 1.018, 27.60E-3,  &
    390.134508D00,  0.0044, 7.318, 19.00E-3,  &
    437.346667D00,  0.0637, 5.015, 13.70E-3,  &
    439.150812D00,  0.9210, 3.561, 16.40E-3,  &
    443.018295D00,  0.1940, 5.015, 14.40E-3,  &
    448.001075D00, 10.6000, 1.370, 23.80E-3,  &
    470.888947D00,  0.3300, 3.561, 18.20E-3,  &
    474.689127D00,  1.2800, 2.342, 19.80E-3,  &
    488.491133D00,  0.2530, 2.814, 24.90E-3,  &
    503.568532D00,  0.0374, 6.693, 11.50E-3,  &
    504.482692D00,  0.0125, 6.693, 11.90E-3/

DATA (nu0(i),b1(i),b2(i),b3(i),i=19,30)/  &
    556.936002D00, 510.000, 0.114, 30.00E-3,  &
    620.700807D00,  5.0900, 2.150, 22.30E-3,  &
    658.006500D00,  0.2740, 7.767, 30.00E-3,  &
    752.033227D00, 250.000, 0.336, 28.60E-3,  &
    841.073593D00,  0.0130, 8.113, 14.10E-3,  &
    859.865000D00,  0.1330, 7.989, 28.60E-3,  &
    899.407000D00,  0.0550, 7.845, 28.60E-3,  &
    902.555000D00,  0.0380, 8.360, 26.40E-3,  &
    906.205524D00,  0.1830, 5.039, 23.40E-3,  &
    916.171582D00,  8.5600, 1.369, 25.30E-3,  &
    970.315022D00,  9.1600, 1.842, 24.00E-3,  &
    987.926764D00, 138.000, 0.178, 28.60E-3/

!  Convert to the units of Liebe.

theta = 300 /  t
e     = 0.001 * pvap
p     = 0.001 * pdry
f     = nu * 1E-9

nr = 2.39 * e * theta + 41.6 * e * theta * theta + 6.47E-6*f**2.05*e*theta**2.4
ni = (0.915*1.40E-6*p + 5.41E-5*e*theta*theta*theta)* f*e*theta**2.5

!  Sum the contributions of the lines.

DO i = 1, nl
   s     = b1(i) * e  * theta**3.5 * EXP (b2(i) * (1 - theta))
   gamma = b3(i) * (p * theta**0.8 + 4.80 * e * theta)
   x     = (nu0(i) - f) * (nu0(i) - f) + gamma * gamma
   y     = (nu0(i) + f) * (nu0(i) + f) + gamma * gamma
   z     = (nu0(i) + gamma * gamma / nu0(i))
   nr    = nr + s * ((z - f) / x + (z + f) / y - 2 / nu0(i))
   ni    = ni + s * ((1/x+1/y)*gamma*f/nu0(i))
END DO

!  Return the result.

refvap = CMPLX (REAL(nr), REAL(ni))

END FUNCTION refvap
