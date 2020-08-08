
C     **************************************************
C     *  getGrowthTime                                 *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calculates the time it takes for a particle to grow
C     from one size to the next by condensation of sulfuric acid (and
C     associated NH3 and water) onto particles.

C     This subroutine assumes that the growth happens entirely in the kinetic
C     regine such that the dDp/dt is not size dependent.  The time for growth 
C     to the first size bin may then be approximated by the time for growth via
C     sulfuric acid (not including nh4 and water) to the size of the first size bin
C     (not including nh4 and water).

C-----INPUTS------------------------------------------------------------

C     d1: intial diameter [m]
C     d2: final diameter [m]
c     h2so4: h2so4 ammount [kg]
c     temp: temperature [K]
c     boxvol: box volume [cm3]

C-----OUTPUTS-----------------------------------------------------------

C     gtime: the time it takes the particle to grow to first size bin [s]

      SUBROUTINE getGrowthTime(d1,d2,h2so4,temp,boxvol,gtime)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision d1,d2    ! initial and final diameters [m]
      double precision h2so4    ! h2so4 ammount [kg]
      double precision temp     ! temperature [K]
      double precision boxvol   ! box volume [cm3]
      double precision gtime    ! the time it will take the particle to grow 
                                ! to first size bin [s]

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision density  ! density of particles in first bin [kg/m3]
      double precision pi, R, MW
      double precision csulf    ! concentration of sulf acid [kmol/m3]
      double precision mspeed   ! mean speed of molecules [m/s]
      double precision alpha    ! accomidation coef

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654d0, R=8.314d0) !pi and gas constant (J/mol K)
      parameter(density=1800.d0,MW=98.d0) ! density [kg/m3], mol wgt sulf [kg/kmol]
      parameter(alpha=0.65)

C-----CODE--------------------------------------------------------------

      csulf = h2so4/MW/(boxvol*1d-6) ! SA conc. [kmol/m3]
      mspeed = sqrt(8.d0*R*temp*1000.d0/(pi*MW))

C     Kinetic regime expression (S&P 11.25) solved for T
      gtime = (d2-d1)/(4.d0*MW/density*mspeed*alpha*csulf)

      RETURN
      END
