
C     **************************************************
C     *  getH2SO4prod                                  *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     Calculates the rate of SO2 + OH -> H2SO4 as in Koch's sulfur
C     model, but for a simple box model.  Assumes typical OH conc.
C     for clear/cloudy and day/night conditions.

C     The output is the average rate of H2SO4 production [kg/s] during
C     the timestep.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE getH2SO4prod(Cso2,ohc,dt,temp,pres,H2SO4rate)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Cso2  !concentrations, arbitrary units
      double precision dt              !chemistry time step
      double precision ohc   !OH conc (molec cm-3)
      double precision temp, pres      !temperature (K) and pressure (Pa)

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision fconv !fraction of so2 converted to h2so4
      double precision rate  !first-order loss rate (s-1)
      double precision H2SO4rate !average rate of H2SO4 production during step (kg s-1)
      double precision Ch2so4f !final concentration, arbittrary units
      double precision RK4, EK4, DMM, TT
      double precision ft, tfac, expo, dfac

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      ft=4.0e-20
      tfac=300.
      expo=3.3
      dfac=1.e-11
      DMM=pres/1.e5/(.082*temp)*6.02E20
      TT=1./temp
      RK4 = ft*((TT*tfac)**(expo))*DMM*dfac
      EK4 = 1./(1. + ((LOG10(RK4/2.0E-12))**2.))
      rate = ohc * (RK4/(1. + RK4/2.0E-12))*(0.45**EK4)
      fconv=exp(-1.*rate*dt)
c      write(*,*) 'boxchem: fconv=',fconv
      Ch2so4f=(1.-fconv)*Cso2*98./64.
      H2SO4rate=Ch2so4f/dt
      Cso2=Cso2*fconv

      RETURN
      END
