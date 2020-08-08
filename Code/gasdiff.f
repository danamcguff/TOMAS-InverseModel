
C     **************************************************
C     *  gasdiff                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, May 2000

C     This function returns the diffusion constant of a species in
C     air (m2/s).  It uses the method of Fuller, Schettler, and
C     Giddings as described in Perry's Handbook for Chemical
C     Engineers.

      double precision FUNCTION gasdiff(temp,pres,mw,Sv)

      IMPLICIT NONE

C     INCLUDE FILES...

C-----ARGUMENT DECLARATIONS------------------------------------------

      double precision temp, pres  !temperature (K) and pressure (Pa) of air
      double precision mw          !molecular weight (g/mol) of diffusing species
      double precision Sv          !sum of atomic diffusion volumes of diffusing species

C-----VARIABLE DECLARATIONS------------------------------------------

      double precision mwair, Svair   !same as above, but for air
      double precision mwf, Svf

C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------

      parameter(mwair=28.9, Svair=20.1)

C-----CODE-----------------------------------------------------------

      mwf=sqrt((mw+mwair)/(mw*mwair))
      Svf=(Sv**(1./3.)+Svair**(1./3.))**2.
      gasdiff=1.0e-7*temp**1.75*mwf/pres*1.0e5/Svf

      RETURN
      END
