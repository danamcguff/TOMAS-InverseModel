
C     **************************************************
C     *  logemiss                                      *
C     **************************************************

C     WRITTEN BY Dana McGuffin, April 2017 FROM Peter Adams, November 1999

C     Initializes aerosol number and mass distributions with 
C     specified species in a lognormal size distribution.

C-----INPUTS------------------------------------------------------------

C     M - total mass of emissions (kg)
C     Dp - lognormal mean diameter (microns)
C     sigma - width parameter
C     bin1,bin2 - range of bins to add mass

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE logemiss(Ntot, MtotFF, MtotBB, Dpbar, sigma, bin1,bin2,
     &           density, idtcomp,FLG, Nkadd, Mkadd)!, fname)


      IMPLICIT NONE

C     INCLUDE FILES

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Ntot, Dpbar, sigma, density
      double precision MtotFF, MtotBB !total mass if FLG = 0
      integer idtcomp, bin1, bin2    !indicates chemical component to be added
      integer FLG  !1=lognormal parameters, 0=emission scaling factors
c                  !                          from Bond et al used in geos chem
      double precision Mkadd(ibins), Nkadd(ibins)
c      character*90 fname

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer b,n     !bin and species counters
      double precision Dl, Dh     !lower and upper bounds of bin 
      double precision Dk         !mean diameter of current bin 
      double precision pi
      double precision No
      double precision ScaleMassDistFF(ibins), ScaleNumDist(ibins)
      double precision ScaleMassDistBB(ibins)
      double precision normalize

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.14159)
      DATA ScaleMassDistFF/
     &   0.0 , 0.0 , 0.0 , 
     &   0.0     , 0.0     , 0.0     , 0.0     , 0.0     ,
     &   0.0     , 0.0     , 0.0     , 0.0     , 0.0     ,
     &   1.04E-03, 2.77E-03, 6.60E-03, 1.41E-02, 2.69E-02,
     &   4.60E-02, 7.06E-02, 9.69E-02, 1.19E-01, 1.31E-01,
     &   1.30E-01, 1.15E-01, 9.07E-02, 6.44E-02, 4.09E-02,
     &   2.33E-02, 1.19E-02, 5.42E-03, 2.22E-03, 8.12E-04,
     &   2.66E-04, 7.83E-05, 2.06E-05, 4.86E-06, 1.03E-06,
     &   1.94E-07, 3.29E-08, 4.99E-09, 6.79E-10, 8.26E-11/

      DATA ScaleMassDistBB/
     &   0.0 , 0.0 , 0.0 , 
     &   0.0     , 0.0     , 0.0     , 0.0     , 0.0     ,
     &   0.0     , 0.0     , 0.0     , 0.0     , 0.0     ,
     &   3.2224e-07, 1.6605e-06, 7.6565e-06, 3.1592e-05, 0.00011664,
     &   0.00038538, 0.0011394 , 0.0030144 , 0.0071362 , 0.015117  ,
     &   0.028657  , 0.048612  , 0.073789  , 0.10023   , 0.12182   ,
     &   0.1325    , 0.12895   , 0.11231   , 0.087525  , 0.061037  ,
     &   0.038089  , 0.02127   , 0.010628  , 0.0047523 , 0.0019015 ,
     &   0.00068081, 0.00021813, 6.2536e-05, 1.6044e-05, 3.6831e-05/

C-----CODE--------------------------------------------------------------
c      print*,'input Diam', Dpbar

      do n=1,ibins
        Nkadd(n)=0d0
        Mkadd(n)=0d0
      enddo

      IF (FLG .eq. 1) then
      !Loop over number of size bins
      do n=bin1,bin2
                  
         !Calculate diameter of this size bin
         Dl=1.0d6*((6.*xk(n))/(density*3.14))**0.3333 ! um
         Dh=1.0d6*((6.*xk(n+1))/(density*3.14))**0.3333 ! um
         Dk=sqrt(Dl*Dh) ! um
           !Calculate number concentration
           Nkadd(n)=Ntot/(sqrt(2.*pi)*Dk*log(sigma))
     &         *exp(-( (log(Dk/Dpbar))**2./(2.*(log(sigma))**2.) )) 
     &         *(Dh-Dl)
           !Calculate component mass
           !Assume spherical
           Mkadd(n)=Nkadd(n)*density*pi/6*(1.0e-6*Dk)**3


      enddo

      ELSE
      !! Use size distribution from GEOS-Chem for 40-bin TOMAS fossil
       !fuel & biomass OC emissions
c         normalize = 0.d0
c      do n=1,ibins
c        Dl=((6.*xk(n))/(density*3.14))**0.3333 ! [m]
c         Dh=((6.*xk(n+1))/(density*3.14))**0.3333 ! [m]
c         Dk=sqrt(Dl*Dh)
c         ScaleNumDist(n) = ScaleMassDist(n)/(Dk**3)
c         normalize = normalize + ScaleNumDist(n)
c      enddo
c      ScaleNumDist(:) = ScaleNumDist(:)/normalize
      do n=1,ibins
          Mkadd(n)=MtotFF*ScaleMassDistFF(n)+MtotBB*ScaleMassDistBB(n)
          Nkadd(n)=Mkadd(n)/sqrt(xk(n)*xk(n+1))
c         Nkadd(n)=Ntot*ScaleNumDist(n)
c         Mkadd(n)=Nkadd(n)*sqrt(xk(n)*xk(n+1))
      enddo

      ENDIF      

1     format( *(G0.4,:,",") )

c      open(unit=36, file=trim("anthro_emis_")//trim(fname)//trim(".txt")
c     &       ,status='old',access='append')
c      write(36,1) (Nkadd(n),n=1,ibins)
c      close(36)

      RETURN
      END
