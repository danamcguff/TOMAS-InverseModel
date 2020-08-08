
C     **************************************************
C     *  getCoagLoss                                   *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calculates the first order loss rate of particles
C     of a given size with respect to coagulation.

C-----INPUTS------------------------------------------------------------

C     d1: diameter in [m] of particle

C-----OUTPUTS-----------------------------------------------------------

C     ltc: first order loss rate with respect to coagulation [s-1]
C     sinkfrac: fraction of sink from each bin

      SUBROUTINE getCoagLoss(d1,ltc,Nko,Mko,sinkfrac)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision d1       ! diameter of the particle [m]
      double precision ltc    ! first order loss rate [s-1]
      double precision Nko(ibins), Mko(ibins, icomp)
      double precision sinkfrac(ibins)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,k,j
      double precision density  ! density of particles [kg/m3]
      double precision density1 ! density of particles in first bin
      double precision pi, R, MW, kB
      double precision mu       !viscosity of air (kg/m s)
      double precision mfp      !mean free path of air molecule (m)
      double precision Kn       !Knudsen number of particle
      double precision Mktot    !total mass of aerosol
      double precision kij(ibins)
      double precision Dpk(ibins) !diameter (m) of particles in bin k
      double precision Dk(ibins),Dk1 !Diffusivity (m2/s) of bin k particles
      double precision ck(ibins),ck1 !Mean velocity (m/2) of bin k particles
      double precision neps
      double precision meps
      double precision mp       ! mass of the particle (kg)
      double precision beta     !correction for coagulation coeff.
      double precision kij_self !coagulation coefficient for self coagulation


C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens
      external aerodens

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654d0, R=8.314d0, kB= 1.38E-23) !pi and gas constant (J/mol K)
      parameter (neps=1E8, meps=1E-8)

C-----CODE--------------------------------------------------------------

      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6

C Calculate particle sizes and diffusivities
      do k=1,ibins
         Mktot = 0.d0
         do j=1, icomp
            Mktot=Mktot+Mko(k,j)
         enddo
         Mktot=Mktot+2.d0*Mko(k,srtso4)/96.d0-Mko(k,srtnh4)/18.d0 ! account for h+

         if (Mktot.gt.meps)then
            density=aerodens(Mko(k,srtso4),0.d0,Mko(k,srtnh4),
     &           Mko(k,srtna),Mko(k,srth2o),Mko(k,srtc) ) !assume bisulfate
         else
            density = 1400.
         endif
         if(Nko(k).gt.neps .and. Mktot.gt.meps)then
            mp=Mktot/Nko(k)
         else
            mp=sqrt(xk(k)*xk(k+1))
         endif
         if (k.eq.1) density1 = density
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
         ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
      enddo
      
      Kn = 2.0*mfp/d1
      Dk1 = kB*temp/(3.0*pi*mu*d1)             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
      mp = 4.d0/3.d0*pi*(d1/2.d0)**3.d0*density1
      ck1=sqrt(8.0*kB*temp/(pi*mp)) !S&P Table 12.1      

      do i=1,ibins
         Kn=4.0*(Dk(i)+Dk1)          
     &        /(sqrt(ck(i)**2+ck1**2)*(Dpk(i)+d1)) !S&P eqn 12.51
         beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
                                !This is S&P eqn 12.46 with non-continuum correction, beta
         kij(i)=2.0*pi*(Dpk(i)+d1)*(Dk(i)+Dk1)*beta
         kij(i)=kij(i)*1.0e6/boxvol !normalize by grid cell volume
c         kij(i)=kij(i)*1.0e6 !cm3/s
      enddo
Cjrp
Cjrp      !self coagulation
Cjrp      Kn=4.0*(Dk1+Dk1)          
Cjrp     &     /(sqrt(ck1**2+ck1**2)*(d1+d1)) !S&P eqn 12.51
Cjrp      beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
Cjrp                                !This is S&P eqn 12.46 with non-continuum correction, beta
Cjrp      kij_self=2.0*pi*(d1+d1)*(Dk1+Dk1)*beta
Cjrp      print*,'kij_self',kij_self*1E6
Cjrp      kij_self=kij_self*1.0e6/boxvol !normalize by grid cell volume      
         
      ltc = 0.d0
c      print*,'ltc',ltc
      do i=1,ibins
         ltc = ltc + kij(i)*Nko(i)
      enddo
      do i=1,ibins
         sinkfrac(i) = kij(i)*Nko(i)/ltc
      enddo
c      print*,'ltc',ltc
      


      RETURN
      END
