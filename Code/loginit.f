
C     **************************************************
C     *  loginit                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     Initializes aerosol number and mass distributions with sulfate
C     only and a lognormal size distribution.

C-----INPUTS------------------------------------------------------------

C     No - total number concentration (#/cm3)
C     Dp - lognormal mean diameter (microns)
C     sigma - width parameter

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE loginit(N1, Dp, sigma, density, idtcomp)

C     INCLUDE FILES

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision N1, Dp, sigma, density
      integer idtcomp    !indicates chemical component to be added

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer b,n     !bin and species counters
      double precision Dl, Dh     !lower and upper bounds of bin (microns)
      double precision Dk         !mean diameter of current bin (microns)
      double precision np         !number of particles in size bin in GCM cell
      double precision pi
      double precision No
      
C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.14159)

C-----CODE--------------------------------------------------------------

C Convert No from #/cm3 to #/box
      No=N1*boxvol

C Loop over grid cells

      do l=1,lm
      do j=1,jm
      do i=1,im

      !Loop over number of size bins   !usually n=14,nbins
      do n=14,NBINS
         
         !Calculate diameter of this size bin
         Dl=1.0d6*((6.*xk(n))/(density*3.14))**0.3333
         Dh=1.0d6*((6.*xk(n+1))/(density*3.14))**0.3333
         Dk=sqrt(Dl*Dh)
c         write(*,*) 'Dk(',n,')= ',Dk,' microns'

         !Calculate number concentration
         np=No/(sqrt(2.*pi)*Dk*log(sigma))*exp(-( (log(Dk/Dp))**2./
     &         (2.*(log(sigma))**2.) )) * (Dh-Dl)
         T0M(i,j,l,IDTNUMD+n-1)=T0M(i,j,l,IDTNUMD+n-1)+np
c         write(*,*)'N=',T0M(i,j,l,IDTNUMD+n-1)

         !Calculate component mass
         T0M(i,j,l,idtcomp+n-1)=T0M(i,j,l,idtcomp+n-1)
     &                         +np*sqrt(xk(n))*sqrt(xk(n+1))
c         write(*,*)'M=',T0M(i,j,l,IDTSO4+n-1)

      enddo
      enddo
      enddo
      enddo

c      print*, 'LOGINIT Nk', (T0M(1,1,1,IDTNUMD+n-1),n=1,nbins)
      RETURN
      END
