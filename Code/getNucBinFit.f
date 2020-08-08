
C     **************************************************
C     *  getNucBinFit                                  *
C     **************************************************

C     WRITTEN BY Jeff Pierce, June 2007

C     This subroutine gets the fit parameters A and B for the sub-bin number distribution
C     of the nucleation pseudo bin given by the following equation.

C     ln(dN/dlnX) = A*lnX + B

C     Where N is the number of particles, X is the mass per particle and A and B
C     are parameters fit using the total mass and total number in the bin.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nnuco - number of nucleation size particles per size bin in grid cell
C     Mtot - mass of nucleation size particle [kg/grid cell]
C     M1 - lower mass limit of bin [kg]
C     M2 - upper mass limit of bin [kg]

C-----OUTPUTS-----------------------------------------------------------

C     Ap, Bp - Fit parameters

      SUBROUTINE getNucBinFit(Nnuco,Mtot,M1,M2,Ap,Bp)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nnuco, Mtot, M1, M2, Ap, Bp

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k
      double precision guess, new !guesses for A
      double precision high, low   !bound for A
      double precision res      !residual

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

C     first, solve for A
      guess = -1.5d0
      call getResidual(guess,Nnuco,Mtot,M1,M2,res)
      if (res .gt. 0.d0) then ! we are below the solution
         low = guess
         do while (res .gt. 0.d0)
            high = low + 1.d0
            call getResidual(high,Nnuco,Mtot,M1,M2,res)
            if (res .gt. 0.d0) then ! still below solution
               low = high
            endif
         enddo
      else ! we are above the solution
         high = guess
         do while (res .lt. 0.d0)
            low = high - 1.d0
            call getResidual(low,Nnuco,Mtot,M1,M2,res)
            if (res .lt. 0.d0) then ! still below solution
               high = low
            endif
         enddo
      endif
c we now have solution bounded within 1
c now bisect 20 times to have solution within 1E-6
      do k=1,40
         new = (high+low)/2.d0
         if (new .eq. -1.d0) then
            new = -1.0000000001d0
         endif
         if (new .eq. 0.d0) then
            new = 0.0000000001d0
         endif
         call getResidual(new,Nnuco,Mtot,M1,M2,res)
         if (res .gt. 0.d0) then ! new is the new low bound
            low = new
         else                   ! new is the new high bound
            high = new
         endif
      enddo

      Ap = (high+low)/2.d0
      if (Ap.gt.130.d0) Ap = 130.d0
      if (Ap.lt.-130.d0) Ap = -130.d0

      Bp = log(Ap*Nnuco/(M2**Ap-M1**Ap))
      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE getResidual(a,Ntot,Mtot,M1,M2,res)

      implicit none

      double precision a,Ntot,Mtot,M1,M2 ! input
      double precision res      ! output

      res = log(Mtot*(a+1)*(M2**(a)-M1**(a))/
     &     (Ntot*(a)*(M2**(a+1)-M1**(a+1))))

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
