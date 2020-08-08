
C     **************************************************
C     *  NH3_GISStoTOMAS                               *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine puts ammonia to the particle phase until 
C     there is 2 moles of ammonium per mole of sulfate and the remainder
C     of ammonia is left in the gas phase.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE NH3_GISStoTOMAS(giss_nh3g,giss_nh4a,Gce,Mke)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Gce(icomp)
      double precision Mke(ibins,icomp)
c      double precision Mnuce(icomp)
      double precision giss_nh3g, giss_nh4a

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k
      double precision tot_nh3  !total kmoles of ammonia
      double precision tot_so4  !total kmoles of so4
      double precision sfrac    !fraction of sulfate that is in that bin

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      ! get the total number of kmol nh3
      tot_nh3 = giss_nh3g/17.d0 + giss_nh4a/18.d0

      ! get the total number of kmol so4
c      tot_so4 = Mnuce(srtso4)/96.d0
       tot_so4=0.d0
      do k=1,ibins
         tot_so4 = tot_so4 + Mke(k,srtso4)/96.d0
      enddo

      ! see if there is free ammonia
      if (tot_nh3/2.d0.lt.tot_so4)then  ! no free ammonia
         Gce(srtnh4) = 0.d0 ! no gas phase ammonia
c         sfrac = Mnuce(srtso4)/96.d0/tot_so4
c         Mnuce(srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
         do k=1,ibins
            sfrac = Mke(k,srtso4)/96.d0/tot_so4
            Mke(k,srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
         enddo
      else ! free ammonia
c         Mnuce(srtnh4) = Mnuce(srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
         do k=1,ibins
            Mke(k,srtnh4) = Mke(k,srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
         enddo
         Gce(srtnh4) = (tot_nh3 - tot_so4*2.d0)*17.d0 ! put whats left over in the gas phase
      endif

c      print*, 'srtc in NH3GISS',( Mke(k,srtc),k=1,ibins)
      RETURN
      END












