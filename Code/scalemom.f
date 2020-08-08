
C     **************************************************
C     *  scalemom                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, September 2000

C     When a tracer concentration decreases, call this routine to
C     decrease T0M and all higher order moments in proportion.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE scalemom(i,j,l,tn,f)

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM' 

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer i,j,l       !grid box
      integer tn          !tracer id number
      double precision f  !factor (0-1) by which to decrease all moments

C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      T0M(i,j,l,tn) = T0M(i,j,l,tn)*f
      TXM(i,j,l,tn) = TXM(i,j,l,tn)*f
      TYM(i,j,l,tn) = TYM(i,j,l,tn)*f
      TZM(i,j,l,tn) = TZM(i,j,l,tn)*f
      TXXM(i,j,l,tn) = TXXM(i,j,l,tn)*f
      TYYM(i,j,l,tn) = TYYM(i,j,l,tn)*f
      TZZM(i,j,l,tn) = TZZM(i,j,l,tn)*f
      TXYM(i,j,l,tn) = TXYM(i,j,l,tn)*f
      TZXM(i,j,l,tn) = TZXM(i,j,l,tn)*f
      TYZM(i,j,l,tn) = TYZM(i,j,l,tn)*f

      RETURN
      END
