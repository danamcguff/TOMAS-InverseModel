
C     **************************************************
C     *  aerodiag                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000

C     Accumulates diagnostics on aerosol microphysical processes.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE aerodiag(ptype,i,j,l)

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer ptype, i,j,l,c

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

C Bulk species
Cjrp      AEROD(j,l,IDTH2SO4,ptype)=AEROD(j,l,IDTH2SO4,ptype)
Cjrp     &     +Gc(srtso4)-Gcd(srtso4)
      AEROD(j,l,IDTH2SO4,ptype)=
     &     Gc(srtso4)-Gcd(srtso4)

C Aerosol number
      do k=1,NBINS
Cjrp         AEROD(j,l,IDTNUMD-1+k,ptype)=AEROD(j,l,IDTNUMD-1+k,ptype)
Cjrp     &        +Nk(k)-Nkd(k)
         AEROD(j,l,IDTNUMD-1+k,ptype)=
     &        Nk(k)-Nkd(k)
      enddo

C Aerosol mass
      do k=1,NBINS
         do c=1,icomp-idiag
C     jrp               AEROD(j,l,IDTSO4-1+(c-1)*NBINS+k,ptype)
C     jrp     &              =AEROD(j,l,IDTSO4-1+(c-1)*NBINS+k,ptype)
C     jrp     &              +Mk(k,c)-Mkd(k,c)
            AEROD(j,l,IDTSO4-1+(c-1)*NBINS+k,ptype)
     &           =Mk(k,c)-Mkd(k,c)
         enddo
      enddo

      RETURN
      END
