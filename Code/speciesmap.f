
C     **************************************************
C     *  speciesmap                                    *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     This routine takes data from the tracer header file (BT263xxxS.COM)
C     that is entered once per species and maps it into arrays that
C     correspond to tracer numbers (i.e. one entry for bulk species and
C     several for size-resolved species).

      SUBROUTINE speciesmap()

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer n,k
      character*8 rootname    !size-resolved species name (bin # gets added)

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

C Map species names into tracer names

      !Bulk species
      do n=1,NBS
         TNAME(n)=SNAME(n)
      enddo

      !Aerosol number distribution
      rootname='ANUMD_00'
      do k=1,NBINS
         write(rootname(7:8),'(I2)') k
         TNAME(NBS+k)=rootname
      enddo

      !Aerosol mass distributions
      do n=1,(NAP+NAD)
         rootname=SNAME(NBS+n)
         do k=1,NBINS
            write(rootname(7:8),'(I2)') k
            TNAME(NBS+k)=rootname
         enddo
      enddo

C Map species properties used by wet deposition code

      do n=1,(NBS+NAP)
         if (n.le.NBS) then
            !bulk species
            NTRTY(n)=NSPTY(n)
            RKD(n)=SPRKD(n)
            DHD(n)=SPDHD(n)
         else
            !size-resolved species
            do k=1,NBINS
               NTRTY(NBS+(n-NBS)*NBINS+k)=NSPTY(n)
               RKD(NBS+(n-NBS)*NBINS+k)=0.0
               DHD(NBS+(n-NBS)*NBINS+k)=0.0
            enddo
         endif
      enddo

      !Aerosol number distribution
      do k=1,NBINS
         NTRTY(NBS+k)=TTAERO
         RKD(NBS+k)=0.0
         DHD(NBS+k)=0.0
      enddo

      RETURN
      END
