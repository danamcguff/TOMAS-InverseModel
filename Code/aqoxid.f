
C     **************************************************
C     *  aqoxid                                        *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000

C     This routine takes an amount of SO4 produced via in-cloud
C     oxidation and condenses it onto an existing aerosol size
C     distribution.  It assumes that only particles larger than the
C     critical activation diameter activate and that all of these have
C     grown to roughly the same size.  Therefore, the mass of SO4 
C     produced by oxidation is partitioned to the various size bins
C     according to the number of particles in that size bin.  
C     Values of tau are calculated for each size bin accordingly and
C     the cond subroutine is called to update Nk and Mk.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE aqoxid(moxid,kmin,i,j,l)

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision moxid !mass of new sulfate from in-cloud oxid.
      integer kmin           !lowest bin number activated
      integer i,j,l          !GCM grid cell

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision Nact, Mact  !#/mass of activated particles
      double precision mpo   !initial particle mass (kg)
      double precision aqtau(ibins)
      integer k,mpnum
      double precision Nko(ibins), Mko(ibins, icomp) !input to cond routine
      double precision Nkf(ibins), Mkf(ibins, icomp) !output from cond routine
      double precision tdt      !the value 2/3
      double precision eps
      integer jc
      integer tracnum
      double precision frac
      double precision moxid_o   !mass of sulfate per particle
      double precision mpw  !wet mass of particle
      double precision WR !wet ratio
      double precision mox(ibins) ! mass added per particle

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (eps=1.d-40)

C-----CODE--------------------------------------------------------------

      tdt=2.d0/3.d0

      !Fix any inconsistencies in M/N distribution (because of advection)
      call mnfix(Nk,Mk)

C Calculate which particles activate

 10   continue  !continue here is kmin has to be lowered

      Nact=0.d0
      Mact=0.d0
      do k=kmin,ibins
         Nact=Nact+Nk(k)
         do jc=1,icomp-idiag
            Mact=Mact+Mk(k,jc)
         enddo
      enddo

      if ((Mact+moxid)/(Nact+eps) .gt. xk(ibins-1)) then
         if (kmin .gt. 8) then
            kmin=kmin-1
            goto 10
         else
c            if (TAU .gt. 8350.) then
            write(*,*) 'ERROR in aqoxidcc: Ave size out of bounds'
            write(*,*) 'Box: ',i,j,l
            write(*,*) 'kmin/Nact: ',kmin,Nact
            write(*,*) 'moxid/Mact: ',moxid,Mact
            do k=1,ibins
               write(*,*) 'k, N, MSO4, MH2O: ',k,Nk(k),
     &              Mk(k,srtso4),Mk(k,srth2o)
            enddo
            Mk(kmin,srtso4)=Mk(kmin,srtso4)+moxid
            Nk(kmin)=Nk(kmin)+moxid/sqrt(xk(kmin)*xk(kmin+1))
            goto 20
c            else
c               !don't worry about the first two weeks
c               goto 20
c            endif
         endif
      endif

C Calculate tau for each size bin
      moxid_o=moxid/(Nact+eps) !now kg H2SO4 per activated particle
      do k=1,ibins
         if (k .lt. kmin) then
            !too small to activate - no sulfate for this particle
            aqtau(k)=0.d0
            mox(k)=0.d0
         else
            mpo=0.0
            mpw=0.0
            !WIN'S CODE MODIFICATION 6/19/06
            !THIS MUST CHANGED WITH THE NEW dmdt_int.f
            do j=1,icomp-idiag
               mpo = mpo+Mk(k,j)  !accumulate dry mass
            enddo
            do j=1,icomp
               mpw = mpw+Mk(k,j)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
            if (Nk(k) .gt. 0.d0) then
               mpw=mpw/Nk(k)
               aqtau(k)=1.5d0*((mpw+moxid_o*WR)**tdt-mpw**tdt)  !added WR to moxid term (win, 5/15/06)
               mox(k)=moxid_o
            else
               !nothing in this bin - set tau to zero
               aqtau(k)=0.0
               mox(k)=0.d0
            endif
         endif
      enddo

C Call condensation algorithm

         !Swap into Nko, Mko
         do k=1,ibins
            Nko(k)=Nk(k)
            do jc=1,icomp
               Mko(k,jc)=Mk(k,jc)
            enddo
         enddo

         call tmcond(aqtau,xk,Mko,Nko,Mkf,Nkf,srtso4,mox)

         !Swap out of Nkf, Mkf
         do k=1,ibins
            Nk(k)=Nkf(k)
            do jc=1,icomp
               Mk(k,jc)=Mkf(k,jc)
            enddo
         enddo

 20      continue   !go here if process is skipped

      RETURN
      END
