
C     **************************************************
C     *  nucleation                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and gets the binary and ternary nucleation rates.
C     The number of particles added to the first size bin is calculated
C     by comparing the growth rate of the new particles to the coagulation
C     sink.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Gcf - same as above, but final values

      SUBROUTINE nucleation(Nki,Mki,Gci,Nkf,Mkf,Gcf,nuc_bin,dt,x,
     &                        time,xmeas,fn_inout,nucInFlg)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      integer nuc_bin, xmeas, nucInFlg
      double precision dt, time
      double precision x ! scaling factor
      logical nflg

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision nh3ppt   ! gas phase ammonia in pptv
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision fn, fn_inout       ! nucleation rate cm-3 s-1
      double precision rnuc     ! critical nucleation radius [nm]
      double precision gtime    ! time to grow to first size bin [s]
      double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      double precision Mktot    ! total mass in bin
      double precision neps
      double precision meps
      double precision density(ibins)  ! density of particle [kg/m3]
      double precision pi
      double precision frac     ! fraction of particles growing into first size bin
      double precision d1,d2    ! diameters of particles [m]
      double precision mp       ! mass of particle [kg]
      double precision mold     ! saved mass in first bin
      double precision mnuc     !mass of nucleation
      double precision dpnm(ibins) ! mean diameters [nm]

C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens
      external aerodens

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159)

C-----CODE--------------------------------------------------------------

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
c      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12
      nh3ppt = NH3ppt_o

      fn = 0.d0
      rnuc = 0.d0
      gtime = 0.d0
      ! Start with zero number and mass of FINAL particles (dlm,04/2017)
      do k=1,ibins
        do j=1,icomp
           Mkf(k,j)=0.d0
        enddo
        Nkf(k)=0.d0
      enddo


C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
c            print*, 'napari'
            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
         elseif (bin_nuc.eq.1) then
c            print*, 'vehk'
            call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
         endif    
c         print*,'srtc bin5 after NUCL routine',Mki(5,srtc),
c     &           Mkf(5,srtc)
      endif
c      fn = 1.d0 ! cm-3 s-1
      If ( NucInFlg .eq. 1 ) then
              fn = x*fn_inout
      Else
              fn = x*fn
      Endif
      fn_inout = fn
C     if nucleation occured, see how many particles grow to join the first size
C     section
      if (fn.gt.0.d0) then
c         print*, 'entered fn'
c         print*, 'fn', fn, 'rnuc',rnuc
         ! get the time it takes to grow to first size bin
Cjrp         d1 = rnuc*2.d0*1E-9
Cjrp         d2 = (6.d0*sqrt(xk(1)*xk(2))/1800.d0/pi)**(1.d0/3.d0)
Cjrp         print*, 'd1', d1, 'd2', d2
Cjrp         call getGrowthTime(d1,d2,Gc(srtso4),temp,
Cjrp     &        boxvol,gtime)
Cjrp         print*, 'gtime', gtime
Cjrp         ! get the first order loss rate of these particles due to coagulation
Cjrp         ! get loss rate for cluster size
Cjrp         d1 = rnuc*2.d0*1E-9
Cjrp         call getCoagLoss(d1,ltc1)
Cjrp         print*, 'd1', d1, 'ltc1', ltc1
Cjrp         ! we want the loss rate for particles at the size where they are added
Cjrp         ! to bin 1
Cjrp         Mktot = 0.d0
Cjrp         do j=1,icomp
Cjrp            Mktot=Mktot+Mk(1,j)
Cjrp         enddo
Cjrp         Mktot=Mktot+2.d0*Mk(1,srtso4)/96.d0-Mk(1,srtnh4)/18.d0 ! account for h+
Cjrp
Cjrp         if (Mktot.gt.meps)then
Cjrp            density=aerodens(Mk(1,srtso4),0.d0,Mk(1,srtnh4),
Cjrp     &           0.d0,Mk(1,srth2o)) !assume bisulfate
Cjrp         else
Cjrp            density = 1400.
Cjrp         endif
Cjrp
Cjrp         if(Nk(1).gt.neps .and. Mktot.gt.meps)then
Cjrp            mp=Mktot/Nk(1)
Cjrp         else
Cjrp            mp=sqrt(xk(1)*xk(2))
Cjrp         endif
Cjrp         
Cjrp         d1 = (6.d0*mp/density/pi)**(1.d0/3.d0)
Cjrp         call getCoagLoss(d1,ltc2)
Cjrp         print*, 'd1', d1, 'ltc2', ltc2 
Cjrp
Cjrp         ! take average of loss rate constants
Cjrp         ltc = (ltc1 + ltc2)/2.d0
Cjrp         print*, 'ltc', ltc
Cjrp
Cjrp	 if (gtime.lt.0.d0)then
Cjrp	    gtime=0.d0
Cjrp         endif
Cjrp
Cjrp         frac = exp(-ltc*gtime)
Cjrp         print*,'frac', frac
Cjrpc         frac = 1.d0

        do k=1,ibins
         density(k)=aerodens(Mki(k,srtso4),0.d0,
     &                0.1875d0*Mki(k,srtso4),Mki(k,srtna),
     &                Mki(k,srth2o),Mki(k,srtc) )
         dpnm(k)=( 6./pi*sqrt( xk(k+1)*xk(k) )/density(k) )**(
     &              1./3.)*1d9 !nm
        enddo
        nuc_bin = 0
        do while (3. .gt. dpnm(nuc_bin+1))
            nuc_bin = nuc_bin + 1
        enddo
c        nuc_bin = 7
         mnuc = (4.d0/3.d0*pi*(dpnm(nuc_bin)/2.*1d-9)**3)*1350.d0
c         print*,'#NUCLEATION J3: nuc_bin',nuc_bin,' J3 = ', fn,
c     &          'DeltaT = ', dt

         mold = Mki(nuc_bin,srtso4)
         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+fn*mnuc*
     &        boxvol*dt
         Nkf(nuc_bin) = Nki(nuc_bin)+fn*boxvol*dt
c         print*,'Nk_NUUUC',Nki(nuc_bin),Nkf(nuc_bin)
         Gcf(srtso4) = Gci(srtso4)! - (Mkf(nuc_bin,srtso4)-mold)
c         print*,'srtc nuc bin', nuc_bin,' in NUCL',Mki(nuc_bin,srtc),
c     &             Mkf(nuc_bin,srtc)
C there is a chance that Gcf will go less than zero because we are artificially growing
C particles into the first size bin.  don't let it go less than zero.
c         if (Gcf(srtso4).lt.0.d0)then
c            Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4) + 
c     &           Gci(srtso4)*96./98.
c            Nkf(nuc_bin) = Nki(1) + Gci(srtso4)*96./98./mnuc
c            Gcf(srtso4) = 0.d0
c            print*,'used up h2so4 in nuc subroutine'
c         endif

      endif

c         print*,'srtc nuc bin', nuc_bin,' in NUCL',Mki(nuc_bin,srtc),
c     &             Mkf(nuc_bin,srtc)
      do k=1,ibins
         if (k .ne. nuc_bin)then
            Nkf(k) = Nki(k)
            do i=1,icomp
               Mkf(k,i) = Mki(k,i)
            enddo
         endif
      enddo

Cdbg      print*,'#nucleation dM',Mkf(nuc_bin,srtso4)-Mki(nuc_bin,srtso4)

      RETURN
      END

