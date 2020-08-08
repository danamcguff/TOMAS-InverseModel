
C     **************************************************
C     *  SOAcond                                        *
C     **************************************************

C     Written by Dana McGuffin, June 2017 
C     Based on EZCOND BY Jeff Pierce, May 2007

C     This subroutine takes a given amount of mass and condenses it
C     across the bins accordingly.  

C     ADD MORE HERE!

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
C     mcond - mass of species to condense [kg/grid cell]
C     spec - the number of the species to condense

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf - same as above, but final values

      SUBROUTINE soacond(Nki,Mki,mcond,spec,time,Nkf,Mkf)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nki(ibins), Mki(ibins, icomp)
      double precision Nkf(ibins), Mkf(ibins, icomp)
      double precision mcond, time
      integer spec

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision CS       ! condensation sink [s^-1]
      double precision sinkfrac(ibins+1) ! fraction of CS in size bin
      double precision Nk1(ibins), Mk1(ibins, icomp)
      double precision Nk2(ibins), Mk2(ibins, icomp)
      double precision madd     ! mass to add to each bin [kg]
      double precision maddp(ibins)    ! mass to add per particle [kg]
      double precision mconds ! mass to add per step [kg]
      integer nsteps            ! number of condensation steps necessary
      integer floor, ceil       ! the floor and ceiling (temporary)
      double precision eps     ! small number
      double precision tdt      !the value 2/3
      double precision mpo,mpw  !dry and "wet" mass of particle
      double precision WR       !wet ratio
      double precision tau(ibins) !driving force for condensation
      double precision totsinkfrac ! total sink fraction not including nuc bin
      double precision CSeps ! lower limit for condensation sink
      double precision fracch(ibins,icomp)
      double precision nfracch(ibins)
      double precision totch
      double precision zeros(icomp)
      double precision thresh
      double precision Mktot(ibins), med(ibins), mtot, medtot
      double precision soaareafrac, partfrac(ibins),avgfrac(ibins)
      integer absall
      double precision num_1, num_2, num_3, num_4, num_5
      double precision so4_1, so4_2, so4_3, so4_4, so4_5

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(eps=1.d-40)
      parameter(CSeps=1.d-20)
      parameter(soaareafrac=1.0d0) !fraction of SOA that goes to SA
c                                  ! instead of mass
      parameter(absall=0)          !1 for SOA to absorb to all species

C-----CODE--------------------------------------------------------------

      tdt=2.d0/3.d0

      ! initialize variables
      num_1=0d0
      so4_1=0d0
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
         num_1=num_1+Nk1(k)
         so4_1=so4_1+Mk1(k,srtso4)
      enddo

      thresh = 1.0

      call mnfix(Nk1,Mk1)
      num_2=0d0
      so4_2=0d0
      do k=1,ibins
         num_2=num_2+Nk1(k)
         so4_2=so4_2+Mk1(k,srtso4)
      enddo

      medtot = 0d0
      med(:) = 0d0
      mtot = 0d0
      Mktot(:) = 0d0
      do k = 1,ibins
        do j = 1,icomp
           Mktot(k) = Mktot(k) + Mk1(k,j)
         enddo
         mtot = mtot + Mktot(k)
         if (absall.eq.1) then !partition to all mass
                 med(k) = Mktot(k)
                 medtot = medtot + Mktot(k)
         else
                 med(k) = Mk1(k,srtc) !partition to just organic Carbon
                 medtot = medtot + Mk1(k,srtc)
         endif
      enddo

      ! Fraction to each bin for mass partitioning
      do k=1,ibins
         partfrac(k) = med(k)/medtot ! MSOA [kg] become [kg SOA per
c                                    !  total absorbing media]
      enddo

      ! get the sink fractions
      call getCondSink(Nk1,Mk1,spec,CS,sinkfrac) ! set Nnuc to zero for this calc
Cdbg      print*,'SOAcond CS=',CS
      num_3=0d0
      so4_3=0d0
      do k=1,ibins
         num_3=num_3+Nk1(k)
         so4_3=so4_3+Mk1(k,srtso4)
      enddo
cdbg      print*,'SOAcond diff Nk & SF after GETCONDSINK', num_3-num_2,
cdbg     &       so4_3-so4_2

      do k = 1,ibins
         avgfrac(k)=soaareafrac*sinkfrac(k)+
     &               (1.d0-soaareafrac)*partfrac(k)
      enddo

cdbg      print*,'CS',CS
      ! make sure that condensation sink isn't too small
       if (CS.lt.CSeps) then ! just make particles in first bin
          print*,'SOAcond CondSink too low'
          Mkf(1,spec) = Mk1(1,spec) + mcond
          Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
          do j=1,icomp
            if (icomp.ne.spec) then
              Mkf(1,j) = Mk1(1,j)
            endif
          enddo
          do k=2,ibins
            Nkf(k) = Nk1(k)
            do j=1,icomp
              Mkf(k,j) = Mk1(k,j)
            enddo
          enddo
          return
       endif
c      print*,'sinkfrac',sinkfrac
c      print*,'mcond',mcond

      ! determine how much mass to add to each size bin
      ! also determine how many condensation steps we need
      totsinkfrac = 0.d0
      do k=1,ibins
        totsinkfrac = totsinkfrac + avgfrac(k) ! get sink frac total not including nuc bin
      enddo
cdbg      print*, 'SOAcond totAVGfrac', totsinkfrac
      nsteps = 1
      do k=1,ibins
         if (avgfrac(k).lt.1.0D-20)then
            madd = 0.d0
         else
            madd = mcond*avgfrac(k)
         endif
         mpo=0.0
         do j=1,icomp-idiag
            mpo=mpo + Mk1(k,j) ! Dry mass
         enddo
         floor = int(madd*2.0/mpo)
         ceil = floor + 1
         nsteps = max(nsteps,ceil) ! don't let the mass increase by more than 10%
      enddo
c      print*,'nsteps',nsteps

      ! mass to condense each step
      mconds = mcond/nsteps

      ! do steps of condensation
      do i=1,nsteps
         if (i.ne.1) then
            call getCondSink(Nk1,Mk1,spec,
     &        CS,sinkfrac)      ! set Nnuc to zero for this calculation
            totsinkfrac = 0.d0
            do k=1,ibins
              totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
            enddo
         endif      

         do k=1,ibins
            mpo=0.0
            mpw=0.0
            !WIN'S CODE MODIFICATION 6/19/06
            !THIS MUST CHANGED WITH THE NEW dmdt_int.f
            do j=1,icomp-idiag
               mpo = mpo+Mk1(k,j)  !accumulate dry mass
            enddo
            do j=1,icomp
               mpw = mpw+Mk1(k,j)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
            if (Nk(k) .gt. 0.d0) then
               maddp(k) = mconds*avgfrac(k)/totsinkfrac/Nk1(k)
               mpw=mpw/Nk1(k)
c               print*,'mpw',mpw,'maddp',maddp(k),'WR',WR
               tau(k)=1.5d0*( ( (mpw+maddp(k) )**tdt )
     &                       -(    mpw**tdt) )  !added WR to moxid term (win, 5/15/06)
c               tau(k)=0.d0
c               maddp(k)=0.d0
               if ( tau(k) .lt. 0d0 ) then
                if (abs(tau(k)) .lt. 1d0 ) then
                    tau(k)=1.d-50
                endif
               endif
            else
               tau(k) = 0.d0
               maddp(k) = 0.d0
            endif
         enddo
c         print*,'tau',tau
c      print*, 'srtc bTmcond in ezcond', (Mk1(k,srtc),k=1,ibins )
         ! do condensation
         call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,maddp)
cdbg      print*, 'SOACOND -- aTMCOND OC', Mk2(41,srtc)
c      print*, 'srtc bMnfix in ezcond', (Mk1(k,srtc),k=1,ibins )
         call mnfix(Nk2,Mk2)
cdbg      print*, 'SOACOND: nk bin 6, 7, 8 aMNFIX',Nk2(6),Nk2(7),Nk2(8)
      num_4=0d0
      so4_4=0d0
      do k=1,ibins
         num_4=num_4+Nk2(k)
         so4_4=so4_4+Mk2(k,srtso4)
      enddo
cdbg      print*,'SOAcond diff Nk & SF after TMCOND+MNFIX2', num_4-num_3,
cdbg     &       so4_4-so4_3
cdbg      print*, 'SOACOND -- aMNFIX OC', Mk1(41,srtc)
c         call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec)
         totch=0.0
         do k=1,ibins
            nfracch(k)=(Nk2(k)-Nk1(k))
            do j=1,icomp
               if (j .ne. spec) then
               fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
               totch = totch + (Mk2(k,j)-Mk1(k,j))
               endif
            enddo
         enddo
cdbg         print*,'num fracch', nfracch
cdbg         print*,'SOACOND totch',totch
cdbg        print*,'so4 mass fracch', (fracch(k,srtso4),k=1,ibins)
cdbg        print*,'Na mass fracch', (fracch(k,srtna),k=1,ibins)
cdbg        print*,'OC mass fracch', (fracch(k,srtc),k=1,ibins)
cdbg        print*,'NH4 mass fracch', (fracch(k,srtnh4),k=1,ibins)
cdbg        print*,'H2O mass fracch', (fracch(k,srth2o),k=1,ibins)

         if (i.ne.nsteps)then
            do k=1,ibins
               Nk1(k)=Nk2(k)
               do j=1,icomp
                  Mk1(k,j)=Mk2(k,j)
               enddo
            enddo            
         endif

      enddo

      do k=1,ibins
         Nkf(k)=Nk2(k)
         do j=1,icomp
            Mkf(k,j)=Mk2(k,j)
         enddo
      enddo
      num_5=0d0
      so4_5=0d0
      do k=1,ibins
         num_5=num_5+Nkf(k)
         so4_5=so4_5+Mkf(k,srtso4)
      enddo
      print*,'SOAcond overall delatN %', 
     &         (num_5-num_1)/num_1*100.
Cdbg     &         (so4_5-so4_1)/so4_1*100.

      return
      end
