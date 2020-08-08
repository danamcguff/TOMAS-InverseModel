
C     **************************************************
C     *  cond_nuc                                      *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     This subroutine calculates the change in the aerosol size distribution
C     due to so4 condensation and binary/ternary nucleation during the
C     overal microphysics timestep.

C     ADD MORE HERE!

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Nnuci - number of nucleation size particles per size bin in grid cell
C     Mnuci - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     H2SO4rate - rate of H2SO4 chemical production [kg s^-1]
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Gcf - same as above, but final values

      SUBROUTINE cond_nuc(Nki,Mki,Gci,Nkf,Mkf,Gcf,H2SO4rate,dt,x,time,
     &                       nuc_bin,nucrate,xmeas,nucInFlg)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      double precision H2SO4rate
      double precision dt, time
      double precision x ! scaling factor
      integer nucInFlg ! 1 = Input nucrate, 0 = Output nucrate

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision CSi,CSa   ! intial and average condensation sinks
      double precision CS1,CS2       ! guesses for condensation sink [s^-1]
      double precision CStest   !guess for condensation sink
      double precision Nk1(ibins), Mk1(ibins, icomp), Gc1(icomp-1)
      double precision Nk2(ibins), Mk2(ibins, icomp), Gc2(icomp-1)
      double precision Nk3(ibins), Mk3(ibins, icomp), Gc3(icomp-1)
c      double precision Nk4(ibins), Mk4(ibins, icomp), Gc4(icomp-1)
      logical nflg ! returned from nucleation, says whether nucleation occurred or not
      double precision mcond,mcond1    !mass to condense [kg]
      double precision tol      !tolerance
      double precision eps      !small number
      double precision sinkfrac(ibins) !fraction of condensation sink coming from bin k
      double precision totmass  !the total mass of H2SO4 generated during the timestep
      double precision tmass
      double precision CSch     !fractional change in condensation sink
      double precision CSch_tol !tolerance in change in condensation sink
      double precision addt     !adaptive timestep time
      double precision time_rem !time remaining
      integer num_iter !number of iteration
      double precision sumH2SO4 !used for finding average H2SO4 conc over timestep
      integer iter ! number of iteration
      double precision fn, rnuc ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      double precision gasConc  !gas concentration [kg]
      double precision mass_change !change in mass during nucleation.f
      double precision total_nh4_1,total_nh4_2
      double precision min_tstep !minimum timestep [s]
      integer nuc_bin ! the nucleation bin
      double precision nucrate, nucrate_add !Overall Nucleation Rate during timestep
      integer xmeas

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(eps=1E-40)
      parameter(CSch_tol=0.01)
      parameter(min_tstep=1.0d0)

C-----CODE--------------------------------------------------------------

C Initialize values of Nkf, Mkf, Gcf, and time
      do j=1,icomp-1
         Gc1(j)=Gci(j)
      enddo
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
      enddo

C     Get initial condensation sink
      CS1 = 0.d0
      call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac)
c      print*,'CONDNUC -- aGetCS OC', Mk1(ibins,srtc)
c      print*,'CS1', CS1
c      CS1 = max(CS1,eps)

C     Get initial H2SO4 concentration guess (assuming no nucleation)
C     Make sure that H2SO4 concentration doesn't exceed the amount generated
C     during that timestep (this will happen when the condensation sink is very low)

C     get the steady state H2SO4 concentration
      call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,x, time)
c      print*,'gasConc',gasConc
      Gc1(srtso4) = gasConc
      addt = min_tstep
c      addt = 3600.d0
      totmass = H2SO4rate*addt*96.d0/98.d0

C     Get change size distribution due to nucleation with initial guess
c      call nucleation(Nk1,Mk1,Gc1,Nnuc1,Mnuc1,totmass,addt,Nk2,Mk2,Gc2,
c     &     Nnuc2,Mnuc2,nflg)    
      fn = 0.d0
      if ( nucInFlg .eq. 1) then
                    fn = nucrate !*addt/dt
      endif
      call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,nuc_bin,addt,x,time,xmeas,
     &                    fn, nucInFlg)

      mass_change = 0.d0
c      print*,'mass_change1',mass_change
      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
      enddo
c      print*,'mass_change2',mass_change
      mcond = totmass-mass_change ! mass of h2so4 to condense
c      print*,'after nucleation'
c      print*,'totmass',totmass,'mass_change',mass_change
c     &     ,'mcond',mcond

      if (mcond.lt.0.d0)then
         tmass = 0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tmass = tmass + Mk2(k,j)
            enddo
         enddo
c	   if (abs(mcond).gt.tmass*1.0D-8) then 
         if (abs(mcond).gt.totmass*1.0d-8) then
             if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
                tmass = 0.d0
                  do j=1,icomp-idiag
                    tmass = tmass + Mk2(nuc_bin,j)
                  enddo
                 Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                 Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                 mcond = 0.d0
             else
               print*,'mcond < 0 in cond_nuc', mcond, totmass
               stop
             endif
         else
           mcond = 0.d0
         endif
      endif
c      print*,'CONDNUC -- aNEGmcond OC', Mk2(ibins,srtc)

c      if (mcond.lt.0.d0)then
c         print*,'mcond < 0 in cond_nuc', mcond
c         stop
c      endif
      tmass = 0.d0
      do k=1,ibins
         do j=1,icomp-idiag
            tmass = tmass + Mk2(k,j)
         enddo
      enddo
c      print*, 'mcond',mcond,'tmass',tmass,'nuc',Nk2(1)-Nk1(1)

C     Get guess for condensation

      call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3)
Cjrp      mcond1 = 0.d0
Cjrp      do k=1,ibins
Cjrp         do j=1,icomp
Cjrp            mcond1 = mcond1 + (Mk3(k,j)-Mk2(k,j))
Cjrp         enddo
Cjrp      enddo
c      print*,'mcond',mcond,'mcond1',mcond1

      Gc3(srtnh4) = Gc1(srtnh4)   

      call eznh3eqm(Gc3,Mk3)
      call ezwatereqm(Mk3)

      ! check to see how much condensation sink changed
      call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac)
      CSch = abs(CS2 - CS1)/CS1    
       
c      if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
         ! get starting adaptive timestep to not allow condensationk sink
         ! to change that much
         addt = addt*CSch_tol/CSch/2
         addt = min(addt,dt)
         addt = max(addt,min_tstep)
         time_rem = dt ! time remaining
         
         num_iter = 0
         sumH2SO4=0.d0
         nucrate_add = 0.d0
         ! do adaptive timesteps
         do while (time_rem .gt. 0.d0)
            num_iter = num_iter + 1
            if ( nucInFlg .eq. 1) then
                    fn = nucrate !*addt/dt
            endif
c            print*, 'iter', num_iter, ' addt', addt
C     get the steady state H2SO4 concentration
            if (num_iter.gt.1)then ! no need to recalculate for first step
               call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,x,
     &                             time)
               Gc1(srtso4) = gasConc
            endif
c            print*,'gasConc',gasConc

            sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
            totmass = H2SO4rate*addt*96.d0/98.d0
c            call nucleation(Nk1,Mk1,Gc1,Nnuc1,Mnuc1,totmass,addt,Nk2,
c     &           Mk2,Gc2,Nnuc2,Mnuc2,nflg) 
            call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,nuc_bin,addt,x,time,
     &                          xmeas, fn, nucInFlg)
            nucrate_add = nucrate_add + fn*addt

            total_nh4_1 = 0. !Mnuc1(srtnh4)
            total_nh4_2 = 0. !Mnuc2(srtnh4)
            do i=1,ibins
               total_nh4_1 = total_nh4_1 + Mk1(i,srtnh4)
               total_nh4_2 = total_nh4_2 + Mk2(i,srtnh4)
            enddo
c            print*,'total_nh4',total_nh4_1,total_nh4_2
c            print*,'nh4 aNucleation',(Mk2(k,srtnh4),k=1,ibins)

            mass_change = 0.d0
c            print*,'mass_change1',mass_change
            do k=1,ibins
               mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
            enddo
c            print*,'mass_change',mass_change
            mcond = totmass-mass_change ! mass of h2so4 to condense
c            print*,'after nucleation'
c            print*,'totmass',totmass,'mass_change',mass_change
c     &           ,'mcond',mcond

            if (mcond.lt.0.d0)then
               tmass = 0.d0
               do k=1,ibins
                  do j=1,icomp-idiag
                     tmass = tmass + Mk2(k,j)
                  enddo
               enddo
c	         if (abs(mcond).gt.tmass*1.0D-8) then
               if (abs(mcond).gt.totmass*1.0D-8) then
                 if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
                    tmass = 0.d0
                    do j=1,icomp-idiag
                       tmass = tmass + Mk2(nuc_bin,j)
                    enddo
                    Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                    Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                    mcond = 0.d0
                 else
                   print*,'mcond < 0 in cond_nuc', mcond, totmass
                   stop
                 endif
                 else
                   mcond = 0.d0
                 endif
               endif
            
c            Gc2(srtnh4) = Gc1(srtnh4)
c            call eznh3eqm(Gc2,Mk2,Mnuc2)
c            call ezwatereqm(Mk2,Mnuc2)

c            call getCondSink(Nk2,Mk2,Nnuc2,Mnuc2,srtso4,CStest,sinkfrac)  

            Gc3(srtnh4) = Gc1(srtnh4)
            call eznh3eqm(Gc3,Mk3)
            call ezwatereqm(Mk3)
            
                                ! check to see how much condensation sink changed
            call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac)  

            time_rem = time_rem - addt
            if (time_rem .gt. 0.d0) then
               CSch = abs(CS2 - CS1)/CS1
Cjrp               if (CSch.lt.0.d0) then
Cjrp                  print*,''
Cjrp                  print*,'CSch LESS THAN ZERO!!!!!', CS1,CStest,CS2
Cjrp                  print*,'Nnuc',Nnuc1,Nnuc2
Cjrp                  print*,''
Cjrp
Cjrp                  addt = min(addt,time_rem)
Cjrp               else
               addt = min(addt*CSch_tol/CSch,addt*1.5d0) ! allow adaptive timestep to change
               addt = min(addt,time_rem) ! allow adaptive timestep to change
               addt = max(addt,min_tstep)
Cjrp               endif
c               print*,'CS1',CS1,'CS2',CS2
               CS1 = CS2
               Gc1(srtnh4)=Gc3(srtnh4)
               do k=1,ibins
                  Nk1(k)=Nk3(k)
                  do j=1,icomp
                     Mk1(k,j)=Mk3(k,j)
                  enddo
               enddo         
            endif
         enddo
         nucrate = nucrate_add/dt
cdbg         print*, '#COND_NUC: J3 = ', nucrate
         Gcf(srtso4)=sumH2SO4/dt
c         print*,'AVERAGE GAS CONC',Gcf(srtso4)
Cjrp      else
Cjrp         num_iter = 1
Cjrp         Gcf(srtso4)=Gc1(srtso4)
Cjrp      endif
      
c      print*, 'cond_nuc num_iter =', num_iter

      do k=1,ibins
         Nkf(k)=Nk3(k)
         do j=1,icomp
            Mkf(k,j)=Mk3(k,j)
         enddo
      enddo      
      Gcf(srtnh4)=Gc3(srtnh4)
c      print*,'CONDNUC -- end OC', Mkf(ibins,srtc)

      return
      end

            
      
