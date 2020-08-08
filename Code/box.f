
C     **************************************************
C     *  box                                           *
C     **************************************************

C     WRITTEN BY Peter Adams, October 1999
C         ADJUSTED by Dana McGuffin for this paper, 2019

C     This program is a box model that is used to test the size
C     resolved aerosol code for the GISS GCM.

      subroutine box(SO2ppt,No,Dp,sigma,Kp,MeasVals,PBO,NzS,NzM,Nz,Np,
     &               x1,x2,x3, xmeas, emis_indx, soa_indx, nuc_indx,
     &               meas_site, NumDays) 

C-----INCLUDE FILES--------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----VARIABLE DECLARATIONS------------------------------------------

      integer i,j,l,n,m        !counters
      integer xmeas !Scenario number 
      double precision No, Dp, sigma     !lognormal parameters (No in #/cm3,Dp in um)
      integer Np  ! Np=number of unknown parameters
      integer Nz  ! Nz=number of inventory variables
      integer NzS, NzM ! number of size cuts and moments to define inventory variables
      integer OBSVR, PBO(Np), FDBKp(Np) !PBO=1 for closed loop and 0 for open loop
      integer MeasVals(Nz), FDBKz(Nz) ! =1 to use measured value for feedback
      double precision Kp(Nz) ! observer gain [hours] for each inventory
      double precision diag(Nz), spinup, NumDays
      integer N3_6indx,N10indx,Vdryindx,Ntotindx
      double precision x1,x2,x3 !scaling factors for three parameters:
c              ! {particle emissions, soacond, nucleation rate}
      double precision x_nuc ! Nucrate scaling factor -- either input
c              ! value (sensitivity runs) or calculated parameter
      double precision time     !elapsed time (seconds)
      double precision hour  !elapsed time (hours)
      double precision adt, adt_pbio   !aerosol microphysics time step (seconds)
      double precision Dpk, vd(ibins)  !particle diameter (um) and drydep velocity (cm/s)
      double precision mnucl    !mass nucleating in a time step (kg)
      integer iseed, getlen
      double precision Ntot,Vtot     !used to check number conservation by condensation
c      double precision area     !box area (m2)
      double precision ohc   !OH conc (molec cm-3)
      logical dayflag, clearflag  !true/false for day/night or clear/cloudy
      double precision cyclet   !how far into 12-hour diurnal cycle we are (s)
      double precision setoh
      double precision u10   !wind speed 10m above ocean [m/s]
      double precision ustar  !friction wind speed for neutral ocean [m/s]
      double precision focean
      double precision SSA(ibins) !sea salt aerosol emissions [#] 
      integer nsteps            !number of timesteps per hour
      double precision endtime
      double precision h2so4diff
      double precision tot_n_1, tot_n_2
      double precision tot_s_1, tot_s_2
      double precision tot_c_1, tot_c_2, total_c, Gntot, Gmtot
      double precision tot_OC11, tot_OC12, tot_num11, tot_num12
      double precision tot_num_1, tot_so4_1, total_num, total_so4
      double precision H2SO4rate
      double precision fn, fnp, rnuc
      logical nflg
      double precision SO2ppt
      double precision CS, sinkfrac(ibins)
      external getlen, setoh
      parameter (u10 = 5.d0, focean = 0.0)
      integer length
      double precision time_prev, soa_inv_rate(Nz,1)

      integer k10, error
      integer inv_num, param_num ! Number of observers turned ON
      integer soa_indx,emis_indx,nuc_indx !Order of PBOs
      double precision tau(NBINS)
      double precision Nnucout, Mnucout(icomp) 
      double precision Nkout(NBINS), Mkout(NBINS,icomp)
      double precision Gcout(icomp-1),Gco(icomp-1)
      double precision Nko(NBINS), Mko(NBINS,icomp)
      double precision coagsink(ibins-1) ! coagulation sink [1/s]
      double precision density(ibins)    ! density [kg/m3]
      double precision total_POA   ! total emissions [Gg/yr] of primary particles
      double precision Nkout1(ibins),Mkout1(ibins)
      double precision Nkout2(ibins),Mkout2(ibins)
      double precision krain  ! loss rate by rainout [1/day]
      double precision mu, dp2, vs(ibins), Ra, g, pi !parameters for dry dep calc
      double precision Cc(ibins), lambda,Rb(ibins),kb, dpact
      double precision weight_emis, weight_dep
      double precision SOA_Nratep(ibins),SOA_Mratep(ibins,icomp)
      double precision SOA_Nratep_tot,SOA_Mratep_tot, total_n
      double precision SOA_Nratef(ibins),SOA_Mratef(ibins,icomp)
      double precision SOA_Nratef_tot,SOA_Mratef_tot
      double precision Epoa1, Epoa2, Epoaf, Evoc1, Evocf,Mpoa1,Mpoa2
      double precision Gmtot_emis, Gmtot_cond, Gmtot_nuc
      integer emis_type
      double precision  f1, dpg1, sigmag1, dpg2, sigmag2
      integer Ni1, Nf1, Ni2, Nf2
      parameter (emis_type=0) ! 1 for lognormal distribution from
c                             ! parameters above, 0 for GEOS-Chem FF,BF,BB input
      parameter (f1=4.07d-1,dpg1=16.2d-3,sigmag1=1.28d0,Ni1=1,Nf1=ibins)
      parameter (dpg2=42.4d-3,sigmag2=1.77d0,Ni2=1,Nf2=ibins)
      parameter (dpact=80.) !activation diamter for particle -> CCN [nm]
      parameter (g=9.8d0, pi=3.14, R=8.314, kb=1.38d-23) !gravitational acceleration [m/s2]
              !pi constant !universal gas constant [Pa-m3/mol/K]
              !!Boltzmann constant [J/K]

      double precision Nkf(ibins),Mkf(ibins,icomp)
      double precision Gcf(icomp-1)
      ! Observer declarations
      double precision Mk_1(ibins,icomp),Nk_1(ibins)
      double precision Mk_2(ibins,icomp),Nk_2(ibins)
      double precision Mk_3(ibins,icomp),Nk_3(ibins)
      double precision Gc_1(icomp-1), Gc_2(icomp-1)
      double precision Mk_prev(ibins,icomp),Nk_prev(ibins)
      double precision Gc_prev(icomp-1), dt, dpnm(ibins)
      double precision dpnm_prev(ibins), dpnm_edge(ibins+1)
      double precision Nnuc, Ntotal, Mdry
c      double precision Nkest(NBINS),Mkest(NBINS,icomp)
c      double precision Gcest(icomp-1)
      double precision F_prev(Nz,1), F(Nz,1)
      double precision con_affine_prev(Nz,Np),con_affine(Nz,Np)
      double precision con_affine_dummy(Nz,Np),con_affine_square(Np,Np)
      double precision theta(Np,1), theta_prev(Np,1) 
      double precision Gmk(ibins,icomp,Np), Gnk(ibins,Np)
      double precision Ggc(icomp-1,Np), Gtotalgc(icomp-1)
      double precision Gtotalmk(ibins,icomp), Gtotalnk(ibins) 
      double precision Fmk(ibins,icomp),Fnk(ibins),Fgc(icomp-1)
      double precision inventory_nominal(Nz,1), Meas(Nz,1)
      double precision con_affineE(Nz,Np)
      double precision inventoryE(Nz,1)
      double precision FE(Nz,1)
      double precision GmkE(ibins,icomp,Np), GnkE(ibins,Np)
      double precision GgcE(icomp-1,Np), GtotalgcE(icomp-1)
      double precision GtotalmkE(ibins,icomp), GtotalnkE(ibins) 
      double precision FmkE(ibins,icomp),FnkE(ibins),FgcE(icomp-1)
      double precision Ztot(NzM,1) 
      integer MeasUsed(Np), chosenMeas(Np)
      double precision, allocatable:: Measmt(:,:),dMeasmtdt(:,:)
      double precision, allocatable:: MeasDist(:,:)
      integer ROW, COL
      character*90 filename_msmt, filename_msmts, print_num

      double precision N2, N3, C2, C3, inventory(Nz,1)
      double  precision Gtot(Nz,1), theta1(Np,1), condnumb
      integer  nuc_bin, reference_date(6)
      double precision state2inventory, inventory2meas, meas2inventory
      external state2inventory, inventory2meas, meas2inventory
      integer TimeLength, NumMeas, flg
c      integer, allocatable :: timestarts(:)
      double precision, allocatable::  meas_diam(:)
      integer count_pbio
      double precision Fadd(Nz,1), Gadd(Nz,Np), Zadd(Nz,1)

      double precision Nkp(ibins),Mkp(ibins,icomp),Gcp(icomp-1)
      double precision NH3Gin,SO2in,H2SO4in,NH4Ain
      double precision NH3Gout,SO2out,H2SO4out,NH4Aout
      integer Bins_meas, Nsteps_inv, window
      parameter (Nsteps_inv = 24)
      double precision avg_rsa(Nz,Np,Nsteps_inv),meas_value(Np)
      double precision avg_rsaFrac(Nz,Np,Nsteps_inv)
      double precision avg_rsaOther(Nz,Np,Nsteps_inv)
      double precision dN_coag(ibins)

      character*50 tau_print, print_scale
      character*90 string, meas_site, fullname_data, fullname_sizeDist
      character*90 runname, datname, diagname, fname, fullname_coag
      character*90 fullname_coagdiag
      character*90 fullname_cond, fullname_param, fullname_NucRate
      character*90 fullname_measured, measurement_site, fullname_soacond
      character*90 fullname_pbio, fname1, fullname_debug_pbio
      common /FLNAMES/ datname, diagname, fullname_coag, fullname_cond, 
     &                  fullname_param, fullname_NucRate,fullname_pbio,
     &                  fullname_measured, fullname_soacond, fname1

C-----CODE-----------------------------------------------------------

      condnumb=0.d0
      iseed=101
      j=0
      MeasUsed(1)=1
      MeasUsed(2)=5
      MeasUsed(3)=16
      meas_value(:)=0.d0

      nuc_bin = 4
      spinup = 3600.*24. !1 day spin-up [seconds]
      time=0.0
      time_prev=0.0
      endtime = NumDays*24. ! hours of run time
      adt_pbio=5.d0*60.!30.d0*60. ! seconds
      window=0

      adt = 5.d0*60. ! seconds
      call boxinit(time,No,Dp,sigma)
c      print*, 'srtc:', 'T0M--',(T0M(1,1,1,idtc+k-1),k=1,nbins),
c     &         'Mk--', (Mk(j,srtc),j=1,ibins)
      hour=time/3600

c      write(*,*) 'Enter name of run'
c      read(*,'(A80)') runname
       fname1='so2_xxxxxx_nh3_xxxxxx_num_xxxxxx'
        write(fname1(5:10),'(E6.1)') SO2_print
        write(fname1(16:21),'(E6.1)') NH3ppt_o
        write(fname1(27:32),'(E6.1)') No
        if (xmeas .ne. 0) then
           fname1=trim(fname1)//trim('_Measxxx')
           write(fname1(38:40),'(I3.3)') xmeas
        endif
        print*,'filenames= ', fname1

        ! Only need 'tau' in Name if at least one estimator (OBSVR) is ON.
        OBSVR=0
        do n=1,Np
          OBSVR=OBSVR+PBO(n)
        enddo
        inv_num=0 ! add to this for the dimension of matrices input to
                ! parameter_estimate()
        do m=1,Nz
           if (MeasVals(m) .eq. 1) then
                   inv_num=inv_num+1
           endif
        enddo
        if (OBSVR .ge. 1) then 
                fname=trim(fname1)//trim('_tau')
              write( tau_print, '(F0.3)' ) Kp(1)
              fname=trim(fname)//trim('_')//trim( tau_print )
        else
                fname = fname1
        endif
        
c        write(fname(47:52),'(E6.1)') x1
c        write(fname(57:62),'(E6.1)') x2
c        write(fname(67:72),'(E6.1)') x3
        datname=trim(fname) // '.dat'
      write(*,*) datname	
	  
	  ! Create files for diagnostics

      open(unit=30,file=datname,status='new')
      close(30)
c      open(unit=31,file=diagname,status='new')
c      close(31)

      fullname_coagdiag=trim('coagDiag_')//trim(fname)//trim('.txt')
      open(unit=36,file=fullname_coagdiag,status='new')
      close(36)

      fullname_coag=trim('coagsink_')//trim(fname)//trim('.txt')
      open(unit=32,file=fullname_coag,status='new')
      close(32)
      fullname_cond=trim('condsink_')//trim(fname)//trim('.txt')
      open(unit=33,file=fullname_cond,status='new')
      close(33)
      fullname_param=trim('parameter_')//trim(fname)//trim('.txt')
      open(unit=34,file=fullname_param,status='new')
      close(34)
      fullname_NucRate=trim('NucRate_')//trim(fname)//trim('.txt')
      open(unit=35,file=fullname_NucRate,status='new')
      close(35)
c      open(unit=36,file=trim("anthro_emis_")//trim(fname)//trim('.txt')
c     &       ,status='new')
c      close(36)
      fullname_measured=trim('measurements_')//trim(fname)//trim('.txt')
      fullname_data=trim('/Users2/danamcg/EUSAAR/')//trim(meas_site)//
     &                         trim('/MeasInfo.txt')
      open(unit=40,file=fullname_measured,status='new')
      open(71, file=fullname_data)
      read(71, *) (reference_date(j),j=1,6)
      read(71, *) Bins_meas
      ALLOCATE( meas_diam(Bins_meas) )
      read(71, *) (meas_diam(j), j = 1,Bins_meas)
      read(71, *) measurement_site
      close(71)
      write(40, *) reference_date
      write(40, *) measurement_site
      write(40, *) Bins_meas
      write(40, *) meas_diam
      close(40)
      fullname_soacond=trim('soacond_')//trim(fname)//trim('.txt')
      open(unit=41,file=fullname_soacond,status='new')
      close(41)
      fullname_pbio=trim('pbio_data')//trim(fname)//trim('.m')
      open(unit=51,file=fullname_pbio,status='new')
      write(51,*) 'i=1;', 'Kc=diag([',(Kp(m),m=1,Nz), ']);'
      close(51)
      fullname_debug_pbio=trim('pbio_debug')//trim(fname)//trim('.m')
      open(unit=91,file=fullname_debug_pbio,status='new')
      write(91,*) 'i=0;'
      close(91)
      fullname_sizeDist=trim('SizeDist_')//trim(fname)//trim('.txt')
      open(unit=92,file=fullname_sizeDist,status='new')
      close(92)

      ! Read & Extract measurements to a large array

      filename_msmt=trim('/Users2/danamcg/EUSAAR/')//
     &                  trim(measurement_site)//trim('/Meas_inventory')
      open(70,file=trim(filename_msmt)//trim('01.txt') )
      read(70,*,IOSTAT=error) ROW
      read(70,*,IOSTAT=error) COL
      allocate( Measmt(ROW,Nz+1) )
      allocate( MeasDist(ROW,COL) )
      allocate( dMeasmtdt(ROW,Nz) )
      close(70)
      DO k=1,Nz     
        write( print_num, '(I2.2)') k
       filename_msmts=trim(filename_msmt)//trim(print_num)//trim('.txt')
        open(70,file=filename_msmts)
        read(70,*,IOSTAT=error) ROW
        read(70,*,IOSTAT=error) COL
        DO I=1,ROW
          read(70,*,IOSTAT=error) Measmt(I,1),(MeasDist(I,J),J=1,COL),
     &                           Measmt(I,k+1), dMeasmtdt(I,k)
        ENDDO
        close(70)
      ENDDO
      T0M(1,1,1,IDTSO2)=SO2ppt/1d12*(boxmass/29.d0)*64.d0

      print*, '#box BOXMASS (kg)', boxmass
      do n=1,NBINS
         Nk(n)=T0M(1,1,1,IDTNUMD-1+n)
         Mk(n,srtso4)=T0M(1,1,1,IDTSO4-1+n)
         Mk(n,srtna)=T0M(1,1,1,IDTNA-1+n) 
         Mk(n,srtc)=T0M(1,1,1,IDTC-1+n) 
         Mk(n,srth2o)=T0M(1,1,1,IDTH2O-1+n)
      enddo

      ! Initialize mass for Ammonia & H2SO4 gas (dlm, 03/2017)
      T0M(1,1,1,IDTNH3G)=NH3ppt_o/1d12*(boxmass/29.d0)*17.d0
      Gc(srtso4)=T0M(1,1,1,IDTH2SO4)
       call NH3_GISStoTOMAS(T0M(1,1,1,IDTNH3G),
     &                    T0M(1,1,1,IDTNH4A),Gc,Mk)
      call eznh3eqm(Gc,Mk)
      call ezwatereqm(Mk)
      call getCondSink(Nk,Mk,srtso4,CS,sinkfrac)

      !Gas-phase SO2 to SO4 during day hours
      if (mod(hour,24.d0) .le. 12.d0) then
         T0M(1,1,1,IDTH2SO4)=T0M(1,1,1,IDTH2SO4)+300.*adt/3600.
      endif
      if (mod(hour,24.d0) .le. 12.d0) then
         dayflag=.true.
      else
         dayflag=.false.
      endif
      if (hour .ge. 192.d0) then
         clearflag=.false.
      else
         clearflag=.true.
      endif
      dayflag=.true.
      clearflag=.true.
      cyclet=dble(mod(hour,24.d0))*3600.d0 !sec past sunup
      ohc=setoh(dayflag,clearflag,cyclet)
      call boxchem(T0M(1,1,1,IDTSO2),T0M(1,1,1,IDTH2SO4),ohc,adt,
     &             temp,pres)
      call getH2SO4prod(T0M(1,1,1,IDTSO2),ohc,adt,
     &             temp,pres,H2SO4rate)
        call getH2SO4conc(H2SO4rate,CS,Gc(srtnh4),gasConc,x2,time)
        T0M(1,1,1,IDTH2SO4)=gasConc
      Gc(srtso4)=T0M(1,1,1,IDTH2SO4)
c      print*, 'aINITIALchem srtc:', (Mk(k,srtc),k=1,nbins)

      do k=1,ibins-1
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &   0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &   Mk(k,srth2o),Mk(k,srtc) )
         coagsink(k)=0.d0
      enddo
         density(ibins)=aerodens(Mk(ibins,srtso4),0.d0,
     &   0.1875d0*Mk(ibins,srtso4),Mk(ibins,srtna), !assume bisulfate
     &   Mk(ibins,srth2o),Mk(ibins,srtc) )
      call getNucRate(Gc, fn, rnuc, nflg,x2,time,xmeas)
      call report(time,CS,coagsink,density,fn)

       do n=1,Np
         if(n.eq.emis_indx) then
             theta1(n,1)=x1
         elseif(n.eq.nuc_indx) then 
             theta1(n,1)=x2
         elseif(n.eq.soa_indx) then
             theta1(n,1)=x3
         endif
       enddo

      Fadd = 0.d0
      Gadd = 0.d0
      Zadd = 0.d0
      count_pbio = 0
 10   continue  !beginning of new time step

      ! Save state-space before TOMAS
      Nko=Nk
      Mko=Mk
      Gco=Gc
      SO2in=T0M(1,1,1,IDTSO2)
      NH3Gin=T0M(1,1,1,IDTNH3G)
      H2SO4in=T0M(1,1,1,IDTH2SO4)
      NH4Ain=T0M(1,1,1,IDTNH4A)

      tot_num11 = 0.d0
      tot_OC11 = 0.d0
      do k = 1,ibins
        tot_num11 = tot_num11 + Nk(k)
        do j = 1,icomp
          tot_OC11 = tot_OC11 + Mk(k,j)
        enddo
      enddo
      soa_nratep(:)=0.d0
      soa_mratep(:,:)=0.d0
      soa_nratep_tot=0.d0
      soa_mratep_tot=0.d0
      fnp = 0.d0
      Epoa1 = 0.d0
      Mpoa1 = 0.d0
      Mpoa2 = 0.d0
	  
	  ! Calculate diameter and density of each size bin
      do k=1,ibins
         density(k)=aerodens(Mko(k,srtso4),0.d0,
     &     0.1875d0*Mko(k,srtso4),Mko(k,srtna), !assume bisulfate
     &     Mko(k,srth2o),Mko(k,srtc)  )
c         print*,'density', density(k)
         dpnm(k) = (  6./pi*sqrt( xk(k+1)*xk(k) )
     &                 /density(k)  )**(1./3.)*1.d9 !nm
      enddo
      dpnm_edge(1) = ( 6./pi*xk(1)/density(1) )**(1./3.)*1.d9
      do k=2,ibins+1
         dpnm_edge(k) = ( 6./pi*xk(k)/density(k-1) )**(1./3.)*1.d9
      enddo
	  
	  ! Get bin number of nucleated ( 3 nm ) particles
      nuc_bin = 1
      do while ( 3.d0 .gt. dpnm_edge(nuc_bin+1))
            nuc_bin = nuc_bin+1
      enddo
      print*, '#box nuc_bin from initial state:', nuc_bin

Cdbg      print*, 'BOX b4 original TOMAS, NH4A', NH4Ain
C Do TOMAS with nominal scaling factors (theta1)
c========================================================
      call tomas(theta1,Np,emis_indx,soa_indx,nuc_indx,xmeas,
     &           Nko,Mko,Gco,SO2in,NH3Gin,H2SO4in,NH4Ain,
     &           Nkf,Mkf,Gcf,SO2out,NH3Gout,H2SO4out,NH4Aout,
     &           soa_nratep,soa_mratep,soa_nratep_tot,
     &           soa_mratep_tot,adt,time,coagsink,fnp,
     &           nuc_bin, Epoa1,Epoa2,emis_type,Mpoa1,Mpoa2,dN_coag,
     &           Evoc1,Nz,soa_inv_rate,Bins_meas,NzM,meas_diam)
c========================================================

C Write new states to T0M
      print*,'#BOX nucBIN dp upper bound = ',dpnm_edge(nuc_bin+1) 
      Nkp=Nkf
      Mkp=Mkf
      Gcp=Gcf
      T0M(1,1,1,IDTSO2)=SO2out
      T0M(1,1,1,IDTNH3G)=NH3Gout
      T0M(1,1,1,IDTH2SO4)=H2SO4out
      T0M(1,1,1,IDTNH4A)=NH4Aout


C Do parameter estimation from passivity based observer
c===========================================================
Cc Multi-Input Multi-Output: Np parameters & Nz outputs
c===========================================================
c      dt=adt
      do k=1,ibins
         density(k)=aerodens(Mko(k,srtso4),0.d0,
     &     0.1875d0*Mko(k,srtso4),Mko(k,srtna), !assume bisulfate
     &     Mko(k,srth2o),Mko(k,srtc)  )

         dpnm(k) = (  6./pi*sqrt( xk(k+1)*xk(k) )
     &                 /density(k)  )**(1./3.)*1.d9 !nm
      enddo
      dpnm_edge(1) = ( 6./pi*xk(1)/density(1) )**(1./3.)*1.d9
      do k=2,ibins+1
         dpnm_edge(k) = ( 6./pi*xk(k)/density(k-1) )**(1./3.)*1.d9
      enddo
      nuc_bin = 1
      do while ( 3.d0 .gt. dpnm_edge(nuc_bin+1))
            nuc_bin = nuc_bin+1
      enddo
	  
	  ! Calculate inventory variables before and after running the nominal TOMAS code
      do n=1,Nz
        inventory(n,1)= state2inventory(Mko, Nko,Gco,dpnm,dpnm_edge,
     &                         meas_diam(Bins_meas),density, nuc_bin,
     &                         adt,n,NzM,0 )
        inventory_nominal(n,1)= state2inventory(Mkp, Nkp,Gcp,dpnm,
     &                               dpnm_edge,meas_diam(Bins_meas),
     &                               density,nuc_bin,adt,n,NzM,0)
      enddo
        if (time .lt. 360d0) dpnm_prev=dpnm
      Gtotalmk(:,:)=0.d0
      Gtotalnk(:)=0.d0
      Gtotalgc(:)=0.d0
        inv_num=0 ! add to this for the dimension of matrices input to
                ! parameter_estimate()
        do m=1,Nz
           if (MeasVals(m) .eq. 1) then
                   inv_num=inv_num+1
           endif
        enddo
	  ! Get rate of change of inventory variables with respect to uncertain processes
	  !    (con_affine)
      param_num=0
      DO n=1,Np
       call affine_dynamics(Gmk(:,:,n),Gnk(:,n),Ggc(:,n),
     &       nuc_bin,fnp,Epoa1,Epoa2,Mpoa1,Mpoa2,dpg1,dpg2,sigmag1,
     &       sigmag2,Ni1,Nf1,Ni2,Nf2,adt,emis_type,SOA_Nratep,
     &       SOA_Mratep,Evoc1,n,emis_indx, nuc_indx, soa_indx, 
     &       theta1,Np, time)
         do m=1,Nz
           con_affine(m,n) =state2inventory(Gmk(:,:,n),Gnk(:,n),Ggc(:,n)
     &                            ,dpnm,dpnm_edge,meas_diam(Bins_meas)
     &                            ,density,nuc_bin,adt,m,NzM,1)
         enddo !loop over inventories

        if (PBO(n) .eq. 1 ) then
          param_num=param_num+1
          Gtotalmk=Gtotalmk+Gmk(:,:,n)
          Gtotalnk=Gtotalnk+Gnk(:,n)
          Gtotalgc=Gtotalgc+Ggc(:,n)
        endif
      ENDDO !loop over parameters

      con_affine_dummy = con_affine
           
          tot_num12 = 0.d0
          tot_OC12 = 0.d0
          Gntot = 0.d0
          Gmtot = 0.d0
          do k = 1,ibins
              Gntot = Gntot + Gtotalnk(k)
              tot_num12 = tot_num12 + Nk(k)

              do j = 1,icomp
                 Gmtot = Gmtot + Gtotalmk(k,j)
                 tot_OC12 = tot_OC12 + Mk(k,j)
              enddo
          enddo
cdbg          print*, '#box all other processes dM= ', tot_OC12 - tot_OC11
cdbg     &              - Gmtot*adt,
cdbg     &          ' dN = ', tot_num12-tot_num11-Gntot*adt
cdbg          print*, '#box uncertain processes dM= ', Gmtot*adt, ' dN = ', 
cdbg     &               Gntot*adt
cdbg      print*, 'control affine matrix', con_affine

      IF( time .ge. spinup ) then 
      ! Calculate full f(x) vector
        call remaining_dynamics(Gtotalmk, Gtotalnk, Gtotalgc,
     &              Mko, Mkp, Nko, Nkp, Gco,
     &              Gcp, adt, Fmk, Fnk, Fgc )
      ! Calculate f(x) for each inventory!
        do n=1,Nz
          F(n,1)=state2inventory(Fmk,Fnk,Fgc,dpnm,dpnm_edge,
     &                           meas_diam(Bins_meas),density,
     &                           nuc_bin,adt,n,NzM,1)
        enddo
        F = (inventory_nominal-inventory)/adt
        do m=1,Np
          F(:,1)=F(:,1)-con_affine(:,m)
        enddo

         ! Ztot : 0-3rd moments integrated over size range measured
		 !  to normalize equations
         Ztot(:,:)=0.d0
         do n=1,NzM
          Ztot(n,1)=state2inventory(Mko,Nko,Gco,dpnm,dpnm_edge,
     &                           meas_diam(Bins_meas),density,
     &                           nuc_bin,adt,Nz+n,NzM,0)
         enddo
         FDBKp=PBO
         FDBKz=MeasVals

        if ( OBSVR .lt. 1) then
           theta = theta1 
        else

                chosenMeas=MeasUsed
            call parameter_estimate(time,FDBKp,FDBKz,Np,Nz,
     &               inventory,inventory_nominal,F_prev,F,
     &               con_affine_prev,con_affine_dummy,Kp,adt_pbio,1,
     &               theta_prev,xmeas,NzM,Ztot,emis_indx,nuc_indx,
     &               soa_indx,spinup,inv_num,param_num,chosenMeas, 
     &               Measmt,dMeasmtdt,MeasDist,COL,ROW,  theta,condnumb)
        endif

         
      ! Calculate new states with parameter estimates
      DO n=1,Np
        theta(n,1)= theta(n,1)*theta1(n,1)
       ! Turn ON/OFF feedback...
       if ( (PBO(n) .eq. 0) .or. (time .lt. 6d0*3600d0) ) then 
        theta(n,1) = theta1(n,1)
       elseif ( mod(time, adt_pbio) .ne. 0 ) then
         theta(n,1)=theta_prev(n,1)
       endif
!      ! Set bounds on parameters
        if ( theta(n,1) .lt. 1d-3 ) theta(n,1)=1d-3
        if ( theta(n,1) .gt. 1d3 ) theta(n,1)=1d3

      ENDDO
      ELSE ! initialize previous values and set previous parameters to 1
               F(:,:)=0d0
               theta = theta1
               condnum=0.d0
       ENDIF

Cdbg      print*, '#box dMtot before Est SOA = ', soa_mratep_tot*adt
Cdbg      print*, '#box dN before Est POA = ', Epoa1
C Do TOMAS with estimated scaling factors (theta)
c========================================================
          call tomas(theta,Np,emis_indx,soa_indx,nuc_indx,xmeas,
     &           Nko,Mko,Gco,SO2in,NH3Gin,H2SO4in,NH4Ain,
     &           Nkf,Mkf,Gcf,SO2out,NH3Gout,H2SO4out,NH4Aout,
     &           soa_nratep,soa_mratep,soa_nratep_tot,
     &           soa_mratep_tot,adt,time,coagsink,fnp,
     &           nuc_bin,Epoa1,Epoa2,emis_type,Mpoa1,Mpoa2,dN_coag,
     &          Evoc1,Nz,soa_inv_rate,Bins_meas,NzM,meas_diam)
c========================================================

C Write new states to T0M
          Nk=Nkf
          Mk=Mkf
          Gc=Gcf
          T0M(1,1,1,IDTSO2)=SO2out
          T0M(1,1,1,IDTNH3G)=NH3Gout
          T0M(1,1,1,IDTH2SO4)=H2SO4out
          T0M(1,1,1,IDTNH4A)=NH4Aout

     ! Save inventory for next timestep
      do k=1,ibins
Cdbg         print*, '#BOX aerodens bin #',k
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
         dpnm(k) = (  6./pi*sqrt( xk(k+1)*xk(k) )
     &                 /density(k)  )**(1./3.)*1.d9 !nm
      enddo
      dpnm_edge(1) = ( 6./pi*xk(1)/density(1) )**(1./3.)*1.d9
      do k=2,ibins+1
         dpnm_edge(k) = ( 6./pi*xk(k)/density(k-1) )**(1./3.)*1.d9
      enddo
      nuc_bin = 1
      do while ( 3.d0 .gt. dpnm_edge(nuc_bin+1))
            nuc_bin = nuc_bin+1
      enddo
      print*,'#estBOX nucBIN dp upper bound = ',dpnm_edge(nuc_bin+1) 
      do n=1,Nz
       inventoryE(n,1)=state2inventory(Mk,Nk,Gc,dpnm,dpnm_edge,
     &                      meas_diam(Bins_meas),density,nuc_bin,dt,
     &                      n,NzM,0)
      enddo
      do n=1,Nz
        Meas(n,1)=inventory2meas( inventoryE(n,1) )
      enddo
      do n=1,Np
        if (PBO(n) .eq. 0)then
           theta(n,1)=theta1(n,1)
        endif
      enddo
      Nk_prev=Nk
      Mk_prev=Mk
      Gc_prev=Gc
      con_affine_prev=con_affine
      F_prev=F
      theta_prev=theta
      dpnm_prev=dpnm
c***********************

      GtotalmkE(:,:)=0.d0
      GtotalnkE(:)=0.d0
      GtotalgcE(:)=0.d0
c      print*, '#box EPOA1: ', Epoa1, 'Nucrate: ', fnp, 'SOA_Mratep: ',
c     &      SOA_Mratep
c      print*, 'box: nucBin =', nuc_bin
      DO n=1,Np
       call affine_dynamics(GmkE(:,:,n),GnkE(:,n),GgcE(:,n),
     &       nuc_bin,fnp,Epoa1,Epoa2,Mpoa1,Mpoa2,dpg1,dpg2,sigmag1,
     &       sigmag2,Ni1,Nf1,Ni2,Nf2,adt,emis_type,SOA_Nratep,
     &      SOA_Mratep,Evoc1,n,emis_indx,nuc_indx,soa_indx,theta,
     &      Np, time)
         do m=1,Nz
           con_affineE(m,n) = state2inventory(GmkE(:,:,n),GnkE(:,n),
     &                          GgcE(:,n),dpnm,dpnm_edge,
     &                          meas_diam(Bins_meas),density,
     &                           nuc_bin,dt,m,NzM,1)
         enddo !loop over inventories

        if (PBO(n) .eq. 1 ) then
          GtotalmkE=GtotalmkE+GmkE(:,:,n)
          GtotalnkE=GtotalnkE+GnkE(:,n)
          GtotalgcE=GtotalgcE+GgcE(:,n)
        endif
      ENDDO !loop over parameters

          tot_num12 = 0.d0
          tot_OC12 = 0.d0
          Gntot = 0.d0
          Gmtot = 0.d0
          Gmtot_emis = 0.d0
          Gmtot_cond = 0.d0
          Gmtot_nuc = 0.d0
          do k = 1,ibins
              Gntot = Gntot + GtotalnkE(k)
              tot_num12 = tot_num12 + Nk(k)
              do j = 1,icomp
                 Gmtot_emis = Gmtot_emis + GmkE(k,j,emis_indx)
                 Gmtot_cond = Gmtot_cond + GmkE(k,j,soa_indx)
                 Gmtot_nuc = Gmtot_nuc + GmkE(k,j,nuc_indx)
                 Gmtot = Gmtot + GtotalmkE(k,j)
                 tot_OC12 = tot_OC12 + Mk(k,j)
              enddo
          enddo
cdbg          print*, '#box all other processes dM= ', tot_OC12 - tot_OC11
cdbg     &              - Gmtot*adt,
cdbg     &          ' dN = ', tot_num12-tot_num11-Gntot*adt
c          print*, '#box dM from SOA, Nuc, POA =',Gmtot_cond*adt,','
c     &             ,Gmtot_nuc*adt,',',Gmtot_emis*adt
cdbg          print*, '#box uncertain processes dM= ', Gmtot*adt, ' dN = ', 
cdbg     &               Gntot*adt

cdbg      print*, 'control affine matrix', con_affine
      IF( time .gt. 360d0) then
      ! Calculate full f(x) vector
        call remaining_dynamics(GtotalmkE, GtotalnkE, GtotalgcE,
     &              Mko, Mkf, Nko, Nkf, Gco,
     &              Gcf, adt, FmkE, FnkE, FgcE )
      ! Calculate f(x) for each inventory!
cdbg       print*, '# of moments=',NzM
       do n=1,Nz
          FE(n,1)=state2inventory(FmkE,FnkE,FgcE,dpnm,dpnm_edge,
     &                           meas_diam(Bins_meas), density, 
     &                           nuc_bin,dt,n,NzM,1)
        enddo
       ELSE
               FE(:,1)=0.d0
       ENDIF


      !Write estiamted inventories and parameters to txt file!
11    format(E15.6, A2,E15.6,A2,*(G0.4,:,","),*(G0.4,:,",") )
      open(unit=34, file=fullname_param, status='old', access='append')
      if ( OBSVR .eq. 0 ) then
      write(34,11) time, ',',condnumb, ',', (Meas(n,1), n=1,Nz),
     &                         (theta1(n,1),n=1,Np) !
      else
      write(34,11) time, ',',condnumb, ',',  (Meas(n,1), n=1,Nz),
     &               (theta(n,1),n=1,Np)!
      endif
      close(34)
12    format(E15.6,A2, G0.4,A2, G0.4)
c12    format(E15.6,A2, *(G0.4,:,","),A2, *(G0.4,:,",") )
      open(unit=41,file=fullname_soacond, status='old', access='append')
c      write(41,12) time,',', (SOA_Nratep(n)*adt,n=1,ibins),',',
c     &                 ((SOA_Mratep(n,j)*adt,n=1,ibins),j=1,icomp-1)
      write(41,12) time, ',', SOA_Mratep_tot, ',',
     &               SOA_Nratep_tot
      close(41)

15    format(E15.6,A2, *(G0.4,:,",") )
      open(unit=41,file=fullname_coagdiag,status='old', access='append')
      write(41,15) time, ',', (dN_coag(i),i=1,ibins)
      close(41)

13    format(E15.6,A2, G0.4)
      open(unit=35,file=fullname_Nucrate, status='old', access='append')
      write(35,13) time, ',', fnp!/x_nuc
      close(35)
14    format(E15.6, A2,*(G0.4,:,",") )
      open(unit=92,file=fullname_sizeDist,status='old', access='append')
      write(92,14) time, ',', (Nk(n),n=1,ibins)
      close(92)
c      Endifi
      ! Debugging file -- is there a difference between the two TOMAS
      ! calls?

20    format(A9)
21    format(A13,E15.6,A3,A13,99(E15.6,","),A3,A13,99(E15.6,","),A3)
22    format(A13,3(E15.6,","),A3,A13,3(E15.6,","),A3,A13,3(E15.6,",")
     &           ,A3)
23    format(A13,3(E15.6,","),A3,3(E15.6,","),A3,3(E15.6,",")
     &           ,A3)

      
      open(unit=51,file=fullname_pbio,status='old',access='append')
      write(51,*) 'i=i+1;'
      close(51)

      !Swap Nk and Mk arrays back to T0M
      T0M(1,1,1,IDTH2SO4)=Gc(srtso4)
      T0M(1,1,1,IDTNH3G)=Gc(srtnh4)
      T0M(1,1,1,IDTNH4A)=0.d0
      do n=1,nbins
        T0M(1,1,1,IDTNH4A)=T0M(1,1,1,IDTNH4A)+Mk(n,srtnh4)
      enddo

      do n=1,NBINS
         T0M(1,1,1,IDTNUMD-1+n)=Nk(n)
         T0M(1,1,1,IDTSO4-1+n)=Mk(n,srtso4)
         T0M(1,1,1,IDTNA-1+n)=Mk(n,srtna)
         T0M(1,1,1,IDTC-1+n)=Mk(n,srtc)
         T0M(1,1,1,IDTH2O-1+n)=Mk(n,srth2o)
      enddo

      call getCondSink(Nk,Mk,srtso4,CS,sinkfrac)
   
C Check for negative tracer problems
      do n=1,ntm
         if (T0M(1,1,1,n) .lt. 0.0) then
            if (abs( T0M(1,1,1,n)/boxvol*1.d18 ) .lt. 1.d-4 ) then
                    T0M(1,1,1,n)=1.d-4*boxvol*1.d-18
            else
            write(*,*) '#boxERROR: Tracer ',n,' < 0',' value:',
     &                 T0M(1,1,1,n)
            STOP
            endif
         endif
      enddo
C End of time step
      time=time+adt
      hour=time/3600
      if (mod(hour,1.d0) .eq. 0.d0) then
         write(*,*) 'Hour ', hour, ' completed'
         do k=1,ibins
           density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
         enddo
         
         call report(time,CS,coagsink,density,fn)
c      print*, 'finalSO2conc: ppt', T0M(1,1,1,IDTSO2)/64.d0/
c     &                          (boxmass/29.d0)*1d12
      print*, '#box finalH2SO4gas', T0M(1,1,1,IDTH2SO4), Gc(srtso4)
      endif


      if (hour .lt. endtime) goto 10

      return
      END   !of main


C     **************************************************
C     *  TOMAS                                         *
C     **************************************************

      subroutine tomas(param, Np, emis_indx, soa_indx, nuc_indx,xmeas,
     &                   Nk1,Mk1,Gc1,SO2,NH3G,H2SO4,NH4A,
     &                   Nkf,Mkf,Gcf,SO2f,NH3Gf,H2SO4f,NH4Af,
     &                   soa_nrate,soa_mrate,soa_nrate_tot,
     &                  soa_mrate_tot,adt,time,coagsink,nucrate,
     &                  nuc_bin, POA1,POA2,emisFLG,POAff,POAbb,dN,
     &                   SOA, Nz,soa_inventory,Bins_meas,NzM,meas_diam)

C-----INCLUDE FILES--------------------------------------------------
      IMPLICIT NONE

      include 'sizecode.COM'

C-----VARIABLE DECLARATIONS------------------------------------------

      integer i,j,l,n,m, Np,k,Nz        !counters
      integer flgtype ! 1 = input process rates, 0= input parameters
      integer emisFLG ! flag for POA emissions scheme
      double precision time     !elapsed time (seconds)
      double precision hour  !elapsed time (hours)
      double precision adt      !aerosol microphysics time step (seconds)
      double precision Dpk, vd(ibins)  !particle diameter (um) and drydep velocity (cm/s)
      double precision mnucl    !mass nucleating in a time step (kg)
      integer iseed, getlen
      double precision Ntot,Vtot     !used to check number conservation by condensation
      double precision ohc   !OH conc (molec cm-3)
      logical dayflag, clearflag  !true/false for day/night or clear/cloudy
      double precision cyclet   !how far into 12-hour diurnal cycle we are (s)
      double precision setoh
      double precision u10   !wind speed 10m above ocean [m/s]
      double precision ustar  !friction wind speed for neutral ocean [m/s]
      double precision focean
      double precision SSA(ibins) !sea salt aerosol emissions [#] 
      integer nsteps            !number of timesteps per hour
      double precision endtime
      double precision h2so4diff
      double precision tot_n_1, tot_n_2, tot_num3
      double precision tot_s_1, tot_s_2, tot_OC3
      double precision tot_c_1, tot_c_2, total_cond
      double precision tot_OC1, tot_OC2, tot_num1, tot_num2
      double precision tot_num_1, tot_so4_1, total_num, total_so4
      double precision H2SO4rate
      double precision fn, rnuc
      logical nflg
      double precision CS, sinkfrac(ibins)
      external getlen, setoh
      parameter (u10 = 5.d0, focean = 0.d0)
      integer length, nuc_bin, xmeas
      double precision Param(Np,1), param_apply
      double precision inventory(Nz,1), inventory_update(Nz,1)
      integer Bins_meas, NzM
      double precision soa_inventory(Nz,1), meas_diam(Bins_meas)

      double precision nucrate, mnuc ! cm-3s-1
      integer soa_indx,emis_indx,nuc_indx !Order of PBOs
      double precision tau(ibins)
      double precision Nnucout, Mnucout(icomp) 
      double precision Nkout(ibins), Mkout(ibins,icomp)
      double precision Gcout(icomp-1)
      double precision coagsink(ibins-1) ! coagulation sink [1/s]
      double precision density(ibins)    ! density [kg/m3]
      double precision total_POA   ! total emissions [Gg/yr] of primary particles
      double precision Nkout1(ibins),Mkout1(ibins)
      double precision Nkout2(ibins),Mkout2(ibins)
      double precision krain  ! loss rate by rainout [1/day]
      double precision mu, dp2, vs(ibins), Ra, g, pi, R !parameters for dry dep calc
      double precision Cc(ibins), lambda,Rb(ibins),kb, dpact
      double precision  OC1, f1, dpg1, sigmag1, total_OC
      double precision OC2, dpg2, sigmag2
      integer Ni1, Nf1, Ni2, Nf2
      double precision POA1, POA2, POAmass, POAnum, POAff,POAbb
      double precision POA, SOA, SOAfrac ![kg/s] & [-]
      double precision SOA_Nrate(ibins),SOA_Mrate(ibins,icomp)
      double precision SOA_Nrate_tot,SOA_Mrate_tot, total_n
      double precision SOA_Nrate_tot1,SOA_Mrate_tot1,condnuc_mass
      double precision SOA_Nrate_tot2,SOA_Mrate_tot2, condnuc_num
      double precision soa_nrate2(ibins), soa_mrate2(ibins,icomp)
      double precision dnumdt, dmassdt, dpnm(ibins), after_cond
      double precision Ntemp(ibins), dN(ibins), dN_frac(ibins)

c      parameter (f1=5.84d-1,dpg1=18.d-3,sigmag1=1.46d0,Ni1=1,Nf1=18)
      parameter (f1=3.3745d-2,dpg1=20.8d-3,sigmag1=1.293d0,Ni1=1,Nf1=18)
      parameter (dpg2=30.d-3,sigmag2=2.2d0,Ni2=19,Nf2=ibins)
c      parameter (dpg1=16.2d-3,sigmag1=1.28d0,Ni1=1,Nf1=ibins)
c      parameter (dpg2=42.4d-3,sigmag2=1.77d0,Ni2=1,Nf2=ibins)
c      parameter (SOAfrac=1.d0, total_OC=10.**(10.5515) )
      parameter (SOAfrac=1.d0, total_OC=1.5134d18*3.d0 ) ! [-] & [1/s]
cBOS      parameter (f1=4.07d-1,SOAfrac=1.d0, total_OC=1.102d21*6.d0 ) ! [-] & [1/s]
c      parameter (f1-1.d0)
c      parameter (SOAfrac=1.d0, total_OC=1.5134d19 ) ! [-] & [1/s]
      parameter (OC1=total_OC*f1,OC2=total_OC*(1.d0-f1) )
c      parameter (OC1=total_OC,OC2=total_OC*(1.d0-f1) )
c      parameter (SOA=SOAfrac*11.4D-2*1.0D9/365.25/24. )
c      parameter (SOA=SOAfrac*3.2387d0 ) ! [kg/s] : 1786 kton/month (Fountoukis et al., 2012)
c      parameter (SOA=SOAfrac*3.2387d-4 ) ! [kg/s] : 1786 kton/month (Fountoukis et al., 2012)
      parameter (dpact=80.) !activation diamter for particle -> CCN [nm]
      parameter (g=9.8d0, pi=3.14, R=8.314, kb=1.38d-23) !gravitational acceleration [m/s2]
              !pi constant !universal gas constant [Pa-m3/mol/K]
              !!Boltzmann constant [J/K]

      double precision Mkf(ibins,icomp),Nkf(ibins)
      double precision Mk1(ibins,icomp),Nk1(ibins)
      double precision Gcf(icomp-1), Gc1(icomp-1)
      double precision SO2 ,H2SO4 ,NH3G ,NH4A
      double precision SO2f,H2SO4f,NH3Gf,NH4Af
      double precision Vol_sulf, Nsulf
      double precision dpnm_edge(ibins+1)
      character*50 tau_print
      character*90 runname, datname, diagname, fname, fullname_coag
      character*90 fullname_cond, fullname_param, fullname_NucRate
      character*90 fullname_measured, measurement_site, fullname_soacond
      character*90 fullname_pbio
      common /FLNAMES/ datname, diagname, fullname_coag, fullname_cond, 
     &                  fullname_param, fullname_NucRate,fullname_pbio,
     &                  fullname_measured, fullname_soacond, fname

      double precision coag_num_change, coag_mass_change
      double precision aerodens, state2inventory
      external aerodens, state2inventory
C-----CODE-----------------------------------------------------------

      hour=time/3600.

      Nk=Nk1
      Mk=Mk1
      Gc=Gc1
      Nkf(:)=0.d0
      Mkf(:,:)=0.d0
      Gcf(:)=0.d0
      SO2f=SO2
      NH3Gf=NH3G
      H2SO4f=H2SO4
      NH4Af=NH4A

      Ntemp = Nk

      tot_num1 = 0.d0
      tot_OC1 = 0.d0
      do k = 1,ibins
        tot_num1 = tot_num1 + Nk(k)
        do j = 1,icomp
          tot_OC1 = tot_OC1 + Mk(k,j)
        enddo
      enddo
Cdbg      print*,'TOMAS nh4a init & final',NH4A,NH4Af
C Emissions

       ! SO2
       !3.57e3 kg/s = 56.3 Tg S/yr [S & P 3rd Edition]
       !1675 kton/month (Fountoukis et al., 2012)
cSMR        SO2f=SO2+3.037d0*6.d1*adt
cPuydeDome       SO2f=SO2+3.037d0*7.d1*adt
        SO2f=SO2+3.037d-1*adt*9.d-1 
		
       do k=20,ibins
        Mk(k,srtso4)=Mk(k,srtso4)+3.037d-2*adt/64.*96./(ibins-19.)
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
         dpnm(k) = (  6./pi*sqrt( xk(k+1)*xk(k) )
     &                 /density(k)  )**(1./3.)*1.d9 !nm
         Vol_sulf=4./3.*pi*(dpnm(k)/2.*1.d-9)**3
         Nsulf=3.037d-2*adt/64.*96./(ibins-19.)/1350.d0/Vol_sulf
         Nk(k)=Nk(k)+Nsulf
c         print*, 'Bin# ',k,' Number sulfate emiss = ', Nsulf
       enddo

      !ammonia
      !501 kton/month (Fountoukis et al., 2012)
cSMR      NH3Gf=NH3G+0.9085d1*6.d1*adt
cPuydeDome      NH3Gf=NH3G+0.9085d1*7.d1*adt
      NH3Gf=NH3G+0.9085d1*6.d1*adt !cBOS

C Sea Spray Emissions
      !sea-salt
      call emisnaclarke(u10,focean,adt,SSA)
c      print*, 'BOX adt', adt,'BOX SSA', (SSA(k),k=1,ibins)
      do k=1,ibins
         Nk(k)=Nk(k)+SSA(k)
         Mk(k,srtna)=Mk(k,srtna)+SSA(k)*sqrt( xk(k)*xk(k+1) )
      enddo  


Cc Do water eqm at appropriate times
      call eznh3eqm(Gc,Mk)
      call ezwatereqm(Mk)

      call NH3_GISStoTOMAS(NH3Gf,NH4Af,Gc,Mk)
Cdbg      print*,'TOMAS nh4a aGISS init & final',NH4A,NH4Af
C Dynamics
      if ((Mk(15,srtso4).lt.0.d0).or.(Mk(16,srtso4).lt.0.d0).or.
     &    (Mk(17,srtso4).lt.0.d0).or.(Mk(18,srtso4).lt.0.d0).or.
     &    (Mk(19,srtso4).lt.0.d0) )  then
              print*, 'Negative tracer after AmSulf',
     &                        (Mk(i,srtso4),i=15,19)
      endif
    
C Chemistry
      !Gas-phase SO2 to SO4 during day hours
      if (mod(hour,24.d0) .le. 12.d0) then
         H2SO4f=H2SO4+300.*adt/3600.
      endif
      if (mod(hour,24.d0) .le. 12.d0) then
         dayflag=.true.
      else
         dayflag=.false.
      endif
      if (hour .ge. 192.d0) then
         clearflag=.false.
      else
         clearflag=.true.
      endif
      dayflag=.true.
      clearflag=.true.
      cyclet=dble(mod(hour,24.d0))*3600.d0 !sec past sunup
      ohc=setoh(dayflag,clearflag,cyclet)
c      print*, 'SO2conc before chem: ppt', T0M(1,1,1,IDTSO2)/64.d0/
c     &                          (boxmass/29.d0)*1d12
      call boxchem(SO2f,H2SO4f,ohc,adt,
     &             temp,pres)
      call getH2SO4prod(SO2f,ohc,adt,
     &             temp,pres,H2SO4rate)
c      H2SO4rate = 5.0d2 ! set this as constant for now
c      print*,'H2SO4rate(kg/s)',H2SO4rate, ', mlc/cm3/s',
c     &        H2SO4rate/boxvol*6.022D23/0.098 
c      print*, 'SO2conc: ppt', T0M(1,1,1,IDTSO2)/64.d0/
c     &                          (boxmass/29.d0)*1d12
c      print*, 'BOXMASS (kg)', boxmass
c      print*, 'aCHEM srtc:', (T0M(1,1,1,idtc+k-1),k=1,nbins)
cdbg      print*, 'aCHEM nk for bins 6,7,8', Nk(6), Nk(7), Nk(8)

C ****************
C Aerosol dynamics
C ****************

      !Swap gas phase data into Gc
      Gc(srtso4)=H2SO4f
c      Gc(srtnh4)=T0M(1,1,1,IDTNH3)
c      Gc(srtno3)=boxmass*HNO3(1,1,1,1)*63./28.9
c      print*, 'H2SO4gas', T0M(1,1,1,IDTH2SO4), Gc(srtso4)

C Rainout
      krain=1/(12.d0*3600.d0) !First order loss rate constant of 12 hours
      call storenm()
      do k=1,iBINS
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
         dp2 = (6./pi*xk(k)/density(k))**(1./3.)*1d9 !nm
c         if (dp2 .gt. dpact ) then
c                 print*, 'bin', k, 'do rainout loss'
           Nk(k)=Nk(k)*exp(-krain*adt)
           do j=1,icomp
             Mk(k,j)=Mk(k,j)*exp(-krain*adt)
           enddo
c         endif
      enddo
      do j=1,icomp-1
          Gc(j)=Gc(j)*exp(-krain*adt)
      enddo
      SO2f=SO2f*exp(-krain*adt)
      call aerodiag(5,1,1,1)

      if ((Mk(15,srtso4).lt.0.d0).or.(Mk(16,srtso4).lt.0.d0).or.
     &    (Mk(17,srtso4).lt.0.d0).or.(Mk(18,srtso4).lt.0.d0).or.
     &    (Mk(19,srtso4).lt.0.d0) )  then
              print*, 'Negative tracer after 1stOrder loss',
     &                        (Mk(i,srtso4),i=15,19)
      endif
C Moist convection

C Large-scale clouds

C Radiation

C Dry deposition
      ! No dry deposition ( assume this is taken care of with the 12-hr loss rate above)
      call storenm()
      mu=2.5277e-7*temp**0.75302 !air viscosity [kg/m s]
      lambda = 2.*mu/pres*( pi*R*temp/(28.9d0*8.)*1d3 )**0.5
      ustar = u10*0.4d0/log(10./5) !neutral central business district
Cdbg      print*,'air viscocity kg/m s', mu
      do k=1,IBINS
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
Cdbg         print*,'density', density(k)
         dp2 = (6./pi*xk(k)/density(k))**(2./3.) !m2
Cdbg         print*,'diam^2 m2', dp2
         Cc(k) = 1. + 2.*lambda/dp2**0.5*(1.257 +
     &            0.4*exp(-1.1*dp2**0.5/lambda) )
         vs(k) = dp2*density(k)*g/(18*mu)*Cc(k)
         Rb(k)=1./( 3.*ustar*(
     &      1./(mu/(pres*28.9d0/R*kb*Cc(k)/3./pi/mu/dp2**0.5))**0.54
     &       + 1./2.*dp2/area**2) ) ![s/m]
Cdbg         print*,'settling velocity m/s',vs(k)
         Ra=1.d4 !very stable conditions [s/m]
         vd(k) = ( 1/( Ra+Rb(k)+Ra*Rb(k)*vs(k) )+vs(k))*100.!*weight_dep ![cm/s]
c         Nk(k)=Nk(k)*exp(-vd(k)/boxvol*area*1d4*adt)
cc         print*, 'bin k',k,'dry dep velocity cm/s',vd(k)
c         do j=1,icomp
c            Mk(k,j)=Mk(k,j)*exp(-vd(k)/boxvol*area*1d4*adt)
c         enddo
      enddo
c      do j=1,icomp-1
c          Gc(j)=Gc(j)*exp(-1.d0/Ra/boxvol*area*1d4*adt) !*weight_dep
c      enddo

      if ((Mk(15,srtso4).lt.0.d0).or.(Mk(16,srtso4).lt.0.d0).or.
     &    (Mk(17,srtso4).lt.0.d0).or.(Mk(18,srtso4).lt.0.d0).or.
     &    (Mk(19,srtso4).lt.0.d0) )  then
              print*, 'Negative tracer after DryDepsn',
     &                        (Mk(i,srtso4),i=15,19)
      endif


C Do water eqm at appropriate times
      call eznh3eqm(Gc,Mk)
      call ezwatereqm(Mk)


      ! get the total mass of N
      tot_n_1 = Gc(srtnh4)/17.
      do k=1,ibins
         tot_n_1 = tot_n_1 + Mk(k,srtnh4)/18.
      enddo

      ! get the total mass of S
      tot_s_1 = SO2f/64.
      tot_s_1 = tot_s_1 + H2SO4rate*adt/98.
      do k=1,ibins
         tot_s_1 = tot_s_1 + Mk(k,srtso4)/96.
      enddo

      ! get the total moles of C
      tot_num_1=0d0
      tot_so4_1=0d0
      do k=1,ibins
         tot_num_1 = tot_num_1+Nk(k)
         tot_so4_1 = tot_so4_1+Mk(k,srtso4)
      enddo

      tot_OC1=0.d0
      do k=1,ibins
        do j=1,icomp
          tot_OC1=tot_OC1+Mk(k,j)
        enddo
      enddo
      !Coagulation
      do k=1,ibins-1
        coagsink(k)=0.d0
      enddo
      call storenm()
      call multicoag(adt, coagsink)

      call mnfix(Nk,Mk)
      do k = 1,ibins
        dN(k) = Nk(k)-Ntemp(k)
        dN_frac(k) = dN(k)/Ntemp(k)
      enddo
c      Ntemp = Nk
c      print*, 'fractional Delta N after coagulation, 3-10 nm', dN_frac
      ! get the total number of kmol nh3
      tot_n_2 = Gc(srtnh4)/17.
      do k=1,ibins
         tot_n_2 = tot_n_2 + Mk(k,srtnh4)/18.
      enddo

      ! get the total mass of S
      tot_s_2 = SO2f/64.
c      tot_s_2 = tot_s_2 + H2*32./98.
      do k=1,ibins
         tot_s_2 = tot_s_2 + Mk(k,srtso4)/96.
      enddo

      ! get the total number of kmol C
      tot_c_2 = Gc(srtc)/12.
      do k=1,ibins
         tot_c_2 = tot_c_2 + Mk(k,srtc)/12.
      enddo

c      print*,'*****Balance Checks [kmoles]*******'
c      print*, 'TOTAL N', tot_n_1, tot_n_2
c      print*, 'TOTAL S', tot_s_1, tot_s_2
c      print*, 'TOTAL C', tot_c_1, tot_c_2
C      print*, 'H2SO4', Gc(srtso4)
c      print*, '************************'

C Do water eqm at appropriate times
      call eznh3eqm(Gc,Mk)
      call ezwatereqm(Mk)
      if ((Mk(15,srtso4).lt.0.d0).or.(Mk(16,srtso4).lt.0.d0).or.
     &    (Mk(17,srtso4).lt.0.d0).or.(Mk(18,srtso4).lt.0.d0).or.
     &    (Mk(19,srtso4).lt.0.d0) )  then
              print*, 'Negative tracer after COAG',
     &                        (Mk(i,srtso4),i=15,19)
      endif

C ***********************
C End of aerosol dynamics
C ***********************
      
      tot_num2 = 0.d0
      tot_OC2 = 0.d0
      do k = 1,ibins
        tot_num2 = tot_num2 + Nk(k)
        do j = 1,icomp
          tot_OC2 = tot_OC2 + Mk(k,j)
        enddo
      enddo

      print*, 'deltaM after Coag %', (tot_OC2-tot_OC1)/tot_OC1*100
      do k = 1,ibins
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
         dpnm(k) = (  6./pi*sqrt( xk(k+1)*xk(k) )
     &                 /density(k)  )**(1./3.)*1.d9 !nm
      enddo

C Secondary organic aerosol production
      Ntemp=Nk
      ! Put VOC in 'gas-phase Carbon' tracer
cc    VOC emissions
c      print*, soa_indx
      SOA=SOAfrac*3.2387d-4*7.2d4*param(soa_indx,1)
c    ![kg/s] : 1786 kton/month (Fountoukis et al., 2012)

      tot_OC1 = 0.d0
      do k = 1,ibins
        do j = 1,icomp
          tot_OC1 = tot_OC1 + Mk(k,j)
        enddo
      enddo
      Gc(srtc)=Gc(srtc)+SOA*adt
      tot_c_1 = Gc(srtc)/12.
      do k=1,ibins
         tot_c_1 = tot_c_1 + Mk(k,srtc)/12.
      enddo
cc    Condensation of VOC into particle OC
Cdbg      print*, 'TOMAS b4 SOACOND Nk(1) & Mk(1)',Nk(1),(Mk(1,j),j=1,icomp)
      call soacond(Nk,Mk,SOA*adt,srtc,time, Nkout,Mkout)
c      call soacond(Nk,Mk,Gc(srtc),srtc, Nkout,Mkout)
       do n=1,Nz
          inventory(n,1)=state2inventory(Mk,Nk,Gc,dpnm,dpnm_edge,
     &                           meas_diam(Bins_meas), density, 
     &                           nuc_bin,adt,n,NzM,0)
        enddo


      total_cond=0d0
      total_num=0d0
      total_so4=0d0
      soa_Nrate_tot=0d0
      soa_Mrate_tot=0d0
      do n=1,iBINS
         soa_Nrate(n)=(Nkout(n)-Nk(n))/adt
         soa_Nrate_tot=soa_Nrate_tot + soa_Nrate(n)
         total_cond = total_cond+Mkout(n,srtc)-Mk(n,srtc)
         do j = 1,icomp
           soa_Mrate(n,j)=(Mkout(n,j)-Mk(n,j))/adt
           soa_Mrate_tot=soa_Mrate_tot + soa_Mrate(n,j)
           Mk(n,j)=Mkout(n,j)
         enddo
         Nk(n)=Nkout(n)
         density(n)=aerodens(Mk(n,srtso4),0.d0,
     &     0.1875d0*Mk(n,srtso4),Mk(n,srtna), !assume bisulfate
     &     Mk(n,srth2o),Mk(n,srtc) )
         dpnm(n) = (  6./pi*sqrt( xk(n+1)*xk(n) )
     &                 /density(n)  )**(1./3.)*1.d9 !nm
      enddo
      do k = 1,ibins
        dN(k) = Nk(k)-Ntemp(k)
        dN_frac(k) = dN(k)/Ntemp(k)
      enddo
       do n=1,Nz
          inventory_update(n,1)=state2inventory(Mk,Nk,Gc,dpnm,dpnm_edge,
     &                           meas_diam(Bins_meas), density, 
     &                           nuc_bin,adt,n,NzM,0)
        enddo
        soa_inventory = (inventory_update-inventory)/adt

      Gc(srtc) = Gc(srtc)-total_cond
      if ((Mk(15,srtso4).lt.0.d0).or.(Mk(16,srtso4).lt.0.d0).or.
     &    (Mk(17,srtso4).lt.0.d0).or.(Mk(18,srtso4).lt.0.d0).or.
     &    (Mk(19,srtso4).lt.0.d0) )  then
              print*, 'Negative tracer after SOAcond',
     &                        (Mk(i,srtso4),i=15,19)
      endif

      after_cond = 0.d0
      do k = 1,ibins
        do j = 1,icomp
          after_cond = after_cond + Mk(k,j)
        enddo
      enddo
      
c      print*, '#TOMAS dM SOA cond = ', after_cond - tot_OC2

      !Condensation/Evaporation

C Do water eqm at appropriate times
      call eznh3eqm(Gc,Mk)
Cdbg      print*, 'TOMAS aEZNH3 Nk(1) & Mk(1)',Nk(1),(Mk(1,j),j=1,icomp)
      call ezwatereqm(Mk)
Cdbg      print*, 'TOMAS aEZWATER Nk(1) & Mk(1)',Nk(1),(Mk(1,j),j=1,icomp)

      !Condensation and nucleation (coupled)
	  !  for this case, we turned off the coupled condensation and nucleation
	  !   and set the baseline nucleation rate to a constant rate of 0.2 #/cm3-s
      call storenm()
	  
	  nucrate = 2.d-1*param(nuc_indx,1) !1/cm3-s

c      tot_c_1 = Nk(nuc_bin)
Cdbg      print*, 'TOMAS b4 CONDNUC Nk(1) & Mk(1)',Nk(1),(Mk(1,j),j=1,icomp)
c      call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,H2SO4rate,adt,
c     &               param(nuc_indx,1),time,nuc_bin,nucrate,xmeas,0)
C!      call nucleation(Nk,Mk,Gc,Nkout,Mkout,Gcout,nuc_bin,adt,x2,time)
      do k = 1,ibins
         density(k)=aerodens(Mk(k,srtso4),0.d0,
     &     0.1875d0*Mk(k,srtso4),Mk(k,srtna), !assume bisulfate
     &     Mk(k,srth2o),Mk(k,srtc) )
         dpnm(k) = (  6./pi*sqrt( xk(k+1)*xk(k) )
     &                 /density(k)  )**(1./3.)*1.d9 !nm
      enddo
      dpnm_edge(1) = ( 6./pi*xk(1)/density(1) )**(1./3.)*1.d9
      do k=2,ibins+1
         dpnm_edge(k) = ( 6./pi*xk(k)/density(k-1) )**(1./3.)*1.d9
      enddo
      nuc_bin = 1
      do while ( 3.d0 .gt. dpnm_edge(nuc_bin+1))
            nuc_bin = nuc_bin+1
      enddo
      print*,'#TOMAS nucBIN dp upper bound = ',dpnm_edge(nuc_bin+1) 
      
      Mkout=Mk
      Nkout=Nk
      Gcout=Gc
c      mnuc = ( 4.d0/3.d0*pi*(dpnm(nuc_bin)/2.*1.d-9)**3 )*1350.d0
      mnuc = ( 4.d0/3.d0*pi*(3./2.*1.d-9)**3 )*1350.d0
      Nkout(nuc_bin) = Nk(nuc_bin)+nucrate*boxvol*adt
      Mkout(nuc_bin,srtso4)=Mkout(nuc_bin,srtso4)+nucrate*boxvol*
     &                                                mnuc*adt
      h2so4diff = 0.d0
      do k=1,ibins
         h2so4diff = h2so4diff + Mkout(k,srtso4) - Mk(k,srtso4)
      enddo

      condnuc_num = -nucrate*adt*boxvol
      condnuc_mass = -nucrate*adt*boxvol*mnuc
      do n=1,iBINS
         CondNuc_num = condnuc_num + Nkout(n) - Nk(n)
         do j=1,icomp
            condnuc_mass = condnuc_mass + Mkout(n,j) - Mk(n,j)
            Mk(n,j)=Mkout(n,j)
         enddo
         Nk(n)=Nkout(n)
      enddo
      Gc(srtnh4) = Gcout(srtnh4)
      Gc(srtso4) = Gcout(srtso4)
      H2SO4f=Gc(srtso4)
c      dN = Nk(nuc_bin:nuc_bin+7)-Ntemp(nuc_bin:nuc_bin+7)
c      Ntemp = Nk
c      print*, 'Delta N after Nucleation, 3-10 nm', dN
Cdbg      print*, '#TOMAS dN Nucleation = ', Nk(nuc_bin) - tot_c_1
      call aerodiag(1,1,1,1)      
c      call getNucRate(Gc,nucrate,rnuc,nflg,param(nuc_indx,1),time)

c      print*, '#TOMAS Nucleation rate = ', nucrate
c      print*, '#TOMAS dM Nucleation = ', nucrate*boxvol*adt*mnuc
C Primary Aerosol Emissions
      POA1 = param(emis_indx,1)*OC1*adt
c      POA2 = 0.d0
      POA2 = param(emis_indx,1)*OC2*adt
      ! emission flux = \mu*( OC*[OrgMatter/OrgC] + BC )
      POAff = param(emis_indx,1)*(1.413d-12*1.8d0+1.74d-12)*adt*area
      POAbb = param(emis_indx,1)*((9.03d-13+1.636d-11)*1.8d0+
     &                               1.664d-13+1.299d-12)*adt*area
      call logemiss(POA1,POAff,POAbb,dpg1,sigmag1,Ni1,Nf1,
     &                ocden,srtc,emisFLG,Nkout1,Mkout1)!,fname)
      do k=1,ibins
       Nk(k)=Nk(k)+Nkout1(k)
       Mk(k,srtc)=Mk(k,srtc)+Mkout1(k)
      enddo


c      call eznh3eqm(Gc,Mk)
c      call ezwatereqm(Mk)

      tot_num3 = 0.d0
      tot_OC3 = 0.d0
      POAnum = 0.d0
      POAmass = 0.d0
      do k = 1,ibins
        POAnum = Nkout1(k) + POAnum
        POAmass = Mkout1(k) + POAmass
        tot_num3 = tot_num3 + Nk(k)
        do j = 1,icomp
          tot_OC3 = tot_OC3 + Mk(k,j)
        enddo
      enddo


      Nkf=Nk
      Mkf=Mk
      Gcf=Gc
      dnumdt = 0.d0
      dmassdt = 0.d0
      do k = 1,ibins
         dnumdt = (Nkf(k)-Nk1(k))/adt +dnumdt
         do j = 1,icomp
            dmassdt = (Mkf(k,j)-Mk1(k,j))/adt + dmassdt
          enddo
       enddo
c       print*, '   dN/dt=', dnumdt, ' dM/dt = ', dmassdt
c      print*, '#TOMAS POA: ', POA, 'Nucrate: ', nucrate, 'SOA_Mrate: ',
c     &      SOA_Mrate


      END

