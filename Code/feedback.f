
C      ****************************************************
C      * FEEDBACK                                         *
C      ****************************************************

C      WRITTEN BY Dana McGuffin, April 2017

C      This routine calculates scaling factors for uncertain processes
C      based on the error between the model and ground-based measurements.
C      The parameters are calculated based on an asymptotically-converging
C      feedback loop.

!====================================================================
      ! MEASUREMENT_DATA !!
!====================================================================


      SUBROUTINE measurement_data(time,dt,Ndt,COL,ROW,
     &                         MEAStime,MEAS,dMEASdt,MEASdist,  Istar, 
     &                         Measured,Measured_numberDist,dMeasuredDt)

      IMPLICIT NONE

C-----INCLUDE FILES---------------------------------------------------

c      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      integer x, j !which inventory-measurement pair?
      INTEGER error, I, ROW, RUNS, Istar, xmeas, COL, Ndt
      double precision MEAStime(ROW),MEASdist(ROW,COL)
      double precision MEAS(ROW), dMEASdt(ROW)
      double precision Measured_numberDist(1,COL)
      double precision time, Measured, dMeasuredDt, dt
      integer MeasureDate(6), k, ROW1
      character*90 filename,MeasureLocation, loc
      character*90 scenario_print,inventory_print

C-----CODE-------------------------------------------------------------

c      print*, 'OBSERVER called at time: ', time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Find measurement for current time stamp from MEAS arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Istar=0
      I=2
        Measured=0.d0
        dMeasuredDt=0.d0
      DO while(Istar .eq. 0) 
       if (I .lt. ROW) then
           IF( (time .le. MEAStime(I)).and.(time .gt. MEAStime(I-1))
     &               ) THEN
             do k = 1,Ndt
                Measured=MEAS(I-k+1)+Measured
                dMeasuredDt=dMEASdt(I-k+1)+dMeasuredDt
                Measured_numberDist(1,:) = MEASdist(I-k+1,:)
             enddo
             Istar=I
           ELSE
                   
             MEasured=0.d0
             dMeasuredDt=0.d0
             Measured_numberDist(1,:) = 0.d0
             Istar=0
             I=I+1

           ENDIF
        else
          GOTO 100
        endif
      ENDDO

c      print*, 'Measured data at time index 10', MEAS(10), dMEASdt(10)
c      print*, 'Measured data for inventory #=',x, 'time',time


         Measured=Measured/Ndt
         dMeasuredDt=dMeasuredDt/Ndt
         Measured_numberDist=Measured_numberDist/Ndt


!!!!!!!!!! Linear Interpolation !!!!!!!!!!!!!!!!!
      if ( MEAS(Istar) .lt. 0.d0 ) then
         Measured = MEAS(Istar)
      elseif ( MEAS(Istar-1) .lt. 0.d0 ) then
         Measured = MEAS(Istar)
      else
      Measured = MEAS(Istar-1) + (MEAS(Istar)-MEAS(Istar-1))/
     &      (MEAStime(Istar)-MEAStime(Istar-1))*(time-MEAStime(Istar-1))
       if ( Measured .lt. 0.d0 ) Measured=0.d0
      endif


100   continue

      RETURN
      END

!====================================================================
      ! STATE2INVENTORY !!
!====================================================================

      double precision FUNCTION state2inventory(xm,xn,xg,dpnm,dpnm_edge,
     &                            dpnm_meas,density,nuc_bin,dt,x,NzM, d)

      IMPLICIT NONE
C-----INCLUDE FILES---------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      integer k, j, first, last, nuc_bin, NzM
      integer x,n !which inventory-measurement pair?
      integer d !0 = state, 1 = derivative
      double precision xm(ibins,icomp) ! particle mass
      double precision xn(ibins) !particle number
      double precision xg(icomp-1) !gas-phase mass
      double precision Z(ibins) !0-3rd moment of particle number
      double precision dpnm(ibins) !To determine nucleation mode cut-off
      double precision dpnm_edge(ibins+1) !To determine nucleation mode cut-off
      double precision dpnm_meas !To determine upper bound cut-off
      double precision density(ibins) !kg/m3 Density to calculate volume
      double precision pi, GR, dt, y
      double precision mass_lb, frac1
      double precision mass_ub, frac2
      parameter (pi=3.14)

C-----CODE-------------------------------------------------------------

      Z(:)=0d0
      IF ( mod(x,NzM) .eq. 1 ) THEN
              !0th moment [Number]
          Z(:) = xn(:)
      ELSEIF ( mod(x,NzM) .eq. 2 ) THEN
              !1st moment [diameter]
         do k=1,ibins
            Z(k) = xn(k)*dpnm(k)*1.d-3
         enddo
      ELSEIF ( mod(x,NzM) .eq. 3 ) THEN
              !2nd moment [surface area]
         do k=1,ibins
            Z(k) = xn(k)*dpnm(k)**2*(1.d-3)**2
         enddo
      ELSEIF ( mod(x,NzM) .eq. 0 ) THEN
              !3rd moment [volume]
        do k=1,ibins
        do j=1,icomp
          if (j .ne. srth2o) then 
             Z(k) = Z(k) + xm(k,j)/density(k)*(1.d6)**3
          endif
        enddo
        enddo

      ENDIF

        State2inventory=0d0

        y = x/dble(NzM)

      IF (ceiling(y) .eq. 1) THEN
cdbg      print*, '#state2inv: 3-6 nm'
      do k=nuc_bin,ibins
         if ( dpnm_edge(k) .le. 6. ) then
               State2inventory=State2inventory+Z(k)/boxvol ! ~Nucleation mode
             if( (dpnm_edge(k+1).ge.3.).and.
     &               (dpnm_edge(k).lt.3.) )then
               mass_lb=density(k)*pi/6.*3.d-9**3.
               frac1=(xk(k+1)-mass_lb)/(xk(k+1)-xk(k))
             elseif( dpnm_edge(k+1).ge.6.)then
               mass_ub=density(k)*pi/6.*6.d-9**3.
               frac2=(mass_ub-xk(k))/(xk(k+1)-xk(k))
c               State2inventory=State2inventory+Z(k)*frac2/boxvol
c             else
c               State2inventory=State2inventory+Z(k)/boxvol ! ~Nucleation mode
cdbg               print*, 'Not first bin of 3-6!', k
             endif
         endif
      enddo

      ELSEIF (ceiling(y) .eq. 2) THEN

      do k=1,ibins
        if  ( (dpnm_edge(k+1) .ge. 10.).and.
     &                (dpnm_edge(k) .le. dpnm_meas) ) then
           State2inventory=State2inventory+Z(k)/boxvol !N_10 ~Aitken mode
        endif
      enddo
      
      ELSEIF ( ceiling(y) .eq. 3) THEN
cdbg      print*, '#state2inv: >10 nm'

      do k=1,ibins
         if ( (dpnm_edge(k+1) .gt. 6.).and.
     &                   (dpnm_edge(k) .le. dpnm_meas) ) then
           State2inventory=State2inventory+Z(k)/boxvol !N_6 ~Total - nucleation
         endif
      enddo

      ELSEIF ( ceiling(y) .eq. 4) THEN
cdbg      print*, '#state2inv: >3 nm'

       do k=nuc_bin,ibins
          if (dpnm_edge(k) .le. dpnm_meas) then 
             State2inventory=State2inventory+Z(k)/boxvol !N_3 ~total number
          endif
       enddo

       ELSE

       do k=1,ibins
          if (dpnm_edge(k) .le. dpnm_meas) then
            State2inventory=State2inventory+Z(k)/boxvol !N_tot
          endif
       enddo


      ENDIF
      RETURN
      END


!====================================================================
      ! MEAS2INVENTORY !!
!====================================================================

      double precision FUNCTION Meas2inventory(y)

C-----INCLUDE FILES---------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      double precision y !measured output
      double precision alpha80 ! Extinction coefficient at 80%RH [m2/g]
      parameter (alpha80=4.06d0)

C-----CODE-------------------------------------------------------------

c     ! Old transformation [for measured AOD in Southern ocean]
c      Meas2inventory=area/alpha80*y*1d-3 !kg of total particles

      Meas2inventory=y

      RETURN
      END

!====================================================================
      ! INVENTORY2MEAS !!
!====================================================================

      double precision FUNCTION Inventory2meas(z)

C-----INCLUDE FILES---------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      double precision z !inventory - kg of dry particles
      double precision alpha80 ! Extinction coefficient at 80%RH [m2/g]
      parameter (alpha80=4.06d0)

C-----CODE-------------------------------------------------------------


      Inventory2meas=z

      RETURN
      END

!====================================================================
      ! PARAMETER_ESTIMATE !!
!====================================================================

      SUBROUTINE parameter_estimate(time,FDBKp,FDBKz,Np,Nz,Z,Z_old,
     &                             F_old,F, G_old,G,tau,dt,Ndt,mu_old,
     &                             xmeas,NumMomt,Ztot,INDemis,INDnuc,
     &                             INDsoa,spinup,NumZ,NumP,indx_include,
     &                             data1Array,data1dt,data2Array,
     &                             Nbin_meas,ROW,    mu, cond)

      IMPLICIT NONE

C-----INCLUDE FILES---------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      INTEGER i, j, k, Np, Nz, n, m, FDBKp(Np),FDBKz(Nz),xmeas
      integer INDemis, INDnuc, INDsoa, Ndt, NumZ, NumP, ROW, Nbin_meas
      double precision F(Nz,1), F_old(Nz,1), spinup 
      double precision G(Nz,NumP), G_old(Nz,NumP), Z(Nz,1)
      double precision Z_old(Nz,1) ! inventory after calling TOMAS with
c                                  !  nominal parameters 
      double precision mu_old(Np,1), mu(Np,1) ! parameter estimate
      double precision, allocatable :: Xm(:,:)
      double precision data1Array(ROW,Nz+1),data2Array(ROW,Nbin_meas)
      double precision data1dt(ROW,Nz)
      integer ipiv(NumP)
      double precision A_o(NumZ,NumP),D(NumZ,NumP),E(NumZ,1)
      double precision A_u(NumZ-1,NumP-1), b_u,Ainit(NumZ,NumP)
c      double precision D(NumZ,NumP), E(NumZ,1),Ainit(NumZ,NumP)
c      double precision, allocatable:: A(:,:),Ainit(:,:)
      integer NumMomt
      double precision Ztot(NumMomt,1), Normalize(Nz)
      integer error, PBIO_on

      integer outflg ! 0 indicates that there are no measurements here
      integer ind1, ind2
      double precision time, Zm(Nz,1), ZmDot(Nz,1), Ym, YmDot
      double precision tau(Nz) ! 1/(proportional gain) [hr]
      double precision dt  ! time step [seconds]
      double precision Gnonzero ! check if all diagonal elements are zero
      double precision Kc(Nz,Nz) ! proportional gain [1/s]
      double precision b1(NumZ,1), b2_o(NumZ,1),b2_u(NumZ-1,1)
      double precision RHSinit(NumZ,1),RHS_o(NumZ,1),RHS_u(NumZ-1,1)
c      double precision, allocatable:: b1(:,:), b2(:,:)
c      double precision,allocatable:: RHSinit(:,:),RHS(:,:)
      integer info, lwork, iwork(NumZ)
      double precision  rcond, anorm, dummy(NumP), cond
      double precision, allocatable:: work(:)
      double precision Adummy(NumZ,NumP),Ainv(NumP,NumZ),rga(NumZ,NumP)
      double precision s(NumP), v(NumP,NumP),u(NumZ,NumZ)

      integer indx_include(NumP), nn, mm, NumMeasd, negFlg(NumZ),Flg_neg
      double precision Measd, Meas(NumZ), Measrate(NumZ)

c      double precision, allocatable :: work(:)
      character*90 string, fullname_soacond, fname1
      character*90 runname, datname, diagname, filename
      character*90 fullname_coag, fullname_cond, fullname_param
      character*90 fullname_NucRate, fullname_measured,fullname_pbio
      common /FLNAMES/ datname, diagname, fullname_coag, fullname_cond,
     &               fullname_param,fullname_NucRate,fullname_pbio,
     &               fullname_measured, fullname_soacond, fname1

!========= EXTERNAL FUNCTIONS  ==============
      double precision meas2inventory
      double precision state2inventory
      external measurement_data, meas2inventory, get_affine_dynamics
      external get_full_dynamics, state2inventory
c      external dgesv
C-----CODE-------------------------------------------------------------

      cond=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get measured output (AOD in this case)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      PBIO_on = 0
      do n=1,Np
        PBIO_on = PBIO_on+FDBKp(n)
      enddo
      IF (PBIO_on .ge. 1) THEN

      ALLOCATE( Xm(1,Nbin_meas) )
      Kc(:,:)=0.d0
c      print*, 'Observer # eqns', NumP
      nn=0
      do n=1,Nz
c       if (FDBKz(n) .eq. 1) then
c        nn=nn+1
        call measurement_data(time-spinup,dt,Ndt,Nbin_meas,ROW,
     &                       data1Array(:,1),data1Array(:,n+1),
     &                       data1dt(:,n),data2Array,outflg,Ym,Xm,YmDot)
c        if (outflg .eq. 0) FDBKz(n)=0
        Zm(n,1)=meas2inventory(Ym) 
        ZmDot(n,1)=meas2inventory(YmDot) 

      enddo

      do n=1,Nz
        Kc(n,n)=1/(3600.*tau(n))
      enddo
c      print*, '##PARAM EST: Feedback off = 0?', FDBK

      !Write measurements to txt file!
13    format(E15.6, A2, *(E9.3,:,",") )
      open(unit=40, file=fullname_measured,status='old',access='append')
      write(40,13) time, ',', Xm, (Zm(n,1),n=1,Nz)
      close(40)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate parameter estimate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( NumP .gt. 0 ) then
      Gnonzero=0.d0
      do n = 1,NumP
        Gnonzero = Gnonzero + abs( G(n,n) )
      enddo

      j=0
      do k=1,Nz
       if (FDBKz(k).eq.1) then
          j=j+1
          Meas(j)=Zm(k,1)
          Measrate(j)=ZmDot(k,1)
          if (Zm(k,1).lt.0.d0) then
              Zm(k,:) = 0.d0
              ZmDot(k,:) = 0.d0
              Meas(j) = 0.d0
              if (time .gt. 7.3d2) then
                mu = mu_old
              else
                mu(:,1) = 1.d0
              endif
              FDBKp(:)=0
          endif
       endif
      enddo

      print*, '#PARAMEST meas & rate= ',Meas, Measrate
      NumMeasd=0
      Measd=0.d0
      do k=1,Np
        Measd=Measd+Meas(k)
        NumMeasd=NumMeasd+FDBKp(k)
      enddo
      IF( time .lt. 7.3d2) then
              mu(:,:)=1d0
      elseif( (Measd .lt. 0.d0 ).or.(NumMeasd.lt.Np) ) then
              Zm(:,:) = 0.d0
              ZmDot(:,:) = 0.d0
              Xm(:,:) = 0.d0
c              mu = mu_old
              mu(:,1)=1.d0
              print*, '#paramEst: Zm ~0'
      else

        if( Gnonzero  .eq. 0d0) then
                mu(:,:)=1d0
                print*, '#paramEst: sensitivity =0'
        else

        ! Nondimensionalize!!
        do n=1,Nz
          m=n
            if ( mod(m,NumMomt) .eq. 1) then
               Normalize(n)=1/Ztot(1,1)
            elseif ( mod(m,NumMomt) .eq. 2) then
               Normalize(n)=1/Ztot(2,1)
            elseif ( mod(m,NumMomt) .eq. 3) then
               Normalize(n)=1/Ztot(3,1)
            else
               Normalize(n)=1/Ztot(4,1)
            endif
        enddo

        RHS_o(:,1)=0.d0
        b2_o(:,1)=0.d0
        A_o(:,:)=0.d0
        ind1=0
        ind2=0
        do n=1,Nz
         if (FDBKz(n).eq.1) then
         ind1=ind1+1
         b2_o(ind1,1)=0.d0
         RHS_o(ind1,1)=0.d0
         do m=1,Nz
          if (FDBKz(m).eq.1) then
           ind2=ind2+1
           b2_o(ind1,1)=b2_o(ind1,1)+Kc(n,m)*(Z(m,1)-Zm(m,1))
          endif
         enddo
         A_o(ind1,:)=G(n,:)*Normalize(n)
         ind2=0
         RHS_o(ind1,1)=(ZmDot(n,1)-b2_o(ind1,1)-F(n,1))*Normalize(n)
        endif
        enddo
        Ainit=A_o
        RHSinit=RHS_o
        lwork = MAX(1,3*MIN(NumZ,NumP)+MAX(NumZ,NumP), 5*MIN(NumZ,NumP))
        Allocate( work(lwork) )
		
        !! Calculate condition number of Normalized matrix
         ! First, normalize by transpose of its inverse
         Adummy=A_o
         call dgesvd('A', 'A', NumZ, NumP, Adummy,NumZ,s,u,NumZ,v,NumP,
     &                 work, lwork, info)
         Ainv(:,:)=0.d0
         do j=1,NumP
          do k = 1,NumZ
           do i = 1,NumP
            Ainv(j,k)=Ainv(j,k)+v(i,j)/s(i)*u(k,i)
           enddo !sum
          enddo !columns
         enddo !rows
         rga(:,:)=0.d0
         ! find if any values negati
         negFlg(:)=0
         do k=1,NumZ
          do j=1,NumP
           rga(k,j) = A_o(k,j)*Ainv(j,k)
           if(rga(k,j).lt.0.d0) negFlg(k)=1
          enddo
         enddo
         ! calculate 1-norm
         anorm=0.d0
         dummy(:)=0.d0
         do n=1,NumP
           do m=1,NumZ
             dummy(n) = ABS( rga(m,n) ) + dummy(n)
cdbg             print*, 'debug..', Ainit(m,n), abs( Ainit(m,n) ), dummy
           enddo
         enddo
         anorm = MAXVAL( dummy )
cdbg         print* ,'PARAMEST debug calc anorm', dummy
        Deallocate( work )
        Allocate( work(4*NumP) )
         ! rcond = reciprocal of condition number
         call dgecon('1',NumP,rga,NumZ,anorm,rcond,work,iwork,info)
         ! condition number
         if (time .lt. spinup) then
            cond=0.d0
         else
            cond=1/rcond
         endif
         Flg_neg=0
         do n=1,NumZ
           Flg_neg = Flg_neg+negFlg(n)
         enddo

         if ( cond .gt. 5.d0 ) then
            FDBKz(16)=0
            FDBKp(INDsoa)=0
c        !use previous estimate of 3rd parameter
        mu(INDsoa,1)=mu_old(INDsoa,1)
        
         ! Parameter Estimate Equation!! :
c        print*, 'PARAM EST normalization=', Normalize
        RHS_u(:,1)=0.d0
        b2_u(:,1)=0.d0
        A_u(:,:)=0.d0
        ind1=0
        ind2=0
        do n=1,Nz
         if (FDBKz(n).eq.1) then
         ind1=ind1+1
         b2_u(ind1,1)=G(n,INDsoa)*mu(INDsoa,1) !0.d0
         RHS_u(ind1,1)=0.d0
         do m=1,Nz
          if (FDBKz(m).eq.1) then
           ind2=ind2+1
           b2_u(ind1,1)=b2_u(ind1,1)+Kc(n,m)*(Z(m,1)-Zm(m,1))
          endif
         enddo
         ind2=0
         do m=1,NumP
          if (FDBKp(m).eq.1) then
             ind2=ind2+1
             A_u(ind1,ind2)=G(n,m)*Normalize(n)
          endif
         enddo
         RHS_u(ind1,1)=(ZmDot(n,1)-b2_u(ind1,1)-F(n,1))*Normalize(n)
        endif
        enddo

c        print*, 'Right hand side=', RHS
        ! Solve for mu : G*mu = b
        info = 1
         !!non-square system:
        lwork= NumP-1 + (NumP-1)*32
        Deallocate( work )
        Allocate( work(lwork) )
        call dgels('N',NumZ-1,Np-1,1,A_u,NumZ-1,RHS_u,NumZ-1,work,lwork,
     &                  info)
       ! Calculated 2 of the 3 parameters
           m=0 
           do n=1,Np
              if (FDBKp(n).eq.1) then
                m=m+1      
                mu(n,1) = RHS_u(m,1)
              endif
            enddo
            print*, '##PARAMEST: 2 Eqns Solved?', info, 'solution', mu

        else !conditional on if we removed a row
         !!non-square system:
        lwork= NumP + NumP*32
        Deallocate( work )
        Allocate( work(lwork) )
        call dgels('N',NumZ,Np,1,A_o,NumZ,RHS_o,NumZ,work,lwork,info)
           do n=1,Np
                mu(n,1) = RHS_o(n,1)
            enddo
            print*, '##PARAMEST: 3 Eqns Solved?', info, 'solution', mu

        endif !row was removed from set of measurements
        endif !Sensitivity > 0
      endif ! time>730sec & Measmts > 0

      ELSE
              mu(:,:)=1.d0
              Zm(:,1)=0.d0
              ZmDot(:,1)=0.d0
      ENDIF ! number parameters estimating >0
      else
              mu(:,:)=1.d0
              Zm(:,1)=0.d0
              ZmDot(:,1)=0.d0
      endif !PBIO_on > 1
c      print*, '##PARAM EST : time (s)', time, 'parameters=',mu
Cdbg      print*, '##PARAM EST end, OC density', ocden

7     format(A13,3(E15.6,","),A3,3(E15.6,","),A3,
     &          3(E15.6,","),A3)
6     format(A13,3(E15.6,","),A3,A13,3(E15.6,","),A3,A13,
     &          3(E15.6,","),A3)
5     format(A13,E15.6,A3,A13,3(E15.6,","),A3,A13,
     &          3(E15.6,","),A3)
8     format(A13,3(E15.6,","),A3,A13,3(E15.6,","),A3)
4     format(A9)
3     format(A13,E15.6,A3,A13,E15.6,A3,A13,E15.6,A3)
9     format(A16,I2,A3)
      open(unit=51, file=fullname_pbio,status='old',access='append')
c      write(51,4) 'i=i+1;'
      if ( NumP .gt. 0 ) then
      write(string, '("(A17,",I3,"(I3,","),A3)")') NumP
      write(51,string) 'Meas_used(i,:)=[', indx_include, '];'
      endif
      write(string, '("(A17,",I3,"(E15.6,","),A3)")') NumZ
c      write(string, '("(A17,",I3,"(E15.6,","),A3)")') NumP
      write(51,string) 'RHS_pbio(i,:)=[', RHSinit, '];'
      write(string, '("(A17,",I3,"(E15.6,","),A3)")') NumP
      write(51,string) 'mu_pbio(i,:)=[', mu, '];'
c      write(string, '("(A13,",I3,"(E15.6,","),A3)")') NumZ
c      write(51,string) 'RS(i,:)=[', rowsum, '];'
      write(string, '("(A13,E15.6,A3,A13,",I3,"(E15.6,","),A3,A13,",
     &               I3,"(E15.6,","),A3)")') Nz,Nz
      write(51,string)'t(i)=',time/3600.,';','Z_old(i,:)=[',Z_old,"]';",
     &     'Z(i,:)=[',Z,"]';"
      write(string, '("(A13,",I3,"(E15.6,","),A3,A13,",
     &               I3,"(E15.6,","),A3)")') Nz,Nz
      write(51,string) 'Zm(i,:)=[',Zm,"]';",'dZmdt(i,:)=[',ZmDot,"]';"!,
c     &     'Kc(i,:)=[',tau,"]';"
      write(string,'("(A13,",I3,"(E15.6,","),A3,A13,",I3,
     &             "(E15.6,","),A3,A13,",I3,"(E15.6,","),A3)")')Np,Nz,Nz
      write(51,string) 'mu_old(i,:)=[',mu_old,"]';",'F(i,:)=[',F,"]';",
     &     'F_old(i,:)=[',F_old,"]';"
      write(string,'("(A2,",I3,"(E15.6,","),A3)")')Nz
c      print*, string
      write(51,*) 'G{i}=[ ...'
      do i=1,Np
        write(51,string) ' ', G(:,i),'...'
        if (i .ne. Np) write(51,*)'; ...'
      enddo
      write(51,*) '];'
      write(string,'("(A2,",I3,"(E15.6,","),A4)")')NumZ
      write(51,*) 'G1_pbio{i}=[ ...'
      do i=1,NumP
        write(51,string) ' ', Ainit(:,i),' ...'
        if (i .ne. NumP) write(51,*)'; ...'
      enddo
      write(51,*) '];'

c      write(51,*) 'G_old{i}=[...'
c      do i=1,Np
c        write(51,string)' ',  G_old(:,i),'...'
c        if (i .ne. Np)    write(51,*)'; ...'
c      enddo
c      write(51,*) '];'
      write(51,*) 'Ztot(i,:)=[',(Ztot(k,1),k=1,NumMomt), '];'
      close(51)

      RETURN
      END 

