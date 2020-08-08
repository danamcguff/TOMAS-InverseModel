
C     **************************************************
C     *  getNucConc                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, June 2007

C     This subroutine calculates the concentration of the nucleation mode
C     at the end of a timestep.  Three processes are involved: 1) Production
C     of new nucleation mode particles, 2) coagulation with larger particles
C     3) modal self coagulation.  The differential equation takes the following
C     form...

C     dNnuc/dt = fn - Kcoag*Nnuc - Kself*Nnuc^2

C     We use a 4th order runge kutta solver that determines its own timestep
C     based on the initial dNnuc/dt and the Nnuc concentration when dNnuc/dt
C     is assumed to be zero (steady state).  The timestep is then Nnuc/(dNnuc/dt).
C     This time step is stable for this system and gives good solutions using
C     the fourth order runge kutta scheme.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nnuci - number of nucleation size particles per size bin in grid cell
C     Mnuci - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     fn - the nucleation production rate [cm-3 s-1]
C     Kcoag - the timescale for coagulation with larger particles [s-1]
C     d1 - diameter of nuc particle [m]
C     mp - average mass of nuc particles [kg]
C     massnuc - mass of a new nucleation cluster [kg]
C     spec - species that is nucleating
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

C     Nnucf - final number of nucleation size particles per size bin in grid cell
C     mass_to_add - mass of each speciesto add to other bins in order to keep 
C         conservation of mass (kg/grid cell
C     grow_to_bin1 - mass of particles that grow through self coagulation into bin1 [kg/box]

      SUBROUTINE getNucConc(Nnuci,Mnuci,fn,Kcoag,d1,mp,massnuc,spec,
     &     dt,Nnucf,Mnucf,mass_to_add,grow_to_bin1)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Nnuci, Nnucf, Mnuci(icomp), Mnucf(icomp)
      double precision fn, Kcoag
      double precision d1
      double precision mp
      double precision massnuc
      double precision dt
      double precision mass_to_add(icomp)
      integer spec
      double precision grow_to_bin1(icomp)

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision Kself    !coagulation kernal for self coagulation [cm3 s-1]
      double precision Nnuco    !number of nucleation mode particles [cm-3]
      double precision Nnuc_init !!number of nucleation mode particles [cm-3]
      double precision Nnucss   !steady state concentration (if time goes to infinity) [cm-3]
      double precision der      !dNnuc/dt [cm-3 s-1]
      double precision maxtstep !maximum time step [s]
      integer nsteps            !number of timesteps to take
      double precision tstep    !actual timstep [s]
      double precision Ntmp,Mtmp     !temperary concentration for RK
      double precision k1,k2,k3,k4 !RK parameters for number
      double precision l1,l2,l3,l4 !RK parameters for mass
      double precision m1,m2,m3,m4 !RK parameters for mass grown into bin 1
      double precision mfp, mu ! mean free path and viscosity
      double precision Dk1,ck1
      double precision Kn       !Knudson number
      double precision pi, R, kB
      double precision beta     !non-continuum correction factor
      double precision mass_move !mass that coagulated
c      double precision mass_to_nuc !mass that coagulated that stays in nuc
      double precision tot_mass, frac_mass(icomp) !total mass and fraction of each species to total mass
      double precision corr_mass !used in correction nuc bin out of bounds
      double precision Mnuc_tot ! the total mass in nucleation mode [kg/cm3]
      double precision Mnuc_toti ! the initial total mass in nucleation mode [kg/cm3]
      double precision grow     !total mass that grows into bin 1 [kg/cm3]
      double precision bin1     !geometric mean mass in bin 1
      double precision mass_part
      double precision fudge ! fudge factor for self coagulation of pseudo bin
                             ! this is necessary because the self coagulation is faster
                             ! when the particles in the nucleation bin are not mono-disperse

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

c      parameter(Kself=2.0D-9) ! for Dp = 3 nm
      parameter(pi=3.141592654d0, R=8.314d0, kB= 1.38E-23)
      parameter(fudge=5.d0)

C-----CODE--------------------------------------------------------------

c      print*,''
c      print*,'Entering getNucConc'
c      print*,'fn',fn
c      print*,'dt',dt
c      print*,'massnuc',massnuc
c	Kcoag=Kcoag/1.0d0

      Nnuc_init = Nnuci/boxvol ! convert number in box to number cm-3
      Nnuco = Nnuc_init
      Mnuc_tot = 0.d0
      do j=1,icomp-idiag
         Mnuc_tot = Mnuc_tot + Mnuci(j)
      enddo
      Mnuc_tot = Mnuc_tot/boxvol
      Mnuc_toti = Mnuc_tot

C     get Kself
      
      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6
      Kn = 2.0*mfp/d1
      Dk1 = kB*temp/(3.0*pi*mu*d1)             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
      ck1=sqrt(8.0*kB*temp/(pi*mp)) !S&P Table 12.1      

      Kn=4.0*(2.d0*Dk1)          
     &     /(sqrt(2.d0*ck1**2)*(2.d0*d1)) !S&P eqn 12.51
      beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
                                !This is S&P eqn 12.46 with non-continuum correction, beta
      Kself=fudge*2.0*pi*(2.d0*d1)*(2.d0*Dk1)*beta*1.0d6
c      Kself = 0.d0
c      print*,'Kself',Kself
      
C     Get steady state solution (for determining timestep)
C     0 = fn - Kcoag*Nnuc - Kself*Nnuc^2
C     0 = C + Bx + Ax^2
C     x = (-b+-sqrt(b^2-4ac))/2a
C     a = -Kself, b = -Kcoag, c = fn
C     The solution with the "-" is the one we want here
c      print*,'Kcoag',Kcoag,'fn',fn
      Nnucss = (Kcoag-sqrt(Kcoag**2+4.d0*Kself*fn))/(-2.d0*Kself)
c      Nnucss = fn/Kcoag
c      print*,'Nnuco',Nnuco,'Nnucss',Nnucss

C     Get intial derivative (for determining timestep)
      der = fn - Kcoag*Nnuco - Kself*Nnuco**2
c      print*,'der',der

C     Get max timestep
      maxtstep = (Nnucss-Nnuco)/der
      if (maxtstep.lt.0.d0)then
         if (Nnucss.lt.1.d0)then
            Nnucf = 0.d0
            return
         else
            print*,'Error in getNucConc, maxtstep < 0', maxtstep
            stop
         endif
      endif

C     Get number of timesteps (round up)
      nsteps =  int(dt/maxtstep)+1
      tstep = dt/nsteps
c      print*,'nsteps',nsteps

      grow = 0.d0
      bin1 = sqrt(xk(1)*xk(2))
C     Start timestepping
      do k=1,nsteps
c         print*,'Mnuc_tot',Mnuc_tot
c         print*,'particle mass',Mnuc_tot/Nnuco
         Ntmp = Nnuco
         Mtmp = Mnuc_tot
c	   mass_part = Mtmp/Ntmp
c         print*,'FRACTION',(Mtmp/Ntmp)/xk(1)
         k1 = tstep*(fn - Kcoag*Ntmp - Kself*Ntmp**2)
c         l1 = tstep*(fn*massnuc - Kcoag*Mtmp)
         l1 = tstep*(fn*massnuc - Kcoag*Mtmp - 
c     &        Kself*(Mtmp/Ntmp)/xk(1)*Ntmp**2*(Mtmp/Ntmp))
     &        Kself*Mtmp**2/bin1)
         m1 = tstep*Kself*Mtmp**2/bin1
         Ntmp = Nnuco+k1/2.d0
         Mtmp = Mnuc_tot+l1/2.d0
c	   mass_part = Mtmp/Ntmp
         k2 = tstep*(fn - Kcoag*Ntmp - Kself*Ntmp**2)
c         l2 = tstep*(fn*massnuc - Kcoag*Mtmp)
         l2 = tstep*(fn*massnuc - Kcoag*Mtmp - 
c     &        Kself*(Mtmp/Ntmp)/xk(1)*Ntmp*Mtmp)
     &        Kself*Mtmp**2/bin1)
         m2 = tstep*Kself*Mtmp**2/bin1
         Ntmp = Nnuco+k2/2.d0
         Mtmp = Mnuc_tot+l2/2.d0
c	   mass_part = Mtmp/Ntmp
         k3 = tstep*(fn - Kcoag*Ntmp - Kself*Ntmp**2)
c         l3 = tstep*(fn*massnuc - Kcoag*Mtmp)
         l3 = tstep*(fn*massnuc - Kcoag*Mtmp -
c     &        Kself*(Mtmp/Ntmp)/xk(1)*Ntmp*Mtmp)
     &        Kself*Mtmp**2/bin1)
         m3 = tstep*Kself*Mtmp**2/bin1
         Ntmp = Nnuco+k3
         Mtmp = Mnuc_tot+l3
c	   mass_part = Mtmp/Ntmp
         k4 = tstep*(fn - Kcoag*Ntmp - Kself*Ntmp**2)
c         l4 = tstep*(fn*massnuc - Kcoag*Mtmp)
         l4 = tstep*(fn*massnuc - Kcoag*Mtmp -
c     &        Kself*(Mtmp/Ntmp)/xk(1)*Ntmp*Mtmp)
     &        Kself*Mtmp**2/bin1)
         m4 = tstep*Kself*Mtmp**2/bin1
         Nnuco = Nnuco + 1.d0/6.d0*(k1 + 2.d0*k2 + 2.d0*k3 + k4)
         Mnuc_tot = Mnuc_tot + 1.d0/6.d0*(l1 + 2.d0*l2 + 2.d0*l3 + l4)
         grow = grow + 1.d0/6.d0*
     &        (m1 + 2.d0*m2 + 2.d0*m3 + m4)
      enddo

      tot_mass = fn*massnuc*boxvol*dt
      do j=1,icomp-idiag
         tot_mass = tot_mass + Mnuci(j) ! total mass at end of step that needs to be
                                                ! distributed
      enddo
c      print*,'tot_mass',tot_mass

      do j=1,icomp-idiag
         if (j.eq.spec) then
            frac_mass(j) = (Mnuci(j)+fn*massnuc*boxvol*dt)/tot_mass 
         else
            frac_mass(j) = Mnuci(j)/tot_mass 
         endif
      enddo

      do j=1,icomp-idiag
         grow_to_bin1(j) = grow*frac_mass(j)*boxvol
c         grow_to_bin1(j) = 0.d0
      enddo
       do j=icomp-idiag+1,icomp
c         grow_to_bin1(j) = grow*frac_mass(j)*boxvol
         grow_to_bin1(j) = 0.d0
       enddo

      
      do j=1,icomp-idiag
         mass_to_add(j) = (tot_mass-Mnuc_tot*boxvol)
     &        *frac_mass(j)-grow_to_bin1(j)     !mass to add to other bins
      enddo
      do j=icomp-idiag+1,icomp
         mass_to_add(j) = 0.d0
      enddo
c      print*,'mass_to_add',mass_to_add
      
      do j=1,icomp-idiag
         Mnucf(j) = Mnuc_tot*boxvol*frac_mass(j)
      enddo
Cjrp      do j=icomp+1,idiag
Cjrp         Mnucf(j) = Mnuci(j)
Cjrp      enddo
      Nnucf = Nnuco*boxvol

      tot_mass = 0.d0
      do j=1,icomp-idiag
         tot_mass = tot_mass + Mnucf(j)
      enddo
      if (tot_mass/Nnucf.gt.xk(1)/1.01d0)then !we need to move some mass out
         corr_mass = Nnucf*xk(1)/1.01d0
         do j=1,icomp-idiag
            Mnucf(j) = Mnucf(j)*corr_mass/tot_mass
         enddo
         mass_to_add(j) = mass_to_add(j) + tot_mass - corr_mass ! spread to larger bins
      endif


c      print*,'Nnucf',Nnucf

      return
      end

