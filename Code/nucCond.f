
C     **************************************************
C     *  nucCond                                       *
C     **************************************************

C     WRITTEN BY Jeff Pierce, June 2007

C     This subroutine does condensational growth of the particles in the nucleation
C     pseudo-bin.  It assumes that the number distribution within the pseudo bin
C     takes the following form...

C     ln(dN/dlnX) = A*lnX + B

C     Where N is the number of particles, X is the mass per particle and A and B
C     are parameters fit using the total mass and total number in the bin.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
C     Nnuci - number of nucleation size particles per size bin in grid cell
C     Mnuci(icomp) - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     mcond - mass of species to condense [kg/grid cell]
C     spec - the number of the species to condense

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Nnuci,Mnuci - same as above, but final values

      SUBROUTINE nucCond(Nki,Mki,Nnuci,Mnuci,mcond,spec,
     &     Nkf,Mkf,Nnucf,Mnucf)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nki(ibins), Mki(ibins, icomp)
      double precision Nkf(ibins), Mkf(ibins, icomp)
      double precision Nnuci, Mnuci(icomp)
      double precision Nnucf, Mnucf(icomp)
      double precision mcond
      integer spec

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision M1i, M2i ! bound of bin (before condensation growth)
      double precision M1f, M2f ! bound of bin (after condensation growth)
      double precision Mtoti,Mtotf !initial and final total masses
      double precision Mtotdryi,Mtotdryf !initial and final total masses
      double precision WR       !ratio of wet mass to dry mass
      double precision mpi,mpf !initial and final average mass per particle
      double precision mfracf(icomp) !final mass fraction
      double precision const    !constant for bin limit growth
      double precision Ap,Bp    !fit parameters of the sub-bin number distribution
      integer top_bin           !the top bin that nuc particles grew into
      integer bot_bin           !the bottom bin that nuc particles grew into
      double precision madd     !mass to add per bin [kg]
      double precision hi,lo    !limits for integration
      double precision scalef   !scale factor to not have double precision overflow
                                !converts units of mass to kg*10^23

      double precision Ntot2,Mtot2
      double precision tot_mass

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654d0, R=8.314159d0) !pi and gas constant (J/mol K)
c      parameter(M1i=1.220703125d-25,M2i=1.0d-21)
      parameter(M1i=2.d-25,M2i=1.0d-21) !M1i=2.d-25
      parameter(scalef=1.0D23)

C-----CODE--------------------------------------------------------------

      ! initialize variables
      do k=1,ibins
         Nkf(k)=Nki(k)
         do j=1,icomp
            Mkf(k,j)=Mki(k,j)
         enddo
      enddo
      Nnucf = 0.d0
      do j=1,icomp
         Mnucf(j) = 0.d0
      enddo

      !determine particle masses before and after condensation
      !and find mass fractions
      Mtoti = 0.d0
      do j=1,icomp
         Mtoti = Mtoti + Mnuci(j)
      enddo
      Mtotf = Mtoti + mcond

      do j=1,icomp
         if (j.eq.spec)then
            mfracf(j) = (Mnuci(j)+mcond)/Mtotf
         else
            mfracf(j) = Mnuci(j)/Mtotf
         endif
      enddo

       !repeat but for not including diagnostic species
      Mtoti = 0.d0
      do j=1,icomp
         Mtoti = Mtoti + Mnuci(j)
      enddo
      Mtotf = Mtoti + mcond

      Mtotdryi = 0.d0
      do j=1,icomp-idiag
         Mtotdryi = Mtotdryi + Mnuci(j)
      enddo
      Mtotdryf = Mtotdryi + mcond

      WR = Mtotf/Mtotdryi

C     determine how much bin limits grow
C     in the kinetic regime particle growth from different sizes follows...
C     m2^(1/3) - m1^(1/3) = Const
C     First determine the constant from the growth in the average particle mass
      mpi = Mtotdryi*WR/Nnuci
      mpf = Mtotdryf*WR/Nnuci
      const = mpf**(1.d0/3.d0)-mpi**(1.d0/3.d0)
C     Now determine new bin limits
      M1f = ((M1i*WR)**(1.d0/3.d0)+const)**3*scalef
      M2f = ((M2i*WR)**(1.d0/3.d0)+const)**3*scalef
      Mtotf = Mtotdryf*WR*scalef
Cjrp      M1f = M1i*(mpf/mpi)
Cjrp      M2f = M2i*(mpf/mpi)
c      print*,''
c      print*,'In nucCond.f'
c      print*,'mpi',mpi,'mpf',mpf
c      print*,'M1f',M1f,'M2f',M2f

C     Now we need to find the A and B parameters of the fit of the sub-bin
C     size distribution
      call getNucBinFit(Nnuci,Mtotf,M1f,M2f,Ap,Bp)
c      print*,'Ap',Ap,'Bp',Bp

      !check to make sure the fit worked
      Ntot2 = exp(Bp)/Ap*(M2f**Ap-M1f**Ap)
      Mtot2 = exp(Bp)/(Ap+1)*(M2f**(Ap+1)-M1f**(Ap+1))
c      print*,'check of nucleation subgrid fit'
c      print*,'Nnucf',Nnuci,'Ntot2',Ntot2
c      print*,'Mtotf',Mtotf,'Mtot2',Mtot2

C     Now we have to remap the mass into various bins
      top_bin = 0
      do while ((M2f .gt. xk(top_bin+1)*WR*scalef).and.
     &         (top_bin.lt.ibins))
         top_bin = top_bin+1
      enddo
      if (top_bin.eq.ibins) then
         tot_mass=0.d0
         do j=1,icomp-idiag
            if (j.eq.spec) then
               tot_mass=tot_mass+Mnuci(j)+mcond
            else
               tot_mass=tot_mass+Mnuci(j)
            endif
         enddo
         Nnucf=0.d0
         Nkf(ibins)=Nkf(ibins)+tot_mass/sqrt(xk(ibins)*xk(ibins+1))
         do j=1,icomp
            Mnucf(j)=0.d0
            if (j.eq.spec) then
               Mkf(ibins,j)=Mkf(ibins,j)+Mnuci(j)+mcond
            else
               Mkf(ibins,j)=Mkf(ibins,j)
            endif
         enddo
      endif
      bot_bin = 0
      do while (M1f .gt. xk(bot_bin+1)*WR*scalef)
         bot_bin = bot_bin+1
      enddo
c      print*,'top_bin',top_bin,'bot_bin',bot_bin


c      print*,'Mkf1',Mkf(1,1)
      
      if (top_bin.eq.0) then ! nothing grew out of bin
         !Nnucf = Nnucf
         do j=1,icomp
            Mnucf(j) = Mtotf*mfracf(j)/scalef
         enddo
      elseif (bot_bin.eq.0) then !the largest particles in the nucleation bin
                                 !grew out, but the smallest ones didn't
         do k=1,top_bin
            if (k.ne.top_bin) then !integrate over entire bin
               hi = xk(k+1)*WR*scalef
               lo = xk(k)*WR*scalef
               Nkf(k)=Nkf(k)+exp(Bp)/Ap*(hi**Ap-lo**Ap)
               madd= exp(Bp)/(Ap+1)*
     &              (hi**(Ap+1)-lo**(Ap+1))/scalef
               do j=1,icomp
                  Mkf(k,j)=Mkf(k,j) + madd*mfracf(j)
               enddo
            else !integrate from bottom of bin to M2f
               hi = M2f
               lo = xk(k)*WR*scalef
               Nkf(k)=Nkf(k)+exp(Bp)/Ap*(hi**Ap-lo**Ap)
               madd = exp(Bp)/(Ap+1)*
     &              (hi**(Ap+1)-lo**(Ap+1))/scalef
               do j=1,icomp
                  Mkf(k,j)=Mkf(k,j) + madd*mfracf(j)
               enddo
            endif
         enddo
         ! leave the rest in nucleation bin
         hi = xk(1)*WR*scalef
         lo = M1f
         Nnucf=exp(Bp)/Ap*(hi**Ap-lo**Ap)
         madd = exp(Bp)/(Ap+1)*(hi**(Ap+1)-lo**(Ap+1))/scalef
         do j=1,icomp
            Mnucf(j)=madd*mfracf(j)
         enddo       
      else ! all particles grew out of nucleation bin
         do k=bot_bin,top_bin
            hi = xk(k+1)*WR*scalef
            lo = M1f            
            if (k.eq.bot_bin) then !integrate from M1f to top of bin
               Nkf(k)=Nkf(k)+exp(Bp)/Ap*(hi**Ap-lo**Ap)
               madd= exp(Bp)/(Ap+1)*
     &              (hi**(Ap+1)-lo**(Ap+1))/scalef
               do j=1,icomp
                  Mkf(k,j)=Mkf(k,j) + madd*mfracf(j)
               enddo               
            elseif (k.ne.top_bin) then !integrate over entire bin
               hi = xk(k+1)*WR*scalef
               lo = xk(k)*WR*scalef              
               Nkf(k)=Nkf(k)+exp(Bp)/Ap*(hi**Ap-lo**Ap)
               madd= exp(Bp)/(Ap+1)*
     &              (hi**(Ap+1)-lo**(Ap+1))/scalef
               do j=1,icomp
                  Mkf(k,j)=Mkf(k,j) + madd*mfracf(j)
               enddo
            else !integrate from bottom of bin to M2f
               hi = M2f
               lo = xk(k)*WR*scalef               
               Nkf(k)=Nkf(k)+exp(Bp)/Ap*(hi**Ap-lo**Ap)
               madd = exp(Bp)/(Ap+1)*
     &              (hi**(Ap+1)-lo**(Ap+1))/scalef
               do j=1,icomp
                  Mkf(k,j)=Mkf(k,j) + madd*mfracf(j)
               enddo
            endif
         enddo         
      endif

c      print*,'Mkf1',Mkf(1,1)

      return
      end
