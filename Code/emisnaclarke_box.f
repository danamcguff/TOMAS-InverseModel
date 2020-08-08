
C     **************************************************
C     *  emisnaclarke_box.f                            *
C     **************************************************

C     WRITTEN BY Jeffrey Pierce, April 2004
C     MODIFIED for box model, September 2005

C     Implements emissions of sea salt based on:

C     Parameters taken from Owens, S.; Clarke, A.D., "Marine aerosol
C     sources part 1: Particle fluxes derived from breaking waves", 
C     submitted to GRL, April 2004

C     This is the same subroutine as emisnaN2clarke.f except that a
C     few changes have been made for the box model.  The wind speed,
C     fraction of grid cell that is ocean, time step and box area
C     must be provided.

C     Like with the emisnamon.f file used previously, this input
C     may cause too high concentrations because this version only 
C     changes T0M, assuming that the flux mixes vertically throughout 
C     the grid cell, rather than being confined to the surface.

C     This program currently assumes an aerosol density of 2.1 g/cc
C     based on Campuzano-Jost et. al., "Near-Real-Time Measurement of
C     Sea-Salt Aerosol during the SEAS Campain", Journal of Atmospheric
C     and Oceanic Technology, 20, 2003.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE emisnaclarke(u10,focean,adt,N)

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision u10    !wind speed 10m above ocean
      double precision focean !fraction of box surface that is ocean
      double precision adt    !timestep

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k !counters
      double precision F100   !number flux at 100% whitecap coverage(m^-2 s^-1)
      double precision W      !fraction whitecap coverage based on U10
      double precision Dbin   !lower limit dimeter for each size bin
      double precision N(NBINS) !number added to bin during timestep
      double precision A      !Flux into each size bin for 100% whitecap coverage

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      dimension A(NBINS)

      data A /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &1719793.975,2887719.254,4086059.079,5222972.121,
     &6172287.155,6789341.855,6954290.435,6647842.508,6030292.47,
     &5411159.039,4920485.633,4467448.678,4379031.834,4180592.479,
     &3836983.331,3328339.218,2675909.44,1972225.823,1384692.112,
     &1062343.821,913194.1118,859176.8257,812688.43,719215.3301,
     &580735.2991,418247.5535,273217.6572,183340.5653,132174.9032,0.0/

C-----CODE--------------------------------------------------------------

c      print*, 'SeaSpray Ocean fraction=', focean
      !Loop over grid cells
      do j=1,JM
         do i=1,IM

           do k=1,NBINS
             N(k)=0d0
           enddo
            !Check if over ocean
c            focean=1.-FDATA(i,j,2)-ODATA(i,j,2)*(1.-FDATA(i,j,2))
            if (focean .gt. 0.5) then  !This assumes only gridcells that
                                       !are at least 50% water are oceans

               !in ocean area - calc wind speed/eqm conc
c               u10=BLDATA(i,j,1)
               !calculate the fraction of whitecap coverage
               W=3.84E-6*u10**(3.41)
c               print*,'EMISCLARKE Whitecap frac:',W, 'area', area,
c     &                'time step',adt,'ocean frac', focean

               !loop over bins
               do k=1,(NBINS)

                  F100 = A(k)

                  !calculate the number of particals added to bin in timestep
                  !DXYP(j) is the area of gridcell (area is a function of latitude)
                  !DT is time in timestep (s) and NDYN is the number of leapfrog
                  !timesteps per hour.
                  N(k) = F100*W*area*adt*focean

c                     T0M(i,j,1,IDTNUMD+k-1)=T0M(i,j,1,IDTNUMD+k-1)+N !Total number
c                     T0M(i,j,1,IDTNA+k-1)=T0M(i,j,1,IDTNA+k-1)
c     &                 +N*sqrt(xk(k)*xk(k+1))  !Total NA mass
cdbg                    print*,'EMISNACLARKE bin', k, 
cdbg     &                 'number emissions [1/s]',N/adt, 
cdbg     &                'Mass emissions [kg/s]',N*sqrt(xk(k)*xk(k+1))/adt
               enddo

            endif

         enddo
      enddo

      RETURN
      END
