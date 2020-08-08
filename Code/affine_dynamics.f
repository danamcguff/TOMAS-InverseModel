!====================================================================
      ! AFFINE_DYNAMICS !!
!====================================================================


      SUBROUTINE affine_dynamics(Gmk, Gnk, Ggc,nuc_bin,fn, 
     &             POA1,POA2,Mff,Mbb,dpg1,dpg2,sigmag1,sigmag2,Ni1,
     &             Nf1,Ni2,Nf2,adt,FLGemis,
     &             soa_nk,soa_mk, voc_emis, FLG, INDemis, INDnuc,
     &              INDsoa,sfact, sfDim, time )

      IMPLICIT NONE

C-----INCLUDE FILES---------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      INTEGER j, k, INDemis, INDnuc, INDsoa, sfDim
      integer FLG !which process?
c    ( 1=NucRate, 2=Poa emissions, 3=SOA production )
      integer FLGemis !use data or lognormal parameters?
      double precision Gmk(ibins, icomp) !input-affine term for particle mass
      double precision Gnk(ibins) !input-affine term for particle number
      double precision Ggc(icomp-1) !input-affine term for gas-phase mass
      double precision sfact(sfDim,1) !scaling factor for each process
!========= VARIABLES RELATED TO DYNAMICS TO BE ESTIAMTED ==============
      double precision  fn, mnuc, pi
      integer nuc_bin
      parameter (pi=3.14159)
      double precision POA1,POA2,dpg1,dpg2,sigmag1,sigmag2,adt,E_dens
      double precision Mff, Mbb
      double precision Nkout1(ibins),Mkout1(ibins)
      double precision Nkout2(ibins),Mkout2(ibins)
      integer Ni1,Nf1,Ni2,Nf2
      double precision soa_nk(ibins), soa_mk(ibins,icomp), condensed
      double precision voc_emis, time
      double precision Nkout(ibins),Mkout(ibins,icomp)
c      double precision Essa(ibins) !sea salt emissions particle number
c      double precision u10 ! 10m-altitude wind speed [m/s]
c      double precision focean ! ocean cover fraction

C-----CODE-------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get (time-varying) portion of dynamics affine wrt input
      ! (parameter estimate)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Initialize Output to zero
      do j=1,icomp
        if( j .ne. icomp) then
         Ggc(j)=0d0
        endif
         do k=1,ibins
            Gmk(k,j)=0d0
            Gnk(k)=0d0
         enddo
      enddo

      IF(FLG.eq. INDemis ) THEN
              call logemiss(POA1,Mff,Mbb,dpg1,sigmag1,Ni1,
     &                        Nf1,ocden,srtc,FLGemis,Nkout1,Mkout1)
              do k=1,ibins
                Gnk(k)=Nkout1(k)/adt
                Gmk(k,srtc)=Mkout1(k)/adt
              enddo
       

      ELSEIF(FLG .eq. INDnuc ) THEN
              mnuc = ( pi/6.*(3.*1d-9)**3)*1350.d0
              Gnk(nuc_bin)=fn*boxvol
              Gmk(nuc_bin,srtso4)=fn*mnuc*boxvol

      ELSEIF(FLG .eq. INDsoa ) THEN
              condensed=0.d0
              do k=1,ibins
                Gnk(k)=soa_nk(k)
                do j=1,icomp
                  Gmk(k,j)=soa_mk(k,j)
                  if (j .eq. srtc ) then
                    condensed=condensed+soa_mk(k,j)
                  endif
                enddo
              enddo
              Ggc(srtc)= -condensed+voc_emis

      ENDIF        

      RETURN
      END
