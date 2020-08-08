C**** BT263box.COM    Common Block for sulfur (QUS)           1/97
C****        Everything double precision*8

C--------------------------------------------------------------------

C----IMPORTANT INTEGER CONSTANTS-------------------------------------

       !constants that have to do with the number of tracers
C    NBS is the number of bulk species (gases and aerosols that don't
C	 have size resolution such as MSA)
C    NAP is the number of size-resolved "prognostic" aerosol species
C	 (ones that undergo transport)
C    NAD is the number of size-resolved "diagnostic" aerosol species
C	 (ones that don't undergo transport or have a complete budget
C	 such as aerosol water and nitrate)
C    NSPECIES is the number of unique chemical species not counting
C	 size-resolved species more than once (the sum of NBS, NAP,
C	 and NAD)
C    NBINS is the number of bins used to resolve the size distribution
C    NTM is the total number of tracer concentrations that the model
C	 tracks.  This counts bulk species, and both prognostic and 
C	 diagnostic aerosols.  Each size-resolved aerosol has a number
C	 of tracers equal to NBINS to resolve its mass distribution.
C	 An additional NBINS are required to resolve the aerosol number
C	 distribution.
C    NTT is the total number of transported tracers.  This is the same
C	 as NTM, but excludes the "diagnostic" aerosol species which
C	 do not undergo transport - ??? will I use this

       !constants that determine the size of diagnostic arrays
C    NXP is the number of transport processes that are tracked
C	 separately
C    NCR is the number of chemical reactions for which data is saved
C    NOPT is the number of aerosol optical properties tracked
C    NFOR is the number of different forcings that are calculated
C    NAERO is the number of aerosol microphysics diagnostics
C    NCONS is the number of conservation quantity diagnostics
C    NUCB is the number of nucleation pseudo-bins

      integer NBS, NAP, NAD, NSPECIES, NBINS, NTM, NTT, 
     &        NXP, NCR, NOPT, NFOR, NAERO, NCONS, NUCB

      PARAMETER (NBS=7, NAP=1, NAD=2,  NBINS=30, NXP=7, NCR=8, 
     &           NOPT=9, NFOR=5, NAERO=6, NCONS=11, NUCB=1,
     &           NSPECIES=NBS+NAP+NAD, 
     &           NTM=NBS+(NBINS+NUCB)*(NAP+NAD+1), 
     &           NTT=NBS+NBINS*(NAP+1))

C    KTACC is the total number of elements in all diagnostic arrays.
C	 This is necessary to give the DIAGT array the correct size
C	 (see below).

      PARAMETER (KTACC = IM*JM*(LM+6)*NTM + JM*LM*11*NTM + IM*JM*NTM 
     &          + JM*LM*9 + JM*NTM*NCONS 
cpja &          + IM*JM*LM*NXP*NTM
cpja &          + IM*JM*LM*NCR + IM*JM*LM*NOPT 
     &          + IM*JM*NFOR
     &          + IM*JM*LM*NTM*NAERO)

C    The following relate chemical species to the appropriate index
C    in the T0M array.  In the case of size-resolved aerosols, the
C    index of the first bin is used.  Bulk species are listed first,
C    followed by a number of tracers for the aerosol number distribution.
C    Prognostic size-resolved aerosols come next, and diagnostic size-
C    resolved aerosols are last.  This ordering is somewhat important
C    because GCM routines that handle transport (advection, convection,
C    etc...) only loop over the first NTT tracers.

      PARAMETER ( IDTH2O2 = 1,                 !Bulk species
     *            IDTSO2  = 2,
     *            IDTMSA  = 3,
     *            IDTDMS  = 4,
     &            IDTH2SO4= 5,
     &            IDTNH3G = 6,
     &            IDTNH4A = 7,
     &            IDTNUC  = 8, ! nucleation number plus mass of each species
     &            IDTNUMD = IDTNUC+NUCB*(NAP+NAD+1),  !12 NBINS for number distribution
     *            IDTSO4  = IDTNUMD+NBINS,     !42;NBINS for sulfate mass dist.
     *            IDTNA   = IDTSO4+NBINS,      !72
     &            IDTH2O  = IDTNH4+NBINS)      !102

C----TRACER CONCENTRATIONS-------------------------------------------

C    The following are the moments of the tracer distributions used
C    for advection.

      REAL*8 T0M,TXM,TYM,TZM,TXXM,TYYM,TZZM,TXYM,TZXM,TYZM

C    The TRAC array is equivalent to the collection of ten moments
C    listed above and is used for input/output.

      REAL*8 TRAC(IM,JM,LM,NTM*10)

C----DIAGNOSTIC ARRAYS-----------------------------------------------

C    See the routines contained in the diag.f file for more details
C    on the following diagnostic arrays.

C    TAIJN contains lat-long (IJ) maps of tracer (N) concentrations, 
C	 deposition and so on.
C    TAJLN contains zonal (JL) average diagnostics for each tracer (N)
C    TCNSRV tracks conservation of tracer mass and is used to
C	 generate budgets.
C    T3DX is an array of 3D transport diagnostics (disabled)
C    T3DC is an array of 3D chemistry diagnostics
C    OPTP constains aerosol optical properties
C    RADF tracks aerosol radiative forcings (W/m2)

      REAL*8 TAIJN, TAJLN, TAIJS, TAJLS, TCNSRV, RADF
cpja  REAL*8 T3DX, T3DC, OPTP
      double precision AEROD

C    DIAGT is an array that is equivalent to all the diagnostic arrays
C	 listed above.  It gets used for input and output.

      REAL*8 DIAGT(KTACC)

C----TRACER COMMON BLOCK---------------------------------------------

      COMMON /TRACER/ T0M(IM,JM,LM,NTM),TXM(IM,JM,LM,NTM),
     *                TYM(IM,JM,LM,NTM),TZM(IM,JM,LM,NTM),
     *  TXXM(IM,JM,LM,NTM),TYYM(IM,JM,LM,NTM),TZZM(IM,JM,LM,NTM),
     *  TXYM(IM,JM,LM,NTM),TZXM(IM,JM,LM,NTM),TYZM(IM,JM,LM,NTM),
     *  TAIJN(IM,JM,LM+6,NTM),TAJLN(JM,LM,11,NTM),TAIJS(IM,JM,NTM),
     *  TAJLS(JM,LM,9),TCNSRV(JM,NTM,NCONS),
cpja *  T3DX(IM,JM,LM,NXP,NTM), 
cpja &  T3DC(IM,JM,LM,NCR), 
cpja &  OPTP(IM,JM,LM,NOPT), 
     &  RADF(IM,JM,NFOR),
     &  AEROD(JM,LM,NTM,NAERO)

      EQUIVALENCE (TRAC,T0M),(DIAGT,TAIJN)

C----AEROSOL CHEMISTRY FIELDS----------------------------------------

C    The following arrays contain data needed to do sulfur chemistry

      COMMON /SCHEM/ OH(IM,JM,LM), DHO2(IM,JM,LM),PERJ(IM,JM,LM),
     &               TNO3(IM,JM,7,12), DMSS(IM,JM,6,12),ICHEMI,JJMON

C    Input nitric acid fields from Harvard.  Assumed to be the total
C    of gas and aerosol nitrate.

      REAL*8 HNO3

C    The following arrays store emission fields

Cpja  NH3OCN contains the oceanic emissions of ammonia for the
Cpja  current month.  NH3AGRO is the sum of all agricultural emissions.
Cpja  These are scaled by number of daylight hours.  NH3SOILS is
Cpja  an array containing soil emissions info.  It is also scaled
Cpja  by number of sunlight hours.

Cpja  DSK(JM) contains time of sunset at the given latitude.  It is
Cpja  used to calculate the number of sunlight hours in a day, and
Cpja  to scale certain ammonia emissions (crops, fertilizers, and
Cpja  domestic animals) accordingly.

      REAL*8 NH3OCN, NH3AGRO, NH3SOIL, NH3SUN, DMSOCN, DSK
      INTEGER NH3MONTH, DMSMONTH

C    These arrays contain a lookup table for aerosol scattering (BNO) 
C    as a function of refractive index (BNORI) and mass ratio (BNOMR).
C    Mass ratio is defined as wet aerosol mass to dry aerosol mass
C    and is used to parameterize the size distribution.  BNO has
C    units of m2 per particle.

      REAL*8 BNO,BNOMR,BNORI

C----AEROSOL COMMON BLOCK--------------------------------------------

      COMMON /AEROSOL/ NH3OCN(IM,JM), NH3AGRO(IM,JM), NH3SOIL(IM,JM), 
     &                 NH3SUN(IM,JM), DMSOCN(IM,JM), DSK(JM),
     &  HNO3(IM,JM,LM,12), BNO(141,101), BNOMR(141), BNORI(101)


C----SPECIES NAMES AND PROPERTIES------------------------------------

C     SNAME contains the name of each unique chemical species.  These
C     get mapped to the TNAME array with bin numbers in the case of
C     size-resolved aerosols.

      CHARACTER*10 SNAME(NSPECIES+NUCB+1), TNAME

      DATA SNAME/
     1           'H2O2    ',
     2           'SO2     ',
     3           'MSA     ',
     4           'DMS     ',
     5           'H2SO4   ',
     6           'NH3gas  ',
	6           'NH4aero ',
     6           'NUC     ',
     6           'ANUM__00',
     7           'ASO4__00',
     8           'ANA___00',
     9           'AH2O__00'/
     
C     TCMASS contains molecular weights
      REAL*8 TCMASS(NTM)
      DATA TCMASS /34., 64., 96., 62., 98., 17., 18.,
     &             1.0, 96., 58.45, 18.,
     &             30*1.0,30*96.,30*58.45,30*18.,30*18./

C     These arrays contains properties used for wet deposition

C     NSPTY indicates whether the species is not scavenged (TTNSCA),
C	 a soluble gas (TTSGAS), or an aerosol (TTAERO).  These get
C	 mapped to the NTRTY array.

      INTEGER NSPTY(NBS+NAP+NAD+1), NTRTY

      INTEGER TTNSCA, TTAERO, TTSGAS
      PARAMETER ( TTNSCA=0,        !not scavenged
     *            TTAERO=1,        !aerosol (infinitely soluble)
     *            TTSGAS=2)        !soluble gas (Henry's law)

      DATA NSPTY /TTSGAS, TTSGAS, TTAERO, TTNSCA, TTAERO, 
     *            TTSGAS, TTAERO, TTAERO, TTAERO, TTAERO,
     *            TTAERO/

C     Note that H2SO4 is being treated as an aerosol (infinitely
C     soluble) for purposes of wet deposition.

C     SPRKD is the Henry's law coefficient for the given species, and
C     SPDHD is the heat of dissolution.  These are only used for
C     soluble gases - set others to zero.  These get mapped to RKD
C     and DHD, respectively.

      REAL*8  SPRKD(NBS), RKD, SPDHD(NBS), DHD

      DATA SPRKD /7.4E4, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0/          !M/atm
      DATA SPDHD /-13.2, -6.27, 0.0, 0.0, 0.0, 0.0, 0.0/        !kcal/mol

Cjrp  Note, NH3 doesn't dissolve yet

C----SPECIES COMMON BLOCK--------------------------------------------

      COMMON /SPECIES/ RKD(NTM), DHD(NTM), NTRTY(NTM), TNAME(NTM)


