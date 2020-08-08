      PROGRAM ARRAY

      !== Input file for box model ==!
      ! this cycles through an array of values to run the program over

      include 'sizecode.COM'

      !! Declare parameters and dimension of arrays to loop over
      integer a, b, c, d, e, f, g, h, i
      integer NH3num, SO2num, Nnum, DISTnum, TAUnum, SCALEnum
      integer SCALE1num, SCALE2num, SCALE3num, INVENTORYnum, MEASnum
      integer INVENTORYtot, soa_indx, nuc_indx, emis_indx, PARAMnum
      integer NumSizes,NumMoments
      parameter(NH3num=1,SO2num=1,Nnum=1,DISTnum=1, TAUnum=1)
      parameter(SCALE1num=1,SCALE2num=1,SCALE3num=1,NumSizes=4)      
      parameter(NumMoments=4) !0-3rd moments define inventory variables
      parameter(INVENTORYnum=NumMoments*NumSizes) !number of inventory variables
      parameter(PARAMnum=3) !Number of unknown parameters
      parameter(soa_indx=3,nuc_indx=1,emis_indx=2) ! index of scaling
c          ! factor corresponding to SOA production, Nucleation, and emissions
      character*90 Site, filename

      !! Declare vectors for multiple initial conditions
      double precision NH3array(NH3num) ! NH4 ppt
      double precision SO2array(SO2num) ! SO2 ppt
      double precision Narray(Nnum)     ! Number cm-3
      double precision Dp(DISTnum), sigma(DISTnum)  ! lognormal parameters (Dp in um)
      integer FDBKp(PARAMnum) !1 turns ON  estimation of each parameter
      integer FDBKz(INVENTORYnum) !1 turns ON feedback of each inventory variable
      double precision NumDays ! number of days to simulate

      !! Declare variables for input to each model run
      double precision SO2ppt ! SO2 ppt
      double precision Number
      double precision Diam, width
      double precision Kp(INVENTORYnum) ! proportional gain [hr]
      double precision ScaleFact1(SCALE1num) ! scaling factor for parameter 1
      double precision ScaleFact2(SCALE2num) ! scaling factor for parameter 2 
      double precision ScaleFact3(SCALE3num) ! scaling factor for parameter 3

      !! Time variables
      real TotTime, timef(2), cntime
	  
      !! SET values for multiple scenarios
       DATA NH3array/1.0D3/
       DATA SO2array/1.0d3/
       DATA Narray/1.0D3/
       DATA Dp/0.1/
       DATA sigma/2.0/
       DATA FDBKp/1, 1, 1/ !nucRate, POA emiss, SOA cond
c Inventory Variables: [0-,1-,2-,3-moments]
c               !Num, Diam, SurfA, Mass
c                                    !diameter range
       DATA FDBKz/1, 0, 0, 0,        !3-6 nm
     &            1, 0, 0, 0,        !>10 nm
     &            0, 0, 0, 0,        !>6 nm
     &            0, 0, 0, 1/        !>3 nm
       DATA Kp/1.5d1,1.d1,1.d1,6.d1, !diagonals of gain matrix for
     &         1.5d1,1.d1,1.d1,6.d1, ! each inventory variables 
     &         6.d1,1.d1,1.d1,6.d1,  !  (minutes)
     &         6.d1,1.5d1,1.5d1,6.d1/
c   ! Nominal value of each scaling factor across the different
c      simulations to run
       DATA ScaleFact1/1.d0/!1.d0,1.d0,1.d0/!8.d-1,1.2d0,3.d0/!3.d0/!POA emiss
       DATA ScaleFact2/1.d0/!1.d0,1.d0,1.d0/!1.d-1,2.d0,1.d1/!1.d1/!NucRate
       DATA ScaleFact3/1.d0/!1.d0,1.d0,1.d0/!2.d-1,2.d0,5.d0/!2.d0/!SOA cond
       Kp=Kp/60.! convert gain to hours
      
      !! BOX conditions
      Site='SanPietroCapofiume' !Measurement station name
      NumDays=370.d0
      boxvol=7.5e19 !cm3
      area=1.5e11 !m2
      temp=298.15 !K
      pres=1.0e5  !Pa
      rh=0.8  !relative humidity (--)

      boxmass=0.0289*pres*boxvol*1.0e-6/8.314/temp

      i=0
      do a=1,NH3num
           do b=1,SO2num
                do c=1,Nnum
                     do d=1,DISTnum
                      do f=1,SCALE1num
                       do g=1,SCALE2num
                        do h=1,SCALE3num
                         NH3ppt_o = NH3array(a)
                         so2_print = SO2array(b)
                         SO2ppt = SO2array(b)
                         Number = Narray(c)
                         Diam = Dp(d)
                         width = sigma(d)
                  TotTime = etime(timef)
                  cntime = timef(1)
                  i=i+1
                  call box(SO2ppt,Number,Diam,width,Kp,FDBKz, 
     &             FDBKp,NumSizes,NumMoments,INVENTORYnum,PARAMnum,
     &             ScaleFact1(f),ScaleFact2(g),ScaleFact3(h),i,
     &             emis_indx,soa_indx,nuc_indx,Site,NumDays )
      TotTime = etime(timef)
      cntime = timef(1) - cntime
      open(93,file='time.txt',status='unknown',access='append')
      write(93,'(E10.3)') cntime
      close(93)
                      enddo
                      enddo
                     enddo
                    enddo
                enddo
           enddo
      enddo

       END ! of main

