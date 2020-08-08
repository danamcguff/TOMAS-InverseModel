    
C     **************************************************
C     *  report                                        *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     Used by box model for testing aerosol microphysics code.  This
C     routine outputs the simulation status.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE report(time, CS, CoagSink, density, nucRate)

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision time, CS, CoagSink(ibins-1), density(ibins), mu
      double precision nucRate

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer n,k

      character*90 runname, datname, diagname, fname
      character*90 fullname_coag, fullname_cond, fullname_param
      character*90 fullname_NucRate, fullname_measured,fullname_pbio
      common /FLNAMES/ datname, diagname,fullname_coag, fullname_cond,
     &                fullname_param, fullname_NucRate,fullname_pbio,
     &                fullname_measured, fname

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

 1    format(A15,43E15.6)
 2    format(325E15.6)
 3    format(E15.6,I4,225E15.6)
 4    format(E15.6,A2, *(G0.4,:,",") )
 5    format(E15.6,A2, 225E15.6) 
 6    format(E15.6,A2, 225E15.6)

C-----CODE--------------------------------------------------------------

C Time

      write(*,1) 'Time (s): ', time
      write(*,1) 'Time (h): ', time/3600.

C Bulk species

      do n=1,NBS
         write(*,1) SNAME(n),T0M(1,1,1,n)
      enddo

C Aerosol number distribution

      write(*,1) 'A#: ',(T0M(1,1,1,IDTNUMD+k-1),k=1,NBINS)

C Aerosol mass distribution

      write(*,1) 'SO4: ',(T0M(1,1,1,IDTSO4+k-1),k=1,NBINS)
      write(*,1) 'NA:  ',(T0M(1,1,1,IDTNA+k-1),k=1,NBINS)
      write(*,1) 'OC:  ',(T0M(1,1,1,IDTC+k-1),k=1,NBINS)
      write(*,1) 'H2O: ',(T0M(1,1,1,IDTH2O+k-1),k=1,NBINS)

C Write everything to a single line in a separate output file
c      print*,'a'

      open(unit=30,file=datname,status='old',access='append')
      write(30,2) time/3600.,(T0M(1,1,1,n),n=1,NTM),
     &           (density(k),k=1,ibins),Gc
      close(30)
c      print*,'b'

c      open(unit=31,file=diagname,status='old',access='append')
c      do k=1,NAERO
c         write(31,3) time/3600.,k,(AEROD(1,1,n,k),n=1,NTM)
c      enddo
c      close(31)

      open(unit=32,file=fullname_coag,status='old',access='append')
      write(32,4) time/3600.,',', (CoagSink(k),k=1,IBINS-1)
      close(32)
c
      open(unit=33,file=fullname_cond,status='old',access='append')
      write(33,5) time/3600.,',',CS
      close(33)

c      open(unit=35,file=fullname_NucRate,status='old',access='append')
c      write(35,5) time/3600.,',',nucRate
c      close(35)

c      print*, 'parameter est', mu
c      open(unit=34,file=fullname_param,status='old',access='append')
c      write(34,6) time/3600.,',',mu
c      close(34)

c      print*,'c'

      RETURN
      END
