
C     **************************************************
C     *  boxinit                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, July 2000

C     Initializes the box model for testing size-resolved microphysics.
C     Starts either with zero concentration or with a restart file if
C     one exists.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE boxinit(time,No,Dp,sigma)

C-----INCLUDE FILES-----------------------------------------------------

      include 'BB192SM9.COM'
      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision time
      double precision No
      double precision Dp
      double precision sigma
      double precision totsulf

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,l,n

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

 2    format(200E15.6)

C-----CODE--------------------------------------------------------------

C Initializations

      !Initialize everything to zero at first
      do i=1,IM
         do j=1,JM
            do l=1,LM
               do n=1,NTM
                  T0M(i,j,l,n)=0.0
               enddo
            enddo
         enddo
      enddo

      call speciesmap()
      call initbounds()
      call loginit(No,Dp,sigma,amsulfden,IDTSO4)

C initallize to ammonium sulfate
      do i=1,IM
         do j=1,JM
            do l=1,LM
                          totsulf = 0.d0
               do k=1,nbins
                          totsulf = totsulf + T0M(i,j,l,IDTSO4+k-1)
               enddo
                          T0M(i,j,l,IDTNH4A)=totsulf/96.*18.*2.
            enddo
         enddo
      enddo

c      do i=1,IM
c         do j=1,JM
c            do l=1,LM
c               T0M(i,j,l,IDTNUC)=10000.d0*boxvol
c               T0M(i,j,l,IDTNUC+1)=10000.d0*boxvol*1.381067932004976D-24
c            enddo
c         enddo
c      enddo  

c      print*,'test1',T0M(1,1,1,IDTNUC)
c      print*,'test1',T0M(1,1,1,IDTNUC+1)

C Read initial conditions from box.rsf if it exists
      open(unit=11,file='box.rsf',status='old',err=100)
      read(11,2) time,(T0M(1,1,1,n),n=1,NTM)
      close(11)
      time=time*3600.
      goto 200

C Otherwise, use these simple initial conditions
 100  continue
      T0M(1,1,1,IDTH2SO4)=0.
      HNO3(1,1,1,1)=100.0e-12    !mixing ratio
 200  continue

      RETURN
      END
