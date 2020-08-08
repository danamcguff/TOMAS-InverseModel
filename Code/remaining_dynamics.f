!====================================================================
      ! REMAINING_DYNAMICS !!
!====================================================================


      SUBROUTINE remaining_dynamics(Gmk,Gnk,Ggc,Mkold,Mknew,Nkold,
     &                        Nknew, Gcold,Gcnew,dt,Fmk,Fnk,Fgc)

      IMPLICIT NONE

C-----INCLUDE FILES---------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS--------------------------------------------

      INTEGER j, k
      double precision Mkold(ibins, icomp) !Mk at previous time-step
      double precision Nkold(ibins) !Nk at previous time-step
      double precision Gcold(icomp-1) !Gc at previous time-step
      double precision Mknew(ibins, icomp) !Mk at current time-step
      double precision Nknew(ibins) !Nk at current time-step
      double precision Gcnew(icomp-1) !Gc at current time-step
      ! Dynamics unrelated to input:
      double precision Fmk(ibins, icomp) ! term for particle mass
      double precision Fnk(ibins) ! term for particle number
      double precision Fgc(icomp-1) ! term for gas-phase mass

      double precision dt
!========= INTERMEDIATE VARIABLES  ==============
      double precision Gmk(ibins, icomp) !input-affine term for particle mass
      double precision Gnk(ibins) !input-affine term for particle number
      double precision Ggc(icomp-1) !input-affine term for gas-phase mass
      double precision XdotN(ibins) !full dynamics- particle number
      double precision XdotM(ibins,icomp) !full dynamics- particle mass
      double precision XdotG(icomp-1) !full dynamics- gas-phase

C-----CODE-------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get rate of change of all states (xdot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=1,ibins
        XdotN(k) = ( Nknew(k) - Nkold(k) )/dt
        do j=1,icomp
          XdotM(k,j)= (Mknew(k,j)-Mkold(k,j) )/dt
          if (j .ne. icomp ) then
            XdotG(j)= (Gcnew(j)-Gcold(j) )/dt
          endif
        enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate f, dynamics unrealated to input (u) :
      !  Xdot = f(x) + g(x,u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c      print*,'XdotG', XdotG
c      print*,'Ggc', Ggc
        Fnk=XdotN-Gnk
        Fmk=XdotM-Gmk
        Fgc=XdotG-Ggc

c        print*,'f(x) for Nk', Fnk, 'xdot for Nk', XdotN, 'g(x) for Nk',
c     &           Gnk
      RETURN
      END
