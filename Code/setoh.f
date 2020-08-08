
C     **************************************************
C     *  setoh                                         *
C     **************************************************

C     WRITTEN BY Peter Adams

C     Returns a typical OH concentration (molec cm-3) depending on
C     clear/cloudy conditions and time of day.

      double precision FUNCTION setoh(day,clear,cyclet)

      IMPLICIT NONE

C     INCLUDE FILES...

C-----ARGUMENT DECLARATIONS------------------------------------------

      logical day,clear
      double precision cyclet   !how far into 12-hour diurnal cycle we are (s)

C-----VARIABLE DECLARATIONS------------------------------------------

      !typical concentrations (#/cm-3) - see Seinfeld and Pandis p. 250
      double precision clearmax, cloudymax, nightconc
      double precision pi

C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------

      parameter(clearmax=4.3d6, cloudymax=8.6d5, nightconc=1.0d5)
      parameter(pi=3.141592654)

C-----CODE-----------------------------------------------------------

      if (day) then
         if (clear) then
            setoh=clearmax*sin(cyclet/(12.*3600.)*pi) + nightconc  ! added this to avoid going to 0
         else
            setoh=cloudymax*sin(cyclet/(12.*3600.)*pi) + nightconc  ! added this to avoid going to 0
         endif
      else
         setoh=nightconc
      endif
      ! Degbugging ...
      setoh=clearmax
c      setoh=cloudymax
c      print*, 'SETOH: hydroxyl conc:',setoh, 'cycletime', cyclet

      RETURN
      END
