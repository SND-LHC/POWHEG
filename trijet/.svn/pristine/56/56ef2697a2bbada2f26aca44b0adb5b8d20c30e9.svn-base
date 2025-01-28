      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.true.)
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
c      PTMIN=100.
      MAXPR=2
      ptrms=2.5d0
      write(*,*)
      write(*,*) '*******************************************'
      write(*,*) '*******************************************'
      write(*,*) ' INITIAL p_T SPREADING OF ',ptrms,' GEV    '
      write(*,*) '*******************************************'
      write(*,*) '*******************************************'
      write(*,*)
      if(.not.uevent) then
         PRSOF=0d0
         write(*,*)
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) ' NO UNDERLYING EVENT  WILL BE GENERATED    '
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*)   
      else
         write(*,*)
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) ' UNDERLYING EVENT  WILL BE GENERATED     '
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*)   
      endif
c----DO NOT USE SOFT ME CORRECTION     
      SOFTME=.FALSE.
c---- Inhibit photon radiation by the shower
c this causes problems due to events that have very small total HT at the Les Houches level,
c and nevertheless shower a very energetic photon.
      VPCUT=1d30
      end

      subroutine getmaxev(numev)
      integer numev
C---  Opens input file and counts number of events, setting MAXEV
      call opencount(numev)
      end
      
      subroutine UPINIT
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      nevhep=0
c first call lhefreadhdr; this sets lprup;
      call lhefreadhdr(97)

C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      if (lprup(1).eq.10015) then
         CALL HWUSTA('TAU-    ')
         CALL HWUSTA('TAU+    ')
      endif
      end
      

      subroutine UPEVNT
      include 'HERWIG65.INC'
      logical ini,changescalup
      data ini/.true./
      save ini,changescalup
      real * 8 powheginput
      external  powheginput
      if(ini) then
         if(powheginput("#changescalup").eq.1) then
            print *,"SCALUP is changed..."
            changescalup=.true.
         else
            changescalup=.false.
         endif
         ini=.false.
      endif      
      maxer=maxev
      call lhefreadev(97)
      if(changescalup) call py_change_scalup
      end


      subroutine hwabeg
      call init_hist
      end


      subroutine hwaend
      implicit none
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      call pwhgsetout
      if(rnd_cwhichseed.eq.'none') then
         call pwhgtopout(pwgprefix(1:lprefix)//
     $                   'POWHEG+HERWIG-output')
      else
         call pwhgtopout(pwgprefix(1:lprefix)//
     $                   '-'//rnd_cwhichseed //'-'//
     $                   'POWHEG+HERWIG-output')
      endif
      end




      subroutine hwanal
      INCLUDE 'HERWIG65.INC'
      include 'LesHouches.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)

      if (ierror.ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call HWUPUP         !hepeup
            call hwuepr         !all the event
         endif
         call pwhgaccumup
         return
      endif
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end



