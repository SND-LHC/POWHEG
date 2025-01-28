      subroutine setup_PYTHIA_tune
      implicit none
      real * 8 powheginput
      external powheginput
      integer tunevar
C   100       A : Rick Field's CDF Tune A                     (Oct 2002)
C   103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
C   320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)

c      tunevar = 350
c      tunevar = 340
      tunevar = 320
      call PYTUNE(tunevar)
      if (tunevar.ne.0) then
        print *,"*****************************************************"
        print *,"*                                                   *"
        print *,"*             Tune is used in PYTHIA:               *"
        if (tunevar.eq.350) then !<
        print *,"*           Perugia2011 tune is in use...           *"
        end if !>
        if (tunevar.eq.340) then !<
        print *,"*             AMBT1 tune is in use...               *"
        end if !>
        if (tunevar.eq.320) then !<
        print *,"*           Perugia2010 tune is in use...           *"
        end if !>
        print *,"*                                                   *"
        print *,"*****************************************************"
      end if
      end

      subroutine setup_PYTHIA_parameters
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      double precision parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      integer nev
      common/cnev/nev
      real * 8 powheginput
      external powheginput
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     multiple interactions
c     (MI can increase a lot the execution time)
c      if(.not.mult_inter) then
c         mstp(81)=20   !No Multiple interactions. Force a call to PYEVNW 
c      else
c         mstp(81)=21   ! MPI on in the PYEVNW MPI scenario
c      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      if(powheginput("#nohad").lt.0) then
c         write(*,*) 'SWITCHING OFF MPI!!!!!!!!!!!!!!!!!!!!!'
c         write(*,*) 'CHANGING value of mstp(81) from ',mstp(81),' to 20'
c         mstp(81)=20
c      endif
c
c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     photon radiation off quarks and leptons
c       mstj(41)=12              
c     No photon radiation off quarks and leptons
c      mstj(41)=11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      mstp(61)=0                !No IS shower
c      mstp(71)=0                !No FS shower
c      mstp(91)=0                !No Primordial kt
c      mstp(131)=0               !No Pile Up
c      mstp(111)=0               !No hadronization
c      mstp(81)=20               !No Multiple interactions. Force a call to PYEVNW 

c      print *,"Printouts: "
c      print *,"IS shower:  ",mstp(61)
c      print *,"FS shower:  ",mstp(71)
c      print *,"Primordial: ",mstp(91)
c      print *,"Pile-up:    ",mstp(131)
c      print *,"Had.:       ",mstp(111)
c      print *,"MPI:        ",mstp(81)
c      print *,"Photons:    ",mstj(41)

c      mstp(64) =3               !use Lambda_MC for IS shower > 6.4.19
c      mstp(64) =1 !use Lambda_MSbar (default)

c     For the case of jet production the following parameter setting
c     limits the transverse momentum of secondary scatterings, due
c     to multiple parton interactions, to be less than that of the
c     primary interaction (see POWHEG Dijet paper arXiv:1012.3380
c     [hep-ph] sec. 4.1 and also the PYTHIA Manual).
      mstp(86)=1

c     number of warnings printed on the shell
      mstu(26)=20
c     call PYLIST(12) to see the PYTHIA decay table
ccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     tolerate 10% of killed events
      mstu(22)=nev/10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(powheginput("#nohad").gt.0) then
         write(*,*)' ****** switching hadronization off ****'
c     Hadronization off
         mstp(111)=0
c     primordial kt off
         mstp(91)=0
c     No multiple parton interactions
         if(mstp(81).eq.1) then
c     Q2 ordered shower
            mstp(81)=0
         elseif(mstp(81).eq.21) then
c     p_t^2 ordered shower
            mstp(81)=20
         endif
      endif
      end

      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end

ccccccc      subroutine UPINIT
      subroutine UPINIT
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      double precision parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp      
      integer hdecaymode,i      
      real * 8 powheginput
      external powheginput
      nevhep=0
c read the header first, so lprup is set
      call lhefreadhdr(97)

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0

      end

ccccccc      subroutine UPEVNT
      subroutine UPEVNT
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      include 'pwhg_physpar.h'
      real * 8 powheginput
      logical ini,changescalup
      data ini/.true./
      save ini,changescalup
      real * 8 weight
      common/weight2cc/weight
      integer oldnevhep
      data oldnevhep/-1/
      save oldnevhep
      if(ini) then
         if(powheginput("#changescalup").eq.1) then
            print *,"SCALUP is changed..."
            changescalup=.true.
         else
            changescalup=.false.
         endif
         ini=.false.
      endif
      if(nevhep.eq.oldnevhep) then
         call pwhgaccumup
c         call lhefwritev(6)
      endif
      oldnevhep = nevhep
      call lhefreadev(97)
      weight=xwgtup
      if(changescalup) call py_change_scalup
      end
 

      subroutine upveto
c pythia routine to abort event
      end

      subroutine pyabeg
      call init_hist
      end

      subroutine pyaend
      implicit none
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
c      if(rnd_cwhichseed.eq.'none') then
c         open(unit=99,file=pwgprefix(1:lprefix)//
c     1     'POWHEG+PYTHIA-output.top',status='unknown')
c      else
c         open(unit=99,file=pwgprefix(1:lprefix)//'-'//
c     1        rnd_cwhichseed //'-'//
c     2     'POWHEG+PYTHIA-output.top',status='unknown')
c      endif
      call pwhgsetout
      if(rnd_cwhichseed.eq.'none') then
         call pwhgtopout(pwgprefix(1:lprefix)//
     >                   'POWHEG+PYTHIA-output')
      else
         call pwhgtopout(pwgprefix(1:lprefix)//
     >                   '-'//rnd_cwhichseed //'-'//
     >                   'POWHEG+PYTHIA-output')
      endif
c      close(99)
      end


      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      integer hdecaymode
      integer mint
      double precision vint
      COMMON/PYINT1/MINT(400),VINT(400)
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      real * 8 powheginput
      external powheginput
      
      if(mint(51).ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call pylist(7)      !hepeup
            call pylist(2)      !all the event
         endif
         return
      endif

      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end


