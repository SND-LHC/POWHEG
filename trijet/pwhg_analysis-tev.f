c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer i,j,k
      real * 8 powheginput
      external powheginput
c
      integer nyrangesIncJet
      parameter (nyrangesIncJet = 6)
      real*8 yrangesIncJet(nyrangesIncJet + 1)
      data yrangesIncJet/0.0d0,0.4d0,0.8d0,1.2d0,1.6d0,2.0d0,2.4d0/
      character*3 cyrangesIncJet(nyrangesIncJet + 1)
      data cyrangesIncJet/'0.0','0.4','0.8','1.2','1.6','2.0','2.4'/
c 0.0 < |y| < 0.4
      integer nptbins1IncJet
      parameter (nptbins1IncJet = 23)
      real*8 ptbins1IncJet(nptbins1IncJet + 1)
      data ptbins1IncJet/50d0,60d0,70d0,80d0,90d0,100d0,110d0,120d0,
     >                   130d0,145d0,160d0,180d0,200d0,220d0,240d0,
     >                   265d0,295d0,325d0,360d0,400d0,445d0,490d0,
     >                   540d0,665d0/
c 0.4 < |y| < 0.8
      integer nptbins2IncJet
      parameter (nptbins2IncJet = 22)
      real*8 ptbins2IncJet(nptbins2IncJet + 1)
      data ptbins2IncJet/50d0,60d0,70d0,80d0,90d0,100d0,110d0,120d0,
     >                   130d0,145d0,160d0,180d0,200d0,220d0,240d0,
     >                   265d0,295d0,325d0,360d0,400d0,445d0,495d0,
     >                   635d0/
c 0.8 < |y| < 1.2
      integer nptbins3IncJet
      parameter (nptbins3IncJet = 20)
      real*8 ptbins3IncJet(nptbins3IncJet + 1)
      data ptbins3IncJet/50d0,60d0,70d0,80d0,90d0,100d0,110d0,125d0,
     >                   140d0,155d0,170d0,190d0,210d0,230d0,250d0,
     >                   270d0,300d0,335d0,375d0,415d0,520d0/
c 1.2 < |y| < 1.6
      integer nptbins4IncJet
      parameter (nptbins4IncJet = 17)
      real*8 ptbins4IncJet(nptbins4IncJet + 1)
      data ptbins4IncJet/50d0,60d0,70d0,80d0,90d0,100d0,110d0,125d0,
     >                   140d0,155d0,170d0,190d0,215d0,240d0,265d0,
     >                   290d0,325d0,415d0/
c 1.6 < |y| < 2.0
      integer nptbins5IncJet
      parameter (nptbins5IncJet = 15)
      real*8 ptbins5IncJet(nptbins5IncJet + 1)
      data ptbins5IncJet/50d0,60d0,70d0,80d0,90d0,100d0,110d0,125d0,
     >                   140d0,160d0,175d0,190d0,210d0,235d0,260d0,
     >                   320d0/
c 2.0 < |y| < 2.4
      integer nptbins6IncJet
      parameter (nptbins6IncJet = 13)
      real*8 ptbins6IncJet(nptbins6IncJet + 1)
      data ptbins6IncJet/50d0,60d0,70d0,80d0,90d0,100d0,110d0,120d0,
     >                   130d0,145d0,160d0,175d0,200d0,230d0/
c ****xxxx****
      integer nymaxrangesMjj
      parameter (nymaxrangesMjj = 6)
      real*8 ymaxrangesMjj(nymaxrangesMjj + 1)
      data ymaxrangesMjj/0.0d0,0.4d0,0.8d0,1.2d0,1.6d0,2.0d0,2.4d0/
      character*3 cymaxrangesMjj(nymaxrangesMjj + 1)
      data cymaxrangesMjj/'0.0','0.4','0.8','1.2','1.6','2.0','2.4'/
c 0.0 < |ymax| < 0.4
      integer nptbins1Mjj
      parameter (nptbins1Mjj = 17)
      real*8 ptbins1Mjj(nptbins1Mjj + 1)
      data ptbins1Mjj/150d0,175d0,200d0,225d0,250d0,300d0,350d0,
     >                400d0,450d0,500d0,560d0,620d0,690d0,770d0,
     >                860d0,950d0,1050d0,1300d0/
c 0.4 < |ymax| < 0.8
      integer nptbins2Mjj
      parameter (nptbins2Mjj = 17)
      real*8 ptbins2Mjj(nptbins2Mjj + 1)
      data ptbins2Mjj/150d0,175d0,200d0,225d0,250d0,300d0,350d0,
     >                400d0,450d0,500d0,560d0,620d0,690d0,770d0,
     >                860d0,950d0,1050d0,1300d0/
c 0.8 < |ymax| < 1.2
      integer nptbins3Mjj
      parameter (nptbins3Mjj = 11)
      real*8 ptbins3Mjj(nptbins3Mjj + 1)
      data ptbins3Mjj/250d0,300d0,350d0,400d0,450d0,500d0,600d0,
     >                700d0,830d0,960d0,1080d0,1400d0/
c 1.2 < |ymax| < 1.6
      integer nptbins4Mjj
      parameter (nptbins4Mjj = 10)
      real*8 ptbins4Mjj(nptbins4Mjj + 1)
      data ptbins4Mjj/300d0,350d0,400d0,450d0,500d0,600d0,700d0,
     >                800d0,960d0,1080d0,1400d0/
c 1.6 < |ymax| < 2.0
      integer nptbins5Mjj
      parameter (nptbins5Mjj = 8)
      real*8 ptbins5Mjj(nptbins5Mjj + 1)
      data ptbins5Mjj/450d0,500d0,600d0,700d0,800d0,920d0,1040d0,
     >                1160d0,1500d0/
c 2.0 < |ymax| < 2.4
      integer nptbins6Mjj
      parameter (nptbins6Mjj = 8)
      real*8 ptbins6Mjj(nptbins6Mjj + 1)
      data ptbins6Mjj/450d0,500d0,600d0,700d0,800d0,920d0,1040d0,
     >                1160d0,1500d0/
c ****xxxx****
      integer nyrangesMjjj
      parameter (nyrangesMjjj = 3)
      real*8 yrangesMjjj(nyrangesMjjj)
      data yrangesMjjj/0.8d0,1.6d0,2.4d0/
      character*3 cyrangesMjjj(nyrangesMjjj)
      data cyrangesMjjj/'0.8','1.6','2.4'/
      integer npt3rangesMjjj
      parameter (npt3rangesMjjj = 3)
      real*8 pt3rangesMjjj(npt3rangesMjjj)
      data pt3rangesMjjj/40d0,70d0,100d0/
      character*3 cpt3rangesMjjj(npt3rangesMjjj)
      data cpt3rangesMjjj/' 40',' 70','100'/
c
c |y| < 0.8, pt3 > 40 GeV
      integer nmjjjbins1Mjjj
      parameter (nmjjjbins1Mjjj = 9)
      real*8 mjjjbins1Mjjj(nmjjjbins1Mjjj + 1)
      data mjjjbins1Mjjj/400d0,450d0,500d0,550d0,610d0,670d0,740d0,
     >                  810d0,900d0,1100d0/
c |y| < 1.6, pt3 > 40 GeV
      integer nmjjjbins2Mjjj
      parameter (nmjjjbins2Mjjj = 10)
      real*8 mjjjbins2Mjjj(nmjjjbins2Mjjj + 1)
      data mjjjbins2Mjjj/400d0,450d0,500d0,550d0,610d0,670d0,740d0,
     >                   810d0,900d0,1100d0,1500d0/
c |y| < 2.4, pt3 > 40 GeV
      integer nmjjjbins3Mjjj
      parameter (nmjjjbins3Mjjj = 10)
      real*8 mjjjbins3Mjjj(nmjjjbins3Mjjj + 1)
      data mjjjbins3Mjjj/400d0,450d0,500d0,550d0,610d0,670d0,740d0,
     >                   810d0,900d0,1100d0,1500d0/
c |y| < 2.4, pt3 > 70 GeV
      integer nmjjjbins4Mjjj
      parameter (nmjjjbins4Mjjj = 10)
      real*8 mjjjbins4Mjjj(nmjjjbins4Mjjj + 1)
      data mjjjbins4Mjjj/400d0,450d0,500d0,550d0,610d0,670d0,740d0,
     >                   810d0,900d0,1100d0,1500d0/
c |y| < 2.4, pt3 > 100 GeV
      integer nmjjjbins5Mjjj
      parameter (nmjjjbins5Mjjj = 10)
      real*8 mjjjbins5Mjjj(nmjjjbins5Mjjj + 1)
      data mjjjbins5Mjjj/400d0,450d0,500d0,550d0,610d0,670d0,740d0,
     >                   810d0,900d0,1100d0,1500d0/
c
      call inihists
c

c vvvvvvvvvvvv This part concerns the inclusive jet pt vvvvvvvvvvvvvvv
c **** arXiv:0802.2400 ****
c
c Inclusive jet pt, 0.0 < |y| < 0.4:
      call bookup('inc. jetpt, '//cyrangesIncJet(1)//' < |y| < '//
     >            cyrangesIncJet(2)//' D0',
     >            nptbins1IncJet,ptbins1IncJet)
c Inclusive jet pt, 0.4 < |y| < 0.8:
      call bookup('inc. jetpt, '//cyrangesIncJet(2)//' < |y| < '//
     >            cyrangesIncJet(3)//' D0',
     >            nptbins2IncJet,ptbins2IncJet)
c Inclusive jet pt, 0.8 < |y| < 1.2:
      call bookup('inc. jetpt, '//cyrangesIncJet(3)//' < |y| < '//
     >            cyrangesIncJet(4)//' D0',
     >            nptbins3IncJet,ptbins3IncJet)
c Inclusive jet pt, 1.2 < |y| < 1.6:
      call bookup('inc. jetpt, '//cyrangesIncJet(4)//' < |y| < '//
     >            cyrangesIncJet(5)//' D0',
     >            nptbins4IncJet,ptbins4IncJet)
c Inclusive jet pt, 1.6 < |y| < 2.0:
      call bookup('inc. jetpt, '//cyrangesIncJet(5)//' < |y| < '//
     >            cyrangesIncJet(6)//' D0',
     >            nptbins5IncJet,ptbins5IncJet)
c Inclusive jet pt, 2.0 < |y| < 2.4:
      call bookup('inc. jetpt, '//cyrangesIncJet(6)//' < |y| < '//
     >            cyrangesIncJet(7)//' D0',
     >            nptbins6IncJet,ptbins6IncJet)
c
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c vvvvvvvvvvvvvvvv This part concerns the dijet mass vvvvvvvvvvvvvvvvv
c **** arXiv:1002.4594 ****
c
c Dijet mass, 0.0 < |ymax| < 0.4
      call bookup('mjj, '//cymaxrangesMjj(1)//' < |ymax| < '//
     >            cymaxrangesMjj(2)//' D0',nptbins1Mjj,ptbins1Mjj)
c Dijet mass, 0.4 < |ymax| < 0.8
      call bookup('mjj, '//cymaxrangesMjj(2)//' < |ymax| < '//
     >            cymaxrangesMjj(3)//' D0',nptbins2Mjj,ptbins2Mjj)
c Dijet mass, 0.8 < |ymax| < 1.2
      call bookup('mjj, '//cymaxrangesMjj(3)//' < |ymax| < '//
     >            cymaxrangesMjj(4)//' D0',nptbins3Mjj,ptbins3Mjj)
c Dijet mass, 1.2 < |ymax| < 1.6
      call bookup('mjj, '//cymaxrangesMjj(4)//' < |ymax| < '//
     >            cymaxrangesMjj(5)//' D0',nptbins4Mjj,ptbins4Mjj)
c Dijet mass, 1.6 < |ymax| < 2.0
      call bookup('mjj, '//cymaxrangesMjj(5)//' < |ymax| < '//
     >            cymaxrangesMjj(6)//' D0',nptbins5Mjj,ptbins5Mjj)
c Dijet mass, 2.0 < |ymax| < 2.4
      call bookup('mjj, '//cymaxrangesMjj(6)//' < |ymax| < '//
     >            cymaxrangesMjj(7)//' D0',nptbins6Mjj,ptbins6Mjj)
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c vvvvvvvvvvvvvvvv This part concerns the trijet mass vvvvvvvvvvvvvvvvv
c **** arXiv:1104.1986 ****
c Trijet mass, |y| < 0.8, 40 GeV < ptj3
      call bookup('mjjj, |y| < '//cyrangesMjjj(1)//', '//
     >            cpt3rangesMjjj(1)//' < ptj3, D0',
     >            nmjjjbins1Mjjj,mjjjbins1Mjjj)
c Trijet mass, |y| < 1.6, 40 GeV < ptj3
      call bookup('mjjj, |y| < '//cyrangesMjjj(2)//', '//
     >            cpt3rangesMjjj(1)//' < ptj3, D0',
     >            nmjjjbins2Mjjj,mjjjbins2Mjjj)
c Trijet mass, |y| < 2.4, 40 GeV < ptj3
      call bookup('mjjj, |y| < '//cyrangesMjjj(3)//', '//
     >            cpt3rangesMjjj(1)//' < ptj3, D0',
     >            nmjjjbins3Mjjj,mjjjbins3Mjjj)
c Trijet mass, |y| < 2.4, 70 GeV < ptj3
      call bookup('mjjj, |y| < '//cyrangesMjjj(3)//', '//
     >            cpt3rangesMjjj(2)//' < ptj3, D0',
     >            nmjjjbins4Mjjj,mjjjbins4Mjjj)
c Trijet mass, |y| < 2.4, 100 GeV < ptj3
      call bookup('mjjj, |y| < '//cyrangesMjjj(3)//', '//
     >            cpt3rangesMjjj(3)//' < ptj3, D0',
     >            nmjjjbins5Mjjj,mjjjbins5Mjjj)
c

      end
     

      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h' 
      include 'pwhg_rad.h' 
      include 'pwhg_weights.h'
      integer   maxjet,mjets,njets,numjets,ntracks
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4),ptot(4),
     2    jetptrel(maxjet)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer ih,il,inu
      real * 8 psum(4)
      real * 8 httot,y,eta,pt,m
      real * 8 dy,deta,delphi,dr
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp
      integer mu,ilep,ivl,nlep,nvl
      real * 8 dsig(7)
      integer nweights
      logical ini
      data ini/.true./
      logical inimulti
      data inimulti/.true./
      integer  minlo
c      data minlo/0/
      save inimulti,minlo,ini
c
      logical pwhg_isfinite
      external pwhg_isfinite
c
      integer iy,ijet,ipt
c
      real*8 p_tmp(4)
      real*8 pj1(4),pj2(4),pj3(4)
      real*8 yj1,yj2,yj3
      real*8 ymax,dymax
      real*8 ptj1,ptj2,ptj3
      real*8 dR12,dR13,dR23
c
      real*8 R
      real*8 MinPtJet,MaxRapJet
c
      real*8 ptmin 
      parameter (ptmin = 8d0)
      real*8 overlapfact
      parameter (overlapfact = 0.5d0)
c
c Declarations for the inclusive jet:
      integer nyrangesIncJet
      parameter (nyrangesIncJet = 6)
      real*8 yrangesIncJet(nyrangesIncJet + 1)
      data yrangesIncJet/0.0d0,0.4d0,0.8d0,1.2d0,1.6d0,2.0d0,2.4d0/
      character*3 cyrangesIncJet(nyrangesIncJet + 1)
      data cyrangesIncJet/'0.0','0.4','0.8','1.2','1.6','2.0','2.4'/
c
      real*8 RconeIncJet
      parameter (RconeIncJet = 0.7d0)
      real*8 MinPtJetIncJet
      parameter (MinPtJetIncJet = 50d0)
      real*8 MaxRapJetIncJet
      parameter (MaxRapJetIncJet = 2.4d0)
c
c Declarations for the dijet mass:
      integer nymaxrangesMjj
      parameter (nymaxrangesMjj = 6)
      real*8 ymaxrangesMjj(nymaxrangesMjj + 1)
      data ymaxrangesMjj/0.0d0,0.4d0,0.8d0,1.2d0,1.6d0,2.0d0,2.4d0/
      character*3 cymaxrangesMjj(nymaxrangesMjj + 1)
      data cymaxrangesMjj/'0.0','0.4','0.8','1.2','1.6','2.0','2.4'/
c
      real*8 RconeMjj
      parameter (RconeMjj = 0.7d0)
      real*8 MinPtJetMjj
      parameter (MinPtJetMjj = 40d0)
      real*8 MaxRapJetMjj
      parameter (MaxRapJetMjj = 1.0d3)
c
c Declarations for the trijet mass:
      integer nyrangesMjjj
      parameter (nyrangesMjjj = 3)
      real*8 yrangesMjjj(nyrangesMjjj)
      data yrangesMjjj/0.8d0,1.6d0,2.4d0/
      character*3 cyrangesMjjj(nyrangesMjjj)
      data cyrangesMjjj/'0.8','1.6','2.4'/
      integer npt3rangesMjjj
      parameter (npt3rangesMjjj = 3)
      real*8 pt3rangesMjjj(npt3rangesMjjj)
      data pt3rangesMjjj/40d0,70d0,100d0/
      character*3 cpt3rangesMjjj(npt3rangesMjjj)
      data cpt3rangesMjjj/' 40',' 70','100'/
c
      real*8 RconeMjjj
      parameter (RconeMjjj = 0.7d0)
      real*8 MinPt1stJetMjjj
      parameter (MinPt1stJetMjjj = 150d0)
      real*8 MinPt3rdJetMjjj
      parameter (MinPt3rdJetMjjj = 40d0)
      real*8 MaxRapJetMjjj
      parameter (MaxRapJetMjjj = 1.0d3)
      real*8 MinRsepMjjj
      parameter (MinRsepMjjj = 1.4d0)
c
c     
      if (ini) then
         minlo=1
         if (powheginput("#minlo").eq.0d0) minlo = 0
         ini = .false.
      end if
c     
      if (.not.pwhg_isfinite(dsig0)) then
         write(*,*) '*** PROBLEMS in subroutine analysis ***'
         return
      endif

      if(inimulti) then
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         inimulti=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
         nweights=1
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
          nweights=weights_num
      endif

      if(sum(abs(dsig)).eq.0) return

      ntracks=0
      mjets=0
c     Loop over final state particles to find jets 
      do ihep=1,nhep
         if (isthep(ihep).eq.1) then
           if (ntracks.eq.maxtrack) then
              write(*,*) 'Too many particles. Increase maxtrack.'//
     #             ' PROGRAM ABORTS'
              call pwhg_exit(-1)
           endif
c     copy momenta to construct jets 
           ntracks=ntracks+1
           do mu=1,4
              ptrack(mu,ntracks)=phep(mu,ihep)
           enddo
        endif
      enddo

c vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c Inclusive jet cross section:
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
        R = RconeIncJet
        MinPtJet = MinPtJetIncJet
        MaxRapJet = MaxRapJetIncJet
        call fastjetd0runiicone(ptrack,ntracks,R,ptmin,overlapfact,
     >                          MinPtJet,MaxRapJet,pj,njets,jetptrel)
c
c We need at least one jet:
        if (njets.lt.1) goto 111
c
c We try to fill the histograms:
        do iy=1,nyrangesIncJet
          do ijet=1,njets
            p_tmp = pj(:,ijet)
            call getyetaptmass(p_tmp,y,eta,pt,m)
            if ((abs(y).gt.yrangesIncJet(iy)).and.
     >          (abs(y).lt.yrangesIncJet(iy+1))) then
              dy = 2d0*(yrangesIncJet(iy+1) - yrangesIncJet(iy))
              call filld('inc. jetpt, '//
     >                   cyrangesIncJet(iy)//' < |y| < '//
     >                   cyrangesIncJet(iy+1)//' D0',pt,dsig/dy)
c              exit
            end if
          end do
        end do
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

111   continue

c vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c Dijet mass:
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
        R = RconeMjj
        MinPtJet = MinPtJetMjj
        MaxRapJet = MaxRapJetMjj
        call fastjetd0runiicone(ptrack,ntracks,R,ptmin,overlapfact,
     >                          MinPtJet,MaxRapJet,pj,njets,jetptrel)
c
c We need at least two jets:
        if (njets.lt.2) goto 222
c
c We select the two hardest jets in the event:
        pj1 = pj(:,1)
        pj2 = pj(:,2)
c
        call pwhg_getrapidity(pj1,yj1)
        call pwhg_getrapidity(pj2,yj2)
        call pwhg_getinvmass(pj1 + pj2,m)
c
        ymax = max(abs(yj1),abs(yj2))
c
c We try to fill the histograms:
        do iy=1,nymaxrangesMjj
          if ((ymax.gt.ymaxrangesMjj(iy)).and.
     >        (ymax.lt.ymaxrangesMjj(iy+1))) then
            dymax = ymaxrangesMjj(iy+1) - ymaxrangesMjj(iy)
            call filld('mjj, '//cymaxrangesMjj(iy)//' < |ymax| < '//
     >            cymaxrangesMjj(iy+1)//' D0',m,dsig/dymax)
          end if
        end do
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

222   continue

c vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c Trijet mass:
      R = RconeMjjj
      MinPtJet = MinPt3rdJetMjjj
      MaxRapJet = MaxRapJetMjjj
      call fastjetd0runiicone(ptrack,ntracks,R,ptmin,overlapfact,
     >                        MinPtJet,MaxRapJet,pj,njets,jetptrel)
c
c We need at least two jets:
      if (njets.lt.3) goto 333
c
c The three hardest jets are selected:
      pj1 = pj(:,1)
      pj2 = pj(:,2)
      pj3 = pj(:,3)
c
      call pwhg_getrapidity(pj1,yj1)
      call pwhg_getrapidity(pj2,yj2)
      call pwhg_getrapidity(pj3,yj3)
c
      call pwhg_getR_phiy(pj1,pj2,dR12)
      call pwhg_getR_phiy(pj1,pj3,dR13)
      call pwhg_getR_phiy(pj2,pj3,dR23)
c
      call pwhg_getinvmass(pj1 + pj2 + pj3,m)
c
      call pwhg_getpt(pj1,ptj1)
      call pwhg_getpt(pj2,ptj2)
      call pwhg_getpt(pj3,ptj3)
c
c The jets should be separated:
      if (min(dR12,dR13,dR23).lt.MinRsepMjjj) goto 333
c The first jet should have a larger minimal pt:
      if (ptj1.lt.MinPt1stJetMjjj) goto 333
c
      do iy=1,nyrangesMjjj
        if (max(abs(yj1),abs(yj2),abs(yj3)).lt.yrangesMjjj(iy)) then
          call filld('mjjj, |y| < '//cyrangesMjjj(iy)//', '//
     >               cpt3rangesMjjj(1)//' < ptj3, D0',
     >             m,dsig)
          if (iy.eq.nyrangesMjjj) then
c Filling up histos with different pt ranges:
            do ipt=2,npt3rangesMjjj
              if (ptj3.gt.pt3rangesMjjj(ipt)) then
                call filld('mjjj, |y| < '//cyrangesMjjj(iy)//', '//
     >                     cpt3rangesMjjj(ipt)//' < ptj3, D0',
     >                     m,dsig)
              end if
            end do
          end if
        end if
      end do
c
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

333   continue

      end

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig(7)
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end

      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass
      call pwhg_getrapidity(p,y)      
      pt=sqrt(p(1)**2+p(2)**2)
      call pwhg_getpseudorapidity(p,eta)
      call pwhg_getinvmass(p,mass)
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end


      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include 'hepevt.h'
      integer j,k
      real * 8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
c bubble sort
      touched=.true.
      do while(touched)
         touched=.false.
         do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
               k=iarr(j)
               iarr(j)=iarr(j+1)
               iarr(j+1)=k
               tmp=pt(j)
               pt(j)=pt(j+1)
               pt(j+1)=tmp
               touched=.true.
            endif
         enddo
      enddo
      end

      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end

      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      subroutine pwhgfinalopshist
      end

      function calcrho(ntrack,njet,ptrack,pjet,jetvec,
     >                 ptmin,ptmax,ymin,ymax,Rpar,dr,r)
      implicit none
c
      real*8 calcrho
c
      integer maxtrack,maxjet
      parameter (maxtrack=2048)
      parameter (maxjet=2048)
c
      integer ntrack,njet
      real*8 ptrack(4,maxtrack),pjet(4,maxjet)
      integer jetvec(maxtrack)
      real*8 ptmin,ptmax
      real*8 ymin,ymax
      real*8 Rpar,dr,r
c
      integer ijet
      integer numjets
      real*8 p_j(4)
      real*8 pt_j,y_j
      real*8 pt_dr,pt_r
c
      real*8 calcptdr
      external calcptdr
c
      calcrho = 0d0
      numjets = 0
c
c we loop over all jets:
      do ijet=1,njet
        p_j = pjet(:,ijet)
        call pwhg_getpt(p_j,pt_j)
        call pwhg_getrapidity(p_j,y_j)
c We only consider those jets for which the pt is in a given range:
        if ((pt_j.lt.ptmin).or.(pt_j.gt.ptmax)) cycle
c We even include the possibility to make cuts on rapidity:
        if ((abs(y_j).lt.ymin).or.(abs(y_j).gt.ymax)) cycle
c We have an extra cut on the jet pt, hence njet is not the number of
c jets considered for rho:
        numjets = numjets + 1
c we calculate the pt contribution of an annulus between r - \Delta r/2
c and r + \Delta r/2
        pt_dr = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,r,dr)
c We calculate the summed pt of all the clusters in the jet:
        pt_r = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,
     >                  Rpar/2d0,Rpar)
C        calcrho = calcrho + pt_dr/pt_j
        calcrho = calcrho + pt_dr/pt_r
      end do
c normalization:
      if (numjets.gt.0) calcrho = calcrho/numjets/dr
c
      end function calcrho
c
c This function calculates the pt contribution coming from an annulus
c between r - dr/2 and r + dr/2. The jet orientation is
c given by the jet momentum (pj). 
      function calcptdr(ijet,pj,ntrack,ptrack,jetvec,r,dr)
      implicit none
c
      real*8 calcptdr
c
      integer maxtrack
      parameter (maxtrack = 2048)
c
      integer ijet
      real*8 pj(4)
      integer ntrack
      real*8 ptrack(4,maxtrack)
      integer jetvec(maxtrack)
      real*8 r,dr
c
      integer itrack
      real*8 dr2
      real*8 deltar,pt_tmp
c
c
      calcptdr = 0d0
c
      dr2 = dr/2d0
c
c We loop over all tracks and select only those which are part
c of the ijet-th jet:
      do itrack=1,ntrack
        if (jetvec(itrack).ne.ijet) cycle
c If we are here the itrack-th track is part of the ijet-th jet:
        call pwhg_getR_phiy(pj,ptrack(:,itrack),deltar)
c is it in the annulus?
        if ((deltar.lt.(r - dr2)).or.(deltar.gt.(r + dr2))) cycle
        call pwhg_getpt(ptrack(:,itrack),pt_tmp)
        calcptdr = calcptdr + pt_tmp
      end do
c
      end function calcptdr
c
c This routine calculates the integrated jet shape. 
c ptmin, ptmax are the cuts on the pt of the considered jets,
c ymin and ymax have the same role but with the rapidities,
c while ptw,yw variables are used to define the pt and y
c window.
      function calcpsi(ntrack,njet,ptrack,pjet,jetvec,
     >                 ptmin,ptmax,ymin,ymax,
     >                 ptwmin,ptwmax,ywmin,ywmax,
     >                 Rpar,r)
      implicit none
c
      real*8 calcpsi
c
      integer maxtrack,maxjet
      parameter (maxtrack=2048)
      parameter (maxjet=2048)
c
      integer ntrack,njet
      real*8 ptrack(4,maxtrack),pjet(4,maxjet)
      integer jetvec(maxtrack)
      real*8 ptmin,ptmax
      real*8 ymin,ymax
      real*8 ptwmin,ptwmax
      real*8 ywmin,ywmax
      real*8 Rpar,r
c
      integer ijet
      integer numjets
      real*8 p_j(4)
      real*8 pt_j,y_j
      real*8 pt_r
c
      real*8 calcptdr
      external calcptdr
c
c We set psi to an unphysical value:
      calcpsi = -1d5
      numjets = 0
c
c we loop over all jets:
      do ijet=1,njet
        p_j = pjet(:,ijet)
        call pwhg_getpt(p_j,pt_j)
        call pwhg_getrapidity(p_j,y_j)
c We only consider those jets for which the pt is in a given range:
        if ((pt_j.lt.ptmin).or.(pt_j.gt.ptmax)) cycle
c We even include the possibility to make cuts on rapidity:
        if ((abs(y_j).lt.ymin).or.(abs(y_j).gt.ymax)) cycle
        numjets = numjets + 1
c If we have at least one jet the psi is initialized to be physical:
        if (numjets.eq.1) calcpsi = 0d0
        if ((pt_j.lt.ptwmin).or.(pt_j.gt.ptwmax)) cycle
        if ((abs(y_j).lt.ywmin).or.(abs(y_j).gt.ywmax)) cycle
c The next few lines are needed for the window:
c we calculate the pt contribution of between 0 and r:
c we can use the previously implemented calcptdr function 
c with parameters: r = r/2 and dr = r, with this
c parameter choice the function calculates the pt 
c contribution from tracks situating in an annulus between 
c zero and r
        pt_r = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,r/2d0,r)
        pt_j = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,Rpar/2d0,Rpar)
        calcpsi = calcpsi + pt_r/pt_j
      end do
c normalization:
      if (numjets.gt.0) calcpsi = calcpsi/numjets
c
      end function calcpsi
