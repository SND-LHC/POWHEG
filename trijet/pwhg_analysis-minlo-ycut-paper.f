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
      integer jy,ipt
      real * 8 dy,dpt,dr
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      integer maxjet
      parameter (maxjet=4)
      integer nptmin2
      parameter (nptmin2=5)
      integer nymax
      parameter (nymax = 2)
      character * 4 cptmin2(nptmin2)
      real * 8 ptminarr2(nptmin2)
      data cptmin2/  '-000','-020','-040','-100','-150'/
      data ptminarr2/ 0d0, 20d0, 40d0,  100d0,  150d0/
      integer nptmin3
      parameter (nptmin3=3)
      character * 4 cptmin3(nptmin3)
      real * 8 ptminarr3(nptmin3)
      data cptmin3/  '-000','-020','-040'/
      data ptminarr3/ 0d0, 20d0, 40d0/
      character*4 cymax(nymax)
      data cymax/ '3.0','inf' /
      real*8 ymaxarr(nymax)
      data ymaxarr/3d0,1d10/
      common/infohist/ptminarr2,ptminarr3,cptmin2,cptmin3,
     1     ymaxarr,cymax,cnum
      save /infohist/
      real * 8 powheginput
      external powheginput
      call inihists

      dy=0.2d0
      dpt=10d0
      dr=0.2d0

      do jy=1,nymax

         call bookupeqbins('j3-pt'//'-ymax-'//cymax(jy),dpt,0d0,400d0)
         call bookupeqbins('j3-ptzoom'//'-ymax-'//cymax(jy),
     1        0.5d0,0d0,20d0)

         do ipt=1,nptmin3
            call bookupeqbins('j3-y'//'-pt3min'//cptmin3(ipt)//
     1           '-ymax-'//cymax(jy),dy,-5d0,5d0)
         enddo

         do ipt=1,nptmin2
            call bookupeqbins('j3-pt'//'-pt2min'//cptmin2(ipt)//
     1           '-ymax-'//cymax(jy),dpt,0d0,400d0)
            call bookupeqbins('j3-ptzoom'//'-pt2min'//cptmin2(ipt)//
     1           '-ymax-'//cymax(jy),0.5d0,0d0,20d0)
            call bookupeqbins('j3-inc-pt'//'-pt2min'//cptmin2(ipt)//
     1           '-ymax-'//cymax(jy),dpt,0d0,400d0)
         enddo

c here ymax is intended inclusivly, i.e. only upon the observed jet
         call bookupeqbins('j-pt'//'-yjmax-'//cymax(jy),
     1        dpt,0d0,400d0)
         call bookupeqbins('j-ptzoom'//'-yjmax-'//cymax(jy),
     1        0.5d0,0d0,20d0)
         
      enddo

      do ipt=1,nptmin2
         call bookupeqbins('j-y'//'-ptjmin'//cptmin2(ipt),
     1        dy,-5d0,5d0)
      enddo
      
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
      include 'pwhg_bookhist-multi.h'
c      include 'pwhg_flg.h'
c      include 'LesHouches.h'
      integer   maxjet,mjets,numjets,ntracks
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),yj(maxjet),maxyj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4),ptot(4)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      character * 1 cnum(9)
      integer nptmin2
      parameter (nptmin2=5)
      integer nymax
      parameter (nymax = 2)
      character * 4 cptmin2(nptmin2)
      real * 8 ptminarr2(nptmin2)      
      integer nptmin3
      parameter (nptmin3=3)
      character * 4 cptmin3(nptmin3)
      real * 8 ptminarr3(nptmin3)
      character*4 cymax(nymax)
      real*8 ymaxarr(nymax)
      common/infohist/ptminarr2,ptminarr3,cptmin2,cptmin3,
     1     ymaxarr,cymax,cnum
      save /infohist/
      integer j,jy,ipt,njets
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp
      real * 8 palg,ptminfastjet,rjet
      real * 8 dsig(maxmulti)
      integer nweights
      logical ini
      data ini/.true./
      logical inimulti
      data inimulti/.true./
      integer  minlo
      data minlo/0/
      save inimulti,minlo,ini
      integer minnumjets,maxnumjets
      logical pwhg_isfinite
      external pwhg_isfinite
      character * 4 prefix
c
      if (.not.pwhg_isfinite(dsig0)) then
         write(*,*) '*** PROBLEMS in subroutine analysis ***'
         return
      endif

c      write(*,*) 'nhep ',nhep

      maxnumjets=4

c      call reweightifneeded(dsig0,dsig)

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

c This routine should be used with minlo for trijet:
      minnumjets = 2

      ntracks=0
      mjets=0
c     Loop over final state particles to find jets 
      do ihep=1,nhep
         if (isthep(ihep).eq.1) then
            if (ntracks.eq.maxtrack) then
               write(*,*) 'Too many particles. Increase maxtrack.'//
     #          ' PROGRAM ABORTS'
               call pwhg_exit(-1)
            endif
c     copy momenta to construct jets 
            ntracks=ntracks+1
            ptrack(1:4,ntracks) = phep(1:4,ihep)
         endif
      enddo

      if (ntracks.eq.0) then
         numjets=0
      else
c     palg=1 is standard kt, -1 is antikt
         palg = -1d0
         rjet = 0.5d0              ! radius parameter
         ptminfastjet = 0
         call fastjetppgenkt(ptrack,ntracks,rjet,palg,ptminfastjet,
     $        pj,numjets,jetvec)
      endif

c      write(*,*) 'numjets ',numjets

      ktj(1:maxnumjets) = 0d0
      yj(1:maxnumjets)  = 0d0
      njets = min(maxnumjets,numjets)
      do j=1,njets
         ktj(j) = sqrt(pj(1,j)**2 + pj(2,j)**2 )            
         call pwhg_getrapidity(pj(:,j),yj(j))
      enddo

      maxyj(1:maxnumjets) = 0d0
      do j=1,njets
         if(j.eq.1) then
            maxyj(j) = abs(yj(j))
         else
            maxyj(j) = max(abs(yj(j)),maxyj(j-1))
         endif
      enddo

      do jy=1,nymax
         if(njets.ge.3.and.maxyj(3).lt.ymaxarr(jy)) then
c fill 3 jet observables
            call filld('j3-pt'//'-ymax-'//cymax(jy),ktj(3),dsig)
            call filld('j3-ptzoom'//'-ymax-'//cymax(jy),ktj(3),dsig)
            do ipt=1,nptmin3
               if(ktj(3).gt.ptminarr3(ipt)) then
                  call filld('j3-y'//'-pt3min'//cptmin3(ipt)//
     1                 '-ymax-'//cymax(jy),yj(3),dsig)
               endif
            enddo

            do ipt=1,nptmin2
               if(ktj(2).gt.ptminarr2(ipt)) then
                  call filld('j3-pt'//'-pt2min'//cptmin2(ipt)//
     1                 '-ymax-'//cymax(jy),ktj(3),dsig)
                  call filld('j3-ptzoom'//'-pt2min'//cptmin2(ipt)//
     1                 '-ymax-'//cymax(jy),ktj(3),dsig)
               endif
            enddo
            
         endif

c inclusive jet plots
         do j=1,njets
            if(abs(yj(j)).lt.ymaxarr(jy)) then
               call filld('j-pt'//'-yjmax-'//cymax(jy),
     1              ktj(j),dsig)
               call filld('j-ptzoom'//'-yjmax-'//cymax(jy),
     1              ktj(j),dsig)
               if(j.ge.3) then
                  do ipt=1,nptmin2
                     if(ktj(2).gt.ptminarr2(ipt)) then
                        call filld('j3-inc-pt'//'-pt2min'//cptmin2(ipt)
     1                       //'-ymax-'//cymax(jy),
     1                       ktj(j),dsig)
                     endif
                  enddo
               endif
            endif
         enddo
      enddo

      do ipt=1,nptmin2
         do j=1,njets
            if(ktj(j).gt.ptminarr2(ipt)) then
               call filld('j-y'//'-ptjmin'//cptmin2(ipt),
     1              yj(j),dsig)
            endif
         enddo
      enddo

      end





c      subroutine yetaptmassplot(p,dsig,prefix)
c      implicit none
c      real * 8 p(4),dsig
c      character *(*) prefix
c      real * 8 y,eta,pt,m
c      call getyetaptmass(p,y,eta,pt,m)
c      call filld(prefix//'-y',y,dsig)
c      call filld(prefix//'-eta',eta,dsig)
c      call filld(prefix//'-pt',pt,dsig)
c      call filld(prefix//'-m',m,dsig)
c      end

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
