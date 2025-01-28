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
      integer j,k,i,l
      real * 8 dy,dpt,dr
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      integer maxjet
      parameter (maxjet=4)
      integer nptmin
      parameter (nptmin=5)
      integer nymax
      parameter (nymax = 3)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)
      data cptmin/  '-000','-005','-020','-040','-060'/
      data ptminarr/ 0d0, 5d0, 20d0,  40d0,  60d0/
      character*4 cymax(nymax)
      data cymax/ '2.0','3.0','inf' /
      real*8 ymaxarr(nymax)
      data ymaxarr/2d0,3d0,1d10/
      common/infohist/ptminarr,ymaxarr,cnum,cptmin,cymax
      save /infohist/
      real * 8 powheginput
      external powheginput
      call inihists

      dy=0.2d0
      dpt=10d0
      dr=0.2d0
      
      do i=1,nptmin
      do j=1,nymax
        call bookupeqbins('sigtot-ptmin'//cptmin(i)
     >                    //'-ymax-'//cymax(j),1d0,0.5d0,1.5d0)

        do k=1,maxjet
          call bookupeqbins('j'//cnum(k)//'-y-ptmin'//cptmin(i)
     >                      //'-ymax-'//cymax(j),dy,-5d0,5d0)
          call bookupeqbins('j'//cnum(k)//'-eta-ptmin'//cptmin(i)
     >                      //'-ymax-'//cymax(j),dy,-5d0,5d0)
          call bookupeqbins('j'//cnum(k)//'-pt-ptmin'//cptmin(i)
     >                      //'-ymax-'//cymax(j),dpt,0d0,400d0)
          call bookupeqbins('j'//cnum(k)//'-ptzoom-ptmin'//cptmin(i)
     >                      //'-ymax-'//cymax(j),2d0,1d0,151d0)
          call bookupeqbins('j'//cnum(k)//'-m-ptmin'//cptmin(i)
     >                      //'-ymax-'//cymax(j),dpt,0d0,400d0) 
          call bookupeqbins('j'//cnum(k)//'-ptzoom2-ptmin'//cptmin(i)
     >                      //'-ymax-'//cymax(j),0.5d0,0d0,20d0)
        enddo

        do k=1,maxjet-1
          do l=k+1,maxjet
            call bookupeqbins('j'//cnum(k)//'j'//cnum(l)//
     1           '-R'//cptmin(i)//'-ymax-'//cymax(j),dr,0d0,5d0)  
c            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
c     1           '-y'//cptmin(i),dy,-5d0,5d0)  
c            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
c     1           '-eta'//cptmin(i),dy,-5d0,5d0)
c            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
c     1           '-pt'//cptmin(i),dpt,0d0,400d0)
c            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
c     1           '-m'//cptmin(i),dpt,0d0,400d0)  
          enddo
        enddo
      enddo

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
      real * 8  ktj(maxjet),etaj(maxjet),yj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4),ptot(4)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      character * 1 cnum(9)
      integer nptmin
      parameter (nptmin=5)
      integer nymax
      parameter (nymax = 3)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)      
      character*4 cymax(nymax)
      real*8 ymaxarr(nymax)
      common/infohist/ptminarr,ymaxarr,cnum,cptmin,cymax
      save /infohist/
      integer j,k,i,jj,ii,l
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
      real * 8 ptmin,mininvmass
      integer maxnumlep
      parameter (maxnumlep=10)
      real * 8 pvl(4,maxnumlep),plep(4,maxnumlep)
      integer mu,ilep,ivl,nlep,nvl
      logical is_W
      real * 8 mV2,ptvb,mvb,ptlep,ptminfastjet,ptvl,R,ylep,yvb,yvl
      real * 8 Wmass,Wwidth,Wmasslow,Wmasshigh
      integer jpart, jjet
      real * 8 palg
c      real * 8 rescfac1,rescfac2
c      common /crescfac/rescfac1,rescfac2
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
      real*8 maxyj
c
      real*8 MinFirstJetPt
      parameter(MinFirstJetPt = 70d0)
      real * 8 rescfac
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

      if (ntracks.eq.0) then
         numjets=0
      else
c     palg=1 is standard kt, -1 is antikt
         palg = -1d0
         R = 0.5d0              ! radius parameter
         ptminfastjet = ptminarr(1)
         call fastjetppgenkt(ptrack,ntracks,R,palg,ptminfastjet,
     $        pj,numjets,jetvec)
      endif

c      write(*,*) 'numjets ',numjets

      ktj(1:maxnumjets) = 0d0
      yj(1:maxnumjets)  = 0d0
      maxyj = 0d0
      do j=1,min(maxnumjets,numjets)
         ktj(j) = sqrt(pj(1,j)**2 + pj(2,j)**2 )            
         call pwhg_getrapidity(pj(:,j),yj(j))
         if (maxyj.lt.abs(yj(j))) maxyj = abs(yj(j))
      enddo

      do i=1,nptmin        
        do j=1,nymax
c The first jet should have MinFirstJetPt pt, while the second
c ptminarr(i):
          if ((ktj(1).le.MinFirstJetPt).or.
     >        (ktj(2).le.ptminarr(i)).or.
     >        (maxyj.gt.ymaxarr(j))) then
c     since the ptminarr array is ordered, if the second jet is not passing the cut, it will never do
            cycle
          endif

          call filld('sigtot-ptmin'//cptmin(i)
     >               //'-ymax-'//cymax(j),1d0,dsig)
                  
c     jets
          mjets=min(maxnumjets,numjets)


c          if(i.eq.1.and.j.eq.1.and.mjets.ge.3.and.ktj(3).gt.300) then
c             write(111,*) '*** checking minlo rescfac formation *** '
c             call setlocalscales(rad_ubornidx,10,rescfac)
c          endif


         
          do k=1,mjets
            call getyetaptmass(pj(:,k),y,eta,pt,m)
            call filld('j'//cnum(k)//'-y-ptmin'//cptmin(i)
     >                 //'-ymax-'//cymax(j), y, dsig)
            call filld('j'//cnum(k)//'-eta-ptmin'//cptmin(i)
     >                 //'-ymax-'//cymax(j), eta, dsig)
            call filld('j'//cnum(k)//'-pt-ptmin'//cptmin(i)
     >                 //'-ymax-'//cymax(j), pt, dsig)
            call filld('j'//cnum(k)//'-ptzoom-ptmin'//cptmin(i)
     >                 //'-ymax-'//cymax(j), pt, dsig)
            call filld('j'//cnum(k)//'-ptzoom2-ptmin'//cptmin(i)
     >                 //'-ymax-'//cymax(j), pt, dsig)
            call filld('j'//cnum(k)//'-m-ptmin'//cptmin(i)
     >                 //'-ymax-'//cymax(j), m, dsig)
          enddo
         


          do k=1,mjets
            do l=k+1,mjets
c              call getyetaptmass(pj(:,k)+pj(:,l),y,eta,pt,m)
              call pwhg_getR_phiy(pj(:,k),pj(:,l),R)
              call filld('j'//cnum(k)//'j'//cnum(l)//'-R'//cptmin(i)
     >                   //'-ymax-'//cymax(j),R,dsig)
c              call filld('j'//cnum(j)//'j'//cnum(k)//'-y'//cptmin(i),
c     $             y, dsig)
c              call filld('j'//cnum(j)//'j'//cnum(k)//'-eta'//cptmin(i),
c     $             eta, dsig)
c              call filld('j'//cnum(j)//'j'//cnum(k)//'-pt'//cptmin(i),
c     $             pt, dsig)
c              call filld('j'//cnum(j)//'j'//cnum(k)//'-m'//cptmin(i), 
c     $             m, dsig)
            enddo
          enddo
         
        enddo
      enddo

c      read(*,*)

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
