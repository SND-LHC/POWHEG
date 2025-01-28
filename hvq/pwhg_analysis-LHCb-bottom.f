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

      integer,parameter :: nybins=5, nybinsnew=18, nptbins=50, nptbinsLHCb = 27
      integer :: i
      character * 10 :: suffix(nybins)
      real * 8 :: ptbin(nptbins+1),ptbinLHCb(nptbinsLHCb+1)
      common/ptarray/ptbin,ptbinLHCb
      common/suffix/suffix
      character *4 tmpsuff
           
      suffix(1) = '2-2.5' 
      suffix(2) = '2.5-3' 
      suffix(3) = '3-3.5' 
      suffix(4) = '3.5-4' 
      suffix(5) = '4-4.5' 

      ptbin = [0.d0,0.25d0,0.5d0,0.75d0,1.d0,1.25d0,1.5d0,1.75d0,2.d0,2.25d0,2.5d0,2.75d0,3.0d0,3.25d0,3.5d0,3.75d0,4.d0,4.25d0,4.5d0,4.75d0,5.d0,5.5d0,6.d0,6.5d0,7.d0,7.5d0,8.d0,8.5d0,9.d0,9.5d0,10.d0,11.d0,12.d0,13.d0,14.d0,15.d0,16.d0,17.d0,18.d0,19d0,20d0,22d0,24d0,26d0,28d0,30d0,32d0,34d0,36d0,38d0,40d0]

      ptbinLHCb = [0d0, 0.5d0, 1d0, 1.5d0, 2d0, 2.5d0, 3d0, 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0, 6d0, 6.5d0, 7d0, 7.5d0, 8d0, 8.5d0, 9d0, 9.5d0, 10d0, 10.5d0, 11.5d0, 12.5d0, 14.0d0, 16.5d0, 23.5d0, 40d0]
      
      call inihists
            
      do i=1,5
         call bookup('pt-bottom'//trim(adjustl(suffix(i))),nptbins,ptbin)
         call bookup('pt-B0'//trim(adjustl(suffix(i))),nptbins,ptbin)
         call bookup('pt-B0b'//trim(adjustl(suffix(i))),nptbins,ptbin)
         call bookup('pt-B+'//trim(adjustl(suffix(i))),nptbins,ptbin)
         call bookup('pt-B-'//trim(adjustl(suffix(i))),nptbins,ptbin)

         call bookup('pt-bottom'//trim(adjustl(suffix(i)))//'LHCb',nptbinsLHCb,ptbinLHCb)
         call bookup('pt-B0'//trim(adjustl(suffix(i)))//'LHCb',nptbinsLHCb,ptbinLHCb)
         call bookup('pt-B0b'//trim(adjustl(suffix(i)))//'LHCb',nptbinsLHCb,ptbinLHCb)
         call bookup('pt-B+'//trim(adjustl(suffix(i)))//'LHCb',nptbinsLHCb,ptbinLHCb)
         call bookup('pt-B-'//trim(adjustl(suffix(i)))//'LHCb',nptbinsLHCb,ptbinLHCb)
      enddo

      do i=1,nybinsnew
         write(tmpsuff,'(I2.0)') i
         call bookupeqbins('pt-bottom-'//trim(adjustl(tmpsuff)),0.25d0,0d0,25d0)
         call bookupeqbins('pt-B0-'//trim(adjustl(tmpsuff)),0.25d0,0d0,25d0)
         call bookupeqbins('pt-Bpm-'//trim(adjustl(tmpsuff)),0.25d0,0d0,25d0)
         call bookupeqbins('pt-Bs0-'//trim(adjustl(tmpsuff)),0.25d0,0d0,25d0)
c         call bookupeqbins('pt-LambdaB-'//trim(adjustl(tmpsuff)),0.25d0,0d0,25d0)
      enddo 
     
      end



      subroutine analysis(dsig0)
      implicit none
      integer, parameter :: dsigdim = 120
      real * 8  dsig0,dsig(dsigdim)
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      logical ini
      data ini/.true./
      save ini
      integer icut,j,mjets
      real * 8 pjet(4),ppairttbar(4),ptt,yt,etat,pttbar,ytbar,etatbar,
     1     ptttbar,mttbar,yttbar,etattbar,ptjet,yjet,etajet,eep,cthep,
     1     phiep,eem,cthem,phiem,dphi
      integer ihep,it,itbar,iem,iep,ijet,mu
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/
      integer   maxjet
      parameter (maxjet=2048)
      real * 8  kt(maxjet),eta(maxjet),rap(maxjet),
     c    phi(maxjet),pj(4,maxjet)
      integer ibtag(maxjet),itags
      logical sonoftop
      external sonoftop
      character * 4 cut

      integer,parameter :: nybins=5, nybinsnew=18, nptbins=50, nptbinsLHCb = 27
      integer :: i
      character * 10 :: suffix(nybins)
      real * 8 :: ptbin(nptbins+1),ptbinLHCb(nptbinsLHCb+1)
      common/ptarray/ptbin,ptbinLHCb
      common/suffix/suffix

      real * 8 :: ptD,yD,etaD
      character *3 tmpsuff
      
      call multi_plot_setup(dsig0,dsig,dsigdim)

      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         do ihep=3,nhep
            
            if (isthep(ihep) == 1 .and. idhep(ihep) == 5 ) then 
               call ptyeta(phep(1,ihep),ptD,yD,etaD)
               if (yD > 2d0 .and. yD<2.5d0) then
                  call filld('pt-bottom'//trim(adjustl(suffix(1))),ptD,dsig)
                  call filld('pt-bottom'//trim(adjustl(suffix(1)))//'LHCb',ptD,dsig)
               endif
               if (yD > 2.5d0 .and. yD<3d0) then
                  call filld('pt-bottom'//trim(adjustl(suffix(2))),ptD,dsig)
                  call filld('pt-bottom'//trim(adjustl(suffix(2)))//'LHCb',ptD,dsig)
               endif
               if (yD > 3d0 .and. yD<3.5d0) then
                  call filld('pt-bottom'//trim(adjustl(suffix(3))),ptD,dsig)
                  call filld('pt-bottom'//trim(adjustl(suffix(3)))//'LHCb',ptD,dsig)
               endif
               if (yD > 3.5d0 .and. yD<4d0) then
                  call filld('pt-bottom'//trim(adjustl(suffix(4))),ptD,dsig)
                  call filld('pt-bottom'//trim(adjustl(suffix(4)))//'LHCb',ptD,dsig)
               endif
               if (yD > 4d0 .and. yD<4.5d0) then
                     call filld('pt-bottom'//trim(adjustl(suffix(5))),ptD,dsig)
                     call filld('pt-bottom'//trim(adjustl(suffix(5)))//'LHCb',ptD,dsig)
               endif               
            endif
         enddo
         
      else 
!     with showers 

         do ihep=1,nhep

            if(isthep(ihep).eq.1) then

               if(idhep(ihep).eq.511) then
                  
                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  
                  if (yD > 2d0 .and. yD<2.5d0) then
                     call filld('pt-B0'//trim(adjustl(suffix(1))),ptD,dsig)
                     call filld('pt-B0'//trim(adjustl(suffix(1)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 2.5d0 .and. yD<3d0) then
                     call filld('pt-B0'//trim(adjustl(suffix(2))),ptD,dsig)
                     call filld('pt-B0'//trim(adjustl(suffix(2)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3d0 .and. yD<3.5d0) then
                     call filld('pt-B0'//trim(adjustl(suffix(3))),ptD,dsig)
                     call filld('pt-B0'//trim(adjustl(suffix(3)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3.5d0 .and. yD<4d0) then
                     call filld('pt-B0'//trim(adjustl(suffix(4))),ptD,dsig)
                     call filld('pt-B0'//trim(adjustl(suffix(4)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 4d0 .and. yD<4.5d0) then
                     call filld('pt-B0'//trim(adjustl(suffix(5))),ptD,dsig)
                     call filld('pt-B0'//trim(adjustl(suffix(5)))//'LHCb',ptD,dsig)
                  endif
                  
                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-B0-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo                  


               else if (idhep(ihep).eq.-511) then

                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  if (yD > 2d0 .and. yD<2.5d0) then
                     call filld('pt-B0b'//trim(adjustl(suffix(1))),ptD,dsig)
                     call filld('pt-B0b'//trim(adjustl(suffix(1)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 2.5d0 .and. yD<3d0) then
                     call filld('pt-B0b'//trim(adjustl(suffix(2))),ptD,dsig)
                     call filld('pt-B0b'//trim(adjustl(suffix(2)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3d0 .and. yD<3.5d0) then
                     call filld('pt-B0b'//trim(adjustl(suffix(3))),ptD,dsig)
                     call filld('pt-B0b'//trim(adjustl(suffix(3)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3.5d0 .and. yD<4d0) then
                     call filld('pt-B0b'//trim(adjustl(suffix(4))),ptD,dsig)
                     call filld('pt-B0b'//trim(adjustl(suffix(4)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 4d0 .and. yD<4.5d0) then
                     call filld('pt-B0b'//trim(adjustl(suffix(5))),ptD,dsig)
                     call filld('pt-B0b'//trim(adjustl(suffix(5)))//'LHCb',ptD,dsig)
                  endif

                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-B0-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo                  

              else if (idhep(ihep).eq. 521) then

                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  if (yD > 2d0 .and. yD<2.5d0) then
                     call filld('pt-B+'//trim(adjustl(suffix(1))),ptD,dsig)
                     call filld('pt-B+'//trim(adjustl(suffix(1)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 2.5d0 .and. yD<3d0) then
                     call filld('pt-B+'//trim(adjustl(suffix(2))),ptD,dsig)
                     call filld('pt-B+'//trim(adjustl(suffix(2)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3d0 .and. yD<3.5d0) then
                     call filld('pt-B+'//trim(adjustl(suffix(3))),ptD,dsig)
                     call filld('pt-B+'//trim(adjustl(suffix(3)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3.5d0 .and. yD<4d0) then
                     call filld('pt-B+'//trim(adjustl(suffix(4))),ptD,dsig)
                     call filld('pt-B+'//trim(adjustl(suffix(4)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 4d0 .and. yD<4.5d0) then
                     call filld('pt-B+'//trim(adjustl(suffix(5))),ptD,dsig)
                     call filld('pt-B+'//trim(adjustl(suffix(5)))//'LHCb',ptD,dsig)
                  endif

                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-Bpm-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo                  
                  
              else if (idhep(ihep).eq. -521) then

                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  if (yD > 2d0 .and. yD<2.5d0) then
                     call filld('pt-B-'//trim(adjustl(suffix(1))),ptD,dsig)
                     call filld('pt-B-'//trim(adjustl(suffix(1)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 2.5d0 .and. yD<3d0) then
                     call filld('pt-B-'//trim(adjustl(suffix(2))),ptD,dsig)
                     call filld('pt-B-'//trim(adjustl(suffix(2)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3d0 .and. yD<3.5d0) then
                     call filld('pt-B-'//trim(adjustl(suffix(3))),ptD,dsig)
                     call filld('pt-B-'//trim(adjustl(suffix(3)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3.5d0 .and. yD<4d0) then
                     call filld('pt-B-'//trim(adjustl(suffix(4))),ptD,dsig)
                     call filld('pt-B-'//trim(adjustl(suffix(4)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 4d0 .and. yD<4.5d0) then
                     call filld('pt-B-'//trim(adjustl(suffix(5))),ptD,dsig)
                     call filld('pt-B-'//trim(adjustl(suffix(5)))//'LHCb',ptD,dsig)
                  endif

                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-Bpm-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo                  

               else if ( abs(idhep(ihep)) .eq. 531) then 

                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  
                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-Bs0-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo                  
                  
               else if ( abs(idhep(ihep)) .eq. 4122 ) then
                  
                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  
                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-LambdaC-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo
                  
               else if ( abs(idhep(ihep)) .eq. 5) then 
              
                  call ptyeta(phep(1,ihep),ptD,yD,etaD)
                  if (yD > 2d0 .and. yD<2.5d0) then
                     call filld('pt-bottom'//trim(adjustl(suffix(1))),ptD,dsig)
                     call filld('pt-bottom'//trim(adjustl(suffix(1)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 2.5d0 .and. yD<3d0) then
                     call filld('pt-bottom'//trim(adjustl(suffix(2))),ptD,dsig)
                     call filld('pt-bottom'//trim(adjustl(suffix(2)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3d0 .and. yD<3.5d0) then
                     call filld('pt-bottom'//trim(adjustl(suffix(3))),ptD,dsig)
                     call filld('pt-bottom'//trim(adjustl(suffix(3)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 3.5d0 .and. yD<4d0) then
                     call filld('pt-bottom'//trim(adjustl(suffix(4))),ptD,dsig)
                     call filld('pt-bottom'//trim(adjustl(suffix(4)))//'LHCb',ptD,dsig)
                  endif
                  if (yD > 4d0 .and. yD<4.5d0) then
                     call filld('pt-bottom'//trim(adjustl(suffix(5))),ptD,dsig)
                     call filld('pt-bottom'//trim(adjustl(suffix(5)))//'LHCb',ptD,dsig)
                  endif

                  do i=1,nybinsnew
                     if ( yD < i*0.5 .and. yD>(i-1)*0.5d0) then 
                        write(tmpsuff,'(I2.0)') i                     
                        call filld('pt-bottom-'//trim(adjustl(tmpsuff)),ptD,dsig)
                     endif
                  enddo                  
               endif
            endif
         enddo
      endif
      end


      subroutine decvariables(pdec0,ppart0,edec,cthdec,phidec)
      implicit none
      include 'pwhg_math.h'
      real * 8 pdec0(4),ppart0(4),edec,cthdec,phidec
      real * 8 pdec(0:3),ppart(0:3)
      real * 8 vec(3),beta,pt
      integer mu
      do mu=1,3
         pdec(mu)=pdec0(mu)
         ppart(mu)=ppart0(mu)
      enddo
      pdec(0)=pdec0(4)
      ppart(0)=ppart0(4)
      vec(1)=0
      vec(2)=0
      vec(3)=-1
      beta=ppart(3)/ppart(0)
      call mboost(1,vec,beta,pdec,pdec)
      call mboost(1,vec,beta,ppart,ppart)
      pt=sqrt(ppart(1)**2+ppart(2)**2)
      vec(1)=ppart(1)/pt
      vec(2)=ppart(2)/pt
      vec(3)=0
      beta=-pt/ppart(0)
      call mboost(1,vec,beta,pdec,pdec)
      call mboost(1,vec,beta,ppart,ppart)
      edec=pdec(0)
      cthdec=pdec(3)/sqrt(pdec(1)**2+pdec(2)**2+pdec(3)**2)
      phidec=atan2(pdec(2),pdec(1))-atan2(vec(2),vec(1))
c bring it back between -pi and pi
      phidec=phidec-nint(phidec/(2*pi))*2*pi
      end

      subroutine ptyeta(p,pt,y,eta)
      implicit none
      real * 8 p(4),pt,y,eta
      real * 8 pp,tiny
      parameter (tiny=1d-12)
      pt=sqrt(p(1)**2+p(2)**2)
      y=log((p(4)+p(3))/(p(4)-p(3)))/2
      pp=sqrt(pt**2+p(3)**2)*(1+tiny)
      eta=log((pp+p(3))/(pp-p(3)))/2
      end





      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ibtag)
c     arrays to reconstruct jets
      implicit none
      integer mjets
      real * 8  kt(mjets),eta(mjets),rap(mjets),phi(mjets),pjet(4,mjets)
      integer ibtag(mjets)
      include   'hepevt.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=20)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8 r,palg,ptmin,pp,tmp
      logical bson
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,mjets
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
C - Extract final state particles to feed to jet finder
      do j=1,nhep
         if (isthep(j).eq.1) then
            if(ntracks.eq.maxtrack) then
               write(*,*) 'analyze: need to increase maxtrack!'
               write(*,*) 'ntracks: ',ntracks
               stop
            endif
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,j)
            enddo
            itrackhep(ntracks)=j
         endif
      enddo
      if (ntracks.eq.0) then
         return
      endif
C --------------------------------------------------------------------- C
C - Inclusive jet pT and Y spectra are to be compared to CDF data:    - C    
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
C     f = 0.75  overlapping fraction
      palg=-1
      r=0.5d0
      ptmin=15
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c Find b decay products among tracks
      do j=1,njets
         ibtag(j)=0
      enddo
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            if(bson(itrackhep(k))) then
               ibtag(jetvec(k))=ibtag(jetvec(k))+1
            endif
         endif
      enddo
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      end

      logical function bson(j)
      implicit none
      integer J
      include   'hepevt.h'
      integer jcurr
      logical bhadr
      jcurr=j
c     This only happens in parton level analysis
      if(abs(idhep(jcurr)).eq.5) then
         bson=.true.
         return
      endif
 1    continue
      bson=.false.
      if(bhadr(idhep(jcurr))) then
         bson=.true.
         return
      endif
      jcurr=jmohep(1,jcurr)
      if(idhep(jcurr).eq.0) then
         bson=.false.
         return
      endif
      goto 1
      end

      logical function bhadr(idhep)
      implicit none
      integer idhep
      integer i1,i2,idigit
      i1=idigit(1,idhep)
      if(i1.eq.1) then
c         is a bottomed meson
         i2=idigit(3,idhep)
      elseif(i1.eq.2) then
c is a bottomed barion
         i2=idigit(5,idhep)
      endif
      if(i2.eq.5) then
         bhadr=.true.
      else
         bhadr=.false.
      endif
      end

      
      logical function sonoftop(j,jtop)
      implicit none
      integer j,jtop
      include   'hepevt.h'
      integer jcurr
      logical bhadr
      jcurr=j
c     This only happens in parton level analysis
 1    continue
      if(abs(idhep(jcurr)).eq.6) then
         sonoftop=.true.
         jtop=jcurr
         return
      endif
      jcurr=jmohep(1,jcurr)
      if(idhep(jcurr).eq.0) then
         sonoftop=.false.
         jtop=0
         return
      endif
      goto 1
      end

      
      function idigit(k,l)
      implicit none
      integer idigit,k,l
      idigit=abs(mod(l,10**k)/10**(k-1))
      end

