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

      integer,parameter :: nybins=6, nenubins=20
      integer :: i
      character * 10 :: suffix(nybins)
      real * 8 :: enubin(nenubins+1)
      common/enuarray/enubin
      common/suffix/suffix
      character *4 tmpsuff

      suffix(1) = '7-7.5'
      suffix(2) = '7.5-8'
      suffix(3) = '8-8.5'
      suffix(4) = '8.5-9'
      suffix(5) = '9-9.5'
      suffix(6) = '9.5-inf'

      enubin = [0.d0,100.d0,200.d0,300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,1000.d0,1100.d0,1200.d0,1300.d0,1400.d0,1500.d0,1600.d0,1700.d0,1800.d0,1900.d0,2000.d0]
      
      call inihists
            
      do i=1,6
         call bookup('e-charm'//trim(adjustl(suffix(i))),nenubins,enubin)
         call bookup('enu-nue'//trim(adjustl(suffix(i))),nenubins,enubin)
         call bookup('enu-numu'//trim(adjustl(suffix(i))),nenubins,enubin)
         call bookup('enu-nutau'//trim(adjustl(suffix(i))),nenubins,enubin)
         call bookup('enu-nueb'//trim(adjustl(suffix(i))),nenubins,enubin)
         call bookup('enu-numub'//trim(adjustl(suffix(i))),nenubins,enubin)
         call bookup('enu-nutaub'//trim(adjustl(suffix(i))),nenubins,enubin)
      enddo
      end

      subroutine analysis(dsig0)
      implicit none
      integer, parameter :: dsigdim = 200
      real * 8  dsig0,dsig(dsigdim)
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      logical ini
      data ini/.true./
      save ini
      integer ihep
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/

      integer,parameter :: nybins=6, nenubins=19
      integer :: i
      character * 10 :: suffix(nybins)
      real * 8 :: enubin(nenubins+1)
      common/enuarray/enubin
      common/suffix/suffix

      real * 8 :: ptD,yD,etaD
      character *3 tmpsuff

      real * 8 :: xextrap,yextrap
      
      call multi_plot_setup(dsig0,dsig,dsigdim)

      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         do ihep=3,nhep
            
            if (isthep(ihep) == 1 .and. idhep(ihep) == 4 ) then 
               call ptyeta(phep(1,ihep),ptD,yD,etaD)
               if (etaD > 7d0 .and. etaD<7.5d0) then
                  call filld('e-charm'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
               endif
               if (etaD > 7.5d0 .and. etaD<8d0) then
                  call filld('e-charm'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
               endif
               if (etaD > 8d0 .and. etaD<8.5d0) then
                  call filld('e-charm'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
               endif
               if (etaD > 8.5d0 .and. etaD<9.d0) then
                  call filld('e-charm'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
               endif
               if (etaD > 9.d0 .and. etaD<9.5d0) then
                  call filld('e-charm'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
               endif
               if (etaD > 9.5d0) then
                  call filld('e-charm'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
               endif
            endif
         enddo
         
      else 
!     with showers 

         do ihep=1,nhep
            
            if(isthep(ihep).eq.1) then
               
               yextrap = phep(2,ihep)/abs(phep(3,ihep))*48000.d0
               if (yextrap > -6.d0) then
                  if (yextrap < 34.d0) then
                     xextrap = phep(1,ihep)/abs(phep(3,ihep))*48000.d0
                     if (xextrap > -30.d0) then
                        if (xextrap < 10.d0) then
                           
                           if(idhep(ihep).eq.12) then
                              
                              call ptyeta(phep(1,ihep),ptD,yD,etaD)

                              if (etaD > 7.d0 .and. etaD<7.5d0) then
                                 call filld('enu-nue'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 7.5d0 .and. etaD<8.d0) then
                                 call filld('enu-nue'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.d0 .and. etaD<8.5d0) then
                                 call filld('enu-nue'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.5d0 .and. etaD<9.d0) then
                                 call filld('enu-nue'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.d0 .and. etaD<9.5d0) then
                                 call filld('enu-nue'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.5d0) then
                                 call filld('enu-nue'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
                              endif
                           endif

                           if(idhep(ihep).eq.14) then
                              
                              call ptyeta(phep(1,ihep),ptD,yD,etaD)

                              if (etaD > 7.d0 .and. etaD<7.5d0) then
                                 call filld('enu-numu'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 7.5d0 .and. etaD<8.d0) then
                                 call filld('enu-numu'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.d0 .and. etaD<8.5d0) then
                                 call filld('enu-numu'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.5d0 .and. etaD<9.d0) then
                                 call filld('enu-numu'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.d0 .and. etaD<9.5d0) then
                                 call filld('enu-numu'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.5d0) then
                                 call filld('enu-numu'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
                              endif
                           endif
                           
                           if(idhep(ihep).eq.16) then
                              
                              call ptyeta(phep(1,ihep),ptD,yD,etaD)

                              if (etaD > 7.d0 .and. etaD<7.5d0) then
                                 call filld('enu-nutau'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 7.5d0 .and. etaD<8.d0) then
                                 call filld('enu-nutau'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.d0 .and. etaD<8.5d0) then
                                 call filld('enu-nutau'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.5d0 .and. etaD<9.d0) then
                                 call filld('enu-nutau'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.d0 .and. etaD<9.5d0) then
                                 call filld('enu-nutau'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.5d0) then
                                 call filld('enu-nutau'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
                              endif
                           endif

                           if(idhep(ihep).eq.-12) then
                              
                              call ptyeta(phep(1,ihep),ptD,yD,etaD)

                              if (etaD > 7.d0 .and. etaD<7.5d0) then
                                 call filld('enu-nueb'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 7.5d0 .and. etaD<8.d0) then
                                 call filld('enu-nueb'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.d0 .and. etaD<8.5d0) then
                                 call filld('enu-nueb'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.5d0 .and. etaD<9.d0) then
                                 call filld('enu-nueb'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.d0 .and. etaD<9.5d0) then
                                 call filld('enu-nueb'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.5d0) then
                                 call filld('enu-nueb'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
                              endif
                           endif

                           if(idhep(ihep).eq.-14) then
                              
                              call ptyeta(phep(1,ihep),ptD,yD,etaD)

                              if (etaD > 7.d0 .and. etaD<7.5d0) then
                                 call filld('enu-numub'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 7.5d0 .and. etaD<8.d0) then
                                 call filld('enu-numub'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.d0 .and. etaD<8.5d0) then
                                 call filld('enu-numub'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.5d0 .and. etaD<9.d0) then
                                 call filld('enu-numub'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.d0 .and. etaD<9.5d0) then
                                 call filld('enu-numub'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.5d0) then
                                 call filld('enu-numub'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
                              endif
                           endif
                           
                           if(idhep(ihep).eq.-16) then
                              
                              call ptyeta(phep(1,ihep),ptD,yD,etaD)

                              if (etaD > 7.d0 .and. etaD<7.5d0) then
                                 call filld('enu-nutaub'//trim(adjustl(suffix(1))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 7.5d0 .and. etaD<8.d0) then
                                 call filld('enu-nutaub'//trim(adjustl(suffix(2))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.d0 .and. etaD<8.5d0) then
                                 call filld('enu-nutaub'//trim(adjustl(suffix(3))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 8.5d0 .and. etaD<9.d0) then
                                 call filld('enu-nutaub'//trim(adjustl(suffix(4))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.d0 .and. etaD<9.5d0) then
                                 call filld('enu-nutaub'//trim(adjustl(suffix(5))), phep(4,ihep), dsig)
                              endif
                              if (etaD > 9.5d0) then
                                 call filld('enu-nutaub'//trim(adjustl(suffix(6))), phep(4,ihep), dsig)
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         enddo
      endif
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
