      subroutine setreal(p,rflav,amp2)
c Wrapper subroutine to call the MadGraph real emission matrix
c elements and set the event-by-event couplings constant
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2,amp2nlojet
      integer proc


      logical initnlojet
      data initnlojet/.true./
      save initnlojet
      integer nsubproc
      parameter (nsubproc=8)
      integer realptrn(nlegreal,nsubproc)
      real * 8 prefactor(nsubproc)
      save realptrn,prefactor
      integer i,iptrn
      real * 8 pmom(0:4,nlegreal)
      real * 8 prefact

      real * 8 amp2_madgraph, ratio

c$$$      goto 10
c$$$      call set_ebe_couplings
c$$$      call sreal_proc(p,rflav,amp2)
c$$$c Cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
c$$$      amp2 = amp2/(st_alpha/(2d0*pi))           
c$$$      amp2_madgraph=amp2
c$$$ 10   continue


      if (initnlojet) then
c     g g g g g g:
         realptrn(1,1) =  0
         realptrn(2,1) =  0
         realptrn(3,1) =  0
         realptrn(4,1) =  0
         realptrn(5,1) =  0
         realptrn(6,1) =  0
         prefactor(1)  =  1536d0
c     q q~ g g g g:
         realptrn(1,2) =  1
         realptrn(2,2) = -1
         realptrn(3,2) =  0
         realptrn(4,2) =  0
         realptrn(5,2) =  0
         realptrn(6,2) =  0
         prefactor(2)  =  216d0
c     q q~ q q~ g g:
         realptrn(1,3) =  1
         realptrn(2,3) = -1
         realptrn(3,3) =  1
         realptrn(4,3) = -1
         realptrn(5,3) =  0
         realptrn(6,3) =  0
         prefactor(3)  =  18d0
c     q q~ Q Q~ g g:
         realptrn(1,4) =  1
         realptrn(2,4) = -1
         realptrn(3,4) =  2
         realptrn(4,4) = -2
         realptrn(5,4) =  0
         realptrn(6,4) =  0
         prefactor(4)  =  18d0
c     q q~ q q~ q q~:
         realptrn(1,5) =  1
         realptrn(2,5) = -1
         realptrn(3,5) =  1
         realptrn(4,5) = -1
         realptrn(5,5) =  1
         realptrn(6,5) = -1
         prefactor(5)  =  36d0
c     q q~ q q~ Q Q~:
         realptrn(1,6) =  1
         realptrn(2,6) = -1
         realptrn(3,6) =  1
         realptrn(4,6) = -1
         realptrn(5,6) =  2
         realptrn(6,6) = -2
         prefactor(6)  =  9d0
c     q q~ Q Q~ R R~:
         realptrn(1,7) =  1
         realptrn(2,7) = -1
         realptrn(3,7) =  2
         realptrn(4,7) = -2
         realptrn(5,7) =  3
         realptrn(6,7) = -3
         prefactor(7)  =  9d0
c     q q~ Q Q~ Q Q~:
         realptrn(1,8) =  1
         realptrn(2,8) = -1
         realptrn(3,8) =  2
         realptrn(4,8) = -2
         realptrn(5,8) =  2
         realptrn(6,8) = -2
         prefactor(8)  =  36d0
c     print out the patterns
         goto 100
         write(*,*) "**************************************************"
         write(*,*) "*                                                *"
         write(*,*) "*       Initializing the real patterns...        *"
         write(*,*) "*                                                *"
         write(*,*) "**************************************************"         
         print *,"The patterns used for the real part: "
         call PrintPatterns(nlegreal,nsubproc,realptrn)
 100     continue
         initnlojet = .false.
      endif

      amp2 = 0d0

c     fill up the pmom array
      do i=1,nlegreal
         pmom(0:3,i) = p(0:3,i)
         pmom(4,i)   = kn_masses(i)
      end do

      
c     reshuffle the momenta and select the nlojet subprocess for the
c     current flavour configuration
      call reshufflemom(nlegreal,rflav,pmom,nsubproc,realptrn,
     $     iptrn,prefact)
      
c     call the nlojet wrapper
      call nlojet3jetreal(iptrn-1,nlegreal,pmom,amp2)

c     fix all missing factors      
      amp2 = amp2*prefactor(iptrn)/prefact
      amp2 = amp2*2.d0*pi*(4.d0*pi*st_alpha)**3*4d0*pi

c      write(*,*) 'real ',amp2

      end


