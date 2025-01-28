      subroutine setborn(p,bflav,born,bornjk,bmunu)
c Wrapper subroutine to call the MadGraph Borns
c and set the event-by-event couplings constant
      implicit none
      include 'nlegborn.h'
      real * 8 p(0:3,nlegborn),bornjk(nlegborn,nlegborn)
      integer bflav(nlegborn),sflav(nlegborn)
      real * 8 bmunu(0:3,0:3,nlegborn),born
      integer j,k,l
      real * 8 born_nlojet,ratio
      real * 8 powheginput
      external powheginput
      logical ini,bornonly
      data ini/.true./
      save bornonly,ini
      if (ini) then
         ini=.false.
         bornonly=(powheginput("#bornonly").eq.1d0)
      endif

c     if bornonly, it's faster to use NLOJET to compute the Born contributions.
c     bornjk and bmunu are NOT needed
      if (bornonly) then
         call setborn_nlojet(p,bflav,born_nlojet)
         born_nlojet=abs(born_nlojet)
         born=born_nlojet
         return
      endif


      call set_ebe_couplings

c convert bflav to an equivalent flavor structure with only u's and d's
      call convert_flav(bflav,sflav)
 
      call sborn_proc(p,sflav,born,bornjk,bmunu)

c      write(*,*) bflav,born,bornjk,bmunu

c      write(*,*) 'Born ',born


      
      end

      subroutine borncolour_lh
c Wrapper subroutine to call the MadGraph code to associate
c a (leading) color structure to an event.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      include 'LesHouches.h'
      integer bflav0(nlegborn),bflav(nlegborn),color(2,nlegborn)
      integer i,j,iborn
      logical samecol,conjcol
c We should reach the madgraph flavour configuration that
c was actually computed, in case smartsig is on
      iborn = rad_ubornidx
      bflav0 = flst_born(:,iborn)
      do while(equivto(iborn).ne.-1)
         iborn=equivto(iborn)
      enddo
      bflav = flst_born(:,iborn)
      call born_color(bflav,color)
c Now we have the colour configuration associated with the
c amplitude that was computed instead of rad_ubornidx.
c However, that amplitude may differ from the original one
c by charge conjugation. Check if this is the case
      call matchcolour(nlegborn,bflav0,color)
      icolup(:,1:nlegborn)=color(:,1:nlegborn)
      end


      subroutine matchcolour(n,flav,colour)
      implicit none
      integer n,flav(n),colour(2,n)
      integer j,k,is,iwhere
c Check that color assignment is compatible with flavour.
c If not, conjugate the colour in the whole colour chain.
c This can be done consistently if the two ends of the colour chain
c have opposite flavour.
c We assumes that the colour array is a self consistent colour assignment.

c first conjugate initial flavours and colours, in order to handle
c all particles as outgoing.
      flav(1)=-flav(1)
      flav(2)=-flav(2)
      call colour_conj(colour(:,1))
      call colour_conj(colour(:,2))

      do k=1,n
         if(flav(k).ne.0) then
            if(flav(k).gt.0) then
               is = 1
            else
               is = 2
            endif
c The following is normal assignment
            if(colour(is,k).gt.0.and.colour(3-is,k).eq.0) cycle
c the following is impossible to handle:
            if(colour(1,k).eq.0.and.colour(2,k).eq.0) goto 999
            if(colour(1,k).ne.0.and.colour(2,k).ne.0) goto 999
c the remaining possibility is colour(is,k)=0 and colour(3-is,k)!=0;
c Correct by exchanging colour and anticolour in the whole colour connected
c chain, provided the last element has opposite sign of flavour
            iwhere = k
 10         continue
            do j=1,n
               if(j.ne.iwhere) then
                  if(colour(is,j).eq.colour(3-is,iwhere)) then
                     if(flav(j).eq.0) then
                        if(colour(1,j).ne.0.and.colour(2,j).ne.0) then
                           call colour_conj(colour(:,iwhere))
                           iwhere = j
                           goto 10
                        else
c this cannot be handled
                           goto 999
                        endif
                     elseif(flav(j)*flav(k).lt.0) then
c found a connected opposite flavour;
c The following can't be, colour(is,j)!=0;
                        if(colour(3-is,j).ne.0) goto 999
                        call colour_conj(colour(:,iwhere))
                        call colour_conj(colour(:,j))
                        exit
                     else
                        goto 999
                     endif
                  endif
               endif
            enddo
         endif
      enddo

c Straighten up colour and flavour
      flav(1)=-flav(1)
      flav(2)=-flav(2)
      call colour_conj(colour(:,1))
      call colour_conj(colour(:,2))

      return

 999  continue
      write(*,*)
     1     ' matchcolour: incompatible colour-flavour configuration',
     2     ' better to switch off smartsig in sigborn.f ...'
      call pwhg_exit(-1)
      end



      subroutine finalize_lh
c Specify here if resonances need be written in the event file.
      implicit none

c     The general reshuffling procedure.
      call lhefinitemasses

      end





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC            NLOJET++         CCCCCCCCCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine setborn_nlojet(p,bflav,born)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn)
      integer bflav(nlegborn)
      real * 8 born
      
      logical initnlojet
      data initnlojet/.true./
      save initnlojet
      integer nsubproc
      parameter (nsubproc=4)
      integer virtptrn(nlegborn,nsubproc)
      real * 8 prefactor(nsubproc)
      save virtptrn,prefactor
      integer iptrn,i
      real * 8 pmom(0:4,nlegborn)
      real * 8 prefact

      if (initnlojet) then
c     g g g g g:
         virtptrn(1,1) =  0
         virtptrn(2,1) =  0
         virtptrn(3,1) =  0
         virtptrn(4,1) =  0
         virtptrn(5,1) =  0
         prefactor(1)  =  384d0
c     q q~ g g g:
         virtptrn(1,2) =  1
         virtptrn(2,2) = -1
         virtptrn(3,2) =  0
         virtptrn(4,2) =  0
         virtptrn(5,2) =  0
         prefactor(2)  =  54d0
c     q q~ q q~ g:
         virtptrn(1,3) =  1
         virtptrn(2,3) = -1
         virtptrn(3,3) =  1
         virtptrn(4,3) = -1
         virtptrn(5,3) =  0
         prefactor(3)  =  9d0
c     q q~ Q Q~ g:
         virtptrn(1,4) =  1
         virtptrn(2,4) = -1
         virtptrn(3,4) =  2
         virtptrn(4,4) = -2
         virtptrn(5,4) =  0
         prefactor(4)  =  9d0        
c     print out the patterns
         goto 100
         write(*,*) "**************************************************"
         write(*,*) "*                                                *"
         write(*,*) "*       Initializing the born patterns...     *"
         write(*,*) "*                                                *"
         write(*,*) "**************************************************"
         print *,"The patterns used for the born part: "
         call PrintPatterns(nlegborn,nsubproc,virtptrn)
 100     continue
         initnlojet = .false.         
      endif
      
      born = 0d0

c     fill up the pmom array
      do i=1,nlegborn
         pmom(0:3,i) = p(0:3,i)
         pmom(4,i)   = kn_masses(i)
      enddo
      
c     reshuffle the momenta and select the nlojet subprocess for the
c     current flavour configuration
      call reshufflemom(nlegborn,bflav,pmom,nsubproc,virtptrn,
     $     iptrn,prefact)
      

c     call the nlojet wrapper
      call nlojet3jetborn(iptrn-1,nlegborn,pmom,born)
      
c     fix all missing factors
      born = born*prefactor(iptrn)/prefact
      born = born*(4.d0*pi*st_alpha)**3

      end




      subroutine convert_flav(vflavin,vflavout)
      implicit none
      include 'nlegborn.h'
      integer vflavin(nlegborn),vflavout(nlegborn)
      integer jj,j,k,tmpind,tmpfl

      vflavout=vflavin
      tmpind=0
c Madgraph wants increasing order in final state flavour
      do jj=1,nlegborn
         if(jj.eq.1) then
            j=3
         elseif(jj.eq.2) then
            j=4
         elseif(jj.eq.3) then
            j=5
         elseif(jj.eq.4) then
            j=1
         elseif(jj.eq.5) then
            j=2
         endif
         tmpfl=abs(vflavout(j))
         if(tmpfl.ne.0.and.tmpfl.lt.1000) then
            tmpind=tmpind+1000
            do k=1,nlegborn
               if(abs(vflavout(k)).eq.tmpfl) then
                  vflavout(k)=sign(tmpind,vflavout(k))
               endif
            enddo
         endif
      enddo

      vflavout = vflavout/1000

      end

