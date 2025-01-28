      subroutine init_processes
      implicit none
      include "brinclude.h"
      include "pwhg_st.h"
      include "pwhg_par.h"
      include "pwhg_flg.h"
      include "coupl.inc"
      integer i
      integer max_num_flav
      integer flst_bornloc(nlegborn,maxprocborn),
     $     flst_realloc(nlegreal,maxprocreal)
      integer flst_nbornloc,flst_nrealloc, j
      logical keep
      real * 8 powheginput
      external powheginput

      flg_doublefsr = .true.
      if(powheginput("#doublefsr").eq.0) flg_doublefsr=.false.
      flg_evenmaxrat = .true.
      if(powheginput("#evenmaxrat").eq.0) flg_evenmaxrat = .false.

      st_bornorder=3
c MiNLO activated by default
      if(powheginput("#minlo").ne.0) then
         flg_minlo=.true.
         flg_minlo_real=.false.
      else
         flg_minlo=.false.
         flg_minlo_real=.false.
      endif

c     the default value 1d-6 could be too small
      par_isrtinycsi = 1d-4
      par_isrtinyy = 1d-4
      par_fsrtinycsi = 1d-4
      par_fsrtinyy = 1d-4

c
c We read in the optimization parameters:
      brpar_diexp  = 1d0
      brpar_dijexp = 1d0         
      if (powheginput("#brpar_diexp").gt.0d0) then
         brpar_diexp  = powheginput("#brpar_diexp")
      endif
      if (powheginput("#brpar_dijexp").gt.0d0) then
         brpar_dijexp = powheginput("#brpar_dijexp")
      endif
c$$$      Print *,"******************************************"
c$$$      Print *,"*                                        *"
c$$$      Print *,"* diexp parameter is:  ",brpar_diexp
c$$$      Print *,"* dijexp parameter is: ",brpar_dijexp
c$$$      Print *,"*                                        *"
c$$$      Print *,"******************************************"

c     The default values of the following parameters is 1 in the POWHEG BOX.
c     For jets, better results are obtained with:
      par_diexp  = 1
      par_dijexp = 1
      par_2gsupp = 1     


      call init_processes_born
      call init_processes_real

      flst_nbornloc=0
      flst_nrealloc=0



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     here we can reduce the number of considered subprocesses!!
c     The list of subprocesses has been generated with 
c     max_num_flav = 5, i.e. all flavors, including the bottom
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      max_num_flav = 5
      if (max_num_flav.lt.5) then        
         
         write(*,*) '********************************************'
         write(*,*) '      NUMBER OF PROCESSES REDUCED!!!!!      '
         write(*,*) '********************************************'

c     do the Born diagrams
         do i=1,flst_nborn
            keep=.true.
            do j=1,nlegborn
               if (abs(flst_born(j,i)).gt.max_num_flav) keep=.false.
            enddo
            if (keep) then
               flst_nbornloc=flst_nbornloc+1
               flst_bornloc(:,flst_nbornloc)=flst_born(:,i)
            endif
         enddo
         flst_nborn=flst_nbornloc
         flst_born(:,:)=flst_bornloc(:,:)

c     do the real diagrams
         do i=1,flst_nreal
            keep=.true.
            do j=1,nlegreal
               if (abs(flst_real(j,i)).gt.max_num_flav) keep=.false.
            enddo
            if (keep) then
               flst_nrealloc=flst_nrealloc+1
               flst_realloc(:,flst_nrealloc)=flst_real(:,i)
            endif
         enddo
         flst_nreal=flst_nrealloc
         flst_real(:,:)=flst_realloc(:,:)
      endif
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Set here the number of light flavours
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      st_nlight=5
      call init_couplings
c     if (tmass.eq.0d0) then
c        st_nlight=6
c     elseif(bmass.eq.0d0) then
c        st_nlight=5
c     elseif(cmass.eq.0d0) then
c        st_nlight=4
c     else
c        st_nlight=3
c     endif
      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.st_nlight) then
            flst_lightpart=i
            exit
         endif
      enddo
 
      end
