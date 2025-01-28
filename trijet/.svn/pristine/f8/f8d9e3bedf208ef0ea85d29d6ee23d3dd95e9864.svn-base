C     returns 2 Re(M_B * M_V)/(as/(2pi)), 
C     where M_B is the Born amplitude and 
C     M_V is the finite parte of the virtual amplitude.
C     The as/(2pi) factor is attached at a later point.
C     The virtual amplitude is generated using GoSam. 
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 p(0:3,nlegborn)
      real * 8 oldp(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual
      
      logical initnlojet
      data initnlojet/.true./
      save initnlojet
      integer nsubproc
      parameter (nsubproc=4)
      integer virtptrn(nlegborn,nsubproc)
      real * 8 prefactor(nsubproc)
      save virtptrn,prefactor
      integer iptrn,i
      real * 8 mur
      real * 8 pmom(0:4,nlegborn)
      real * 8 virt,tree
      real * 8 prefact
      real * 8 det, ecmsq, Q(4,4)
      real * 8 virtual_gosam, ratio
      integer j,k,l,m
      real * 8 dotp,sss(10)
      data oldp/20*0d0/
      save oldp


c Look for almost identical invariants (they cause troubles in
c the virtual routines). Set the virtual to 0 if this is the case.
      if(.not.flg_in_smartsig) then
         l=1
         do j=1,5
            do k=j+1,5
               sss(l) = dotp(p(:,j),p(:,k))
               do m=l-1,1,-1
                  if(abs(abs(sss(l))-abs(sss(m)))/
     1                 (abs(sss(l))+abs(sss(m))).lt.1d-6) then
                     virtual = 0
                     if(abs(sum(oldp-p)/sum(oldp+p)).gt.0d0) then
c                        write(*,*)
c     1                  'Skipped risky kinematic point in setvirtual'
                        call increasecnt
     1                       ("Risky kinematic points in setvirtual")
                        oldp = p
                     endif
                     return
                  endif
               enddo
               l=l+1
            enddo
         enddo

c  Check if the 5 momenta are independent. Otherwise return virtual equal 0
c  Construct a matrix Q (4x4) from invariants such that
c  Q_ij = dot(p_i,p_j).
c  If det(Q)=0, then p3 or p4 are parallel to the beam axis or the event is planar         
         do i=1,3
            do j=i+1,4
               Q(i,j) = dotp(p(:,i),p(:,j))
            enddo
         enddo         
         det = Q(1,2)**2*Q(3,4)**2 - 2*Q(1,2)*Q(2,3)*Q(3,4)*Q(1,4) 
     $        - 2*Q(1,2)*Q(2,4)*Q(1,3)*Q(3,4) + Q(1,3)**2*Q(2,4)**2 
     $        - 2*Q(1,3)*Q(2,4)*Q(2,3)*Q(1,4) + Q(2,3)**2*Q(1,4)**2         
         ecmsq = 2*Q(1,2)
c         write(*,*) det/ecmsq**4
         if (abs(det/ecmsq**4).lt.1d-14) then
            virtual=0
            if(abs(sum(oldp-p)/sum(oldp+p)).gt.0d0) then
c               write(*,*)
c     1              ' Skipped risky kinematic point. Planar config'
               call increasecnt
     1              ("Risky kinematic points. Planar config")   
               oldp = p
            endif
            return
         endif
      endif

      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC              GOSAM          CCCCCCCCCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      call setvirtual_gosam(p,vflav,virtual_gosam)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC            NLOJET++         CCCCCCCCCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
         write(*,*) "*       Initializing the virtual patterns...     *"
         write(*,*) "*                                                *"
         write(*,*) "**************************************************"
         print *,"The patterns used for the virtual part: "
         call PrintPatterns(nlegborn,nsubproc,virtptrn)
 100     continue
         initnlojet = .false.         
      endif
      
      virtual = 0d0

c     fill up the pmom array
      do i=1,nlegborn
         pmom(0:3,i) = p(0:3,i)
         pmom(4,i)   = kn_masses(i)
      enddo
      
c     reshuffle the momenta and select the nlojet subprocess for the
c     current flavour configuration
      call reshufflemom(nlegborn,vflav,pmom,nsubproc,virtptrn,
     $     iptrn,prefact)
      
c     set the renormalization scale
      mur = sqrt(st_muren2)

c     call the nlojet wrapper
      call nlojet3jetvirt(iptrn-1,nlegborn,pmom,mur,st_nlight,virt)
      
c     fix all missing factors
      virt = virt*prefactor(iptrn)/prefact
      virtual = virt*2.d0*pi*(4.d0*pi*st_alpha)**3*4d0*pi

      end




! module reshuffle_data
! implicit none
! !
!   integer :: firstlight
!   integer :: nlight
!   integer , parameter :: maxpart = 20
!   integer , parameter :: maxpartons = 9
! !
! contains
! !
      
      subroutine PrintPatterns(n,numptrn,ptrn)
      implicit none
      integer n,numptrn
      integer ptrn(n,numptrn)
      integer iptrn,ipart      
!     We print out the patterns:
      print 100,"The number of patterns we have: ",numptrn
 100  format(A,I2)
      do iptrn=1,numptrn
         print 101,"The ",iptrn,"th pattern is: ",ptrn(:,iptrn)
 101     format(A,I2,A,20(I3,1X))
      end do
      end subroutine PrintPatterns




      subroutine reshufflemom(n,flav,p,numptrn,ptrn,
     $     returnptrn,prefactor)      
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'      
c      include 'reshuffle_data.h'     
      integer n,numptrn
      integer flav(n)
      real * 8 p(0:4,n)
      integer ptrn(n,numptrn)
      integer returnptrn
      real * 8  prefactor
!     
c      integer flavin(maxpartons) 
c      integer flav_ptrn(maxpartons)
c      integer  io(maxpartons),ordering(maxpartons)

      integer flavin(n) 
      integer flav_ptrn(n)
      integer io(n),ordering(n)


      integer ipart,iparton,jparton,iptrn
      integer nparton
      integer nextparton,partonid
c      integer partonlist(9)
      integer partonlist(n)
      logical found
      real * 8 prefact
c      real * 8 pmom_tmp(0:4,maxpart)
      real * 8 pmom_tmp(0:4,n)
      integer firstlight,nlight

      firstlight=flst_lightpart
      nlight= st_nlight
      nextparton = 1
      iparton = 0
      partonlist = 0
      returnptrn = -1000
      prefactor = 0d0
!     print *,flav(:)
!     We encode the input flavour configuration with our scheme:
      do ipart=1,n
!     print *,"flav(ipart) is: ",flav(ipart)
!     skip heavy and/or non-QCD stuff:
         if ((ipart.gt.2).and.(ipart.lt.firstlight)) cycle
!     print *,"iparton is increased..."
         iparton = iparton + 1
!     we select light, non-gluon partons:
         if ((flav(ipart).ne.0).and.(abs(flav(ipart)).le.nlight)) then
!      print *,"quark!"
!     We presume that the partontype is not yet stored in partonlist:
            partonid = nextparton
!     We search for it in the partonlist:
            do jparton=1,nextparton-1
! The presumption was wrong we found it, so use jparton for id:
               if (partonlist(jparton).eq.abs(flav(ipart))) then
!     print *,"We found the quark in the list!"
                  partonid = jparton
                  exit
               end if
            end do
!     If we did not find the parton in the list we have to insert it:
            if (partonid.eq.nextparton) then
               nextparton = nextparton + 1
               partonlist(partonid) = abs(flav(ipart))
            end if
!     In any case we insert it into the flavour list:
            flavin(iparton) = sign(partonid,flav(ipart))
!     If it is a gluon we simply include a 'g':
         elseif (flav(ipart).eq.0) then
!     print *,"gluon!"
            flavin(iparton) = 0
         end if
      end do
!     In the following we only book-keep the massless partons:
      nparton = iparton
!     print 200,flavin(1:nparton)
!     
!     At this point we bear the flavor configuration in an internal 
!     representation.
!     We have to take all the available patterns and compare flavin with 
!     them:
      do iptrn=1,numptrn
         flav_ptrn(1:n) = ptrn(1:nparton,iptrn)
         call CompFlavors(nparton,flav_ptrn,flavin,io,ordering,found)
!     print *,"iptrn, found: ",iptrn,found
         if (found) then
!     print 200,ordering(1:nparton)
!     print 200,io(1:nparton)
            exit
         end if
      end do
!     
! For safety measures we check wether we indeed found the pattern:
      if (.not.found) then
         print*,"We failed to found the pattern for: ",flavin(1:nparton)
         stop
      end if
!     
!     Since we can have non-QCD and/or massive partons we have to
!     modify the numbers in ordering if they are larger than 2:
      where (ordering > 2)
         ordering = ordering + firstlight - 3
      end where
!     print *,"The ordering when nonQCD part is considered: "
!     print 200,ordering(1:nparton)
!     
!     We have to take the original pmom array and reshuffle the momenta
!     using a previously specified convention:
!     We can have non-QCD and/or heavy particles, for these we have to 
!     simply copy them:
      if (firstlight.ne.3) then
         pmom_tmp(:,3:firstlight-1) = p(:,3:firstlight-1)
      end if
!     For all the massless QCD parton we reshuffle:
!     -io adds a minus sign if the particle is in the initial state:
      do ipart=1,nparton
         if (ipart.lt.3) 
     $        pmom_tmp(:,ipart) = -io(ipart)*p(:,ordering(ipart))
         if (ipart.ge.3) pmom_tmp(:,ipart + firstlight - 3) = 
     $        -io(ipart)*p(:,ordering(ipart))
      end do
!
!     We copy back pmom_tmp:
      p(:,1:2) = -pmom_tmp(:,1:2)
      p(:,3:n) = pmom_tmp(:,3:n)
!     
! We calculate the needed new prefactor as well:
      call CalcPrefactor(nparton,flavin,prefact)
      prefactor = prefact

      returnptrn = iptrn

 200  format(20(I2,1x))
     
      end subroutine reshufflemom




      subroutine CompFlavors(npart,flavA,flavB,io,ordering,found)      
      implicit none
c      include 'reshuffle_data.h'
      integer npart
c      integer flavA(maxpartons),flavB(maxpartons)
c      integer io(maxpartons),ordering(maxpartons)
      integer flavA(npart),flavB(npart)
      integer io(npart),ordering(npart)
      logical found

      integer ipart,jpart
      integer flv1,flv2
c      integer flav_tmp(maxpartons)
      integer flav_tmp(npart)
      
      found = .false.
      io(1:2) = 1
      io(3:npart) = -1
      ordering = 0
!     
!     print *,"The input flavor configs are: "
!     print 110,"flavA : ",flavA(1:npart)
!     print 110,"flavB : ",flavB(1:npart)
!     
!     We keep flavA and we try to shuffle the flavors in flavB to arrive at
!     flavA:
      flav_tmp = flavB
!     
!     We go through all partons in flavA and try to find the same in flav_tmp:
      do ipart=1,npart
         do jpart=1,npart+1
            if (jpart.eq.npart+1) return
            if (flav_tmp(jpart).eq.99) cycle
            if (((ipart.lt.3).and.(jpart.gt.2)).or. 
     $           ((ipart.gt.2).and.(jpart.lt.3))) then
            flv1 = -flavA(ipart)
         else
            flv1 =  flavA(ipart)
         end if
         flv2 = flav_tmp(jpart)
         if (flv1.eq.flv2) then
            ordering(ipart) = jpart
            flav_tmp(jpart) = 99
            exit
         end if
      end do
      end do
!     
!     When we found a good pattern we have to save the io array as well,
!     which is used in HELAC but can be useful even in NJETS too:
      do ipart=1,npart
         if (((ipart.lt.3).and.(ordering(ipart).ge.3)).or.
     $        ((ipart.ge.3).and.(ordering(ipart).lt.3))) then
         io(ipart) = -io(ipart)
      end if
      end do
!     
      found = .true.
!     
 100  format(20(I2,1x))
 110  format(A,20(I2,1x))
!     
      end subroutine CompFlavors



!
!     This routine presents a simple way how to recalculate the symmetry
!     factors associated to the colour factor of initial states and the
!     factor associated to identical partons in the final state. The 
!     factor is calculated for each evaluation, so the efficiency can
!     be at least argued, but this routine is intended to be used for
!     complex processes, that is, the CPU time needed to perform these
!     calculations is negligible to calculations related to the dynamics.
      subroutine CalcPrefactor(npart,flav,prefact)
      implicit none
c      include 'reshuffle_data.h'
      integer npart
c      integer flav(maxpartons)
      integer flav(npart)
      real * 8  prefact
!
      integer ipart
c      integer numparton(-maxpartons:maxpartons)
      integer numparton(-npart:npart)
!
      real * 8 factorial
      external factorial
!
      prefact = 1d0
      numparton = 0
!     We loop through all final state particles and count the
!     number of each flavour.
      do ipart=3,npart
         numparton(flav(ipart)) = numparton(flav(ipart)) + 1
      end do
!     
!     print *,"Number of different partons: "
!     print *,numparton
!     
c      do ipart=-maxpartons,maxpartons
      do ipart=-npart,npart
         if (numparton(ipart).ne.0) then
            prefact = prefact*factorial(numparton(ipart))
         endif
      end do
!     
!     We also include the colour factor associated to the initial state:
      do ipart=1,2
         if (flav(ipart).ne.0) prefact = prefact*3d0
         if (flav(ipart).eq.0) prefact = prefact*8d0
      end do
!     
!     print *,"The prefactor is: ",prefact
!     
      end subroutine CalcPrefactor
!

      function factorial(n) 
      implicit none
      integer  n
      real * 8 factorial    
      real ans
      integer i
      ans = 1d0
      if (n .gt. 1) then
         do i=n,2,-1
            ans=ans*i
         enddo
      endif
      factorial=ans
      end 
      
