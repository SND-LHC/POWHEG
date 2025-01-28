      subroutine py_change_scalup
      implicit none
      include 'LesHouches.h'
c      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 ptmin,ptmin2,pcm(0:3,maxnup),beta,vec(3),ptkj2
      integer k,mu,j
      logical ini
      save ini
      data ini/.true./
      real * 8 dotp
      external dotp
      integer npart
      npart=0
      do k=3,nup
c     only light partons
         if (idup(k).eq.21.or.abs(idup(k)).le.5) then
            npart=npart+1
            do mu=1,3
               pcm(mu,npart)=pup(mu,k)
            enddo
            pcm(0,npart)=pup(4,k)
         endif
      enddo
c     compute min pt of light partons with respect to the incoming beam
      ptmin2=1d30
      do k=1,npart
         ptmin2=min(ptmin2,pcm(1,k)**2+pcm(2,k)**2)
      enddo

c     compute pt's of the final state partons with respect to each other
      beta=-(pup(3,1)+pup(3,2))/(pup(4,1)+pup(4,2))
      vec(1)=0
      vec(2)=0
      vec(3)=1
c     go in the CM frame   
      call mboost(npart,vec,beta,pcm,pcm)
      do k=1,npart-1
         do j=k+1,npart
            ptkj2 = 2*dotp(pcm(0,k),pcm(0,j))*
     $           pcm(0,k)*pcm(0,j)/(pcm(0,k)+pcm(0,j))**2
            ptmin2=min(ptmin2,ptkj2)
         enddo
      enddo
      ptmin=sqrt(ptmin2)
      if(scalup.gt.ptmin) then
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) 'scalup set to the min pt in the event'
            write(*,*) '*************************************'
            ini=.false.
         endif         
         scalup = ptmin
      endif
      end




      subroutine py_change_scalup_bi
      implicit none
      include 'LesHouches.h'
c      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 ptmin,ptmin2,pcm(0:3,maxnup),beta,vec(3),ptkj2
      integer k,mu,j
      logical ini
      save ini
      data ini/.true./
      real * 8 dotp,r,ptj1,ptj2
      external dotp
      integer npart
      npart=0
      do k=3,nup
c     only light partons
         if (idup(k).eq.21.or.abs(idup(k)).le.5) then
            npart=npart+1
            do mu=1,3
               pcm(mu,npart)=pup(mu,k)
            enddo
            pcm(0,npart)=pup(4,k)
         endif
      enddo
c     compute min pt of light partons with respect to the incoming beam
      ptmin2=1d30
      do k=1,npart
         ptmin2=min(ptmin2,pcm(1,k)**2+pcm(2,k)**2)
      enddo

      do k=1,npart-1
         do j=k+1,npart
            ptj1 = pcm(1,k)**2+pcm(2,k)**2
            ptj2 = pcm(1,j)**2+pcm(2,j)**2
            call pwhg_getR_phiy(pcm(:,k),pcm(:,j),R)
            ptmin2=min(ptmin2,r**2*min(ptj1,ptj2))
         enddo
      enddo
      ptmin=sqrt(ptmin2)
      if(scalup.gt.ptmin) then
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) 'scalup set to the min pt in the event'
            write(*,*) '*************************************'
            ini=.false.
         endif         
         scalup = ptmin
      endif
      end
