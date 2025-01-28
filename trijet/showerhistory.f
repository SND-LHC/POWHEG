      subroutine showerhistory(ntr,clusters,m)
      implicit none
      integer ntr,m
      real * 8 clusters(4,ntr)
c Constructs a clusterization of the input clusters. Clusters that are
c given zero 4 momentum are ignored. Recombination is as in POWHEG.
c ntr: number of initial clusters
c clusters(4,ntr): clusters four momenta: x y z t
c m: desired final number of clusters
      real * 8 ptmin
      parameter (ptmin=0.1d-10)
      real * 8 dst,dsttmp
      integer ip0,ip1,ip2,l,j,k
      real * 8 distance, distance_beams,r,r2
      parameter (r=0.65d0,r2=r*r)

      do j=1,ntr
         if(abs(clusters(1,j)) < ptmin.and.abs(clusters(2,j)) < ptmin)
     1        clusters(:,j) = 0
      enddo

c compress the array, putting zero entries at the end
      call compress_tracks(ntr,clusters)

      do l=1,ntr-m
         dst=-1
         ip1=0
         ip2=0
         ip0=0
         do j=1,ntr
            if(clusters(4,j) == 0) cycle
            if(dst.lt.0) then
               dst =  distance_beams(clusters(:,j))
               ip0=j
               ip1=0
               ip2=0
            else
               dsttmp =  distance_beams(clusters(:,j))
               if(dsttmp.lt.dst) then
                  dst=dsttmp
                  ip0=j
                  ip1=0
                  ip2=0
               endif
            endif
            do k=j+1,ntr
               if(clusters(4,k) == 0) cycle
               if(dst.lt.0) then
                  dst = distance(clusters(:,j),clusters(:,k))/r2
                  ip1=j
                  ip2=k
                  ip0=0
               else
                  dsttmp = distance(clusters(:,j),clusters(:,k))/r2
                  if(dsttmp.lt.dst) then
                     dst=dsttmp
                     ip1=j
                     ip2=k
                     ip0=0
                  endif
               endif
            enddo
         enddo
         
         if(ip1.ne.0) then
            call mergefs(ntr,ip1,ip2,clusters)
c            write(*,*) ' merge fsr '
         elseif(ip0.ne.0) then
            call merge_isr(ntr,ip0,clusters)
c            write(*,*) ' merge isr '
         else
            write(*,*)' no pair found!'
            stop
         endif
      enddo

c compress the array, putting zero entries at the end
      call compress_tracks(ntr,clusters)

      if(m.ge.ntr) then
         m=ntr
         return
      endif         

      end


      subroutine compress_tracks(ntr,clusters)
      implicit none
      integer ntr
      double precision clusters(4,ntr)
      integer j,m
c count the non-zero elements
      m=0
      do j=1,ntr
         if(clusters(4,j) /= 0) then
            m=m+1
         endif
      enddo

      j=0
 1    continue
      j=j+1
c backup ntr until we find a non-zero element
      do while(clusters(4,ntr) == 0 .and. ntr > 0)
         ntr = ntr - 1
      enddo
c now either cluster(4,ntr) /= 0 or ntr = 0 
      if(j >= ntr) goto 3
c now cluster(4,ntr) /= 0 and ntr > 0 (because j<ntr) 
      if(clusters(4,j) == 0) then
         clusters(:,j) = clusters(:,ntr)
         clusters(:,ntr) = 0
      endif
      goto 1

 3    continue

c sanity check:
      if(ntr /= m) then
         write(*,*) ' compress_tracks error!'
         stop
      endif
      do j=1,ntr
         if(clusters(4,j) == 0) then
            write(*,*) ' compress_tracks error!'
            stop
         endif
      enddo

      end

      function distance(p1,p2)
      implicit none
      real * 8 distance
      real * 8 p1(4),p2(4)
      real * 8 y1,y2,dy,phi1,phi2,dphi,ptsq1,ptsq2
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      y1=0.5d0*log((p1(4)+p1(3))/(p1(4)-p1(3)))
      y2=0.5d0*log((p2(4)+p2(3))/(p2(4)-p2(3)))
      dy=abs(y1-y2)
      phi1=atan2(p1(2),p1(1))
      phi2=atan2(p2(2),p2(1))
      dphi=abs(phi1-phi2)
      if(dphi.gt.pi) dphi=abs(2*pi-dphi)
      ptsq1=p1(1)**2+p1(2)**2
      ptsq2=p2(1)**2+p2(2)**2
      distance = (min(ptsq1,ptsq2)) * (dy**2+dphi**2)
      end


      function distance_beams(p)
      implicit none
      real * 8 distance_beams
      real * 8 p(4)
      distance_beams = p(1)**2+p(2)**2
      end


      subroutine mergefs(ntr,j,k,clusters)
      implicit none
      integer ntr,j,k
      real * 8 clusters(4,ntr)
      real * 8 ptot(4),krec(4),q2,ptotu2,ptotu,beta,vec(3),
     1     k0rec,kurec2,kurec,vec1(3),beta1,ptsum(2)
      integer l
      
      ptot = sum(clusters(:,1:ntr),dim=2)

c boost to cm frame

      ptotu2 = ptot(1)**2+ptot(2)**2+ptot(3)**2

      ptotu = sqrt(ptotu2)
      beta = - ptotu/ptot(4)
      vec = ptot(1:3)/ptotu

      call mboost14(ntr,vec,beta,clusters,clusters)
c verify ptot timelike
      call mboost14(1,vec,beta,ptot,ptot)
c      q2 = ptot(4)**2 -ptotu2
      q2 = ptot(4)**2
c      ptot(4) = sqrt(q2)
c      ptot(1:3) = 0
      
c Start recombination

      clusters(:,j) = clusters(:,j) + clusters(:,k)
      clusters(:,k) = 0

c now compute boost on recoil system

      krec = ptot - clusters(:,j)
      k0rec = krec(4)
      kurec2 = krec(1)**2+krec(2)**2+krec(3)**2
      kurec = sqrt(kurec2)

      vec1 = krec(1:3)/kurec
      beta1 = (q2 - (k0rec+kurec)**2)/(q2 + (k0rec+kurec)**2)

      call mboost14(ntr,vec1,beta1,clusters,clusters)
      call mboost14(1,vec1,beta1,krec,krec)

c now fix the jth momentum

      clusters(1:3,j) = -krec(1:3)
      clusters(4,j) = sqrt(krec(1)**2+krec(2)**2+krec(3)**2)

c      write(*,*) sqrt(q2)/(clusters(4,j)+krec(4))

c OK; now boost back

      beta = -beta
      call mboost14(ntr,vec,beta,clusters,clusters)

c Verify pt balance
      ptsum=sum(clusters(1:2,:),dim=2)

      if(abs(ptsum(1))+abs(ptsum(2)).gt.1d-3) then
         write(*,*) ' unbalanced pt: ',ptsum(1),ptsum(2)
      endif

      end



      subroutine merge_isr(ntr,j,clusters)
      implicit none
      integer ntr,j
      double precision clusters(4,ntr)
      real * 8 ptot(4),vec(3),beta,pt2,pt,vec1(3),beta1,ptsum(2)
      clusters(:,j)=0
      ptot = sum(clusters(:,1:ntr), dim = 2)
c longitudinall boost
      beta = -ptot(3)/ptot(4)
      vec = (/ 0d0,0d0,1d0 /)
c 
      call mboost14(ntr,vec,beta,clusters,clusters)
      call mboost14(1,vec,beta,ptot,ptot)
c transverse boost
      pt2 = ptot(1)**2+ptot(2)**2
      pt=sqrt(pt2)
      vec1 = ptot(1:3)/pt
      beta1 = - pt/ptot(4)
      call mboost14(ntr,vec1,beta1,clusters,clusters)
      call mboost14(1,vec1,beta1,ptot,ptot)
c boost back longitudinal
      beta = -beta
      call mboost14(ntr,vec,beta,clusters,clusters)
c Verify pt balance
      ptsum=sum(clusters(1:2,1:ntr),dim=2)

      if(abs(ptsum(1))+abs(ptsum(2)).gt.1d-3) then
         write(*,*) ' unbalanced pt: ',ptsum(1),ptsum(2)
      endif

      end

      
      subroutine mboost14(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(1:4,m),vout(1:4,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
      enddo
      end

