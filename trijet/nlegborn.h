c -*- Fortran -*-      
      integer nlegborn,nlegreal
      parameter (nlegborn=           5 )
      parameter (nlegreal=nlegborn+1)
      
      integer ndiminteg
      parameter (ndiminteg=(nlegreal-2)*3-4+2-1
c     one extra variable to separate ISR and FSR for second jet
     $ +1 )
     
      integer maxprocborn,maxprocreal
      parameter (maxprocborn=300,maxprocreal=900)

      integer maxalr
      parameter (maxalr=maxprocreal*nlegreal*(nlegreal-1)/2)


      
