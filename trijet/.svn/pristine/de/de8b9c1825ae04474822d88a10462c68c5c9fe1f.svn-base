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
      integer i,j,k
      real * 8 powheginput
      external powheginput
c
      integer nnum
      parameter (nnum = 9)
      character*1 cnum(nnum)
      data cnum/'1','2','3','4','5','6','7','8','9'/
c
      integer ndphibins1
      parameter (ndphibins1 = 15)
      integer ndphibins2
      parameter (ndphibins2 = 10)
      integer ndphibins3
      parameter (ndphibins3 = 9)
      integer ndphibins4
      parameter (ndphibins4 = 4)
      integer nptmaxbins
      parameter (nptmaxbins = 9)
      real*8 dphibins1(ndphibins1 + 1)
      real*8 dphibins2(ndphibins2 + 1)
      real*8 dphibins3(ndphibins3 + 1)
      real*8 dphibins4(ndphibins4 + 1)
      real*8 dptmaxbins(nptmaxbins + 1)
      data dphibins1/0.50d0,0.59d0,0.66d0,0.71d0,0.75d0,0.78d0,0.81d0,
     1               0.84d0,0.86d0,0.88d0,0.90d0,0.92d0,0.94d0,0.96d0,
     2               0.98d0,1.00d0/
      data dphibins2/0.50d0,0.66d0,0.75d0,0.81d0,0.86d0,0.90d0,0.92d0,
     1               0.94d0,0.96d0,0.98d0,1.00d0/
      data dphibins3/0.50d0,0.66d0,0.75d0,0.81d0,0.86d0,0.90d0,0.94d0,
     1               0.96d0,0.98d0,1.00d0/
      data dphibins4/0.75d0,0.86d0,0.94d0,0.98d0,1.00d0/
      data dptmaxbins/110d0,160d0,210d0,260d0,310d0,400d0,500d0,600d0,
     1                800d0,10000d0/
      character*3 cptmax(9)
      data cptmax/'110','160','210','260','310','400','500',
     1            '600','800'/
c
      integer nptrangebins
      parameter (nptrangebins = 11)
      real*8 ptrangebins(nptrangebins + 1)
      data ptrangebins/30d0,40d0,60d0,80d0,110d0,160d0,210d0,260d0,
     1                 310d0,400d0,500d0,600d0/
      character*3 cptrange(nptrangebins + 1)
      data cptrange/' 30',' 40',' 60',' 80','110','160','210','260',
     1              '310','400','500','600'/
      integer nyrangebins
      parameter (nyrangebins = 5)
      real*8 yrangebins(nyrangebins + 1)
      data yrangebins/0.0d0,0.3d0,0.8d0,1.2d0,2.1d0,2.8d0/
      character*3 cyrange(nyrangebins + 1)
      data cyrange/'0.0','0.3','0.8','1.2','2.1','2.8'/
      real*8 deltar
      parameter (deltar = 0.1d0)
      real*8 Rparshapes
      parameter (Rparshapes = 0.6d0)
      real*8 MaxJetRapshapes
      parameter (MaxJetRapshapes = 2.8d0)
c
      integer nptbins
      parameter (nptbins = 4)
      real*8 ptbins(nptbins + 1)
      data ptbins/200d0,300d0,400d0,500d0,600d0/
      character*3 cptrangesub(nptbins + 1)
      data cptrangesub/'200','300','400','500','600'/
c plot limits used with CA without filtering:
      real*8 jetmassCAminsub(nptbins)
      data jetmassCAminsub/10d0,10d0,20d0,20d0/
      real*8 jetmassCAmaxsub(nptbins)
      data jetmassCAmaxsub/190d0,190d0,290d0,290d0/
c plot limits used with CA with filtering:
      real*8 jetmassCASFminsub(nptbins)
      data jetmassCASFminsub/20d0,30d0,40d0,50d0/
      real*8 jetmassCASFmaxsub(nptbins)
      data jetmassCASFmaxsub/160d0,190d0,290d0,290d0/
c plot limits used with antikt:
      real*8 jetmassantiktminsub(nptbins)
      data jetmassantiktminsub/10d0,10d0,20d0,20d0/
      real*8 jetmassantiktmaxsub(nptbins)
      data jetmassantiktmaxsub/140d0,190d0,230d0,290d0/
c
c Ratio of 3-jet xs to 2-jet, arXiv:1304.7498, CMS
      integer nptrangeR32
      parameter (nptrangeR32 = 27)
      real*8 ptrangeR32(nptrangeR32 + 1)
      data ptrangeR32/250d0,275d0,300d0,330d0,360d0,390d0,420d0,450d0,
     >               480d0,510d0,540d0,570d0,600d0,640d0,680d0,720d0,
     >               760d0,800d0,850d0,900d0,950d0,1000d0,1060d0,
     >               1120d0,1180d0,1250d0,1320d0,1390d0/
c
      integer nptrangesveto
      parameter (nptrangesveto = 7)
      real*8 ptrangesveto(nptrangesveto + 1)
      data ptrangesveto/70d0,90d0,120d0,150d0,180d0,210d0,240d0,270d0/
      character*3 cptrangesveto(nptrangesveto + 1)
      data cptrangesveto/' 70',' 90','120','150',
     >                   '180','210','240','270'/
      integer nyrangesveto
      parameter (nyrangesveto = 5)
      real*8 yrangesveto(nyrangesveto + 1)
      data yrangesveto/1d0,2d0,3d0,4d0,5d0,6d0/
      character*1 cyrangesveto(nyrangesveto + 1)
      data cyrangesveto/'1','2','3','4','5','6'/
      integer nptbinsveto
      parameter (nptbinsveto = 18)
      real*8 ptbinsveto(nptbinsveto + 1)
      data ptbinsveto/50d0,60d0,70d0,80d0,90d0,105d0,120d0,135d0,150d0,
     >                180d0,210d0,240d0,270d0,300d0,340d0,380d0,420d0,
     >                460d0,500d0/
c
      integer nyrangesdiffxs
      parameter (nyrangesdiffxs = 5)
      real*8 yrangesdiffxs(nyrangesdiffxs + 1)
      data yrangesdiffxs/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0/
      character*3 cyrangesdiffxs(nyrangesdiffxs + 1)
      data cyrangesdiffxs/'0.0','0.5','1.0','1.5','2.0','2.5'/
c Inclusive jet 0.0 < |y| < 0.5:
      integer nptbinsdiffxs1
      parameter (nptbinsdiffxs1 = 33)
      real*8 ptbinsdiffxs1(nptbinsdiffxs1 + 1)
      data ptbinsdiffxs1/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                   272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                   507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                   846d0,905d0,967d0,1032d0,1101d0,1172d0,
     >                   1248d0,1327d0,1410d0,1497d0,1588d0,1784d0,
     >                   2116d0/
c Inclusive jet, 0.5 < |y| < 1.0:
      integer nptbinsdiffxs2
      parameter (nptbinsdiffxs2 = 30)
      real*8 ptbinsdiffxs2(nptbinsdiffxs2 + 1)
      data ptbinsdiffxs2/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                   272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                   507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                   846d0,905d0,967d0,1032d0,1101d0,1172d0,
     >                   1248d0,1327d0,1410d0,1784d0/
c Inclusive jet, 1.0 < |y| < 1.5:
      integer nptbinsdiffxs3
      parameter (nptbinsdiffxs3 = 27)
      real*8 ptbinsdiffxs3(nptbinsdiffxs3 + 1)
      data ptbinsdiffxs3/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                   272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                   507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                   846d0,905d0,967d0,1032d0,1101d0,1172d0,
     >                   1684d0/
c Inclusive jet, 1.5 < |y| < 2.0:
      integer nptbinsdiffxs4
      parameter (nptbinsdiffxs4 = 24)
      real*8 ptbinsdiffxs4(nptbinsdiffxs4 + 1)
      data ptbinsdiffxs4/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                   272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                   507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                   846d0,905d0,967d0,1248d0/
c Inclusive jet, 2.0 < |y| < 2.5:
      integer nptbinsdiffxs5
      parameter (nptbinsdiffxs5 = 19)
      real*8 ptbinsdiffxs5(nptbinsdiffxs5 + 1)
      data ptbinsdiffxs5/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                   272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                   507d0,548d0,592d0,638d0,686d0,905d0/
c Dijet mass, 0.0 < |y| < 0.5:
      integer nmbinsdiffxs6
      parameter (nmbinsdiffxs6 = 13)
      real*8 mbinsdiffxs6(nmbinsdiffxs6 + 1)
      data mbinsdiffxs6/197d0,296d0,419d0,565d0,740d0,944d0,1181d0,
     >                   1455d0,1770d0,2132d0,2546d0,3019d0,3416d0,
     >                   4010d0/
c Dijet mass, 0.5 < |y| < 1.0:
      integer nmbinsdiffxs7
      parameter (nmbinsdiffxs7 = 12)
      real*8 mbinsdiffxs7(nmbinsdiffxs7 + 1)
      data mbinsdiffxs7/270d0,386d0,526d0,693d0,890d0,1118d0,1383d0,
     >                   1687d0,2037d0,2438d0,2895d0,3416d0,4010d0/
c Dijet mass, 1.0 < |y| < 1.5:
      integer nmbinsdiffxs8
      parameter (nmbinsdiffxs8 = 11)
      real*8 mbinsdiffxs8(nmbinsdiffxs8 + 1)
      data mbinsdiffxs8/419d0,565d0,740d0,944d0,1181d0,1455d0,1770d0,
     >                   2132d0,2546d0,3147d0,3854d0,4509d0/
c Dijet mass, 1.5 < |y| < 2.0:
      integer nmbinsdiffxs9
      parameter (nmbinsdiffxs9 = 10)
      real*8 mbinsdiffxs9(nmbinsdiffxs9 + 1)
      data mbinsdiffxs9/565d0,740d0,944d0,1181d0,1455d0,1770d0,2132d0,
     >                   2546d0,3019d0,3558d0,5058d0/
c Dijet mass, 2.0 < |y| < 2.5:
      integer nmbinsdiffxs10
      parameter (nmbinsdiffxs10 = 8)
      real*8 mbinsdiffxs10(nmbinsdiffxs10 + 1)
      data mbinsdiffxs10/1000d0,1246d0,1530d0,1856d0,2231d0,2659d0,
     >                    3147d0,3704d0,5058d0/
c
c Differential dijet cross section, Atlas:
      integer nystarrangesdiffxsAtlas
      parameter (nystarrangesdiffxsAtlas = 6)
      real*8 ystarrangesdiffxsAtlas(nystarrangesdiffxsAtlas + 1)
      data ystarrangesdiffxsAtlas/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0,3.0d0/
      character*3 cystarrangesdiffxsAtlas(nystarrangesdiffxsAtlas + 1)
      data cystarrangesdiffxsAtlas/'0.0','0.5','1.0','1.5','2.0','2.5',
     >                             '3.0'/
      integer nRpardiffxsAtlas
      parameter (nRpardiffxsAtlas = 2)
      real*8 RpardiffxsAtlas(nRpardiffxsAtlas)
      data RpardiffxsAtlas/0.4d0,0.6d0/
      character*3 cRpardiffxsAtlas(nRpardiffxsAtlas)
      data cRpardiffxsAtlas/'0.4','0.6'/
c Dijet mass, ystar < 0.5:
      integer nmbinsdiffxsAtlas1
      parameter (nmbinsdiffxsAtlas1 = 21)
      real*8 mbinsdiffxsAtlas1(nmbinsdiffxsAtlas1 + 1)
      data mbinsdiffxsAtlas1/260d0,310d0,370d0,440d0,510d0,590d0,
     >                        670d0,760d0,850d0,950d0,1060d0,1180d0,
     >                        1310d0,1450d0,1600d0,1760d0,1940d0,2120d0,
     >                        2330d0,2550d0,2780d0,4270d0/
c Dijet mass, 0.5 < ystar < 1.0:
      integer nmbinsdiffxsAtlas2
      parameter (nmbinsdiffxsAtlas2 = 21)
      real*8 mbinsdiffxsAtlas2(nmbinsdiffxsAtlas2 + 1)
      data mbinsdiffxsAtlas2/310d0,370d0,440d0,510d0,590d0,670d0,760d0,
     >                        850d0,950d0,1060d0,1180d0,1310d0,1450d0,
     >                        1600d0,1760d0,1940d0,2120d0,2330d0,2550d0,
     >                        2780d0,3040d0,4270d0/
c Dijet mass, 1.0 < ystar < 1.5:
      integer nmbinsdiffxsAtlas3
      parameter (nmbinsdiffxsAtlas3 = 19)
      real*8 mbinsdiffxsAtlas3(nmbinsdiffxsAtlas3 + 1)
      data mbinsdiffxsAtlas3/510d0,590d0,670d0,760d0,850d0,950d0,
     >                        1060d0,1180d0,1310d0,1450d0,1600d0,1760d0,
     >                        1940d0,2120d0,2330d0,2550d0,2780d0,3040d0,
     >                        3310d0,4640d0/
c Dijet mass, 1.5 < ystar < 2.0:
      integer nmbinsdiffxsAtlas4
      parameter (nmbinsdiffxsAtlas4 = 17)
      real*8 mbinsdiffxsAtlas4(nmbinsdiffxsAtlas4 + 1)
      data mbinsdiffxsAtlas4/760d0,850d0,950d0,1060d0,1180d0,1310d0,
     >                        1450d0,1600d0,1760d0,1940d0,2120d0,2330d0,
     >                        2550d0,2780d0,3040d0,3310d0,3610d0,4640d0/
c Dijet mass, 2.0 < ystar < 2.5:
      integer nmbinsdiffxsAtlas5
      parameter (nmbinsdiffxsAtlas5 = 8)
      real*8 mbinsdiffxsAtlas5(nmbinsdiffxsAtlas5 + 1)
      data mbinsdiffxsAtlas5/1310d0,1450d0,1600d0,1760d0,2120d0,
     >                        2550d0,3040d0,3610d0,5040d0/
c Dijet mass, 2.5 < ystar < 3.0:
      integer nmbinsdiffxsAtlas6
      parameter (nmbinsdiffxsAtlas6 = 4)
      real*8 mbinsdiffxsAtlas6(nmbinsdiffxsAtlas6 + 1)
      data mbinsdiffxsAtlas6/2120d0,2550d0,3040d0,3930d0,5040d0/
c
c Trijet mass following CMS PAS
      integer nymaxrangesdiffxsCMSPAS
      parameter (nymaxrangesdiffxsCMSPAS = 2)
      real*8 ymaxrangesdiffxsCMSPAS(nymaxrangesdiffxsCMSPAS + 1)
      data ymaxrangesdiffxsCMSPAS/0d0,1d0,2d0/
      character*3 cymaxrangesdiffxsCMSPAS(nymaxrangesdiffxsCMSPAS + 1)
      data cymaxrangesdiffxsCMSPAS/'0.0','1.0','2.0'/
c
      integer nRpardiffxsCMSPAS
      parameter (nRpardiffxsCMSPAS = 2)
      real*8 RpardiffxsCMSPAS(nRpardiffxsCMSPAS)
      data RpardiffxsCMSPAS/0.5d0,0.7d0/
      character*3 cRpardiffxsCMSPAS(nRpardiffxsCMSPAS)
      data cRpardiffxsCMSPAS/'0.5','0.7'/
c
      integer nm123binsdiffxsCMSPAS
      parameter (nm123binsdiffxsCMSPAS = 42)
      real*8 m123binsdiffxsCMSPAS(nm123binsdiffxsCMSPAS + 1)
      data m123binsdiffxsCMSPAS/
     >     150d0,177d0,207d0,240d0,275d0,313d0,354d0,398d0,445d0,
     >     495d0,548d0,604d0,664d0,727d0,794d0,864d0,938d0,1016d0,
     >     1098d0,1184d0,1274d0,1369d0,1469d0,1573d0,1682d0,1796d0,
     >     1916d0,2041d0,2172d0,2309d0,2452d0,2602d0,2758d0,2921d0,
     >     3092d0,3270d0,3456d0,3650d0,3852d0,4063d0,4283d0,4513d0,
     >     4753d0/
c
c Inclusive jet and dijet production, ATLAS
c arXiv:1112.6297
      integer nyrangesIncjetAtlas
      parameter (nyrangesIncjetAtlas = 7)
      real*8 yrangesIncjetAtlas(nyrangesIncjetAtlas + 1)
      data yrangesIncjetAtlas/0d0,0.3d0,0.8d0,1.2d0,2.1d0,
     >                        2.8d0,3.6d0,4.4d0/
      character*3 cyrangesIncjetAtlas(nyrangesIncjetAtlas + 1)
      data cyrangesIncjetAtlas/'0.0','0.3','0.8','1.2','2.1',
     >                         '2.8','3.6','4.4'/
      integer nRparIncjetAtlas
      parameter (nRparIncjetAtlas = 2)
      real*8 RparIncjetAtlas(nRparIncjetAtlas)
      data RparIncjetAtlas/0.4d0,0.6d0/
      character*3 cRparIncjetAtlas(nRparIncjetAtlas)
      data cRparIncjetAtlas/'0.4','0.6'/
c pt bins for 0.0 < |y| < 0.3
      integer nptbinsIncjetAtlas1
      parameter (nptbinsIncjetAtlas1 = 16)
      real*8 ptbinsIncjetAtlas1(nptbinsIncJetAtlas1 + 1)
      data ptbinsIncjetAtlas1/
     >     20d0,30d0,45d0,60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >     400d0,500d0,600d0,800d0,1000d0,1200d0,1500d0/
c pt bins for 0.3 < |y| < 0.8
      integer nptbinsIncjetAtlas2
      parameter (nptbinsIncjetAtlas2 = 16)
      real*8 ptbinsIncjetAtlas2(nptbinsIncJetAtlas2 + 1)
      data ptbinsIncjetAtlas2/
     >     20d0,30d0,45d0,60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >     400d0,500d0,600d0,800d0,1000d0,1200d0,1500d0/
c pt bins for 0.8 < |y| < 1.2
      integer nptbinsIncjetAtlas3
      parameter (nptbinsIncjetAtlas3 = 16)
      real*8 ptbinsIncjetAtlas3(nptbinsIncJetAtlas3 + 1)
      data ptbinsIncjetAtlas3/
     >     20d0,30d0,45d0,60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >     400d0,500d0,600d0,800d0,1000d0,1200d0,1500d0/
c pt bins for 1.2 < |y| < 2.1
      integer nptbinsIncjetAtlas4
      parameter (nptbinsIncjetAtlas4 = 15)
      real*8 ptbinsIncjetAtlas4(nptbinsIncJetAtlas4 + 1)
      data ptbinsIncjetAtlas4/
     >     20d0,30d0,45d0,60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >     400d0,500d0,600d0,800d0,1000d0,1200d0/
c pt bins for 2.1 < |y| < 2.8
      integer nptbinsIncjetAtlas5
      parameter (nptbinsIncjetAtlas5 = 12)
      real*8 ptbinsIncjetAtlas5(nptbinsIncJetAtlas5 + 1)
      data ptbinsIncjetAtlas5/
     >     20d0,30d0,45d0,60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >     400d0,500d0,600d0/
c pt bins for 2.8 < |y| < 3.6
      integer nptbinsIncjetAtlas6
      parameter (nptbinsIncjetAtlas6 = 9)
      real*8 ptbinsIncjetAtlas6(nptbinsIncJetAtlas6 + 1)
      data ptbinsIncjetAtlas6/
     >     20d0,30d0,45d0,60d0,80d0,110d0,160d0,210d0,260d0,310d0/
c pt bins for 3.6 < |y| < 4.4
      integer nptbinsIncjetAtlas7
      parameter (nptbinsIncjetAtlas7 = 6)
      real*8 ptbinsIncjetAtlas7(nptbinsIncJetAtlas7 + 1)
      data ptbinsIncjetAtlas7/20d0,30d0,45d0,60d0,80d0,110d0,160d0/
c These declarations are needed for the dijet mass:
      integer nystarrangesm12Atlas
      parameter (nystarrangesm12Atlas = 9)
      real*8 ystarrangesm12Atlas(nystarrangesm12Atlas + 1)
      data ystarrangesm12Atlas/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0,3d0,
     >                         3.5d0,4d0,4.4d0/
      character*3 cystarrangesm12Atlas(nystarrangesm12Atlas + 1)
      data cystarrangesm12Atlas/'0.0','0.5','1.0','1.5','2.0','2.5',
     >                          '3.0','3.5','4.0','4.4'/
c m12 bins for 0.0 < ystar < 0.5:
      integer nm12binsm12Atlas1
      parameter (nm12binsm12Atlas1 = 20)
      real*8 m12binsm12Atlas1(nm12binsm12Atlas1 + 1)
      data m12binsm12Atlas1/70d0,110d0,160d0,210d0,260d0,310d0,370d0,
     >                      440d0,510d0,590d0,670d0,760d0,850d0,950d0,
     >                      1060d0,1180d0,1310d0,1450d0,1600d0,1940d0,
     >                      2780d0/
c m12 bins for 0.5 < ystar < 1.0:
      integer nm12binsm12Atlas2
      parameter (nm12binsm12Atlas2 = 20)
      real*8 m12binsm12Atlas2(nm12binsm12Atlas2 + 1)
      data m12binsm12Atlas2/110d0,160d0,210d0,260d0,310d0,370d0,
     >                      440d0,510d0,590d0,670d0,760d0,850d0,950d0,
     >                      1060d0,1180d0,1310d0,1450d0,1600d0,1760d0,
     >                      2120d0,3310d0/
c m12 bins for 1.0 < ystar < 1.5:
      integer nm12binsm12Atlas3
      parameter (nm12binsm12Atlas3 = 21)
      real*8 m12binsm12Atlas3(nm12binsm12Atlas3 + 1)
      data m12binsm12Atlas3/160d0,210d0,260d0,310d0,370d0,
     >                      440d0,510d0,590d0,670d0,760d0,850d0,950d0,
     >                      1060d0,1180d0,1310d0,1450d0,1600d0,1760d0,
     >                      1940d0,2120d0,2550d0,3610d0/
c m12 bins for 1.5 < ystar < 2.0:
      integer nm12binsm12Atlas4
      parameter (nm12binsm12Atlas4 = 20)
      real*8 m12binsm12Atlas4(nm12binsm12Atlas4 + 1)
      data m12binsm12Atlas4/260d0,310d0,370d0,
     >                      440d0,510d0,590d0,670d0,760d0,850d0,950d0,
     >                      1060d0,1180d0,1310d0,1450d0,1600d0,1760d0,
     >                      1940d0,2120d0,2330d0,2780d0,3930d0/
c m12 bins for 2.0 < ystar < 2.5:
      integer nm12binsm12Atlas5
      parameter (nm12binsm12Atlas5 = 19)
      real*8 m12binsm12Atlas5(nm12binsm12Atlas5 + 1)
      data m12binsm12Atlas5/370d0,
     >                      440d0,510d0,590d0,670d0,760d0,850d0,950d0,
     >                      1060d0,1180d0,1310d0,1450d0,1600d0,1760d0,
     >                      1940d0,2120d0,2330d0,2550d0,3040d0,4270d0/
c m12 bins for 2.5 < ystar < 3.0:
      integer nm12binsm12Atlas6
      parameter (nm12binsm12Atlas6 = 16)
      real*8 m12binsm12Atlas6(nm12binsm12Atlas6 + 1)
      data m12binsm12Atlas6/670d0,760d0,850d0,950d0,1060d0,1180d0,
     >                      1310d0,1450d0,1600d0,1760d0,1940d0,2120d0,
     >                      2330d0,2550d0,2780d0,3310d0,4640d0/
c m12 bins for 3.0 < ystar < 3.5:
      integer nm12binsm12Atlas7
      parameter (nm12binsm12Atlas7 = 12)
      real*8 m12binsm12Atlas7(nm12binsm12Atlas7 + 1)
      data m12binsm12Atlas7/1180d0,1310d0,1450d0,1600d0,1760d0,1940d0,
     >                      2120d0,2330d0,2550d0,2780d0,3040d0,3610d0,
     >                      5040d0/
c m12 bins for 3.5 < ystar < 4.0:
      integer nm12binsm12Atlas8
      parameter (nm12binsm12Atlas8 = 9)
      real*8 m12binsm12Atlas8(nm12binsm12Atlas8 + 1)
      data m12binsm12Atlas8/1760d0,1940d0,2120d0,2330d0,2550d0,2780d0,
     >                      3040d0,3310d0,3930d0,5470d0/
c m12 bins for 4.0 < ystar < 4.4:
      integer nm12binsm12Atlas9
      parameter (nm12binsm12Atlas9 = 2)
      real*8 m12binsm12Atlas9(nm12binsm12Atlas9 + 1)
      data m12binsm12Atlas9/2550d0,3040d0,4270d0/
c
c Declarations for multi-jet cross section analysis of Atlas:
      integer nptbinsMultiJetAtlas1
      parameter (nptbinsMultiJetAtlas1 = 9)
      real*8 ptbinsMultiJetAtlas1(nptbinsMultiJetAtlas1 + 1)
      data ptbinsMultiJetAtlas1/80d0,110d0,160d0,210d0,260d0,310d0,
     >                          400d0,500d0,600d0,800d0/
      integer nptbinsMultiJetAtlas2
      parameter (nptbinsMultiJetAtlas2 = 10)
      real*8 ptbinsMultiJetAtlas2(nptbinsMultiJetAtlas2 + 1)
      data ptbinsMultiJetAtlas2/60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >                          400d0,500d0,600d0,800d0/
      integer nptbinsMultiJetAtlas3
      parameter (nptbinsMultiJetAtlas3 = 8)
      real*8 ptbinsMultiJetAtlas3(nptbinsMultiJetAtlas3 + 1)
      data ptbinsMultiJetAtlas3/60d0,80d0,110d0,160d0,210d0,260d0,310d0,
     >                          400d0,500d0/
      integer nptbinsMultiJetAtlas4
      parameter (nptbinsMultiJetAtlas4 = 4)
      real*8 ptbinsMultiJetAtlas4(nptbinsMultiJetAtlas4 + 1)
      data ptbinsMultiJetAtlas4/60d0,80d0,110d0,160d0,210d0/
c Ht binsizes:
      integer nhtbinsMultiJetAtlas1
      parameter (nhtbinsMultijetAtlas1 = 9)
      real*8 htbinsMultijetAtlas1(nhtbinsMultiJetAtlas1 + 1)
      data htbinsMultiJetAtlas1/140d0,200d0,260d0,350d0,450d0,600d0,
     >                          790d0,1030d0,1300d0,1600d0/
c
      integer nhtbinsMultiJetAtlas2
      parameter (nhtbinsMultijetAtlas2 = 8)
      real*8 htbinsMultijetAtlas2(nhtbinsMultiJetAtlas2 + 1)
      data htbinsMultiJetAtlas2/200d0,260d0,350d0,450d0,600d0,790d0,
     >                          1030d0,1300d0,1600d0/
c
      integer nhtbinsMultiJetAtlas3
      parameter (nhtbinsMultijetAtlas3 = 6)
      real*8 htbinsMultijetAtlas3(nhtbinsMultiJetAtlas3 + 1)
      data htbinsMultiJetAtlas3/260d0,350d0,450d0,600d0,790d0,1030d0,
     >                          1300d0/
c Declarations for the three-to-two ratios:
      integer nRpar3to2MultiJetAtlas
      parameter (nRpar3to2MultiJetAtlas = 2)
      real*8 Rpar3to2MultiJetAtlas(nRpar3to2MultiJetAtlas)
      data Rpar3to2MultiJetAtlas/0.4d0,0.6d0/
      character*3 cRpar3to2MultiJetAtlas(nRpar3to2MultiJetAtlas)
      data cRpar3to2MultiJetAtlas/'0.4','0.6'/
c
      integer nptcut1st3to2MultiJetAtlas
      parameter (nptcut1st3to2MultiJetAtlas = 3)
      real*8 ptcut1st3to2MultiJetAtlas(nptcut1st3to2MultiJetAtlas)
      data ptcut1st3to2MultiJetAtlas/80d0,110d0,160d0/
      character*3 cptcut1st3to2MultiJetAtlas(nptcut1st3to2MultiJetAtlas)
      data cptcut1st3to2MultiJetAtlas/' 80','110','160'/
      integer nptcut2nd3to2MultiJetAtlas
      parameter (nptcut2nd3to2MultiJetAtlas = 3)
      real*8 ptcut2nd3to2MultiJetAtlas(nptcut2nd3to2MultiJetAtlas)
      data ptcut2nd3to2MultiJetAtlas/60d0,80d0,110d0/
      character*3 cptcut2nd3to2MultiJetAtlas(nptcut2nd3to2MultiJetAtlas)
      data cptcut2nd3to2MultiJetAtlas/' 60',' 80','110'/
c
      integer nptbins3to2MultiJetAtlas1
      parameter (nptbins3to2MultiJetAtlas1 = 9)
      real*8 ptbins3to2MultiJetAtlas1(nptbins3to2MultiJetAtlas1 + 1)
      data ptbins3to2MultiJetAtlas1/80d0,110d0,160d0,210d0,260d0,310d0,
     >                              400d0,500d0,600d0,800d0/
      integer nptbins3to2MultiJetAtlas2
      parameter (nptbins3to2MultiJetAtlas2 = 8)
      real*8 ptbins3to2MultiJetAtlas2(nptbins3to2MultiJetAtlas2 + 1)
      data ptbins3to2MultiJetAtlas2/110d0,160d0,210d0,260d0,310d0,400d0,
     >                              500d0,600d0,800d0/
      integer nptbins3to2MultiJetAtlas3
      parameter (nptbins3to2MultiJetAtlas3 = 7)
      real*8 ptbins3to2MultiJetAtlas3(nptbins3to2MultiJetAtlas3 + 1)
      data ptbins3to2MultiJetAtlas3/160d0,210d0,260d0,310d0,400d0,500d0,
     >                              600d0,800d0/
c
c pt sum of the first two hardest jet:
      integer nht2bins3to2MultiJetAtlas
      parameter (nht2bins3to2MultiJetAtlas = 8)
      real*8 ht2bins3to2MultiJetAtlas(nht2bins3to2MultiJetAtlas + 1)
      data ht2bins3to2MultiJetAtlas/140d0,200d0,260d0,350d0,450d0,600d0,
     >                              790d0,1030d0,1300d0/
c
c Inclusive jet and dijet differential cross section, CMS, arXiv:1212.6660
c |y| < 0.5
      integer nyrangesIncJetCMS
      parameter (nyrangesIncJetCMS = 5)
      real*8 yrangesIncJetCMS(nyrangesIncJetCMS + 1)
      data yrangesIncJetCMS/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0/
      character*3 cyrangesIncJetCMS(nyrangesIncJetCMS + 1)
      data cyrangesIncJetCMS/'0.0','0.5','1.0','1.5','2.0','2.5'/
c
      integer nptbinsIncJetCMS1
      parameter (nptbinsIncJetCMS1 = 33)
      real*8 ptbinsIncJetCMS1(nptbinsIncJetCMS1 + 1)
      data ptbinsIncJetCMS1/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                      272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                      507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                      846d0,905d0,967d0,1032d0,1101d0,1172d0,
     >                      1248d0,1327d0,1410d0,1497d0,1588d0,1784d0,
     >                      2116d0/
c 0.5 < |y| < 1.0
      integer nptbinsIncJetCMS2
      parameter (nptbinsIncJetCMS2 = 30)
      real*8 ptbinsIncJetCMS2(nptbinsIncJetCMS2 + 1)
      data ptbinsIncJetCMS2/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                      272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                      507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                      846d0,905d0,967d0,1032d0,1101d0,1172d0,
     >                      1248d0,1327d0,1410d0,1784d0/
c 1.0 < |y| < 1.5
      integer nptbinsIncJetCMS3
      parameter (nptbinsIncJetCMS3 = 27)
      real*8 ptbinsIncJetCMS3(nptbinsIncJetCMS3 + 1)
      data ptbinsIncJetCMS3/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                      272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                      507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                      846d0,905d0,967d0,1032d0,1101d0,1172d0,
     >                      1684d0/
c 1.5 < |y| < 2.0
      integer nptbinsIncJetCMS4
      parameter (nptbinsIncJetCMS4 = 24)
      real*8 ptbinsIncJetCMS4(nptbinsIncJetCMS4 + 1)
      data ptbinsIncJetCMS4/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                      272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                      507d0,548d0,592d0,638d0,686d0,737d0,790d0,
     >                      846d0,905d0,967d0,1248d0/
c 2.0 < |y| < 2.5
      integer nptbinsIncJetCMS5
      parameter (nptbinsIncJetCMS5 = 19)
      real*8 ptbinsIncJetCMS5(nptbinsIncJetCMS5 + 1)
      data ptbinsIncJetCMS5/114d0,133d0,153d0,174d0,196d0,220d0,245d0,
     >                      272d0,300d0,330d0,362d0,395d0,430d0,468d0,
     >                      507d0,548d0,592d0,638d0,686d0,905d0/
c
c |ystar| < 0.5 
      integer nm12binsCMS1
      parameter (nm12binsCMS1 = 13)
      real*8 m12binsCMS1(nm12binsCMS1 + 1)
      data m12binsCMS1/197d0,296d0,419d0,565d0,740d0,944d0,1181d0,
     >                 1455d0,1770d0,2132d0,2546d0,3019d0,3416d0,
     >                 4010d0/
c 0.5 < |ystar| < 1.0 
      integer nm12binsCMS2
      parameter (nm12binsCMS2 = 12)
      real*8 m12binsCMS2(nm12binsCMS2 + 1)
      data m12binsCMS2/270d0,386d0,526d0,693d0,890d0,1118d0,1383d0,
     >                 1687d0,2037d0,2438d0,2895d0,3416d0,4010d0/
c 1.0 < |ystar| < 1.5 
      integer nm12binsCMS3
      parameter (nm12binsCMS3 = 11)
      real*8 m12binsCMS3(nm12binsCMS3 + 1)
      data m12binsCMS3/419d0,565d0,740d0,944d0,1181d0,1455d0,1770d0,
     >                 2132d0,2546d0,3147d0,3854d0,4509d0/
c 1.5 < |ystar| < 2.0 
      integer nm12binsCMS4
      parameter (nm12binsCMS4 = 10)
      real*8 m12binsCMS4(nm12binsCMS4 + 1)
      data m12binsCMS4/565d0,740d0,944d0,1181d0,1455d0,1770d0,2132d0,
     >                 2546d0,3019d0,3558d0,5058d0/
c 2.0 < |ystar| < 2.5 
      integer nm12binsCMS5
      parameter (nm12binsCMS5 = 8)
      real*8 m12binsCMS5(nm12binsCMS5 + 1)
      data m12binsCMS5/1000d0,1246d0,1530d0,1856d0,2231d0,2659d0,3147d0,
     >                 3704d0,5058d0/
c
c Hadronic event shapes with CMS:
c arXiv:1102.0068:
      integer nptrangesEvntShapeCMS
      parameter (nptrangesEvntShapeCMS = 3)
      real*8 ptrangesEvntShapeCMS(nptrangesEvntShapeCMS + 1)
      data ptrangesEvntShapeCMS/90d0,125d0,200d0,1d10/
      character*3 cptrangesEvntShapeCMS(nptrangesEvntShapeCMS + 1)
      data cptrangesEvntShapeCMS/' 90','125','200','inf'/
c Transverse, central thrust for 90 < pt < 125:
      integer nTauPerpCbinsEvntShapeCMS1
      parameter (nTauPerpCbinsEvntShapeCMS1 = 13)
      real*8 TauPerpCbinsEvntShapeCMS1(nTauPerpCbinsEvntShapeCMS1 + 1)
      data TauPerpCbinsEvntShapeCMS1/
     >     -13.75d0,-11.15d0,-9.2d0,-7.9d0,-6.6d0,-5.95d0,-5.3d0,
     >     -4.65d0,-4.0d0,-3.35d0,-2.7d0,-2.05d0,-1.4d0,-0.75d0/
c Transverse, central thrust for 125 < pt < 200:
      integer nTauPerpCbinsEvntShapeCMS2
      parameter (nTauPerpCbinsEvntShapeCMS2 = 14)
      real*8 TauPerpCbinsEvntShapeCMS2(nTauPerpCbinsEvntShapeCMS2 + 1)
      data TauPerpCbinsEvntShapeCMS2/
     >     -13.75d0,-11.15d0,-9.2d0,-7.9d0,-7.25d0,-6.6d0,-5.95d0,
     >     -5.3d0,-4.65d0,-4.0d0,-3.35d0,-2.7d0,-2.05d0,-1.4d0,-0.75d0/
c Transverse, central thrust for 200 < pt <:
      integer nTauPerpCbinsEvntShapeCMS3
      parameter (nTauPerpCbinsEvntShapeCMS3 = 12)
      real*8 TauPerpCbinsEvntShapeCMS3(nTauPerpCbinsEvntShapeCMS3 + 1)
      data TauPerpCbinsEvntShapeCMS3/
     >     -13.75d0,-11.3125d0,-8.875d0,-8.0625d0,-7.25d0,-6.4375d0,
     >     -5.625d0,-4.8125d0,-4.0d0,-3.1875d0,-2.375d0,-1.5625d0,
     >     -0.75d0/
c Minor, central thrust for 90 < pt < 125:
      integer nTMinorCbinsEvntShapeCMS1
      parameter (nTMinorCbinsEvntShapeCMS1 = 13)
      real*8 TMinorCbinsEvntShapeCMS1(nTMinorCbinsEvntShapeCMS1 + 1)
      data TMinorCbinsEvntShapeCMS1/
     >     -6.25d0,-5.05d0,-4.15d0,-3.55d0,-2.95d0,-2.65d0,-2.35d0,
     >     -2.05d0,-1.75d0,-1.45d0,-1.15d0,-0.85d0,-0.55d0,-0.25d0/
c Minor, central thrust for 125 < pt < 200:
      integer nTMinorCbinsEvntShapeCMS2
      parameter (nTMinorCbinsEvntShapeCMS2 = 14)
      real*8 TMinorCbinsEvntShapeCMS2(nTMinorCbinsEvntShapeCMS2 + 1)
      data TMinorCbinsEvntShapeCMS2/
     >     -6.25d0,-5.05d0,-4.15d0,-3.55d0,-3.25d0,-2.95d0,-2.65d0,
     >     -2.35d0,-2.05d0,-1.75d0,-1.45d0,-1.15d0,-0.85d0,-0.55d0,
     >     -0.25d0/
c Minor, central thrust for 200 < pt:
      integer nTMinorCbinsEvntShapeCMS3
      parameter (nTMinorCbinsEvntShapeCMS3 = 12)
      real*8 TMinorCbinsEvntShapeCMS3(nTMinorCbinsEvntShapeCMS3 + 1)
      data TMinorCbinsEvntShapeCMS3/
     >     -6.25d0,-5.125d0,-4.0d0,-3.625d0,-3.25d0,-2.875d0,-2.5d0,
     >     -2.125d0,-1.75d0,-1.375d0,-1.0d0,-0.625d0,-0.25d0/
c
      call inihists
c

c vvvvvvvvvvvv This part concerns the azimuthal decorrelation vvvvvvvv
c **** arXiv:1102.2696 ****
c for normalization purposes we store the cross sections:
      do j=1,nptmaxbins
        if (j.le.5) then
          call bookupeqbins('sigma '//cptmax(j)//' < ptmax < '//
     >                         cptmax(j+1),1d0,0.5d0,1.5d0)
        else if ((j.gt.5).and.(j.le.7)) then
          call bookupeqbins('sigma '//cptmax(j)//' < ptmax < '//
     >                         cptmax(j+1),1d0,0.5d0,1.5d0)
        else if (j.eq.8) then
          call bookupeqbins('sigma '//cptmax(j)//' < ptmax < '//
     >                         cptmax(j+1),1d0,0.5d0,1.5d0)
        else if (j.eq.9) then
          call bookupeqbins('sigma '//cptmax(j)//' < ptmax',
     >                1d0,0.5d0,1.5d0)
        end if
      end do
      do j=1,nptmaxbins
        if (j.le.5) then
          call bookup('dphi '//cptmax(j)//' < ptmax < '//
     >                         cptmax(j+1),ndphibins1,dphibins1)
        else if ((j.gt.5).and.(j.le.7)) then
          call bookup('dphi '//cptmax(j)//' < ptmax < '//
     >                         cptmax(j+1),ndphibins2,dphibins2)
        else if (j.eq.8) then
          call bookup('dphi '//cptmax(j)//' < ptmax < '//
     >                         cptmax(j+1),ndphibins3,dphibins3)
        else if (j.eq.9) then
          call bookup('dphi '//cptmax(j)//' < ptmax',
     >                ndphibins4,dphibins4)
        end if
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvvvvvvvv This part concerns the jet shapes vvvvvvvvvvvvvvvvvvvvv
c **** arXiv:1101.0070 ****
c Differential jet shapes:
      do j=1,nptrangebins
        call bookupeqbins('sigma rho(r) '//cptrange(j)//' < pt < '//
     >                    cptrange(j+1),1d0,0.5d0,1.5d0)
        call bookupeqbins('rho(r) '//cptrange(j)//' < pt < '//
     >                    cptrange(j+1),deltar,0d0,Rparshapes)
      end do
      do j=1,nptrangebins
        do k=1,nyrangebins
          call bookupeqbins('sigma rho(r) '//cptrange(j)//' < pt < '//
     >                      cptrange(j+1)//' '//cyrange(k)
     >                      //' < |y| < '//cyrange(k+1),1d0,0.5d0,1.5d0)
          call bookupeqbins('rho(r) '//cptrange(j)//' < pt < '//
     >                      cptrange(j+1)//' '//cyrange(k)
     >                      //' < |y| < '//cyrange(k+1),
     >                      deltar,0d0,Rparshapes)
        end do
      end do
c Integrated jet shapes:
      call bookup('sigma 1 - psi(r=0.3)',nptrangebins,ptrangebins)
      call bookup('1 - psi(r=0.3)',nptrangebins,ptrangebins)
      do k=1,nyrangebins
        call bookup('sigma 1 - psi(r=0.3) '//cyrange(k)//' < |y| < '//
     >              cyrange(k+1),nptrangebins,ptrangebins)
        call bookup('1 - psi(r=0.3) '//cyrange(k)//' < |y| < '//
     >              cyrange(k+1),nptrangebins,ptrangebins)
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvvvvvvvv This part concerns the jet substructures vvvvvvvvvvvvvv
c **** arXiv:1203.4606 ****
c Jet masses with CA and without filtering:
      do j=1,nptbins
        call bookupeqbins('sigma mass, CA, R=1.2, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//
     >                    cptrangesub(j+1),1d0,0.5d0,1.5d0)
        call bookupeqbins('mass, CA, R=1.2, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    10d0,jetmassCAminsub(j),jetmassCAmaxsub(j))
      end do
c Jet masses with CA and with filtering:
      do j=1,nptbins
        call bookupeqbins('sigma mass, CA SF, R=1.2, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//
     >                    cptrangesub(j+1),1d0,0.5d0,1.5d0)
        call bookupeqbins('mass, CA SF, R=1.2, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    10d0,jetmassCASFminsub(j),
     >                    jetmassCASFmaxsub(j))
      end do
c Jet masses with antikt:
      do j=1,nptbins
        call bookupeqbins('sigma mass, antikt, R=1.0, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//
     >                    cptrangesub(j+1),1d0,0.5d0,1.5d0)
        call bookupeqbins('mass, antikt, R=1.0, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    10d0,jetmassantiktminsub(j),
     >                    jetmassantiktmaxsub(j))
      end do
c kt splitting scales with antikt:
      do j=1,nptbins
        call bookupeqbins('sd12, antikt, R=1.0, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    5d0,0d0,100d0)
      end do
      do j=1,nptbins
        call bookupeqbins('sd23, antikt, R=1.0, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    2d0,0d0,36d0)
      end do
c subjettiness with CA without filtering:
      do j=1,nptbins
        call bookupeqbins('tau21, CA, R=1.2, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    0.08d0,0.08d0,1.2d0)
      end do
      do j=1,nptbins
        call bookupeqbins('tau32, CA, R=1.2, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    0.08d0,0.32d0,1.12d0)
      end do
c subjettiness with antikt:
      do j=1,nptbins
        call bookupeqbins('tau21, antikt, R=1.0, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    0.08d0,0.08d0,1.2d0)
      end do
      do j=1,nptbins
        call bookupeqbins('tau32, antikt, R=1.0, |y| < 2, '//
     >                    cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                    0.08d0,0.32d0,1.12d0)
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvvv This part concerns the incl. 3-jet and 2-jet ratio, ATLAS vvvvvvvvv
c**** arXiv:1304.7498 ****
      call bookup('R32, incl. 2jet',nptrangeR32,ptrangeR32)
      call bookup('R32, incl. 3jet',nptrangeR32,ptrangeR32)
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvvv This part concerns the dijet production with a veto vvvvvvvv
c**** arXiv:1107.1641 ****
c The gap fraction as a function of rapidity difference when the two
c hardest jets are selected:
      do j=1,nptrangesveto
        call bookupeqbins('gap tot leadpt, '//cptrangesveto(j)//
     >                    ' < avg. pt < '//cptrangesveto(j+1),
     >                    0.5d0,0d0,6d0)
        call bookupeqbins('gap nojet leadpt, '//cptrangesveto(j)//
     >                    ' < avg. pt < '//cptrangesveto(j+1),
     >                    0.5d0,0d0,6d0)
      end do
c The gap fraction as a function of rapidity difference when the most
c forward/backward jets are selected:
      do j=1,nptrangesveto
        call bookupeqbins('gap tot maxdy, '//cptrangesveto(j)//
     >                    ' < avg. pt < '//cptrangesveto(j+1),
     >                    0.5d0,0d0,6d0)
        call bookupeqbins('gap nojet maxdy, '//cptrangesveto(j)//
     >                    ' < avg. pt < '//cptrangesveto(j+1),
     >                    0.5d0,0d0,6d0)
      end do
c The gap function as the average leading pt for various rapidity
c difference bins when the two leading jets are selected:
      do j=1,nyrangesveto
        call bookup('gap tot leadpt, '//cyrangesveto(j)//
     >              ' < dy < '//cyrangesveto(j+1),
     >                    nptbinsveto,ptbinsveto)
        call bookup('gap nojet leadpt, '//cyrangesveto(j)//
     >              ' < dy < '//cyrangesveto(j+1),
     >                    nptbinsveto,ptbinsveto)
c These histos are needed for the mean number of jets in the gap:
        call bookup('gap njets leadpt, '//cyrangesveto(j)//
     >              ' < dy < '//cyrangesveto(j+1),
     >                    nptbinsveto,ptbinsveto)
      end do
c The gap function as the average leading pt for various rapidity
c difference bins when the most further jets are selected:
      do j=1,nyrangesveto
        call bookup('gap tot maxdy, '//cyrangesveto(j)//
     >              ' < dy < '//cyrangesveto(j+1),
     >                    nptbinsveto,ptbinsveto)
        call bookup('gap nojet maxdy, '//cyrangesveto(j)//
     >              ' < dy < '//cyrangesveto(j+1),
     >                    nptbinsveto,ptbinsveto)
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvv This part concerns the differential jet cross sections vvvvvv
c**** arXiv:1212.6660 ****
c
c Inclusive jet, 0.0 < |y| < 0.5:
      call bookup('jetpt CMS, R=0.7, '//cyrangesdiffxs(1)//' < |y| < '//
     >            cyrangesdiffxs(2),nptbinsdiffxs1,ptbinsdiffxs1)
c Inclusive jet, 0.5 < |y| < 1.0:
      call bookup('jetpt CMS, R=0.7, '//cyrangesdiffxs(2)//' < |y| < '//
     >            cyrangesdiffxs(3),nptbinsdiffxs2,ptbinsdiffxs2)
c Inclusive jet, 1.0 < |y| < 1.5:
      call bookup('jetpt CMS, R=0.7, '//cyrangesdiffxs(3)//' < |y| < '//
     >            cyrangesdiffxs(4),nptbinsdiffxs3,ptbinsdiffxs3)
c Inclusive jet, 1.5 < |y| < 2.0:
      call bookup('jetpt CMS, R=0.7, '//cyrangesdiffxs(4)//' < |y| < '//
     >            cyrangesdiffxs(5),nptbinsdiffxs4,ptbinsdiffxs4)
c Inclusive jet, 2.0 < |y| < 2.5:
      call bookup('jetpt CMS, R=0.7, '//cyrangesdiffxs(5)//' < |y| < '//
     >            cyrangesdiffxs(6),nptbinsdiffxs5,ptbinsdiffxs5)
c
c Dijet mass, 0.0 < |y|_max < 0.5:
      call bookup('m12 CMS, R=0.7, '//cyrangesdiffxs(1)//' < |y| < '//
     >            cyrangesdiffxs(2),nmbinsdiffxs6,mbinsdiffxs6)
c Dijet mass, 0.5 < |y|_max < 1.0:
      call bookup('m12 CMS, R=0.7, '//cyrangesdiffxs(2)//' < |y| < '//
     >            cyrangesdiffxs(3),nmbinsdiffxs7,mbinsdiffxs7)
c Dijet mass, 1.0 < |y|_max < 1.5:
      call bookup('m12 CMS, R=0.7, '//cyrangesdiffxs(3)//' < |y| < '//
     >            cyrangesdiffxs(4),nmbinsdiffxs8,mbinsdiffxs8)
c Dijet mass, 1.5 < |y|_max < 2.0:
      call bookup('m12 CMS, R=0.7, '//cyrangesdiffxs(4)//' < |y| < '//
     >            cyrangesdiffxs(5),nmbinsdiffxs9,mbinsdiffxs9)
c Dijet mass, 2.0 < |y|_max < 2.5:
      call bookup('m12 CMS, R=0.7, '//cyrangesdiffxs(5)//' < |y| < '//
     >            cyrangesdiffxs(6),nmbinsdiffxs10,mbinsdiffxs10)
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvv This part concerns the differential dijet cross sections vvvvvv
c**** arXiv:1312.3524 ****
c 
      do i=1,nRpardiffxsAtlas
c Dijet mass, y* < 0.5:
        call bookup('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >              cystarrangesdiffxsAtlas(1)//' < ystar < '//
     >              cystarrangesdiffxsAtlas(2),
     >              nmbinsdiffxsAtlas1,mbinsdiffxsAtlas1)
c Dijet mass, 0.5 < y* < 1.0:
        call bookup('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >              cystarrangesdiffxsAtlas(2)//' < ystar < '//
     >              cystarrangesdiffxsAtlas(3),
     >              nmbinsdiffxsAtlas2,mbinsdiffxsAtlas2)
c Dijet mass, 1.0 < y* < 1.5:
        call bookup('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >              cystarrangesdiffxsAtlas(3)//' < ystar < '//
     >              cystarrangesdiffxsAtlas(4),
     >              nmbinsdiffxsAtlas3,mbinsdiffxsAtlas3)
c Dijet mass, 1.5 < y* < 2.0:
        call bookup('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >              cystarrangesdiffxsAtlas(4)//' < ystar < '//
     >              cystarrangesdiffxsAtlas(5),
     >              nmbinsdiffxsAtlas4,mbinsdiffxsAtlas4)
c Dijet mass, 2.0 < y* < 2.5:
        call bookup('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >              cystarrangesdiffxsAtlas(5)//' < ystar < '//
     >              cystarrangesdiffxsAtlas(6),
     >              nmbinsdiffxsAtlas5,mbinsdiffxsAtlas5)
c Dijet mass, 2.5 < y* < 3.0:
        call bookup('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >              cystarrangesdiffxsAtlas(6)//' < ystar < '//
     >              cystarrangesdiffxsAtlas(7),
     >              nmbinsdiffxsAtlas6,mbinsdiffxsAtlas6)
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvv This part concerns the trijet invariant mass distro. vvvvvv
c**** CMS PAS SMP-12-027 ****
      do i=1,nRpardiffxsCMSPAS
        do j=1,nymaxrangesdiffxsCMSPAS
          call bookup('m123, R='//cRpardiffxsCMSPAS(i)//', '//
     >                cymaxrangesdiffxsCMSPAS(j)//' < |ymax| < '//
     >                cymaxrangesdiffxsCMSPAS(j+1),
     >                nm123binsdiffxsCMSPAS,m123binsdiffxsCMSPAS)
        end do
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvv Inclusive jet and dijet production, ATLAS vvvvvv
c**** arXiv:1112.6297 ****
      do i=1,nRparIncjetAtlas
c 0.0 < |y| < 0.3:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(1)//' < |y| < '//
     >              cyrangesIncjetAtlas(2),nptbinsIncjetAtlas1,
     >              ptbinsIncjetAtlas1)
c 0.3 < |y| < 0.8:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(2)//' < |y| < '//
     >              cyrangesIncjetAtlas(3),nptbinsIncjetAtlas2,
     >              ptbinsIncjetAtlas2)
c 0.8 < |y| < 1.2:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(3)//' < |y| < '//
     >              cyrangesIncjetAtlas(4),nptbinsIncjetAtlas3,
     >              ptbinsIncjetAtlas3)
c 1.2 < |y| < 2.1:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(4)//' < |y| < '//
     >              cyrangesIncjetAtlas(5),nptbinsIncjetAtlas4,
     >              ptbinsIncjetAtlas4)
c 2.1 < |y| < 2.8:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(5)//' < |y| < '//
     >              cyrangesIncjetAtlas(6),nptbinsIncjetAtlas5,
     >              ptbinsIncjetAtlas5)
c 2.8 < |y| < 3.6:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(6)//' < |y| < '//
     >              cyrangesIncjetAtlas(7),nptbinsIncjetAtlas6,
     >              ptbinsIncjetAtlas6)
c 3.6 < |y| < 4.4:
        call bookup('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >              cyrangesIncjetAtlas(7)//' < |y| < '//
     >              cyrangesIncjetAtlas(8),nptbinsIncjetAtlas7,
     >              ptbinsIncjetAtlas7)
      end do
c dijet cross sections:
      do i=1,nRparIncjetAtlas
c 0.0 < ystar < 0.5
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(1)//' < ystar < '//
     >              cystarrangesm12Atlas(2),
     >              nm12binsm12Atlas1,m12binsm12Atlas1)
c 0.5 < ystar < 1.0
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(2)//' < ystar < '//
     >              cystarrangesm12Atlas(3),
     >              nm12binsm12Atlas2,m12binsm12Atlas2)
c 1.0 < ystar < 1.5
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(3)//' < ystar < '//
     >              cystarrangesm12Atlas(4),
     >              nm12binsm12Atlas3,m12binsm12Atlas3)
c 1.5 < ystar < 2.0
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(4)//' < ystar < '//
     >              cystarrangesm12Atlas(5),
     >              nm12binsm12Atlas4,m12binsm12Atlas4)
c 2.0 < ystar < 2.5
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(5)//' < ystar < '//
     >              cystarrangesm12Atlas(6),
     >              nm12binsm12Atlas5,m12binsm12Atlas5)
c 2.5 < ystar < 3.0
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(6)//' < ystar < '//
     >              cystarrangesm12Atlas(7),
     >              nm12binsm12Atlas6,m12binsm12Atlas6)
c 3.0 < ystar < 3.5
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(7)//' < ystar < '//
     >              cystarrangesm12Atlas(8),
     >              nm12binsm12Atlas7,m12binsm12Atlas7)
c 3.5 < ystar < 4.0
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(8)//' < ystar < '//
     >              cystarrangesm12Atlas(9),
     >              nm12binsm12Atlas8,m12binsm12Atlas8)
c 4.0 < ystar < 4.4
        call bookup('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >              cystarrangesm12Atlas(9)//' < ystar < '//
     >              cystarrangesm12Atlas(10),
     >              nm12binsm12Atlas9,m12binsm12Atlas9)
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvvvv Multi-jet cross sections with Atlas vvvvvvvvvv
c**** arXiv:1107.2092 ****
c Inclusive jet multiplicity:
      call bookupeqbins('Inclusive Jet Multiplicity, Atlas',
     >                  1d0,1.5d0,6.5d0)
c pt 1st leading jet:
      call bookup('pt 1st leading jet, Atlas',nptbinsMultiJetAtlas1,
     >            ptbinsMultiJetAtlas1)
c pt 2nd leading jet:
      call bookup('pt 2nd leading jet, Atlas',nptbinsMultiJetAtlas2,
     >            ptbinsMultiJetAtlas2)
c pt 3rd leading jet:
      call bookup('pt 3rd leading jet, Atlas',nptbinsMultiJetAtlas3,
     >            ptbinsMultiJetAtlas3)
c pt 4th leading jet:
      call bookup('pt 4th leading jet, Atlas',nptbinsMultiJetAtlas4,
     >            ptbinsMultiJetAtlas4)
c     Ht for events with njet >= 2:
      call bookup('Ht, njet >= 2, Atlas',nhtbinsMultiJetAtlas1,
     >            htbinsMultiJetAtlas1)
c     Ht for events with njet >= 3:
      call bookup('Ht, njet >= 3, Atlas',nhtbinsMultiJetAtlas2,
     >            htbinsMultiJetAtlas2)
c     Ht for events with njet >= 4:
      call bookup('Ht, njet >= 4, Atlas',nhtbinsMultiJetAtlas3,
     >            htbinsMultiJetAtlas3)
c three-to-two ratios:
      do i=1,nRpar3to2MultiJetAtlas ! R-parameter for jet recon.
        do j=2,3 ! Number of jets
          call bookup('ptj1, Atlas, njet>='//cnum(j)//
     >                ', R='//cRpar3to2MultiJetAtlas(i)//
     >                ', ptj1>'//cptcut1st3to2MultiJetAtlas(1)//
     >                ', ptj2>'//cptcut2nd3to2MultiJetAtlas(1),
     >                nptbins3to2MultiJetAtlas1,
     >                ptbins3to2MultiJetAtlas1)
          call bookup('ptj1, Atlas, njet>='//cnum(j)//
     >                ', R='//cRpar3to2MultiJetAtlas(i)//
     >                ', ptj1>'//cptcut1st3to2MultiJetAtlas(2)//
     >                ', ptj2>'//cptcut2nd3to2MultiJetAtlas(2),
     >                nptbins3to2MultiJetAtlas2,
     >                ptbins3to2MultiJetAtlas2)
          call bookup('ptj1, Atlas, njet>='//cnum(j)//
     >                ', R='//cRpar3to2MultiJetAtlas(i)//
     >                ', ptj1>'//cptcut1st3to2MultiJetAtlas(3)//
     >                ', ptj2>'//cptcut2nd3to2MultiJetAtlas(3),
     >                nptbins3to2MultiJetAtlas3,
     >                ptbins3to2MultiJetAtlas3)
        end do
      end do
      do i=1,nRpar3to2MultiJetAtlas ! R-parameter for jet recon.
        do j=2,3 ! Number of jets
          call bookup('Ht2, Atlas, njet>='//cnum(j)//
     >                ', R='//cRpar3to2MultiJetAtlas(i)//
     >                ', ptj1>'//cptcut1st3to2MultiJetAtlas(1)//
     >                ', ptj2>'//cptcut2nd3to2MultiJetAtlas(1),
     >                nHt2bins3to2MultiJetAtlas,
     >                Ht2bins3to2MultiJetAtlas)
          call bookup('Ht2, Atlas, njet>='//cnum(j)//
     >                ', R='//cRpar3to2MultiJetAtlas(i)//
     >                ', ptj1>'//cptcut1st3to2MultiJetAtlas(2)//
     >                ', ptj2>'//cptcut2nd3to2MultiJetAtlas(2),
     >                nHt2bins3to2MultiJetAtlas,
     >                Ht2bins3to2MultiJetAtlas)
          call bookup('Ht2, Atlas, njet>='//cnum(j)//
     >                ', R='//cRpar3to2MultiJetAtlas(i)//
     >                ', ptj1>'//cptcut1st3to2MultiJetAtlas(3)//
     >                ', ptj2>'//cptcut2nd3to2MultiJetAtlas(3),
     >                nHt2bins3to2MultiJetAtlas,
     >                Ht2bins3to2MultiJetAtlas)
        end do
      end do
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c vvvvvvvvvvv Hadronic event shapes with CMS vvvvvvvvvvvv
c**** arXiv:1102.0068 ****
c Corresponding cross sections:
      do i=1,nptrangesEvntShapeCMS
c The cross section should be the same in the following two
c histos, both are kept for better readibility and transparency:
        call bookupeqbins('sigma logTauPerpC CMS, '//
     >         cptrangesEvntShapeCMS(i)//' < pt < '//
     >         cptrangesEvntShapeCMS(i+1),1d0,0.5d0,1.5d0)
        call bookupeqbins('sigma logTMinorC CMS, '//
     >         cptrangesEvntShapeCMS(i)//' < pt < '//
     >         cptrangesEvntShapeCMS(i+1),1d0,0.5d0,1.5d0)
      end do
c log of transverse thrust, central, 90 < pt < 125:
      call bookup('logTauPerpC CMS, '//
     >            cptrangesEvntShapeCMS(1)//' < pt < '//
     >            cptrangesEvntShapeCMS(2),
     >            nTauPerpCbinsEvntShapeCMS1,
     >            TauPerpCbinsEvntShapeCMS1)
c log of minor thrust, central, 90 < pt < 125:
      call bookup('logTMinorC CMS, '//
     >            cptrangesEvntShapeCMS(1)//' < pt < '//
     >            cptrangesEvntShapeCMS(2),
     >            nTMinorCbinsEvntShapeCMS1,
     >            TMinorCbinsEvntShapeCMS1)
c log of transverse thrust, central, 125 < pt < 200:
      call bookup('logTauPerpC CMS, '//
     >            cptrangesEvntShapeCMS(2)//' < pt < '//
     >            cptrangesEvntShapeCMS(3),
     >            nTauPerpCbinsEvntShapeCMS2,
     >            TauPerpCbinsEvntShapeCMS2)
c log of minor thrust, central, 125 < pt < 200:
      call bookup('logTMinorC CMS, '//
     >            cptrangesEvntShapeCMS(2)//' < pt < '//
     >            cptrangesEvntShapeCMS(3),
     >            nTMinorCbinsEvntShapeCMS2,
     >            TMinorCbinsEvntShapeCMS2)
c log of transverse thrust, central, 200 < pt:
      call bookup('logTauPerpC CMS, '//
     >            cptrangesEvntShapeCMS(3)//' < pt < '//
     >            cptrangesEvntShapeCMS(4),
     >            nTauPerpCbinsEvntShapeCMS3,
     >            TauPerpCbinsEvntShapeCMS3)
c log of minor thrust, central, 200 < pt:
      call bookup('logTMinorC CMS, '//
     >            cptrangesEvntShapeCMS(3)//' < pt < '//
     >            cptrangesEvntShapeCMS(4),
     >            nTMinorCbinsEvntShapeCMS3,
     >            TMinorCbinsEvntShapeCMS3)
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
c      include 'pwhg_kn.h' 
c      include 'pwhg_flg.h'
c      include 'LesHouches.h'
      integer   maxjet,mjets,njets,numjets,ntracks
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4),ptot(4),pjnew(4,maxjet)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer j,k,i,jj,ii
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
      real * 8 dsig(7)
      integer nweights
      logical ini
      data ini/.true./
      logical inimulti
      data inimulti/.true./
      integer  minlo
c      data minlo/1/
      save inimulti,minlo,ini,minnumjets
      integer minnumjets,maxnumjets
      logical pwhg_isfinite
      external pwhg_isfinite
      character * 4 prefix
c
      real*8 MinJetPt,MaxJetRap,MaxJetPseudorap
c
      real*8 dy12,deta12,dr12,dphi12
c
      real*8 pj1(4),pj2(4),pj3(4)
      real*8 ptj1,ptj2,ptj3
      real*8 etaj1,etaj2,etaj3
      real*8 yj1,yj2,yj3
      real*8 mj1,mj2,mj3
c
      integer l
      real*8 rho,psi
c
      real*8 calcrho,calcpsi
      external calcrho,calcpsi
c
      integer ijet
      real*8 p_j(4)
      real*8 pt_j,eta_j,y_j,m_j
c
      integer ijetymax,ijetymin
      real*8 ptavg12
c
      real*8 ystar
      real*8 dystar
c
      real*8 ymax
      real*8 dymax
c
      real*8 Ht
c
      integer njets_gap
c
      real*8 tauperpC,TminorC
c
      integer nnum
      parameter (nnum = 9)
      character*1 cnum(nnum)
      data cnum/'1','2','3','4','5','6','7','8','9'/
c
c These declarations are needed for azimuthal decorrelations:
      real*8 MinJetPtdphi
      parameter (MinJetPtdphi = 100d0)
      real*8 MaxJetRapdphi
      parameter (MaxJetRapdphi = 2.8d0)
      real*8 MaxRapJ1J2dphi
      parameter (MaxRapJ1J2dphi = 0.8d0)
      real*8 Rpardphi
      parameter (Rpardphi = 0.6d0)
c
      integer nptmaxbins
      parameter (nptmaxbins = 9)
      real*8 dptmaxbins(nptmaxbins + 1)
      data dptmaxbins/110d0,160d0,210d0,260d0,310d0,400d0,500d0,600d0,
     1                800d0,10000d0/
      character*3 cptmax(9)
      data cptmax/'110','160','210','260','310','400','500',
     1            '600','800'/
c
c These declarations for the jet shapes:
      integer nptrangebins
      parameter (nptrangebins = 11)
      real*8 ptrangebins(nptrangebins + 1)
      data ptrangebins/30d0,40d0,60d0,80d0,110d0,160d0,210d0,260d0,
     1                 310d0,400d0,500d0,600d0/
      character*3 cptrange(nptrangebins + 1)
      data cptrange/' 30',' 40',' 60',' 80','110','160','210','260',
     1              '310','400','500','600'/
      integer nyrangebins
      parameter (nyrangebins = 5)
      real*8 yrangebins(nyrangebins + 1)
      data yrangebins/0.0d0,0.3d0,0.8d0,1.2d0,2.1d0,2.8d0/
      character*3 cyrange(nyrangebins + 1)
      data cyrange/'0.0','0.3','0.8','1.2','2.1','2.8'/
      real*8 deltar
      parameter (deltar = 0.1d0)
      real*8 Rparshapes
      parameter (Rparshapes = 0.6d0)
      real*8 MaxJetRapshapes
      parameter (MaxJetRapshapes = 2.8d0)
      real*8 MinJetPtshapes
      parameter (MinJetPtshapes = 30d0)
c
c These declarations for the jet subtructure and jet mass:
      integer nptbins
      parameter (nptbins = 4)
      real*8 ptbins(nptbins + 1)
      data ptbins/200d0,300d0,400d0,500d0,600d0/
      character*3 cptrangesub(nptbins + 1)
      data cptrangesub/'200','300','400','500','600'/
      real*8 jetmassmaxsub(nptbins)
      data jetmassmaxsub/220d0,220d0,310d0,310d0/
c
      real*8 RparCAsub
      parameter (RparCAsub = 1.2d0)
      real*8 RparAntiKTsub
      parameter (RparAntiKTsub = 1.0d0)
      real*8 MaxJetRapsub
      parameter (MaxJetRapsub = 2.0d0)
      real*8 MinJetPtsub
      parameter (MinJetPtsub = 200d0)
      real*8 sd12(maxjet),sd23(maxjet)
      real*8 tau21(maxjet),tau32(maxjet)
c
c These declarations for the R32 ratio:
      real*8 MinJetPtR32
      parameter (MinJetPtR32 = 150d0)
      real*8 MaxJetRapR32
      parameter (MaxJetRapR32 = 2.5d0)
      real*8 RparR32
      parameter (RparR32 = 0.7d0)
c
c These declarations for the dijet veto:
      integer nptrangesveto
      parameter (nptrangesveto = 7)
      real*8 ptrangesveto(nptrangesveto + 1)
      data ptrangesveto/70d0,90d0,120d0,150d0,180d0,210d0,240d0,270d0/
      character*3 cptrangesveto(nptrangesveto + 1)
      data cptrangesveto/' 70',' 90','120','150',
     >                   '180','210','240','270'/
      integer nyrangesveto
      parameter (nyrangesveto = 5)
      real*8 yrangesveto(nyrangesveto + 1)
      data yrangesveto/1d0,2d0,3d0,4d0,5d0,6d0/
      character*1 cyrangesveto(nyrangesveto + 1)
      data cyrangesveto/'1','2','3','4','5','6'/
c
      real*8 MinJetPtVeto
      parameter (MinJetPtVeto = 20d0)
      real*8 MaxJetRapVeto
      parameter (MaxJetRapVeto = 4.4d0)
      real*8 MinAvgPt12Veto
      parameter (MinAvgPt12Veto = 50d0)
      real*8 RparVeto
      parameter (RparVeto = 0.6d0)
      real*8 Q0Veto
      parameter (Q0Veto = 20d0)
      real*8 MaxRapDiff12Veto
      parameter (MaxRapDiff12Veto = 6d0)
c
c These declarations for the differential cross sections:
      real*8 Min1stJetPtdiffxs
      parameter (Min1stJetPtdiffxs = 60d0)
      real*8 Min2ndJetPtdiffxs
      parameter (Min2ndJetPtdiffxs = 30d0)
      real*8 MaxJetRapdiffxs
      parameter (MaxJetRapdiffxs = 2.5d0)
      real*8 Rpardiffxs
      parameter (Rpardiffxs = 0.7d0)
c
      integer nyrangesdiffxs
      parameter (nyrangesdiffxs = 5)
      real*8 yrangesdiffxs(nyrangesdiffxs + 1)
      data yrangesdiffxs/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0/
      character*3 cyrangesdiffxs(nyrangesdiffxs + 1)
      data cyrangesdiffxs/'0.0','0.5','1.0','1.5','2.0','2.5'/
c
c These declarations for the differential cross section of Atlas:
      real*8 Min1stJetPtdiffxsAtlas
      parameter (Min1stJetPtdiffxsAtlas = 100d0)
      real*8 Min2ndJetPtdiffxsAtlas
      parameter (Min2ndJetPtdiffxsAtlas = 50d0)
      real*8 MaxJetRapdiffxsAtlas
      parameter (MaxJetRapdiffxsAtlas = 3.0d0)
c
      integer nystarrangesdiffxsAtlas
      parameter (nystarrangesdiffxsAtlas = 6)
      real*8 ystarrangesdiffxsAtlas(nystarrangesdiffxsAtlas + 1)
      data ystarrangesdiffxsAtlas/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0,3.0d0/
      character*3 cystarrangesdiffxsAtlas(nystarrangesdiffxsAtlas + 1)
      data cystarrangesdiffxsAtlas/'0.0','0.5','1.0','1.5','2.0','2.5',
     >                             '3.0'/
      integer nRpardiffxsAtlas
      parameter (nRpardiffxsAtlas = 2)
      real*8 RpardiffxsAtlas(nRpardiffxsAtlas)
      data RpardiffxsAtlas/0.4d0,0.6d0/
      character*3 cRpardiffxsAtlas(nRpardiffxsAtlas)
      data cRpardiffxsAtlas/'0.4','0.6'/
c
c These declarations are needed for the three-jet invariant mass, CMS PAS:
      real*8 MinJetPtdiffxsCMSPAS
      parameter (MinJetPtdiffxsCMSPAS = 100d0)
      real*8 MaxJetRapdiffxsCMSPAS
      parameter (MaxJetRapdiffxsCMSPAS = 3.0d0)
c
      integer nRpardiffxsCMSPAS
      parameter (nRpardiffxsCMSPAS = 2)
      real*8 RpardiffxsCMSPAS(nRpardiffxsCMSPAS)
      data RpardiffxsCMSPAS/0.5d0,0.7d0/
      character*3 cRpardiffxsCMSPAS(nRpardiffxsCMSPAS)
      data cRpardiffxsCMSPAS/'0.5','0.7'/
c
      integer nymaxrangesdiffxsCMSPAS
      parameter (nymaxrangesdiffxsCMSPAS = 2)
      real*8 ymaxrangesdiffxsCMSPAS(nymaxrangesdiffxsCMSPAS + 1)
      data ymaxrangesdiffxsCMSPAS/0d0,1d0,2d0/
      character*3 cymaxrangesdiffxsCMSPAS(nymaxrangesdiffxsCMSPAS + 1)
      data cymaxrangesdiffxsCMSPAS/'0.0','1.0','2.0'/
c These declarations for inclusive jet and dijet cross sections
c at Atlas
      integer nyrangesIncjetAtlas
      parameter (nyrangesIncjetAtlas = 7)
      real*8 yrangesIncjetAtlas(nyrangesIncjetAtlas + 1)
      data yrangesIncjetAtlas/0d0,0.3d0,0.8d0,1.2d0,2.1d0,
     >                        2.8d0,3.6d0,4.4d0/
      character*3 cyrangesIncjetAtlas(nyrangesIncjetAtlas + 1)
      data cyrangesIncjetAtlas/'0.0','0.3','0.8','1.2','2.1',
     >                         '2.8','3.6','4.4'/
      integer nRparIncjetAtlas
      parameter (nRparIncjetAtlas = 2)
      real*8 RparIncjetAtlas(nRparIncjetAtlas)
      data RparIncjetAtlas/0.4d0,0.6d0/
      character*3 cRparIncjetAtlas(nRparIncjetAtlas)
      data cRparIncjetAtlas/'0.4','0.6'/
c
      real*8 MinJetPtIncjetAtlas
      parameter (MinJetPtIncjetAtlas = 20d0)
      real*8 MaxJetRapIncjetAtlas
      parameter (MaxJetRapIncjetAtlas = 4.4d0)
c
      real*8 Min1stJetPtm12Atlas
      parameter (Min1stJetPtm12Atlas = 30d0)
c
      integer nystarrangesm12Atlas
      parameter (nystarrangesm12Atlas = 9)
      real*8 ystarrangesm12Atlas(nystarrangesm12Atlas + 1)
      data ystarrangesm12Atlas/0d0,0.5d0,1d0,1.5d0,2d0,2.5d0,3d0,
     >                         3.5d0,4d0,4.4d0/
      character*3 cystarrangesm12Atlas(nystarrangesm12Atlas + 1)
      data cystarrangesm12Atlas/'0.0','0.5','1.0','1.5','2.0','2.5',
     >                          '3.0','3.5','4.0','4.4'/
c 
c Declarations for multi-jet cross section with Atlas:
      real*8 Min1stJetPtMultiJetAtlas
      parameter (Min1stJetPtMultiJetAtlas = 80d0)
      real*8 Min2ndJetPtMultiJetAtlas
      parameter (Min2ndJetPtMultiJetAtlas = 60d0)
      real*8 MaxJetRapMultiJetAtlas
      parameter (MaxJetRapMultiJetAtlas = 2.8d0)
c
      real*8 RparMultiJetAtlas
      parameter (RparMultiJetAtlas = 0.4d0)
c
      integer nRpar3to2MultiJetAtlas
      parameter (nRpar3to2MultiJetAtlas = 2)
      real*8 Rpar3to2MultiJetAtlas(nRpar3to2MultiJetAtlas)
      data Rpar3to2MultiJetAtlas/0.4d0,0.6d0/
      character*3 cRpar3to2MultiJetAtlas(nRpar3to2MultiJetAtlas)
      data cRpar3to2MultiJetAtlas/'0.4','0.6'/
c
      integer nptcut1st3to2MultiJetAtlas
      parameter (nptcut1st3to2MultiJetAtlas = 3)
      real*8 ptcut1st3to2MultiJetAtlas(nptcut1st3to2MultiJetAtlas)
      data ptcut1st3to2MultiJetAtlas/80d0,110d0,160d0/
      character*3 cptcut1st3to2MultiJetAtlas(nptcut1st3to2MultiJetAtlas)
      data cptcut1st3to2MultiJetAtlas/' 80','110','160'/
      integer nptcut2nd3to2MultiJetAtlas
      parameter (nptcut2nd3to2MultiJetAtlas = 3)
      real*8 ptcut2nd3to2MultiJetAtlas(nptcut2nd3to2MultiJetAtlas)
      data ptcut2nd3to2MultiJetAtlas/60d0,80d0,110d0/
      character*3 cptcut2nd3to2MultiJetAtlas(nptcut2nd3to2MultiJetAtlas)
      data cptcut2nd3to2MultiJetAtlas/' 60',' 80','110'/
c Hadronic event shapes with CMS:
c arXiv:1102.0068:
      integer nptrangesEvntShapeCMS
      parameter (nptrangesEvntShapeCMS = 3)
      real*8 ptrangesEvntShapeCMS(nptrangesEvntShapeCMS + 1)
      data ptrangesEvntShapeCMS/90d0,125d0,200d0,1d10/
      character*3 cptrangesEvntShapeCMS(nptrangesEvntShapeCMS + 1)
      data cptrangesEvntShapeCMS/' 90','125','200','inf'/
c
      real*8 MaxJetPseudorapEvntShapeCMS
      parameter (MaxJetPseudorapEvntShapeCMS = 5d0)
      real*8 MaxJetPseudorapCentralEvntShapeCMS
      parameter (MaxJetPseudorapCentralEvntShapeCMS = 1.3d0)
      real*8 MinJetPtEvntShapeCMS
      parameter (MinJetPtEvntShapeCMS = 30d0)
      real*8 RparEvntShapeCMS
      parameter (RparEvntShapeCMS = 0.5d0)
c


      if (.not.pwhg_isfinite(dsig0)) then
         write(*,*) '*** dsig0 is NaN in analysis subroutine ***'
         return
      endif

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

      if (ini) then
         minlo=powheginput('#minlo')                           
         if (minlo.eq.0) then
            minnumjets = 3
            minlo=0
         else
            minnumjets = 2
            minlo=1
         endif
         ini=.false.
      endif

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
      
      jetvec=0

c vvvvvvvvvvvv This part concerns the azimuthal decorrelation vvvvvvvv
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
      R = Rpardphi
      palg = -1d0
      call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     $                        MinJetPtdphi,MaxJetRapdphi,
     $                        pj,njets,jetvec)
c
      if (njets.lt.2) goto 111
c
      pj1 = pj(:,1)
      pj2 = pj(:,2)
c
      call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
      call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
      call getdydetadphidr(pj1,pj2,dy12,deta12,dphi12,dr12)
c
      dphi12 = dphi12/pi
c
c The two leading jets should be very central:
      if ((abs(yj1).gt.MaxRapJ1J2dphi).or.
     >    (abs(yj2).gt.MaxRapJ1J2dphi)) goto 111
c
      do j=1,nptmaxbins
        if ((ptj1.gt.dptmaxbins(j)).and.(ptj1.lt.dptmaxbins(j+1))) then
          if (j.ne.nptmaxbins) then
            if ((j.lt.3).and.
     >          (dphi12.lt.23d0/24d0).and.
     >          (dphi12.gt.0.5d0)
     >         ) then
              call filld('sigma '//cptmax(j)//' < ptmax < '//
     >                 cptmax(j+1),1.0d0,dsig)
              call filld('dphi '//cptmax(j)//' < ptmax < '//
     >                   cptmax(j+1),dphi12,dsig)
            end if
            if ((j.gt.2).and.
     >          (dphi12.lt.47d0/48d0).and.
     >          (dphi12.gt.0.5d0)
     >         ) then
              call filld('sigma '//cptmax(j)//' < ptmax < '//
     >                 cptmax(j+1),1.0d0,dsig)
              call filld('dphi '//cptmax(j)//' < ptmax < '//
     >                   cptmax(j+1),dphi12,dsig)
            end if
          else if (j.eq.nptmaxbins) then
            if ((dphi12.lt.47d0/48d0).and.
     >          (dphi12.gt.0.75d0)
     >         ) then
              call filld('dphi '//cptmax(j)//' < ptmax',dphi12,dsig)
              call filld('sigma '//cptmax(j)//' < ptmax',1.0d0,dsig)
            end if
          end if
        end if
      end do
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 111  continue
      jetvec=0
c vvvvvvvvvvvv This part concerns the jet shapes vvvvvvvvvvvvvvvvvvvvv
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
c
      R = Rparshapes
      palg = -1d0
      call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     $                        MinJetPtshapes,MaxJetRapshapes,
     $                        pj,njets,jetvec)
c
      if (njets.lt.1) goto 333

c Differential in pt:
      do j=1,nptrangebins
        do l=0,5
          rho = calcrho(ntracks,njets,ptrack,pj,jetvec,
     >                  ptrangebins(j),ptrangebins(j+1),
     >                  0d0,yrangebins(nyrangebins+1),
     >                  Rparshapes,deltar,deltar/2d0 + deltar*l)
          call filld('sigma rho(r) '//cptrange(j)//' < pt < '//
     >               cptrange(j+1),1d0,rho*dsig)
          call filld('rho(r) '//cptrange(j)//' < pt < '//
     >               cptrange(j+1),deltar/2d0 + deltar*l,rho*dsig)
        end do
      end do
c
c Differential in pt and y:
c Loop over the whole pt range:
      do j=1,nptrangebins
c Loop over the total rapidity range:
        do k=1,nyrangebins
          do l=0,5
            rho = calcrho(ntracks,njets,ptrack,pj,jetvec,
     >                    ptrangebins(j),ptrangebins(j+1),
     >                    yrangebins(k),yrangebins(k+1),
     >                    Rparshapes,deltar,deltar/2d0 + deltar*l)
            call filld('sigma rho(r) '//cptrange(j)//' < pt < '//
     >                 cptrange(j+1)//' '//cyrange(k)//' < |y| < '//
     >                 cyrange(k+1),1d0,rho*dsig)
            call filld('rho(r) '//cptrange(j)//' < pt < '//
     >                 cptrange(j+1)//' '//cyrange(k)//' < |y| < '//
     >                 cyrange(k+1),deltar/2d0 + deltar*l,rho*dsig)
          end do
        end do
      end do
c
c We calculate psi(r = 0.3) for different pt ranges:
c We plot 1 - psi(r = 0.3) as a function of jet pt:
      do j=1,nptrangebins
        psi = calcpsi(ntracks,njets,ptrack,pj,jetvec,
     >                ptrangebins(j),ptrangebins(j + 1),
     >                yrangebins(1),yrangebins(nyrangebins + 1),
     >                ptrangebins(j),ptrangebins(j + 1),
     >                yrangebins(1),yrangebins(nyrangebins + 1),
     >                Rparshapes,0.3d0)
        if (psi.ge.0d0) then
          call filld('sigma 1 - psi(r=0.3)',(ptrangebins(j+1) + 
     >               ptrangebins(j))/2d0,dsig)
          call filld('1 - psi(r=0.3)',(ptrangebins(j+1) + 
     >               ptrangebins(j))/2d0,(1d0 - psi)*dsig)
        end if
      end do
c we plot 1 - psi(r = 0.3) as a function of jet pt, in fixed y range:
      do j=1,nptrangebins
        do k=1,nyrangebins
          psi = calcpsi(ntracks,njets,ptrack,pj,jetvec,
     >                  ptrangebins(j),ptrangebins(j + 1),
     >                  yrangebins(k),yrangebins(k + 1),
     >                  ptrangebins(j),ptrangebins(j + 1),
     >                  yrangebins(k),yrangebins(k + 1),
     >                  Rparshapes,0.3d0)
          if (psi.ge.0d0) then
            call filld('sigma 1 - psi(r=0.3) '//
     >                 cyrange(k)//' < |y| < '//
     >                 cyrange(k+1),(ptrangebins(j+1) + 
     >                 ptrangebins(j))/2d0,dsig)
            call filld('1 - psi(r=0.3) '//cyrange(k)//' < |y| < '//
     >                 cyrange(k+1),(ptrangebins(j+1) + 
     >                 ptrangebins(j))/2d0,(1d0 - psi)*dsig)
          end if
        end do
      end do
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

333   continue

c vvvvvvvvv This part concerns the jet substructure and mass vvvvvvvvv
c
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
c First we do the analysis with CA, R=1.2:
      R = RparCAsub
      palg = 0d0
      call fastjetalg(ptrack,ntracks,R,palg,.false.,.false.,.true.,
     >                MinJetPtsub,MaxJetRapsub,
     >                pj,njets,sd12,sd23,tau21,tau32)
c The Jet mass should be plotted for various pt windows:
c We go through the jets, if it qualifies for actual ptrange
c we plot it:
      do ijet=1,njets
        p_j = pj(:,ijet)
        call getyetaptmass(p_j,y,eta,pt,m)
c loop over pt ranges:
        do j=1,nptbins
          if ((pt.gt.ptbins(j)).and.(pt.lt.ptbins(j+1))) then
            call filld('sigma mass, CA, R=1.2, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 1d0,dsig)
            call filld('mass, CA, R=1.2, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 m,dsig)
            call filld('tau21, CA, R=1.2, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 tau21(ijet),dsig)
            call filld('tau32, CA, R=1.2, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 tau32(ijet),dsig)
          end if
        end do
      end do
c We do the analysis with CA with Split and Filtering, R=1.2:
      R = RparCAsub
      palg = 0d0
      call fastjetalg(ptrack,ntracks,R,palg,.true.,.false.,.false.,
     >                MinJetPtsub,MaxJetRapsub,
     >                pj,njets,sd12,sd23,tau21,tau32)
c The Jet mass should be plotted for various pt windows:
c We go through the jets, if it qualifies for actual ptrange
c we plot it:
c      print *,"Number of jets is: ",njets
      do ijet=1,njets
        p_j = pj(:,ijet)
        call getyetaptmass(p_j,y,eta,pt,m)
c        print *,"ijet: ",ijet
c        print *,"m: ",m
c loop over pt ranges:
        do j=1,nptbins
          if ((pt.gt.ptbins(j)).and.(pt.lt.ptbins(j+1))) then
            call filld('sigma mass, CA SF, R=1.2, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 1d0,dsig)
            call filld('mass, CA SF, R=1.2, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 m,dsig)
          end if
        end do
      end do
c We do the analysis with antikt:
      R = RparAntiKTsub
      palg = -1d0
      call fastjetalg(ptrack,ntracks,R,palg,.false.,.true.,.false.,
     >                MinJetPtsub,MaxJetRapsub,
     >                pj,njets,sd12,sd23,tau21,tau32)
c The Jet mass should be plotted for various pt windows:
c We go through the jets, if it qualifies for actual ptrange
c we plot it:
c      print *,"Number of jets is: ",njets
      do ijet=1,njets
        p_j = pj(:,ijet)
        call getyetaptmass(p_j,y,eta,pt,m)
c        print *,"ijet: ",ijet
c        print *,"m: ",m
c        print *,"sd12: ",sd12(ijet)
c        print *,"sd23: ",sd23(ijet)
c loop over pt ranges:
        do j=1,nptbins
          if ((pt.gt.ptbins(j)).and.(pt.lt.ptbins(j+1))) then
            call filld('sigma mass, antikt, R=1.0, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 1d0,dsig)
            call filld('mass, antikt, R=1.0, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 m,dsig)
c we use the same loop to plot the kt splitting scales:
c for these we did not define a cross section histo, we use the one
c set up for the jet mass:
            call filld('sd12, antikt, R=1.0, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 sd12(ijet),dsig)
            call filld('sd23, antikt, R=1.0, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 sd23(ijet),dsig)
            call filld('tau21, antikt, R=1.0, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 tau21(ijet),dsig)
            call filld('tau32, antikt, R=1.0, |y| < 2, '//
     >                 cptrangesub(j)//' < pt < '//cptrangesub(j+1),
     >                 tau32(ijet),dsig)
          end if
        end do
      end do
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      jetvec=0
c vvvvvvvvvvv This part concerns the R32 ratio calculation vvvvvvvvvvv
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
      R = RparR32
      palg = -1d0
      MinJetPt = MinJetPtR32
      MaxJetRap = 1d10
      call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     $                        MinJetPt,MaxJetRap,
     $                        pj,njets,jetvec)
c
      if (njets.lt.2) goto 444
c
      pj1 = pj(:,1)
      pj2 = pj(:,2)
      if (njets.gt.2) pj3 = pj(:,3)
c
c yj3 initially put to zero:
      yj3 = 0d0
      call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
      call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
      if (njets.gt.2) call getyetaptmass(pj3,yj3,etaj3,ptj3,mj3)
c
c The first two(three) jets should be central:
      if (max(abs(yj1),abs(yj2),abs(yj3)).gt.MaxJetRapR32) goto 444
c
      ptavg12 = (ptj1 + ptj2)/2d0
c
      if (njets.ge.2) call filld('R32, incl. 2jet',ptavg12,dsig)
      if (njets.ge.3) call filld('R32, incl. 3jet',ptavg12,dsig)
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

444   continue
      jetvec=0
c vvvvvvvvvvvvvvvv This part concerns the dijet veto vvvvvvvvvvvvvvvvv
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
      R = RparVeto
      palg = -1d0
      call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     $                        MinJetPtVeto,MaxJetRapVeto,
     $                        pj,njets,jetvec)
c
      if (njets.lt.2) goto 555
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **** In this part we start with the two hardest jets: ****
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c We select the dijet system by their pt which is required to be the
c two hardest:
      pj1 = pj(:,1)
      pj2 = pj(:,2)
c
      call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
      call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
c
      ptavg12 = (ptj1 + ptj2)/2d0
      dy12 = abs(yj1 - yj2)
c
c We have an additional cut on the average pt of the dijet system
      if (ptavg12.lt.MinAvgPt12Veto) goto 566
c We have a maximal difference in rapidity too:
      if (dy12.gt.MaxRapDiff12Veto) goto 566
c To calculate the gap fraction we need the cross section for all
c jets passing the cuts:
c For different pt ranges:
      do j=1,nptrangesveto
        if ((ptavg12.gt.ptrangesveto(j)).and.
     >      (ptavg12.lt.ptrangesveto(j+1))) then
          call filld('gap tot leadpt, '//cptrangesveto(j)//
     >               ' < avg. pt < '//cptrangesveto(j+1),dy12,dsig)
        end if
      end do
c For different \Delta y ranges:
      do j=1,nyrangesveto
        if ((dy12.gt.yrangesveto(j)).and.
     >      (dy12.lt.yrangesveto(j+1))) then
          call filld('gap tot leadpt, '//cyrangesveto(j)//
     >               ' < dy < '//cyrangesveto(j+1),ptavg12,dsig)
        end if
      end do
c We count the number of jets in the gap:
      njets_gap = 0
      do j=3,njets
        p_j = pj(:,j)
        call getyetaptmass(p_j,y_j,eta_j,pt_j,m_j)
        if ((pt_j.gt.Q0Veto).and.
     >      (y_j.gt.min(yj1,yj2)).and.
     >      (y_j.lt.max(yj1,yj2))) njets_gap = njets_gap + 1
      end do
c We store the number of jets contained by the gap:
      do j=1,nyrangesveto
        if ((dy12.gt.yrangesveto(j)).and.
     >      (dy12.lt.yrangesveto(j+1))) then
          call filld('gap njets leadpt, '//cyrangesveto(j)//
     >               ' < dy < '//cyrangesveto(j+1),
     >               ptavg12,njets_gap*dsig)
        end if
      end do
c We should have no jet between the two hardest:
      do j=3,njets
        p_j = pj(:,j)
        call getyetaptmass(p_j,y_j,eta_j,pt_j,m_j)
        if ((pt_j.gt.Q0Veto).and.
     >      (y_j.gt.min(yj1,yj2)).and.
     >      (y_j.lt.max(yj1,yj2))) goto 566
      end do
c We go through all pt ranges and we try to fill in the adequate one:
      do j=1,nptrangesveto
        if ((ptavg12.gt.ptrangesveto(j)).and.
     >      (ptavg12.lt.ptrangesveto(j+1))) then
          call filld('gap nojet leadpt, '//cptrangesveto(j)//
     >               ' < avg. pt < '//cptrangesveto(j+1),dy12,dsig)
        end if
      end do
c We do the same but with \Delta y ranges:
      do j=1,nyrangesveto
        if ((dy12.gt.yrangesveto(j)).and.
     >      (dy12.lt.yrangesveto(j+1))) then
c We have to divide by the width of the rapidity window to obtain the
c correct differential cross section:
          call filld('gap nojet leadpt, '//cyrangesveto(j)//
     >               ' < dy < '//cyrangesveto(j+1),ptavg12,dsig)
        end if
      end do
c
566   continue
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **** In this part we start with the two furthest jets: ****
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c We go through the jets and select the two most further ones:
c j1 should have the most positive rapidity while j2 the opposite:
      yj1 = -1d5
      yj2 = 1d5
      do j=1,njets
        p_j = pj(:,j)
        call getyetaptmass(p_j,y_j,eta_j,pt_j,m_j)
        if (y_j.gt.yj1) then
          yj1 = y_j
          ijetymax = j
        end if
        if (y_j.lt.yj2) then
          yj2 = y_j
          ijetymin = j
        end if
      end do

      pj1 = pj(:,ijetymax)
      pj2 = pj(:,ijetymin)
c
      call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
      call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
c
      ptavg12 = (ptj1 + ptj2)/2d0
      dy12 = abs(yj1 - yj2)
c
c We have an additional cut on the average pt of the dijet system
      if (ptavg12.lt.MinAvgPt12Veto) goto 577
c We have a maximal difference in rapidity too:
      if (dy12.gt.MaxRapDiff12Veto) goto 577
c To calculate the gap fraction we need the cross section for all
c jets passing the cuts:
      do j=1,nptrangesveto
        if ((ptavg12.gt.ptrangesveto(j)).and.
     >      (ptavg12.lt.ptrangesveto(j+1))) then
          call filld('gap tot maxdy, '//cptrangesveto(j)//
     >               ' < avg. pt < '//cptrangesveto(j+1),dy12,dsig)
        end if
      end do
c For different \Delta y ranges:
      do j=1,nyrangesveto
        if ((dy12.gt.yrangesveto(j)).and.
     >      (dy12.lt.yrangesveto(j+1))) then
          call filld('gap tot maxdy, '//cyrangesveto(j)//
     >               ' < dy < '//cyrangesveto(j+1),ptavg12,dsig)
        end if
      end do
c We should have no jet between the two hardest:
      do j=1,njets
c We go through all the jets and skip the furthest ones:
        if ((j.eq.ijetymax).or.(j.eq.ijetymin)) cycle
        p_j = pj(:,j)
        call getyetaptmass(p_j,y_j,eta_j,pt_j,m_j)
        if ((pt_j.gt.Q0Veto).and.
     >      (y_j.gt.min(yj1,yj2)).and.
     >      (y_j.lt.max(yj1,yj2))) goto 577
      end do
c We go through all pt ranges and we try to fill in the adequate one:
      do j=1,nptrangesveto
        if ((ptavg12.gt.ptrangesveto(j)).and.
     >      (ptavg12.lt.ptrangesveto(j+1))) then
          call filld('gap nojet maxdy, '//cptrangesveto(j)//
     >               ' < avg. pt < '//cptrangesveto(j+1),dy12,dsig)
        end if
      end do
c We do the same but with \Delta y ranges:
      do j=1,nyrangesveto
        if ((dy12.gt.yrangesveto(j)).and.
     >      (dy12.lt.yrangesveto(j+1))) then
c We have to divide by the width of the rapidity window to obtain the
c correct differential cross section:
          call filld('gap nojet maxdy, '//cyrangesveto(j)//
     >               ' < dy < '//cyrangesveto(j+1),ptavg12,dsig)
        end if
      end do
c
577   continue
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

555   continue
      jetvec=0
c vvvvvvvv This part concerns the differential cross sections vvvvvvvv
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
      R = Rpardiffxs
      palg = -1d0
      call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     $                        Min2ndJetPtdiffxs,MaxJetRapdiffxs,
     $                        pj,njets,jetvec)
c
c For the inclusive jet study we need at least one reconstructed jet:
      if (njets.lt.1) goto 666
c
c We go through all the jets and plot their pt's:
      do j=1,njets
        pj1   = pj(:,j)
        call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
        do k=1,nyrangesdiffxs
          if ((abs(yj1).gt.yrangesdiffxs(k)).and.
     >        (abs(yj1).lt.yrangesdiffxs(k+1))) then
            call filld('jetpt CMS, R=0.7, '//cyrangesdiffxs(k)//
     >                 ' < |y| < '//cyrangesdiffxs(k+1),ptj1,dsig)
          end if
        end do
      end do
c For dijet mass we need at least two jets to be present:
      if (njets.lt.2) goto 666
c We only interested in the first two jets:
      pj1   = pj(:,1)
      pj2   = pj(:,2)
      call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
      call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
c The first jet has to have a larger transverse momentum:
      if (ptj1.lt.Min1stJetPtdiffxs) goto 666
c
      call pwhg_getinvmass(pj1 + pj2,m)
c We go through all the y ranges and plot the mass where it belongs:
        do k=1,nyrangesdiffxs
          if ((max(abs(yj1),abs(yj2)).gt.yrangesdiffxs(k)).and.
     >        (max(abs(yj1),abs(yj2)).lt.yrangesdiffxs(k+1))) then
            call filld('m12 CMS, R=0.7, '//cyrangesdiffxs(k)//
     >                 ' < |y| < '//cyrangesdiffxs(k+1),m,dsig)
          end if
        end do
c 
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

666   continue
      jetvec=0
c vvvv This part concerns the differential cross sections at Atlas vvvv
c We can only do this analysis if minlo is activated:
      if (minlo.eq.1) then !<minlo>
      do i=1,nRpardiffxsAtlas !<sweep through R's>
        palg = -1d0
        R = RpardiffxsAtlas(i)
        MinJetPt = Min2ndJetPtdiffxsAtlas
        MaxJetRap = MaxJetRapdiffxsAtlas
        call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     $                          MinJetPt,MaxJetRap,
     $                          pj,njets,jetvec)
c
c For dijet mass we need at least two jets to be present:
        if (njets.lt.2) cycle
c We only interested in the first two jets:
        pj1   = pj(:,1)
        pj2   = pj(:,2)
        call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
        call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
c The first jet has to have a larger transverse momentum:
        if (ptj1.lt.Min1stJetPtdiffxsAtlas) cycle
c
        call pwhg_getinvmass(pj1 + pj2,m)
c We go through all the ystar ranges and plot the mass where it belongs:
        ystar = abs(yj1 - yj2)/2d0
        do k=1,nystarrangesdiffxsAtlas
          if ((ystar.gt.ystarrangesdiffxsAtlas(k)).and.
     >        (ystar.lt.ystarrangesdiffxsAtlas(k+1))) then
c We calculate the width of the ystar bin:
            dystar = ystarrangesdiffxsAtlas(k+1)
     >             - ystarrangesdiffxsAtlas(k)
            call filld('m12 Atlas2013, R='//cRpardiffxsAtlas(i)//', '//
     >                 cystarrangesdiffxsAtlas(k)//' < ystar < '//
     >                 cystarrangesdiffxsAtlas(k+1),m,dsig/dystar)
          end if
        end do
c 
      end do !<sweep through R's>
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      jetvec=0
c vvvv This part concerns the trijet invariant mass, CMS PAS vvvv
      do i=1,nRpardiffxsCMSPAS !<sweep through R's>
        jetvec=0
        palg = -1d0
        R = RpardiffxsCMSPAS(i)
        call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     >                          MinJetPtdiffxsCMSPAS,
     >                          MaxJetRapdiffxsCMSPAS,
     >                          pj,njets,jetvec)
c
c For trijet mass we need at least three jets to be present:
        if (njets.lt.3) cycle
c We only interested in the first two jets:
        pj1   = pj(:,1)
        pj2   = pj(:,2)
        pj3   = pj(:,3)
c
        call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
        call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
        call getyetaptmass(pj3,yj3,etaj3,ptj3,mj3)
c
        call pwhg_getinvmass(pj1 + pj2 + pj3,m)
c We go through all the ymax ranges and plot the mass where it belongs:
        ymax = sign(max(abs(yj1),abs(yj2),abs(yj3)),
     >              abs(max(yj1,yj2,yj3)) - abs(min(yj1,yj2,yj3)))
        do k=1,nymaxrangesdiffxsCMSPAS
          if ((abs(ymax).gt.ymaxrangesdiffxsCMSPAS(k)).and.
     >        (abs(ymax).lt.ymaxrangesdiffxsCMSPAS(k+1))) then
c We calculate the width of the ystar bin:
            dymax = 2d0*(ymaxrangesdiffxsCMSPAS(k+1)
     >             - ymaxrangesdiffxsCMSPAS(k))
            call filld('m123, R='//cRpardiffxsCMSPAS(i)//', '//
     >                 cymaxrangesdiffxsCMSPAS(k)//' < |ymax| < '//
     >                 cymaxrangesdiffxsCMSPAS(k+1),m,dsig/dymax)
          end if
        end do
c 
      end do !<sweep through R's>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cvvvvvvv Inclusive jet and dijet xs at Atlas vvvvvvvvv
      if (minlo.eq.1) then !<minlo>
      jetvec = 0
      do i=1,nRparIncjetAtlas
        R = RparIncjetAtlas(i)
        palg = -1d0
        MinJetPt = MinJetPtIncjetAtlas
        MaxJetRap = MaxJetRapIncjetAtlas
        jetvec = 0
        call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     >                          MinJetPt,MaxJetRap,
     >                          pj,njets,jetvec)
c For the inclusive jet we need at least one jet:
        if (njets.lt.1) cycle
        do ijet=1,njets
          p_j = pj(:,ijet)
          call getyetaptmass(p_j,y,eta,pt,m)
          do j=1,nyrangesIncjetAtlas
c     use abs(y)
            y = abs(y)
            if ((y.gt.yrangesIncjetAtlas(j)).and.
     >          (y.lt.yrangesIncjetAtlas(j+1))) then
              dy = yrangesIncjetAtlas(j+1) - yrangesIncjetAtlas(j)
              call filld('jetpt Atlas, R='//cRparIncjetAtlas(i)//', '//
     >                   cyrangesIncjetAtlas(j)//' < |y| < '//
     >                   cyrangesIncjetAtlas(j+1),pt,dsig/dy/2)
            end if
          end do
        end do
c For the dijet cross section we need at least two jets:
        if (njets.gt.1) then
          pj1   = pj(:,1)
          pj2   = pj(:,2)
c
          call getyetaptmass(pj1,yj1,etaj1,ptj1,mj1)
          call getyetaptmass(pj2,yj2,etaj2,ptj2,mj2)
c
c For the dijet cross section we have a higher minimal pt for the
c hardest jet: 
          if (ptj1.lt.Min1stJetPtm12Atlas) cycle
          call pwhg_getinvmass(pj1+pj2,m)
          ystar = abs(yj1 - yj2)/2d0
          do j=1,nystarrangesm12Atlas
            if ((ystar.gt.ystarrangesm12Atlas(j)).and.
     >          (ystar.lt.ystarrangesm12Atlas(j+1))) then
              dystar = (ystarrangesm12Atlas(j+1) 
     >               -  ystarrangesm12Atlas(j))
              call filld('m12 Atlas, R='//cRparIncJetAtlas(i)//', '//
     >                   cystarrangesm12Atlas(j)//' < ystar < '//
     >                   cystarrangesm12Atlas(j+1),m,dsig/dystar)
            end if
          end do
        end if
      end do
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



cvvvvvvvvvv Multi-Jet inclusive cross section Atlas. arXiv:1107.2092 vvvvvvvvvvvvvv
      if (minlo.eq.1) then !<minlo>
      do i=1,nRpar3to2MultiJetAtlas
        jetvec = 0
        R = Rpar3to2MultiJetAtlas(i)
        palg = -1d0
        MinJetPt = Min2ndJetPtMultiJetAtlas
        MaxJetRap = MaxJetRapMultiJetAtlas
        jetvec = 0
        call fastjetppgenkt_pty(ptrack,ntracks,R,palg,
     >                          MinJetPt,MaxJetRap,
     >                          pj,njets,jetvec)
        do j=2,3 ! Least number of jets
c Only the number of different pt cuts is important for the 
c do loop:
          do k=1,nptcut2nd3to2MultiJetAtlas ! pt cuts
c We apply a new pt cut on the jets:
            MinJetPt = ptcut2nd3to2MultiJetAtlas(k)
            call posterjetcuts(njets,pj,MinJetPt,mjets,pjnew)
c We need at least j jets:
            if (mjets.lt.j) cycle
c For the first jet a higher cut is needed, ptj1 is calculated:
            ptj1 = sqrt(pjnew(1,1)**2 + pjnew(2,1)**2)
c for later purposes we calculate the pt of the second hardest jet:
            ptj2 = sqrt(pjnew(1,2)**2 + pjnew(2,2)**2)
c For the first jet we have a higher minimal pt:
            if (ptj1.lt.ptcut1st3to2MultiJetAtlas(k)) cycle
c We bin the multiplicity, pt's and Ht if the conditions are
c least severe:
c For these observables we need:
c R = 0.4
c ptj1 > 80 GeV, ptj2 > 60 GeV
c mjets >= 2
            if ((i.eq.1).and.(j.eq.2).and.(k.eq.1)) then
c binning the jet multiplicity:
              call filld('Inclusive Jet Multiplicity, Atlas',
     >                   dble(mjets),dsig)
c We calcute the Ht too:
              Ht = 0d0
              do ijet=1,mjets
                p_j = pj(:,ijet)
                call getyetaptmass(p_j,y,eta,pt,m)
c Ht is the sum of jet pts:
                Ht = Ht + pt
c binning the jet pts
                if (ijet.eq.1) call filld('pt 1st leading jet, Atlas',
     >                                    pt,dsig)
                if (ijet.eq.2) call filld('pt 2nd leading jet, Atlas',
     >                                    pt,dsig)
                if (ijet.eq.3) call filld('pt 3rd leading jet, Atlas',
     >                                    pt,dsig)
                if (ijet.eq.4) call filld('pt 4th leading jet, Atlas',
     >                                    pt,dsig)
              end do
c Binning Ht
              if (mjets.ge.2) call filld('Ht, njet >= 2, Atlas',Ht,dsig)
              if (mjets.ge.3) call filld('Ht, njet >= 3, Atlas',Ht,dsig)
              if (mjets.ge.4) call filld('Ht, njet >= 4, Atlas',Ht,dsig)
            endif
c We bin the leading jet pt for the 3-to-2 ratio
            call filld('ptj1, Atlas, njet>='//cnum(j)//
     >                 ', R='//cRpar3to2MultiJetAtlas(i)//
     >                 ', ptj1>'//cptcut1st3to2MultiJetAtlas(k)//
     >                 ', ptj2>'//cptcut2nd3to2MultiJetAtlas(k),
     >                  ptj1,dsig)
c We also bin the sum of the pt of the first two leading jets:
            call filld('Ht2, Atlas, njet>='//cnum(j)//
     >                 ', R='//cRpar3to2MultiJetAtlas(i)//
     >                 ', ptj1>'//cptcut1st3to2MultiJetAtlas(k)//
     >                 ', ptj2>'//cptcut2nd3to2MultiJetAtlas(k),
     >                  ptj1 + ptj2,dsig)
          end do
        end do
      end do
c
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
777   continue 

cvvvvvvvvvvvvvvvvvvvvvv Event shapes, CMS  arXiv:1102.0068  vvvvvvvvvvvvvvvvvvvvvvvvvvvv
      if (minlo.eq.1) then !<minlo>
        R = RparEvntShapeCMS
        palg = -1d0
        MinJetPt = MinJetPtEvntShapeCMS
        MaxJetPseudorap = MaxJetPseudorapEvntShapeCMS
        jetvec = 0
        call fastjetppgenkt_pteta(ptrack,ntracks,R,palg,
     >                            MinJetPt,MaxJetPseudorap,
     >                            pj,njets,jetvec)
c The first jet should have a higher pt:
         p_j = pj(:,1)
         call getyetaptmass(p_j,y,eta,pt,m)
         if (pt.lt.ptrangesEvntShapeCMS(1)) goto 888
         if (abs(eta).gt.MaxJetPseudorapCentralEvntShapeCMS) goto 888
c We go through all the jets and consider only those which are central
c enough:
         mjets = 0
         do i=1,njets
           p_j = pj(:,i)
           call getyetaptmass(p_j,y,eta,pt,m)
           if (abs(eta).gt.MaxJetPseudorapCentralEvntShapeCMS) cycle
c The jet satisfies the central criteria, hence added to the list:
           mjets = mjets + 1
           pjnew(:,mjets) = pj(:,i)
         end do
c We need at least two jests in the central region:
         if (mjets.lt.2) goto 888
c We calculate the transverse thrust and thrust minor for the selected,
c central jets:
         call CalcThrustPerp(mjets,pjnew,tauperpC,TminorC)
c The first jet should be the same so the original can be used to 
c trigger the binning:
         p_j = pj(:,1)
         call getyetaptmass(p_j,y,eta,pt,m)
         do i=1,nptrangesEvntShapeCMS
           if ((pt.gt.ptrangesEvntShapeCMS(i)).and.
     >         (pt.lt.ptrangesEvntShapeCMS(i+1))) then
             call filld('sigma logTauPerpC CMS, '//
     >                  cptrangesEvntShapeCMS(i)//' < pt < '//
     >                  cptrangesEvntShapeCMS(i+1),1d0,dsig)
             call filld('sigma logTMinorC CMS, '//
     >                  cptrangesEvntShapeCMS(i)//' < pt < '//
     >                  cptrangesEvntShapeCMS(i+1),1d0,dsig)
             call filld('logTauPerpC CMS, '//
     >                  cptrangesEvntShapeCMS(i)//' < pt < '//
     >                  cptrangesEvntShapeCMS(i+1),log(tauperpC),dsig)
             call filld('logTMinorC CMS, '//
     >                  cptrangesEvntShapeCMS(i)//' < pt < '//
     >                  cptrangesEvntShapeCMS(i+1),log(TminorC),dsig)
           end if
         end do
      end if !<minlo>
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
888   continue 
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

      function calcrho(ntrack,njet,ptrack,pjet,jetvec,
     >                 ptmin,ptmax,ymin,ymax,Rpar,dr,r)
      implicit none
c
      real*8 calcrho
c
      integer maxtrack,maxjet
      parameter (maxtrack=2048)
      parameter (maxjet=2048)
c
      integer ntrack,njet
      real*8 ptrack(4,maxtrack),pjet(4,maxjet)
      integer jetvec(maxtrack)
      real*8 ptmin,ptmax
      real*8 ymin,ymax
      real*8 Rpar,dr,r
c
      integer ijet
      integer numjets
      real*8 p_j(4)
      real*8 pt_j,y_j
      real*8 pt_dr,pt_r
c
      real*8 calcptdr
      external calcptdr
c
      calcrho = 0d0
      numjets = 0
c
c we loop over all jets:
      do ijet=1,njet
        p_j = pjet(:,ijet)
        call pwhg_getpt(p_j,pt_j)
        call pwhg_getrapidity(p_j,y_j)
c We only consider those jets for which the pt is in a given range:
        if ((pt_j.lt.ptmin).or.(pt_j.gt.ptmax)) cycle
c We even include the possibility to make cuts on rapidity:
        if ((abs(y_j).lt.ymin).or.(abs(y_j).gt.ymax)) cycle
c We have an extra cut on the jet pt, hence njet is not the number of
c jets considered for rho:
        numjets = numjets + 1
c we calculate the pt contribution of an annulus between r - \Delta r/2
c and r + \Delta r/2
        pt_dr = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,r,dr)
c We calculate the summed pt of all the clusters in the jet:
        pt_r = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,
     >                  Rpar/2d0,Rpar)
C        calcrho = calcrho + pt_dr/pt_j
        calcrho = calcrho + pt_dr/pt_r
      end do
c normalization:
      if (numjets.gt.0) calcrho = calcrho/numjets/dr
c
      end function calcrho
c
c This function calculates the pt contribution coming from an annulus
c between r - dr/2 and r + dr/2. The jet orientation is
c given by the jet momentum (pj). 
      function calcptdr(ijet,pj,ntrack,ptrack,jetvec,r,dr)
      implicit none
c
      real*8 calcptdr
c
      integer maxtrack
      parameter (maxtrack = 2048)
c
      integer ijet
      real*8 pj(4)
      integer ntrack
      real*8 ptrack(4,maxtrack)
      integer jetvec(maxtrack)
      real*8 r,dr
c
      integer itrack
      real*8 dr2
      real*8 deltar,pt_tmp
c
c
      calcptdr = 0d0
c
      dr2 = dr/2d0
c
c We loop over all tracks and select only those which are part
c of the ijet-th jet:
      do itrack=1,ntrack
        if (jetvec(itrack).ne.ijet) cycle
c If we are here the itrack-th track is part of the ijet-th jet:
        call pwhg_getR_phiy(pj,ptrack(:,itrack),deltar)
c is it in the annulus?
        if ((deltar.lt.(r - dr2)).or.(deltar.gt.(r + dr2))) cycle
        call pwhg_getpt(ptrack(:,itrack),pt_tmp)
        calcptdr = calcptdr + pt_tmp
      end do
c
      end function calcptdr
c
c This routine calculates the integrated jet shape. 
c ptmin, ptmax are the cuts on the pt of the considered jets,
c ymin and ymax have the same role but with the rapidities,
c while ptw,yw variables are used to define the pt and y
c window.
      function calcpsi(ntrack,njet,ptrack,pjet,jetvec,
     >                 ptmin,ptmax,ymin,ymax,
     >                 ptwmin,ptwmax,ywmin,ywmax,
     >                 Rpar,r)
      implicit none
c
      real*8 calcpsi
c
      integer maxtrack,maxjet
      parameter (maxtrack=2048)
      parameter (maxjet=2048)
c
      integer ntrack,njet
      real*8 ptrack(4,maxtrack),pjet(4,maxjet)
      integer jetvec(maxtrack)
      real*8 ptmin,ptmax
      real*8 ymin,ymax
      real*8 ptwmin,ptwmax
      real*8 ywmin,ywmax
      real*8 Rpar,r
c
      integer ijet
      integer numjets
      real*8 p_j(4)
      real*8 pt_j,y_j
      real*8 pt_r
c
      real*8 calcptdr
      external calcptdr
c
c We set psi to an unphysical value:
      calcpsi = -1d5
      numjets = 0
c
c we loop over all jets:
      do ijet=1,njet
        p_j = pjet(:,ijet)
        call pwhg_getpt(p_j,pt_j)
        call pwhg_getrapidity(p_j,y_j)
c We only consider those jets for which the pt is in a given range:
        if ((pt_j.lt.ptmin).or.(pt_j.gt.ptmax)) cycle
c We even include the possibility to make cuts on rapidity:
        if ((abs(y_j).lt.ymin).or.(abs(y_j).gt.ymax)) cycle
        numjets = numjets + 1
c If we have at least one jet the psi is initialized to be physical:
        if (numjets.eq.1) calcpsi = 0d0
        if ((pt_j.lt.ptwmin).or.(pt_j.gt.ptwmax)) cycle
        if ((abs(y_j).lt.ywmin).or.(abs(y_j).gt.ywmax)) cycle
c The next few lines are needed for the window:
c we calculate the pt contribution of between 0 and r:
c we can use the previously implemented calcptdr function 
c with parameters: r = r/2 and dr = r, with this
c parameter choice the function calculates the pt 
c contribution from tracks situating in an annulus between 
c zero and r
        pt_r = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,r/2d0,r)
        pt_j = calcptdr(ijet,p_j,ntrack,ptrack,jetvec,Rpar/2d0,Rpar)
        calcpsi = calcpsi + pt_r/pt_j
      end do
c normalization:
      if (numjets.gt.0) calcpsi = calcpsi/numjets
c
      end function calcpsi
c
c This routine applies a pt cut on jet momenta a posteriori:
      subroutine posterjetcuts(njets,pjet,minpt,njetsout,pjetout)
      implicit none
c
      integer njets,njetsout
      real*8 pjet(4,njets),pjetout(4,njets)
      real*8 minpt
c
      integer ijet
      real*8 pt_tmp
c
c
      njetsout = 0
      do ijet=1,njets
        pt_tmp = sqrt(pjet(1,ijet)**2 + pjet(2,ijet)**2)
        if (pt_tmp.gt.minpt) then
          njetsout = njetsout + 1
          pjetout(:,njetsout) = pjet(:,ijet)
        end if
      end do
c
      end subroutine posterjetcuts
c
c This routine calculates the Thrust for a set of jets:
c This routine uses the algorithm published by:
c Yamamoto, J. Comput. Phys. 52, 597-601, (1983)
      subroutine CalcThrust(njet,pjet,Thrust)
      implicit none
c
      integer maxjet
      parameter (maxjet = 2048)
c
      integer njet
      real*8 pjet(4,maxjet)
      real*8 Thrust
c
      integer i,j,k
      integer si,sj,sk
      real*8 pkpixpj,maxabssp,abssp,sumabsmom
      real*8 nvec(3),pixpj(3),sumskpk(3),sumsp(3),maxsp(3)
c
c
      maxabssp = 0d0
      nvec = 0d0
      do i=1,njet-1
        do j=i+1,njet
          pixpj(1) = pjet(2,i)*pjet(3,j) - pjet(3,i)*pjet(2,j)
          pixpj(2) = pjet(3,i)*pjet(1,j) - pjet(1,i)*pjet(3,j)
          pixpj(3) = pjet(1,i)*pjet(2,j) - pjet(2,i)*pjet(1,j)
c For each (i,j) pair we calculate the sum of s_k*p_k: 
          sumskpk = 0d0
          do k=1,njet
            if ((k.eq.i).or.(k.eq.j)) cycle
            pkpixpj = pjet(1,k)*pixpj(1) + pjet(2,k)*pixpj(2)
     >              + pjet(3,k)*pixpj(3)
c This fixes sk:
            sk = sign(1d0,pkpixpj)
            sumskpk = sumskpk + sk*pjet(1:3,k)
          end do
          do si=-1,1,2
            do sj=-1,1,2
c We calculate the total sum of products: \sum_{l=1}^{n}s_l p_l for
c each (si,sj) pair:
              sumsp = sumskpk + si*pjet(1:3,i) + sj*pjet(1:3,j)
              abssp = sqrt(sumsp(1)**2 + sumsp(2)**2 + sumsp(3)**2)
              if (maxabssp.lt.abssp) then
                maxabssp = abssp
                nvec = sumsp/abssp
              end if
            end do
          end do
        end do
      end do
c
      sumabsmom = 0d0
      do i=1,njet
        sumabsmom = sumabsmom + sqrt(pjet(1,i)**2 + pjet(2,i)**2
     >            + pjet(3,i)**2)
      end do
c
      Thrust = maxabssp/sumabsmom
c
      end subroutine CalcThrust
c
c This routine calculates the transverse thrust, and 
c and heavily relies upon the previous routine:
      subroutine CalcThrustperp(njet,pjet,Tau,Tminor)
      implicit none
c
      integer maxjet
      parameter (maxjet = 2048)
c
      integer njet
      real*8 pjet(4,maxjet)
      real*8 Tau,Tminor
c
      integer i,j
      integer si,sj
      real*8 Tperp
      real*8 ptiptj,abssumspt,abssumspt_max,ptixnT,sumabspt
      real*8 sumsjptj(2),sumspt(2),sumspt_max(2),nvec(2)
c
c
      nvec = 0d0
      abssumspt_max = 0d0
      sumspt_max = 0d0
c      print *,"number of jets: ",njet
      do i=1,njet
c This will contain the \sum_{j\ne i} s_j p_{\perp\,,j} sum:
        sumsjptj = 0d0
        do j=1,njet
c When the two jets coincide simply cycle:
          if (i.eq.j) cycle
c The dotproduct of ptiptj fixes the sign for momentum pj: 
          ptiptj = pjet(1,i)*pjet(1,j) + pjet(2,i)*pjet(2,j)
c          print *,"pti.ptj: ",ptiptj
          sj = sign(1d0,ptiptj)
c          print *,"sj: ",sj
c The sum of the product of transverse momenta with appropriate signs:
          sumsjptj = sumsjptj + sj*pjet(1:2,j)
        end do
c The sign of pti can take two values:
        do si=-1,1,2
          sumspt = sumsjptj + si*pjet(1:2,i)
          abssumspt = sqrt(sumspt(1)**2 + sumspt(2)**2)
          if (abssumspt_max.lt.abssumspt) then
            abssumspt_max = abssumspt
            sumspt_max = sumspt
          end if
        end do
      end do
c
c We have to calculate the integer sum of pt's:
      sumabspt = 0d0
      do i=1,njet
        sumabspt = sumabspt + sqrt(pjet(1,i)**2 + pjet(2,i)**2)
      end do
c The normal vector can be optained from \sum_{i=1}^{n}s_i p_{\perp\,,i}
c by norming it: 
      nvec = sumspt_max / sqrt(sumspt_max(1)**2 + sumspt_max(2)**2)
c
      Tperp = abssumspt_max/sumabspt
c
c We give back 1 - Tperp, instead:
      Tau = 1d0 - Tperp
c
c We also calculate the transverse Thrust minor:
      Tminor = 0d0
      do i=1,njet
c We have to calculate the cross product of the ith jet pt and nvec:
        Tminor = Tminor + abs(pjet(1,i)*nvec(2) - pjet(2,i)*nvec(1))
      end do
      Tminor = Tminor / sumabspt
c
c      read(*,*)
c
      end subroutine CalcThrustperp
