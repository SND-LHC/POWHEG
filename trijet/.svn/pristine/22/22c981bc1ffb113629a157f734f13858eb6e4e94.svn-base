      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include "coupl.inc"
      integer i
 
      call init_processes_born
      call init_processes_real
      call init_couplings
      if (tmass.eq.0d0) then
         st_nlight=6
      elseif(bmass.eq.0d0) then
         st_nlight=5
      elseif(cmass.eq.0d0) then
         st_nlight=4
      else
         st_nlight=3
      endif
      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.st_nlight) then
            flst_lightpart=i
            exit
         endif
      enddo
 
      return
      end
 
 
 
      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=          -1
      flst_born(   3,   1)=          -1
      flst_born(   4,   1)=          -1
      flst_born(   5,   1)=           0
 
      flst_born(   1,   2)=          -1
      flst_born(   2,   2)=           1
      flst_born(   3,   2)=           1
      flst_born(   4,   2)=          -1
      flst_born(   5,   2)=           0
 
      flst_born(   1,   3)=          -1
      flst_born(   2,   3)=           1
      flst_born(   3,   3)=           2
      flst_born(   4,   3)=          -2
      flst_born(   5,   3)=           0
 
      flst_born(   1,   4)=          -1
      flst_born(   2,   4)=           1
      flst_born(   3,   4)=           0
      flst_born(   4,   4)=           0
      flst_born(   5,   4)=           0
 
      flst_born(   1,   5)=          -1
      flst_born(   2,   5)=          -2
      flst_born(   3,   5)=          -1
      flst_born(   4,   5)=          -2
      flst_born(   5,   5)=           0
 
      flst_born(   1,   6)=          -1
      flst_born(   2,   6)=           2
      flst_born(   3,   6)=          -1
      flst_born(   4,   6)=           2
      flst_born(   5,   6)=           0
 
      flst_born(   1,   7)=          -1
      flst_born(   2,   7)=           0
      flst_born(   3,   7)=           1
      flst_born(   4,   7)=          -1
      flst_born(   5,   7)=          -1
 
      flst_born(   1,   8)=          -1
      flst_born(   2,   8)=           0
      flst_born(   3,   8)=          -1
      flst_born(   4,   8)=           2
      flst_born(   5,   8)=          -2
 
      flst_born(   1,   9)=          -1
      flst_born(   2,   9)=           0
      flst_born(   3,   9)=          -1
      flst_born(   4,   9)=           0
      flst_born(   5,   9)=           0
 
      flst_born(   1,  10)=           1
      flst_born(   2,  10)=          -1
      flst_born(   3,  10)=           1
      flst_born(   4,  10)=          -1
      flst_born(   5,  10)=           0
 
      flst_born(   1,  11)=           1
      flst_born(   2,  11)=          -1
      flst_born(   3,  11)=           2
      flst_born(   4,  11)=          -2
      flst_born(   5,  11)=           0
 
      flst_born(   1,  12)=           1
      flst_born(   2,  12)=          -1
      flst_born(   3,  12)=           0
      flst_born(   4,  12)=           0
      flst_born(   5,  12)=           0
 
      flst_born(   1,  13)=           1
      flst_born(   2,  13)=           1
      flst_born(   3,  13)=           1
      flst_born(   4,  13)=           1
      flst_born(   5,  13)=           0
 
      flst_born(   1,  14)=           1
      flst_born(   2,  14)=          -2
      flst_born(   3,  14)=           1
      flst_born(   4,  14)=          -2
      flst_born(   5,  14)=           0
 
      flst_born(   1,  15)=           1
      flst_born(   2,  15)=           2
      flst_born(   3,  15)=           1
      flst_born(   4,  15)=           2
      flst_born(   5,  15)=           0
 
      flst_born(   1,  16)=           1
      flst_born(   2,  16)=           0
      flst_born(   3,  16)=           1
      flst_born(   4,  16)=           1
      flst_born(   5,  16)=          -1
 
      flst_born(   1,  17)=           1
      flst_born(   2,  17)=           0
      flst_born(   3,  17)=           1
      flst_born(   4,  17)=           2
      flst_born(   5,  17)=          -2
 
      flst_born(   1,  18)=           1
      flst_born(   2,  18)=           0
      flst_born(   3,  18)=           1
      flst_born(   4,  18)=           0
      flst_born(   5,  18)=           0
 
      flst_born(   1,  19)=          -2
      flst_born(   2,  19)=          -1
      flst_born(   3,  19)=          -1
      flst_born(   4,  19)=          -2
      flst_born(   5,  19)=           0
 
      flst_born(   1,  20)=          -2
      flst_born(   2,  20)=           1
      flst_born(   3,  20)=           1
      flst_born(   4,  20)=          -2
      flst_born(   5,  20)=           0
 
      flst_born(   1,  21)=          -2
      flst_born(   2,  21)=          -2
      flst_born(   3,  21)=          -2
      flst_born(   4,  21)=          -2
      flst_born(   5,  21)=           0
 
      flst_born(   1,  22)=          -2
      flst_born(   2,  22)=           2
      flst_born(   3,  22)=           1
      flst_born(   4,  22)=          -1
      flst_born(   5,  22)=           0
 
      flst_born(   1,  23)=          -2
      flst_born(   2,  23)=           2
      flst_born(   3,  23)=           2
      flst_born(   4,  23)=          -2
      flst_born(   5,  23)=           0
 
      flst_born(   1,  24)=          -2
      flst_born(   2,  24)=           2
      flst_born(   3,  24)=           0
      flst_born(   4,  24)=           0
      flst_born(   5,  24)=           0
 
      flst_born(   1,  25)=          -2
      flst_born(   2,  25)=           0
      flst_born(   3,  25)=           1
      flst_born(   4,  25)=          -1
      flst_born(   5,  25)=          -2
 
      flst_born(   1,  26)=          -2
      flst_born(   2,  26)=           0
      flst_born(   3,  26)=           2
      flst_born(   4,  26)=          -2
      flst_born(   5,  26)=          -2
 
      flst_born(   1,  27)=          -2
      flst_born(   2,  27)=           0
      flst_born(   3,  27)=          -2
      flst_born(   4,  27)=           0
      flst_born(   5,  27)=           0
 
      flst_born(   1,  28)=           2
      flst_born(   2,  28)=          -1
      flst_born(   3,  28)=          -1
      flst_born(   4,  28)=           2
      flst_born(   5,  28)=           0
 
      flst_born(   1,  29)=           2
      flst_born(   2,  29)=           1
      flst_born(   3,  29)=           1
      flst_born(   4,  29)=           2
      flst_born(   5,  29)=           0
 
      flst_born(   1,  30)=           2
      flst_born(   2,  30)=          -2
      flst_born(   3,  30)=           1
      flst_born(   4,  30)=          -1
      flst_born(   5,  30)=           0
 
      flst_born(   1,  31)=           2
      flst_born(   2,  31)=          -2
      flst_born(   3,  31)=           2
      flst_born(   4,  31)=          -2
      flst_born(   5,  31)=           0
 
      flst_born(   1,  32)=           2
      flst_born(   2,  32)=          -2
      flst_born(   3,  32)=           0
      flst_born(   4,  32)=           0
      flst_born(   5,  32)=           0
 
      flst_born(   1,  33)=           2
      flst_born(   2,  33)=           2
      flst_born(   3,  33)=           2
      flst_born(   4,  33)=           2
      flst_born(   5,  33)=           0
 
      flst_born(   1,  34)=           2
      flst_born(   2,  34)=           0
      flst_born(   3,  34)=           1
      flst_born(   4,  34)=          -1
      flst_born(   5,  34)=           2
 
      flst_born(   1,  35)=           2
      flst_born(   2,  35)=           0
      flst_born(   3,  35)=           2
      flst_born(   4,  35)=           2
      flst_born(   5,  35)=          -2
 
      flst_born(   1,  36)=           2
      flst_born(   2,  36)=           0
      flst_born(   3,  36)=           2
      flst_born(   4,  36)=           0
      flst_born(   5,  36)=           0
 
      flst_born(   1,  37)=           0
      flst_born(   2,  37)=          -1
      flst_born(   3,  37)=           1
      flst_born(   4,  37)=          -1
      flst_born(   5,  37)=          -1
 
      flst_born(   1,  38)=           0
      flst_born(   2,  38)=          -1
      flst_born(   3,  38)=          -1
      flst_born(   4,  38)=           2
      flst_born(   5,  38)=          -2
 
      flst_born(   1,  39)=           0
      flst_born(   2,  39)=          -1
      flst_born(   3,  39)=          -1
      flst_born(   4,  39)=           0
      flst_born(   5,  39)=           0
 
      flst_born(   1,  40)=           0
      flst_born(   2,  40)=           1
      flst_born(   3,  40)=           1
      flst_born(   4,  40)=           1
      flst_born(   5,  40)=          -1
 
      flst_born(   1,  41)=           0
      flst_born(   2,  41)=           1
      flst_born(   3,  41)=           1
      flst_born(   4,  41)=           2
      flst_born(   5,  41)=          -2
 
      flst_born(   1,  42)=           0
      flst_born(   2,  42)=           1
      flst_born(   3,  42)=           1
      flst_born(   4,  42)=           0
      flst_born(   5,  42)=           0
 
      flst_born(   1,  43)=           0
      flst_born(   2,  43)=          -2
      flst_born(   3,  43)=           1
      flst_born(   4,  43)=          -1
      flst_born(   5,  43)=          -2
 
      flst_born(   1,  44)=           0
      flst_born(   2,  44)=          -2
      flst_born(   3,  44)=           2
      flst_born(   4,  44)=          -2
      flst_born(   5,  44)=          -2
 
      flst_born(   1,  45)=           0
      flst_born(   2,  45)=          -2
      flst_born(   3,  45)=          -2
      flst_born(   4,  45)=           0
      flst_born(   5,  45)=           0
 
      flst_born(   1,  46)=           0
      flst_born(   2,  46)=           2
      flst_born(   3,  46)=           1
      flst_born(   4,  46)=          -1
      flst_born(   5,  46)=           2
 
      flst_born(   1,  47)=           0
      flst_born(   2,  47)=           2
      flst_born(   3,  47)=           2
      flst_born(   4,  47)=           2
      flst_born(   5,  47)=          -2
 
      flst_born(   1,  48)=           0
      flst_born(   2,  48)=           2
      flst_born(   3,  48)=           2
      flst_born(   4,  48)=           0
      flst_born(   5,  48)=           0
 
      flst_born(   1,  49)=           0
      flst_born(   2,  49)=           0
      flst_born(   3,  49)=           1
      flst_born(   4,  49)=          -1
      flst_born(   5,  49)=           0
 
      flst_born(   1,  50)=           0
      flst_born(   2,  50)=           0
      flst_born(   3,  50)=           2
      flst_born(   4,  50)=          -2
      flst_born(   5,  50)=           0
 
      flst_born(   1,  51)=           0
      flst_born(   2,  51)=           0
      flst_born(   3,  51)=           0
      flst_born(   4,  51)=           0
      flst_born(   5,  51)=           0
 
      flst_nborn=          51
 
      return
      end
 
 
 
      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_real(   1,   1)=          -1
      flst_real(   2,   1)=          -1
      flst_real(   3,   1)=           1
      flst_real(   4,   1)=          -1
      flst_real(   5,   1)=          -1
      flst_real(   6,   1)=          -1
 
      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=          -1
      flst_real(   3,   2)=          -1
      flst_real(   4,   2)=          -1
      flst_real(   5,   2)=           2
      flst_real(   6,   2)=          -2
 
      flst_real(   1,   3)=          -1
      flst_real(   2,   3)=          -1
      flst_real(   3,   3)=          -1
      flst_real(   4,   3)=          -1
      flst_real(   5,   3)=           0
      flst_real(   6,   3)=           0
 
      flst_real(   1,   4)=          -1
      flst_real(   2,   4)=           1
      flst_real(   3,   4)=           1
      flst_real(   4,   4)=           1
      flst_real(   5,   4)=          -1
      flst_real(   6,   4)=          -1
 
      flst_real(   1,   5)=          -1
      flst_real(   2,   5)=           1
      flst_real(   3,   5)=           1
      flst_real(   4,   5)=          -1
      flst_real(   5,   5)=           2
      flst_real(   6,   5)=          -2
 
      flst_real(   1,   6)=          -1
      flst_real(   2,   6)=           1
      flst_real(   3,   6)=           1
      flst_real(   4,   6)=          -1
      flst_real(   5,   6)=           0
      flst_real(   6,   6)=           0
 
      flst_real(   1,   7)=          -1
      flst_real(   2,   7)=           1
      flst_real(   3,   7)=           2
      flst_real(   4,   7)=           2
      flst_real(   5,   7)=          -2
      flst_real(   6,   7)=          -2
 
      flst_real(   1,   8)=          -1
      flst_real(   2,   8)=           1
      flst_real(   3,   8)=           2
      flst_real(   4,   8)=          -2
      flst_real(   5,   8)=           0
      flst_real(   6,   8)=           0
 
      flst_real(   1,   9)=          -1
      flst_real(   2,   9)=           1
      flst_real(   3,   9)=           0
      flst_real(   4,   9)=           0
      flst_real(   5,   9)=           0
      flst_real(   6,   9)=           0
 
      flst_real(   1,  10)=          -1
      flst_real(   2,  10)=          -2
      flst_real(   3,  10)=           1
      flst_real(   4,  10)=          -1
      flst_real(   5,  10)=          -1
      flst_real(   6,  10)=          -2
 
      flst_real(   1,  11)=          -1
      flst_real(   2,  11)=          -2
      flst_real(   3,  11)=          -1
      flst_real(   4,  11)=           2
      flst_real(   5,  11)=          -2
      flst_real(   6,  11)=          -2
 
      flst_real(   1,  12)=          -1
      flst_real(   2,  12)=          -2
      flst_real(   3,  12)=          -1
      flst_real(   4,  12)=          -2
      flst_real(   5,  12)=           0
      flst_real(   6,  12)=           0
 
      flst_real(   1,  13)=          -1
      flst_real(   2,  13)=           2
      flst_real(   3,  13)=           1
      flst_real(   4,  13)=          -1
      flst_real(   5,  13)=          -1
      flst_real(   6,  13)=           2
 
      flst_real(   1,  14)=          -1
      flst_real(   2,  14)=           2
      flst_real(   3,  14)=          -1
      flst_real(   4,  14)=           2
      flst_real(   5,  14)=           2
      flst_real(   6,  14)=          -2
 
      flst_real(   1,  15)=          -1
      flst_real(   2,  15)=           2
      flst_real(   3,  15)=          -1
      flst_real(   4,  15)=           2
      flst_real(   5,  15)=           0
      flst_real(   6,  15)=           0
 
      flst_real(   1,  16)=          -1
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=           1
      flst_real(   4,  16)=          -1
      flst_real(   5,  16)=          -1
      flst_real(   6,  16)=           0
 
      flst_real(   1,  17)=          -1
      flst_real(   2,  17)=           0
      flst_real(   3,  17)=          -1
      flst_real(   4,  17)=           2
      flst_real(   5,  17)=          -2
      flst_real(   6,  17)=           0
 
      flst_real(   1,  18)=          -1
      flst_real(   2,  18)=           0
      flst_real(   3,  18)=          -1
      flst_real(   4,  18)=           0
      flst_real(   5,  18)=           0
      flst_real(   6,  18)=           0
 
      flst_real(   1,  19)=           1
      flst_real(   2,  19)=          -1
      flst_real(   3,  19)=           1
      flst_real(   4,  19)=           1
      flst_real(   5,  19)=          -1
      flst_real(   6,  19)=          -1
 
      flst_real(   1,  20)=           1
      flst_real(   2,  20)=          -1
      flst_real(   3,  20)=           1
      flst_real(   4,  20)=          -1
      flst_real(   5,  20)=           2
      flst_real(   6,  20)=          -2
 
      flst_real(   1,  21)=           1
      flst_real(   2,  21)=          -1
      flst_real(   3,  21)=           1
      flst_real(   4,  21)=          -1
      flst_real(   5,  21)=           0
      flst_real(   6,  21)=           0
 
      flst_real(   1,  22)=           1
      flst_real(   2,  22)=          -1
      flst_real(   3,  22)=           2
      flst_real(   4,  22)=           2
      flst_real(   5,  22)=          -2
      flst_real(   6,  22)=          -2
 
      flst_real(   1,  23)=           1
      flst_real(   2,  23)=          -1
      flst_real(   3,  23)=           2
      flst_real(   4,  23)=          -2
      flst_real(   5,  23)=           0
      flst_real(   6,  23)=           0
 
      flst_real(   1,  24)=           1
      flst_real(   2,  24)=          -1
      flst_real(   3,  24)=           0
      flst_real(   4,  24)=           0
      flst_real(   5,  24)=           0
      flst_real(   6,  24)=           0
 
      flst_real(   1,  25)=           1
      flst_real(   2,  25)=           1
      flst_real(   3,  25)=           1
      flst_real(   4,  25)=           1
      flst_real(   5,  25)=           1
      flst_real(   6,  25)=          -1
 
      flst_real(   1,  26)=           1
      flst_real(   2,  26)=           1
      flst_real(   3,  26)=           1
      flst_real(   4,  26)=           1
      flst_real(   5,  26)=           2
      flst_real(   6,  26)=          -2
 
      flst_real(   1,  27)=           1
      flst_real(   2,  27)=           1
      flst_real(   3,  27)=           1
      flst_real(   4,  27)=           1
      flst_real(   5,  27)=           0
      flst_real(   6,  27)=           0
 
      flst_real(   1,  28)=           1
      flst_real(   2,  28)=          -2
      flst_real(   3,  28)=           1
      flst_real(   4,  28)=           1
      flst_real(   5,  28)=          -1
      flst_real(   6,  28)=          -2
 
      flst_real(   1,  29)=           1
      flst_real(   2,  29)=          -2
      flst_real(   3,  29)=           1
      flst_real(   4,  29)=           2
      flst_real(   5,  29)=          -2
      flst_real(   6,  29)=          -2
 
      flst_real(   1,  30)=           1
      flst_real(   2,  30)=          -2
      flst_real(   3,  30)=           1
      flst_real(   4,  30)=          -2
      flst_real(   5,  30)=           0
      flst_real(   6,  30)=           0
 
      flst_real(   1,  31)=           1
      flst_real(   2,  31)=           2
      flst_real(   3,  31)=           1
      flst_real(   4,  31)=           1
      flst_real(   5,  31)=          -1
      flst_real(   6,  31)=           2
 
      flst_real(   1,  32)=           1
      flst_real(   2,  32)=           2
      flst_real(   3,  32)=           1
      flst_real(   4,  32)=           2
      flst_real(   5,  32)=           2
      flst_real(   6,  32)=          -2
 
      flst_real(   1,  33)=           1
      flst_real(   2,  33)=           2
      flst_real(   3,  33)=           1
      flst_real(   4,  33)=           2
      flst_real(   5,  33)=           0
      flst_real(   6,  33)=           0
 
      flst_real(   1,  34)=           1
      flst_real(   2,  34)=           0
      flst_real(   3,  34)=           1
      flst_real(   4,  34)=           1
      flst_real(   5,  34)=          -1
      flst_real(   6,  34)=           0
 
      flst_real(   1,  35)=           1
      flst_real(   2,  35)=           0
      flst_real(   3,  35)=           1
      flst_real(   4,  35)=           2
      flst_real(   5,  35)=          -2
      flst_real(   6,  35)=           0
 
      flst_real(   1,  36)=           1
      flst_real(   2,  36)=           0
      flst_real(   3,  36)=           1
      flst_real(   4,  36)=           0
      flst_real(   5,  36)=           0
      flst_real(   6,  36)=           0
 
      flst_real(   1,  37)=          -2
      flst_real(   2,  37)=          -1
      flst_real(   3,  37)=           1
      flst_real(   4,  37)=          -1
      flst_real(   5,  37)=          -1
      flst_real(   6,  37)=          -2
 
      flst_real(   1,  38)=          -2
      flst_real(   2,  38)=          -1
      flst_real(   3,  38)=          -1
      flst_real(   4,  38)=           2
      flst_real(   5,  38)=          -2
      flst_real(   6,  38)=          -2
 
      flst_real(   1,  39)=          -2
      flst_real(   2,  39)=          -1
      flst_real(   3,  39)=          -1
      flst_real(   4,  39)=          -2
      flst_real(   5,  39)=           0
      flst_real(   6,  39)=           0
 
      flst_real(   1,  40)=          -2
      flst_real(   2,  40)=           1
      flst_real(   3,  40)=           1
      flst_real(   4,  40)=           1
      flst_real(   5,  40)=          -1
      flst_real(   6,  40)=          -2
 
      flst_real(   1,  41)=          -2
      flst_real(   2,  41)=           1
      flst_real(   3,  41)=           1
      flst_real(   4,  41)=           2
      flst_real(   5,  41)=          -2
      flst_real(   6,  41)=          -2
 
      flst_real(   1,  42)=          -2
      flst_real(   2,  42)=           1
      flst_real(   3,  42)=           1
      flst_real(   4,  42)=          -2
      flst_real(   5,  42)=           0
      flst_real(   6,  42)=           0
 
      flst_real(   1,  43)=          -2
      flst_real(   2,  43)=          -2
      flst_real(   3,  43)=           1
      flst_real(   4,  43)=          -1
      flst_real(   5,  43)=          -2
      flst_real(   6,  43)=          -2
 
      flst_real(   1,  44)=          -2
      flst_real(   2,  44)=          -2
      flst_real(   3,  44)=           2
      flst_real(   4,  44)=          -2
      flst_real(   5,  44)=          -2
      flst_real(   6,  44)=          -2
 
      flst_real(   1,  45)=          -2
      flst_real(   2,  45)=          -2
      flst_real(   3,  45)=          -2
      flst_real(   4,  45)=          -2
      flst_real(   5,  45)=           0
      flst_real(   6,  45)=           0
 
      flst_real(   1,  46)=          -2
      flst_real(   2,  46)=           2
      flst_real(   3,  46)=           1
      flst_real(   4,  46)=           1
      flst_real(   5,  46)=          -1
      flst_real(   6,  46)=          -1
 
      flst_real(   1,  47)=          -2
      flst_real(   2,  47)=           2
      flst_real(   3,  47)=           1
      flst_real(   4,  47)=          -1
      flst_real(   5,  47)=           2
      flst_real(   6,  47)=          -2
 
      flst_real(   1,  48)=          -2
      flst_real(   2,  48)=           2
      flst_real(   3,  48)=           1
      flst_real(   4,  48)=          -1
      flst_real(   5,  48)=           0
      flst_real(   6,  48)=           0
 
      flst_real(   1,  49)=          -2
      flst_real(   2,  49)=           2
      flst_real(   3,  49)=           2
      flst_real(   4,  49)=           2
      flst_real(   5,  49)=          -2
      flst_real(   6,  49)=          -2
 
      flst_real(   1,  50)=          -2
      flst_real(   2,  50)=           2
      flst_real(   3,  50)=           2
      flst_real(   4,  50)=          -2
      flst_real(   5,  50)=           0
      flst_real(   6,  50)=           0
 
      flst_real(   1,  51)=          -2
      flst_real(   2,  51)=           2
      flst_real(   3,  51)=           0
      flst_real(   4,  51)=           0
      flst_real(   5,  51)=           0
      flst_real(   6,  51)=           0
 
      flst_real(   1,  52)=          -2
      flst_real(   2,  52)=           0
      flst_real(   3,  52)=           1
      flst_real(   4,  52)=          -1
      flst_real(   5,  52)=          -2
      flst_real(   6,  52)=           0
 
      flst_real(   1,  53)=          -2
      flst_real(   2,  53)=           0
      flst_real(   3,  53)=           2
      flst_real(   4,  53)=          -2
      flst_real(   5,  53)=          -2
      flst_real(   6,  53)=           0
 
      flst_real(   1,  54)=          -2
      flst_real(   2,  54)=           0
      flst_real(   3,  54)=          -2
      flst_real(   4,  54)=           0
      flst_real(   5,  54)=           0
      flst_real(   6,  54)=           0
 
      flst_real(   1,  55)=           2
      flst_real(   2,  55)=          -1
      flst_real(   3,  55)=           1
      flst_real(   4,  55)=          -1
      flst_real(   5,  55)=          -1
      flst_real(   6,  55)=           2
 
      flst_real(   1,  56)=           2
      flst_real(   2,  56)=          -1
      flst_real(   3,  56)=          -1
      flst_real(   4,  56)=           2
      flst_real(   5,  56)=           2
      flst_real(   6,  56)=          -2
 
      flst_real(   1,  57)=           2
      flst_real(   2,  57)=          -1
      flst_real(   3,  57)=          -1
      flst_real(   4,  57)=           2
      flst_real(   5,  57)=           0
      flst_real(   6,  57)=           0
 
      flst_real(   1,  58)=           2
      flst_real(   2,  58)=           1
      flst_real(   3,  58)=           1
      flst_real(   4,  58)=           1
      flst_real(   5,  58)=          -1
      flst_real(   6,  58)=           2
 
      flst_real(   1,  59)=           2
      flst_real(   2,  59)=           1
      flst_real(   3,  59)=           1
      flst_real(   4,  59)=           2
      flst_real(   5,  59)=           2
      flst_real(   6,  59)=          -2
 
      flst_real(   1,  60)=           2
      flst_real(   2,  60)=           1
      flst_real(   3,  60)=           1
      flst_real(   4,  60)=           2
      flst_real(   5,  60)=           0
      flst_real(   6,  60)=           0
 
      flst_real(   1,  61)=           2
      flst_real(   2,  61)=          -2
      flst_real(   3,  61)=           1
      flst_real(   4,  61)=           1
      flst_real(   5,  61)=          -1
      flst_real(   6,  61)=          -1
 
      flst_real(   1,  62)=           2
      flst_real(   2,  62)=          -2
      flst_real(   3,  62)=           1
      flst_real(   4,  62)=          -1
      flst_real(   5,  62)=           2
      flst_real(   6,  62)=          -2
 
      flst_real(   1,  63)=           2
      flst_real(   2,  63)=          -2
      flst_real(   3,  63)=           1
      flst_real(   4,  63)=          -1
      flst_real(   5,  63)=           0
      flst_real(   6,  63)=           0
 
      flst_real(   1,  64)=           2
      flst_real(   2,  64)=          -2
      flst_real(   3,  64)=           2
      flst_real(   4,  64)=           2
      flst_real(   5,  64)=          -2
      flst_real(   6,  64)=          -2
 
      flst_real(   1,  65)=           2
      flst_real(   2,  65)=          -2
      flst_real(   3,  65)=           2
      flst_real(   4,  65)=          -2
      flst_real(   5,  65)=           0
      flst_real(   6,  65)=           0
 
      flst_real(   1,  66)=           2
      flst_real(   2,  66)=          -2
      flst_real(   3,  66)=           0
      flst_real(   4,  66)=           0
      flst_real(   5,  66)=           0
      flst_real(   6,  66)=           0
 
      flst_real(   1,  67)=           2
      flst_real(   2,  67)=           2
      flst_real(   3,  67)=           1
      flst_real(   4,  67)=          -1
      flst_real(   5,  67)=           2
      flst_real(   6,  67)=           2
 
      flst_real(   1,  68)=           2
      flst_real(   2,  68)=           2
      flst_real(   3,  68)=           2
      flst_real(   4,  68)=           2
      flst_real(   5,  68)=           2
      flst_real(   6,  68)=          -2
 
      flst_real(   1,  69)=           2
      flst_real(   2,  69)=           2
      flst_real(   3,  69)=           2
      flst_real(   4,  69)=           2
      flst_real(   5,  69)=           0
      flst_real(   6,  69)=           0
 
      flst_real(   1,  70)=           2
      flst_real(   2,  70)=           0
      flst_real(   3,  70)=           1
      flst_real(   4,  70)=          -1
      flst_real(   5,  70)=           2
      flst_real(   6,  70)=           0
 
      flst_real(   1,  71)=           2
      flst_real(   2,  71)=           0
      flst_real(   3,  71)=           2
      flst_real(   4,  71)=           2
      flst_real(   5,  71)=          -2
      flst_real(   6,  71)=           0
 
      flst_real(   1,  72)=           2
      flst_real(   2,  72)=           0
      flst_real(   3,  72)=           2
      flst_real(   4,  72)=           0
      flst_real(   5,  72)=           0
      flst_real(   6,  72)=           0
 
      flst_real(   1,  73)=           0
      flst_real(   2,  73)=          -1
      flst_real(   3,  73)=           1
      flst_real(   4,  73)=          -1
      flst_real(   5,  73)=          -1
      flst_real(   6,  73)=           0
 
      flst_real(   1,  74)=           0
      flst_real(   2,  74)=          -1
      flst_real(   3,  74)=          -1
      flst_real(   4,  74)=           2
      flst_real(   5,  74)=          -2
      flst_real(   6,  74)=           0
 
      flst_real(   1,  75)=           0
      flst_real(   2,  75)=          -1
      flst_real(   3,  75)=          -1
      flst_real(   4,  75)=           0
      flst_real(   5,  75)=           0
      flst_real(   6,  75)=           0
 
      flst_real(   1,  76)=           0
      flst_real(   2,  76)=           1
      flst_real(   3,  76)=           1
      flst_real(   4,  76)=           1
      flst_real(   5,  76)=          -1
      flst_real(   6,  76)=           0
 
      flst_real(   1,  77)=           0
      flst_real(   2,  77)=           1
      flst_real(   3,  77)=           1
      flst_real(   4,  77)=           2
      flst_real(   5,  77)=          -2
      flst_real(   6,  77)=           0
 
      flst_real(   1,  78)=           0
      flst_real(   2,  78)=           1
      flst_real(   3,  78)=           1
      flst_real(   4,  78)=           0
      flst_real(   5,  78)=           0
      flst_real(   6,  78)=           0
 
      flst_real(   1,  79)=           0
      flst_real(   2,  79)=          -2
      flst_real(   3,  79)=           1
      flst_real(   4,  79)=          -1
      flst_real(   5,  79)=          -2
      flst_real(   6,  79)=           0
 
      flst_real(   1,  80)=           0
      flst_real(   2,  80)=          -2
      flst_real(   3,  80)=           2
      flst_real(   4,  80)=          -2
      flst_real(   5,  80)=          -2
      flst_real(   6,  80)=           0
 
      flst_real(   1,  81)=           0
      flst_real(   2,  81)=          -2
      flst_real(   3,  81)=          -2
      flst_real(   4,  81)=           0
      flst_real(   5,  81)=           0
      flst_real(   6,  81)=           0
 
      flst_real(   1,  82)=           0
      flst_real(   2,  82)=           2
      flst_real(   3,  82)=           1
      flst_real(   4,  82)=          -1
      flst_real(   5,  82)=           2
      flst_real(   6,  82)=           0
 
      flst_real(   1,  83)=           0
      flst_real(   2,  83)=           2
      flst_real(   3,  83)=           2
      flst_real(   4,  83)=           2
      flst_real(   5,  83)=          -2
      flst_real(   6,  83)=           0
 
      flst_real(   1,  84)=           0
      flst_real(   2,  84)=           2
      flst_real(   3,  84)=           2
      flst_real(   4,  84)=           0
      flst_real(   5,  84)=           0
      flst_real(   6,  84)=           0
 
      flst_real(   1,  85)=           0
      flst_real(   2,  85)=           0
      flst_real(   3,  85)=           1
      flst_real(   4,  85)=           1
      flst_real(   5,  85)=          -1
      flst_real(   6,  85)=          -1
 
      flst_real(   1,  86)=           0
      flst_real(   2,  86)=           0
      flst_real(   3,  86)=           1
      flst_real(   4,  86)=          -1
      flst_real(   5,  86)=           2
      flst_real(   6,  86)=          -2
 
      flst_real(   1,  87)=           0
      flst_real(   2,  87)=           0
      flst_real(   3,  87)=           1
      flst_real(   4,  87)=          -1
      flst_real(   5,  87)=           0
      flst_real(   6,  87)=           0
 
      flst_real(   1,  88)=           0
      flst_real(   2,  88)=           0
      flst_real(   3,  88)=           2
      flst_real(   4,  88)=           2
      flst_real(   5,  88)=          -2
      flst_real(   6,  88)=          -2
 
      flst_real(   1,  89)=           0
      flst_real(   2,  89)=           0
      flst_real(   3,  89)=           2
      flst_real(   4,  89)=          -2
      flst_real(   5,  89)=           0
      flst_real(   6,  89)=           0
 
      flst_real(   1,  90)=           0
      flst_real(   2,  90)=           0
      flst_real(   3,  90)=           0
      flst_real(   4,  90)=           0
      flst_real(   5,  90)=           0
      flst_real(   6,  90)=           0
 
      flst_nreal=          90
 
      return
      end
 
