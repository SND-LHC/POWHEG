 1    subroutine born_phsp(xborn)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'coupl.inc'
      real * 8 xborn(ndiminteg-3),jac,tmpvec(0:3),bwcutoff
      integer i,j,k
      logical ini,fullphsp
      data ini/.true./
      save ini,fullphsp
      real * 8 powheginput
c
      real*8 pt1,pt2,pt3,m12,m13,m23
      real*8 ptsqrel
      real*8 dotp
      external dotp
      ptsqrel(i,j)=2*dotp(kn_cmpborn(0,i),kn_cmpborn(0,j))
     1     *  kn_cmpborn(0,i)*kn_cmpborn(0,j)
     2     /(kn_cmpborn(0,i)**2+kn_cmpborn(0,j)**2)
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         fullphsp = .true.
         if(powheginput("#fullphsp").eq.0) then
            fullphsp = .false.
         endif
c We read in the min kt: 
         brkn_ktmin = powheginput("bornktmin")
         Print *,"******************************************"
         Print *,"*                                        *"
         Print *,"* kT min for Born phase-space generation *"
         Print *,"*   ktmin = ",brkn_ktmin
         Print *,"*                                        *"
         Print *,"******************************************"
         ini=.false.
      endif      

c      brkn_ktmin=0
      call born_phsp_jj(xborn)

      if(fullphsp) then
c We construct the 3-parton final state with final state radiation if
c xi is greater than 1/3, since we have 9 different regions, and 3 is
c associated to the initial states:
c         print *,xborn(ndiminteg-3)
c         print *,int(xborn(ndiminteg-3)*9)+1
         if (xborn(ndiminteg-3).gt.1d0/3d0) then
c            print *,"Region 3 -> 35"
c            print *,"4th region"
c We always consider the 3 parton (the 1st one in the final state) as
c the emitter:
            brkn_emitter=3
c The emitted parton will be inserted as the fifth one:
            call br_real_phsp_fsr(xborn(ndiminteg-6),jac)
            jac = jac * 1/brkn_dijterm(3,5) 
     >          / (1/brkn_dijterm(0,3) 
     >          +  1/brkn_dijterm(0,4)
     >          +  1/brkn_dijterm(0,5)
     >          +  1/brkn_dijterm(3,5)
     >          +  1/brkn_dijterm(4,5)
     >          +  1/brkn_dijterm(3,4))
     >        * brkn_cmpreal(0,5)/(brkn_cmpreal(0,3)+brkn_cmpreal(0,5))
            kn_cmpborn=brkn_cmpreal
            kn_pborn=brkn_preal
            kn_xb1=brkn_x1
            kn_xb2=brkn_x2
c In the original splitting parton n. 3 was the emitter, after swapping
c it is the n. 5:
            if ((xborn(ndiminteg-3).lt.5d0/9d0).and.      
     >          (xborn(ndiminteg-3).gt.4d0/9d0)) then
c              print *,"Region 5 -> 53"
c              print *,"5th region"
c The change needed: 3 <-> 5
              call swap_momenta(kn_cmpborn(:,3),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,3),kn_pborn(:,5))
c This time we consider parton n. 4 as the emitter and n. 5 as the
c emitted one:
            elseif ((xborn(ndiminteg-3).lt.6d0/9d0).and.
     >              (xborn(ndiminteg-3).gt.5d0/9d0)) then
c              print *,"Region 4 -> 45"
c              print *,"6th region"
c The change needed: 3 <-> 4
              call swap_momenta(kn_cmpborn(:,3),kn_cmpborn(:,4))
              call swap_momenta(kn_pborn(:,3),kn_pborn(:,4))
c This time we consider parton n. 5 as the emitter and n. 4 as the
c emitted one:
            elseif ((xborn(ndiminteg-3).lt.7d0/9d0).and.
     >              (xborn(ndiminteg-3).gt.6d0/9d0)) then
c              print *,"Region 5 -> 54"
c              print *,"7th region"
c The change needed: 3 <-> 5
              call swap_momenta(kn_cmpborn(:,3),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,3),kn_pborn(:,5))
c Now we make the changes: 3 <-> 4, 
              call swap_momenta(kn_cmpborn(:,3),kn_cmpborn(:,4))
              call swap_momenta(kn_pborn(:,3),kn_pborn(:,4))
c This time we consider parton n. 3 as the emitter and n. 4 as the
c emitted one:
            elseif ((xborn(ndiminteg-3).lt.8d0/9d0).and.
     >              (xborn(ndiminteg-3).gt.7d0/9d0)) then
c              print *,"Region 3 -> 34"
c              print *,"8th region"
c The change needed: 4 <-> 5
              call swap_momenta(kn_cmpborn(:,4),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,4),kn_pborn(:,5))
c This time we consider parton n. 4 as the emitter and n. 3 as the
c emitted one:
            elseif (xborn(ndiminteg-3).gt.8d0/9d0) then
c              print *,"Region 4 -> 43"
c              print *,"9th region"
c The change needed: 3 <-> 5
              call swap_momenta(kn_cmpborn(:,3),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,3),kn_pborn(:,5))
c Now we make the changes: 4 <-> 5, 
              call swap_momenta(kn_cmpborn(:,4),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,4),kn_pborn(:,5))
            endif
          else
c We consider the case when parton n. 5 is emitted:
            brkn_emitter=0
c            print *,"1st region"
            call br_real_phsp_isr(xborn(ndiminteg-6),jac)
            jac = jac * 1/brkn_dijterm(0,5)
     >          / (1/brkn_dijterm(0,3) 
     >          +  1/brkn_dijterm(0,4)
     >          +  1/brkn_dijterm(0,5)
     >          +  1/brkn_dijterm(3,5)
     >          +  1/brkn_dijterm(4,5)
     >          +  1/brkn_dijterm(3,4))
            kn_cmpborn=brkn_cmpreal
            kn_pborn=brkn_preal
            kn_xb1=brkn_x1
            kn_xb2=brkn_x2
c This time the emitted parton is the 3rd one:
            if ((xborn(ndiminteg-3).lt.2d0/9d0).and. 
     >          (xborn(ndiminteg-3).gt.1d0/9d0)) then
c              print *,"2nd region"
c The change needed: 3 <-> 5
              call swap_momenta(kn_cmpborn(:,3),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,3),kn_pborn(:,5))
c This time the emitted parton is the 4th one:
            elseif (xborn(ndiminteg-3).gt.2d0/9d0) then
c              print *,"3rd region"
c The change needed: 4 <-> 5
              call swap_momenta(kn_cmpborn(:,4),kn_cmpborn(:,5))
              call swap_momenta(kn_pborn(:,4),kn_pborn(:,5))
            endif
         endif
      else
         brkn_emitter=0
         call br_real_phsp_isr(xborn(ndiminteg-6),jac)
         kn_cmpborn=brkn_cmpreal
         kn_pborn=brkn_preal
         kn_xb1=brkn_x1
         kn_xb2=brkn_x2
         jac=jac/9
      endif
c factor of 9 to compensate for the 1/9 range for each
c of the above contributions.
      kn_jacborn=9*brkn_jacborn*jac

c     set the CMS energy 
      kn_sborn=brkn_sreal

c     minimal final state mass 
c     kn_minmass=sqrt(ph_Hmass2low)
c      kn_minmass=kn_ktmin + sqrt(kn_ktmin**2 + kn_masses(3)**2)
      kn_minmass=2*kn_ktmin
c
c We include an extra cut on the pts and invariant masses of
c final-state particles
      pt1 = sqrt(kn_cmpborn(1,3)**2 + kn_cmpborn(2,3)**2)
      pt2 = sqrt(kn_cmpborn(1,4)**2 + kn_cmpborn(2,4)**2)
      pt3 = sqrt(kn_cmpborn(1,5)**2 + kn_cmpborn(2,5)**2)
      m12 = sqrt(abs(ptsqrel(3,4)))
      m13 = sqrt(abs(ptsqrel(3,5)))
      m23 = sqrt(abs(ptsqrel(4,5)))
c We set the jacobian to zero for sufficiently small pts and
c invariant masses
      if (min(pt1,pt2,pt3,m12,m13,m23).lt.0.3d0) then
        kn_jacborn = 0d0
        call increasecnt("Small pt or inv. mass")
      end if
      end

         

      subroutine born_phsp_jj(xborn)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'

      real * 8 xborn(ndiminteg-4)
      real * 8 vec(3),phidec
      integer mu,k
      logical ini
      data ini/.true./
      save ini

      real*8 xx(ndiminteg-3)
      real*8 tau,tau_min,tau_max
      real*8 ycm,ycm_min,ycm_max
      real*8 cth1,cth1_min,cth1_max,lcthl,lcthh,lcth,tcth
      real*8 powheginput,beta,tmp
      integer ixx_tau,ixx_ycm,ixx_cth1
      parameter(
     $     ixx_tau  =1,
     $     ixx_ycm  =2,
     $     ixx_cth1 =3)

C -   N.B: To check the phase space volume is correct,
C -   replace res(j) in sigborn.f by res(j)=kn_jacborn
C -   then check that the cross section is equal to,
C -   hc^2/(8*pi) = 1.54929*10^7.

C -   Parameter to select phase space importance sampling:
C -   psgen=0:     flat in 1/tau, flat in 1/(1-cth1)
C -   psgen=1:     flat in tau, flat in cth1
c -   psgen=2:     flat in log tau, flat in log((1+cth1)/(1-cth1))
      integer psgen
      parameter (psgen=2)

C -   Phase space: 1 /(16 pi) d tau d y d cth 
C -   4*kn_ktmin**2 / S < tau < 1
C -   |y| < ln(tau)/2
C -   |cth1| < sqrt(1-4*kn_ktmin**2 / S /tau)

      if (ini) then
         write(*,*) '************************************'
         write(*,*) '  Minimum kT of jets in underlying  '
         write(*,*) '  Born is limited to',kn_ktmin,' GeV'
         write(*,*) '************************************'
      endif

C -   Set minimum sqrt(\hat{s}) value
      kn_minmass = 2.*kn_ktmin

C -   Set initial- and final-state masses for Born and real
      if(ini) then
         do k=1,br_nlegborn
            brkn_masses(k)=0
         enddo
         ini=.false.
      endif

C -   Making a local copy of the unit random variables (?)
      do k=1,ndiminteg-4
         xx(k)=xborn(k)
      enddo

C -   Initial Jacobian prefactor
      brkn_jacborn = 1.0d0/16.0d0/pi
      
C -   Set minimum and maximum tau values then get tau:
      tau_min = 4.*brkn_ktmin**2 / kn_sbeams
      tau_max = 1.0d0
      if(psgen.eq.0)then
C -      Sampling flat in 1/tau
         tmp  = 1.0d0/ tau_max + xx(ixx_tau) * (1d0/tau_min-1d0/tau_max)
         tau  = 1.0d0/ tmp
         brkn_jacborn = brkn_jacborn * tau**2 *(1d0/tau_min-1d0/tau_max)
      elseif(psgen.eq.1) then
C -      Sampling flat in tau
         tau  = tau_min + xx(ixx_tau) * (tau_max-tau_min)
         brkn_jacborn = brkn_jacborn * (tau_max-tau_min)
      elseif(psgen.eq.2) then
C -      Flat in log(tau)
         tau = tau_min*exp(xx(ixx_tau) * log(tau_max/tau_min))
         brkn_jacborn = brkn_jacborn * tau * log(tau_max/tau_min)
      else
         write(*,*) 'Wrong psgen in gen_born_vars'
         call exit(1)
      endif

C -   Now we can work out \hat{s}
      brkn_sborn = tau * kn_sbeams

C -   Set minimum and maximum y values 
      ycm_min =  log(tau)/2.
      ycm_max = -log(tau)/2.
C -   Now get a y uniformly in the allowed range
      ycm  = ycm_min + xx(ixx_ycm) * (ycm_max-ycm_min)
      brkn_jacborn = brkn_jacborn * (ycm_max-ycm_min)

C -   Feynman x's can be computed and saved now:
      brkn_xb1 = sqrt(tau) * exp( ycm)
      brkn_xb2 = sqrt(tau) * exp(-ycm)

C -   Set minimum and maximum cos(theta) values (theta is the 
C -   angle of 1st outgoing particle wrt +z axis).
      tmp = sqrt(brkn_sborn)
      cth1_max =  sqrt(1.-2.d0*brkn_ktmin/tmp)
     1     *sqrt(1.+2.d0*brkn_ktmin/tmp)
      cth1_min = -cth1_max
      if(psgen.eq.0) then
C -      Sample cos(theta) as 1/(1-cos(theta)) ~ 1/sin^4(theta/2)
         tmp  = 2.0D0*cth1_max/(1.0D0-cth1_max)/(1.0D0+cth1_max) 
         cth1 = 1.0D0-1.0D0/(1.0D0/(1.0D0-cth1_min) + xx(ixx_cth1)*tmp)
         brkn_jacborn = brkn_jacborn * (1D0-cth1)**2 * tmp
      elseif(psgen.eq.1) then
C -      Sample cos(theta) uniformly
         cth1 = cth1_min + xx(ixx_cth1) * (cth1_max-cth1_min)
         brkn_jacborn = brkn_jacborn * (cth1_max-cth1_min)
      elseif(psgen.eq.2) then
c -      Sample cos(theta) uniformly in log((1-cos)/(1+cos))
         lcthl=log((1-cth1_max)/(1+cth1_max))
         lcthh=-lcthl
         lcth=xx(ixx_cth1)*(lcthh-lcthl)+lcthl
         brkn_jacborn = brkn_jacborn * (lcthh-lcthl)
         tcth=exp(lcth)
         brkn_jacborn = brkn_jacborn * tcth
         cth1=(tcth-1)/(tcth+1)
         brkn_jacborn = brkn_jacborn * 2/(tcth+1)**2
      else
         write(*,*) 'Wrong psgen in born_phsp_jj'
         call pwhg_exit(-1)
      endif
 
      if(brkn_jacborn.lt.0.0) then
         write(6,*) 'born_phsp: Error brkn_jacborn < 0. !'
         write(*,*) 'brkn_jacborn        = ', brkn_jacborn
         write(*,*) 'sqrt(brkn_sborn)    = ', sqrt(brkn_sborn)
         write(*,*) 'cth1_max-cth1_min = ', (cth1_max-cth1_min)
         call pwhg_exit(-1)
      endif

C -   Assign born jacobian and default kinematics variables
c      brkn_born_pt2 = (1d0-cth1**2)*brkn_sborn/4.0d0
c      kn_cthdec   = cth1
      phidec      = 0d0

C -   With this choice px_3 is always positive, but at the
C -   end the whole event will be randomly rotated around z-axis

C -   Initial state momenta:
      brkn_cmpborn(0,1) =  sqrt(brkn_sborn)/2        ! E
      brkn_cmpborn(0,2) =  brkn_cmpborn(0,1)         !
      brkn_cmpborn(3,1) =  brkn_cmpborn(0,1)         ! Pz
      brkn_cmpborn(3,2) = -brkn_cmpborn(0,2)         !
      brkn_cmpborn(1,1) = 0                        ! Px
      brkn_cmpborn(1,2) = 0                        !
      brkn_cmpborn(2,1) = 0                        ! Py
      brkn_cmpborn(2,2) = 0                        !

C -   Final state momenta:
      brkn_cmpborn(0,3) = sqrt(brkn_sborn)/2            ! E
      brkn_cmpborn(0,4) =  brkn_cmpborn(0,3)            !
      brkn_cmpborn(3,3) =  cth1*brkn_cmpborn(0,3)  ! Pz
      brkn_cmpborn(3,4) = -brkn_cmpborn(3,3)            !
      brkn_cmpborn(1,3) = sqrt(1-cth1**2)*sin(phidec)
     1     *brkn_cmpborn(0,3)   ! Px
      brkn_cmpborn(1,4) = -brkn_cmpborn(1,3)                                  !
      brkn_cmpborn(2,3) = sqrt(1-cth1**2)*cos(phidec)
     1     *brkn_cmpborn(0,3)   ! Py
      brkn_cmpborn(2,4) = -brkn_cmpborn(2,3)                                  !

C -   Boost to the LAB frame:
      beta=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
c We use the br_nlegborn variable which contains one less leg as
c compared to nlegborn:
c note that br_nlegborn-2 specifies how many 4-momenta are needed to
c be boosted:
      call mboost(br_nlegborn-2,vec,beta,brkn_cmpborn(0,3),
     >            brkn_pborn(0,3))
      do mu=0,3
         brkn_pborn(mu,1)=brkn_xb1*kn_beams(mu,1)
         brkn_pborn(mu,2)=brkn_xb2*kn_beams(mu,2)
      enddo
c We have to use the br_nlegborn to check for momenta:
      call checkmomzero(br_nlegborn,brkn_pborn)
C -      call checkmass(2,brkn_pborn(0,3))
      call br_compute_csimax_fsr

      end


c Mappings of the underlying born configuration in
c brkn_cmpborn(0:3,br_nlegborn), and the xrad(1:3) variables
c in the unit cube, into brkn_real(0:3,nlegreal).
c The factor jac_over_csi*csi*brkn_csimax, multiplied
c by the Born phase space jacobian, yields the real phase
c space jacobian.
c More explicitly:
c d Phi_n = d^3 xrad jac_over_csi csi csimax d Phi_{n-1}
c Since
c  d Phi_n = d phi d y d csi Jrad d Phi_{n-1}
c (where Jrad is given in FNO2006) we get
c                                  d phi d y d csi
c csimax csi jac_over_csi = Jrad  ----------------
c                                    d^3 xrad
c Notice that using d csi=d csitilde csimax the csimax
c factor cancels, and jac_over_csi is as given in the
c code below (see notes on xscaled.tm).
c br_real_phsp_fsr: provides the mapping for the final state
c radiation, assuming that the emitter is the brkn_emitter-th
c particle, and the emitted particle is the nlegreal-th particle
c br_real_phsp_isr: mapping for the initial state radiation
      subroutine br_real_phsp_fsr(xrad,jac)
      implicit none
      real * 8 xrad(3),jac
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 q0,q2,xjac
c Boost the underlying Born variables to their cm frame
      q0=2*brkn_cmpborn(0,1)
      q2=brkn_sborn
      brkn_csitilde=xrad(1)
      xjac=1
      brkn_y=1-2*xrad(2)
      xjac=xjac*2
c importance sampling for brkn_y
      xjac=xjac*1.5d0*(1-brkn_y**2)
      brkn_y=1.5d0*(brkn_y-brkn_y**3/3)
      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      brkn_csimax=brkn_csimax_arr(brkn_emitter)
      brkn_csi=brkn_csitilde*brkn_csimax     
c remember: no csimax in the jacobian factor, we are integrating in csitilde 
      call br_real_phsp_fsr_rad
      jac=xjac*brkn_jacreal*brkn_csimax
      end

c br_real_phsp_fsr_rad: provides the mapping for the final state
c radiation, assuming that we are considering the region rad_kinreg
c and the emitted particle is the nlegreal-th particle,
c for given brkn_csi, brkn_y, brkn_azi. Sets the jacobian
c brkn_jacreal so that brkn_jacreal d brkn_csi d brkn_y d brkn_azi times
c the underlying Born jacobian is the phase space volume
      subroutine br_real_phsp_fsr_rad
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 vec(3),q0,beta
      integer i
      data vec/0d0,0d0,1d0/
      save vec
      q0=2*brkn_cmpborn(0,1)
c remember: no csimax factor, we are integrating in csitilde 
      call barradmap(br_nlegborn-2,brkn_emitter-2,q0,brkn_cmpborn(0,3),
     1    brkn_csi,brkn_y,brkn_azi,brkn_preal(0,3),brkn_jacreal)
      beta=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      call mboost(br_nlegreal-2,vec,beta,
     1    brkn_preal(0,3),brkn_preal(0,3))
      do i=0,3
         brkn_preal(i,1)=brkn_pborn(i,1)
         brkn_preal(i,2)=brkn_pborn(i,2)
      enddo
      brkn_x1=brkn_xb1
      brkn_x2=brkn_xb2
      brkn_sreal=brkn_sborn
c      call checkmomzero(br_nlegreal,brkn_preal)
      call br_compcmkin
      call br_compdij
      end

      subroutine br_real_phsp_isr(xrad,jac)
      implicit none
      real * 8 xrad(3),jac
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 xjac
      brkn_csitilde=(3-2*xrad(1))*xrad(1)**2
      xjac=6*(1-xrad(1))*xrad(1)
      brkn_y=1-2*xrad(2)
      xjac=xjac*2
      xjac=xjac*1.5d0*(1-brkn_y**2)
      brkn_y=1.5d0*(brkn_y-brkn_y**3/3)
      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      call br_compcsimax
      brkn_csi=brkn_csitilde*brkn_csimax
      call br_real_phsp_isr_rad
      jac=xjac*brkn_jacreal*brkn_csimax
      end

      subroutine br_compcsimax
      implicit none
      include 'brinclude.h'
      real * 8 y,xb1,xb2
      xb1=brkn_xb1
      xb2=brkn_xb2
      y=brkn_y
      brkn_csimax=1-max(2*(1+y)*xb1**2/
     1    (sqrt((1+xb1**2)**2*(1-y)**2+16*y*xb1**2)+(1-y)*(1-xb1**2)),
     1            2*(1-y)*xb2**2/
     1    (sqrt((1+xb2**2)**2*(1+y)**2-16*y*xb2**2)+(1+y)*(1-xb2**2)))
      end

      subroutine br_real_phsp_isr_rad
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      include 'pwhg_kn.h'
      real * 8 y,xb1,xb2,x1,x2,betal,betat,vecl(3),vect(3),
     1         cth,sth,cph,sph,csi,pt2
      integer i,mu
      real * 8 dotp
      external dotp
c the following call sets brkn_csimax, brkn_csimaxp, brkn_csimaxm
c also when br_real_phsp_isr_rad is called directly
c (i.e. not through br_real_phsp_isr_rad0)
      call br_compcsimax
      y=brkn_y
      xb1=brkn_xb1
      xb2=brkn_xb2
      csi=brkn_csi
      cth=y
      sth=sqrt(1-cth**2)
      cph=cos(brkn_azi)
      sph=sin(brkn_azi)
      x1=xb1/sqrt(1-csi)*sqrt((2-csi*(1-y))/(2-csi*(1+y)))
      x2=xb2/sqrt(1-csi)*sqrt((2-csi*(1+y))/(2-csi*(1-y)))
      brkn_x1=x1
      brkn_x2=x2
      do mu=0,3
         brkn_preal(mu,1)=kn_beams(mu,1)*x1
         brkn_preal(mu,2)=kn_beams(mu,2)*x2
      enddo
      brkn_sreal=brkn_sborn/(1-csi)
c Build k_n+1 in the rest frame of brkn_preal
c      write(*,*) ' br_nlegreal ',br_nlegreal
      brkn_preal(0,br_nlegreal)=sqrt(brkn_sreal)*csi/2
      brkn_preal(1,br_nlegreal)=brkn_preal(0,br_nlegreal)*sth*sph
      brkn_preal(2,br_nlegreal)=brkn_preal(0,br_nlegreal)*sth*cph
      brkn_preal(3,br_nlegreal)=brkn_preal(0,br_nlegreal)*cth
c boost it to the frame of brkn_preal
      do i=1,3
         vecl(i)=(brkn_preal(i,1)+brkn_preal(i,2))
     1          /(brkn_preal(0,1)+brkn_preal(0,2))
      enddo      
      betal=sqrt(vecl(1)**2+vecl(2)**2+vecl(3)**2)
      if(betal.gt.0) then
         do i=1,3
            vecl(i)=vecl(i)/betal
         enddo
      else
         vecl(1)=1
         vecl(2)=0
         vecl(3)=0
      endif
      call mboost(1,vecl,betal,
     1    brkn_preal(0,br_nlegreal),brkn_preal(0,br_nlegreal))
c longitudinal boost of underlying Born to zero rapidity frame
      do i=1,3
         vecl(i)=(brkn_pborn(i,1)+brkn_pborn(i,2))
     1          /(brkn_pborn(0,1)+brkn_pborn(0,2))
      enddo
      betal=sqrt(vecl(1)**2+vecl(2)**2+vecl(3)**2)
      if(betal.gt.0) then
         do i=1,3
            vecl(i)=vecl(i)/betal
         enddo
      else
         vecl(1)=1
         vecl(2)=0
         vecl(3)=0
      endif
      call mboost(br_nlegborn-2,vecl,-betal,
     1 brkn_pborn(0,3),brkn_preal(0,3))
c      call printtot(br_nlegborn,brkn_preal(0,1))
c construct transverse boost velocity
      vect(3)=0
      vect(1)=brkn_preal(1,br_nlegreal)
      vect(2)=brkn_preal(2,br_nlegreal)
      pt2=vect(1)**2+vect(2)**2
c      betat=1/sqrt(1+(brkn_sreal*(1-csi))/pt2)
      betat=sqrt(pt2/(pt2+brkn_sreal*(1-csi)))
      if(pt2.eq.0) then
         vect(1)=1
         vect(2)=0
      else
         vect(1)=vect(1)/sqrt(pt2)
         vect(2)=vect(2)/sqrt(pt2)
      endif
c     write(*,*) ' k+1: ',(brkn_preal(mu,br_nlegreal),mu=0,3)
         call mboost(br_nlegborn-2,vect,-betat,
     1        brkn_preal(0,3),brkn_preal(0,3))
c      call printtot(nlegborn,brkn_preal(0,1))
c longitudinal boost in opposite direction
      call mboost(br_nlegborn-2,vecl,betal,
     1 brkn_preal(0,3),brkn_preal(0,3))
c      call printtot(br_nlegreal,brkn_preal(0,1))
      brkn_jacreal=brkn_sreal/(4*pi)**3*csi/(1-csi)
      call br_compcmkin
      call br_compdij
      end


      subroutine br_compcmkin
      implicit none
      include 'brinclude.h'
      real * 8 vecl(3),betal
      data vecl/0d0,0d0,1d0/
      save vecl
      betal=-(brkn_preal(3,1)+brkn_preal(3,2))/
     1 (brkn_preal(0,1)+brkn_preal(0,2))
      call mboost(br_nlegreal,vecl,betal,brkn_preal,brkn_cmpreal)
      end

      subroutine br_compdij
      implicit none
      include 'brinclude.h'
      integer j,k
      real * 8 y
      real * 8 crossp,dotp
      external crossp,dotp
      do j=3,br_nlegreal
         y=1-dotp(brkn_cmpreal(0,1),brkn_cmpreal(0,j))
     1 /(brkn_cmpreal(0,1)*brkn_cmpreal(0,j))
         brkn_dijterm(0,j)=(brkn_cmpreal(0,j)**2
     1 *(1-y**2))**brpar_diexp
         brkn_dijterm(1,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1-y))**brpar_diexp
         brkn_dijterm(2,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1+y))**brpar_diexp
      enddo
      do j=3,br_nlegreal
         do k=j+1,br_nlegreal
            brkn_dijterm(j,k)=
     1(2*dotp(brkn_cmpreal(0,k),brkn_cmpreal(0,j))*
     1       brkn_cmpreal(0,k)*brkn_cmpreal(0,j)
     2    /  (brkn_cmpreal(0,k)+brkn_cmpreal(0,j))**2)**brpar_dijexp
c     2    /  ((brkn_cmpreal(1,k)+brkn_cmpreal(1,j))**2+
c     3        (brkn_cmpreal(2,k)+brkn_cmpreal(2,j))**2+
c     4        (brkn_cmpreal(3,k)+brkn_cmpreal(3,j))**2))**brpar_dijexp
         enddo
      enddo
      end

      subroutine br_compute_csimax_fsr
      implicit none
c Compute csimax for all possible final state emitters;
c for initial state emitters it is not possible, since
c csimax depends upon y in this case.
      include 'brinclude.h'
      integer j
      real * 8 q0,mrec2
      logical valid_emitter
      external valid_emitter
c For later convenience we calculate csimax for all possible final
c state emitters, due at the present we only use particle n. 3 as
c an emitter and obtain all the other regions by momentum swapping:
      do j=3,br_nlegborn
        q0=2*brkn_cmpborn(0,1)
        mrec2=(q0-brkn_cmpborn(0,j))**2
     1       -brkn_cmpborn(1,j)**2-brkn_cmpborn(2,j)**2
     1       -brkn_cmpborn(3,j)**2
        brkn_csimax_arr(j)=1-mrec2/q0**2
      end do
      end



      subroutine born_suppression(fact)
      implicit none
      real * 8 fact
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      real * 8 bornsuppfact,powheginput,dotp,scale,
     1     pt1sq,pt2sq,pt3sq,msq12,msq23,msq31
      logical ini
      data ini/.true./
      external powheginput,dotp
      integer i,j
      real * 8 ptsqrel
      real*8 pt1,pt2,pt3,Ht
      real * 8 p,pp,rapsuppfact,y1,y2,y3
      parameter (p=2d0)
      parameter (pp=1d0)
      real * 8 pwhg_rapidity
      external pwhg_rapidity
      real * 8 scale2,scale3
      save ini,scale,scale2,scale3,rapsuppfact
      ptsqrel(i,j)=2*dotp(kn_cmpborn(0,i),kn_cmpborn(0,j))
     1     *  kn_cmpborn(0,i)*kn_cmpborn(0,j)
     2     /(kn_cmpborn(0,i)**2+kn_cmpborn(0,j)**2)
      if (ini) then
         if (.not.flg_minlo) then

            scale2 = powheginput("#suppnominlo")
            if(scale2.lt.0) scale2 = 0
            scale2 = scale2**2        
            
            if(scale2.eq.0) then
               scale3 = powheginput("#otherscale")
               if(scale3.lt.0) scale3 = 0
               if(scale3.eq.0) then
                  write(*,*) ' you should specify either '//
     1                 'suppnominlo or otherscale'//
     2                 ' in the powheg.input file!'
                  call pwhg_exit(-1)
               endif
               scale3 = scale3**2
            else
               scale3 = 0
            endif
         endif
         scale  = powheginput("bornsuppfact")**2
         rapsuppfact = powheginput("#rapsuppfact")
         if (rapsuppfact.lt.0d0) rapsuppfact=0d0
         ini = .false.
      endif
      pt1sq = kn_cmpborn(1,3)**2+kn_cmpborn(2,3)**2
      pt2sq = kn_cmpborn(1,4)**2+kn_cmpborn(2,4)**2
      pt3sq = kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2
      msq12 = ptsqrel(3,4)
      msq23 = ptsqrel(4,5)
      msq31 = ptsqrel(5,3)
      if (.not.flg_minlo) then
         if(scale2.gt.0) then
            fact = exp(-scale2**p * (
     1           1/pt1sq**p
     2           +1/pt2sq**p
     3           +1/pt3sq**p
     4           +1/msq12**p
     5           +1/msq23**p
     6           +1/msq31**p ))
            pt1 = sqrt(pt1sq)
            pt2 = sqrt(pt2sq)
            pt3 = sqrt(pt3sq)
            Ht  = pt1 + pt2 + pt3
            fact = fact/scale**p/(1d0/scale + 1d0/Ht**2)**p
         else
            fact = 1/scale3**pp/(1/scale3+1/pt1sq)**pp
     1           * 1/scale3**pp/(1/scale3+1/pt2sq)**pp
     2           * 1/scale3**pp/(1/scale3+1/pt3sq)**pp
     3           * 1/scale3**pp/(1/scale3+1/msq12)**pp
     4           * 1/scale3**pp/(1/scale3+1/msq23)**pp
     5           * 1/scale3**pp/(1/scale3+1/msq31)**pp
            pt1 = sqrt(pt1sq)
            pt2 = sqrt(pt2sq)
            pt3 = sqrt(pt3sq)
            Ht  = pt1 + pt2 + pt3
            fact = fact/scale**p/(1d0/scale + 1d0/Ht**2)**p
         endif
      else
c     In the case of MINLO evaluation sum of pt's are used:
         pt1 = sqrt(pt1sq)
         pt2 = sqrt(pt2sq)
         pt3 = sqrt(pt3sq)
         Ht  = pt1 + pt2 + pt3
         fact = 1d0/scale**p/(1d0/scale + 1d0/Ht**2)**p
      endif 
c     suppress regions with high rapidities
      if (rapsuppfact.ne.0) then
         y1=pwhg_rapidity(kn_cmpborn(0,3))
         y2=pwhg_rapidity(kn_cmpborn(0,4))
         y3=pwhg_rapidity(kn_cmpborn(0,5))
         fact=fact*exp(-(y1**2+y2**2+y3**2)/rapsuppfact**2)
      endif
      end
      

      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'PhysPars.h'
      logical ini,fixedscale
      data ini/.true./
      real * 8 muf,mur
      real * 8 muref
      real * 8 pt1,pt2,pt3,pt4
      real * 8 Ht
      real * 8 powheginput
      external powheginput
      save ini,fixedscale,muref
      if (ini) then
         if(flg_minlo) then
            write(*,*) '****************************************'
            write(*,*) '*******          MINLO ACTIVE    *******'
            write(*,*) '****************************************'
            write(*,*) '*******     FIXED SCALES!          *****'
            fixedscale=.true.
            muref=ph_Zmass
         else        
            muref=powheginput("#fixedscale")                        
            print *,"********************************************"
            print *,"*                                          *"
            if(muref.lt.0) then
               print *,"*         Dynamical scale is used          *"
               fixedscale=.false.
            else
               print *,"*       Fixed scale is used with: ",muref
               fixedscale=.true.
            endif
            print *,"*                                          *"
            print *,"********************************************"
         endif
         ini=.false.
      endif
      if (.not.fixedscale) then
         if ((flg_btildepart.eq.'b').or.(flg_btildepart.eq.'c')) then
            pt1 = sqrt(kn_cmpborn(1,3)**2+kn_cmpborn(2,3)**2)
            pt2 = sqrt(kn_cmpborn(1,4)**2+kn_cmpborn(2,4)**2)
            pt3 = sqrt(kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2)
            Ht = pt1 + pt2 + pt3
         elseif ((flg_btildepart.eq.'r')) then
            pt1 = sqrt(kn_preal(1,3)**2+kn_preal(2,3)**2)
            pt2 = sqrt(kn_preal(1,4)**2+kn_preal(2,4)**2)
            pt3 = sqrt(kn_preal(1,5)**2+kn_preal(2,5)**2)
            pt4 = sqrt(kn_preal(1,6)**2+kn_preal(2,6)**2)
            Ht = pt1 + pt2 + pt3 + pt4
         else
            print *,"Problem occured in set_fac_ren_scales"
            print *,"flg_btildepart: ",flg_btildepart
            call pwhg_exit(-1)
         end if
c     We use the Ht/2:
         muref = Ht/2
      endif
c freeze alpha to avoid NaN's from alpha_s
      muref = max(muref,0.4d0)
      muf=muref
      mur=muref
      end


      subroutine swap_momenta(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3),tmp(0:3)
      tmp=p1
      p1=p2
      p2=tmp
      end
