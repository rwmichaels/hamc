        double precision function cross_section(Z,A,e,th,ef)

        implicit none

        double precision Z,A,alpha,e,ef
        double precision mn,mp,mott,nbarn,nu
        double precision s2,t2,theta,th
        double precision xbj,qmu2,wmm
        double precision f1psfun,f2psfun,f1nsfun,f2nsfun
        double precision f1dsfun,f2dsfun,f1hesfun,f2hesfun
        double precision f1,f2,w1,w2
c
        if (th.le.0) then
           theta=-th
        else
           theta=th
        endif
        mp=0.9382727
        mn=0.9395653

        alpha=7.3e-3 !  alpha = 1/137
        nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn
c
        s2=sin(abs(theta)/2.)**2

c       Calculate kinematic variables
        nu = e - ef             !energy transfer
        qmu2 = 4*e*ef*s2        !Q2
        xbj=qmu2/2/nu/mp        !Bjorken x

        wmm=mp*mp+2*mp*nu-qmu2  !invariant mass

c       Calculate Mott cross section in nbarn/(sr)
        mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn
        if (A==1.) then    ! for proton
           f1 = f1psfun(xbj,qmu2)
           f2 = f2psfun(xbj,qmu2)
        else if (A==2) then  ! for deuteron
           f1 = 2.*f1dsfun(xbj,qmu2)
           f2 = 2.*f2dsfun(xbj,qmu2)
c           print *,'for deuteron',f1,f2
        else if (A==3) then     ! for 3He.
           f1 = 3.*f1hesfun(xbj,qmu2)
           f2 = 3.*f2hesfun(xbj,qmu2)
        else                    ! for heavy nuclei
           f1 = Z*f1psfun(xbj,qmu2)+(A-Z)*f1nsfun(xbj,qmu2)
           f2 = Z*f2psfun(xbj,qmu2)+(A-Z)*f2nsfun(xbj,qmu2)
        endif
c        print *,'in NMC: Z=',Z,' A=',A,' f1=',f1,' f2=',f2

        w1 = f1/mp
        w2 = f2/nu
c
        t2 = tan(abs(theta)/2.)**2
c
c       DIS cross section in nbarn/sr-GeV
c
        cross_section = mott * (w2 + 2 * w1 * t2)
c        print *,'mott=',mott,' f1=',f1,' f2=',f2,' w1=',w1,' w2=',w2,
c     > ' xsec=', cross_section
 103    continue
        end

********************** f1psfun ***********************************
* nucleon averaged structure functions

      double precision function f1psfun(aks,t)
      implicit double precision(a-h,o-z)
      f2p=d95f2h8(t,aks)
      amp=0.9382727
      dnu=t/(2*amp*aks)
      f1psfun=f2p*(1+t/dnu**2)/(2.*aks*(1.+r1998(aks,t)))
      end


********************** f2psfun ***********************************
* nucleon averaged structure functions

      double precision function f2psfun(aks,t)
      implicit double precision(a-h,o-z)
      f2psfun=d95f2h8(t,aks)
      end


********************** f1nsfun ***********************************
      double precision function f1nsfun(aks,t)
      implicit double precision(a-h,o-z)
      df1d=f1dsfun(aks,t)
      df1p=f1psfun(aks,t)
      f1nsfun=2*df1d-df1p
      end



********************** f2nsfun ***********************************
      double precision function f2nsfun(aks,t)
      implicit double precision(a-h,o-z)
      df2d=f2dsfun(aks,t)
      df2p=f2psfun(aks,t)
      f2nsfun=2*df2d-df2p
      end


********************** f1dsfun ***********************************
* note: per nucleon definition
*
      double precision function f1dsfun(aks,t)
      implicit double precision(a-h,o-z)
      f2=d95f2d8(t,aks)
      amp=0.9382727
      dnu=t/(2*amp*aks)
      f1dsfun=f2*(1+t/dnu**2)/(2.*aks*(1.+r1998(aks,t)))
c      print *,'f1dsfun:',f1dsfun,f2,r1998(aks,t)
      end

********************** f2dsfun ***********************************
* note: per nucleon definition
*
      double precision function f2dsfun(aks,t)
      implicit double precision(a-h,o-z)
      f2dsfun=d95f2d8(t,aks)
      end


********************** f1hesfun ***********************************
      double precision function f1hesfun(aks,t)
      implicit double precision(a-h,o-z)
      f2=f2hesfun(aks,t)
      amp=0.9382727
      dnu=t/(2*amp*aks)
      f1hesfun=f2*(1+t/dnu**2)/(2.*aks*(1.+r1998(aks,t)))
c      print *,'nmc95:',aks,t,f1hesfun
      end

********************** f2hesfun ***********************************
      double precision function f2hesfun(aks,t)
      implicit double precision(a-h,o-z)
      f2p=d95f2h8(t,aks)
      f2d=d95f2d8(t,aks)
      f2hesfun=(f2p+2d0*f2d)/3d0
c      print *,'nmc95:',f2hesfun
      end

********************** d95f2d8  ***********************************
      double precision function d95f2d8(dq2,dx)
*:=====================================================================:
*:                                                                     :
*:      author:    x.zheng        last update: 11.25.2001              :
*:                                 tested: xxx                         :
*:                                                                     :
*:      arguments: dq2,dx: double prec. input xbj,q2                   :
*:                 d95f2h8* double prec f2  output                     :
*:                                                                     :
*:      called by: mkf2                                                :
*:                                                                     :
*:      action:    calculate f2 structure function of the deuteron     :
*:                 nmc fit of dis-region with 15 parameters fit        :
*:                 kinematics range: 0.5<q2<75 GeV^2, 0.006<x<0.9      :
*:                                                                     :
*:                 parametrized with a1-7,b1-4,c1-4  as                :
*:                                                                     :
*:                 f2_dis(x,q2) ~prop.                                 :
*:                   A(x)*(ln(q2/dl2)/ln(q20/dl2))**(B(x))*(1+C(x)/q2) :
*:                   with x = (q2+m_a)/(2m*nu + m_b**2)                :
*:                        dl2 = (0.250 GeV)^2, (so-called lambda^2)    :
*:                        q20 = 20 GeV^2                               :
*:                        A(x) = x**a1*(1-x)**a2*(a3+a4*(1-x)          :
*:                               +a5*(1-x)**2+a6*(1-x)**3+a7*(1-x)**4) :
*:                        B(x) = b1+b2*x+b3/(x+b4)                     :
*:                        C(x) = c1*x+c2*x**2+c3*x**3+c4*x**4          :
*:                 reference:                                          :
*:                 the new muon collaboration                          :
*:                 Phys.Lett.B364(1995)107~115                         :
*:                 hep-ph/9509406                                      :
*:                                                                     :
*:                 resonance contribution is calculated in the same    :
*:                 way as df2d8(dq2,dx)                                :
*:                                                                     :
*:       comments: The proton and deuteron structure functions F2p and :
*:                 F2d were measured in the kinematic range 0.006<x<0.6:
*:                 and 0.5<q2<75 GeV^2, by inclusive deep inelastic    :
*:                 muon scattering at 90,120,300 and 280 GeV.  The     :
*:                 measurements are in good agreement with earlier high:
*:                 precision results.  The present and earlier results :
*:                 together have been parametrised to give descriptions:
*:                 of the proton and deuteron structure functions F2   :
*:                 and their uncertainties over the range 0.006<x<0.9  :
*:=====================================================================:
c
      implicit double precision (d)
c
c
c *** a1,..a7,b1..b4,c1..c4 = 15 param of nmc, slac, bcdms (95)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c            values: daw=1.512(gev2), dbw=0.351(gev2)
c            ref:  bodek et al., p.r.d20(1979)1427.
c            see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.2**2 = 0.04 (gev2)
c *** q0**2 = 2 gev2 ... (2+0.351)/0.04 = 58.771
c *** fit by y.m.(25-nov-88 19h43m14s)
c
      data a1,a2,a3,a4,a5,a6
     :     ,a7,b1,b2,b3,b4
     :     ,c1,c2,c3,c4
     :     ,d9,d10
     :     ,daw,dbw
c
c     f2 from nmc phys.lett.b364(1995)107
     :     /-0.04858,2.863,0.8367,-2.532,9.145,-12.504
     :     ,5.473,-0.008,-2.227,0.0551,0.0570
     :     ,-1.509,8.553,-31.20,39.98
c     resonance-region
     :     ,.89456,.16452
     :     ,1.512,.351 /
c
c
      d95f2d8=1.d-30
      dl2 = 0.25**2
      dq20 = 20
      d95f2d8 = (dx**a1*(1-dx)**a2*(a3+a4*(1-dx)
     :      +a5*(1-dx)**2+a6*(1-dx)**3+a7*(1-dx)**4))
     :      *(log(dq2/dl2)/log(dq20/dl2))
     :      **(b1+b2*dx+b3/(dx+b4))
     :      *(1+(c1*dx+c2*dx**2+c3*dx**3+c4*dx**4)/dq2)
C
      if (d95f2d8.gt.0d0) return
      d95f2d8 = 1.d-30
      return
      end


********************** d95f2h8  ***********************************
      double precision function d95f2h8(dq2,dx)
*:=====================================================================:
*:                                                                     :
*:      author:    x.zheng        last update: 11.25.2001              :
*:                                tested: xxx                          :
*:                                                                     :
*:      arguments: dq2,dx: double prec. input xbj,q2                   :
*:                 d95f2h8* double prec f2  output                   :
*:                                                                     :
*:      called by: mkf2                                                :
*:                                                                     :
*:      action:    calculate f2 structure function of the proton       :
*:                 nmc fit of dis-region with 15 parameters fit        :
*:                 kinematics range: 0.5<q2<75 GeV^2, 0.006<x<0.9      :
*:                                                                     :
*:                 parametrized with a1-7,b1-4,c1-4  as                :
*:                                                                     :
*:                 f2_dis(x,q2) ~prop.                                 :
*:                   A(x)*(ln(q2/dl2)/ln(q20/dl2))**(B(x))*(1+C(x)/q2) :
*:                   with x = (q2+m_a)/(2m*nu + m_b**2)                :
*:                        dl2 = (0.250 GeV)^2, (so-called lambda^2)    :
*:                        q20 = 20 GeV^2                               :
*:                        A(x) = x**a1*(1-x)**a2*(a3+a4*(1-x)          :
*:                               +a5*(1-x)**2+a6*(1-x)**3+a7*(1-x)**4) :
*:                        B(x) = b1+b2*x+b3/(x+b4)                     :
*:                        C(x) = c1*x+c2*x**2+c3*x**3+c4*x**4          :
*:                 reference:                                          :
*:                 the new muon collaboration                          :
*:                 Phys.Lett.B364(1995)107~115                         :
*:                 hep-ph/9509406                                      :
*:                                                                     :
*:                 resonance contribution is calculated in the same    :
*:                 way as df2d8(dq2,dx)                                :
*:                                                                     :
*:       comments: The proton and deuteron structure functions F2p and :
*:                 F2d were measured in the kinematic range 0.006<x<0.6:
*:                 and 0.5<q2<75 GeV^2, by inclusive deep inelastic    :
*:                 muon scattering at 90,120,300 and 280 GeV.  The     :
*:                 measurements are in good agreement with earlier high:
*:                 precision results.  The present and earlier results :
*:                 together have been parametrised to give descriptions:
*:                 of the proton and deuteron structure functions F2   :
*:                 and their uncertainties over the range 0.006<x<0.9  :
*:=====================================================================:
c
      implicit double precision (d)
*
c *** a1,..a7,b1..b4,c1..c4 = 15 param of nmc, slac, bcdms (95)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c            values: daw=1.512(gev2), dbw=0.351(gev2)
c            ref:  bodek et al., p.r.d20(1979)1427.
c            see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.25**2 = 0.0625 (gev2)
c *** q0**2 = 20 gev2
*
      data
     .amm/2.7928456d0/,amn/-1.913148d0/,pi/3.1415926d0/
     .,alfa/.729735d-2/,amh/.938272d0/,ampi/.104/
c
      data a1,a2,a3,a4,a5,a6
     :     ,a7,b1,b2,b3,b4
     :     ,c1,c2,c3,c4
     :     ,d9,d10,d11,d12,d13,d14
     :     ,d15,d16
     :     ,daw,dbw
c
c     f2 from nmc phys.lett.b364(1995)107
     :     /-.02778,2.926,1.0362,-1.840,8.123,-13.074
     :     ,6.215,0.285,-2.694,0.0188,0.0274
     :     ,-1.413,9.366,-37.79,47.10
c     resonance-region:
     :     ,.1179, .044735, .038445, .27921, 8.8228d-5, 6.2099d-5
     :     ,1.421,1.2582
     :     ,1.642, .376/
c
      d95f2h8=1.d-30
      dl2 = 0.25**2
      dq20 = 20.
c
      d95f2h8=1.d-30
      dl2 = 0.25**2
      dq20 = 20.
      d95f2h8 = (dx**a1*(1-dx)**a2*(a3+a4*(1-dx)
     :      +a5*(1-dx)**2+a6*(1-dx)**3+a7*(1-dx)**4))
     :      *(log(dq2/dl2)/log(dq20/dl2))
     :      **(b1+b2*dx+b3/(dx+b4))
     :      *(1+(c1*dx+c2*dx**2+c3*dx**3+c4*dx**4)/dq2)
c
      dres = 0.d0
c
c *** (total) = (qcd part) + (resonance region)
c
      d95f2h8 = d95f2h8 + dres
c
      if(d95f2h8 .gt. 0.d0) return
      d95f2h8=1.d-30
c
      return
      end
c
