
********************** r1998 ************************************
      double precision function r1998(aks,tc)
      implicit double precision(a-h,o-z)
c
c the following part is r1998, Dec. 2001
c
c this is r1998 from E143
c ref.: hep-ex/9808028, Phys.Lett.B452(1999)194~200
c teta and zn are the same as r1990
c   for r1990, q2 set to >=0.35
c   for r1998, fits work for 0.005<=x<=0.86, 0.5<=q2<=130
c              here set q2 >= 0.5
c coefficients for r1998:
c fit   1      2      3       4        5       6
c  a 0.0485 0.5470 2.0621  -0.3804  0.5090  -0.0285 
c  b 0.0481 0.6114 -0.3509 -0.4611  0.7172  -0.0317
c  c 0.0577 0.4644 1.8288  12.3708 -43.1043 41.7415
c --------------------------------------------------
      t=tc
      if(tc.lt..5d0)t=0.50
      teta=1.+12.*t/(1.+t)*(.125**2/(aks**2+.125**2))
      zn=teta/log(t/.04)
      ra=.0485*zn+(.5470/(t**4+2.0621**4)**(.25d0))
     .    *(1-0.3804*aks+0.5090*aks**2)*aks**(-0.0285)
      rb=.0481*zn+(.6114/t-.3509/(t**2+.09))
     .    *(1-0.4611*aks+0.7172*aks**2)*aks**(-0.0317)
      rc=.0577*zn+.4644/sqrt((t-(12.3708*aks-43.1043
     .    *aks**2+41.7415*aks**3))**2+1.8288**2)
c 513  format('r1998: '5(f7.4))
c      print 513,aks,t,ra,rb,rc
      rrr=(ra+rb+rc)/3.
      r1998=rrr
c
c end of r1998
c
      return
      end

      double precision function dr1998(aks,tc)
      implicit double precision(a-h,o-z)
      x=aks
      Q2=tc
      dr1998=0.0078-0.013*x+(0.070-0.39*x+0.70*x*x)/(1.7+Q2);
      return
      end

