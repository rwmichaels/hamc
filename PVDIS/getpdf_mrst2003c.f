      SUBROUTINE GETPDF_MRST2003C(Eb,xb,q2,Ad)
C
c template code to print the MRST2004c PDFs
C
      implicit double precision (a-h,o-z)

      double precision r1998
      integer mode  
      CHARACTER*50 filename
      CHARACTER*80 buf
      logical STAT ! .false. = table created
                   ! .true. = value found
      logical MRSTDEBUG
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/ ! MRST limit

      DATA PI/3.1416/
      DATA SWnominal/0.235/
      DATA AKAPPA/539.4/         ! 3*GF/2/sqrt(2)/pi/alpha, in ppm GeV-2

      MRSTDEBUG=.false.
c      MRSTDEBUG=.true.
C
      STAT=.true.
C
C
C     clear all outputs
C
      uv=0
      dv=0
      u=0
      d=0
      
      s=0
      c=0
      ubar=0
      dbar=0
      Rs=0
      Rv=0
      Rc=0
      
      Ad=0
      Ap=0
      a2=0
      b2=0
      c2=0
      
C
C     now start the program
C
c      if (MRSTDEBUG) print *,'in MRST2003c',xb,q2,Eb

      if(q2.lt.qsqmin.or.q2.gt.qsqmax) then
c         print 99,q2
         return
      endif
      if(xb.lt.xmin.or.xb.gt.xmax) then
c         print 98,xb
         return
      endif
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
C
C check if table is available
C

      if (abs(Eb-2.2*nint(Eb/2.2)).lt.0.1) then
         Ebmin=2.2
         Ebmax=11.
         Ebin=2.2
         nEbin=(Ebmax-Ebmin)/Ebin+1
         iEb=nint((Eb-Ebmin)/Ebin)+1
         if (mode.eq.1) then
            WRITE(filename,
     >           '("pdffiles_11gev/mrstPDF2003c_nlo.dat.",I1)') iEb
         else
         WRITE(filename,
     >           '("pdffiles_11gev/mrstPDF2003c_nnlo.dat.",I1)') iEb
         endif
      elseif (abs(Eb-1.2*nint(Eb/1.2)).lt.0.1) then
         Ebmin=1.2
         Ebmax=6.0
         Ebin=1.2
         nEbin=(Ebmax-Ebmin)/Ebin+1
         iEb=nint((Eb-Ebmin)/Ebin)+1
         if (mode.eq.1) then
            WRITE(filename,
     >           '("pdffiles_6gev/mrstPDF2003c_nlo.dat.",I1)') iEb
         else
         WRITE(filename,
     >           '("pdffiles_6gev/mrstPDF2003c_nnlo.dat.",I1)') iEb
         endif
      else
         CALL CRPDF_MRST2003c_SINGLE(Eb,xb,q2,mode,uv,dv,u,d,s,c,
     >        ubar,dbar,Rs,Rv,Rc,Ad,Ap,a2,b2,c2,STAT)
         return
      endif
c      if (MRSTDEBUG) PRINT *,'open ',filename
      
      itemp=NextUn()
      OPEN(itemp,FILE=filename,STATUS='OLD',err=101)
C
c      if (MRSTDEBUG) PRINT *,'read MRST2003c from ',filename,xb,q2
C
      CLOSE(itemp)
      CALL READPDF_SINGLE(filename,xb,q2,uv,dv,u,d,s,c,ubar,dbar,
     >Rs,Rv,Rc,Ad,Ap,a2,b2,c2)


      dnu=q2/2/0.938/xb
      y=dnu/Eb
      dR=r1998(xb,q2)
c      if (MRSTDEBUG) print *,'MRST2003c check dR=',dR,r1998(xb,q2)
      dY=(1.-(1.-y)*(1.-y))
     >     /(1.+(1.-y)*(1.-y)-y*y*dR/(1.+dR))

      C1u=-0.5+4./3.*SWnominal
      C1d=0.5-2./3.*SWnominal
      C2u=-0.5+2.*SWnominal
      C2d=0.5-2.*SWnominal

      a2_new=1/(AKAPPA*q2)*(5+Rs+4*Rc)/(dY*Rv)
      b2_new=-(2*C1u*(1+Rc)-C1d*(1+Rs))/(dY*Rv)

      if (abs(a2_new/a2-1).gt.0.04) then
c         print *,'MRST2003c check a2:',a2_new,a2
      endif
      if (abs(b2_new/b2-1).gt.0.04) then
c         print *,'MRST2003c check b2:',b2_new,b2
      endif

      RETURN
C
C     calculate MRST PDF and save to table
C
 101  if (MRSTDEBUG) print *,'creat MRST2003c table'
      CALL CRPDF_MRST2003C(Eb)
      STAT=.false.

      return
      end

      SUBROUTINE CRPDF_MRST2003C(Eb_i)

      implicit double precision (a-h,o-z)
c      implicit real*8 (a-h,o-z)

      real*8 xb_r,q_r
      real*8 pdfs(9),dpdfs(9,15),delpdf(9)
      real*8 upv1,dnv1,usea1,dsea1,str1,chm1,bot1,glu1,phot1
      real*8 upv2,dnv2,usea2,dsea2,str2,chm2,bot2,glu2,phot2
      double precision r1998
      CHARACTER*50 filename

      parameter (NEB=10)
      dimension iutil(NEB)
      dimension vAd(NEB),vAp(NEB),vdAd(NEB),vdAp(NEB)
      dimension vAd1(NEB),vAp1(NEB),vAd2(NEB),vAp2(NEB)
      dimension vAd0(NEB),vAp0(NEB)
      CHARACTER*80 buf1
      dimension va2(NEB),vb2(NEB),vc2(NEB),vda2(NEB),vdb2(NEB),vdc2(NEB)
      dimension va21(NEB),vb21(NEB),vc21(NEB),
     >     va22(NEB),vb22(NEB),vc22(NEB),va20(NEB),vb20(NEB),vc20(NEB)

      integer mode

      logical MRSTDEBUG1

      DATA PI/3.1416/
      DATA SWnominal/0.235/
      DATA AKAPPA/539.4/         ! 3*GF/2/sqrt(2)/pi/alpha, in ppm GeV-2


      MRSTDEBUG1=.false.
c      MRSTDEBUG1=.true.


      xbmin=0.01
      xbmax=0.9
      nxbin=90
      xbin=(xbmax-xbmin)/(nxbin-1)
      q2min=1.25
c      q2max=10
c      q2bin=0.25
      q2max=12.
      q2bin=0.1
      nq2bin=(q2max-q2min)/q2bin+1

      C1u=-0.5+4./3.*SWnominal
      C1d=0.5-2./3.*SWnominal
      C2u=-0.5+2.*SWnominal
      C2d=0.5-2.*SWnominal
      
      Eb=Eb_i
      if (abs(Eb-2.2*nint(Eb/2.2)).lt.0.1) then
         Ebmin=2.2
         Ebmax=11.
         Ebin=2.2
         nEbin=(Ebmax-Ebmin)/Ebin+1
      else
         Ebmin=1.2
         Ebmax=6.0
         Ebin=1.2
         nEbin=(Ebmax-Ebmin)/Ebin+1
      endif

      do mode=1,2
         do iEb=1,5
            if (abs(Eb-2.2*nint(Eb/2.2)).lt.0.1) then
               if (mode.eq.1) then
                  WRITE(filename,
     >               '("pdffiles_11gev/mrstPDF2003c_nlo.dat.",I1)') iEb
               else
                  WRITE(filename,
     >               '("pdffiles_11gev/mrstPDF2003c_nnlo.dat.",I1)') iEb
               endif
            else
               if (mode.eq.1) then
                  WRITE(filename,
     >                 '("pdffiles_6gev/mrstPDF2003c_nlo.dat.",I1)') iEb
               else
                  WRITE(filename,
     >               '("pdffiles_6gev/mrstPDF2003c_nnlo.dat.",I1)') iEb
               endif
            endif
            iutil(iEb)=NextUn()
            if (MRSTDEBUG1) 
     >           print *,' creating file ',filename,' unit=',iutil(iEb)

            OPEN(iutil(iEb),FILE=filename,STATUS='UNKNOWN')
            WRITE(iutil(iEb),108) xbmin,xbmax,nxbin,xbin
            WRITE(iutil(iEb),109) q2min,q2max,nq2bin,q2bin
            Eb=Ebmin+(iEb-1)*Ebin
            WRITE(iutil(iEb),110) Eb
            
            WRITE(iutil(iEb),117)
         enddo
C
C     now start calculating for a grid of x and Q2
C
         do ix=1,nxbin
            if (MRSTDEBUG1) WRITE(*,116)

            do iq2=1,nq2bin
               xb=xbmin+(ix-1)*xbin
               q2=q2min+(iq2-1)*q2bin            
C     
               q=q2**0.5
               xb_r=xb
               q_r=q
               
               if (MRSTDEBUG1) print *,
     >              'calling mrst2003c, x=',xb_r,' Q=',q_r,' mode=',mode            
               call mrst2003c(xb_r,q_r,mode,upv1,dnv1,usea1,dsea1,
     >              str1,chm1,bot1,glu1,phot1)
               
               upv=upv1
               dnv=dnv1
               usea=usea1
               dsea=dsea1
               str=str1
               chm=chm1
               bot=bot1
               glu=glu1
               phot=phot1
               
               uv=upv/xb
               dv=dnv/xb
               u=(upv+usea)/xb
               d=(dnv+dsea)/xb
               c=chm/xb
               cbar=chm/xb
               s=str/xb
               ubar=usea/xb
               dbar=dsea/xb
C     
               u0=u
               d0=d
               s0=s
               c0=c
               ubar0=ubar
               dbar0=dbar
               uv0=uv
               dv0=dv
C     
               Rc=2*(chm+chm)/(upv+dnv+usea+dsea+usea+dsea)
               Rs=2*(str+str)/(upv+dnv+usea+dsea+usea+dsea)
               Rv=(upv+dnv)/(upv+dnv+usea+dsea+usea+dsea)
C     
               Rc0=Rc
               Rs0=Rs
               Rv0=Rv
C     
               do iEb=1,5
                  Eb=Ebmin+(iEb-1)*Ebin
                  dnu=q2/2./0.938272/xb
                  y=dnu/Eb
                  dR=r1998(xb,q2)
                  dY=(1.-(1.-y)*(1.-y))
     >                 /(1.+(1.-y)*(1.-y)-y*y*dR/(1.+dR)) ! bigY
C     
                  if (MRSTDEBUG1) print *,'check Eb',
     >                 Eb,dnu,y,dR,r1998(xb,q2),dY
                  vAd0(iEb)=(q2*AKAPPA*(2.*C1u*(1.+Rc0)-C1d*(1.+Rs0)
     >                 +dY*(2.*C2u-C2d)*Rv0))/(5.+Rs0+4.*Rc0)
                  vAp0(iEb)=q2*AKAPPA*(2*C1u*(u+ubar)-C1d*(d+dbar+s+s)
     >                 +dY*(2*C2u*uv-C2d*dv))/(4*(u+ubar)+d+dbar+s+s)
                  if (MRSTDEBUG1) print *,'vA0=',vAd0(iEb),vAp0(iEb)
C     
C     calculate coefficient for Ad->(2C2u+C2d)
C     
                  va20(iEb)=1/(AKAPPA*q2)*(5+Rs+4*Rc)/(dY*Rv)
                  vb20(iEb)=-(2*C1u*(1+Rc)-C1d*(1+Rs))/(dY*Rv)   
                  vc20(iEb)=va20(iEb)*vAd0(iEb)+vb20(iEb)
               enddo
               
               do iEb=1,5
                  Eb=Ebmin+(iEb-1)*Ebin
                  sAd=vAd0(iEb)
                  sAp=vAp0(iEb)
                  
                  a2=va20(iEb)
                  b2=vb20(iEb)
                  c2=vc20(iEb)
C     
                  WRITE(iutil(iEb),107) xb,q2,uv,dv,u,d,s,c,
     >                 ubar,dbar,Rs,Rv,Rc,sAd,sAp,a2,b2,c2
                  if (MRSTDEBUG1) WRITE(*,106)
     >                 Eb,xb,q2,uv,dv,u,d,s,c,ubar,dbar,
     >                 Rs,Rv,Rc,sAd,sAp,a2,b2,c2
               enddo
            enddo
         enddo
C
C     close all files
C
         do iEb=1,5
            CLOSE(iutil(iEb))
         enddo
      
      enddo  ! mode=1,2

 108  FORMAT(2(2X,F5.3),1X,I5,1(2X,F5.3),' :  xbj min, max, #bins, bin')
 109  FORMAT(2(1X,F6.3),1X,I5,1(1X,F6.3),' :  q2  min, max, #bins, bin')
 110  FORMAT((1X,F6.3),20X,' :  Eb')

 106  FORMAT(2(1X,F5.2),(1X,F6.2),11(2X,F7.4),"/",2(2X,F10.2),
     >  "/",(2X,F11.7),2(2X,F10.5))
 107  FORMAT((1X,F5.2),(1X,F6.2),11(2X,F7.4),2(2X,F10.2),
     > (2X,F11.7),2(2X,F10.5))

 116  FORMAT("   Eb",4X,"x",6X,"Q2",5X,"uv",7X,"dv",7X,"u"   
     >     ,8X,"d",8X,"s",8X,"c",8X,"ubar",6X,"dbar",6X,"Rs",7X
     >     ,"Rv",7X,"Rc",10X,"Ad",10X,"Ap",8X,"a2",12X,"b2",9X,"c2")
 117  FORMAT("   x",6X,"Q2",5X,"uv",7X,"dv",7X,"u"   
     >     ,8X,"d",8X,"s",8X,"c",8X,"ubar",6X,"dbar",6X,"Rs",7X
     >     ,"Rv",7X,"Rc",10X,"Ad",10X,"Ap",8X,"a2",12X,"b2",9X,"c2")
    
      return
      end



      SUBROUTINE CRPDF_MRST2003C_SINGLE(Eb,xb,q2,mode,uv,dv,u,d,s,c,
     >ubar,dbar,Rs,Rv,Rc,Ad,Ap,a2,b2,c2,STAT)

      implicit double precision (a-h,o-z)
c      implicit real*8 (a-h,o-z)

      real*8 xb_r,q_r
      real*8 pdfs(9),dpdfs(9,15),delpdf(9)
      real*8 upv1,dnv1,usea1,dsea1,str1,chm1,bot1,glu1,phot1
      real*8 upv2,dnv2,usea2,dsea2,str2,chm2,bot2,glu2,phot2
      double precision r1998
      CHARACTER*50 filename

      parameter (NEB=1)
      dimension vAd(NEB),vAp(NEB),vdAd(NEB),vdAp(NEB)
      dimension vAd1(NEB),vAp1(NEB),vAd2(NEB),vAp2(NEB)
      dimension vAd0(NEB),vAp0(NEB)
      CHARACTER*80 buf1
      dimension va2(NEB),vb2(NEB),vc2(NEB),vda2(NEB),vdb2(NEB),vdc2(NEB)
      dimension va21(NEB),vb21(NEB),vc21(NEB),
     >     va22(NEB),vb22(NEB),vc22(NEB),va20(NEB),vb20(NEB),vc20(NEB)
      logical STAT ! not used

      integer mode

      logical MRSTDEBUG1

      DATA PI/3.1416/
      DATA SWnominal/0.235/
      DATA AKAPPA/539.4/         ! 3*GF/2/sqrt(2)/pi/alpha, in ppm GeV-2


      MRSTDEBUG1=.false.
c      MRSTDEBUG1=.true.

      C1u=-0.5+4./3.*SWnominal
      C1d=0.5-2./3.*SWnominal
      C2u=-0.5+2.*SWnominal
      C2d=0.5-2.*SWnominal

C
C     now start calculating for the given x and q2           
C     
      q=q2**0.5
      xb_r=xb
      q_r=q
               
      if (MRSTDEBUG1) print *,
     >     'calling mrst2003c, x=',xb_r,' Q=',q_r,' mode=',mode            
      call mrst2003c(xb_r,q_r,mode,upv1,dnv1,usea1,dsea1,
     >     str1,chm1,bot1,glu1,phot1)
               
      upv=upv1
      dnv=dnv1
      usea=usea1
      dsea=dsea1
      str=str1
      chm=chm1
      bot=bot1
      glu=glu1
      phot=phot1
      
      uv=upv/xb
      dv=dnv/xb
      u=(upv+usea)/xb
      d=(dnv+dsea)/xb
      c=chm/xb
      cbar=chm/xb
      s=str/xb
      ubar=usea/xb
      dbar=dsea/xb
C     
      u0=u
      d0=d
      s0=s
      c0=c
      ubar0=ubar
      dbar0=dbar
      uv0=uv
      dv0=dv
C     
      Rc=2*(chm+chm)/(upv+dnv+usea+dsea+usea+dsea)
      Rs=2*(str+str)/(upv+dnv+usea+dsea+usea+dsea)
      Rv=(upv+dnv)/(upv+dnv+usea+dsea+usea+dsea)
C     
      Rc0=Rc
      Rs0=Rs
      Rv0=Rv
C     
      iEb=1
      dnu=q2/2./0.938272/xb
      y=dnu/Eb
      dR=r1998(xb,q2)
      dY=(1.-(1.-y)*(1.-y))
     >     /(1.+(1.-y)*(1.-y)-y*y*dR/(1.+dR)) ! bigY
C     
      vAd0(iEb)=(q2*AKAPPA*(2.*C1u*(1.+Rc0)-C1d*(1.+Rs0)
     >     +dY*(2.*C2u-C2d)*Rv0))/(5.+Rs0+4.*Rc0)
      vAp0(iEb)=q2*AKAPPA*(2*C1u*(u+ubar)-C1d*(d+dbar+s+s)
     >     +dY*(2*C2u*uv-C2d*dv))/(4*(u+ubar)+d+dbar+s+s)
C     
C     calculate coefficient for Ad->(2C2u+C2d)
C     
      va20(iEb)=1/(AKAPPA*q2)*(5+Rs+4*Rc)/(dY*Rv)
      vb20(iEb)=-(2*C1u*(1+Rc)-C1d*(1+Rs))/(dY*Rv)   
      vc20(iEb)=va20(iEb)*vAd0(iEb)+vb20(iEb)

      Ad=vAd0(iEb)
      Ap=vAp0(iEb)
      
      a2=va20(iEb)
      b2=vb20(iEb)
      c2=vc20(iEb)
C     
      return
      end

