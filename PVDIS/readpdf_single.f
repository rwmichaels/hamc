
      SUBROUTINE READPDF_SINGLE(filename,xb,q2,uv,dv,u,d,
     >s,c,ubar,dbar,Rs,Rv,Rc,Ad,Ap,a2,b2,c2)

c     subroutine to read PDFs from look-up tables without error bars

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*50 filename
      CHARACTER*300 buf
      parameter(VDIM=16)
      dimension data1(VDIM),data2(VDIM),data3(VDIM),data4(VDIM)
      dimension data12(VDIM),data34(VDIM),data0(VDIM)
      logical READDEBUG

c      READDEBUG=.true.
      READDEBUG=.false.

      if (READDEBUG) PRINT *,'read PDF_SINGLE from ',filename,xb,q2
      OPEN(26,FILE=filename,STATUS='OLD',err=101)
      READ(26,108,err=102,end=102) xbmin,xbmax,nxbin,xbin
      READ(26,109,err=102,end=102) q2min,q2max,nq2bin,q2bin
 108  FORMAT(2(2X,F5.3),1X,I5,1(2X,F5.3))
 109  FORMAT(2(1X,F6.3),1X,I5,1(1X,F6.3))
      READ(26,*,err=102,end=102) Eb
      READ(26,*,err=102,end=102) buf
      if(READDEBUG) print *,'xb binning=',xbmin,xbmax,nxbin,xbin
      if(READDEBUG) print *,'q2 binning=',q2min,q2max,nq2bin,q2bin
      if(READDEBUG) print *,'Eb=',Eb
      xb_r=0.
      ixb=(xb-xbmin)/xbin+1
      iq2=(q2-q2min)/q2bin+1
      if (ixb.ne.0.and.iq2.ne.0) then
         do i=1,(ixb-1)*nq2bin+(iq2-1)
            READ (26,*,end=102) buf
         enddo
      endif
      READ(26,*,end=102) xb_r1,q2_r1,data1(1),data1(2),data1(3),
     >     data1(4),data1(5),data1(6),data1(7),data1(8),data1(9),
     >     data1(10),data1(11),data1(12),data1(13),data1(14),data1(15),
     >     data1(16)
      READ(26,*,end=102) xb_r2,q2_r2,data2(1),data2(2),data2(3),
     >     data2(4),data2(5),data2(6),data2(7),data2(8),data2(9),
     >     data2(10),data2(11),data2(12),data2(13),data2(14),data2(15),
     >     data2(16)
      if (q2.gt.q2_r2 .or. q2.lt.q2_r1) then
         print *,'READPDF_SINGLE q2 interpolation doesnot match ',
     >        q2,q2_r2,q2_r1
         return
      endif
      do i=1,nq2bin-2
         READ (26,*,end=102) buf
      enddo
      READ(26,*,end=102) xb_r3,q2_r3,data3(1),data3(2),data3(3),
     >     data3(4),data3(5),data3(6),data3(7),data3(8),data3(9),
     >     data3(10),data3(11),data3(12),data3(13),data3(14),data3(15),
     >     data3(16)
      READ(26,*,end=102) xb_r4,q2_r4,data4(1),data4(2),data4(3),
     >     data4(4),data4(5),data4(6),data4(7),data4(8),data4(9),
     >     data4(10),data4(11),data4(12),data4(13),data4(14),data4(15),
     >     data4(16)

      if (q2.gt.q2_r4 .or. q2.lt.q2_r3) then
        print *,'READPDF_SINGLE q2 interpolation doesnot match ',
     >        q2,q2_r4,q2_r3
         return
      endif

      CLOSE(26)

      if (xb.gt.xb_r3 .or. xb.lt.xb_r1) then
        print *,'READPDF_SINGLE x interpolation doesnot match ',
     >        xb,xb_r3,xb_r1
         return
      endif

      do i=1,VDIM
         data12(i)=data2(i)*(q2-q2_r1)/(q2_r2-q2_r1)
     >        +data1(i)*(q2-q2_r2)/(q2_r1-q2_r2)
         data34(i)=data4(i)*(q2-q2_r3)/(q2_r4-q2_r3)
     >        +data3(i)*(q2-q2_r4)/(q2_r3-q2_r4)
         data0(i)=data12(i)*(xb-xb_r3)/(xb_r1-xb_r3)
     >        +data34(i)*(xb-xb_r1)/(xb_r3-xb_r1)
      enddo
      if (READDEBUG) print *,'READPDF_SINGLE: xb=',xb_r1,xb_r2,xb_r3,xb_r4
      if (READDEBUG) print *,'READPDF_SINGLE:',
     >     q2,q2_r1,q2_r2,data1(1),data2(1),data12(1)
      if (READDEBUG) print *,'READPDF_SINGLE:',
     >     q2,q2_r3,q2_r4,data3(1),data4(1),data34(1)
      if (READDEBUG) print *,'READPDF_SINGLE:',
     >     xb,xb_r1,xb_r3,data12(1),data34(1),data0(1)
      uv=data0(1)
      dv=data0(2)
      u=data0(3)
      d=data0(4)
      s=data0(5)
      c=data0(6)
      ubar=data0(7)
      dbar=data0(8)
      Rs=data0(9)
      Rv=data0(10)
      Rc=data0(11)
      Ad=data0(12)
      Ap=data0(13)
      a2=data0(14)
      b2=data0(15)
      c2=data0(16)

      if (READDEBUG) print *,'check READPDF',filename,a2,Ad,b2,a2*Ad+b2

      RETURN
 101  PRINT *,'cannot open ',filename,' READPDF_SINGLE failed'
      RETURN
 102  PRINT *,'error reading ',filename,' READPDF_SINGLE failed'
      RETURN
      END
