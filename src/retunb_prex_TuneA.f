C Forward transfer functions for right hrs with septum based on prex_retunb_dir.dat
c  setup for (x|theta)=(y|phi)=0 143 cm downstream of the 1st VDC
c HRS + PREX room temperature septum (right side)
c                     -JJL 3/29/2010

c This is tune A

c typical call: answer = function(x,5)
c INPUTS: x = 5 or more element array 
c              x(1)=x0  (meters)
c              x(2)=theta0 (really tan(theta0))
c              x(3)=y0   (meters)
c              x(4)=phi0 (really tan(phi0))
c              x(5)=delta (fractional value NOT percent)
c         M=5
c
c OUTPUT: units are the same as inputs
c 
c NOMENCLATURE: function name = prefix + _sp_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means path length
c     
c           suffixes:     fp means target to focus
c                         q1ex means target to Q1 exit
c                         q2ex means target to Q2 exit
c                         den  means target to Dipole entrance
c                         dex  means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3m  means target to Q3 middle
c                         q3ex means target to Q3 exit
c                         sen  means septum entrance
c                         sm   means halfway through the septum
c                         sex  means septum exit
c                         col  means target to Q1 collimator
c                         cq1x means collimator to q1 exit
c                         cden means collimator to dipole entrance
c                         cdex means collimator to dipole exit
c                         cq3e means collimator to q3 entrance
c                         cq3x means collimator to q3 exit
c                         cfp  means collimator to focus
c                               i.e. The plane perpendicular to the optic axis
c                                    that crosses the 1st VDC
c
c          _sp_ is for septum PREX
c
c APERTURES:
c    Coordinate systems are different in the spectrometer with septum model compared to the spectrometer
c     without septum. So some apertures even in the body of the spectrometer appear to be different. Numbers
c    given here supercede any numbers given in regard to transfer functions for the spectrometers alone.
c
c     sen: +0.088 < x < +0.382
c          -0.120 < y < 0.120
c      sm: +0.088 < x < +0.382
c          -0.120 < y < 0.120
c     sex: +0.088 < x < +0.382
c          -0.120 < y < 0.120
c     q1ex is a circle of radius 0.1492 m
c     q2ex is a circle of radius 0.30 m
c     den is a trapazoid:
c                                   -5.22008 < x < -4.98099
c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784  
c     dex is also a trapazoid: 
c                                   -0.46188 < x < 0.46188
c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
c     q3en is a circle of radius 0.3 m
c     q3m  is a circle of radius 0.3 m
c     q3ex is a circle of radius 0.3 m
c
 
      function x_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.8693400E-03/
      data xmin/
     1 -0.49945E-02,-0.51289E-01,-0.19995E-01,-0.28039E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.51728E-01, 0.19983E-01, 0.24966E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.70582121E-03, 0.12114862E+00, 0.46429415E-02,-0.46390048E-02,
     +  0.17536914E-02, 0.13469047E-02,-0.23920091E-02, 0.20589498E-02,
     +  0.32736890E-02,-0.80348284E-03, 0.24426731E-02,-0.42569268E-03,
     +  0.23728753E-03, 0.21928694E-03, 0.13193269E-02,-0.22622824E-02,
     + -0.21281645E-02, 0.67342375E-03,-0.62352454E-03,-0.31560531E-03,
     + -0.63657010E-03, 0.75510266E-03,-0.14704954E-02, 0.21381530E-02,
     +  0.20923070E-02, 0.33175373E-04,-0.16786551E-03, 0.26453927E-03,
     + -0.15351095E-03,-0.21072716E-03, 0.47645246E-03,-0.80300629E-03,
     +  0.62889059E-03,-0.11341627E-03,-0.38804760E-03,-0.35038491E-03,
     +  0.17100145E-04,-0.74109179E-04, 0.11945469E-03, 0.25634223E-03,
     +  0.16216822E-02, 0.15439111E-02,-0.30509630E-03,-0.28238414E-03,
     + -0.17734958E-03,-0.27791192E-03, 0.14806384E-04,-0.11377947E-03,
     + -0.32581336E-05, 0.14916761E-03, 0.45318110E-03, 0.89588248E-05,
     + -0.46678155E-03,-0.28656764E-03,-0.14966431E-04, 0.23982251E-04,
     +  0.29382605E-04,-0.21109119E-04, 0.64763204E-04,-0.40411625E-04,
     +  0.62747873E-04, 0.12576716E-04,-0.72294060E-04,-0.23513522E-03,
     +  0.32174931E-03, 0.26920531E-03, 0.56431742E-04,-0.27349273E-04,
     + -0.10696945E-02,-0.11544245E-03,-0.10531541E-02, 0.11700497E-03,
     +  0.59739233E-03, 0.41355105E-03, 0.77787023E-04,-0.24988090E-04,
     + -0.81546132E-04, 0.59645314E-04,-0.26735602E-03,-0.33874829E-04,
     + -0.31434619E-03,-0.62081293E-04,-0.16579621E-03, 0.13061732E-03,
     + -0.13902964E-03, 0.24936648E-03, 0.12220122E-03, 0.11066852E-03,
     +  0.29549975E-03,-0.74160249E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x23*x31        
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x23            
      x_sp_col    =x_sp_col    
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x31        
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x42    
      x_sp_col    =x_sp_col    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x24            
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21*x32        
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)                *x51
      x_sp_col    =x_sp_col    
     9  +coeff( 27)*x11    *x31        
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11*x22*x31        
     5  +coeff( 32)    *x24*x31        
     6  +coeff( 33)    *x23*x32        
     7  +coeff( 34)    *x23*x31    *x51
     8  +coeff( 35)        *x31*x41    
      x_sp_col    =x_sp_col    
     9  +coeff( 36)            *x42    
     1  +coeff( 37)            *x41*x51
     2  +coeff( 38)    *x22        *x51
     3  +coeff( 39)    *x21*x31    *x51
     4  +coeff( 40)*x11*x21    *x41    
     5  +coeff( 41)    *x22*x31*x41    
     6  +coeff( 42)    *x22    *x42    
     7  +coeff( 43)    *x23    *x41*x51
     8  +coeff( 44)*x11*x24            
      x_sp_col    =x_sp_col    
     9  +coeff( 45)*x11*x23*x31        
     1  +coeff( 46)    *x24*x32        
     2  +coeff( 47)        *x32*x42*x52
     3  +coeff( 48)        *x32        
     4  +coeff( 49)            *x41*x52
     5  +coeff( 50)*x11*x21*x31        
     6  +coeff( 51)    *x22*x32        
     7  +coeff( 52)*x12*x21            
     8  +coeff( 53)    *x21*x31*x42    
      x_sp_col    =x_sp_col    
     9  +coeff( 54)    *x21    *x43    
     1  +coeff( 55)    *x21*x31*x41*x51
     2  +coeff( 56)    *x21    *x42*x51
     3  +coeff( 57)    *x21        *x53
     4  +coeff( 58)*x11*x22        *x51
     5  +coeff( 59)    *x24        *x51
     6  +coeff( 60)            *x42*x52
     7  +coeff( 61)*x12*x22            
     8  +coeff( 62)    *x22*x31    *x52
      x_sp_col    =x_sp_col    
     9  +coeff( 63)    *x22    *x41*x52
     1  +coeff( 64)*x11*x23    *x41    
     2  +coeff( 65)    *x21*x32*x42    
     3  +coeff( 66)    *x21*x31*x43    
     4  +coeff( 67)    *x21*x33    *x51
     5  +coeff( 68)*x11*x22*x32        
     6  +coeff( 69)    *x24*x31*x41    
     7  +coeff( 70)*x11*x22    *x42    
     8  +coeff( 71)    *x24    *x42    
      x_sp_col    =x_sp_col    
     9  +coeff( 72)    *x23*x32*x41    
     1  +coeff( 73)    *x23*x31*x42    
     2  +coeff( 74)    *x23    *x43    
     3  +coeff( 75)    *x23*x32    *x51
     4  +coeff( 76)    *x23*x31    *x52
     5  +coeff( 77)*x11*x21    *x41*x52
     6  +coeff( 78)    *x23    *x41*x52
     7  +coeff( 79)*x11*x24*x31        
     8  +coeff( 80)*x12*x22    *x41    
      x_sp_col    =x_sp_col    
     9  +coeff( 81)*x11*x24    *x41    
     1  +coeff( 82)*x12*x21*x31*x41    
     2  +coeff( 83)    *x21*x32*x43    
     3  +coeff( 84)    *x21*x33*x41*x51
     4  +coeff( 85)*x12*x24            
     5  +coeff( 86)    *x23*x33*x41    
     6  +coeff( 87)    *x23    *x41*x53
     7  +coeff( 88)*x11*x24*x31*x41    
     8  +coeff( 89)*x11*x24    *x42    
      x_sp_col    =x_sp_col    
     9  +coeff( 90)*x11*x24    *x41*x51
c
      return
      end
      function t_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.6932339E-03/
      data xmin/
     1 -0.49945E-02,-0.51289E-01,-0.19995E-01,-0.28039E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.51728E-01, 0.19983E-01, 0.24966E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.83620696E-04, 0.44951227E-01,-0.56202146E-02,-0.28576415E-02,
     +  0.19809101E-02,-0.40304285E-03, 0.86165749E-03, 0.38275805E-02,
     + -0.71496045E-03, 0.30154779E-03, 0.23360619E-04, 0.12845623E-02,
     +  0.42706993E-03, 0.19738008E-02,-0.54499556E-04,-0.38641863E-03,
     +  0.16357985E-03, 0.71001245E-03,-0.19747554E-02,-0.18264749E-02,
     +  0.46463785E-03, 0.30408088E-04,-0.31700774E-03,-0.48892276E-03,
     +  0.17445222E-03, 0.19094801E-04, 0.56998723E-03, 0.19131680E-02,
     +  0.17269318E-02, 0.13455478E-03,-0.86988075E-04,-0.24374411E-03,
     + -0.21933242E-03,-0.16406791E-03, 0.11735272E-03, 0.57516736E-03,
     +  0.53985004E-03,-0.99901437E-04, 0.31331589E-03, 0.69703630E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x23            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23    *x41    
      t_sp_col    =t_sp_col    
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x22*x33        
     7  +coeff( 16)        *x31        
     8  +coeff( 17)*x11*x21            
      t_sp_col    =t_sp_col    
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x23*x32        
     4  +coeff( 22)                *x51
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x21*x32        
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)*x11    *x33        
      t_sp_col    =t_sp_col    
     9  +coeff( 27)*x11*x22    *x41    
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)    *x21*x33    *x51
     4  +coeff( 31)*x11*x22*x33        
     5  +coeff( 32)        *x31*x41    
     6  +coeff( 33)            *x42    
     7  +coeff( 34)*x11    *x31        
     8  +coeff( 35)*x11*x21    *x41    
      t_sp_col    =t_sp_col    
     9  +coeff( 36)    *x22*x31*x41    
     1  +coeff( 37)    *x22    *x42    
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)*x11*x22*x31        
     4  +coeff( 40)*x11*x23*x31        
c
      return
      end
      function y_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.9400545E-02/
      data xmin/
     1 -0.49945E-02,-0.51289E-01,-0.19995E-01,-0.28039E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.51728E-01, 0.19983E-01, 0.24966E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56388811E-02, 0.20872107E-01, 0.66277854E-01, 0.46263677E-02,
     + -0.13928600E-02,-0.40211123E-02,-0.26263238E-02,-0.41340516E-03,
     + -0.14659781E-02,-0.23572592E-03,-0.52651780E-03,-0.42151217E-03,
     +  0.92696626E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)*x11*x21            
      y_sp_col    =y_sp_col    
     9  +coeff(  9)    *x22*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21*x33        
     4  +coeff( 13)    *x23            
c
      return
      end
      function p_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4096258E-02/
      data xmin/
     1 -0.49945E-02,-0.51289E-01,-0.19995E-01,-0.28039E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.51728E-01, 0.19983E-01, 0.24966E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31141571E-02,-0.17120824E-02, 0.12601212E-02, 0.29153345E-01,
     +  0.66118971E-02,-0.45902161E-02,-0.68394491E-03, 0.13487670E-02,
     + -0.28711697E-02,-0.63863998E-04,-0.11437887E-02,-0.35106958E-03,
     + -0.15614325E-02, 0.51002018E-03,-0.61430060E-03, 0.12888956E-02,
     +  0.12120340E-02, 0.26729456E-03,-0.55863743E-03, 0.10180962E-02,
     + -0.25376137E-02,-0.24333524E-02, 0.46913254E-04,-0.98617544E-04,
     +  0.27577573E-03, 0.23445646E-03,-0.29499296E-03, 0.92352981E-04,
     +  0.22172107E-03,-0.98105666E-04,-0.15158854E-03, 0.29800841E-03,
     +  0.20716446E-03, 0.27940964E-03,-0.80057391E-03, 0.16865933E-03,
     + -0.13526005E-03,-0.28347599E-03,-0.33452365E-03,-0.11706002E-04,
     +  0.74946831E-04,-0.68629583E-04,-0.90253860E-04,-0.26793597E-03,
     +  0.63042593E-04,-0.66949047E-04,-0.24106317E-03, 0.17311072E-03,
     + -0.98844110E-04,-0.12853389E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x23            
      p_sp_col    =p_sp_col    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x22*x33        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)            *x42    
      p_sp_col    =p_sp_col    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x22*x32        
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)        *x32    *x52
     6  +coeff( 24)*x11                
     7  +coeff( 25)        *x32        
     8  +coeff( 26)*x11*x22            
      p_sp_col    =p_sp_col    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x12*x21            
     2  +coeff( 29)    *x21*x31*x42    
     3  +coeff( 30)        *x31    *x53
     4  +coeff( 31)            *x41*x53
     5  +coeff( 32)*x11*x23            
     6  +coeff( 33)*x11*x23*x31        
     7  +coeff( 34)*x12*x21    *x43    
     8  +coeff( 35)    *x23*x31*x43*x52
      p_sp_col    =p_sp_col    
     9  +coeff( 36)    *x21        *x51
     1  +coeff( 37)    *x21*x32        
     2  +coeff( 38)    *x21*x31*x41    
     3  +coeff( 39)    *x21    *x42    
     4  +coeff( 40)        *x31*x41*x51
     5  +coeff( 41)            *x42*x51
     6  +coeff( 42)            *x41*x52
     7  +coeff( 43)                *x53
     8  +coeff( 44)*x11*x21*x31        
      p_sp_col    =p_sp_col    
     9  +coeff( 45)*x11*x21        *x51
     1  +coeff( 46)*x12    *x31        
     2  +coeff( 47)    *x23        *x51
     3  +coeff( 48)    *x22    *x41*x51
     4  +coeff( 49)    *x22        *x52
     5  +coeff( 50)    *x21*x31    *x52
c
      return
      end
      function l_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3719136E-03/
      data xmin/
     1 -0.49945E-02,-0.51289E-01,-0.19995E-01,-0.28039E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.51728E-01, 0.19983E-01, 0.24966E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73199393E-03,-0.26609884E-02,-0.60257334E-02,-0.13234274E-03,
     + -0.29994326E-02,-0.69621849E-04,-0.93470572E-03,-0.12899714E-03,
     +  0.34094718E-03, 0.13093649E-03, 0.10985035E-03, 0.40343948E-04,
     +  0.26697844E-04,-0.10956508E-03, 0.15521559E-03, 0.16174931E-03,
     +  0.43854423E-04,-0.67162197E-04, 0.12729209E-04, 0.42255026E-04,
     + -0.45822024E-04,-0.13926605E-04, 0.43055866E-05,-0.14414318E-05,
     +  0.35196048E-04,-0.11544387E-04, 0.45109218E-05,-0.60053317E-05,
     +  0.17775828E-05, 0.22600252E-04,-0.11862204E-04,-0.11119253E-04,
     + -0.19146062E-04, 0.45894444E-05,-0.89871610E-05,-0.15069944E-04,
     +  0.11465439E-04,-0.62675040E-05,-0.39040355E-04, 0.16200491E-04,
     +  0.11103752E-04, 0.19782448E-04, 0.91448437E-05,-0.27022932E-04,
     +  0.17121938E-04,-0.73766209E-05, 0.89125060E-05, 0.18446306E-04,
     + -0.21434069E-04, 0.27596718E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x41*x51
      l_sp_col    =l_sp_col    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x22*x31*x41    
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)    *x21            
      l_sp_col    =l_sp_col    
     9  +coeff( 18)    *x23            
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)*x11*x23*x32        
     5  +coeff( 23)    *x21        *x51
     6  +coeff( 24)                *x52
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)            *x42*x51
      l_sp_col    =l_sp_col    
     9  +coeff( 27)                *x53
     1  +coeff( 28)*x11*x22            
     2  +coeff( 29)*x11*x21*x31        
     3  +coeff( 30)    *x22*x32        
     4  +coeff( 31)    *x21*x31*x42    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)    *x21*x31    *x52
     8  +coeff( 35)        *x32    *x52
      l_sp_col    =l_sp_col    
     9  +coeff( 36)*x11*x23            
     1  +coeff( 37)*x11*x21    *x42    
     2  +coeff( 38)*x12    *x32        
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)*x12        *x41*x51
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)    *x21*x31*x41*x52
     7  +coeff( 43)    *x22        *x53
     8  +coeff( 44)*x11*x23    *x41    
      l_sp_col    =l_sp_col    
     9  +coeff( 45)*x11*x21*x31    *x52
     1  +coeff( 46)*x12*x23            
     2  +coeff( 47)*x12    *x31*x42    
     3  +coeff( 48)*x12        *x42*x51
     4  +coeff( 49)    *x23    *x41*x52
     5  +coeff( 50)    *x22    *x42*x52
c
      return
      end
      function x_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1663723E-03/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22172378E-02, 0.14830737E+00, 0.25612949E-02,-0.11674454E-01,
     +  0.27715706E-02, 0.42169252E-02, 0.26199063E-02,-0.60055922E-02,
     +  0.38189800E-02,-0.17085378E-02, 0.52677542E-02,-0.71749707E-06,
     +  0.78984685E-02,-0.89876784E-03, 0.44159850E-03, 0.27356660E-02,
     +  0.13449113E-02,-0.12908381E-02,-0.74454362E-03,-0.54028328E-02,
     + -0.52231513E-02, 0.13569556E-02, 0.12861788E-04,-0.73950214E-03,
     + -0.70699828E-03,-0.39421167E-03, 0.10364457E-03,-0.13917184E-02,
     + -0.11029733E-03,-0.58784246E-06, 0.20816911E-02, 0.23254594E-02,
     + -0.43034132E-03, 0.76518371E-03, 0.14699483E-02,-0.35680798E-02,
     +  0.57095345E-02, 0.57090800E-02, 0.49976876E-04,-0.24742362E-03,
     +  0.17887302E-03, 0.31954175E-03, 0.61806000E-03, 0.67972700E-03,
     + -0.18316898E-02,-0.56418363E-03,-0.34213689E-03, 0.12127443E-03,
     + -0.43609631E-03,-0.15445771E-02,-0.33717003E-03, 0.62135652E-04,
     + -0.57362460E-04, 0.35148216E-03,-0.23808090E-03, 0.72020113E-04,
     +  0.98793680E-04,-0.22359402E-03,-0.74599090E-03, 0.22707946E-03,
     + -0.20231959E-04, 0.98499666E-04, 0.12057597E-02,-0.36606024E-03,
     + -0.10968489E-03,-0.53137116E-03,-0.26446485E-03, 0.29378227E-03,
     + -0.90293499E-04,-0.12656212E-03, 0.19859835E-02, 0.20617119E-02,
     +  0.21525900E-03, 0.19756092E-02, 0.98769460E-03,-0.62975887E-03,
     +  0.71084680E-03, 0.14046850E-02, 0.15198187E-02, 0.45728768E-03,
     +  0.63132978E-03,-0.13080082E-02,-0.70923980E-03,-0.25919991E-03,
     + -0.32882910E-03, 0.98466955E-03, 0.16495107E-04,-0.93253017E-04,
     + -0.70251939E-04, 0.18860113E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x33        
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)        *x31        
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)*x11*x22            
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 18)    *x24            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)    *x23*x32        
     5  +coeff( 23)*x11    *x33        
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)            *x42    
     8  +coeff( 26)*x11    *x31        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)    *x21*x32        
     2  +coeff( 29)    *x21        *x52
     3  +coeff( 30)        *x32    *x51
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)*x11*x22    *x41    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 36)    *x24    *x41    
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)                *x51
     4  +coeff( 40)        *x32        
     5  +coeff( 41)    *x21*x31    *x51
     6  +coeff( 42)    *x21    *x41*x51
     7  +coeff( 43)*x11*x21    *x41    
     8  +coeff( 44)    *x22*x32        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 45)    *x24*x31        
     1  +coeff( 46)*x11*x24            
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x11*x21*x33        
     4  +coeff( 49)*x11*x23*x33        
     5  +coeff( 50)    *x24    *x42*x52
     6  +coeff( 51)*x12*x24*x32        
     7  +coeff( 52)            *x41*x51
     8  +coeff( 53)    *x22        *x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 54)*x11*x21*x31        
     1  +coeff( 55)*x11        *x42    
     2  +coeff( 56)*x12*x21            
     3  +coeff( 57)    *x21*x31    *x52
     4  +coeff( 58)    *x22    *x43    
     5  +coeff( 59)*x11*x23    *x41    
     6  +coeff( 60)    *x21*x31*x43    
     7  +coeff( 61)*x11*x23        *x51
     8  +coeff( 62)*x12    *x32        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 63)    *x23    *x43    
     1  +coeff( 64)    *x23*x31    *x52
     2  +coeff( 65)*x11    *x33*x41    
     3  +coeff( 66)*x11    *x31*x43    
     4  +coeff( 67)    *x22    *x43*x51
     5  +coeff( 68)    *x22*x32    *x52
     6  +coeff( 69)*x11    *x31*x41*x52
     7  +coeff( 70)*x11*x22        *x53
     8  +coeff( 71)    *x23*x32*x42    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 72)    *x23*x31*x43    
     1  +coeff( 73)*x11*x21*x31*x41*x52
     2  +coeff( 74)*x11*x24*x31*x41    
     3  +coeff( 75)*x11*x24    *x42    
     4  +coeff( 76)    *x24*x32*x42    
     5  +coeff( 77)    *x23*x33*x42    
     6  +coeff( 78)    *x23*x32*x41*x52
     7  +coeff( 79)    *x23*x31*x42*x52
     8  +coeff( 80)*x11    *x33*x43    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 81)*x12*x23    *x43    
     1  +coeff( 82)*x11*x24*x33*x41    
     2  +coeff( 83)*x12*x22*x32    *x52
     3  +coeff( 84)*x12*x22    *x41*x53
     4  +coeff( 85)*x11*x23*x32    *x53
     5  +coeff( 86)*x12*x22*x31*x42*x52
     6  +coeff( 87)        *x31    *x51
     7  +coeff( 88)        *x31*x42    
     8  +coeff( 89)            *x43    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 90)        *x31    *x52
c
      return
      end
      function t_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3012097E-03/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15726645E-03, 0.16899521E-02,-0.63567513E-05,-0.16546772E-02,
     +  0.37814278E-03,-0.15443350E-02,-0.30544649E-02, 0.21753672E-02,
     +  0.55026889E-04, 0.73498004E-03, 0.34057823E-04,-0.45092405E-04,
     +  0.44102212E-05,-0.36571240E-06, 0.21766247E-02,-0.11518413E-04,
     +  0.51488023E-03,-0.44457276E-04,-0.18780041E-03,-0.36140264E-03,
     + -0.16643292E-03, 0.68438858E-04, 0.31634147E-03, 0.60220633E-03,
     + -0.13050295E-02,-0.11166438E-03, 0.18906929E-03, 0.10566910E-02,
     +  0.27005688E-04, 0.10924365E-02, 0.19362044E-04, 0.28535658E-04,
     + -0.29080035E-03,-0.12218689E-02,-0.79022800E-04, 0.31730274E-03,
     +  0.86967571E-03, 0.40543196E-03,-0.13592499E-03,-0.12720408E-03,
     + -0.13992978E-03,-0.99414472E-04, 0.29397859E-04, 0.42890089E-04,
     +  0.37621593E-03, 0.22221951E-04, 0.48741011E-03, 0.74530493E-04,
     +  0.24206722E-03, 0.29089305E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)        *x32*x41    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)        *x33*x42    
     8  +coeff( 17)    *x22    *x43    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 18)    *x23*x33        
     1  +coeff( 19)        *x31        
     2  +coeff( 20)            *x41    
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x22*x31        
     6  +coeff( 24)    *x22    *x41    
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)    *x21        *x52
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)    *x23*x31        
     2  +coeff( 29)*x11    *x33        
     3  +coeff( 30)    *x23*x31*x41    
     4  +coeff( 31)    *x21*x33*x41    
     5  +coeff( 32)    *x23*x33*x41    
     6  +coeff( 33)    *x21*x32        
     7  +coeff( 34)    *x21*x31*x41    
     8  +coeff( 35)    *x21    *x41*x51
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)*x11*x22*x33        
     4  +coeff( 40)        *x31*x41    
     5  +coeff( 41)            *x42    
     6  +coeff( 42)*x11    *x31        
     7  +coeff( 43)    *x22        *x51
     8  +coeff( 44)*x11*x21    *x41    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 45)    *x22*x31*x41    
     1  +coeff( 46)    *x21*x32*x41    
     2  +coeff( 47)    *x22    *x42    
     3  +coeff( 48)    *x23        *x51
     4  +coeff( 49)*x11*x22*x31        
     5  +coeff( 50)    *x23*x32        
c
      return
      end
      function y_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.1062998E-01/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24409750E-02, 0.30434271E-01, 0.14350954E+00, 0.19196080E-01,
     + -0.49238042E-02, 0.46034413E-02,-0.13327494E-01,-0.16491811E-02,
     + -0.11495489E-02,-0.17656893E-02,-0.47069797E-02,-0.83389888E-02,
     +  0.40177163E-02,-0.15085756E-01, 0.44188122E-02,-0.52653038E-03,
     + -0.13217647E-02,-0.24320013E-02, 0.10935171E-02,-0.12648975E-01,
     + -0.30318890E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)            *x41*x51
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)    *x21*x31        
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)        *x34        
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22*x34        
c
      return
      end
      function p_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4036574E-02/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12018534E-02, 0.80871787E-02, 0.55243935E-01, 0.94823958E-02,
     + -0.63075302E-02, 0.14855913E-02,-0.21497600E-02,-0.14824192E-02,
     + -0.90969872E-03,-0.18648012E-02,-0.36160808E-03,-0.64055744E-03,
     + -0.23411585E-02,-0.37523888E-02, 0.10257369E-02,-0.10342359E-02,
     +  0.22513731E-02, 0.21179600E-02, 0.32870620E-03,-0.64936499E-02,
     + -0.75523946E-02,-0.16339276E-03, 0.47693923E-03, 0.11947934E-03,
     +  0.41307937E-03,-0.94765739E-03,-0.13136382E-02, 0.16182216E-02,
     +  0.29277013E-03,-0.94644830E-03,-0.10745311E-02, 0.20627916E-03,
     +  0.60178171E-03,-0.69089164E-03,-0.12339437E-02,-0.44542740E-03,
     + -0.78887657E-04,-0.17210640E-03,-0.36259592E-03,-0.14364174E-03,
     + -0.18964025E-03,-0.34256352E-03, 0.15593519E-02,-0.96443499E-03,
     +  0.16277068E-03, 0.77906158E-03, 0.48949645E-03,-0.61531493E-03,
     +  0.62883511E-03, 0.73015742E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x21            
     8  +coeff(  8)            *x41*x51
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)        *x31*x41    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 18)            *x42    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11                
     5  +coeff( 23)        *x32        
     6  +coeff( 24)    *x21        *x51
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)*x11*x21    *x41    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 27)    *x22*x32        
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x21*x31*x41    
     4  +coeff( 31)*x11*x21    *x42    
     5  +coeff( 32)*x11*x23*x31        
     6  +coeff( 33)    *x22*x33    *x52
     7  +coeff( 34)    *x21*x31*x41    
     8  +coeff( 35)    *x21    *x42    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 36)*x11*x21*x31        
     1  +coeff( 37)*x11    *x31    *x51
     2  +coeff( 38)*x12*x21            
     3  +coeff( 39)        *x31*x43    
     4  +coeff( 40)        *x31    *x53
     5  +coeff( 41)*x11*x21    *x41*x51
     6  +coeff( 42)    *x23*x32        
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)    *x22    *x43    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 45)        *x32*x43    
     1  +coeff( 46)*x11*x23    *x41    
     2  +coeff( 47)*x12*x23            
     3  +coeff( 48)    *x22*x32*x41*x51
     4  +coeff( 49)    *x23    *x41*x52
     5  +coeff( 50)    *x22*x31*x41*x52
c
      return
      end
      function l_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.9150847E-03/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13235277E-02, 0.11099711E-03,-0.26360319E-02,-0.57485751E-02,
     + -0.33821820E-02,-0.33135733E-04,-0.55089651E-03,-0.27562783E-02,
     + -0.81661332E-04,-0.78567344E-03, 0.81000169E-03,-0.10757502E-03,
     +  0.27937276E-03, 0.20933240E-03,-0.68893467E-04,-0.93366420E-04,
     +  0.67476853E-04, 0.29973475E-04, 0.38910941E-04, 0.89828180E-04,
     +  0.69958580E-04, 0.95041185E-04,-0.31751103E-03, 0.66728576E-03,
     + -0.84290014E-05,-0.18125898E-03, 0.28277543E-05, 0.97437187E-04,
     + -0.18697758E-04, 0.22961533E-04, 0.17586484E-04,-0.85352789E-04,
     +  0.61027787E-03, 0.10121156E-04, 0.50870562E-03,-0.31437554E-04,
     +  0.19589565E-04, 0.39709568E-04, 0.38571429E-03, 0.52937961E-04,
     + -0.13031682E-03,-0.14403355E-03, 0.66150583E-05, 0.72493291E-04,
     +  0.12585365E-04, 0.18928486E-04, 0.51012277E-04, 0.26222289E-04,
     + -0.23124805E-04, 0.23994848E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)                *x51
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x21*x31        
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 18)    *x21        *x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)            *x42*x53
     8  +coeff( 26)    *x22*x33*x41    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 27)        *x32*x41    
     1  +coeff( 28)    *x21    *x42    
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)    *x23*x31        
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x22    *x43    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 36)        *x32*x43    
     1  +coeff( 37)        *x31*x41*x53
     2  +coeff( 38)        *x33*x43    
     3  +coeff( 39)    *x22*x33*x42    
     4  +coeff( 40)    *x21*x31*x41    
     5  +coeff( 41)        *x31*x42    
     6  +coeff( 42)            *x43    
     7  +coeff( 43)        *x31*x41*x51
     8  +coeff( 44)            *x42*x51
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 45)        *x31    *x52
     1  +coeff( 46)                *x53
     2  +coeff( 47)    *x22*x32        
     3  +coeff( 48)    *x21*x32*x41    
     4  +coeff( 49)    *x23        *x51
     5  +coeff( 50)    *x22*x31    *x51
c
      return
      end
      function x_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1451013E-02/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43754159E-02, 0.25693557E+00,-0.30711483E-01, 0.16788099E-03,
     +  0.66469372E-02,-0.15859794E-01, 0.72349948E-02, 0.92936000E-02,
     + -0.44609648E-02,-0.19673225E-02,-0.23213080E-02, 0.13496816E-01,
     +  0.21280661E-01,-0.44047632E-02,-0.52143697E-03, 0.10981939E-02,
     +  0.68787797E-02,-0.20035673E-02,-0.15153597E-01, 0.34610508E-02,
     + -0.31181816E-02, 0.11170815E-01, 0.16423933E-01,-0.10239078E-04,
     +  0.30575630E-02,-0.61479851E-03, 0.13278943E-03,-0.34554747E-02,
     +  0.17832842E-03, 0.40544014E-03,-0.37229713E-02,-0.15827598E-01,
     +  0.32583126E-03,-0.10904504E-02, 0.41312077E-02,-0.84018856E-02,
     +  0.16064551E-01,-0.15802348E-05,-0.18202494E-02,-0.16604566E-02,
     + -0.99092722E-03, 0.49214071E-03,-0.31089829E-03, 0.15568759E-02,
     +  0.53473460E-02, 0.50636935E-02, 0.20642544E-02, 0.40481552E-02,
     + -0.14683430E-02,-0.74854982E-03,-0.96364594E-04,-0.53188769E-03,
     +  0.11250368E-03,-0.13573981E-03, 0.10799296E-02, 0.19563693E-02,
     +  0.19611372E-03, 0.11077136E-02,-0.14163664E-02,-0.18834122E-02,
     +  0.31739338E-02,-0.11080729E-02, 0.10239020E-03,-0.54723694E-03,
     + -0.34526743E-02,-0.29018873E-03,-0.54929620E-02,-0.30185208E-02,
     +  0.17237671E-02, 0.46643587E-02, 0.37858917E-02, 0.30068760E-02,
     + -0.42570333E-03, 0.45554654E-04, 0.68353271E-04,-0.39700459E-03,
     + -0.38568734E-03,-0.12475663E-02,-0.11673815E-02, 0.27937102E-03,
     +  0.54514740E-03, 0.76892553E-03, 0.35571420E-03,-0.24970801E-03,
     +  0.15790065E-03,-0.59641211E-03, 0.50212466E-02,-0.82638313E-03,
     +  0.44822553E-03, 0.42253110E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x33        
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff(  9)            *x41    
     1  +coeff( 10)*x11                
     2  +coeff( 11)        *x31        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)    *x24*x31        
     6  +coeff( 15)    *x23*x33        
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x22*x31        
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)*x11    *x33        
     7  +coeff( 25)    *x21*x33*x41    
     8  +coeff( 26)    *x21*x32    *x52
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 27)    *x21*x31*x41*x52
     1  +coeff( 28)    *x23*x33*x41    
     2  +coeff( 29)                *x51
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)    *x21*x32        
     5  +coeff( 32)    *x21*x31*x41    
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)*x11*x23            
     8  +coeff( 35)*x11*x22    *x41    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 36)    *x24    *x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)*x11*x24*x31        
     3  +coeff( 39)        *x31*x41    
     4  +coeff( 40)            *x42    
     5  +coeff( 41)*x11    *x31        
     6  +coeff( 42)    *x21    *x41*x51
     7  +coeff( 43)    *x21        *x52
     8  +coeff( 44)*x11*x21    *x41    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 45)    *x22*x31*x41    
     1  +coeff( 46)    *x22    *x42    
     2  +coeff( 47)*x11*x22*x31        
     3  +coeff( 48)    *x23*x32        
     4  +coeff( 49)*x11*x24            
     5  +coeff( 50)    *x24*x32        
     6  +coeff( 51)*x11*x21*x33        
     7  +coeff( 52)        *x32        
     8  +coeff( 53)            *x41*x51
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 54)    *x22        *x51
     1  +coeff( 55)*x11*x21*x31        
     2  +coeff( 56)    *x22*x32        
     3  +coeff( 57)*x12*x21            
     4  +coeff( 58)    *x22    *x41*x52
     5  +coeff( 59)*x11*x23*x31        
     6  +coeff( 60)*x11*x23    *x41    
     7  +coeff( 61)    *x21*x31*x43    
     8  +coeff( 62)    *x23*x31    *x52
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 63)    *x21    *x42*x53
     1  +coeff( 64)*x12*x24            
     2  +coeff( 65)    *x24    *x41*x52
     3  +coeff( 66)*x11*x22        *x53
     4  +coeff( 67)    *x23*x32*x42    
     5  +coeff( 68)    *x23    *x43*x51
     6  +coeff( 69)*x11*x24*x31*x41    
     7  +coeff( 70)    *x21*x32*x42*x52
     8  +coeff( 71)    *x21*x31*x43*x52
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 72)    *x23    *x43*x52
     1  +coeff( 73)*x11    *x31*x43*x53
     2  +coeff( 74)        *x31    *x51
     3  +coeff( 75)*x12                
     4  +coeff( 76)*x11    *x31*x41    
     5  +coeff( 77)*x11        *x42    
     6  +coeff( 78)    *x21*x31*x42    
     7  +coeff( 79)    *x21    *x43    
     8  +coeff( 80)    *x21*x32    *x51
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 81)    *x21*x31*x41*x51
     1  +coeff( 82)    *x21    *x42*x51
     2  +coeff( 83)    *x21*x31    *x52
     3  +coeff( 84)    *x23*x31    *x51
     4  +coeff( 85)*x11    *x32*x41    
     5  +coeff( 86)    *x22    *x43    
     6  +coeff( 87)    *x21*x32*x42    
     7  +coeff( 88)    *x21    *x42*x52
     8  +coeff( 89)    *x21*x31    *x53
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 90)    *x21    *x41*x53
c
      return
      end
      function t_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.6105365E-03/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13028642E-02, 0.75126745E-01,-0.15271656E-02,-0.10655079E-01,
     +  0.37749135E-02, 0.13818472E-02,-0.55080005E-02, 0.29447670E-02,
     + -0.13571661E-02, 0.88599866E-03, 0.64357278E-04, 0.76812478E-02,
     + -0.74798759E-03, 0.24432370E-02,-0.53605549E-02, 0.75045973E-03,
     + -0.25791320E-03, 0.44579264E-02, 0.92507405E-04, 0.31455411E-03,
     + -0.61401527E-03, 0.14230994E-02,-0.48181945E-02,-0.64065375E-05,
     +  0.10771607E-02,-0.59769087E-03,-0.49206754E-03,-0.30200661E-03,
     +  0.14198141E-03,-0.11441448E-02, 0.21112573E-03,-0.19194762E-03,
     +  0.23862277E-02, 0.57788927E-03, 0.11784706E-02, 0.39734528E-02,
     + -0.22379090E-05, 0.79383688E-04, 0.26111040E-03, 0.27548592E-03,
     +  0.19845350E-02,-0.25924919E-04,-0.26103435E-03, 0.63247877E-04,
     +  0.21399428E-03,-0.85236774E-04,-0.20279257E-03, 0.16339347E-03,
     +  0.49148826E-03,-0.49532618E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x23*x31        
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x23            
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)        *x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x22*x33        
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)                *x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)    *x22*x31        
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)*x11    *x33        
     7  +coeff( 25)    *x23*x32        
     8  +coeff( 26)        *x31*x41    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 27)            *x42    
     1  +coeff( 28)*x11    *x31        
     2  +coeff( 29)*x11            *x51
     3  +coeff( 30)    *x21*x32        
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)    *x22    *x42    
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)*x11*x22    *x41    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x22*x33*x41    
     2  +coeff( 38)            *x41*x51
     3  +coeff( 39)*x11*x21    *x41    
     4  +coeff( 40)    *x22*x32        
     5  +coeff( 41)    *x22*x31*x41    
     6  +coeff( 42)        *x33*x41    
     7  +coeff( 43)        *x32*x42    
     8  +coeff( 44)        *x33    *x51
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 45)    *x21    *x42*x51
     1  +coeff( 46)        *x32    *x52
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)    *x21    *x41*x53
     4  +coeff( 49)*x11*x23*x31        
     5  +coeff( 50)*x11*x23*x33        
c
      return
      end
      function y_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.1315558E-01/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35219593E-02, 0.30649807E-01, 0.18212245E+00, 0.28592531E-01,
     + -0.72472468E-02, 0.69052652E-02,-0.19639716E-01,-0.25764548E-02,
     +  0.16397507E-02,-0.10554438E-02,-0.68422505E-02,-0.12850484E-01,
     +  0.59131617E-02,-0.22482000E-01, 0.64920210E-02,-0.18880155E-02,
     + -0.34506551E-02, 0.15852231E-02,-0.18534021E-01,-0.43636509E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)*x11*x21            
      y_sp_q2ex   =y_sp_q2ex   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)    *x21    *x41    
      y_sp_q2ex   =y_sp_q2ex   
     9  +coeff( 18)        *x34        
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)    *x22*x34        
c
      return
      end
      function p_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2362812E-02/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.51079778E-03,-0.78085400E-02,-0.29419733E-01,-0.32013967E-02,
     +  0.22524025E-02, 0.36686361E-02, 0.72468811E-03,-0.60741976E-03,
     +  0.81619649E-03, 0.67995535E-03,-0.67103759E-03, 0.34093604E-03,
     +  0.62639487E-03, 0.69875590E-03,-0.78121276E-03,-0.18143308E-03,
     +  0.29238011E-03,-0.31824192E-03, 0.27015104E-03, 0.37479709E-03,
     + -0.80776581E-03,-0.20576823E-03, 0.26252074E-03,-0.62746223E-03,
     +  0.22490730E-02, 0.26507473E-02,-0.25988498E-04, 0.54271641E-04,
     + -0.20438139E-03,-0.36373246E-03,-0.51398051E-03, 0.12414175E-04,
     + -0.10551811E-03,-0.19336680E-04,-0.13128293E-03, 0.17210346E-03,
     +  0.55471278E-03, 0.14485240E-03,-0.15171607E-03, 0.33017789E-03,
     +  0.34996285E-03, 0.17789718E-02, 0.19835012E-02,-0.21447323E-03,
     + -0.18104233E-03, 0.31477152E-03, 0.27460908E-03,-0.87448920E-04,
     + -0.65239612E-04, 0.27815092E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x23            
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff(  9)    *x21            
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)    *x22        *x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x21        *x51
     8  +coeff( 17)    *x22    *x41    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)        *x32*x42    
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)            *x42    
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)    *x22    *x42    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 27)*x11*x23*x31        
     1  +coeff( 28)*x11                
     2  +coeff( 29)        *x32        
     3  +coeff( 30)        *x31*x42    
     4  +coeff( 31)            *x43    
     5  +coeff( 32)        *x31*x41*x51
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)        *x31    *x52
     8  +coeff( 35)*x11*x22            
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 36)*x11*x21*x31        
     1  +coeff( 37)    *x22*x32        
     2  +coeff( 38)    *x23        *x51
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)*x11*x21*x31*x41    
     5  +coeff( 41)*x11*x21    *x42    
     6  +coeff( 42)    *x22*x31*x42    
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)    *x21    *x43*x51
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 45)    *x22*x31    *x52
     1  +coeff( 46)    *x21*x31*x41*x52
     2  +coeff( 47)    *x21    *x42*x52
     3  +coeff( 48)    *x21*x31    *x53
     4  +coeff( 49)*x11*x21        *x53
     5  +coeff( 50)        *x33*x43    
c
      return
      end
      function l_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2020909E-02/
      data xmin/
     1 -0.49945E-02,-0.49894E-01,-0.19959E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.49454E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24425848E-02, 0.25053124E-03,-0.25989290E-02,-0.55761812E-02,
     + -0.37166931E-05,-0.67816088E-02,-0.72154900E-04,-0.11080862E-02,
     + -0.45260210E-02,-0.15694587E-03,-0.13704542E-02, 0.16243160E-02,
     + -0.24280263E-04, 0.50252600E-03,-0.12056854E-03, 0.19737362E-03,
     + -0.17779680E-03, 0.50101947E-03, 0.13814283E-03,-0.57812216E-03,
     +  0.11431341E-02,-0.17925815E-03,-0.68901289E-04, 0.12101118E-03,
     +  0.27545055E-04, 0.22637309E-03, 0.16446442E-03,-0.18451874E-03,
     +  0.87461894E-03, 0.15005658E-02, 0.10596973E-04, 0.78791374E-04,
     +  0.43675356E-03, 0.82338825E-04, 0.19128619E-03,-0.13252220E-03,
     +  0.15923277E-04, 0.98285534E-04, 0.23146673E-04,-0.38021288E-04,
     +  0.42084288E-04, 0.91938055E-04, 0.45925575E-04,-0.51231229E-04,
     +  0.44580425E-04, 0.41826428E-04,-0.29339199E-03, 0.74844004E-03,
     + -0.27995170E-03,-0.24611401E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22*x33        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x22*x33*x41    
     5  +coeff( 23)                *x51
     6  +coeff( 24)    *x21*x31        
     7  +coeff( 25)    *x21        *x51
     8  +coeff( 26)            *x42*x51
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x23*x31        
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)    *x22    *x43    
     4  +coeff( 31)        *x33*x41*x51
     5  +coeff( 32)    *x22        *x53
     6  +coeff( 33)    *x22*x33*x42    
     7  +coeff( 34)        *x32*x41    
     8  +coeff( 35)    *x21    *x42    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 36)            *x43    
     1  +coeff( 37)        *x32    *x51
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)        *x31    *x52
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)*x11*x21*x31        
     6  +coeff( 42)    *x22*x32        
     7  +coeff( 43)        *x33*x41    
     8  +coeff( 44)*x11*x23            
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)    *x23*x31*x41    
     2  +coeff( 47)    *x23    *x42    
     3  +coeff( 48)    *x22*x31*x42    
     4  +coeff( 49)        *x33*x42    
     5  +coeff( 50)        *x32*x43    
c
      return
      end
      function x_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5099132E+01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.44253320E-02,-0.12225901E+00, 0.15517200E-01,-0.33401868E-02,
     +  0.39317943E-02, 0.23040159E-02, 0.78555588E-02,-0.23216021E-02,
     + -0.30154833E-02, 0.29447645E-04,-0.75768880E-02,-0.53209249E-04,
     +  0.20568150E-02,-0.33811177E-03, 0.37631457E-03,-0.93782484E-03,
     + -0.39172624E-02, 0.14250119E-02, 0.73244902E-02,-0.12919597E-02,
     + -0.69796555E-02,-0.44325562E-02,-0.49913400E-04, 0.17993834E-02,
     +  0.16563631E-02, 0.76464284E-03,-0.31073109E-03, 0.17960328E-02,
     +  0.75482894E-02, 0.13371488E-02, 0.49680634E-03,-0.20788494E-03,
     +  0.26327177E-03,-0.10616280E-02,-0.34956511E-02,-0.32693234E-02,
     +  0.52776828E-03,-0.17475937E-02, 0.25834923E-02, 0.42580641E-03,
     + -0.24305478E-04, 0.14993871E-03,-0.16240640E-03, 0.22012614E-03,
     + -0.46683391E-03,-0.61905565E-03, 0.15455272E-03, 0.13979011E-02,
     + -0.30348257E-02,-0.10682349E-03,-0.96881442E-03, 0.45388186E-03,
     +  0.43719533E-03,-0.57516016E-04, 0.26184120E-03,-0.10595947E-03,
     +  0.11973055E-02, 0.24271894E-02, 0.14061780E-02,-0.54181699E-03,
     + -0.92695293E-04,-0.10138775E-02,-0.10676851E-02, 0.52746566E-03,
     +  0.49130915E-03,-0.10605027E-03,-0.54315835E-04, 0.64390962E-03,
     + -0.11143882E-02,-0.17895451E-03,-0.26834811E-03, 0.10915480E-03,
     + -0.10256864E-03,-0.18699460E-03,-0.85170730E-04,-0.10714490E-02,
     + -0.10727458E-02,-0.62397320E-03,-0.23998593E-04, 0.39442792E-03,
     +  0.61732280E-03, 0.48391201E-03,-0.33263736E-04, 0.30783427E-04,
     +  0.13246144E-03, 0.20902583E-03,-0.16224227E-03,-0.35650967E-03,
     + -0.77632903E-04,-0.11523155E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23*x31        
     5  +coeff(  5)            *x41    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21        *x51
      x_sp_den  =x_sp_den  
     9  +coeff(  9)    *x23            
     1  +coeff( 10)        *x33        
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22*x33        
     4  +coeff( 13)        *x31        
     5  +coeff( 14)                *x51
     6  +coeff( 15)    *x22            
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x22*x31        
      x_sp_den  =x_sp_den  
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)*x11    *x33        
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)            *x42    
     8  +coeff( 26)*x11    *x31        
      x_sp_den  =x_sp_den  
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)    *x21*x32        
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)    *x24            
     4  +coeff( 31)        *x32        
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)    *x21        *x52
     7  +coeff( 34)*x11*x21    *x41    
     8  +coeff( 35)    *x22*x31*x41    
      x_sp_den  =x_sp_den  
     9  +coeff( 36)    *x22    *x42    
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)    *x24    *x41    
     4  +coeff( 40)*x11*x23*x31        
     5  +coeff( 41)    *x24*x32        
     6  +coeff( 42)*x11*x24*x31        
     7  +coeff( 43)            *x41*x51
     8  +coeff( 44)    *x22        *x51
      x_sp_den  =x_sp_den  
     9  +coeff( 45)    *x21    *x41*x51
     1  +coeff( 46)*x11*x21*x31        
     2  +coeff( 47)    *x22*x31    *x51
     3  +coeff( 48)    *x24*x31        
     4  +coeff( 49)    *x23    *x42    
     5  +coeff( 50)        *x31    *x51
     6  +coeff( 51)    *x22*x32        
     7  +coeff( 52)*x11    *x31*x41    
     8  +coeff( 53)*x11        *x42    
      x_sp_den  =x_sp_den  
     9  +coeff( 54)*x11        *x41*x51
     1  +coeff( 55)    *x22    *x41*x51
     2  +coeff( 56)*x12*x21            
     3  +coeff( 57)    *x21*x32*x41    
     4  +coeff( 58)    *x21*x31*x42    
     5  +coeff( 59)    *x21    *x43    
     6  +coeff( 60)    *x21    *x42*x51
     7  +coeff( 61)    *x21    *x41*x52
     8  +coeff( 62)*x11*x22*x31        
      x_sp_den  =x_sp_den  
     9  +coeff( 63)    *x23*x32        
     1  +coeff( 64)    *x23    *x41*x51
     2  +coeff( 65)*x11*x24            
     3  +coeff( 66)*x11    *x32*x41    
     4  +coeff( 67)*x11    *x31*x42    
     5  +coeff( 68)*x11*x23    *x41    
     6  +coeff( 69)    *x21*x31*x43    
     7  +coeff( 70)    *x21*x33    *x51
     8  +coeff( 71)*x11*x22*x32        
      x_sp_den  =x_sp_den  
     9  +coeff( 72)        *x33*x42    
     1  +coeff( 73)        *x32*x43    
     2  +coeff( 74)*x11*x22*x31    *x51
     3  +coeff( 75)    *x23*x33        
     4  +coeff( 76)    *x23*x32*x42    
     5  +coeff( 77)*x11*x24*x31*x41    
     6  +coeff( 78)*x11*x24    *x42    
     7  +coeff( 79)*x12                
     8  +coeff( 80)        *x32*x41    
      x_sp_den  =x_sp_den  
     9  +coeff( 81)        *x31*x42    
     1  +coeff( 82)            *x43    
     2  +coeff( 83)                *x53
     3  +coeff( 84)*x11*x21        *x51
     4  +coeff( 85)*x11    *x32        
     5  +coeff( 86)    *x21*x33        
     6  +coeff( 87)    *x21*x32    *x51
     7  +coeff( 88)    *x21*x31*x41*x51
     8  +coeff( 89)            *x41*x53
      x_sp_den  =x_sp_den  
     9  +coeff( 90)*x11*x21*x31*x41    
c
      return
      end
      function t_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1299281E+01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.38278807E-03, 0.57920333E-01,-0.15726447E-02, 0.52839830E-02,
     + -0.84911594E-02, 0.20430775E-02,-0.20745825E-02, 0.69967035E-03,
     + -0.43193316E-02, 0.17049501E-02, 0.15406407E-04,-0.10702632E-02,
     +  0.63118146E-03,-0.88353769E-03, 0.32775849E-03,-0.64303656E-03,
     +  0.17296628E-02,-0.40986235E-02, 0.53668139E-03, 0.41082497E-02,
     +  0.34785680E-04, 0.11220398E-04, 0.13954350E-02,-0.23777314E-03,
     +  0.14596934E-03, 0.16276147E-03, 0.96933916E-03,-0.64239599E-03,
     + -0.34259607E-02, 0.77963353E-03,-0.34256579E-03,-0.13932049E-03,
     +  0.19341843E-03, 0.19995958E-03,-0.12513054E-03, 0.26289807E-03,
     +  0.11367031E-02, 0.39426182E-03, 0.15864321E-02,-0.10637827E-04,
     + -0.29010963E-03, 0.47017162E-04, 0.25596419E-04, 0.66645596E-04,
     +  0.79506339E-04, 0.17325315E-03, 0.35616336E-03, 0.94257714E-03,
     + -0.41734797E-03,-0.30359224E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)            *x41    
     8  +coeff(  8)                *x51
      t_sp_den  =t_sp_den  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)        *x33        
     3  +coeff( 12)        *x31        
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x22    *x41    
      t_sp_den  =t_sp_den  
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)*x11    *x33        
     4  +coeff( 22)    *x22*x33        
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)        *x32        
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)*x11            *x51
      t_sp_den  =t_sp_den  
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)    *x21*x32        
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)    *x21        *x52
      t_sp_den  =t_sp_den  
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)    *x22*x31*x41    
     2  +coeff( 38)*x11*x22*x31        
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)*x11*x23*x31        
     5  +coeff( 41)            *x42    
     6  +coeff( 42)        *x31    *x51
     7  +coeff( 43)            *x43    
     8  +coeff( 44)    *x21    *x41*x51
      t_sp_den  =t_sp_den  
     9  +coeff( 45)            *x41*x52
     1  +coeff( 46)*x11*x21*x31        
     2  +coeff( 47)    *x22*x32        
     3  +coeff( 48)    *x22    *x42    
     4  +coeff( 49)    *x21*x31*x42    
     5  +coeff( 50)    *x21    *x43    
c
      return
      end
      function y_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.4588855E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19634292E-02, 0.38058439E-03, 0.68885721E-01, 0.16194781E-01,
     + -0.49838391E-02, 0.34862482E-02, 0.25289957E-02, 0.15941383E-01,
     +  0.17506860E-02, 0.58757880E-03, 0.41725058E-02,-0.71181399E-02,
     + -0.11395161E-02,-0.73854318E-02, 0.21156354E-02,-0.13827837E-03,
     +  0.33861708E-02,-0.82869432E-04,-0.10762605E-02,-0.30207934E-02,
     + -0.14043419E-02, 0.77270181E-03,-0.37248584E-03, 0.34029465E-03,
     + -0.96791954E-03,-0.65261562E-03,-0.39137793E-02,-0.33783889E-02,
     +  0.62367221E-03,-0.13264221E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sp_den  =y_sp_den  
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x22            
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22*x33        
     8  +coeff( 17)        *x31*x41    
      y_sp_den  =y_sp_den  
     9  +coeff( 18)        *x33        
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)        *x32        
     5  +coeff( 23)*x11                
     6  +coeff( 24)    *x21        *x51
     7  +coeff( 25)            *x43    
     8  +coeff( 26)    *x21    *x41*x51
      y_sp_den  =y_sp_den  
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)    *x22*x34        
c
      return
      end
      function p_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1018158E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38273025E-02, 0.42293309E-02,-0.29086664E-01,-0.99486127E-01,
     + -0.94531458E-02, 0.63561564E-02,-0.46381112E-02,-0.26626118E-01,
     + -0.47526490E-02, 0.15281579E-01, 0.27955591E-02,-0.11294666E-01,
     +  0.29239554E-02,-0.36923937E-02, 0.72236243E-03,-0.81865670E-03,
     +  0.75622724E-03, 0.46823695E-02, 0.14754985E-02,-0.10489315E-02,
     + -0.20131361E-03, 0.25625320E-03,-0.29286457E-03, 0.74808148E-03,
     +  0.34989545E-02, 0.30268068E-03, 0.31204644E-03, 0.26285637E-03,
     +  0.46175084E-03, 0.95121539E-03,-0.31224482E-02, 0.34947053E-02,
     +  0.16506153E-03,-0.69265882E-03,-0.20743592E-03, 0.34201427E-02,
     + -0.50781405E-03,-0.53505506E-03, 0.19434841E-03, 0.67640445E-04,
     +  0.50387805E-03, 0.32985708E-03,-0.72482170E-03, 0.14488838E-03,
     + -0.30852563E-03,-0.19936361E-03,-0.21586053E-03, 0.53773896E-03,
     + -0.95952541E-03, 0.71829936E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sp_den  =p_sp_den  
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)*x11        *x41    
      p_sp_den  =p_sp_den  
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)*x11                
     5  +coeff( 23)        *x32        
     6  +coeff( 24)    *x23            
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)        *x31*x42    
      p_sp_den  =p_sp_den  
     9  +coeff( 27)    *x21*x31    *x51
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)    *x21        *x52
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)    *x23*x32        
     7  +coeff( 34)    *x22*x31    *x52
     8  +coeff( 35)*x11*x23*x31        
      p_sp_den  =p_sp_den  
     9  +coeff( 36)    *x22*x33*x41    
     1  +coeff( 37)    *x23*x31    *x52
     2  +coeff( 38)            *x42    
     3  +coeff( 39)*x11    *x31        
     4  +coeff( 40)*x11            *x51
     5  +coeff( 41)    *x21*x32        
     6  +coeff( 42)        *x32*x41    
     7  +coeff( 43)            *x43    
     8  +coeff( 44)        *x32    *x51
      p_sp_den  =p_sp_den  
     9  +coeff( 45)            *x42*x51
     1  +coeff( 46)                *x53
     2  +coeff( 47)*x11*x22            
     3  +coeff( 48)*x11*x21*x31        
     4  +coeff( 49)    *x23*x31        
     5  +coeff( 50)    *x22*x32        
c
      return
      end
      function l_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1440277E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12519561E-01, 0.23541974E+00,-0.12954248E-01,-0.45322455E-02,
     + -0.12643890E-01,-0.30001398E-01, 0.64689242E-02,-0.65018674E-02,
     + -0.15220565E-01, 0.46881535E-02, 0.64432900E-02,-0.89311898E-02,
     +  0.22714017E-02, 0.13248716E-01, 0.19947471E-03, 0.65031543E-03,
     + -0.49766623E-02,-0.15063304E-02,-0.23969377E-02, 0.65365508E-02,
     + -0.13292078E-01, 0.12920807E-01,-0.39641676E-03, 0.42138421E-02,
     + -0.10221617E-02, 0.61659812E-03,-0.23591486E-02,-0.11921554E-01,
     +  0.17774695E-02,-0.96566119E-03,-0.43052729E-03, 0.64236310E-03,
     +  0.41653932E-03, 0.15949880E-02, 0.56485836E-02, 0.56194742E-02,
     + -0.11931218E-02, 0.24762515E-02, 0.41926284E-02,-0.12827913E-03,
     +  0.19442629E-02,-0.14401085E-03,-0.19344901E-03, 0.46430682E-03,
     + -0.43564744E-03, 0.25562034E-03, 0.88568270E-03, 0.16394302E-02,
     + -0.94241032E-03,-0.39855688E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)        *x31        
      l_sp_den  =l_sp_den  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x23            
     3  +coeff( 12)            *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22*x33        
     7  +coeff( 16)                *x51
     8  +coeff( 17)        *x31*x41    
      l_sp_den  =l_sp_den  
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)*x11    *x33        
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)        *x32        
     8  +coeff( 26)*x11            *x51
      l_sp_den  =l_sp_den  
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11    *x31        
     4  +coeff( 31)    *x22        *x51
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)    *x21    *x41*x51
     7  +coeff( 34)*x11*x21    *x41    
     8  +coeff( 35)    *x22*x31*x41    
      l_sp_den  =l_sp_den  
     9  +coeff( 36)    *x22    *x42    
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)*x11*x23*x31        
     5  +coeff( 41)*x11*x22*x33        
     6  +coeff( 42)        *x31    *x51
     7  +coeff( 43)                *x52
     8  +coeff( 44)            *x42*x51
      l_sp_den  =l_sp_den  
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)            *x41*x52
     2  +coeff( 47)*x11*x21*x31        
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x21*x31*x42    
     5  +coeff( 50)    *x21*x31    *x52
c
      return
      end
      function x_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2593375E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12533084E-01, 0.35400668E+00,-0.12047050E-01, 0.13457318E+00,
     + -0.94145797E-02, 0.22166895E-01,-0.49522705E-01, 0.30144000E-01,
     +  0.11024419E-01,-0.62445058E-02,-0.24964465E-01, 0.11850312E-01,
     + -0.35780850E-02, 0.19982955E-02,-0.42639030E-02, 0.18022319E-01,
     + -0.26111925E-01, 0.43587927E-02, 0.23181992E-01,-0.27734179E-02,
     +  0.12123093E-01,-0.25110787E-04,-0.75568701E-02,-0.22505692E-02,
     +  0.91248490E-02,-0.63715405E-02,-0.24354829E-01,-0.26714909E-02,
     + -0.45374250E-02, 0.18327229E-03, 0.31075656E-03, 0.30940250E-03,
     + -0.18299581E-02,-0.66836718E-02,-0.67948207E-03, 0.43443975E-03,
     + -0.92995662E-03, 0.29900763E-02, 0.95731327E-02, 0.85916044E-02,
     +  0.55638184E-02, 0.98303817E-02,-0.98995364E-03, 0.49574472E-04,
     + -0.71690656E-03,-0.35447738E-03, 0.16941425E-02,-0.11532099E-02,
     + -0.70092352E-02, 0.33394150E-02,-0.50176838E-02, 0.42687971E-02,
     + -0.64539636E-04,-0.15352995E-02, 0.68900394E-02,-0.13822166E-02,
     +  0.18618964E-02,-0.21982244E-02,-0.12738537E-02, 0.24629669E-03,
     +  0.75453060E-03, 0.27682995E-02,-0.14605648E-02,-0.13566180E-02,
     +  0.31584810E-03,-0.65742788E-03,-0.27991191E-02,-0.31363883E-02,
     +  0.46133070E-03,-0.15342018E-02, 0.55164341E-02, 0.20512522E-03,
     +  0.19934829E-02, 0.12737322E-02, 0.78372762E-03, 0.46687643E-03,
     +  0.37178986E-02, 0.17483124E-02,-0.60366339E-03,-0.13062058E-03,
     + -0.10641675E-02, 0.17757030E-03, 0.34064133E-03, 0.21871127E-03,
     + -0.23176147E-03, 0.27357636E-03, 0.44613660E-03, 0.29194387E-03,
     +  0.18331841E-03, 0.36441212E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)        *x31        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x23            
     4  +coeff( 13)                *x52
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x24*x31        
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)*x11    *x33        
     5  +coeff( 23)            *x42    
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x21*x32        
      x_sp_dex    =x_sp_dex    
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)    *x24            
     3  +coeff( 30)        *x33*x41    
     4  +coeff( 31)        *x32*x42    
     5  +coeff( 32)        *x31*x43    
     6  +coeff( 33)        *x32        
     7  +coeff( 34)        *x31*x41    
     8  +coeff( 35)            *x41*x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 36)*x11            *x51
     1  +coeff( 37)    *x21        *x52
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)*x11*x22    *x41    
     6  +coeff( 42)    *x23    *x42    
     7  +coeff( 43)*x11*x23*x31        
     8  +coeff( 44)    *x24*x32        
      x_sp_dex    =x_sp_dex    
     9  +coeff( 45)*x11*x24*x31        
     1  +coeff( 46)    *x22        *x51
     2  +coeff( 47)*x11*x21*x31        
     3  +coeff( 48)*x11*x23            
     4  +coeff( 49)    *x21*x31*x42    
     5  +coeff( 50)*x11*x22*x31        
     6  +coeff( 51)    *x24    *x41    
     7  +coeff( 52)    *x23*x32        
     8  +coeff( 53)    *x23*x31    *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 54)*x11*x24            
     1  +coeff( 55)    *x21*x31*x43    
     2  +coeff( 56)    *x21*x32*x43    
     3  +coeff( 57)    *x23*x32*x42    
     4  +coeff( 58)        *x31*x42    
     5  +coeff( 59)            *x43    
     6  +coeff( 60)                *x53
     7  +coeff( 61)    *x23        *x51
     8  +coeff( 62)    *x22*x32        
      x_sp_dex    =x_sp_dex    
     9  +coeff( 63)*x11    *x31*x41    
     1  +coeff( 64)*x11        *x42    
     2  +coeff( 65)*x12*x21            
     3  +coeff( 66)    *x21*x33        
     4  +coeff( 67)    *x21*x32*x41    
     5  +coeff( 68)    *x21    *x43    
     6  +coeff( 69)    *x21    *x41*x52
     7  +coeff( 70)*x11*x23    *x41    
     8  +coeff( 71)    *x21*x32*x42    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 72)*x12    *x32        
     1  +coeff( 73)*x11*x22*x31*x41    
     2  +coeff( 74)*x11*x22    *x42    
     3  +coeff( 75)        *x33*x41*x51
     4  +coeff( 76)    *x23*x32    *x51
     5  +coeff( 77)    *x23*x33*x41    
     6  +coeff( 78)    *x23*x31*x42*x51
     7  +coeff( 79)    *x21*x31    *x51
     8  +coeff( 80)        *x33        
      x_sp_dex    =x_sp_dex    
     9  +coeff( 81)        *x32*x41    
     1  +coeff( 82)        *x32    *x51
     2  +coeff( 83)            *x42*x51
     3  +coeff( 84)*x11*x21        *x51
     4  +coeff( 85)*x11    *x32        
     5  +coeff( 86)    *x22        *x52
     6  +coeff( 87)    *x21    *x42*x51
     7  +coeff( 88)    *x21        *x53
     8  +coeff( 89)*x11*x22        *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 90)*x11*x21*x31*x41    
c
      return
      end
      function t_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.5329418E+00/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.51955367E-03,-0.76543473E-01, 0.22676529E-01,-0.58806604E-02,
     +  0.10692037E-01,-0.22910980E-02, 0.27529830E-02, 0.14577387E-02,
     +  0.53592608E-02,-0.31887377E-02,-0.62952357E-04, 0.15092926E-02,
     + -0.11465288E-02,-0.12639518E-02, 0.75296982E-03,-0.22023206E-02,
     + -0.47079958E-02, 0.93584153E-04,-0.16359766E-02, 0.43698477E-04,
     + -0.27937215E-03,-0.20179787E-03,-0.19869623E-03, 0.36532752E-03,
     + -0.19832890E-03,-0.12453323E-02, 0.67388767E-03, 0.32475844E-02,
     +  0.35331622E-02, 0.10948117E-02, 0.22873213E-03,-0.63739461E-03,
     + -0.63252803E-04,-0.12733768E-03,-0.11472839E-03,-0.33504737E-03,
     + -0.35879406E-03, 0.88495249E-03,-0.40099627E-03,-0.18689181E-02,
     + -0.19156084E-02,-0.85684779E-03,-0.16891557E-02, 0.98947086E-04,
     + -0.53681922E-03, 0.11032104E-03,-0.26501785E-03,-0.45491473E-03,
     +  0.30541603E-03, 0.17046610E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)            *x41    
     8  +coeff(  8)*x11                
      t_sp_dex    =t_sp_dex    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)        *x33        
     3  +coeff( 12)        *x31        
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x23    *x41    
      t_sp_dex    =t_sp_dex    
     9  +coeff( 18)    *x22*x33        
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)    *x23*x33*x41    
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)            *x42    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x22*x31        
      t_sp_dex    =t_sp_dex    
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)    *x21    *x42    
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)    *x21        *x52
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)        *x33*x41*x51
     7  +coeff( 34)        *x32        
     8  +coeff( 35)*x11*x21            
      t_sp_dex    =t_sp_dex    
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)    *x21    *x41*x51
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)*x11*x21    *x41    
     4  +coeff( 40)    *x22*x31*x41    
     5  +coeff( 41)    *x22    *x42    
     6  +coeff( 42)*x11*x22    *x41    
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)*x11*x23*x31        
      t_sp_dex    =t_sp_dex    
     9  +coeff( 45)*x11*x22*x33        
     1  +coeff( 46)        *x32    *x51
     2  +coeff( 47)*x11*x21*x31        
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x21*x31*x42    
     5  +coeff( 50)    *x23        *x51
c
      return
      end
      function y_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 35)
      data ncoeff/ 34/
      data avdat/  0.1162915E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.34511033E-02,-0.54175854E-01,-0.77403724E-01, 0.89946790E-02,
     + -0.13845073E-03, 0.31306406E-02, 0.87730559E-02, 0.51222723E-01,
     +  0.74799657E-02,-0.32950060E-02,-0.36808524E-01,-0.85110748E-02,
     + -0.70555496E-03,-0.37135528E-02,-0.43588812E-02,-0.27240838E-02,
     + -0.22273779E-01,-0.76789870E-02, 0.42967056E-02, 0.21441057E-02,
     + -0.38331610E-02,-0.13359918E-02, 0.35893626E-02,-0.14184045E-02,
     +  0.62705792E-03,-0.87201659E-03,-0.10744992E-03, 0.18014018E-02,
     + -0.88856812E-03,-0.93057274E-03,-0.45545620E-03, 0.15947591E-02,
     + -0.58531447E-03,-0.22478884E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sp_dex    =y_sp_dex    
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)    *x22            
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x22    *x41    
      y_sp_dex    =y_sp_dex    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x23            
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)        *x31*x42    
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)        *x31    *x52
      y_sp_dex    =y_sp_dex    
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)                *x53
     3  +coeff( 30)    *x21*x31    *x51
     4  +coeff( 31)*x11*x21            
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)    *x22    *x41*x51
c
      return
      end
      function p_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.3732368E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12361675E-02,-0.13964322E-01,-0.31178396E-01,-0.36213265E-03,
     + -0.18492162E-03,-0.73283818E-02, 0.19769500E-03,-0.21486566E-02,
     +  0.14341930E-02, 0.98537700E-02, 0.16732384E-02, 0.69547544E-03,
     + -0.28554145E-02,-0.16423237E-02, 0.48538900E-03, 0.78543706E-03,
     +  0.43963359E-03, 0.10569214E-02, 0.26702663E-03,-0.77938894E-03,
     + -0.10583723E-02, 0.17803127E-03, 0.48962561E-03,-0.85494318E-03,
     + -0.16976449E-03, 0.38993625E-04,-0.18384440E-04, 0.10246831E-04,
     +  0.69549809E-04, 0.13019024E-03,-0.12726318E-03,-0.11469852E-03,
     + -0.10111983E-03, 0.54863503E-03, 0.76395477E-03, 0.44555299E-03,
     + -0.35117657E-03, 0.33928492E-04,-0.10264327E-03,-0.50745287E-03,
     +  0.50581672E-04,-0.14627352E-03, 0.73113137E-04, 0.39521090E-03,
     +  0.19848310E-03,-0.18130564E-03,-0.38538731E-03, 0.23658929E-03,
     + -0.19414879E-03, 0.22249417E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_dex    =p_sp_dex    
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x21            
     7  +coeff( 16)    *x22            
     8  +coeff( 17)    *x22*x31        
      p_sp_dex    =p_sp_dex    
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)    *x21*x31    *x51
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)*x11                
      p_sp_dex    =p_sp_dex    
     9  +coeff( 27)        *x32        
     1  +coeff( 28)        *x31*x41    
     2  +coeff( 29)*x11*x21            
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)            *x41*x52
     5  +coeff( 32)                *x53
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)    *x22    *x42    
      p_sp_dex    =p_sp_dex    
     9  +coeff( 36)    *x23        *x51
     1  +coeff( 37)    *x22*x31    *x51
     2  +coeff( 38)*x11    *x31        
     3  +coeff( 39)        *x32*x41    
     4  +coeff( 40)        *x31*x42    
     5  +coeff( 41)        *x32    *x51
     6  +coeff( 42)            *x42*x51
     7  +coeff( 43)*x11*x21    *x41    
     8  +coeff( 44)    *x22*x31*x41    
      p_sp_dex    =p_sp_dex    
     9  +coeff( 45)    *x21    *x43    
     1  +coeff( 46)    *x21    *x41*x52
     2  +coeff( 47)    *x23    *x42    
     3  +coeff( 48)    *x22*x31*x42    
     4  +coeff( 49)    *x23*x31    *x51
     5  +coeff( 50)*x11*x22    *x42*x51
c
      return
      end
      function l_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1857960E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.93412381E-02,-0.49589071E+00,-0.10015509E+00, 0.11955075E-01,
     + -0.43356311E-01, 0.70522562E-01,-0.25491158E-01,-0.16438061E-01,
     +  0.11725446E-01, 0.36080863E-01,-0.19132515E-01, 0.63467883E-02,
     +  0.43975124E-02, 0.35890320E-03, 0.54006758E-02,-0.14685604E-01,
     +  0.34312293E-01,-0.43255184E-02,-0.32835707E-01,-0.17042223E-04,
     + -0.36254787E-03,-0.15175828E-01,-0.19408267E-02, 0.17771372E-02,
     + -0.18216291E-02,-0.77036540E-02, 0.72720987E-02, 0.30120850E-01,
     +  0.16233501E-02,-0.59474621E-03, 0.27618138E-02,-0.94660238E-03,
     + -0.45887727E-03,-0.90486305E-02,-0.79539912E-02,-0.65670796E-02,
     + -0.13893206E-01, 0.66418224E-03,-0.81212551E-03, 0.12713921E-02,
     +  0.88754401E-03, 0.85500139E-03,-0.15224514E-02,-0.29483205E-02,
     + -0.28145709E-02, 0.37303290E-02, 0.28035131E-02, 0.13892052E-02,
     + -0.29032968E-02,-0.37701731E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23*x31        
      l_sp_dex    =l_sp_dex    
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)        *x31        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)*x11    *x33        
     3  +coeff( 21)    *x22*x33        
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)                *x52
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x22*x31        
      l_sp_dex    =l_sp_dex    
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)        *x32        
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)    *x22*x31*x41    
     8  +coeff( 35)    *x22    *x42    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)*x11*x23    *x41    
     3  +coeff( 39)*x11*x22*x33        
     4  +coeff( 40)    *x21    *x41*x51
     5  +coeff( 41)            *x42*x51
     6  +coeff( 42)    *x21        *x52
     7  +coeff( 43)*x11*x21*x31        
     8  +coeff( 44)*x11*x21    *x41    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 45)    *x22*x32        
     1  +coeff( 46)    *x21*x31*x42    
     2  +coeff( 47)    *x21    *x43    
     3  +coeff( 48)*x11*x23            
     4  +coeff( 49)*x11*x22*x31        
     5  +coeff( 50)    *x23*x32        
c
      return
      end
      function x_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5309471E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12893287E-01, 0.22888893E+00,-0.40939990E-02,-0.78634582E-02,
     +  0.13971300E+00,-0.67108935E-02, 0.25533039E-01,-0.33371106E-01,
     +  0.25367245E-01, 0.75834501E-02,-0.16837811E-01,-0.59653004E-02,
     +  0.97943339E-02,-0.28111837E-02, 0.93054036E-02,-0.18094990E-01,
     +  0.30165475E-02, 0.15801737E-01,-0.11530914E-02, 0.28161285E-03,
     +  0.89845108E-02, 0.27581042E-04,-0.57172128E-02, 0.94462623E-03,
     + -0.14730068E-02,-0.39741560E-02,-0.16679250E-01,-0.24546406E-02,
     + -0.24973282E-02, 0.56807948E-02,-0.10263722E-03, 0.14367092E-04,
     +  0.36488753E-04, 0.16078127E-03,-0.12614701E-02,-0.46415869E-02,
     + -0.97680918E-03, 0.44619581E-02, 0.18287767E-02,-0.52959885E-03,
     +  0.49065039E-02, 0.39529111E-02,-0.17016407E-03,-0.43547028E-03,
     +  0.23154361E-03,-0.77340129E-03, 0.41020283E-03, 0.18132540E-02,
     +  0.10069213E-02, 0.21917964E-02, 0.59514097E-02,-0.57236722E-03,
     +  0.27903463E-03, 0.74543245E-03,-0.13994939E-03, 0.87892741E-03,
     +  0.41646756E-04, 0.17919824E-02,-0.71879942E-03,-0.68231946E-03,
     + -0.75810991E-03,-0.53639226E-02,-0.21101958E-02, 0.63884858E-03,
     +  0.64281229E-03,-0.22110683E-02, 0.21433260E-02,-0.13252313E-02,
     + -0.70492405E-03, 0.20051065E-02,-0.35731014E-03, 0.10210664E-03,
     + -0.26660322E-03,-0.15932268E-02, 0.33178943E-03, 0.49397261E-02,
     +  0.47429712E-03, 0.35482224E-02, 0.69347589E-04,-0.82839317E-04,
     + -0.86562219E-03,-0.14651857E-02,-0.10040699E-02,-0.16755347E-03,
     +  0.15064316E-03,-0.37720590E-03,-0.16968350E-02, 0.27664044E-03,
     +  0.81931852E-04,-0.48958481E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x23*x31        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x23            
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x24*x31        
     2  +coeff( 20)    *x24        *x51
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)*x11    *x33        
     5  +coeff( 23)            *x42    
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)    *x21*x32        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)    *x24            
     3  +coeff( 30)    *x22*x31*x41    
     4  +coeff( 31)        *x33*x41    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)        *x31*x43    
     7  +coeff( 34)            *x41*x53
     8  +coeff( 35)        *x32        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 36)        *x31*x41    
     1  +coeff( 37)            *x41*x51
     2  +coeff( 38)    *x22*x31        
     3  +coeff( 39)    *x22        *x51
     4  +coeff( 40)    *x21        *x52
     5  +coeff( 41)    *x22    *x42    
     6  +coeff( 42)*x11*x22    *x41    
     7  +coeff( 43)    *x24*x32        
     8  +coeff( 44)*x11*x24*x31        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 45)*x11            *x51
     1  +coeff( 46)    *x21*x31    *x51
     2  +coeff( 47)                *x53
     3  +coeff( 48)*x11*x21    *x41    
     4  +coeff( 49)    *x23        *x51
     5  +coeff( 50)*x11*x22*x31        
     6  +coeff( 51)    *x23    *x42    
     7  +coeff( 52)*x11*x23*x31        
     8  +coeff( 53)*x11*x21*x33        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 54)    *x21*x31*x42*x52
     1  +coeff( 55)        *x31    *x51
     2  +coeff( 56)*x11*x21*x31        
     3  +coeff( 57)*x11*x21        *x51
     4  +coeff( 58)    *x22*x32        
     5  +coeff( 59)*x11    *x31*x41    
     6  +coeff( 60)*x11        *x42    
     7  +coeff( 61)*x11*x23            
     8  +coeff( 62)    *x21*x31*x42    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 63)    *x21    *x43    
     1  +coeff( 64)    *x21    *x42*x51
     2  +coeff( 65)    *x21    *x41*x52
     3  +coeff( 66)    *x24    *x41    
     4  +coeff( 67)    *x23*x32        
     5  +coeff( 68)*x11*x24            
     6  +coeff( 69)*x11*x23    *x41    
     7  +coeff( 70)    *x21*x31*x43    
     8  +coeff( 71)        *x33*x42    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 72)*x12*x23            
     1  +coeff( 73)    *x23*x32*x41    
     2  +coeff( 74)    *x21*x32*x43    
     3  +coeff( 75)    *x24*x33        
     4  +coeff( 76)    *x23*x32*x42    
     5  +coeff( 77)    *x23*x33    *x51
     6  +coeff( 78)    *x23*x33*x43    
     7  +coeff( 79)*x12                
     8  +coeff( 80)        *x33        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 81)        *x32*x41    
     1  +coeff( 82)        *x31*x42    
     2  +coeff( 83)            *x43    
     3  +coeff( 84)*x11    *x32        
     4  +coeff( 85)*x12*x21            
     5  +coeff( 86)    *x21*x33        
     6  +coeff( 87)    *x21*x32*x41    
     7  +coeff( 88)    *x21*x31*x41*x51
     8  +coeff( 89)*x12        *x41    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 90)    *x23    *x41*x51
c
      return
      end
      function t_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 50)
      data ncoeff/ 49/
      data avdat/ -0.1464176E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13733167E-02,-0.77285156E-01, 0.22559708E-01,-0.75042802E-02,
     +  0.10628263E-01,-0.25504604E-02, 0.25034426E-02, 0.14745987E-02,
     +  0.54326188E-02,-0.23983854E-02,-0.13735149E-02,-0.29246791E-02,
     + -0.11311965E-04, 0.12862440E-02, 0.80054620E-03, 0.46061138E-02,
     + -0.66689512E-03,-0.52072974E-02, 0.15928570E-03,-0.12901046E-02,
     +  0.13919073E-02, 0.11607127E-02,-0.11139541E-02, 0.78321836E-03,
     + -0.21985364E-02, 0.40533496E-02,-0.23245903E-03, 0.41303801E-03,
     +  0.29893365E-03,-0.15928413E-03, 0.34148496E-03,-0.38251985E-03,
     + -0.85669570E-03,-0.98237838E-03,-0.31689779E-03,-0.75883605E-03,
     + -0.18297185E-03, 0.13165358E-03,-0.13141455E-03,-0.27394592E-03,
     +  0.87654218E-04,-0.54463604E-03,-0.22274972E-03,-0.25798552E-03,
     +  0.19032022E-03, 0.24075124E-03,-0.10937754E-02,-0.30060930E-03,
     +  0.48147776E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)            *x41    
     8  +coeff(  8)*x11                
      t_sp_q3en   =t_sp_q3en   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23            
     4  +coeff( 13)        *x33        
     5  +coeff( 14)        *x31        
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)*x11    *x33        
     2  +coeff( 20)    *x23*x31*x41    
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)            *x42    
     5  +coeff( 23)    *x22*x31        
     6  +coeff( 24)    *x21*x32        
     7  +coeff( 25)    *x22    *x41    
     8  +coeff( 26)    *x21*x31*x41    
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 27)    *x22*x32        
     1  +coeff( 28)        *x32        
     2  +coeff( 29)*x11    *x31        
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)    *x21        *x52
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)*x11*x22    *x41    
     8  +coeff( 35)*x11*x21*x33        
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 36)*x11*x22*x33        
     1  +coeff( 37)*x11*x21            
     2  +coeff( 38)    *x22        *x51
     3  +coeff( 39)    *x21*x31    *x51
     4  +coeff( 40)            *x42*x51
     5  +coeff( 41)                *x53
     6  +coeff( 42)    *x22    *x42    
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)    *x22    *x41*x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 45)    *x21*x31*x41*x51
     1  +coeff( 46)*x11*x23            
     2  +coeff( 47)    *x23    *x42    
     3  +coeff( 48)        *x33*x41*x51
     4  +coeff( 49)    *x21*x31*x42*x52
c
      return
      end
      function y_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/  0.1532255E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.46881558E-02,-0.68816967E-01,-0.10893060E+00, 0.88319303E-02,
     + -0.13194732E-03, 0.39489758E-02, 0.11256397E-01, 0.64057887E-01,
     +  0.92552807E-02,-0.16867259E-02,-0.40107857E-01,-0.10595278E-01,
     + -0.18105975E-02,-0.23740327E-01,-0.93180258E-02, 0.49082730E-02,
     + -0.13632134E-02, 0.17139000E-03, 0.22068913E-02, 0.26863255E-03,
     + -0.11376198E-02, 0.76174433E-03,-0.96561312E-03,-0.40210653E-02,
     +  0.34937847E-02,-0.11034206E-02,-0.50869682E-02,-0.44515001E-03,
     + -0.23537860E-02,-0.85318030E-03,-0.31716824E-02, 0.21682014E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sp_q3en   =y_sp_q3en   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23            
     8  +coeff( 17)            *x43*x52
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 18)    *x21    *x43*x51
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)        *x31*x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)        *x31    *x52
     6  +coeff( 24)            *x41*x52
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)                *x53
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)*x11*x21            
     2  +coeff( 29)            *x45    
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)    *x21*x33*x41    
c
      return
      end
      function p_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.3825038E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12571997E-02,-0.14869470E-01,-0.31019928E-01, 0.12046281E-03,
     + -0.13227164E-02,-0.99744927E-02,-0.22747149E-02, 0.15301822E-02,
     +  0.11005373E-01, 0.18874052E-02, 0.77894435E-03,-0.44258712E-02,
     + -0.18055072E-02, 0.61042723E-03, 0.74921933E-03, 0.39515819E-03,
     + -0.18439409E-02,-0.29246049E-03, 0.16258193E-03, 0.22247669E-03,
     +  0.11431240E-02,-0.26677863E-03,-0.77551650E-03,-0.54332538E-03,
     + -0.66297717E-03, 0.56343852E-03,-0.41100398E-05, 0.30655166E-04,
     +  0.13693696E-03, 0.46361881E-03,-0.12142026E-03,-0.12982765E-03,
     +  0.69125748E-03, 0.40600667E-03, 0.32914781E-04, 0.54892604E-04,
     + -0.10486385E-03, 0.38092687E-04, 0.76659882E-04, 0.22458042E-03,
     + -0.42135372E-04, 0.62470484E-04, 0.89781766E-04,-0.14497824E-03,
     +  0.34284350E-03, 0.38303458E-03,-0.22543782E-03, 0.38281889E-03,
     + -0.31545511E-03, 0.25389605E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x31    *x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x21            
     6  +coeff( 15)    *x22            
     7  +coeff( 16)            *x42    
     8  +coeff( 17)    *x22    *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 18)    *x23    *x41*x51
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)        *x31*x42    
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)    *x22*x31    *x51
     8  +coeff( 26)    *x23*x31*x41    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 27)        *x31*x43*x51
     1  +coeff( 28)*x11                
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)    *x21*x31*x41    
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x11*x21        *x51
     6  +coeff( 33)    *x22    *x42    
     7  +coeff( 34)    *x23        *x51
     8  +coeff( 35)*x11*x21            
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 36)*x11    *x31        
     1  +coeff( 37)    *x21*x32        
     2  +coeff( 38)    *x21*x31    *x51
     3  +coeff( 39)        *x32    *x51
     4  +coeff( 40)        *x31*x41*x51
     5  +coeff( 41)    *x21        *x52
     6  +coeff( 42)        *x31    *x52
     7  +coeff( 43)*x11*x21    *x41    
     8  +coeff( 44)    *x23    *x41    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 45)    *x22*x31*x41    
     1  +coeff( 46)    *x21    *x43    
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)    *x23*x32        
     4  +coeff( 49)    *x22    *x42*x51
     5  +coeff( 50)    *x21*x33*x42    
c
      return
      end
      function l_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1020145E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.75171120E-02,-0.31906408E+00,-0.32790456E-01, 0.72003817E-02,
     + -0.37029110E-01, 0.45699514E-01,-0.10744084E-01, 0.23233010E-01,
     + -0.81914719E-02,-0.13014623E-01, 0.56859106E-02,-0.36630982E-02,
     + -0.10000300E-03, 0.32819579E-02,-0.23952932E-02, 0.34743743E-02,
     +  0.21099683E-01,-0.27746423E-02,-0.21611307E-01,-0.52023889E-02,
     + -0.43673315E-02,-0.18399906E-03,-0.99761318E-02,-0.71894709E-03,
     + -0.41211327E-02, 0.47456399E-02,-0.73281261E-02, 0.19096041E-01,
     + -0.43694354E-02,-0.62448252E-03,-0.28131891E-03,-0.76130993E-03,
     +  0.18549102E-02,-0.79713255E-03, 0.16840110E-02,-0.14474940E-02,
     + -0.20988530E-02,-0.88951411E-02,-0.13263919E-02, 0.72944246E-03,
     +  0.93801116E-03,-0.34567164E-03, 0.75092877E-03, 0.37625068E-03,
     +  0.47432727E-03,-0.16217309E-02, 0.26627802E-02, 0.21153332E-02,
     +  0.79080299E-03,-0.25078501E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
      l_sp_q3en   =l_sp_q3en   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)        *x33        
     5  +coeff( 14)        *x31        
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x21    *x42    
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11    *x33        
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x21*x32        
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 27)    *x22    *x41    
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)                *x52
     5  +coeff( 32)*x11*x21            
     6  +coeff( 33)*x11    *x31        
     7  +coeff( 34)    *x21*x31    *x51
     8  +coeff( 35)            *x42*x51
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)*x11*x22*x31        
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)*x11*x21*x33        
     4  +coeff( 40)        *x32        
     5  +coeff( 41)        *x31*x41    
     6  +coeff( 42)    *x21    *x41*x51
     7  +coeff( 43)        *x31*x41*x51
     8  +coeff( 44)    *x21        *x52
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 45)            *x41*x52
     1  +coeff( 46)    *x22*x32        
     2  +coeff( 47)    *x21*x31*x42    
     3  +coeff( 48)    *x21    *x43    
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)    *x23*x32        
c
      return
      end
      function x_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.8403291E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12273602E-01, 0.13044804E+00,-0.26176963E-02,-0.49084434E-02,
     +  0.18948768E+00, 0.16502773E-01,-0.20302564E-01, 0.23399411E-01,
     + -0.95903343E-02,-0.51259468E-02,-0.10442187E-01, 0.58470205E-02,
     + -0.42325621E-02, 0.62578772E-02, 0.21250530E-02,-0.12972950E-01,
     +  0.96369367E-02,-0.37109616E-03, 0.56716022E-02,-0.26004252E-02,
     + -0.11394504E-02, 0.73628873E-03,-0.17092583E-02,-0.11426355E-01,
     + -0.27503371E-02, 0.18337390E-02,-0.30304473E-02, 0.45810398E-02,
     +  0.11235618E-02, 0.17043514E-02,-0.60077064E-03,-0.86314254E-03,
     +  0.27750784E-02,-0.29417558E-02,-0.10534186E-02,-0.63619827E-03,
     +  0.61497581E-03, 0.24211842E-02,-0.23265371E-03,-0.77648193E-03,
     +  0.11787325E-02, 0.37684438E-02, 0.33442245E-02, 0.42439373E-02,
     +  0.64359978E-03,-0.36574056E-03,-0.20500182E-03, 0.12133745E-03,
     + -0.14689969E-03,-0.81439869E-03,-0.65681648E-04, 0.66151272E-03,
     +  0.35442036E-03, 0.12443725E-02,-0.45963819E-03,-0.44034340E-03,
     + -0.22931595E-02, 0.13744447E-02, 0.41602523E-03,-0.12347990E-02,
     +  0.46412728E-03,-0.81792875E-03, 0.18921080E-02, 0.30661150E-03,
     + -0.34486915E-04,-0.57947863E-03,-0.10189399E-03,-0.57792681E-03,
     +  0.27294651E-04, 0.45233817E-04,-0.54313010E-03,-0.45149741E-03,
     + -0.29987207E-03,-0.33053386E-03, 0.11717447E-03,-0.79293823E-03,
     + -0.85797632E-03, 0.28442676E-03, 0.10314799E-02, 0.19015568E-03,
     + -0.56003570E-03, 0.37380599E-03, 0.18476132E-03,-0.31518724E-03,
     + -0.52119768E-03, 0.38377833E-03, 0.66047709E-03,-0.52220677E-03,
     +  0.12761797E-02,-0.33626918E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x23            
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x23    *x41    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 18)    *x24*x31        
     1  +coeff( 19)    *x23*x31*x41    
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)*x11*x22            
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 27)    *x24            
     1  +coeff( 28)    *x23*x31        
     2  +coeff( 29)*x11*x22*x31        
     3  +coeff( 30)    *x23*x32        
     4  +coeff( 31)        *x32        
     5  +coeff( 32)*x11    *x31        
     6  +coeff( 33)    *x22*x31        
     7  +coeff( 34)    *x21*x32        
     8  +coeff( 35)    *x21*x31    *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 36)    *x21        *x52
     1  +coeff( 37)                *x53
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)        *x31    *x51
     4  +coeff( 40)            *x42*x51
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x22*x31*x41    
     7  +coeff( 43)    *x22    *x42    
     8  +coeff( 44)    *x23    *x42    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 45)    *x22*x31*x41*x51
     1  +coeff( 46)*x11*x23*x31        
     2  +coeff( 47)    *x24*x32        
     3  +coeff( 48)        *x33*x41*x51
     4  +coeff( 49)        *x32    *x51
     5  +coeff( 50)        *x31*x41*x51
     6  +coeff( 51)            *x41*x52
     7  +coeff( 52)*x11*x21*x31        
     8  +coeff( 53)    *x23        *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 54)    *x22*x32        
     1  +coeff( 55)    *x22        *x52
     2  +coeff( 56)*x11*x23            
     3  +coeff( 57)    *x21*x31*x42    
     4  +coeff( 58)    *x21    *x42*x51
     5  +coeff( 59)    *x21    *x41*x52
     6  +coeff( 60)    *x24    *x41    
     7  +coeff( 61)    *x23*x31    *x51
     8  +coeff( 62)*x11*x24            
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 63)    *x21*x31*x43    
     1  +coeff( 64)    *x23*x31*x41*x51
     2  +coeff( 65)    *x21*x33*x42    
     3  +coeff( 66)    *x21*x32*x43    
     4  +coeff( 67)    *x23*x32*x42    
     5  +coeff( 68)*x12    *x31*x42*x52
     6  +coeff( 69)*x11            *x51
     7  +coeff( 70)        *x33        
     8  +coeff( 71)        *x31*x42    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 72)            *x43    
     1  +coeff( 73)*x11    *x31*x41    
     2  +coeff( 74)*x11        *x42    
     3  +coeff( 75)*x12*x21            
     4  +coeff( 76)    *x21*x32*x41    
     5  +coeff( 77)    *x21    *x43    
     6  +coeff( 78)    *x21*x32    *x51
     7  +coeff( 79)    *x21*x31*x41*x51
     8  +coeff( 80)    *x21*x31    *x52
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 81)    *x24        *x51
     1  +coeff( 82)            *x42*x52
     2  +coeff( 83)            *x41*x53
     3  +coeff( 84)    *x23        *x52
     4  +coeff( 85)    *x22*x32*x41    
     5  +coeff( 86)    *x22    *x43    
     6  +coeff( 87)    *x22    *x42*x51
     7  +coeff( 88)*x11*x23    *x41    
     8  +coeff( 89)    *x21*x32*x42    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 90)    *x21    *x43*x51
c
      return
      end
      function t_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 50)
      data ncoeff/ 49/
      data avdat/ -0.3508774E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17663747E-02,-0.44094801E-01, 0.62799529E-03, 0.13000639E-02,
     +  0.58323272E-01,-0.40222388E-02,-0.49306406E-02, 0.56224796E-02,
     + -0.13415082E-02, 0.28682426E-02, 0.95667021E-03,-0.14941937E-02,
     +  0.32100070E-03, 0.63207152E-03, 0.25206909E-03, 0.39493822E-03,
     +  0.39895141E-03,-0.38544793E-03,-0.27807599E-02, 0.18868002E-03,
     + -0.10166099E-03,-0.47577309E-03,-0.96374721E-03, 0.12349308E-02,
     + -0.42585118E-03,-0.61352638E-03,-0.64278423E-03,-0.24572629E-03,
     + -0.21378204E-04, 0.50766401E-04, 0.62245020E-03,-0.48178926E-03,
     + -0.87929649E-04,-0.26285677E-03, 0.26113042E-03, 0.75771517E-04,
     + -0.81139915E-04,-0.24852986E-03, 0.39750357E-04, 0.22681794E-03,
     +  0.13374377E-02, 0.10461683E-03,-0.13320941E-03,-0.79800724E-04,
     +  0.51000324E-03, 0.49800897E-03, 0.44390385E-03,-0.37207382E-03,
     + -0.29792404E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      t_sp_q3m    =t_sp_q3m    
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)*x11                
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)        *x32        
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)                *x53
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)*x11    *x31        
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x22*x31        
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)        *x31*x41*x51
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)    *x21*x33    *x51
     8  +coeff( 35)    *x23*x33*x41    
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 36)            *x42    
     1  +coeff( 37)        *x31    *x51
     2  +coeff( 38)            *x41*x51
     3  +coeff( 39)*x11*x21            
     4  +coeff( 40)    *x21*x32        
     5  +coeff( 41)    *x21*x31*x41    
     6  +coeff( 42)    *x22        *x51
     7  +coeff( 43)        *x32    *x51
     8  +coeff( 44)    *x21        *x52
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 45)    *x21*x31*x41*x51
     1  +coeff( 46)    *x21    *x42*x51
     2  +coeff( 47)        *x31*x41*x52
     3  +coeff( 48)*x11*x22    *x41    
     4  +coeff( 49)*x11*x22*x33        
c
      return
      end
      function y_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/  0.1921629E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.60658357E-02,-0.83353668E-01,-0.14145458E+00, 0.76251365E-02,
     +  0.74676709E-03, 0.43711658E-02, 0.12176132E-01, 0.72774135E-01,
     +  0.11078024E-01,-0.33841706E-02,-0.50121110E-01,-0.12370532E-01,
     + -0.13160720E-02,-0.27205605E-01,-0.10954374E-01, 0.54672887E-02,
     + -0.34789492E-02,-0.37155934E-02, 0.50716302E-02,-0.55823894E-02,
     + -0.64336816E-02, 0.23754812E-02,-0.15850604E-02, 0.11240854E-02,
     +  0.91702677E-03, 0.25806639E-02,-0.23471620E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sp_q3m    =y_sp_q3m    
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23            
     8  +coeff( 17)            *x43    
      y_sp_q3m    =y_sp_q3m    
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)        *x31*x41    
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)    *x22            
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)    *x21*x31*x41    
      y_sp_q3m    =y_sp_q3m    
     9  +coeff( 27)    *x22*x31    *x51
c
      return
      end
      function p_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1131117E-03/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.75422058E-05, 0.45963706E-03, 0.23604159E-02,-0.23313195E-02,
     + -0.14510915E-02, 0.49551803E-03,-0.43638371E-03,-0.90506743E-04,
     +  0.89831662E-03,-0.36980529E-03,-0.43936842E-03, 0.54135645E-03,
     + -0.58107981E-02,-0.39715593E-03, 0.22353579E-02, 0.42041653E-03,
     + -0.21924348E-04, 0.59144979E-03, 0.96889744E-05, 0.14985830E-03,
     +  0.16951758E-02,-0.26954673E-03,-0.22586751E-03,-0.29342224E-04,
     + -0.84432861E-03, 0.45320476E-05, 0.14906946E-04,-0.85688560E-04,
     + -0.69377034E-04,-0.22461424E-03, 0.12402377E-04, 0.11455990E-04,
     + -0.77109917E-06, 0.16839844E-03, 0.61953790E-03, 0.75554723E-04,
     + -0.77217199E-04, 0.19240644E-03, 0.10067515E-04,-0.19884710E-02,
     + -0.38903664E-03, 0.86705311E-03, 0.37040751E-03, 0.33692471E-03,
     + -0.43185847E-03, 0.11862550E-03, 0.10368069E-03, 0.18075616E-03,
     +  0.94037583E-04,-0.37812148E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_q3m    =p_sp_q3m    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21*x31    *x51
     8  +coeff( 17)        *x32    *x51
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)        *x31*x41*x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)    *x22*x31    *x51
     6  +coeff( 24)        *x33    *x51
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)        *x32*x41*x51
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 27)        *x31*x42*x51
     1  +coeff( 28)            *x43*x51
     2  +coeff( 29)        *x31    *x53
     3  +coeff( 30)            *x41*x53
     4  +coeff( 31)*x11*x21        *x52
     5  +coeff( 32)    *x22*x33        
     6  +coeff( 33)    *x23        *x52
     7  +coeff( 34)    *x22*x31    *x52
     8  +coeff( 35)    *x22        *x53
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 36)    *x22*x33    *x51
     1  +coeff( 37)        *x33*x42*x51
     2  +coeff( 38)    *x22*x31    *x53
     3  +coeff( 39)        *x33    *x53
     4  +coeff( 40)        *x31    *x51
     5  +coeff( 41)    *x23            
     6  +coeff( 42)    *x22*x31        
     7  +coeff( 43)            *x43    
     8  +coeff( 44)        *x31    *x52
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 45)    *x21    *x41*x52
     1  +coeff( 46)*x11*x21            
     2  +coeff( 47)    *x21*x31*x41    
     3  +coeff( 48)        *x31*x42    
     4  +coeff( 49)                *x53
     5  +coeff( 50)    *x23    *x41    
c
      return
      end
      function l_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1159556E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49924E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.60213422E-02,-0.31925532E+00,-0.32573752E-01, 0.72010113E-02,
     + -0.40832613E-01, 0.46016116E-01,-0.10884196E-01, 0.23447765E-01,
     + -0.13712756E-01, 0.57513826E-02,-0.44970149E-02,-0.50706188E-02,
     + -0.13209859E-03, 0.33177894E-02,-0.25547813E-02,-0.12534508E-02,
     +  0.34683079E-02, 0.20860638E-01,-0.27645621E-02,-0.21697763E-01,
     + -0.39658328E-02,-0.34407524E-02,-0.11444580E-03,-0.99001583E-02,
     + -0.77429024E-03,-0.37733994E-02, 0.47934894E-02,-0.66977073E-02,
     +  0.18972412E-01,-0.17634231E-03,-0.44428553E-02, 0.18164770E-02,
     + -0.99064666E-03,-0.76049584E-03, 0.19707887E-02,-0.14846558E-02,
     + -0.21891538E-02,-0.88612279E-02, 0.23486481E-03, 0.46898800E-03,
     + -0.57188625E-03,-0.68722526E-03, 0.90630539E-03, 0.53049013E-03,
     + -0.10189697E-02,-0.12126836E-02, 0.26121391E-02, 0.21304642E-02,
     +  0.81311603E-03,-0.25773954E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
      l_sp_q3m    =l_sp_q3m    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x33        
     5  +coeff( 14)        *x31        
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)                *x52
     8  +coeff( 17)*x11        *x41    
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11    *x33        
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x22*x31        
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)        *x33    *x51
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)*x11    *x31        
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)            *x42*x51
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)*x11*x22*x31        
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)*x11*x23*x31        
     4  +coeff( 40)        *x32        
     5  +coeff( 41)        *x31    *x51
     6  +coeff( 42)*x11*x21            
     7  +coeff( 43)        *x31*x41*x51
     8  +coeff( 44)            *x41*x52
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 45)*x11*x21*x31        
     1  +coeff( 46)    *x22*x32        
     2  +coeff( 47)    *x21*x31*x42    
     3  +coeff( 48)    *x21    *x43    
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)    *x23*x32        
c
      return
      end
      function x_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2075658E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.42139608E-02, 0.89980438E-01,-0.22277089E-02,-0.40495647E-02,
     +  0.32268015E+00, 0.29269811E-01,-0.21804249E-01,-0.57938439E-02,
     +  0.14937590E-01,-0.16104972E-01,-0.83348351E-02, 0.48202649E-02,
     + -0.49079233E-02, 0.59039323E-02, 0.32627431E-02,-0.13080277E-01,
     + -0.17854536E-02, 0.84251974E-03, 0.24469811E-02,-0.13450335E-02,
     + -0.11307781E-01,-0.43073855E-02,-0.11023757E-02,-0.49954788E-02,
     +  0.22482059E-02, 0.78146271E-02,-0.18724386E-05, 0.16255222E-02,
     + -0.21728585E-02,-0.52767806E-03,-0.70959248E-03,-0.15466418E-02,
     +  0.14499387E-02, 0.36006114E-02,-0.28416084E-02,-0.15517152E-02,
     + -0.15354119E-02, 0.85646426E-03,-0.13774537E-02, 0.18253918E-02,
     + -0.18628601E-02, 0.11554157E-03, 0.19250551E-03,-0.33367015E-03,
     + -0.40944541E-03, 0.35681625E-03,-0.67171024E-03, 0.38837807E-02,
     +  0.35671510E-02,-0.17142274E-02, 0.40437799E-03, 0.22941292E-02,
     +  0.78184856E-03, 0.11223175E-02, 0.57462822E-02, 0.33911683E-02,
     +  0.79644274E-03, 0.93109888E-03,-0.87029045E-03,-0.42011001E-03,
     +  0.12346659E-02, 0.16581341E-03, 0.24081633E-03, 0.29875722E-03,
     + -0.89170979E-04, 0.14582820E-02, 0.31683149E-03,-0.19984171E-03,
     + -0.75824122E-03,-0.46795156E-03, 0.14831836E-02,-0.69109129E-03,
     +  0.23335210E-03, 0.12754354E-03, 0.60650532E-03,-0.58073405E-03,
     + -0.30686791E-03, 0.94776251E-03,-0.34508607E-03, 0.10455906E-02,
     + -0.31976899E-03,-0.45516514E-03, 0.13611853E-03, 0.21185954E-02,
     +  0.31260252E-02,-0.87914924E-03, 0.42171916E-03, 0.25071390E-02,
     +  0.20942844E-04,-0.53207776E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x23            
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)            *x41*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)    *x24            
     7  +coeff( 25)                *x53
     8  +coeff( 26)    *x23    *x41    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 27)        *x33*x41    
     1  +coeff( 28)    *x23*x32        
     2  +coeff( 29)        *x31*x41    
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)    *x23*x31        
     8  +coeff( 35)    *x21*x32        
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 36)        *x31*x41*x51
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x22        *x52
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)    *x24        *x51
     6  +coeff( 42)        *x32    *x53
     7  +coeff( 43)*x11*x22*x33        
     8  +coeff( 44)        *x32        
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 45)        *x32    *x51
     1  +coeff( 46)*x11*x21*x31        
     2  +coeff( 47)    *x23        *x51
     3  +coeff( 48)    *x22*x31*x41    
     4  +coeff( 49)    *x22    *x42    
     5  +coeff( 50)    *x21*x31*x42    
     6  +coeff( 51)    *x21*x32    *x51
     7  +coeff( 52)    *x21    *x42*x51
     8  +coeff( 53)*x11*x22*x31        
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 54)            *x42*x52
     1  +coeff( 55)    *x23*x31*x41    
     2  +coeff( 56)    *x23    *x42    
     3  +coeff( 57)    *x23*x31    *x51
     4  +coeff( 58)    *x23    *x41*x51
     5  +coeff( 59)    *x23        *x52
     6  +coeff( 60)    *x24*x32        
     7  +coeff( 61)    *x23*x31*x41*x51
     8  +coeff( 62)    *x21*x33*x42    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 63)    *x21*x33*x41*x51
     1  +coeff( 64)        *x33*x41*x52
     2  +coeff( 65)*x11            *x51
     3  +coeff( 66)    *x22*x32        
     4  +coeff( 67)    *x22    *x41*x51
     5  +coeff( 68)*x11*x23            
     6  +coeff( 69)    *x21*x32*x41    
     7  +coeff( 70)    *x21    *x43    
     8  +coeff( 71)    *x21*x31*x41*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 72)    *x24    *x41    
     1  +coeff( 73)*x11*x22        *x51
     2  +coeff( 74)        *x33    *x51
     3  +coeff( 75)        *x31*x41*x52
     4  +coeff( 76)*x11*x24            
     5  +coeff( 77)    *x22*x31*x41*x51
     6  +coeff( 78)    *x22    *x42*x51
     7  +coeff( 79)*x11*x23    *x41    
     8  +coeff( 80)    *x21*x31*x43    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 81)*x11*x22    *x41*x51
     1  +coeff( 82)        *x31*x42*x52
     2  +coeff( 83)*x12*x23            
     3  +coeff( 84)    *x24*x31*x41*x51
     4  +coeff( 85)    *x23*x32*x42    
     5  +coeff( 86)    *x23*x32    *x52
     6  +coeff( 87)*x11*x23    *x42*x51
     7  +coeff( 88)    *x23*x33*x43    
     8  +coeff( 89)*x12                
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 90)        *x32*x41    
c
      return
      end
      function t_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/ -0.1790809E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80691493E-03,-0.20141762E-01, 0.12556766E-03, 0.38378502E-03,
     +  0.10567356E+00, 0.50261701E-02,-0.98503092E-02,-0.77012519E-03,
     + -0.12158006E-02, 0.17037536E-02, 0.27469979E-03,-0.67430094E-03,
     + -0.41637840E-03, 0.81528473E-03, 0.20065240E-03,-0.71514642E-03,
     + -0.41204452E-03, 0.17509110E-03,-0.32068870E-04,-0.10946806E-02,
     + -0.80432778E-03,-0.10874136E-02, 0.49973640E-03,-0.94117998E-03,
     +  0.12294644E-03,-0.65313434E-04, 0.31213555E-03,-0.97977558E-04,
     + -0.67913439E-03, 0.17499617E-03,-0.19675968E-03,-0.28928765E-03,
     + -0.10513500E-02,-0.14723407E-03, 0.39189649E-03,-0.25768270E-03,
     +  0.57898310E-03,-0.25768863E-03, 0.48656148E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)        *x32        
     7  +coeff( 16)            *x42    
     8  +coeff( 17)            *x41*x51
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22*x32    *x51
     8  +coeff( 26)        *x33*x41*x51
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 27)    *x22        *x53
     1  +coeff( 28)        *x31    *x51
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)    *x22        *x51
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)        *x32    *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x22    *x42    
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 36)    *x22        *x52
     1  +coeff( 37)            *x42*x52
     2  +coeff( 38)    *x23*x32        
     3  +coeff( 39)        *x33*x41*x52
c
      return
      end
      function y_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1389974E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.36890572E-02,-0.62007837E-01,-0.11466108E+00, 0.43258271E-02,
     +  0.13858142E-02, 0.29345441E-02, 0.59061074E-02, 0.47140896E-01,
     +  0.79286220E-02,-0.35253770E-02,-0.38663127E-01,-0.90262387E-02,
     + -0.66326797E-03,-0.63255509E-04,-0.18547993E-01,-0.82708504E-02,
     + -0.23771771E-02, 0.15403978E-02, 0.87503099E-03, 0.39548627E-02,
     + -0.41948915E-02, 0.36893964E-02,-0.60936157E-02, 0.14534460E-02,
     + -0.10639309E-02, 0.70148101E-03, 0.21299832E-02,-0.21693760E-02,
     +  0.53545088E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)            *x43    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 18)    *x22            
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x23            
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)*x11        *x41    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x22*x31    *x51
     2  +coeff( 29)*x11*x22            
c
      return
      end
      function p_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3799357E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10185717E-02, 0.23024881E-03, 0.19001111E-01, 0.26572345E-01,
     + -0.29566197E-02, 0.40470198E-03, 0.11489336E-01,-0.86474430E-03,
     + -0.11815556E-02, 0.30839490E-02,-0.51671951E-02,-0.21605507E-01,
     + -0.27896108E-02, 0.11819835E-02, 0.76766140E-02, 0.26541941E-02,
     +  0.46522855E-02,-0.14931314E-02, 0.16810346E-02, 0.14531875E-02,
     +  0.25326360E-03,-0.90779940E-03, 0.15811183E-02, 0.72428264E-03,
     +  0.53986232E-03,-0.18688942E-03, 0.16374787E-03,-0.20550133E-03,
     +  0.90165435E-04,-0.10583092E-03, 0.50050864E-03,-0.20233673E-03,
     +  0.17674344E-03,-0.83482789E-03,-0.12133226E-02,-0.60027489E-03,
     +  0.46868441E-04, 0.84449374E-03,-0.29694245E-03,-0.42111584E-03,
     +  0.52986940E-03, 0.57326548E-03,-0.32058070E-03,-0.21423641E-04,
     +  0.20509874E-03,-0.87125496E-04,-0.13084772E-03, 0.14919516E-03,
     + -0.13613485E-03,-0.39601963E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)            *x41*x52
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)    *x22            
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)                *x53
     8  +coeff( 26)        *x32        
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)            *x41*x53
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)*x11                
     2  +coeff( 38)        *x31*x42    
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)    *x22*x31    *x51
     6  +coeff( 42)    *x22    *x41*x51
     7  +coeff( 43)        *x31    *x53
     8  +coeff( 44)*x11    *x31        
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 45)        *x32*x41    
     1  +coeff( 46)        *x32    *x51
     2  +coeff( 47)    *x21        *x52
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x21*x31*x42    
     5  +coeff( 50)    *x21    *x43    
c
      return
      end
      function l_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1733875E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20571414E-03,-0.31879425E+00,-0.32285608E-01, 0.71961000E-02,
     + -0.41669268E-01, 0.46736799E-01,-0.10679879E-01, 0.23553604E-01,
     + -0.76795043E-02,-0.14816926E-01, 0.60371682E-02,-0.44701854E-02,
     + -0.17441188E-03, 0.34181096E-02,-0.14449595E-02, 0.34649372E-02,
     +  0.19395325E-01,-0.21886462E-01,-0.47372016E-02,-0.42258902E-02,
     + -0.22093549E-03,-0.75465147E-02, 0.17620169E-02,-0.37407249E-02,
     +  0.34965819E-02,-0.67908848E-02, 0.17823942E-01,-0.27273847E-02,
     + -0.45401533E-02,-0.22748769E-03,-0.24118652E-02,-0.67220832E-03,
     +  0.11238625E-02,-0.96499885E-03,-0.13342447E-02,-0.19400465E-03,
     + -0.22570128E-02,-0.62655387E-02,-0.20925235E-03,-0.73951500E-03,
     +  0.22923051E-03,-0.53414452E-03,-0.24851135E-03,-0.11073326E-02,
     + -0.98859624E-03, 0.23657898E-02, 0.11795970E-02,-0.93682355E-03,
     + -0.53325924E-03, 0.99828793E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)        *x33        
     5  +coeff( 14)        *x31        
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x21    *x42    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)            *x41*x53
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)*x11    *x31        
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)    *x21*x32        
     8  +coeff( 26)    *x22    *x41    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)*x11*x22            
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)            *x42*x53
     4  +coeff( 31)            *x41*x51
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)    *x21        *x52
     8  +coeff( 35)*x11*x21    *x41    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 36)        *x33    *x51
     1  +coeff( 37)*x11*x22*x31        
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x21*x33    *x51
     4  +coeff( 40)    *x21    *x41*x53
     5  +coeff( 41)*x11*x23*x31        
     6  +coeff( 42)        *x31    *x51
     7  +coeff( 43)*x11*x21            
     8  +coeff( 44)    *x21*x31    *x51
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 45)    *x21    *x41*x51
     1  +coeff( 46)            *x42*x51
     2  +coeff( 47)                *x53
     3  +coeff( 48)*x11*x21*x31        
     4  +coeff( 49)    *x22*x32        
     5  +coeff( 50)    *x21*x31*x42    
c
      return
      end
      function x_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/  0.1191852E+00/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53854575E-02,-0.20122120E-01,-0.34139078E-01,-0.19767629E-04,
     + -0.44410535E-04, 0.69190661E-04,-0.14711030E-04, 0.34262364E-04,
     + -0.81985609E-05, 0.30482637E-04, 0.20698608E-04, 0.72317953E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      x_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)        *x31*x41    
      x_sp_sen    =x_sp_sen    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x24*x31        
     3  +coeff( 12)    *x21    *x41    
c
      return
      end
      function t_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.9635380E-01/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.26645276E-02,-0.61499566E-03,-0.27497200E-01,-0.23944663E-03,
     + -0.46063139E-03,-0.19531260E-04, 0.55534218E-03,-0.29690634E-03,
     + -0.23944273E-03, 0.34463781E-03, 0.23971358E-03,-0.34391440E-04,
     + -0.18750248E-03,-0.88656489E-04,-0.37312122E-05, 0.13751909E-03,
     + -0.82598999E-04, 0.35949764E-04,-0.17418295E-03, 0.26077038E-03,
     +  0.21575636E-03, 0.39077870E-04, 0.46370416E-04, 0.38886163E-04,
     + -0.32305634E-04, 0.48912538E-04, 0.67059787E-04,-0.64894986E-04,
     +  0.42631425E-04, 0.32830692E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)        *x31*x41    
      t_sp_sen    =t_sp_sen    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)*x11                
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)        *x32        
      t_sp_sen    =t_sp_sen    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)        *x33    *x51
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)*x11*x21    *x41    
      t_sp_sen    =t_sp_sen    
     9  +coeff( 27)    *x22*x32        
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)    *x22*x32    *x51
     3  +coeff( 30)*x11*x23*x31        
c
      return
      end
      function y_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/  0.5446269E-03/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.57491125E-03,-0.37380003E-05, 0.66430971E-01, 0.49871812E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x51 = x5
c
c                  function
c
      y_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)                *x51
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x11                
c
      return
      end
      function p_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.2639816E-03/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53959811E-03, 0.50686941E-01,-0.12283202E-02, 0.18280916E-03,
     + -0.17048563E-03,-0.12519833E-03, 0.20325009E-03,-0.63584105E-03,
     +  0.40009830E-03,-0.96435237E-04, 0.61855353E-04, 0.15948413E-03,
     + -0.22801322E-03,-0.18712114E-03, 0.28682430E-03,-0.56707431E-04,
     +  0.96763972E-04,-0.64342792E-04, 0.78036377E-04,-0.41379710E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      p_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23*x31        
     5  +coeff(  5)            *x41    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sp_sen    =p_sp_sen    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)        *x31        
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x22*x31        
      p_sp_sen    =p_sp_sen    
     9  +coeff( 18)    *x21*x32        
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11    *x33        
c
      return
      end
      function l_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1355398E-02/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11827123E-02, 0.17550093E-02, 0.30392946E-02,-0.17190408E-02,
     + -0.43993982E-05,-0.45777764E-03, 0.29875305E-05, 0.17915127E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)                *x51
c
      return
      end
      function x_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1800067E+00/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.70729782E-02,-0.20640416E-01,-0.47987439E-01,-0.84801111E-03,
     +  0.98780240E-03, 0.15237500E-02, 0.29527966E-03, 0.91019855E-03,
     +  0.42768457E-03,-0.24675543E-03,-0.44390236E-03, 0.25655574E-03,
     + -0.40504584E-03, 0.39847291E-04, 0.14348759E-03,-0.39392943E-03,
     +  0.15415340E-02, 0.13765736E-02,-0.78057044E-03,-0.12096538E-03,
     + -0.22980086E-03, 0.46705463E-04,-0.40195009E-04, 0.31347771E-03,
     +  0.23710394E-03,-0.23712542E-03, 0.23315474E-03, 0.41137281E-03,
     + -0.11574220E-03,-0.49482554E-03,-0.92103815E-03,-0.80595980E-03,
     +  0.30848427E-04, 0.31887026E-04, 0.86270717E-04,-0.13406840E-03,
     +  0.13020070E-03, 0.14021304E-03,-0.33864411E-03, 0.14702878E-03,
     + -0.20789393E-03, 0.10520906E-03, 0.90698735E-03, 0.52560697E-03,
     + -0.23053217E-03,-0.28831794E-03,-0.12918395E-04,-0.67327994E-04,
     + -0.38370406E-04,-0.65574044E-04,-0.40804498E-04, 0.12126722E-03,
     +  0.11631755E-03,-0.75121388E-04, 0.45011620E-03,-0.11345084E-03,
     + -0.36708119E-04,-0.10754010E-03,-0.24091236E-03, 0.15335279E-04,
     + -0.14505122E-03,-0.15385079E-03, 0.51398158E-04,-0.55618404E-03,
     + -0.28747256E-03, 0.16968651E-03, 0.30066798E-04,-0.57720949E-05,
     + -0.27158831E-05, 0.88660154E-05, 0.66781477E-05, 0.29897887E-04,
     +  0.28688413E-04, 0.17017359E-04, 0.76769611E-05,-0.15053235E-03,
     +  0.55600402E-04, 0.14276126E-04, 0.11973261E-04, 0.88564320E-05,
     +  0.10412546E-04, 0.47198333E-04, 0.45423996E-04, 0.77060804E-05,
     +  0.12186261E-04, 0.82979768E-05, 0.19957122E-04,-0.16617509E-04,
     +  0.12517995E-03, 0.91424818E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)    *x22*x31        
      x_sp_sm     =x_sp_sm     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x24            
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x22*x31*x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x24    *x41    
     2  +coeff( 20)        *x32        
     3  +coeff( 21)            *x42    
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)    *x23*x31        
      x_sp_sm     =x_sp_sm     
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x22*x32        
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)    *x24*x31        
     4  +coeff( 31)    *x24*x31*x41    
     5  +coeff( 32)    *x24    *x42    
     6  +coeff( 33)*x11                
     7  +coeff( 34)        *x31    *x51
     8  +coeff( 35)    *x21*x32        
      x_sp_sm     =x_sp_sm     
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)*x11*x21*x31        
     2  +coeff( 38)*x11*x21*x31*x41    
     3  +coeff( 39)    *x23*x31*x41    
     4  +coeff( 40)*x11*x21    *x42    
     5  +coeff( 41)    *x23    *x42    
     6  +coeff( 42)*x11*x24            
     7  +coeff( 43)    *x22*x31*x42    
     8  +coeff( 44)    *x22    *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 45)*x11*x23    *x41    
     1  +coeff( 46)    *x24*x32*x41    
     2  +coeff( 47)    *x21        *x51
     3  +coeff( 48)        *x31*x42    
     4  +coeff( 49)            *x43    
     5  +coeff( 50)    *x22*x31    *x51
     6  +coeff( 51)    *x22    *x41*x51
     7  +coeff( 52)        *x32*x42    
     8  +coeff( 53)        *x31*x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 54)    *x23*x32        
     1  +coeff( 55)    *x22*x32*x41    
     2  +coeff( 56)*x11*x23*x31        
     3  +coeff( 57)    *x21*x32*x42    
     4  +coeff( 58)    *x21*x31*x43    
     5  +coeff( 59)    *x24*x32        
     6  +coeff( 60)        *x32*x43    
     7  +coeff( 61)*x11*x23*x31*x41    
     8  +coeff( 62)*x11*x23    *x42    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 63)    *x24*x33        
     1  +coeff( 64)    *x24*x31*x42    
     2  +coeff( 65)    *x24    *x43    
     3  +coeff( 66)    *x23*x31*x43    
     4  +coeff( 67)*x11        *x41    
     5  +coeff( 68)    *x21*x31    *x51
     6  +coeff( 69)        *x32*x41    
     7  +coeff( 70)        *x31*x41*x51
     8  +coeff( 71)            *x41*x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 72)    *x21*x31*x42    
     1  +coeff( 73)    *x21    *x43    
     2  +coeff( 74)    *x21    *x42*x51
     3  +coeff( 75)*x12        *x41    
     4  +coeff( 76)*x11*x22    *x41    
     5  +coeff( 77)        *x33*x41    
     6  +coeff( 78)            *x42*x52
     7  +coeff( 79)*x11*x21*x32        
     8  +coeff( 80)*x11    *x33        
      x_sp_sm     =x_sp_sm     
     9  +coeff( 81)*x12        *x42    
     1  +coeff( 82)        *x33*x42    
     2  +coeff( 83)    *x24*x31    *x51
     3  +coeff( 84)*x12        *x41*x51
     4  +coeff( 85)    *x24        *x52
     5  +coeff( 86)            *x42*x53
     6  +coeff( 87)*x11*x21    *x43    
     7  +coeff( 88)*x11*x24*x31        
     8  +coeff( 89)*x11*x24    *x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 90)*x11    *x33    *x51
c
      return
      end
      function t_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1604923E+00/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38543399E-02, 0.10386721E-02,-0.12099126E-02,-0.28883358E-01,
     + -0.34110737E-02, 0.25739477E-02, 0.27118667E-02, 0.14681452E-02,
     +  0.12214561E-02, 0.59192616E-03,-0.74682769E-03, 0.32036311E-04,
     + -0.91855548E-03,-0.65080984E-03, 0.19378909E-03,-0.12324257E-02,
     +  0.20415077E-02, 0.16805300E-02,-0.48191904E-03, 0.62834343E-03,
     + -0.16232405E-03, 0.48418868E-04, 0.87669752E-04, 0.53456682E-03,
     +  0.67179056E-03,-0.13140582E-03, 0.55391050E-03,-0.59390673E-03,
     +  0.33088730E-03,-0.78999018E-03,-0.60029636E-03,-0.59417280E-03,
     +  0.62138475E-04, 0.98039825E-04,-0.59028051E-04, 0.43809876E-04,
     + -0.10572325E-03, 0.25671345E-03, 0.28698618E-03,-0.12569074E-03,
     + -0.52064791E-03, 0.50161075E-03, 0.42825486E-03, 0.85771142E-04,
     + -0.25513693E-03, 0.14909086E-03,-0.27281520E-03,-0.25017379E-03,
     + -0.19758285E-04, 0.57110679E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x22*x31        
      t_sp_sm     =t_sp_sm     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x33        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x22*x31*x41    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)*x11*x23            
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)        *x32        
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)    *x22        *x51
      t_sp_sm     =t_sp_sm     
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x23*x31        
     2  +coeff( 29)    *x22*x32        
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)*x11*x23    *x41    
     5  +coeff( 32)    *x22*x31*x43    
     6  +coeff( 33)*x11                
     7  +coeff( 34)    *x21*x32        
     8  +coeff( 35)            *x42*x51
      t_sp_sm     =t_sp_sm     
     9  +coeff( 36)            *x41*x52
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)*x11*x21*x31        
     3  +coeff( 39)        *x31*x43    
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)    *x23*x31*x41    
     6  +coeff( 42)    *x22*x31*x42    
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)        *x32    *x53
      t_sp_sm     =t_sp_sm     
     9  +coeff( 45)*x11*x23*x31        
     1  +coeff( 46)    *x23*x32*x41    
     2  +coeff( 47)    *x23    *x42*x52
     3  +coeff( 48)    *x23    *x43*x52
     4  +coeff( 49)    *x21        *x51
     5  +coeff( 50)        *x32*x41    
c
      return
      end
      function y_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  8)
      data ncoeff/  7/
      data avdat/  0.6688596E-03/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.79294143E-03,-0.29369961E-04, 0.13579262E-05, 0.90666190E-01,
     +  0.49539772E-02,-0.10769847E-02,-0.50613756E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21*x31        
c
      return
      end
      function p_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/  0.2295544E-03/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.77782694E-03, 0.49143501E-01,-0.45296857E-02,-0.22390520E-02,
     +  0.20944686E-02, 0.82695682E-03, 0.29733144E-02,-0.54925954E-03,
     +  0.10136085E-02, 0.16479676E-02,-0.19043684E-03,-0.29306018E-03,
     + -0.26312904E-03, 0.54891419E-03,-0.17523753E-02,-0.17246724E-02,
     +  0.38948722E-03, 0.15618079E-02, 0.14168711E-02, 0.89284833E-04,
     +  0.17805706E-03,-0.22720823E-03,-0.43644576E-03, 0.11493733E-03,
     + -0.19342029E-03, 0.33928122E-03,-0.10144227E-03, 0.66817345E-04,
     +  0.16487326E-03,-0.26799669E-03, 0.22571572E-03,-0.15474636E-03,
     + -0.17349143E-03,-0.78372927E-04, 0.15951828E-03, 0.31303184E-03,
     + -0.99986777E-04, 0.18874275E-03, 0.51454834E-04,-0.10536281E-03,
     +  0.43261368E-03, 0.31087128E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)            *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23*x31        
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)        *x31        
     4  +coeff( 13)*x11                
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      p_sp_sm     =p_sp_sm     
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)    *x21        *x51
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x32        
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)*x11*x23            
     8  +coeff( 26)*x11*x22    *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff( 27)*x11    *x31        
     1  +coeff( 28)    *x21*x31    *x51
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)        *x33*x41    
     4  +coeff( 31)    *x22    *x42    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)            *x42*x52
     8  +coeff( 35)*x11*x22*x31        
      p_sp_sm     =p_sp_sm     
     9  +coeff( 36)    *x23*x32        
     1  +coeff( 37)*x12        *x42    
     2  +coeff( 38)    *x21*x32*x42    
     3  +coeff( 39)*x12            *x52
     4  +coeff( 40)*x12*x22    *x41    
     5  +coeff( 41)    *x22*x33*x41    
     6  +coeff( 42)    *x22*x32*x42    
c
      return
      end
      function l_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/ -0.1920060E-02/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16516329E-02, 0.18158886E-02, 0.47300989E-02, 0.11195156E-03,
     + -0.24329007E-02, 0.22384856E-04,-0.62287936E-03,-0.34125518E-04,
     + -0.28805203E-04,-0.31465301E-04,-0.46970359E-04,-0.56644381E-04,
     +  0.14785050E-04,-0.22177954E-04, 0.10551767E-04,-0.70174456E-05,
     + -0.15222380E-04, 0.24731615E-04, 0.22676324E-04,-0.59909471E-04,
     + -0.33046839E-04,-0.26607927E-05,-0.37738773E-05,-0.87052440E-05,
     + -0.57192324E-05, 0.43645887E-05, 0.41292155E-05,-0.59305794E-05,
     + -0.13849220E-04, 0.10154054E-04, 0.19292536E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21            
      l_sp_sm     =l_sp_sm     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)        *x32        
     7  +coeff( 16)                *x52
     8  +coeff( 17)*x11*x21            
      l_sp_sm     =l_sp_sm     
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11                
     5  +coeff( 23)        *x31    *x51
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)    *x22        *x51
      l_sp_sm     =l_sp_sm     
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)    *x22*x32        
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)    *x22*x31*x43    
c
      return
      end
      function x_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2748169E+00/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.93835527E-02,-0.21352610E-01,-0.62645614E-01,-0.34025181E-02,
     +  0.39081224E-02, 0.96443930E-03, 0.40984633E-02,-0.12072674E-02,
     +  0.97283476E-03, 0.48426515E-03,-0.64146816E-03,-0.13791157E-02,
     +  0.59156626E-03,-0.10150005E-02, 0.16961611E-03, 0.22922698E-02,
     + -0.85785467E-03, 0.35958730E-02, 0.32327005E-02,-0.22001062E-02,
     +  0.38713586E-03,-0.65423251E-03, 0.13020333E-03,-0.16649772E-03,
     +  0.78848947E-03, 0.72070514E-03,-0.53323142E-03, 0.55538845E-03,
     + -0.31876969E-03,-0.16613321E-02, 0.11199249E-02, 0.84186853E-04,
     + -0.59942533E-04,-0.27387615E-03, 0.70236609E-04, 0.25684241E-03,
     + -0.28659048E-03, 0.32501295E-03, 0.91310969E-03, 0.38134758E-03,
     + -0.79372403E-03,-0.69127104E-03, 0.19130635E-03, 0.96697983E-03,
     +  0.63256413E-03,-0.53305842E-03,-0.20420603E-02,-0.81025058E-03,
     + -0.95915125E-03, 0.11744379E-03, 0.10905149E-03,-0.15719779E-03,
     + -0.10296891E-03, 0.32935926E-03, 0.20592513E-03,-0.25456335E-03,
     +  0.33228437E-03, 0.37155920E-03,-0.29236844E-03,-0.34531058E-04,
     + -0.47308634E-03,-0.31656356E-03,-0.37730220E-03, 0.80650451E-03,
     +  0.45270717E-04,-0.14482561E-04, 0.18116858E-04, 0.42034175E-04,
     +  0.31066360E-04, 0.31234460E-04, 0.22566883E-03, 0.37270984E-04,
     +  0.12885936E-04,-0.10625226E-03, 0.12448395E-03, 0.23589773E-04,
     +  0.35748031E-04,-0.22709777E-04, 0.28761657E-03,-0.48681664E-04,
     +  0.27466393E-04, 0.67762725E-04, 0.98788980E-04, 0.22960099E-04,
     +  0.30592317E-04,-0.26930927E-03,-0.32654876E-03, 0.55271707E-04,
     + -0.18441286E-03,-0.37575759E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x24*x31        
      x_sp_sex    =x_sp_sex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x23    *x41    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x24    *x41    
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)            *x42    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)    *x21    *x42    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)    *x24    *x42    
     4  +coeff( 31)    *x24*x31*x43    
     5  +coeff( 32)*x11                
     6  +coeff( 33)    *x21        *x51
     7  +coeff( 34)        *x32        
     8  +coeff( 35)        *x31    *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 36)    *x21*x32        
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)*x11*x21*x31        
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)        *x31*x43    
     5  +coeff( 41)    *x23*x31*x41    
     6  +coeff( 42)    *x23    *x42    
     7  +coeff( 43)*x11*x24            
     8  +coeff( 44)    *x22*x31*x42    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 45)    *x22    *x43    
     1  +coeff( 46)*x11*x23    *x41    
     2  +coeff( 47)    *x24*x31*x41    
     3  +coeff( 48)    *x22*x32*x42    
     4  +coeff( 49)    *x22*x31*x43    
     5  +coeff( 50)    *x24*x32*x41    
     6  +coeff( 51)        *x32*x41    
     7  +coeff( 52)    *x22*x31    *x51
     8  +coeff( 53)    *x22    *x41*x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 54)    *x21*x31*x42    
     1  +coeff( 55)    *x21    *x43    
     2  +coeff( 56)    *x23*x32        
     3  +coeff( 57)*x11*x21*x31*x41    
     4  +coeff( 58)*x11*x21    *x42    
     5  +coeff( 59)*x11*x23*x31        
     6  +coeff( 60)    *x21*x32*x42    
     7  +coeff( 61)    *x24*x32        
     8  +coeff( 62)*x11*x23*x31*x41    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 63)*x11*x23    *x42    
     1  +coeff( 64)    *x24*x32*x42    
     2  +coeff( 65)*x11        *x41    
     3  +coeff( 66)    *x21*x31    *x51
     4  +coeff( 67)*x12                
     5  +coeff( 68)        *x33        
     6  +coeff( 69)        *x31*x41*x51
     7  +coeff( 70)    *x23        *x51
     8  +coeff( 71)    *x21*x32*x41    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 72)    *x21    *x42*x51
     1  +coeff( 73)*x12        *x41    
     2  +coeff( 74)*x11*x22    *x41    
     3  +coeff( 75)        *x33*x41    
     4  +coeff( 76)            *x42*x52
     5  +coeff( 77)*x11*x21*x32        
     6  +coeff( 78)*x12*x22            
     7  +coeff( 79)    *x22*x32*x41    
     8  +coeff( 80)    *x22*x31*x41*x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 81)*x11    *x31    *x52
     1  +coeff( 82)        *x33*x42    
     2  +coeff( 83)    *x24*x31    *x51
     3  +coeff( 84)*x12        *x41*x51
     4  +coeff( 85)    *x24        *x52
     5  +coeff( 86)    *x23*x32*x41    
     6  +coeff( 87)    *x23*x31*x42    
     7  +coeff( 88)*x11*x21    *x43    
     8  +coeff( 89)    *x23    *x43    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 90)*x11*x24*x31        
c
      return
      end
      function t_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2216374E+00/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.50754487E-02, 0.18630525E-02,-0.12998981E-02,-0.29351339E-01,
     + -0.64644297E-02, 0.60313907E-02, 0.30508831E-02, 0.10336877E-02,
     + -0.13680841E-02, 0.17452129E-02, 0.13116965E-02, 0.35926144E-03,
     +  0.62216789E-03,-0.12067389E-02,-0.96219685E-03,-0.27651465E-03,
     + -0.12221916E-02, 0.24560166E-02, 0.19720190E-02,-0.82598510E-03,
     + -0.21704558E-03,-0.81740882E-05, 0.72747607E-04, 0.63778926E-03,
     +  0.88346505E-03,-0.25143032E-03, 0.70262086E-05, 0.54034329E-03,
     +  0.38653779E-05,-0.53031993E-03, 0.34253090E-03,-0.88935753E-03,
     + -0.57688344E-03,-0.10941322E-02,-0.74633415E-03, 0.11110228E-03,
     + -0.16774413E-03, 0.90414585E-04, 0.68182271E-04, 0.44772049E-04,
     + -0.79776029E-04, 0.10019255E-03, 0.67729197E-04, 0.96218144E-04,
     +  0.32874823E-03, 0.21869285E-03,-0.17170170E-03, 0.11079470E-03,
     + -0.16803823E-03,-0.61693619E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)*x11*x21            
      t_sp_sex    =t_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x23    *x41    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)*x11*x23            
     3  +coeff( 21)        *x32        
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)*x11*x22            
      t_sp_sex    =t_sp_sex    
     9  +coeff( 27)*x11    *x32        
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11            *x52
     3  +coeff( 30)    *x23*x31        
     4  +coeff( 31)    *x22*x32        
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)*x11*x23    *x41    
     7  +coeff( 34)    *x22*x31*x43    
     8  +coeff( 35)    *x23    *x42*x52
      t_sp_sex    =t_sp_sex    
     9  +coeff( 36)*x11                
     1  +coeff( 37)    *x21        *x51
     2  +coeff( 38)    *x21*x32        
     3  +coeff( 39)        *x32    *x51
     4  +coeff( 40)        *x31*x41*x51
     5  +coeff( 41)            *x42*x51
     6  +coeff( 42)    *x21        *x52
     7  +coeff( 43)            *x41*x52
     8  +coeff( 44)*x11*x21*x31        
      t_sp_sex    =t_sp_sex    
     9  +coeff( 45)        *x31*x43    
     1  +coeff( 46)    *x23        *x51
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)        *x31    *x53
     4  +coeff( 49)    *x22*x33        
     5  +coeff( 50)    *x23*x31*x41    
c
      return
      end
      function y_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/  0.7532500E-03/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12396590E-02,-0.49101665E-04, 0.12713833E-04, 0.11339740E+00,
     +  0.48828041E-02,-0.35431818E-02,-0.10657887E-02, 0.65507652E-03,
     +  0.19593232E-02, 0.18431770E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x22            
      y_sp_sex    =y_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x23    *x41    
c
      return
      end
      function p_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/ -0.3026909E-04/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.75633620E-03, 0.46953391E-01,-0.54503419E-02,-0.26563348E-02,
     +  0.21117015E-02, 0.84008626E-03, 0.36788150E-02,-0.66342443E-03,
     + -0.39345905E-03, 0.23185261E-03, 0.27978613E-04, 0.12128494E-02,
     + -0.20522345E-02, 0.18583145E-02,-0.10023457E-03, 0.20334001E-02,
     + -0.27936796E-03,-0.21623827E-02, 0.40746137E-03, 0.45350808E-03,
     +  0.16204854E-02,-0.36749706E-03, 0.20073417E-03,-0.12770681E-03,
     +  0.69367042E-03,-0.50378515E-03, 0.13944459E-03,-0.97465742E-03,
     + -0.71208848E-03,-0.21008297E-03, 0.20307896E-03, 0.27147167E-04,
     + -0.21529425E-03,-0.19130039E-03, 0.76810938E-04, 0.20681070E-03,
     +  0.49215398E-03,-0.22659967E-03, 0.45120690E-03, 0.37898324E-03,
     + -0.25521600E-03, 0.59979211E-03, 0.43636034E-03,-0.43332443E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      p_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)            *x41    
      p_sp_sex    =p_sp_sex    
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x22*x33        
     7  +coeff( 16)    *x23*x31*x41    
     8  +coeff( 17)*x11        *x41    
      p_sp_sex    =p_sp_sex    
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11*x22    *x41    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)        *x31        
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x21*x32        
      p_sp_sex    =p_sp_sex    
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11*x22*x31        
     5  +coeff( 32)                *x51
     6  +coeff( 33)        *x31*x41    
     7  +coeff( 34)            *x42    
     8  +coeff( 35)    *x21*x31    *x51
      p_sp_sex    =p_sp_sex    
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)    *x22*x31*x41    
     2  +coeff( 38)    *x21*x32*x41    
     3  +coeff( 39)    *x22    *x42    
     4  +coeff( 40)    *x23*x32        
     5  +coeff( 41)*x11*x23    *x41    
     6  +coeff( 42)    *x23*x31*x42    
     7  +coeff( 43)    *x23    *x43    
     8  +coeff( 44)    *x23*x31*x43    
c
      return
      end
      function l_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/ -0.2777880E-02/
      data xmin/
     1 -0.49950E-02,-0.52029E-01,-0.19995E-01,-0.28039E-01,-0.49946E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49970E-02, 0.52011E-01, 0.19983E-01, 0.25006E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23474235E-02,-0.15900336E-03, 0.19391376E-02, 0.74882507E-02,
     +  0.60497381E-03,-0.34614997E-02, 0.95667114E-04,-0.76192524E-03,
     + -0.11744224E-03,-0.21667345E-03,-0.12533421E-03,-0.41111649E-04,
     + -0.86353335E-04, 0.93488983E-04,-0.16673206E-03, 0.66808374E-04,
     + -0.22411904E-03,-0.19055614E-03,-0.84223531E-04, 0.37233225E-04,
     +  0.11197008E-03,-0.11245149E-04, 0.31639618E-04, 0.15284641E-04,
     + -0.10563765E-04,-0.74624899E-04,-0.74998148E-04, 0.17116981E-04,
     + -0.28894203E-04,-0.55243599E-04,-0.57381327E-04, 0.54478609E-04,
     + -0.13740240E-04, 0.83092655E-05,-0.13496938E-04, 0.13328943E-04,
     +  0.66421882E-04, 0.86899818E-04,-0.74512827E-04,-0.31830688E-04,
     + -0.10322759E-04, 0.55659792E-04,-0.37390339E-04, 0.36659367E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_sp_sex    =l_sp_sex    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)    *x22*x31*x41    
      l_sp_sex    =l_sp_sex    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21*x31        
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)*x11                
     5  +coeff( 23)        *x32        
     6  +coeff( 24)    *x21        *x51
     7  +coeff( 25)        *x31    *x51
     8  +coeff( 26)    *x21*x31*x41    
      l_sp_sex    =l_sp_sex    
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)*x11*x22            
     2  +coeff( 29)*x11*x21*x31        
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)    *x22*x32        
     5  +coeff( 32)*x11*x23            
     6  +coeff( 33)    *x21*x32        
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)    *x23        *x51
      l_sp_sex    =l_sp_sex    
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x22*x31*x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)        *x32    *x53
     6  +coeff( 42)*x11*x23    *x41    
     7  +coeff( 43)    *x22*x32*x43    
     8  +coeff( 44)*x11*x23*x33        
c
      return
      end
      function x_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.9852073E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.74351467E-02, 0.40053315E-02,-0.16517057E-02,-0.25754757E-02,
     +  0.77773601E+00, 0.51842157E-01,-0.64694665E-01,-0.90816366E-02,
     +  0.15325301E-01,-0.94868811E-02,-0.14252036E-02,-0.48361802E-02,
     + -0.51240874E-02,-0.81081623E-02,-0.39725890E-02, 0.13979734E-02,
     +  0.34188768E-02, 0.67796996E-02,-0.18273866E-01,-0.85881259E-02,
     + -0.12466125E-01, 0.65167728E-02, 0.18229784E-02, 0.70218137E-02,
     + -0.13214859E-01,-0.51309681E-02,-0.48640822E-02, 0.18174677E-02,
     +  0.26826831E-03,-0.83751604E-03,-0.37679926E-02,-0.31303500E-02,
     + -0.16435430E-02, 0.48567508E-02,-0.46687284E-02,-0.66288062E-02,
     +  0.28425446E-03,-0.69812252E-02,-0.87759143E-03, 0.48522101E-03,
     +  0.21641827E-02, 0.41728470E-03,-0.42369096E-02, 0.56280964E-03,
     +  0.57315864E-02, 0.26986911E-02,-0.30955535E-02, 0.11803952E-04,
     + -0.52028066E-02, 0.38800971E-02, 0.26825522E-02, 0.43689128E-03,
     +  0.43567995E-03,-0.34098100E-03,-0.50476211E-03,-0.33472059E-03,
     +  0.16139200E-03, 0.28408335E-02, 0.19325003E-02, 0.19704672E-02,
     +  0.71813550E-03, 0.46189320E-02, 0.41585024E-02, 0.14241335E-02,
     +  0.29650242E-04, 0.64634456E-03, 0.67572342E-03, 0.42305035E-02,
     +  0.50536883E-02, 0.48888582E-02, 0.22153258E-02,-0.89074794E-03,
     +  0.32005380E-03, 0.21651816E-02, 0.10852581E-02,-0.35194501E-02,
     + -0.35121860E-02, 0.22935444E-02,-0.75735926E-03, 0.52686148E-04,
     +  0.28538369E-02, 0.13725280E-02,-0.23475528E-03, 0.23677181E-03,
     +  0.10647619E-02, 0.83057693E-03, 0.14907090E-02, 0.35945483E-03,
     +  0.90310810E-03, 0.71798760E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      x_sp_fp     =x_sp_fp     
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)            *x42    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      x_sp_fp     =x_sp_fp     
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x24            
     4  +coeff( 22)                *x53
     5  +coeff( 23)    *x22*x31        
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)        *x31*x41*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)    *x23*x32        
     2  +coeff( 29)        *x32        
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)    *x21*x32        
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)        *x32    *x51
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)    *x22        *x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 36)    *x24        *x51
     1  +coeff( 37)    *x23*x32    *x51
     2  +coeff( 38)    *x23    *x42*x52
     3  +coeff( 39)*x11        *x41    
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)    *x23*x31        
     6  +coeff( 42)*x11*x21    *x41    
     7  +coeff( 43)    *x23        *x51
     8  +coeff( 44)*x11*x22*x31        
      x_sp_fp     =x_sp_fp     
     9  +coeff( 45)            *x42*x52
     1  +coeff( 46)    *x23    *x41*x51
     2  +coeff( 47)    *x23        *x52
     3  +coeff( 48)*x11    *x33        
     4  +coeff( 49)    *x21*x31*x41*x52
     5  +coeff( 50)    *x24    *x42    
     6  +coeff( 51)    *x23    *x42*x51
     7  +coeff( 52)*x11*x24    *x41    
     8  +coeff( 53)        *x33*x41*x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 54)    *x23*x33*x41*x51
     1  +coeff( 55)*x11    *x31        
     2  +coeff( 56)*x11            *x51
     3  +coeff( 57)*x11*x21*x31        
     4  +coeff( 58)    *x22*x31*x41    
     5  +coeff( 59)    *x22    *x42    
     6  +coeff( 60)    *x22    *x41*x51
     7  +coeff( 61)    *x21    *x43    
     8  +coeff( 62)    *x21*x31*x41*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 63)    *x21    *x42*x51
     1  +coeff( 64)    *x24*x31        
     2  +coeff( 65)        *x33*x41    
     3  +coeff( 66)*x11*x22        *x51
     4  +coeff( 67)        *x32    *x52
     5  +coeff( 68)        *x31*x41*x52
     6  +coeff( 69)    *x23*x31*x41    
     7  +coeff( 70)    *x23    *x42    
     8  +coeff( 71)    *x23*x31    *x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 72)    *x21*x32    *x52
     1  +coeff( 73)    *x24*x32        
     2  +coeff( 74)    *x24*x31*x41    
     3  +coeff( 75)    *x24*x31    *x51
     4  +coeff( 76)        *x31*x41*x53
     5  +coeff( 77)            *x42*x53
     6  +coeff( 78)    *x23*x31*x41*x51
     7  +coeff( 79)*x11        *x43*x51
     8  +coeff( 80)*x11    *x31    *x53
      x_sp_fp     =x_sp_fp     
     9  +coeff( 81)    *x21*x31*x43*x52
     1  +coeff( 82)    *x22*x32        
     2  +coeff( 83)*x11    *x31    *x51
     3  +coeff( 84)*x11*x23            
     4  +coeff( 85)    *x21*x32    *x51
     5  +coeff( 86)    *x21*x31    *x52
     6  +coeff( 87)    *x21    *x41*x52
     7  +coeff( 88)    *x21        *x53
     8  +coeff( 89)*x11*x22    *x41    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 90)            *x41*x53
c
      return
      end
      function t_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/ -0.1790684E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80665067E-03,-0.20141661E-01, 0.12558750E-03, 0.38379387E-03,
     +  0.10567504E+00, 0.50257533E-02,-0.98496545E-02,-0.77008683E-03,
     + -0.12157555E-02, 0.17041595E-02, 0.27494080E-03,-0.67610125E-03,
     + -0.41685873E-03, 0.81571453E-03, 0.20081832E-03,-0.71512995E-03,
     + -0.41112246E-03, 0.17510440E-03,-0.32036693E-04,-0.10953841E-02,
     + -0.80362603E-03,-0.10917467E-02, 0.49598847E-03,-0.94151578E-03,
     +  0.12302806E-03,-0.65460124E-04, 0.31307922E-03,-0.97514428E-04,
     + -0.68024499E-03, 0.17395050E-03,-0.19648527E-03,-0.29064342E-03,
     + -0.10557306E-02,-0.14724440E-03, 0.39173165E-03,-0.25849292E-03,
     +  0.58300130E-03,-0.25823287E-03, 0.48971031E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_sp_fp     =t_sp_fp     
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)        *x32        
     7  +coeff( 16)            *x42    
     8  +coeff( 17)            *x41*x51
      t_sp_fp     =t_sp_fp     
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22*x32    *x51
     8  +coeff( 26)        *x33*x41*x51
      t_sp_fp     =t_sp_fp     
     9  +coeff( 27)    *x22        *x53
     1  +coeff( 28)        *x31    *x51
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)    *x22        *x51
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)        *x32    *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x22    *x42    
      t_sp_fp     =t_sp_fp     
     9  +coeff( 36)    *x22        *x52
     1  +coeff( 37)            *x42*x52
     2  +coeff( 38)    *x23*x32        
     3  +coeff( 39)        *x33*x41*x52
c
      return
      end
      function y_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2472616E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56944136E-03, 0.19783553E-01,-0.36136975E-03,-0.83306842E-02,
     +  0.25294558E-02,-0.54302614E-03,-0.23446153E-02,-0.28259342E-02,
     + -0.16170947E-01,-0.44874314E-01, 0.19441938E-03,-0.27996064E-02,
     + -0.21287960E-02, 0.84196050E-02, 0.37491480E-02, 0.13208697E-02,
     +  0.15043268E-02, 0.27567835E-02, 0.73606702E-02, 0.22120556E-01,
     +  0.19125799E-02, 0.34784796E-02, 0.57147152E-03,-0.18116403E-02,
     +  0.14857551E-01, 0.29428569E-02,-0.69977599E-02,-0.28307007E-02,
     + -0.74458291E-03,-0.30018825E-02,-0.50126226E-03, 0.10367483E-02,
     +  0.27199537E-02,-0.60155871E-03, 0.35359606E-02,-0.87219762E-03,
     +  0.47293473E-02,-0.16875969E-02, 0.39370505E-04,-0.12471529E-02,
     +  0.34140691E-03,-0.16344236E-02, 0.24668558E-02, 0.28944184E-03,
     +  0.33542307E-03, 0.21952977E-02,-0.90485153E-03,-0.87978708E-03,
     + -0.44753155E-03,-0.68028504E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_fp     =y_sp_fp     
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)        *x32*x41    
     8  +coeff( 17)            *x42*x51
      y_sp_fp     =y_sp_fp     
     9  +coeff( 18)    *x22            
     1  +coeff( 19)        *x31    *x52
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x21*x31    *x51
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)    *x22    *x41    
     8  +coeff( 26)    *x22        *x51
      y_sp_fp     =y_sp_fp     
     9  +coeff( 27)            *x41*x53
     1  +coeff( 28)    *x23            
     2  +coeff( 29)        *x32*x43    
     3  +coeff( 30)    *x22    *x41*x51
     4  +coeff( 31)    *x21    *x43*x51
     5  +coeff( 32)    *x22*x33        
     6  +coeff( 33)    *x23    *x41*x51
     7  +coeff( 34)    *x22*x35        
     8  +coeff( 35)            *x43    
      y_sp_fp     =y_sp_fp     
     9  +coeff( 36)    *x21    *x42    
     1  +coeff( 37)    *x22*x31        
     2  +coeff( 38)        *x31    *x53
     3  +coeff( 39)        *x31*x44    
     4  +coeff( 40)                *x54
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)    *x22        *x52
     7  +coeff( 43)        *x31*x42    
     8  +coeff( 44)        *x31*x41*x51
      y_sp_fp     =y_sp_fp     
     9  +coeff( 45)    *x21*x32        
     1  +coeff( 46)    *x21    *x41*x51
     2  +coeff( 47)            *x43*x51
     3  +coeff( 48)    *x21*x31*x41*x51
     4  +coeff( 49)*x11*x22            
     5  +coeff( 50)    *x23*x31        
c
      return
      end
      function p_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3799727E-02/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10192617E-02, 0.23049899E-03, 0.19000314E-01, 0.26570668E-01,
     + -0.29560388E-02, 0.40447304E-03, 0.11489868E-01,-0.86469873E-03,
     + -0.11814415E-02, 0.30834733E-02,-0.51671560E-02,-0.21602303E-01,
     + -0.27922005E-02, 0.11817635E-02, 0.76754973E-02, 0.26537615E-02,
     +  0.46720202E-02,-0.14931844E-02, 0.16847776E-02, 0.14620366E-02,
     +  0.25313487E-03,-0.90781838E-03, 0.15809712E-02, 0.72726927E-03,
     +  0.53943880E-03,-0.18726724E-03, 0.16343591E-03,-0.20530383E-03,
     +  0.89980240E-04,-0.10736001E-03, 0.49983547E-03,-0.20271541E-03,
     +  0.17707584E-03,-0.83582109E-03,-0.12313814E-02,-0.60138875E-03,
     +  0.46969235E-04, 0.84480824E-03,-0.29625141E-03,-0.42035955E-03,
     +  0.53066638E-03, 0.57564850E-03,-0.32366675E-03,-0.21535319E-04,
     +  0.20531853E-03,-0.87737535E-04,-0.13257824E-03, 0.14979701E-03,
     + -0.13447221E-03,-0.39451465E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_fp     =p_sp_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)            *x41*x52
      p_sp_fp     =p_sp_fp     
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)    *x22            
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)                *x53
     8  +coeff( 26)        *x32        
      p_sp_fp     =p_sp_fp     
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)            *x41*x53
      p_sp_fp     =p_sp_fp     
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)*x11                
     2  +coeff( 38)        *x31*x42    
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)    *x22*x31    *x51
     6  +coeff( 42)    *x22    *x41*x51
     7  +coeff( 43)        *x31    *x53
     8  +coeff( 44)*x11    *x31        
      p_sp_fp     =p_sp_fp     
     9  +coeff( 45)        *x32*x41    
     1  +coeff( 46)        *x32    *x51
     2  +coeff( 47)    *x21        *x52
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x21*x31*x42    
     5  +coeff( 50)    *x21    *x43    
c
      return
      end
      function l_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2514068E-01/
      data xmin/
     1 -0.49945E-02,-0.35030E-01,-0.19928E-01,-0.27864E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49899E-02, 0.39151E-01, 0.19983E-01, 0.24342E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.77875969E-02,-0.31936467E+00,-0.31890001E-01, 0.71883583E-02,
     + -0.43328125E-01, 0.47062978E-01, 0.84584123E-02,-0.32505229E-01,
     + -0.10779140E-01, 0.62895771E-02, 0.23826374E-01,-0.14176809E-01,
     +  0.36232267E-02,-0.75405170E-02,-0.37324934E-02,-0.31201099E-02,
     + -0.30093144E-02, 0.33958755E-02, 0.56882990E-02,-0.27184316E-02,
     + -0.22271253E-01,-0.74573320E-04,-0.66190475E-03,-0.39264304E-02,
     + -0.69366447E-02, 0.17775405E-01, 0.18270573E-01,-0.23681694E-02,
     +  0.53079207E-02,-0.42485674E-02,-0.21629543E-02,-0.26483077E-03,
     + -0.54542511E-03,-0.64586150E-03, 0.18151776E-02, 0.46275961E-02,
     +  0.16549381E-02,-0.17959921E-02, 0.35917433E-02,-0.14586715E-02,
     + -0.21261999E-02,-0.82746278E-02,-0.58882907E-02, 0.25532034E-03,
     + -0.28543556E-03,-0.35034769E-03, 0.72176871E-03,-0.10589282E-02,
     + -0.38124307E-03,-0.13585492E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      l_sp_fp     =l_sp_fp     
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x23            
     4  +coeff( 13)        *x31        
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)            *x41*x51
      l_sp_fp     =l_sp_fp     
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)*x11    *x33        
     5  +coeff( 23)        *x32        
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)    *x22    *x41    
     8  +coeff( 26)    *x21*x31*x41    
      l_sp_fp     =l_sp_fp     
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)            *x42*x51
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x23*x31    *x51
     6  +coeff( 33)        *x33*x41*x51
     7  +coeff( 34)        *x31    *x51
     8  +coeff( 35)*x11    *x31        
      l_sp_fp     =l_sp_fp     
     9  +coeff( 36)    *x21*x32        
     1  +coeff( 37)    *x22        *x51
     2  +coeff( 38)    *x21*x31    *x51
     3  +coeff( 39)        *x31*x41*x51
     4  +coeff( 40)*x11*x21    *x41    
     5  +coeff( 41)*x11*x22*x31        
     6  +coeff( 42)    *x23*x31*x41    
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)*x11*x23*x31        
      l_sp_fp     =l_sp_fp     
     9  +coeff( 45)*x11*x21            
     1  +coeff( 46)*x11            *x51
     2  +coeff( 47)            *x41*x52
     3  +coeff( 48)*x11*x21*x31        
     4  +coeff( 49)*x11*x21        *x51
     5  +coeff( 50)    *x22*x31*x41    
c
      return
      end
      function x_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.1663723E-03/
      data xmin/
     1 -0.11739E+00,-0.45385E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11507E+00, 0.44097E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16379266E-02, 0.70080958E-01, 0.79911970E-01, 0.41134373E-03,
     +  0.17471267E-02, 0.41969411E-05,-0.23358156E-04,-0.20551852E-03,
     + -0.50118080E-04, 0.24175430E-04,-0.13105776E-03,-0.10753057E-03,
     +  0.47557121E-04,-0.16341259E-03, 0.64471483E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)*x11            *x51
     6  +coeff(  6)    *x23        *x52
     7  +coeff(  7)                *x51
     8  +coeff(  8)    *x23            
      x_sp_cq1x   =x_sp_cq1x   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)*x12*x21    *x42    
     6  +coeff( 15)    *x22            
c
      return
      end
      function t_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 36)
      data ncoeff/ 35/
      data avdat/ -0.3012097E-03/
      data xmin/
     1 -0.11739E+00,-0.45385E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11507E+00, 0.44097E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11302172E-03, 0.36012869E-01,-0.22066250E-04,-0.33097420E-01,
     +  0.99608524E-05, 0.72497828E-05, 0.41641961E-03, 0.15824307E-02,
     + -0.29691972E-04,-0.24371034E-03,-0.94713454E-04, 0.26074631E-03,
     + -0.76382006E-04,-0.27823859E-04, 0.63144871E-04, 0.29014638E-04,
     + -0.97200478E-04,-0.43529722E-04, 0.11632624E-03, 0.56811255E-05,
     +  0.15418060E-03,-0.25699320E-03,-0.15075007E-04, 0.10585733E-03,
     + -0.20085190E-04,-0.18821491E-04, 0.91501119E-06, 0.35463841E-05,
     +  0.77781186E-06, 0.10235804E-04,-0.16232286E-03,-0.19309773E-03,
     + -0.18830873E-04, 0.46318619E-04,-0.26511846E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)*x11            *x51
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)*x11*x22            
     2  +coeff( 11)*x11            *x52
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x23*x32        
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x23*x31*x41    
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff( 18)    *x21*x33*x41    
     1  +coeff( 19)*x12*x23    *x41    
     2  +coeff( 20)    *x22            
     3  +coeff( 21)    *x23            
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)    *x21*x33        
     6  +coeff( 24)*x11    *x32*x41    
     7  +coeff( 25)*x12*x23            
     8  +coeff( 26)*x12*x23*x32        
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff( 27)        *x32        
     1  +coeff( 28)            *x42    
     2  +coeff( 29)    *x21*x32    *x51
     3  +coeff( 30)    *x21    *x42*x51
     4  +coeff( 31)    *x21*x32*x42    
     5  +coeff( 32)*x11*x22    *x42    
     6  +coeff( 33)*x11    *x31*x42*x51
     7  +coeff( 34)    *x23*x31*x41*x51
     8  +coeff( 35)*x12*x23*x31*x41    
c
      return
      end
      function y_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1062998E-01/
      data xmin/
     1 -0.11739E+00,-0.45385E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11507E+00, 0.44097E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.78096059E-02, 0.66251449E-01, 0.62145881E-01,-0.10592953E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)        *x31    *x51
c
      return
      end
      function p_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.4036574E-02/
      data xmin/
     1 -0.11739E+00,-0.45385E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11507E+00, 0.44097E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.28291901E-02, 0.15636479E-01, 0.38245264E-01,-0.81308704E-03,
     +  0.18249585E-03,-0.35130879E-03,-0.16006817E-03, 0.29618894E-04,
     +  0.57216461E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)        *x31    *x51
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)    *x22            
      p_sp_cq1x   =p_sp_cq1x   
     9  +coeff(  9)        *x31    *x52
c
      return
      end
      function l_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/ -0.4478042E-03/
      data xmin/
     1 -0.11739E+00,-0.45385E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11507E+00, 0.44097E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.42107634E-03, 0.10516183E-03, 0.34393300E-03, 0.10474348E-05,
     + -0.14086750E-02, 0.10040980E-07,-0.11082583E-03,-0.61690691E-03,
     + -0.10966265E-02,-0.83617614E-07, 0.12884881E-04, 0.26260066E-04,
     + -0.81417757E-05, 0.13796326E-02,-0.65573760E-04, 0.33543885E-04,
     + -0.90344010E-05,-0.51944837E-03,-0.89562900E-05, 0.87574796E-06,
     + -0.49501369E-06, 0.14475892E-04, 0.45709246E-04,-0.11327957E-05,
     +  0.90198796E-06,-0.34020222E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)        *x32    *x51
     3  +coeff( 12)    *x21            
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x41*x51
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 18)*x12                
     1  +coeff( 19)*x11                
     2  +coeff( 20)    *x21        *x51
     3  +coeff( 21)                *x52
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)*x11*x21        *x51
     6  +coeff( 24)        *x32    *x52
     7  +coeff( 25)        *x31    *x52
     8  +coeff( 26)        *x32*x42*x52
c
      return
      end
      function x_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 54)
      data ncoeff/ 53/
      data avdat/ -0.5099132E+01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.56942073E-02,-0.13214388E+00, 0.48385886E-02, 0.12433766E-01,
     +  0.24268748E-02, 0.14061846E-03,-0.81837308E-04,-0.13907054E-03,
     +  0.64220010E-04,-0.39660241E-02,-0.27972669E-03, 0.23370103E-05,
     + -0.11308641E-03,-0.12380078E-02,-0.16622081E-03, 0.24484895E-03,
     +  0.25489385E-03, 0.15908161E-04,-0.21215089E-03, 0.18796338E-03,
     +  0.45895178E-03,-0.11371356E-03, 0.61420363E-03, 0.24927565E-03,
     + -0.13733386E-03,-0.23810619E-04,-0.67385508E-05, 0.31988635E-04,
     +  0.95553451E-05,-0.48547785E-03,-0.16088797E-03,-0.18423982E-03,
     +  0.32223332E-04, 0.55207747E-05,-0.32799773E-04, 0.15486439E-03,
     +  0.30938193E-03, 0.33762404E-04, 0.19808078E-04,-0.34798804E-03,
     +  0.25369276E-03,-0.24593019E-04, 0.48433958E-05,-0.34570839E-04,
     +  0.48733113E-03,-0.14479734E-04,-0.97625380E-04, 0.60964574E-03,
     +  0.29226919E-03,-0.26037893E-03,-0.12746087E-03,-0.20073199E-03,
     +  0.41380027E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21*x32        
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)                *x51
      x_sp_cden   =x_sp_cden   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)*x11            *x51
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x23        *x52
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x23*x32        
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x21        *x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)    *x21*x32*x41    
     2  +coeff( 20)    *x21*x33*x41    
     3  +coeff( 21)*x12*x21*x31*x41    
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)    *x21*x32*x42    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 27)        *x31        
     1  +coeff( 28)        *x32        
     2  +coeff( 29)                *x52
     3  +coeff( 30)    *x23            
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)        *x32    *x51
     8  +coeff( 35)    *x23        *x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 36)*x11    *x32        
     1  +coeff( 37)*x11            *x52
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)    *x21*x31*x42    
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)*x11*x24*x32        
     6  +coeff( 42)        *x31*x41    
     7  +coeff( 43)        *x31    *x51
     8  +coeff( 44)        *x31*x41*x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 45)*x11    *x31*x41    
     1  +coeff( 46)    *x21        *x53
     2  +coeff( 47)    *x23*x31*x41    
     3  +coeff( 48)    *x21*x31*x43    
     4  +coeff( 49)*x11*x22    *x42    
     5  +coeff( 50)    *x21*x32*x43    
     6  +coeff( 51)*x11*x22*x33        
     7  +coeff( 52)*x11*x22    *x43    
     8  +coeff( 53)    *x23*x33*x41    
c
      return
      end
      function t_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.1299281E+01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.93272777E-03, 0.69825202E-01,-0.13163514E-01, 0.42994688E-02,
     +  0.59107394E-03,-0.64020736E-04, 0.32704538E-04,-0.72819408E-03,
     + -0.26547511E-05, 0.82950683E-04, 0.18959078E-03, 0.20349973E-02,
     + -0.30909316E-03,-0.75575651E-03, 0.70645206E-03,-0.18236740E-02,
     + -0.10990469E-03,-0.86191008E-04,-0.11261005E-03,-0.42598371E-04,
     + -0.33818622E-04,-0.29599853E-03, 0.31739770E-03,-0.86280597E-04,
     + -0.17283307E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)    *x21*x31*x41    
      t_sp_cden   =t_sp_cden   
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x23            
     5  +coeff( 14)*x11*x23            
     6  +coeff( 15)            *x42    
     7  +coeff( 16)    *x21        *x51
     8  +coeff( 17)    *x22        *x51
      t_sp_cden   =t_sp_cden   
     9  +coeff( 18)        *x32    *x51
     1  +coeff( 19)            *x41    
     2  +coeff( 20)        *x31    *x51
     3  +coeff( 21)                *x52
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)        *x31*x41*x51
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)    *x22*x31*x41    
c
      return
      end
      function y_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/ -0.4588855E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.33370752E-02,-0.69191759E-02, 0.85821651E-01,-0.20507821E-02,
     + -0.25325489E-04, 0.47896570E-03, 0.47516685E-04, 0.55826749E-02,
     +  0.98088244E-02, 0.40029627E-02, 0.17263817E-04,-0.11643120E-02,
     + -0.84926945E-03,-0.10914000E-02, 0.10561523E-02,-0.99224120E-03,
     +  0.37871781E-03,-0.62178163E-03,-0.15046430E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sp_cden   =y_sp_cden   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x21            
     5  +coeff( 14)        *x33    *x52
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x22            
      y_sp_cden   =y_sp_cden   
     9  +coeff( 18)    *x21*x31    *x51
     1  +coeff( 19)    *x22    *x41    
c
      return
      end
      function p_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.1018158E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31051603E-02, 0.44998927E-02,-0.67737974E-01,-0.11540750E-01,
     + -0.19285496E-02,-0.81145309E-03,-0.16059630E-01,-0.15304365E-01,
     +  0.11874623E-03, 0.23640410E-03, 0.60480116E-02,-0.63494886E-02,
     +  0.16809114E-02, 0.83091421E-04, 0.79611922E-02,-0.89448094E-02,
     +  0.30722567E-02, 0.51593606E-03, 0.13961890E-02,-0.42038949E-03,
     + -0.19861736E-03, 0.12718274E-03,-0.25431023E-03,-0.13970567E-02,
     +  0.95586752E-03,-0.47053926E-03,-0.12104395E-02,-0.38148003E-03,
     +  0.22097775E-02, 0.13199232E-03, 0.19465941E-02,-0.11268014E-02,
     +  0.97911526E-03, 0.30447909E-03, 0.29966931E-03, 0.15546411E-03,
     + -0.79691055E-03,-0.20633748E-03, 0.37322659E-03, 0.14377913E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sp_cden   =p_sp_cden   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22            
     5  +coeff( 14)        *x32        
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)*x11    *x31        
      p_sp_cden   =p_sp_cden   
     9  +coeff( 18)        *x32*x41    
     1  +coeff( 19)    *x21*x31    *x51
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)                *x52
     5  +coeff( 23)        *x33        
     6  +coeff( 24)        *x31*x42    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)            *x41*x52
      p_sp_cden   =p_sp_cden   
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)*x11*x23*x31        
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)    *x23            
     4  +coeff( 31)*x11*x21*x31        
     5  +coeff( 32)*x11    *x31    *x51
     6  +coeff( 33)    *x23    *x41    
     7  +coeff( 34)    *x22*x31*x41    
     8  +coeff( 35)    *x21*x32*x41    
      p_sp_cden   =p_sp_cden   
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)    *x22*x32*x41    
     2  +coeff( 38)    *x21*x32*x42    
     3  +coeff( 39)    *x23*x31    *x51
     4  +coeff( 40)*x11*x23    *x41    
c
      return
      end
      function l_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.1533945E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14215362E-01, 0.25556234E+00,-0.21277634E-01,-0.24743851E-01,
     + -0.10999225E-02,-0.47750007E-02, 0.58291596E-04, 0.71241241E-03,
     + -0.14720244E-02, 0.66592772E-02, 0.14719604E-03,-0.48966482E-04,
     +  0.29310238E-03, 0.23670794E-03, 0.78414166E-02,-0.29238833E-02,
     +  0.53042552E-03, 0.60429936E-03,-0.20302716E-02,-0.15543170E-03,
     + -0.19689131E-03, 0.35972151E-03,-0.34999553E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)*x11                
     5  +coeff(  5)        *x32        
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)        *x32*x41    
     8  +coeff(  8)            *x41    
      l_sp_cden   =l_sp_cden   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x23*x31        
     4  +coeff( 13)                *x51
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)        *x31        
      l_sp_cden   =l_sp_cden   
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)        *x31    *x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)    *x21        *x52
c
      return
      end
      function x_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 80)
      data ncoeff/ 79/
      data avdat/ -0.5309471E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15285586E-01, 0.28345057E+00, 0.32121749E-03, 0.17726225E-03,
     +  0.13951607E+00,-0.60192190E-01, 0.21401839E-01, 0.18145628E-01,
     + -0.58611897E-02,-0.12466680E-02, 0.36749193E-02,-0.12120471E-02,
     +  0.12589335E-03, 0.28967654E-03, 0.84003661E-03,-0.65142690E-03,
     +  0.49566990E-03, 0.12168040E-02, 0.51551149E-02,-0.21248562E-02,
     +  0.53795418E-04,-0.19019977E-02,-0.31339878E-02,-0.15862423E-02,
     +  0.36004046E-03, 0.14975699E-03, 0.25985439E-03,-0.10787777E-05,
     +  0.74104035E-04, 0.76728390E-03,-0.49286924E-03,-0.11592360E-02,
     + -0.14666485E-02, 0.44356431E-04, 0.70612674E-04, 0.41022527E-03,
     + -0.37396880E-04, 0.29665173E-03, 0.94063178E-03, 0.39942828E-04,
     +  0.83479937E-03,-0.12475435E-02,-0.44681612E-04, 0.30887267E-03,
     + -0.44665311E-03,-0.53646398E-03, 0.49130170E-03,-0.12247264E-02,
     +  0.69475491E-05,-0.25037793E-03,-0.14771586E-04,-0.40088443E-03,
     + -0.65420364E-03,-0.30590862E-03,-0.25232811E-03, 0.22353391E-03,
     +  0.91347138E-04, 0.12377929E-03,-0.12614432E-03,-0.64785639E-03,
     + -0.12156190E-03,-0.50580635E-03, 0.52447844E-03, 0.49485912E-03,
     + -0.58561178E-04,-0.89273562E-04,-0.85056847E-04,-0.81403137E-04,
     + -0.85910702E-04, 0.45144118E-03, 0.34738454E-03, 0.86359546E-03,
     +  0.16940688E-03, 0.14711162E-02,-0.61116071E-03,-0.88248064E-03,
     + -0.20866771E-03,-0.29575936E-02,-0.29397555E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21        *x51
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff(  9)                *x52
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x23*x31        
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x41    
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x24            
     3  +coeff( 21)*x11            *x52
     4  +coeff( 22)            *x42    
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)                *x53
     8  +coeff( 26)    *x22*x32        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 27)    *x21*x32*x41    
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)*x12    *x31        
     3  +coeff( 30)    *x23*x32        
     4  +coeff( 31)    *x21*x33*x41    
     5  +coeff( 32)*x11*x22*x31*x41    
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)        *x32*x41    
     8  +coeff( 35)        *x32    *x51
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 36)        *x31*x41*x51
     1  +coeff( 37)    *x22*x31    *x51
     2  +coeff( 38)    *x21*x32    *x51
     3  +coeff( 39)*x11*x22    *x41    
     4  +coeff( 40)    *x23*x31    *x51
     5  +coeff( 41)    *x21*x32*x42    
     6  +coeff( 42)*x11    *x31*x43    
     7  +coeff( 43)        *x31    *x51
     8  +coeff( 44)*x11    *x31        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)            *x42*x51
     2  +coeff( 47)*x11*x21        *x51
     3  +coeff( 48)    *x22*x31*x41    
     4  +coeff( 49)*x11    *x31    *x51
     5  +coeff( 50)    *x22        *x52
     6  +coeff( 51)        *x33*x41    
     7  +coeff( 52)    *x24        *x51
     8  +coeff( 53)    *x23    *x42    
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 54)    *x23        *x52
     1  +coeff( 55)*x11*x24            
     2  +coeff( 56)    *x22*x32*x41    
     3  +coeff( 57)    *x22*x32    *x51
     4  +coeff( 58)    *x22    *x41    
     5  +coeff( 59)    *x21*x31    *x51
     6  +coeff( 60)    *x23    *x41    
     7  +coeff( 61)    *x23        *x51
     8  +coeff( 62)*x11    *x32        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 63)    *x22    *x42    
     1  +coeff( 64)*x11*x23            
     2  +coeff( 65)    *x21        *x53
     3  +coeff( 66)        *x31*x43    
     4  +coeff( 67)            *x42*x52
     5  +coeff( 68)    *x22        *x53
     6  +coeff( 69)    *x21*x32    *x52
     7  +coeff( 70)    *x21*x33*x42    
     8  +coeff( 71)    *x21*x32*x43    
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 72)*x11*x22*x32*x41    
     1  +coeff( 73)*x12        *x42*x51
     2  +coeff( 74)    *x23*x33*x41    
     3  +coeff( 75)*x11*x24*x32        
     4  +coeff( 76)    *x21*x33*x43    
     5  +coeff( 77)*x12    *x33*x41    
     6  +coeff( 78)*x11*x22*x33*x41    
     7  +coeff( 79)    *x24*x33*x41    
c
      return
      end
      function t_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 34)
      data ncoeff/ 33/
      data avdat/ -0.1464176E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21504872E-02,-0.83882228E-01, 0.16781517E-03, 0.22745371E-01,
     + -0.77687297E-02, 0.37730468E-03,-0.14523474E-02, 0.79887444E-02,
     + -0.35458503E-03,-0.48475391E-04, 0.89689094E-03,-0.41902375E-04,
     +  0.33580905E-03,-0.17618206E-03,-0.52222946E-04,-0.27558827E-03,
     +  0.26435687E-02,-0.20861458E-02,-0.44194839E-03, 0.30141370E-03,
     +  0.16480000E-03,-0.50006951E-04, 0.12485670E-04, 0.26862920E-03,
     +  0.82581845E-03,-0.14912919E-02, 0.60908473E-03, 0.84331128E-04,
     + -0.97789729E-04, 0.29494916E-03,-0.20302049E-03,-0.19213662E-03,
     +  0.49538154E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x23*x31        
     4  +coeff( 13)*x12    *x32        
     5  +coeff( 14)        *x31        
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)    *x21    *x41    
     8  +coeff( 17)*x11*x21            
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)        *x32    *x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)        *x32    *x52
     4  +coeff( 22)        *x31*x41*x52
     5  +coeff( 23)            *x41*x53
     6  +coeff( 24)    *x22        *x53
     7  +coeff( 25)        *x32        
     8  +coeff( 26)        *x31*x41    
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 27)            *x42    
     1  +coeff( 28)        *x31    *x51
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)    *x23        *x51
     5  +coeff( 32)*x11*x23            
     6  +coeff( 33)*x11*x22*x31*x41    
c
      return
      end
      function y_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.1532255E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14408195E-02,-0.17508604E+00, 0.13452946E+00,-0.81634531E-02,
     +  0.45433268E-02, 0.65182737E-03,-0.29651416E-03, 0.17617790E-02,
     +  0.24368899E-01, 0.34784023E-01, 0.26663646E-03, 0.59483532E-03,
     + -0.27202419E-02,-0.41595764E-01,-0.37435153E-02,-0.19106567E-01,
     + -0.52086120E-02, 0.28099872E-02,-0.68496435E-03,-0.47714412E-03,
     + -0.26473599E-02, 0.59388479E-03,-0.41041044E-02, 0.19154307E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)*x11        *x41*x51
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff( 18)*x12                
     1  +coeff( 19)        *x33    *x52
     2  +coeff( 20)        *x32*x41    
     3  +coeff( 21)        *x31    *x52
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)        *x31*x42    
     6  +coeff( 24)    *x23*x31        
c
      return
      end
      function p_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/  0.3825038E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.68731565E-03,-0.37111606E-01, 0.19509820E-01,-0.13631538E-02,
     + -0.24719682E-03,-0.37747363E-02,-0.74216747E-02, 0.23763796E-03,
     +  0.55662866E-04, 0.28584457E-02, 0.80085816E-02, 0.34690430E-03,
     +  0.52295922E-03,-0.47847964E-02,-0.12320579E-02, 0.14730814E-02,
     +  0.11072259E-02,-0.12344064E-02,-0.11676901E-02, 0.13362992E-02,
     +  0.21908840E-03, 0.22396218E-03, 0.10083312E-03,-0.51790872E-03,
     + -0.41230658E-04,-0.67458772E-04, 0.11075645E-02, 0.75731616E-04,
     + -0.62959596E-04, 0.30156857E-03,-0.96546537E-04,-0.79763331E-03,
     + -0.49471145E-03, 0.58267557E-04,-0.81948255E-03,-0.31382737E-04,
     + -0.53414889E-03,-0.68267889E-03, 0.12279153E-03, 0.86045556E-03,
     +  0.57893936E-04, 0.26941934E-03,-0.56967483E-03, 0.54467452E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)        *x32*x41    
     7  +coeff( 16)    *x21            
     8  +coeff( 17)    *x21*x31    *x51
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)*x11    *x31        
     3  +coeff( 21)        *x33        
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)    *x21*x31    *x52
     7  +coeff( 25)*x11*x23*x31*x42    
     8  +coeff( 26)    *x21*x32        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 27)        *x31*x42    
     1  +coeff( 28)    *x21        *x52
     2  +coeff( 29)    *x22*x32        
     3  +coeff( 30)    *x22*x31*x41    
     4  +coeff( 31)        *x31    *x53
     5  +coeff( 32)    *x23*x31    *x51
     6  +coeff( 33)*x11*x23*x31        
     7  +coeff( 34)    *x23            
     8  +coeff( 35)            *x43    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)*x11    *x31    *x51
     2  +coeff( 38)    *x23*x31        
     3  +coeff( 39)    *x21*x33        
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x21*x31*x41*x51
     6  +coeff( 42)    *x21    *x41*x52
     7  +coeff( 43)    *x22*x31*x42    
     8  +coeff( 44)    *x23    *x41*x51
c
      return
      end
      function l_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 23)
      data ncoeff/ 22/
      data avdat/ -0.1012984E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11349961E-01,-0.36261192E+00, 0.44899489E-03,-0.32205462E-01,
     + -0.34129266E-01, 0.50263461E-01, 0.35769292E-02,-0.28896206E-02,
     + -0.23837311E-02, 0.10643180E-02,-0.12500651E-03, 0.10192256E-02,
     +  0.28603172E-02,-0.70621609E-03,-0.40577829E-03, 0.81829857E-02,
     + -0.95106900E-03,-0.92942426E-02, 0.11293673E-02,-0.19188165E-02,
     + -0.20445623E-03, 0.22575455E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x32        
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)        *x31        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21*x33        
     6  +coeff( 15)        *x31    *x51
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)*x11        *x41    
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)        *x32    *x51
     2  +coeff( 20)            *x42    
     3  +coeff( 21)                *x52
     4  +coeff( 22)    *x23*x33*x41    
c
      return
      end
      function x_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 61)
      data ncoeff/ 60/
      data avdat/ -0.8403291E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13680696E-01, 0.18053657E+00, 0.13913662E-03, 0.44452224E-03,
     +  0.18936408E+00, 0.19579949E-01,-0.95107052E-02, 0.11583416E-02,
     + -0.54186102E-01,-0.19863865E-03,-0.27907242E-02, 0.49906026E-04,
     +  0.11089029E-01, 0.74125570E-03,-0.81784819E-03, 0.19034108E-02,
     +  0.18339657E-02,-0.23663824E-02,-0.30778401E-03,-0.24490198E-02,
     + -0.81219873E-03, 0.69927773E-03, 0.88465589E-04,-0.13047850E-02,
     +  0.16767246E-03, 0.26906014E-02, 0.73800044E-03, 0.16137430E-02,
     +  0.63831219E-03,-0.20375306E-04, 0.10423200E-02,-0.14610460E-03,
     + -0.95131592E-03,-0.86373510E-03,-0.13545927E-02,-0.12519254E-02,
     +  0.17908579E-02,-0.95550023E-03,-0.62773674E-03,-0.13737427E-02,
     +  0.11151702E-03, 0.20895742E-02,-0.38579677E-03, 0.88372844E-03,
     + -0.24211954E-02, 0.13006750E-02, 0.17131070E-03,-0.20068903E-03,
     +  0.11969599E-03,-0.30752615E-03,-0.38790662E-03,-0.34608226E-02,
     + -0.26672665E-03,-0.12250707E-03, 0.22813614E-03,-0.32747182E-03,
     +  0.27798864E-03, 0.11484745E-03,-0.36805408E-03, 0.10886419E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11*x21            
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff(  9)*x11                
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x23*x31        
     4  +coeff( 13)    *x22            
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22        *x51
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x24            
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)                *x53
     5  +coeff( 23)    *x21    *x41    
     6  +coeff( 24)            *x42    
     7  +coeff( 25)        *x31    *x51
     8  +coeff( 26)*x11            *x51
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 27)    *x21*x31*x42    
     1  +coeff( 28)    *x21*x32    *x51
     2  +coeff( 29)        *x32    *x52
     3  +coeff( 30)    *x21*x33    *x51
     4  +coeff( 31)        *x31*x41*x51
     5  +coeff( 32)        *x31    *x52
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)    *x22        *x52
     8  +coeff( 35)    *x21*x31*x41*x51
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 36)    *x24        *x51
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x23        *x52
     3  +coeff( 39)*x11*x24            
     4  +coeff( 40)    *x21*x32*x42    
     5  +coeff( 41)    *x22    *x41    
     6  +coeff( 42)    *x21*x31*x41    
     7  +coeff( 43)    *x21*x31    *x51
     8  +coeff( 44)    *x22*x32        
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 45)    *x22*x31*x41    
     1  +coeff( 46)    *x22    *x42    
     2  +coeff( 47)    *x21*x31    *x52
     3  +coeff( 48)    *x21        *x53
     4  +coeff( 49)    *x22*x32    *x51
     5  +coeff( 50)    *x22        *x53
     6  +coeff( 51)    *x21*x32    *x52
     7  +coeff( 52)*x11*x22    *x42    
     8  +coeff( 53)        *x32    *x53
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 54)            *x41*x51
     1  +coeff( 55)    *x21    *x41*x51
     2  +coeff( 56)            *x42*x51
     3  +coeff( 57)    *x23    *x41    
     4  +coeff( 58)    *x21*x33        
     5  +coeff( 59)        *x31*x41*x52
     6  +coeff( 60)        *x31    *x53
c
      return
      end
      function t_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.3508774E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12643007E-02,-0.40095020E-01,-0.13797733E-03, 0.26000498E-03,
     +  0.58371820E-01,-0.52379398E-02,-0.49052495E-02, 0.13583234E-02,
     + -0.35413161E-02,-0.15844788E-02, 0.26870833E-03, 0.31245390E-02,
     + -0.78153447E-03, 0.43216973E-03, 0.80408424E-03, 0.13249143E-03,
     +  0.25888401E-03,-0.72014530E-03,-0.37733684E-03,-0.53188601E-03,
     + -0.15467541E-03, 0.58443501E-03,-0.12048005E-03, 0.59594115E-03,
     + -0.71166427E-03,-0.19175861E-03, 0.22798558E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21        *x51
      t_sp_cq3m   =t_sp_cq3m   
     9  +coeff(  9)*x11                
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)        *x32    *x51
     5  +coeff( 14)                *x53
     6  +coeff( 15)        *x32        
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)        *x32    *x52
      t_sp_cq3m   =t_sp_cq3m   
     9  +coeff( 18)*x11*x23            
     1  +coeff( 19)*x12*x21        *x51
     2  +coeff( 20)    *x21*x32    *x52
     3  +coeff( 21)    *x23            
     4  +coeff( 22)        *x31*x41*x51
     5  +coeff( 23)        *x31    *x52
     6  +coeff( 24)    *x21*x32    *x51
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)    *x22        *x52
      t_sp_cq3m   =t_sp_cq3m   
     9  +coeff( 27)    *x21*x31    *x52
c
      return
      end
      function y_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.1921629E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22251876E-02,-0.20998147E+00, 0.14865577E+00,-0.90893256E-02,
     +  0.61956150E-02, 0.86880807E-03,-0.49452344E-03, 0.20843165E-02,
     +  0.25076257E-01, 0.43330681E-01,-0.16191007E-03, 0.52283902E-03,
     + -0.58434661E-02,-0.48530683E-01,-0.38905889E-02,-0.23309220E-01,
     +  0.33203124E-02, 0.20561542E-02,-0.70891487E-02,-0.26280335E-02,
     + -0.92377956E-03, 0.59696549E-03,-0.47985353E-02,-0.22503077E-02,
     +  0.30890554E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cq3m   =y_sp_cq3m   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)*x12                
      y_sp_cq3m   =y_sp_cq3m   
     9  +coeff( 18)        *x33        
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)    *x22*x31    *x51
     7  +coeff( 25)*x11*x22    *x41    
c
      return
      end
      function p_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1131117E-03/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23200795E-03,-0.50459214E-04, 0.67351838E-02,-0.12281295E-01,
     +  0.80699526E-03,-0.25739250E-03,-0.15472700E-02,-0.83581122E-04,
     +  0.29626908E-02, 0.68472458E-04,-0.47588707E-02, 0.39859246E-04,
     +  0.14795795E-02, 0.16132859E-03, 0.79619419E-03, 0.23369979E-03,
     + -0.14902165E-03, 0.71111106E-03,-0.10761673E-02,-0.21207286E-03,
     + -0.74893127E-04,-0.67539193E-04,-0.81391173E-03, 0.74905163E-03,
     + -0.35302958E-03,-0.96742267E-04, 0.49706484E-03, 0.97077958E-04,
     + -0.31190881E-03,-0.20619525E-03, 0.14573590E-03,-0.14157845E-03,
     +  0.83793631E-04, 0.26322840E-03,-0.59157443E-04, 0.25875773E-03,
     +  0.11914200E-03, 0.15864529E-03, 0.31079366E-03,-0.37561473E-03,
     +  0.41725940E-04,-0.47850993E-03,-0.12250723E-03,-0.38270846E-04,
     + -0.19470474E-03, 0.47216323E-03, 0.20027603E-03, 0.84069863E-04,
     + -0.14426054E-04, 0.25451765E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 18)        *x31    *x52
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)                *x52
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)            *x41*x52
     7  +coeff( 25)    *x21*x31    *x52
     8  +coeff( 26)            *x42    
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 27)            *x41*x51
     1  +coeff( 28)        *x33        
     2  +coeff( 29)    *x23*x31        
     3  +coeff( 30)*x12        *x41    
     4  +coeff( 31)*x12            *x51
     5  +coeff( 32)        *x31    *x53
     6  +coeff( 33)    *x21        *x52
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)        *x33    *x51
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)    *x22*x31*x41*x51
     3  +coeff( 39)*x11*x23*x31    *x52
     4  +coeff( 40)    *x22*x33*x42*x51
     5  +coeff( 41)    *x23            
     6  +coeff( 42)    *x22    *x41    
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)    *x21*x32    *x51
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 45)    *x22    *x41*x51
     1  +coeff( 46)*x11*x21*x31    *x51
     2  +coeff( 47)    *x23*x31    *x51
     3  +coeff( 48)    *x21*x31    *x53
     4  +coeff( 49)*x11*x23*x31        
     5  +coeff( 50)*x11*x23    *x41    
c
      return
      end
      function l_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.1152398E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10029471E-01,-0.36317816E+00, 0.35006937E-03,-0.31909894E-01,
     + -0.38123038E-01, 0.50472748E-01,-0.25360398E-02, 0.70619546E-02,
     + -0.28486783E-02,-0.81850518E-03, 0.11892524E-02,-0.11393810E-02,
     +  0.88607827E-02, 0.68167731E-03, 0.20975081E-03, 0.53290324E-03,
     + -0.95791125E-03,-0.22467209E-02,-0.97622592E-02,-0.27409883E-03,
     + -0.45897582E-03,-0.29078571E-03, 0.27774859E-02, 0.99188229E-03,
     +  0.23237099E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)*x11                
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21        *x51
      l_sp_cq3m   =l_sp_cq3m   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x43*x51
     2  +coeff( 11)        *x31        
     3  +coeff( 12)                *x52
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)    *x21*x32*x41    
      l_sp_cq3m   =l_sp_cq3m   
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)    *x21    *x41    
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)    *x23*x33*x41    
c
      return
      end
      function x_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 77)
      data ncoeff/ 76/
      data avdat/ -0.2075658E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.52556000E-02, 0.15762660E+00, 0.73039555E-03, 0.32252622E+00,
     +  0.25687737E-01,-0.21662755E-01,-0.71911030E-01, 0.72156233E-02,
     +  0.13097656E-02,-0.42166589E-02,-0.55141128E-02,-0.44834206E-03,
     +  0.13921899E-02,-0.64290725E-02, 0.24793854E-02, 0.46925070E-02,
     +  0.25067225E-02, 0.30021528E-02, 0.68634101E-02,-0.41602864E-02,
     + -0.18080722E-02, 0.53101720E-03, 0.25004793E-02, 0.19838081E-02,
     + -0.27719259E-02,-0.17951146E-02, 0.28471034E-02,-0.25711916E-02,
     + -0.20409597E-02,-0.16121153E-03,-0.38257471E-03,-0.89143537E-03,
     + -0.12581992E-02, 0.97128452E-03,-0.11797666E-02,-0.11116089E-02,
     + -0.71119081E-04,-0.70395315E-03,-0.66667929E-03,-0.35074231E-03,
     +  0.69297108E-04,-0.37146383E-03, 0.32075483E-03, 0.25296595E-03,
     + -0.40922736E-04, 0.15574612E-03,-0.10880398E-02,-0.46192645E-03,
     +  0.56868827E-03,-0.62284217E-03, 0.31965357E-03,-0.65728644E-04,
     + -0.11733171E-02, 0.20138335E-02, 0.38293013E-03,-0.40000366E-03,
     +  0.37452483E-03, 0.10615574E-02,-0.19677393E-02,-0.43801835E-03,
     + -0.97411551E-03, 0.70210203E-03, 0.96350070E-03,-0.39055016E-04,
     +  0.50648628E-03, 0.23358665E-04,-0.40887104E-03,-0.31236582E-03,
     +  0.26213515E-03, 0.56872651E-03,-0.41978646E-03, 0.13556422E-03,
     + -0.11256621E-02, 0.45450742E-03, 0.65308064E-04,-0.21832764E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x24            
     6  +coeff( 15)                *x53
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x22        *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 27)    *x21*x32    *x51
     1  +coeff( 28)    *x24        *x51
     2  +coeff( 29)    *x23        *x52
     3  +coeff( 30)        *x31        
     4  +coeff( 31)            *x42    
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)    *x21*x31*x41*x51
     7  +coeff( 34)        *x32    *x52
     8  +coeff( 35)*x11*x24            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 36)*x11        *x42*x51
     1  +coeff( 37)        *x33    *x52
     2  +coeff( 38)    *x23*x32*x42    
     3  +coeff( 39)    *x21    *x41    
     4  +coeff( 40)            *x41*x51
     5  +coeff( 41)        *x33        
     6  +coeff( 42)        *x31    *x52
     7  +coeff( 43)    *x21*x33        
     8  +coeff( 44)    *x21    *x43    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 45)*x11*x22*x31        
     1  +coeff( 46)    *x22*x31*x41*x51
     2  +coeff( 47)    *x21*x32    *x52
     3  +coeff( 48)*x12    *x31    *x51
     4  +coeff( 49)*x12        *x41*x51
     5  +coeff( 50)        *x32    *x53
     6  +coeff( 51)    *x21    *x41*x51
     7  +coeff( 52)        *x32*x41    
     8  +coeff( 53)*x11        *x42    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 54)*x11*x23            
     1  +coeff( 55)    *x21*x31    *x52
     2  +coeff( 56)            *x42*x52
     3  +coeff( 57)        *x31    *x53
     4  +coeff( 58)    *x23*x32        
     5  +coeff( 59)    *x23*x31*x41    
     6  +coeff( 60)    *x22        *x53
     7  +coeff( 61)    *x21*x33*x41    
     8  +coeff( 62)*x12*x21        *x51
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 63)    *x24*x32        
     1  +coeff( 64)*x12    *x31*x41    
     2  +coeff( 65)    *x23*x31*x42    
     3  +coeff( 66)*x12*x22*x31        
     4  +coeff( 67)*x12*x22    *x41    
     5  +coeff( 68)    *x22*x31        
     6  +coeff( 69)*x11        *x41    
     7  +coeff( 70)    *x22    *x41    
     8  +coeff( 71)            *x42*x51
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 72)            *x41*x52
     1  +coeff( 73)    *x22*x31*x41    
     2  +coeff( 74)    *x22    *x42    
     3  +coeff( 75)        *x33    *x51
     4  +coeff( 76)            *x41*x53
c
      return
      end
      function t_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.1790809E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.57672290E-03,-0.43640202E-02, 0.36220989E-03, 0.10569038E+00,
     +  0.50981683E-02,-0.98832948E-02,-0.15841939E-01,-0.36454524E-02,
     + -0.17908234E-02,-0.88765420E-03, 0.80183410E-03, 0.39817565E-02,
     + -0.93408761E-03,-0.15024613E-02,-0.13039074E-02,-0.20708365E-03,
     +  0.31166727E-03, 0.30463107E-03, 0.11995251E-02, 0.52636722E-03,
     + -0.18921959E-03,-0.64106198E-03,-0.48465448E-03, 0.53979020E-03,
     +  0.37428871E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21        *x52
     2  +coeff( 11)        *x32        
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)*x11*x23            
     7  +coeff( 16)        *x31        
     8  +coeff( 17)    *x21*x31        
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff( 18)        *x31    *x51
     1  +coeff( 19)        *x31*x41*x51
     2  +coeff( 20)                *x53
     3  +coeff( 21)        *x31    *x52
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)        *x32    *x52
     7  +coeff( 25)*x11    *x32    *x51
c
      return
      end
      function y_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/  0.1389974E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80125947E-02,-0.15640648E+00, 0.95509313E-01,-0.39564483E-02,
     +  0.31910669E-02, 0.34524876E-03, 0.98755502E-03, 0.10934507E-01,
     +  0.33935204E-01,-0.11587657E-03,-0.71364194E-02,-0.32001011E-01,
     + -0.71570068E-03,-0.18015521E-01,-0.44225841E-02, 0.14897709E-02,
     + -0.47365800E-02, 0.11387543E-02,-0.28007501E-02, 0.12552880E-02,
     +  0.22095246E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)*x11        *x41*x51
     7  +coeff( 16)        *x33        
     8  +coeff( 17)        *x32*x41    
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff( 18)    *x22            
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)*x12*x21    *x41    
c
      return
      end
      function p_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3799357E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24806496E-02,-0.73934870E-03, 0.48939902E-01,-0.39513644E-01,
     +  0.24665503E-02, 0.86887208E-04,-0.54094085E-03,-0.16216209E-03,
     +  0.13382467E-01, 0.35941481E-03,-0.37648735E-03,-0.25057801E-03,
     + -0.11408830E-01,-0.65590083E-02, 0.21893517E-02, 0.34227627E-02,
     +  0.40344372E-02, 0.27352350E-02,-0.59091335E-03, 0.80752093E-03,
     + -0.11002240E-02,-0.74307295E-03,-0.60992455E-03, 0.45378928E-03,
     + -0.27434083E-03, 0.13798194E-02,-0.17736766E-02,-0.19260746E-03,
     +  0.10473082E-02, 0.58518734E-03,-0.10938349E-03,-0.35193824E-03,
     + -0.24035831E-04, 0.77926750E-04,-0.16598213E-03,-0.18539542E-03,
     + -0.40143833E-03, 0.28911344E-03,-0.11602186E-02, 0.69520996E-04,
     +  0.26770312E-03,-0.29702834E-03,-0.56874403E-03,-0.13462363E-02,
     +  0.12124167E-02,-0.12497714E-03,-0.17795947E-03,-0.85298503E-04,
     + -0.10482050E-03,-0.11743526E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x22    *x41    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 18)    *x21*x31    *x51
     1  +coeff( 19)                *x52
     2  +coeff( 20)        *x33        
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)        *x31    *x53
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)            *x42*x53
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)        *x31*x42    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)        *x33    *x51
     2  +coeff( 29)    *x22*x31*x42    
     3  +coeff( 30)*x11*x23*x31        
     4  +coeff( 31)    *x22*x32*x42    
     5  +coeff( 32)        *x32        
     6  +coeff( 33)    *x21        *x52
     7  +coeff( 34)                *x53
     8  +coeff( 35)    *x22    *x42    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)            *x41*x53
     2  +coeff( 38)    *x23*x31    *x51
     3  +coeff( 39)*x11        *x41    
     4  +coeff( 40)    *x21    *x42    
     5  +coeff( 41)        *x32    *x51
     6  +coeff( 42)        *x31*x41*x51
     7  +coeff( 43)*x11*x21*x31        
     8  +coeff( 44)*x11    *x31    *x51
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 45)*x11        *x41*x51
     1  +coeff( 46)    *x21*x32*x41    
     2  +coeff( 47)        *x32*x41*x51
     3  +coeff( 48)        *x32    *x52
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)*x11*x21*x32        
c
      return
      end
      function l_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/ -0.1721881E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.42729550E-02,-0.36326334E+00,-0.31879708E-01,-0.38643457E-01,
     +  0.50481968E-01,-0.75634904E-02,-0.45089764E-02, 0.71726145E-05,
     +  0.14696424E-02,-0.28804876E-02, 0.19291020E-02, 0.89283176E-02,
     +  0.87597137E-02, 0.13444221E-02,-0.13071486E-02,-0.90846408E-03,
     + -0.34051071E-03,-0.59299904E-03,-0.79640048E-02, 0.15716908E-02,
     +  0.13509012E-02,-0.24431527E-02,-0.39741481E-03, 0.29559068E-02,
     +  0.54186076E-03, 0.26059414E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)*x11                
     6  +coeff(  6)                *x52
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x33        
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)    *x21*x31        
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff( 18)        *x31    *x51
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)                *x53
     4  +coeff( 22)            *x42    
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)    *x23*x32*x42    
c
      return
      end
      function x_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 71)
      data ncoeff/ 70/
      data avdat/ -0.9852073E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78008897E-02, 0.13899027E+00, 0.21112531E-02, 0.77743626E+00,
     +  0.46513982E-01,-0.64383641E-01,-0.14100747E+00,-0.39255251E-02,
     + -0.12007952E-01,-0.32299224E-02, 0.34223285E-02, 0.36100456E-02,
     + -0.14686529E-01, 0.71755992E-02, 0.44853599E-02, 0.22415966E-02,
     +  0.17596750E-01, 0.64763832E-02, 0.16051156E-02,-0.78447051E-02,
     +  0.86752353E-02,-0.88205223E-03,-0.13514300E-01,-0.61351876E-02,
     + -0.43247370E-02,-0.52707428E-02,-0.70776530E-02, 0.43321932E-02,
     + -0.57126023E-02,-0.29524483E-02,-0.97014697E-03,-0.16605679E-02,
     + -0.14664575E-02,-0.13217140E-02,-0.91180559E-02, 0.15469552E-03,
     + -0.40792669E-02, 0.33473147E-02, 0.43867780E-02,-0.18194327E-01,
     +  0.25344141E-01,-0.14377772E-02, 0.32973683E-02, 0.53172740E-02,
     +  0.48043396E-03, 0.88125430E-02,-0.27309258E-02,-0.18779301E-02,
     + -0.50214697E-02,-0.16506338E-02,-0.49311803E-02,-0.29018938E-02,
     + -0.28866480E-03,-0.14598754E-02, 0.41859905E-04,-0.54500269E-03,
     +  0.10551164E-02,-0.24124880E-02, 0.76723669E-03, 0.41030715E-02,
     +  0.14826715E-02,-0.19591654E-03, 0.20035554E-03, 0.10722420E-02,
     +  0.32635464E-02, 0.56565722E-03,-0.43130267E-03,-0.80515660E-03,
     +  0.18718097E-02, 0.19080127E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      x_sp_cfp    =x_sp_cfp    
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21        *x52
     2  +coeff( 11)*x11    *x32        
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x24            
     5  +coeff( 14)                *x53
     6  +coeff( 15)        *x32        
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)*x11*x21            
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)        *x31        
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)    *x22        *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 27)    *x24        *x51
     1  +coeff( 28)        *x32    *x52
     2  +coeff( 29)    *x23        *x52
     3  +coeff( 30)*x11*x24            
     4  +coeff( 31)*x11    *x32    *x51
     5  +coeff( 32)    *x21    *x41    
     6  +coeff( 33)            *x41*x51
     7  +coeff( 34)        *x31    *x52
     8  +coeff( 35)    *x22*x31*x41    
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 36)*x11    *x31    *x51
     1  +coeff( 37)    *x21*x31*x41*x51
     2  +coeff( 38)    *x23            
     3  +coeff( 39)*x11            *x51
     4  +coeff( 40)    *x21*x32        
     5  +coeff( 41)    *x21*x31*x41    
     6  +coeff( 42)            *x42*x51
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)    *x22    *x42    
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 45)    *x21*x33        
     1  +coeff( 46)    *x21*x32    *x51
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)            *x42*x52
     4  +coeff( 49)    *x23*x31*x41    
     5  +coeff( 50)    *x21*x33*x41    
     6  +coeff( 51)    *x21*x32    *x52
     7  +coeff( 52)        *x32    *x53
     8  +coeff( 53)    *x23*x31    *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 54)*x11*x21        *x53
     1  +coeff( 55)        *x33    *x53
     2  +coeff( 56)    *x22*x31        
     3  +coeff( 57)    *x22    *x41    
     4  +coeff( 58)    *x21*x31    *x51
     5  +coeff( 59)    *x21    *x41*x51
     6  +coeff( 60)*x11*x23            
     7  +coeff( 61)    *x21*x31    *x52
     8  +coeff( 62)*x11*x22*x31        
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 63)        *x33    *x51
     1  +coeff( 64)        *x31    *x53
     2  +coeff( 65)    *x23*x32        
     3  +coeff( 66)*x11    *x31*x42    
     4  +coeff( 67)    *x22*x32    *x51
     5  +coeff( 68)*x11*x23*x31        
     6  +coeff( 69)    *x21*x31*x41*x52
     7  +coeff( 70)    *x24*x32        
c
      return
      end
      function t_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.1790684E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.57657878E-03,-0.43644155E-02, 0.36205805E-03, 0.10569151E+00,
     +  0.50975014E-02,-0.98830983E-02,-0.15841449E-01,-0.36459249E-02,
     + -0.17911042E-02,-0.88922156E-03, 0.80181326E-03, 0.39818492E-02,
     + -0.93435205E-03,-0.15082682E-02,-0.13034560E-02,-0.20691560E-03,
     +  0.31179571E-03, 0.30579101E-03, 0.12038780E-02, 0.52286946E-03,
     + -0.19034573E-03,-0.64040138E-03,-0.48483707E-03, 0.54312736E-03,
     +  0.37444496E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_sp_cfp    =t_sp_cfp    
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21        *x52
     2  +coeff( 11)        *x32        
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)*x11*x23            
     7  +coeff( 16)        *x31        
     8  +coeff( 17)    *x21*x31        
      t_sp_cfp    =t_sp_cfp    
     9  +coeff( 18)        *x31    *x51
     1  +coeff( 19)        *x31*x41*x51
     2  +coeff( 20)                *x53
     3  +coeff( 21)        *x31    *x52
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)        *x32    *x52
     7  +coeff( 25)*x11    *x32    *x51
c
      return
      end
      function y_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 49)
      data ncoeff/ 48/
      data avdat/ -0.2472616E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.25324011E-02, 0.54874353E-01,-0.74935153E-01, 0.65161306E-02,
     +  0.20403274E-03, 0.15379499E-02,-0.11169595E-02,-0.38188878E-01,
     +  0.73053334E-02, 0.11033195E-03,-0.28280267E-02,-0.79136696E-02,
     +  0.19509597E-01, 0.92939759E-03, 0.70203480E-03,-0.10907941E-02,
     +  0.42914487E-02,-0.24149022E-02, 0.11470974E-01,-0.86742695E-02,
     +  0.16978629E-03, 0.23175560E-04, 0.84622558E-02,-0.58212876E-03,
     + -0.40500718E-02,-0.42555882E-02,-0.93718932E-04,-0.15745849E-03,
     + -0.13611107E-02,-0.90437272E-03, 0.10176701E-01,-0.16275435E-02,
     + -0.46384535E-02,-0.82876289E-03, 0.60880007E-02,-0.13649367E-02,
     + -0.15912807E-02, 0.65527204E-03,-0.79216884E-03, 0.46465849E-03,
     + -0.33368983E-02,-0.22231361E-03, 0.14839014E-02, 0.46806070E-02,
     +  0.28310288E-02, 0.14846050E-02, 0.27236410E-02, 0.16819990E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sp_cfp    =y_sp_cfp    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11                
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)        *x33        
     6  +coeff( 15)            *x43    
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)        *x31    *x52
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)        *x34        
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x22*x31        
     6  +coeff( 24)    *x22    *x41    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)        *x31    *x53
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 27)*x12                
     1  +coeff( 28)    *x21        *x53
     2  +coeff( 29)    *x22*x31    *x51
     3  +coeff( 30)    *x21        *x51
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)        *x33    *x51
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)    *x23*x31        
     8  +coeff( 35)        *x31    *x54
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 36)        *x32        
     1  +coeff( 37)    *x22            
     2  +coeff( 38)                *x53
     3  +coeff( 39)        *x32    *x52
     4  +coeff( 40)    *x22        *x51
     5  +coeff( 41)            *x41*x53
     6  +coeff( 42)    *x21*x31    *x52
     7  +coeff( 43)        *x33    *x52
     8  +coeff( 44)    *x21*x31    *x53
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 45)    *x22*x31    *x52
     1  +coeff( 46)    *x23*x31    *x51
     2  +coeff( 47)    *x22    *x41*x53
     3  +coeff( 48)*x11*x23*x32*x41    
c
      return
      end
      function p_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3799727E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24673766E-02,-0.76424633E-03, 0.48941325E-01,-0.39530519E-01,
     +  0.24721315E-02, 0.91065966E-04,-0.64128340E-03,-0.16824077E-03,
     +  0.13319228E-01, 0.31137269E-03,-0.37404575E-03,-0.24818364E-03,
     + -0.11412320E-01,-0.65512350E-02, 0.22818085E-02, 0.34377377E-02,
     +  0.41385768E-02, 0.27433899E-02,-0.59215759E-03, 0.79823635E-03,
     + -0.10972304E-02,-0.74051670E-03,-0.59768913E-03, 0.47687325E-03,
     + -0.27231063E-03, 0.14184703E-02,-0.17479007E-02,-0.21550404E-03,
     +  0.48868160E-03, 0.40263907E-03,-0.55388687E-03,-0.38206903E-03,
     + -0.25837104E-04, 0.69037720E-04, 0.97102042E-04,-0.19523528E-03,
     + -0.43063174E-03, 0.27371047E-03, 0.11452144E-02,-0.11011391E-02,
     +  0.35303583E-04, 0.83095103E-04, 0.27791795E-03,-0.32119919E-03,
     + -0.55853132E-03,-0.13559237E-02, 0.11897534E-02,-0.11779962E-03,
     + -0.13597155E-03,-0.84116138E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sp_cfp    =p_sp_cfp    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x22    *x41    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 18)    *x21*x31    *x51
     1  +coeff( 19)                *x52
     2  +coeff( 20)        *x33        
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)        *x31    *x53
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)            *x42*x53
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)        *x31*x42    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)        *x33    *x51
     2  +coeff( 29)    *x22*x31*x42    
     3  +coeff( 30)*x11*x23*x31        
     4  +coeff( 31)    *x22*x32*x42    
     5  +coeff( 32)        *x32        
     6  +coeff( 33)    *x21        *x52
     7  +coeff( 34)                *x53
     8  +coeff( 35)        *x33*x41    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)            *x41*x53
     2  +coeff( 38)    *x23*x31    *x51
     3  +coeff( 39)*x11*x23*x32*x43    
     4  +coeff( 40)*x11        *x41    
     5  +coeff( 41)    *x23            
     6  +coeff( 42)    *x21    *x42    
     7  +coeff( 43)        *x32    *x51
     8  +coeff( 44)        *x31*x41*x51
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 45)*x11*x21*x31        
     1  +coeff( 46)*x11    *x31    *x51
     2  +coeff( 47)*x11        *x41*x51
     3  +coeff( 48)    *x21*x32*x41    
     4  +coeff( 49)        *x32*x41*x51
     5  +coeff( 50)        *x32    *x52
c
      return
      end
      function l_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/ -0.2502072E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.32006E-01,-0.49852E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34876864E-02,-0.36202925E+00,-0.31350438E-01,-0.38877252E-01,
     +  0.99461097E-02,-0.32131806E-01, 0.49064424E-01, 0.29717379E-02,
     + -0.96776579E-02,-0.41558086E-02,-0.75049051E-02, 0.59020650E-02,
     + -0.13057715E-02, 0.29121260E-02, 0.34062632E-02,-0.12320568E-02,
     +  0.97379042E-02,-0.55274935E-02, 0.80361478E-02, 0.15053276E-02,
     + -0.25946300E-02,-0.27996779E-03,-0.54483209E-03, 0.48753670E-02,
     +  0.13352958E-02, 0.87431556E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)        *x31        
      l_sp_cfp    =l_sp_cfp    
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21        *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)                *x53
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)        *x32    *x51
     7  +coeff( 16)            *x41    
     8  +coeff( 17)        *x31*x41    
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 18)            *x42    
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x22*x31*x41*x52
     4  +coeff( 22)    *x21*x31        
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)*x11*x22            
     7  +coeff( 25)    *x21        *x53
     8  +coeff( 26)*x11    *x32    *x51
c
      return
      end
      function x_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.2593375E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16211534E-01, 0.42377925E+00, 0.19120258E-03, 0.13430269E+00,
     +  0.19436345E-01,-0.78342572E-01, 0.12735782E-01,-0.34477543E-02,
     + -0.90385228E-03, 0.73709548E-03,-0.12531520E-02,-0.66312624E-03,
     +  0.11328941E-02, 0.20930378E-02, 0.71578228E-03, 0.77924407E-02,
     + -0.10017459E-02,-0.83100196E-03,-0.17542307E-02, 0.32961165E-03,
     + -0.21103034E-02,-0.23834235E-02, 0.71445754E-03,-0.20737355E-02,
     +  0.16861954E-02,-0.12099701E-02, 0.50826243E-03,-0.48775150E-03,
     +  0.10240532E-03, 0.51404862E-03,-0.47068324E-03,-0.18258423E-02,
     + -0.81337616E-03, 0.19996399E-03,-0.10018842E-03,-0.79615781E-03,
     + -0.53553563E-03, 0.52978512E-03,-0.14252191E-02,-0.14542080E-02,
     +  0.54757624E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)                *x52
      x_sp_cdex   =x_sp_cdex   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x22        *x51
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)    *x24            
     2  +coeff( 20)        *x31        
     3  +coeff( 21)            *x42    
     4  +coeff( 22)    *x23*x31*x43    
     5  +coeff( 23)*x12*x21    *x43    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x21*x32*x43    
     8  +coeff( 26)*x11*x21            
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 27)        *x32    *x51
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)                *x53
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)    *x22*x32        
     5  +coeff( 32)*x11    *x31*x41    
     6  +coeff( 33)*x12*x21            
     7  +coeff( 34)*x12    *x31        
     8  +coeff( 35)        *x33    *x51
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 36)    *x21*x33*x41    
     1  +coeff( 37)*x11*x22*x32        
     2  +coeff( 38)    *x23*x33        
     3  +coeff( 39)*x12*x23*x31*x41    
     4  +coeff( 40)*x11    *x33*x43    
     5  +coeff( 41)        *x33        
c
      return
      end
      function t_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/  0.5329418E+00/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13298964E-02,-0.82902692E-01,-0.44232816E-03, 0.22777388E-01,
     + -0.59928573E-02,-0.11907567E-02, 0.76731918E-02,-0.24074095E-02,
     +  0.54766372E-03, 0.36357308E-02, 0.18770897E-02, 0.91196829E-03,
     + -0.13619513E-02,-0.25682573E-03,-0.25908719E-02,-0.16186273E-02,
     +  0.93298027E-03, 0.21310372E-02,-0.43290312E-03, 0.20543138E-03,
     +  0.14292329E-03,-0.38892918E-03, 0.90025773E-04, 0.11038687E-02,
     + -0.39486375E-03, 0.28184123E-03,-0.87754997E-04,-0.22574209E-03,
     +  0.21615354E-03, 0.26667197E-03, 0.18394212E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)        *x32        
      t_sp_cdex   =t_sp_cdex   
     9  +coeff(  9)        *x31        
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x32    *x51
     4  +coeff( 13)            *x42    
     5  +coeff( 14)        *x31    *x51
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x21    *x42    
      t_sp_cdex   =t_sp_cdex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x21*x32        
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)        *x31*x41*x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)*x11*x22            
     7  +coeff( 25)    *x22*x32        
     8  +coeff( 26)    *x21*x31*x41*x51
      t_sp_cdex   =t_sp_cdex   
     9  +coeff( 27)            *x43*x51
     1  +coeff( 28)        *x32    *x52
     2  +coeff( 29)*x11*x22        *x51
     3  +coeff( 30)    *x22*x32    *x51
     4  +coeff( 31)*x11*x23        *x52
c
      return
      end
      function y_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.1162915E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.70422021E-03,-0.13913377E+00, 0.11614332E+00,-0.65059168E-02,
     +  0.49692844E-02, 0.59906585E-03,-0.23220663E-03, 0.13768341E-02,
     +  0.18951179E-01, 0.28676642E-01,-0.56734501E-03, 0.46529403E-03,
     + -0.55417325E-02,-0.31921882E-01,-0.25514341E-02,-0.17611361E-02,
     + -0.47506848E-02,-0.15246538E-01,-0.22345500E-03, 0.14917030E-02,
     + -0.51063728E-02, 0.51955006E-03, 0.26279925E-02,-0.18820587E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cdex   =y_sp_cdex   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)    *x21*x31    *x51
     8  +coeff( 17)    *x22*x31        
      y_sp_cdex   =y_sp_cdex   
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)*x12                
     2  +coeff( 20)        *x33        
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)    *x22            
     6  +coeff( 24)    *x21    *x41*x51
c
      return
      end
      function p_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 50)
      data ncoeff/ 49/
      data avdat/  0.3732368E-02/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.74530346E-03,-0.34607064E-01, 0.15730698E-01,-0.12183005E-02,
     + -0.10791004E-02,-0.81352219E-02,-0.10528385E-04, 0.16240068E-03,
     +  0.27110886E-02, 0.68900776E-02,-0.80436155E-04, 0.82913938E-03,
     +  0.10660802E-02,-0.42207655E-02, 0.34077250E-03,-0.36675230E-03,
     +  0.11922464E-02,-0.19248786E-02, 0.31787640E-03,-0.20250844E-03,
     + -0.34230645E-03, 0.98112144E-03, 0.98891382E-03,-0.17287077E-03,
     +  0.63106592E-04,-0.35433861E-03, 0.10021487E-03,-0.31032178E-03,
     + -0.31009357E-03,-0.64588251E-03, 0.45052173E-04, 0.51383855E-03,
     + -0.28940657E-03, 0.42702704E-04, 0.11911245E-03,-0.23961350E-03,
     + -0.11524865E-03,-0.34511904E-03, 0.21265962E-03, 0.13674652E-03,
     +  0.62894011E-04,-0.72516654E-04,-0.69854541E-04,-0.52278512E-03,
     +  0.57422329E-03,-0.53515618E-04, 0.45604666E-03,-0.92038978E-03,
     +  0.19186428E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x21            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22            
     7  +coeff( 16)        *x33        
     8  +coeff( 17)    *x21*x31    *x51
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x22*x31    *x51
     4  +coeff( 22)*x11    *x31        
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)*x11                
     7  +coeff( 25)    *x21        *x51
     8  +coeff( 26)        *x31*x42    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)    *x21*x31    *x52
     3  +coeff( 30)    *x22*x31*x42    
     4  +coeff( 31)    *x23        *x52
     5  +coeff( 32)        *x31*x42*x52
     6  +coeff( 33)*x11*x23*x31        
     7  +coeff( 34)    *x22*x33*x41    
     8  +coeff( 35)        *x32        
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 36)            *x43    
     1  +coeff( 37)            *x41*x52
     2  +coeff( 38)*x11    *x31    *x51
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x21*x31*x42    
     5  +coeff( 41)        *x33    *x51
     6  +coeff( 42)        *x31    *x53
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)    *x23*x31    *x51
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)        *x33*x41*x51
     2  +coeff( 47)        *x33    *x52
     3  +coeff( 48)        *x32*x41*x52
     4  +coeff( 49)*x11*x22*x31    *x52
c
      return
      end
      function l_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.1849175E-01/
      data xmin/
     1 -0.84036E-01,-0.30564E-01,-0.58615E-01,-0.35058E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.93675E-01, 0.33188E-01, 0.40902E-01, 0.28301E-01, 0.49972E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.14899017E-01,-0.57504714E+00, 0.29868135E-03,-0.99463649E-01,
     + -0.36236912E-01,-0.77394894E-02, 0.89767598E-01,-0.88783045E-03,
     +  0.25940367E-02,-0.81184454E-03,-0.28627813E-02, 0.17931744E-02,
     + -0.27835988E-02, 0.58519086E-02, 0.72518981E-03, 0.72211958E-02,
     + -0.13373333E-01,-0.14197760E-02, 0.42637167E-03, 0.73230662E-03,
     +  0.71067957E-03, 0.70519571E-03, 0.11215328E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11                
     8  +coeff(  8)        *x32        
      l_sp_cdex   =l_sp_cdex   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)        *x31        
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)*x11            *x51
      l_sp_cdex   =l_sp_cdex   
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)            *x43    
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)*x11*x23            
c
      return
      end
