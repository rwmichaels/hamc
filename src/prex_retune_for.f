C Forward transfer functions for right hrs with septum based on prex_retune_dir.dat
c  setup for (x|theta)=(y|phi)=0 60 cm downstream of the 1st VDC
c HRS + PREX room temperature septum (right side)
c                     -JJL 4/22/09


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
      function x_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1444055E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35427525E-02, 0.62848264E-02,-0.11982995E-02, 0.73272508E+00,
     +  0.58588479E-01,-0.59608590E-01,-0.10824529E-01, 0.27201604E-01,
     + -0.15164130E-01, 0.89448327E-02,-0.73907617E-02,-0.16616251E-01,
     + -0.31673145E-02,-0.30547222E-02, 0.80417832E-02,-0.84229968E-02,
     + -0.99614197E-02,-0.23861548E-02, 0.64678136E-02, 0.12564475E-01,
     +  0.41055956E-03, 0.66389688E-02,-0.96120866E-03,-0.51378994E-03,
     +  0.86318096E-03,-0.68984972E-03,-0.40025050E-02, 0.36356123E-02,
     + -0.69083311E-02, 0.32943534E-02, 0.49244375E-02, 0.10001199E-01,
     +  0.76920941E-03, 0.37337267E-02,-0.29966110E-04,-0.18583935E-02,
     +  0.39941585E-02,-0.31726749E-02,-0.81190892E-03, 0.56779431E-03,
     + -0.34478249E-02, 0.21966072E-02, 0.29632787E-02, 0.28074055E-02,
     +  0.79924765E-03,-0.14093436E-02, 0.37411253E-02, 0.40158461E-03,
     + -0.14745044E-03,-0.25321747E-03, 0.12500750E-03, 0.78833802E-03,
     +  0.22336876E-03, 0.91965921E-03,-0.56247006E-03,-0.15412881E-02,
     + -0.18561464E-02, 0.97659603E-03, 0.35473965E-02, 0.61088428E-03,
     +  0.26258692E-03, 0.65346691E-03,-0.32651077E-02,-0.26235231E-02,
     + -0.57763653E-03,-0.48457668E-03,-0.12915778E-02, 0.70606414E-02,
     + -0.30658722E-02,-0.16675411E-02, 0.14221204E-02,-0.33942919E-03,
     +  0.14071040E-03, 0.35094717E-03, 0.54471107E-03,-0.20139008E-02,
     + -0.68820489E-03, 0.30805604E-03,-0.41480103E-03, 0.54229947E-03,
     +  0.24512534E-02, 0.99767058E-03,-0.12589167E-02, 0.40767979E-03,
     +  0.12759182E-02, 0.35615481E-03, 0.16666901E-02,-0.59010269E-03,
     + -0.28896451E-02, 0.23844285E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      x_sp_fp     =x_sp_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x24            
      x_sp_fp     =x_sp_fp     
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)    *x23        *x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)*x11*x22    *x41    
     2  +coeff( 29)    *x24        *x51
     3  +coeff( 30)            *x42*x52
     4  +coeff( 31)    *x23    *x41*x51
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)    *x22    *x42    
     8  +coeff( 35)    *x22    *x41*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x23        *x52
     3  +coeff( 39)            *x43    
     4  +coeff( 40)    *x21    *x42*x51
     5  +coeff( 41)    *x24    *x41    
     6  +coeff( 42)*x11*x22        *x51
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)    *x22    *x42*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 45)*x11*x22    *x42    
     1  +coeff( 46)            *x42*x53
     2  +coeff( 47)    *x23    *x42*x51
     3  +coeff( 48)*x11*x21        *x51
     4  +coeff( 49)*x11        *x42    
     5  +coeff( 50)*x11        *x41*x51
     6  +coeff( 51)*x12*x21            
     7  +coeff( 52)*x11*x23            
     8  +coeff( 53)    *x21    *x41*x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 54)            *x41*x53
     1  +coeff( 55)*x11*x23    *x41    
     2  +coeff( 56)    *x21    *x43*x51
     3  +coeff( 57)    *x21    *x42*x52
     4  +coeff( 58)*x11*x22    *x41*x51
     5  +coeff( 59)    *x24    *x41*x51
     6  +coeff( 60)    *x24        *x52
     7  +coeff( 61)            *x43*x52
     8  +coeff( 62)    *x23        *x53
      x_sp_fp     =x_sp_fp     
     9  +coeff( 63)*x11*x24    *x41    
     1  +coeff( 64)*x11*x24        *x51
     2  +coeff( 65)    *x22    *x42*x52
     3  +coeff( 66)*x12*x24            
     4  +coeff( 67)*x11*x22    *x42*x51
     5  +coeff( 68)    *x24    *x41*x52
     6  +coeff( 69)*x11*x23    *x41*x52
     7  +coeff( 70)*x11*x24    *x41*x52
     8  +coeff( 71)*x12*x24        *x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 72)            *x41*x52
     1  +coeff( 73)*x11            *x52
     2  +coeff( 74)*x11*x21    *x42    
     3  +coeff( 75)*x11*x21        *x52
     4  +coeff( 76)    *x22    *x41*x52
     5  +coeff( 77)*x11*x23        *x51
     6  +coeff( 78)    *x21    *x41*x53
     7  +coeff( 79)*x11*x21    *x42*x51
     8  +coeff( 80)*x11*x21    *x41*x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 81)    *x23    *x41*x52
     1  +coeff( 82)    *x22    *x43*x51
     2  +coeff( 83)*x11*x23        *x52
     3  +coeff( 84)    *x21    *x43*x52
     4  +coeff( 85)    *x24        *x53
     5  +coeff( 86)*x11*x21    *x43*x51
     6  +coeff( 87)*x11*x24    *x42    
     7  +coeff( 88)*x11*x24        *x52
     8  +coeff( 89)    *x22    *x42*x53
      x_sp_fp     =x_sp_fp     
     9  +coeff( 90)*x12*x21        *x53
c
      return
      end
      function t_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.2373933E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.68322918E-03,-0.18201835E-01, 0.15316324E-03, 0.11588956E+00,
     +  0.66362740E-02,-0.10697904E-01,-0.11075734E-02, 0.13132613E-02,
     + -0.10398206E-02,-0.31648777E-03, 0.42700252E-03,-0.13551463E-02,
     +  0.20865427E-03,-0.11498965E-02,-0.74272306E-03, 0.56239375E-03,
     +  0.22735087E-03, 0.16611398E-03, 0.28656831E-03,-0.49690588E-03,
     +  0.92302985E-03, 0.55087556E-03,-0.39714118E-03, 0.87799341E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_sp_fp     =t_sp_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x21    *x41    
      t_sp_fp     =t_sp_fp     
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)            *x42*x52
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)    *x23    *x41*x51
c
      return
      end
      function y_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 46)
      data ncoeff/ 45/
      data avdat/ -0.1689464E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.59307308E-03,-0.17907152E-03,-0.70862663E-02, 0.13493743E-02,
     + -0.28802773E-02,-0.36622718E-01,-0.35172896E-02, 0.13328217E-01,
     +  0.52113170E-02, 0.19971046E-02, 0.11572900E-02, 0.48930766E-02,
     + -0.27787997E-03, 0.12394955E-01, 0.21243270E-02, 0.65945199E-03,
     + -0.10960817E-02, 0.13593372E-01, 0.47669895E-02,-0.61762556E-02,
     + -0.46853004E-02, 0.32635869E-02,-0.18108850E-02, 0.25343394E-02,
     + -0.18176943E-02, 0.52492315E-03,-0.22372082E-02, 0.53048199E-02,
     + -0.99855871E-03,-0.11188334E-02,-0.97918138E-03,-0.19532444E-03,
     +  0.10246574E-02, 0.23553103E-03,-0.21258872E-02,-0.13338678E-02,
     +  0.63553912E-03, 0.24179127E-02, 0.45077997E-03, 0.22245350E-02,
     + -0.14420411E-02, 0.17891169E-02, 0.18036237E-02, 0.12184870E-02,
     +  0.33920491E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_fp     =y_sp_fp     
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x43    
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)    *x22            
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)                *x53
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x21        *x52
      y_sp_fp     =y_sp_fp     
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)            *x41*x53
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)    *x21    *x41*x53
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)            *x44    
      y_sp_fp     =y_sp_fp     
     9  +coeff( 27)                *x54
     1  +coeff( 28)            *x41*x54
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)*x11        *x41*x51
     6  +coeff( 33)    *x21    *x41*x52
     7  +coeff( 34)*x11*x21        *x51
     8  +coeff( 35)    *x21        *x53
      y_sp_fp     =y_sp_fp     
     9  +coeff( 36)    *x23    *x41    
     1  +coeff( 37)*x11*x21    *x42    
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)*x11*x21    *x41*x51
     4  +coeff( 40)    *x22    *x42*x51
     5  +coeff( 41)    *x21        *x54
     6  +coeff( 42)    *x22    *x41*x52
     7  +coeff( 43)    *x23    *x41*x51
     8  +coeff( 44)    *x23        *x52
      y_sp_fp     =y_sp_fp     
     9  +coeff( 45)    *x23        *x53
c
      return
      end
      function p_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3011403E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10179573E-02, 0.19467058E-01,-0.26463510E-02, 0.96691074E-03,
     +  0.12601298E-01,-0.13886375E-02, 0.32001238E-02,-0.21297211E-01,
     + -0.35653694E-02,-0.30540728E-02, 0.89558875E-02, 0.35977152E-02,
     +  0.41679987E-02,-0.28908029E-03,-0.10268216E-02, 0.12439637E-02,
     +  0.16789044E-02, 0.42445189E-03,-0.36998640E-03, 0.37447235E-03,
     +  0.47138505E-03, 0.14733167E-02,-0.11673782E-02, 0.37198249E-03,
     +  0.94683457E-03, 0.20965922E-02,-0.72146836E-03,-0.31269342E-03,
     +  0.17803954E-03, 0.36583759E-03, 0.10015860E-02, 0.39640957E-03,
     + -0.31118063E-03, 0.64677806E-04,-0.64438168E-03,-0.22541733E-03,
     + -0.83676283E-03, 0.26811313E-03,-0.25585345E-04,-0.22386400E-03,
     +  0.26292310E-03,-0.22817276E-03, 0.23544571E-03,-0.20970423E-03,
     +  0.60336333E-04, 0.13338634E-02, 0.54325827E-03, 0.10310099E-02,
     + -0.17394059E-03,-0.53056027E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41*x51
      p_sp_fp     =p_sp_fp     
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x21            
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x21    *x41*x51
      p_sp_fp     =p_sp_fp     
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)            *x41*x53
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x22    *x42*x51
      p_sp_fp     =p_sp_fp     
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)            *x43*x51
     2  +coeff( 29)    *x21    *x41*x52
     3  +coeff( 30)    *x21        *x53
     4  +coeff( 31)    *x22    *x43    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)*x11*x22    *x42*x51
     7  +coeff( 34)*x11                
     8  +coeff( 35)    *x21        *x52
      p_sp_fp     =p_sp_fp     
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)    *x22        *x52
     4  +coeff( 40)*x11*x23            
     5  +coeff( 41)*x11*x22    *x41    
     6  +coeff( 42)*x11*x22        *x51
     7  +coeff( 43)*x11*x21    *x41*x51
     8  +coeff( 44)*x11*x21        *x52
      p_sp_fp     =p_sp_fp     
     9  +coeff( 45)*x12        *x42    
     1  +coeff( 46)    *x23        *x52
     2  +coeff( 47)    *x21    *x42*x52
     3  +coeff( 48)    *x22        *x53
     4  +coeff( 49)*x11*x22    *x42    
     5  +coeff( 50)*x11*x23        *x51
c
      return
      end
      function l_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/ -0.1455364E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10633183E-01,-0.29081002E+00,-0.33077829E-01, 0.89993486E-02,
     + -0.38628705E-01, 0.38061101E-01,-0.32576084E-01,-0.16670607E-01,
     + -0.49122409E-02, 0.27858100E-02,-0.51958254E-02, 0.17165311E-01,
     +  0.60553229E-02,-0.23913218E-01,-0.25816543E-02,-0.38734644E-02,
     +  0.38523837E-02,-0.20214869E-02, 0.78065228E-03,-0.31097766E-02,
     + -0.57951934E-02,-0.10871630E-01,-0.38344825E-02, 0.16814813E-02,
     + -0.30197506E-03,-0.19029408E-02, 0.90532156E-03,-0.22932959E-02,
     +  0.10391430E-02,-0.27015470E-02,-0.29224477E-03,-0.10909261E-02,
     + -0.28315337E-01, 0.93457147E-05,-0.10876826E-02, 0.85272826E-02,
     + -0.12012064E-02,-0.12631958E-02, 0.11270032E-02,-0.52333623E-02,
     +  0.15817633E-02,-0.51185037E-02, 0.36900791E-02, 0.33895415E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x23            
      l_sp_fp     =l_sp_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)                *x53
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)            *x42*x51
      l_sp_fp     =l_sp_fp     
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)*x11        *x43    
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)*x11*x22    *x42    
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x22        *x51
      l_sp_fp     =l_sp_fp     
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)*x11*x21    *x42    
     5  +coeff( 32)*x11*x23    *x41    
     6  +coeff( 33)    *x23    *x43    
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)*x11*x21    *x41    
      l_sp_fp     =l_sp_fp     
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x22    *x41*x51
     2  +coeff( 38)            *x42*x52
     3  +coeff( 39)    *x21        *x53
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)            *x43*x52
     6  +coeff( 42)*x11*x22    *x43    
     7  +coeff( 43)*x11*x23    *x41*x52
     8  +coeff( 44)*x11*x22    *x42*x52
c
      return
      end
      function x_sp_col     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1742675E-02/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51386E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14639602E-02, 0.12083183E+00, 0.37593034E-02,-0.34698253E-02,
     +  0.18215601E-02, 0.11529680E-02, 0.22849545E-02,-0.55656658E-03,
     +  0.14193709E-02,-0.13998586E-02, 0.16592095E-03, 0.29949626E-03,
     + -0.52654702E-03, 0.13425618E-02, 0.15349494E-03,-0.22380290E-03,
     +  0.53985475E-03,-0.79711695E-03,-0.48243388E-03, 0.34070163E-04,
     + -0.25050470E-03, 0.15173425E-03, 0.20741727E-03, 0.90932421E-03,
     + -0.16847918E-03, 0.11282061E-04,-0.97753713E-04,-0.79200854E-05,
     +  0.34259243E-04,-0.30850916E-03, 0.84795618E-04,-0.21982715E-03,
     +  0.36001974E-03, 0.42870632E-03,-0.13749819E-03,-0.21210271E-03,
     +  0.15138500E-03,-0.11947990E-03, 0.71531672E-05,-0.69963447E-04,
     +  0.33384958E-05, 0.64325890E-04, 0.33440214E-04, 0.55385125E-03,
     + -0.53681688E-04,-0.90506655E-04,-0.14935229E-03,-0.31222808E-03,
     +  0.20837912E-03, 0.63782500E-04,-0.26946445E-03,-0.56803707E-04,
     + -0.32788891E-03,-0.27858748E-03, 0.30577858E-03,-0.75122975E-04,
     +  0.37752237E-03,-0.73560259E-05,-0.30778872E-05,-0.44252215E-05,
     +  0.49585201E-05, 0.15011069E-03,-0.20584346E-04, 0.74927821E-05,
     +  0.85072941E-04, 0.12868551E-03, 0.10040370E-03,-0.43384272E-04,
     +  0.80925283E-05, 0.48844908E-04, 0.23204320E-04,-0.94990210E-04,
     + -0.23691184E-05,-0.14468087E-03,-0.19603853E-04,-0.33645720E-04,
     + -0.20219291E-04, 0.29624362E-04, 0.41808489E-04, 0.10906631E-03,
     + -0.17599456E-03, 0.13526402E-04,-0.19403124E-03,-0.45396137E-04,
     +  0.93772040E-04, 0.11191967E-03,-0.25223728E-04, 0.72991206E-04,
     +  0.51440835E-04,-0.14634305E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_col     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)            *x41    
      x_sp_col     =x_sp_col     
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)*x11*x22    *x41    
      x_sp_col     =x_sp_col     
     9  +coeff( 18)    *x24    *x41    
     1  +coeff( 19)    *x24    *x42    
     2  +coeff( 20)                *x51
     3  +coeff( 21)            *x42    
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)*x11*x23            
     8  +coeff( 26)            *x41*x51
      x_sp_col     =x_sp_col     
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x12*x21            
     3  +coeff( 30)    *x21    *x43    
     4  +coeff( 31)    *x24        *x51
     5  +coeff( 32)*x11*x23    *x41    
     6  +coeff( 33)*x11*x22    *x42    
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)*x11*x21    *x41*x52
      x_sp_col     =x_sp_col     
     9  +coeff( 36)*x11*x24    *x41    
     1  +coeff( 37)*x11*x21    *x42*x52
     2  +coeff( 38)            *x43    
     3  +coeff( 39)            *x42*x51
     4  +coeff( 40)*x11        *x42    
     5  +coeff( 41)*x12        *x41    
     6  +coeff( 42)    *x23    *x41*x51
     7  +coeff( 43)*x11*x24            
     8  +coeff( 44)    *x22    *x43    
      x_sp_col     =x_sp_col     
     9  +coeff( 45)*x11        *x42*x52
     1  +coeff( 46)*x12*x21    *x41*x51
     2  +coeff( 47)*x12*x24            
     3  +coeff( 48)    *x24    *x43    
     4  +coeff( 49)    *x24    *x41*x52
     5  +coeff( 50)*x11*x21    *x43*x51
     6  +coeff( 51)    *x23    *x43*x51
     7  +coeff( 52)*x11*x21    *x41*x53
     8  +coeff( 53)*x11*x24    *x42    
      x_sp_col     =x_sp_col     
     9  +coeff( 54)    *x22    *x43*x52
     1  +coeff( 55)*x11*x22    *x42*x52
     2  +coeff( 56)*x12*x23    *x42    
     3  +coeff( 57)*x11*x21    *x43*x52
     4  +coeff( 58)    *x21        *x52
     5  +coeff( 59)*x12                
     6  +coeff( 60)*x11        *x41*x51
     7  +coeff( 61)*x11            *x52
     8  +coeff( 62)    *x21    *x42*x51
      x_sp_col     =x_sp_col     
     9  +coeff( 63)    *x21    *x41*x52
     1  +coeff( 64)            *x43*x51
     2  +coeff( 65)*x12*x22            
     3  +coeff( 66)*x12*x21    *x41    
     4  +coeff( 67)    *x21    *x43*x51
     5  +coeff( 68)    *x21    *x42*x52
     6  +coeff( 69)*x11*x22    *x41*x51
     7  +coeff( 70)*x11*x22        *x52
     8  +coeff( 71)            *x43*x52
      x_sp_col     =x_sp_col     
     9  +coeff( 72)*x11*x21    *x43    
     1  +coeff( 73)*x11*x21    *x42*x51
     2  +coeff( 74)    *x23    *x42*x51
     3  +coeff( 75)*x11*x21        *x53
     4  +coeff( 76)*x12*x22    *x41    
     5  +coeff( 77)*x11*x23    *x42    
     6  +coeff( 78)    *x24    *x42*x51
     7  +coeff( 79)*x11*x22    *x41*x52
     8  +coeff( 80)*x11*x24    *x41*x51
      x_sp_col     =x_sp_col     
     9  +coeff( 81)*x11*x24        *x52
     1  +coeff( 82)*x11        *x43*x52
     2  +coeff( 83)*x12*x21    *x43    
     3  +coeff( 84)*x12*x21    *x42*x51
     4  +coeff( 85)*x11*x23    *x42*x51
     5  +coeff( 86)*x11*x23    *x41*x52
     6  +coeff( 87)*x12*x24        *x51
     7  +coeff( 88)    *x24    *x43*x51
     8  +coeff( 89)    *x24    *x42*x52
      x_sp_col     =x_sp_col     
     9  +coeff( 90)*x11*x22    *x41*x53
c
      return
      end
      function t_sp_col     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.9818955E-03/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51386E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.39164562E-03, 0.45146391E-01,-0.43272334E-02, 0.68456784E-03,
     +  0.16372181E-02, 0.29426806E-02,-0.56328875E-03,-0.31039360E-03,
     +  0.97406015E-03, 0.29067631E-03,-0.13504670E-02, 0.30837976E-03,
     +  0.13359722E-02, 0.34260232E-04, 0.56982888E-04,-0.23415334E-03,
     +  0.17754294E-03, 0.85170475E-04, 0.36669069E-03,-0.23499948E-03,
     + -0.84339503E-04, 0.42721210E-03, 0.16355139E-03,-0.15379590E-03,
     +  0.32453150E-04, 0.25334847E-03, 0.25664826E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_col     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x23    *x41    
     7  +coeff(  7)            *x41    
     8  +coeff(  8)*x11                
      t_sp_col     =t_sp_col     
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)                *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x21    *x41*x51
      t_sp_col     =t_sp_col     
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)            *x42*x52
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x22    *x42*x52
     6  +coeff( 24)            *x42    
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)*x11*x21    *x42*x52
      t_sp_col     =t_sp_col     
     9  +coeff( 27)*x11*x22    *x42*x52
c
      return
      end
      function y_sp_col     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.8833224E-02/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51386E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.19098083E-03, 0.55173006E-01, 0.46393820E-02,-0.13540791E-02,
     + -0.34628164E-02,-0.21711951E-02,-0.47640497E-03, 0.95080351E-03,
     +  0.67514722E-03,-0.37292662E-03,-0.11770079E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      y_sp_col     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      y_sp_col     =y_sp_col     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22    *x42    
c
      return
      end
      function p_sp_col     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4133915E-02/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51386E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.95886510E-03,-0.15253081E-02, 0.24188867E-01, 0.66525913E-02,
     + -0.43738089E-02,-0.20828727E-02,-0.84126234E-03,-0.30171633E-03,
     +  0.11672101E-02, 0.63826027E-03,-0.37063658E-03,-0.12548680E-02,
     + -0.50550140E-03,-0.53737120E-04,-0.31538118E-03, 0.60138159E-03,
     + -0.24139194E-04,-0.27859426E-03, 0.18631388E-03,-0.97172931E-04,
     + -0.91655216E-04,-0.78338030E-03, 0.46515340E-03,-0.47494835E-03,
     +  0.98707282E-03,-0.14480040E-03,-0.38746886E-04, 0.11743413E-03,
     +  0.57419813E-04, 0.20085499E-03,-0.78635210E-04,-0.21753152E-03,
     +  0.30275271E-03,-0.12518071E-03,-0.97025011E-04, 0.16574826E-03,
     +  0.30064778E-03, 0.18498460E-03,-0.26406872E-03, 0.25453273E-03,
     + -0.11093831E-03, 0.18558576E-04,-0.53360913E-03, 0.10312002E-03,
     +  0.15703362E-03, 0.46844903E-03,-0.21197395E-03, 0.30684622E-03,
     +  0.28936728E-03, 0.58840052E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_col     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      p_sp_col     =p_sp_col     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)*x11*x21    *x41    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)            *x41*x51
      p_sp_col     =p_sp_col     
     9  +coeff( 18)            *x43    
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11            *x52
     3  +coeff( 21)*x11*x23            
     4  +coeff( 22)*x11*x21    *x42    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)    *x23    *x41*x51
     7  +coeff( 25)*x11*x23    *x42    
     8  +coeff( 26)*x11        *x42*x53
      p_sp_col     =p_sp_col     
     9  +coeff( 27)*x11                
     1  +coeff( 28)    *x21        *x51
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)            *x41*x52
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)            *x42*x52
     7  +coeff( 34)            *x41*x53
     8  +coeff( 35)*x11        *x41*x52
      p_sp_col     =p_sp_col     
     9  +coeff( 36)    *x22    *x42*x51
     1  +coeff( 37)            *x43*x52
     2  +coeff( 38)    *x22        *x53
     3  +coeff( 39)            *x42*x53
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x21    *x43    
     6  +coeff( 42)    *x22    *x43*x51
     7  +coeff( 43)    *x22    *x42*x52
     8  +coeff( 44)*x12*x23        *x51
      p_sp_col     =p_sp_col     
     9  +coeff( 45)*x12        *x43*x51
     1  +coeff( 46)    *x22    *x42*x53
     2  +coeff( 47)*x11*x21    *x43*x52
     3  +coeff( 48)*x12*x22    *x43    
     4  +coeff( 49)*x12*x23    *x43    
     5  +coeff( 50)*x12*x23    *x42*x53
c
      return
      end
      function l_sp_col     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.4697743E-03/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51386E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12799997E-02, 0.39470913E-04,-0.48939385E-02,-0.12428229E-03,
     + -0.29927765E-02,-0.63647685E-03,-0.10367396E-03, 0.25641185E-03,
     +  0.84177751E-04,-0.55865210E-04,-0.74732008E-04, 0.89804395E-04,
     +  0.64400319E-05, 0.34818720E-04, 0.30278038E-05,-0.34325994E-05,
     + -0.88787165E-05,-0.71713785E-05, 0.40362393E-04, 0.15588685E-04,
     + -0.10016148E-04, 0.48764759E-05,-0.22143991E-04, 0.63728648E-05,
     +  0.51723218E-04,-0.99648878E-05, 0.10849750E-04,-0.48059734E-04,
     +  0.24949039E-04,-0.21083006E-04, 0.88711031E-05,-0.20629530E-04,
     +  0.57205525E-05,-0.11739577E-04, 0.23620520E-04,-0.11869878E-04,
     +  0.25374491E-04,-0.19519834E-04,-0.57373971E-04, 0.28471799E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_col     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_col     =l_sp_col     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)                *x52
     8  +coeff( 17)*x11*x22            
      l_sp_col     =l_sp_col     
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)*x11            *x52
     5  +coeff( 23)            *x42*x52
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)*x11*x21    *x42    
     8  +coeff( 26)*x11        *x42*x51
      l_sp_col     =l_sp_col     
     9  +coeff( 27)*x11        *x41*x52
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)    *x23    *x41*x51
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)            *x42*x53
     5  +coeff( 32)*x11*x23    *x41    
     6  +coeff( 33)*x11        *x43*x51
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)    *x23    *x42*x51
      l_sp_col     =l_sp_col     
     9  +coeff( 36)    *x22    *x43*x51
     1  +coeff( 37)    *x22    *x42*x52
     2  +coeff( 38)    *x21    *x42*x53
     3  +coeff( 39)*x11*x23    *x42    
     4  +coeff( 40)*x11        *x42*x53
c
      return
      end
      function x_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2221237E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19774621E-02, 0.13825060E+00,-0.77310065E-02, 0.17664997E-02,
     +  0.34362825E-02, 0.33739172E-02,-0.12433161E-02, 0.24179632E-02,
     +  0.30902829E-02, 0.41137300E-02,-0.33189578E-02, 0.74537122E-03,
     + -0.11197731E-02, 0.29257268E-02,-0.11261661E-03, 0.77529090E-04,
     + -0.48271799E-03, 0.27578216E-03, 0.21236483E-02,-0.13732993E-02,
     +  0.55209437E-03,-0.43480299E-03, 0.97307457E-04,-0.16851713E-03,
     +  0.38777577E-03,-0.16856990E-03,-0.27114505E-03, 0.10679500E-02,
     +  0.66189743E-02, 0.44070279E-04, 0.27186863E-03,-0.21104547E-02,
     + -0.20661163E-03,-0.74329058E-03,-0.19769667E-03, 0.35538444E-04,
     +  0.77702927E-04, 0.14870352E-03, 0.19484639E-03, 0.23171093E-03,
     +  0.88727550E-03,-0.39829192E-03,-0.12766054E-02,-0.34018388E-03,
     + -0.10192574E-03, 0.11208366E-03,-0.20040867E-03,-0.76593581E-03,
     + -0.17918331E-02,-0.19024438E-03, 0.19796416E-02, 0.25889075E-02,
     +  0.20370202E-02,-0.43835660E-03,-0.34343022E-04,-0.46445242E-04,
     + -0.46716034E-03, 0.10497288E-03,-0.17954038E-04,-0.38102444E-03,
     + -0.43159659E-03,-0.97133299E-04, 0.78840391E-03, 0.34375564E-03,
     +  0.96457085E-03, 0.85522392E-03, 0.49063208E-04, 0.11002002E-03,
     + -0.39408685E-03,-0.35474732E-03,-0.45832599E-03,-0.15595952E-02,
     +  0.51459621E-04,-0.35908823E-04, 0.32748361E-03, 0.26422480E-03,
     + -0.38189365E-03,-0.49992342E-03,-0.18158002E-03,-0.67163550E-04,
     + -0.82852901E-03,-0.18650660E-03,-0.73676615E-03, 0.94457489E-03,
     + -0.15162744E-02, 0.68561482E-03,-0.36224051E-05, 0.30621904E-04,
     +  0.31985037E-04, 0.85873769E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23            
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22            
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)*x11        *x43    
     7  +coeff( 16)                *x51
     8  +coeff( 17)            *x42    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x24    *x41    
     3  +coeff( 21)*x11*x22    *x42    
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11        *x42    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 27)*x11*x23            
     1  +coeff( 28)*x11*x22    *x41    
     2  +coeff( 29)    *x23    *x43    
     3  +coeff( 30)            *x41*x51
     4  +coeff( 31)    *x21    *x41*x51
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x24            
     7  +coeff( 34)*x11*x24    *x41    
     8  +coeff( 35)            *x43    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)*x12*x21            
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x11*x21    *x42    
     4  +coeff( 40)*x12*x22            
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)*x11*x23    *x41    
     7  +coeff( 43)    *x24    *x42    
     8  +coeff( 44)*x12*x24            
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 45)*x11*x22    *x42*x51
     1  +coeff( 46)*x11*x22    *x41*x52
     2  +coeff( 47)*x12*x23        *x51
     3  +coeff( 48)*x12*x22    *x41*x51
     4  +coeff( 49)*x12*x24    *x41    
     5  +coeff( 50)*x12*x21    *x42*x52
     6  +coeff( 51)*x12*x24    *x41*x51
     7  +coeff( 52)*x12*x23    *x42*x51
     8  +coeff( 53)*x11*x23    *x42*x53
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 54)*x12*x23    *x41*x53
     1  +coeff( 55)*x12        *x41    
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)    *x22    *x41*x52
     4  +coeff( 58)*x12*x21    *x41    
     5  +coeff( 59)*x11*x23        *x51
     6  +coeff( 60)    *x21    *x43*x51
     7  +coeff( 61)    *x24    *x41*x51
     8  +coeff( 62)*x11*x21        *x53
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 63)*x12*x22    *x41    
     1  +coeff( 64)*x11*x23    *x41*x51
     2  +coeff( 65)*x11*x22    *x43    
     3  +coeff( 66)    *x24    *x41*x52
     4  +coeff( 67)*x12*x23    *x41    
     5  +coeff( 68)*x11*x21    *x41*x53
     6  +coeff( 69)*x12*x22    *x42    
     7  +coeff( 70)    *x22    *x42*x53
     8  +coeff( 71)*x12*x21    *x42*x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 72)*x11*x23    *x42*x51
     1  +coeff( 73)*x12*x21        *x53
     2  +coeff( 74)*x11*x23        *x53
     3  +coeff( 75)    *x21    *x43*x53
     4  +coeff( 76)    *x24    *x41*x53
     5  +coeff( 77)*x12*x23    *x42    
     6  +coeff( 78)*x12*x23    *x41*x51
     7  +coeff( 79)*x11*x21    *x43*x52
     8  +coeff( 80)    *x23    *x42*x53
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 81)*x11*x23    *x41*x53
     1  +coeff( 82)*x11*x22    *x43*x52
     2  +coeff( 83)    *x23    *x43*x53
     3  +coeff( 84)*x12*x24    *x43    
     4  +coeff( 85)*x11*x24    *x43*x52
     5  +coeff( 86)*x11*x24    *x42*x53
     6  +coeff( 87)                *x52
     7  +coeff( 88)            *x42*x51
     8  +coeff( 89)            *x41*x52
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 90)*x11*x21        *x51
c
      return
      end
      function t_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 39)
      data ncoeff/ 38/
      data avdat/ -0.1493079E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.27774577E-03,-0.11465000E-01,-0.16236467E-02,-0.16903701E-02,
     +  0.27960106E-02,-0.20526171E-03, 0.23204116E-03, 0.33754535E-03,
     +  0.32174072E-03,-0.77359501E-03, 0.12047648E-02, 0.19596491E-04,
     +  0.69787173E-04,-0.16568693E-03,-0.14465212E-03, 0.79551268E-04,
     +  0.11031738E-04,-0.79712612E-04, 0.45308137E-04, 0.26035082E-03,
     +  0.10972311E-03, 0.20807538E-03, 0.47707304E-03, 0.12694637E-03,
     + -0.87783352E-04,-0.64136133E-04, 0.28302915E-04, 0.67047447E-04,
     + -0.32770029E-04,-0.36109850E-04, 0.21642254E-03, 0.10649214E-02,
     + -0.31788458E-03, 0.35027559E-04, 0.46855668E-04,-0.10888973E-03,
     + -0.10791676E-03,-0.67416528E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11        *x43    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 18)            *x42    
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x23        *x51
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)            *x43    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11        *x42    
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)    *x21    *x42*x51
     8  +coeff( 35)    *x22    *x42*x51
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 36)*x11*x22    *x41*x52
     1  +coeff( 37)*x11*x23    *x41*x52
     2  +coeff( 38)*x12*x22    *x41*x52
c
      return
      end
      function y_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.9565258E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.66634500E-02, 0.12371802E+00, 0.20155212E-01,-0.50948975E-02,
     +  0.31855351E-02,-0.13784861E-01,-0.43805414E-02,-0.18860981E-02,
     + -0.11843657E-02,-0.35310406E-02,-0.17503272E-02,-0.13640367E-02,
     +  0.37178318E-02,-0.98765716E-02, 0.10356137E-02, 0.35718910E-02,
     + -0.67380182E-02,-0.21785197E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)            *x41*x51
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x22    *x43    
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff( 18)*x11*x21    *x43    
c
      return
      end
      function p_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4042409E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24644516E-02,-0.25704440E-02, 0.51857434E-01, 0.10633092E-01,
     + -0.69846367E-02,-0.16075280E-02,-0.62156777E-03, 0.18267636E-02,
     + -0.17094245E-02,-0.96238358E-03, 0.15553412E-02,-0.77690848E-03,
     +  0.45848294E-03,-0.49180361E-02,-0.11874987E-02, 0.40639023E-03,
     + -0.64579194E-03,-0.45301379E-02,-0.65017179E-04,-0.11248686E-02,
     +  0.12131345E-03,-0.55222312E-03, 0.47254536E-03, 0.11853987E-02,
     +  0.10677902E-02,-0.46144170E-02, 0.24302285E-03,-0.15347204E-03,
     +  0.15747009E-03, 0.14355585E-03,-0.34155758E-03, 0.19030578E-02,
     + -0.22145147E-02,-0.35375440E-04,-0.18223419E-03,-0.36118960E-03,
     +  0.27739041E-03,-0.24874960E-04, 0.10572172E-02, 0.12160421E-03,
     + -0.42087614E-03, 0.34029491E-03, 0.24550839E-03, 0.98045180E-04,
     + -0.31243436E-03, 0.31161704E-03,-0.10476743E-02, 0.46878058E-03,
     + -0.40827107E-03,-0.23956064E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x23            
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x43    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)                *x52
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x21    *x42    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)*x11                
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)*x11*x23    *x41    
     8  +coeff( 26)*x11*x23    *x43*x53
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 27)    *x21        *x51
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)    *x23        *x51
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)*x12*x22    *x42*x51
     7  +coeff( 34)*x12                
     8  +coeff( 35)*x11            *x52
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)            *x42*x52
     2  +coeff( 38)*x11        *x43    
     3  +coeff( 39)    *x22    *x42*x51
     4  +coeff( 40)    *x21    *x41*x53
     5  +coeff( 41)*x11*x21    *x43    
     6  +coeff( 42)*x11*x22        *x52
     7  +coeff( 43)*x12*x23            
     8  +coeff( 44)*x12*x22    *x41    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 45)*x12*x21    *x42    
     1  +coeff( 46)*x12*x22        *x51
     2  +coeff( 47)    *x22    *x42*x52
     3  +coeff( 48)*x11*x22    *x41*x52
     4  +coeff( 49)*x11*x21    *x42*x52
     5  +coeff( 50)*x11        *x43*x52
c
      return
      end
      function l_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1101841E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17708678E-02,-0.18907538E-04,-0.43344027E-02,-0.35633126E-02,
     +  0.19218533E-03,-0.21375674E-02,-0.77067706E-03, 0.55269775E-03,
     + -0.78796307E-04, 0.39054622E-03,-0.65052744E-04, 0.11748135E-03,
     + -0.48324418E-04, 0.39774401E-04, 0.59979639E-05,-0.12452735E-03,
     +  0.57919828E-04, 0.71297960E-04, 0.77994191E-04,-0.95922842E-04,
     +  0.18956879E-03, 0.71369734E-03, 0.16308588E-03, 0.10527851E-04,
     + -0.72269281E-05,-0.15160431E-03, 0.12322361E-03,-0.37561238E-03,
     +  0.32657717E-04,-0.12977326E-05, 0.16877120E-04,-0.22286860E-04,
     +  0.60131537E-04,-0.33231143E-04,-0.86226828E-04, 0.45284189E-03,
     +  0.41468261E-05, 0.10655326E-05, 0.84303301E-05,-0.47956746E-05,
     + -0.28065704E-04,-0.18580226E-04,-0.16690397E-04,-0.17502973E-04,
     + -0.19473071E-04, 0.88579181E-05, 0.10175439E-04, 0.36291927E-04,
     + -0.19158215E-04, 0.46902314E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)                *x51
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)            *x43    
     8  +coeff( 17)            *x42*x51
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)*x11*x23            
     8  +coeff( 26)    *x23    *x42    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 27)    *x22    *x42*x51
     1  +coeff( 28)    *x23    *x43    
     2  +coeff( 29)*x11*x23    *x42*x53
     3  +coeff( 30)*x11                
     4  +coeff( 31)                *x53
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)            *x43*x52
     8  +coeff( 35)*x11*x23    *x41    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 36)*x11*x23    *x43*x53
     1  +coeff( 37)*x11        *x41    
     2  +coeff( 38)*x11            *x51
     3  +coeff( 39)    *x21    *x41*x51
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)    *x22        *x52
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)*x11*x22    *x41    
     8  +coeff( 44)*x11*x22        *x51
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)*x11        *x41*x52
     2  +coeff( 47)*x11            *x53
     3  +coeff( 48)    *x22    *x41*x52
     4  +coeff( 49)    *x22        *x53
     5  +coeff( 50)*x11*x21    *x43    
c
      return
      end
      function x_sp_q2ex  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4344766E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22455694E-02, 0.19236839E+00,-0.31389650E-02,-0.19501496E-01,
     +  0.11639248E-01, 0.56751673E-02, 0.70305276E-02,-0.52016799E-03,
     +  0.11787278E-01, 0.46038614E-02,-0.29710850E-02,-0.11221578E-02,
     +  0.72524557E-02,-0.89356462E-02, 0.16415343E-02,-0.26832414E-02,
     +  0.26889693E-03, 0.65187033E-03,-0.10310895E-02, 0.45505349E-03,
     + -0.81378268E-03, 0.24553731E-02,-0.28747392E-02, 0.51970989E-02,
     +  0.93372265E-03,-0.66306890E-03, 0.14961407E-01,-0.47314265E-02,
     + -0.48039301E-03, 0.27568310E-02,-0.31933186E-02, 0.13644910E-03,
     +  0.25577755E-02, 0.16050268E-03, 0.53624233E-03,-0.48982666E-03,
     +  0.14568951E-03, 0.61670278E-03, 0.40775546E-03,-0.91834762E-03,
     + -0.14215877E-02,-0.92289736E-03,-0.67720073E-03, 0.28580648E-03,
     + -0.62332593E-03,-0.41023479E-02, 0.27194546E-03,-0.45944205E-02,
     +  0.12127555E-03, 0.45168385E-03,-0.45555894E-03,-0.79010038E-04,
     +  0.43379111E-03, 0.26357215E-03,-0.12507328E-02,-0.30645475E-03,
     + -0.11739549E-02,-0.38016142E-03, 0.15611858E-02, 0.56792889E-03,
     + -0.11743586E-03, 0.19731449E-02,-0.53768605E-03, 0.98963072E-04,
     +  0.41662759E-03, 0.19478479E-02,-0.84089983E-03,-0.62813453E-03,
     + -0.57469610E-04, 0.39138980E-03,-0.89914154E-03,-0.93704093E-05,
     + -0.14258932E-02, 0.16610981E-02,-0.59108320E-03, 0.23130409E-02,
     +  0.53626634E-02, 0.10555582E-02, 0.18167197E-02,-0.22194492E-04,
     + -0.33701730E-04,-0.46033547E-04,-0.12488752E-03,-0.14378960E-03,
     + -0.69041816E-04,-0.91334754E-04,-0.86659915E-03,-0.18777154E-03,
     +  0.18539462E-03, 0.25928131E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q2ex  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)            *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x24            
     8  +coeff( 17)                *x51
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x23            
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 27)    *x23    *x43    
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)*x11*x24            
     3  +coeff( 30)    *x22    *x43    
     4  +coeff( 31)    *x24    *x42    
     5  +coeff( 32)*x11        *x42*x52
     6  +coeff( 33)*x11*x22    *x43    
     7  +coeff( 34)            *x41*x51
     8  +coeff( 35)    *x21    *x41*x51
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 36)*x11        *x42    
     1  +coeff( 37)*x12*x21            
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x12*x22            
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x24    *x41    
     6  +coeff( 42)*x11*x23    *x41*x51
     7  +coeff( 43)*x12*x24            
     8  +coeff( 44)*x11*x22    *x42*x51
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 45)*x11*x24    *x42    
     1  +coeff( 46)*x12*x24    *x41    
     2  +coeff( 47)*x12*x24        *x51
     3  +coeff( 48)*x11*x24    *x43*x52
     4  +coeff( 49)    *x22        *x51
     5  +coeff( 50)    *x23        *x51
     6  +coeff( 51)    *x22    *x41*x51
     7  +coeff( 52)*x12        *x41    
     8  +coeff( 53)*x11*x21    *x42    
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 54)*x11*x21    *x41*x51
     1  +coeff( 55)    *x23    *x41*x51
     2  +coeff( 56)*x11        *x43    
     3  +coeff( 57)    *x21    *x43*x51
     4  +coeff( 58)    *x21    *x42*x52
     5  +coeff( 59)*x11*x22    *x42    
     6  +coeff( 60)    *x23    *x42*x51
     7  +coeff( 61)*x11*x21    *x41*x52
     8  +coeff( 62)*x12*x22    *x41    
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 63)*x12*x22        *x51
     1  +coeff( 64)*x12        *x42*x51
     2  +coeff( 65)*x11*x22    *x41*x52
     3  +coeff( 66)    *x24    *x41*x52
     4  +coeff( 67)*x12*x23        *x51
     5  +coeff( 68)*x12*x22    *x41*x51
     6  +coeff( 69)*x12*x22        *x52
     7  +coeff( 70)    *x22    *x43*x52
     8  +coeff( 71)*x12*x21    *x42*x51
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 72)*x12        *x41*x53
     1  +coeff( 73)*x12*x23    *x41*x51
     2  +coeff( 74)    *x23    *x43*x52
     3  +coeff( 75)*x12*x21    *x42*x52
     4  +coeff( 76)*x12*x24    *x41*x51
     5  +coeff( 77)*x12*x23    *x42*x51
     6  +coeff( 78)*x12*x24        *x53
     7  +coeff( 79)*x11*x24    *x42*x53
     8  +coeff( 80)                *x53
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 81)*x11        *x41*x51
     1  +coeff( 82)*x11            *x52
     2  +coeff( 83)    *x22        *x52
     3  +coeff( 84)*x11*x22        *x51
     4  +coeff( 85)            *x41*x53
     5  +coeff( 86)*x11        *x42*x51
     6  +coeff( 87)    *x22    *x41*x52
     7  +coeff( 88)    *x22        *x53
     8  +coeff( 89)*x12*x21    *x41    
      x_sp_q2ex  =x_sp_q2ex  
     9  +coeff( 90)*x12*x21        *x51
c
      return
      end
      function t_sp_q2ex  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/ -0.1352613E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54975506E-03, 0.49597032E-01,-0.18997565E-02,-0.66241901E-02,
     +  0.30865707E-02, 0.21939178E-02,-0.87880349E-03, 0.11435297E-02,
     +  0.12456782E-02,-0.35148747E-02, 0.42363266E-02, 0.89334164E-04,
     +  0.17429247E-03,-0.36733504E-03, 0.47979018E-03,-0.47382477E-04,
     +  0.20196210E-02, 0.20676993E-02,-0.40257873E-03, 0.15549608E-03,
     + -0.35445386E-03,-0.29557844E-03, 0.14069846E-02, 0.65223424E-03,
     +  0.43612032E-04,-0.14362801E-03, 0.26296967E-03,-0.87245353E-04,
     +  0.51099770E-02, 0.25557072E-03, 0.13614201E-03,-0.15191090E-02,
     +  0.20112442E-03, 0.43575317E-04,-0.64969368E-04, 0.19984080E-03,
     +  0.21183814E-04,-0.11723672E-03,-0.21002487E-03, 0.74827526E-03,
     + -0.34463592E-03,-0.51387580E-03,-0.32679961E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q2ex  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23            
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22            
      t_sp_q2ex  =t_sp_q2ex  
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)*x11        *x43    
     8  +coeff( 17)    *x23    *x42    
      t_sp_q2ex  =t_sp_q2ex  
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)    *x21    *x41*x51
      t_sp_q2ex  =t_sp_q2ex  
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)    *x23    *x43    
     3  +coeff( 30)    *x21    *x42*x53
     4  +coeff( 31)*x12*x21            
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)            *x42*x53
     8  +coeff( 35)*x11*x23    *x41    
      t_sp_q2ex  =t_sp_q2ex  
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)*x11*x22        *x52
     2  +coeff( 38)*x11        *x42*x52
     3  +coeff( 39)*x12*x23            
     4  +coeff( 40)*x11*x22    *x43    
     5  +coeff( 41)*x11*x22    *x41*x52
     6  +coeff( 42)*x11*x23    *x41*x52
     7  +coeff( 43)*x12*x22    *x41*x52
c
      return
      end
      function y_sp_q2ex  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/ -0.1255690E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.86893057E-02, 0.16437545E+00, 0.30927019E-01,-0.86207306E-02,
     +  0.76957135E-02,-0.20065106E-01,-0.16776239E-02,-0.89972671E-02,
     + -0.30937952E-02, 0.70535857E-02,-0.21273898E-01, 0.64245408E-03,
     + -0.10702156E-02,-0.25301829E-02,-0.31363743E-02,-0.31252610E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q2ex  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22    *x41    
      y_sp_q2ex  =y_sp_q2ex  
     9  +coeff(  9)    *x21    *x43    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)            *x44    
     7  +coeff( 16)*x11*x21    *x42    
c
      return
      end
      function p_sp_q2ex  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 48)
      data ncoeff/ 47/
      data avdat/  0.2066354E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13244015E-02, 0.84210705E-03,-0.24872359E-01,-0.36383627E-02,
     +  0.23949414E-02, 0.33792909E-02, 0.81041036E-03,-0.60195540E-03,
     +  0.40643945E-03, 0.23935171E-03,-0.57797477E-03, 0.56365428E-04,
     +  0.15893631E-02,-0.44618431E-03,-0.21043691E-03,-0.29839147E-03,
     + -0.21519675E-03, 0.16445002E-02, 0.46598574E-03, 0.16658405E-03,
     + -0.12427395E-03, 0.19016658E-03,-0.59988146E-04, 0.43061111E-04,
     + -0.49367154E-04,-0.90241818E-04,-0.77729303E-04,-0.14690584E-04,
     +  0.21874244E-03,-0.47865469E-03,-0.60190196E-03,-0.15888859E-03,
     + -0.26677773E-03,-0.48308437E-04, 0.10874084E-02, 0.39810225E-03,
     + -0.13306741E-04,-0.11798259E-03,-0.16810039E-03,-0.11331667E-03,
     +  0.17573079E-03,-0.52693558E-04,-0.84431580E-03, 0.50621643E-03,
     +  0.65397040E-03,-0.14924785E-03, 0.41297948E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q2ex  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22        *x51
      p_sp_q2ex  =p_sp_q2ex  
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)            *x43    
     8  +coeff( 17)            *x41*x52
      p_sp_q2ex  =p_sp_q2ex  
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)*x11                
     7  +coeff( 25)                *x53
     8  +coeff( 26)*x11*x22            
      p_sp_q2ex  =p_sp_q2ex  
     9  +coeff( 27)*x11*x21        *x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)*x11*x21    *x42    
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)    *x22    *x42*x51
     5  +coeff( 32)    *x21    *x43*x51
     6  +coeff( 33)*x11*x23    *x41    
     7  +coeff( 34)*x11*x23    *x41*x51
     8  +coeff( 35)*x11*x23    *x43*x51
      p_sp_q2ex  =p_sp_q2ex  
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x22        *x52
     2  +coeff( 38)            *x42*x52
     3  +coeff( 39)*x11*x21    *x41*x51
     4  +coeff( 40)            *x43*x52
     5  +coeff( 41)*x11*x21    *x43    
     6  +coeff( 42)*x12*x23            
     7  +coeff( 43)    *x23    *x43    
     8  +coeff( 44)    *x22    *x42*x52
      p_sp_q2ex  =p_sp_q2ex  
     9  +coeff( 45)*x12*x22    *x42*x51
     1  +coeff( 46)*x12*x22        *x53
     2  +coeff( 47)*x11*x23    *x41*x53
c
      return
      end
      function l_sp_q2ex  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1817624E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24707757E-02,-0.95761270E-05,-0.39369981E-02,-0.49924357E-02,
     +  0.33065991E-03,-0.36673981E-02,-0.14008803E-02,-0.14072334E-03,
     +  0.10444118E-02, 0.16399984E-03, 0.26587004E-04, 0.89168032E-04,
     + -0.97320750E-04, 0.58400159E-03, 0.73890486E-04, 0.17622103E-03,
     +  0.10218564E-03,-0.19062281E-03, 0.12504532E-02,-0.10627510E-03,
     +  0.26154425E-03,-0.20599984E-03, 0.48870817E-04, 0.64213500E-05,
     + -0.24085153E-03,-0.72834082E-05, 0.31275212E-03,-0.51993622E-04,
     + -0.13448969E-04, 0.71941977E-04, 0.21562724E-03,-0.10778074E-03,
     +  0.19747027E-04,-0.63733582E-03, 0.32696468E-03, 0.12533603E-04,
     +  0.15559726E-04, 0.32175925E-04,-0.66042259E-04,-0.69832211E-04,
     + -0.32156237E-04,-0.33813052E-04,-0.48321563E-04, 0.11051204E-03,
     +  0.18811390E-04,-0.17421224E-04, 0.25017804E-03,-0.11397631E-03,
     +  0.56876914E-04, 0.53799693E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q2ex  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)                *x52
      l_sp_q2ex  =l_sp_q2ex  
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)                *x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11*x21    *x41    
      l_sp_q2ex  =l_sp_q2ex  
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)            *x43*x52
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)*x11                
      l_sp_q2ex  =l_sp_q2ex  
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x21    *x42    
     4  +coeff( 31)    *x22    *x42*x51
     5  +coeff( 32)*x11*x23    *x41    
     6  +coeff( 33)*x11*x21    *x42*x51
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)*x11*x23    *x42*x51
      l_sp_q2ex  =l_sp_q2ex  
     9  +coeff( 36)*x11            *x51
     1  +coeff( 37)    *x21    *x41*x51
     2  +coeff( 38)                *x53
     3  +coeff( 39)    *x22        *x52
     4  +coeff( 40)            *x42*x52
     5  +coeff( 41)*x11*x22        *x51
     6  +coeff( 42)    *x22    *x41*x52
     7  +coeff( 43)    *x22        *x53
     8  +coeff( 44)*x11*x21    *x43    
      l_sp_q2ex  =l_sp_q2ex  
     9  +coeff( 45)*x11        *x41*x53
     1  +coeff( 46)*x12*x23            
     2  +coeff( 47)    *x22    *x42*x52
     3  +coeff( 48)*x11*x22    *x41*x52
     4  +coeff( 49)*x11        *x43*x52
     5  +coeff( 50)*x12*x22    *x41*x51
c
      return
      end
      function x_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5096715E+01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16981035E-02,-0.11521436E+00, 0.31246699E-02, 0.13506039E-01,
     + -0.69051618E-02,-0.46062935E-02, 0.20845260E-02,-0.63098497E-02,
     +  0.66511915E-02,-0.83984444E-02,-0.64626266E-03,-0.10985623E-02,
     +  0.21396268E-02, 0.32996051E-02, 0.19442633E-03,-0.20981077E-03,
     + -0.74405211E-03, 0.71802689E-03,-0.32466976E-03, 0.27429161E-03,
     +  0.70010894E-03,-0.33849804E-02,-0.34230638E-02,-0.87713957E-03,
     + -0.12558028E-03, 0.75290905E-03,-0.65974775E-03, 0.29659248E-03,
     +  0.49159321E-03,-0.21381399E-02,-0.10129706E-01,-0.23533861E-03,
     +  0.31058372E-02,-0.77749009E-03, 0.32054700E-03, 0.15184367E-02,
     +  0.34542356E-03,-0.81409002E-04,-0.35819455E-03,-0.15552107E-02,
     +  0.67867147E-03, 0.15830037E-02, 0.60240773E-03, 0.11853307E-05,
     + -0.25622157E-03, 0.11194149E-03, 0.47052471E-03,-0.32191896E-04,
     + -0.11570938E-02, 0.16058946E-02, 0.35357245E-03, 0.17922436E-02,
     + -0.31546862E-02, 0.24559586E-02,-0.37266003E-03, 0.19673014E-03,
     +  0.77716795E-05, 0.27477107E-03,-0.19845925E-03, 0.58250816E-03,
     +  0.45461868E-03,-0.21850229E-03,-0.15199352E-03, 0.99464774E-03,
     + -0.28384336E-04,-0.62044004E-04,-0.13502338E-03,-0.16014196E-02,
     + -0.10845889E-03,-0.10985750E-02, 0.40152392E-03,-0.15893058E-03,
     +  0.12065711E-03,-0.18516650E-03, 0.65309345E-03, 0.10596196E-02,
     + -0.68055472E-03, 0.69201493E-03, 0.25283778E-03, 0.11728761E-03,
     + -0.10277742E-02,-0.22965875E-03,-0.10315423E-02, 0.50406950E-03,
     +  0.18227343E-02,-0.57265215E-03,-0.29656192E-03, 0.98349305E-03,
     + -0.15570035E-02, 0.14493680E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23            
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22    *x41    
      x_sp_den  =x_sp_den  
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x24    *x41    
     6  +coeff( 15)*x11        *x43    
     7  +coeff( 16)                *x51
     8  +coeff( 17)    *x22            
      x_sp_den  =x_sp_den  
     9  +coeff( 18)            *x42    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)*x11        *x41    
      x_sp_den  =x_sp_den  
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)    *x21    *x41*x51
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)    *x21    *x42*x51
     8  +coeff( 35)*x11*x24            
      x_sp_den  =x_sp_den  
     9  +coeff( 36)*x11*x24    *x41    
     1  +coeff( 37)            *x43    
     2  +coeff( 38)*x12*x21            
     3  +coeff( 39)*x12*x22            
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)*x11*x23    *x41    
     6  +coeff( 42)    *x24    *x42    
     7  +coeff( 43)*x12*x22        *x51
     8  +coeff( 44)*x11        *x42*x52
      x_sp_den  =x_sp_den  
     9  +coeff( 45)*x11*x23    *x42    
     1  +coeff( 46)*x11*x23    *x41*x51
     2  +coeff( 47)*x12*x24            
     3  +coeff( 48)*x12*x22    *x42    
     4  +coeff( 49)*x12*x24        *x51
     5  +coeff( 50)*x12*x23    *x41*x51
     6  +coeff( 51)*x12*x21    *x42*x52
     7  +coeff( 52)*x11*x22    *x43*x52
     8  +coeff( 53)*x12*x23    *x42*x51
      x_sp_den  =x_sp_den  
     9  +coeff( 54)    *x23    *x43*x53
     1  +coeff( 55)    *x23        *x51
     2  +coeff( 56)    *x22    *x41*x51
     3  +coeff( 57)    *x21    *x41*x52
     4  +coeff( 58)    *x24        *x51
     5  +coeff( 59)*x11*x21    *x42    
     6  +coeff( 60)    *x23    *x41*x51
     7  +coeff( 61)    *x22    *x41*x52
     8  +coeff( 62)*x12*x21    *x41    
      x_sp_den  =x_sp_den  
     9  +coeff( 63)*x12*x21        *x51
     1  +coeff( 64)    *x21    *x43*x51
     2  +coeff( 65)            *x42*x53
     3  +coeff( 66)*x12*x21    *x41*x51
     4  +coeff( 67)*x12*x21        *x52
     5  +coeff( 68)*x11*x22    *x43    
     6  +coeff( 69)*x12        *x42*x51
     7  +coeff( 70)    *x24    *x41*x52
     8  +coeff( 71)*x12*x23        *x51
      x_sp_den  =x_sp_den  
     9  +coeff( 72)*x11*x21    *x42*x52
     1  +coeff( 73)*x12*x22        *x52
     2  +coeff( 74)*x11        *x43*x52
     3  +coeff( 75)*x12*x21    *x42*x51
     4  +coeff( 76)*x11*x23    *x42*x51
     5  +coeff( 77)    *x21    *x43*x53
     6  +coeff( 78)*x12*x24    *x41    
     7  +coeff( 79)*x12*x23    *x42    
     8  +coeff( 80)*x11*x21    *x43*x52
      x_sp_den  =x_sp_den  
     9  +coeff( 81)    *x23    *x43*x52
     1  +coeff( 82)*x11*x21    *x42*x53
     2  +coeff( 83)*x12*x22    *x43    
     3  +coeff( 84)    *x22    *x43*x53
     4  +coeff( 85)*x12*x24    *x42    
     5  +coeff( 86)*x12*x24    *x41*x51
     6  +coeff( 87)*x11*x21    *x43*x53
     7  +coeff( 88)*x11*x23    *x43*x52
     8  +coeff( 89)*x12*x23    *x43*x51
      x_sp_den  =x_sp_den  
     9  +coeff( 90)*x12*x23    *x42*x52
c
      return
      end
      function t_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/  0.1297845E+01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.58682117E-03, 0.51097244E-01,-0.19126425E-02, 0.49371142E-02,
     + -0.68945331E-02, 0.30999614E-02, 0.23256394E-02,-0.24112206E-03,
     +  0.61964296E-03,-0.10552508E-02, 0.10667635E-03, 0.71490800E-03,
     + -0.38128598E-02, 0.44102105E-02, 0.30545253E-03, 0.17922612E-03,
     + -0.38797007E-03, 0.46140654E-03, 0.13413153E-02,-0.56323552E-04,
     +  0.23902082E-02, 0.16886522E-02, 0.40229145E-03,-0.35770956E-03,
     + -0.20620700E-03, 0.25328985E-03, 0.16095322E-03,-0.13888617E-03,
     +  0.75153448E-03, 0.56599965E-02, 0.15267635E-03, 0.72946132E-04,
     + -0.16790567E-02,-0.14394557E-03, 0.20110911E-03,-0.83469800E-04,
     +  0.22272936E-03,-0.14148229E-03, 0.76005515E-03,-0.46708214E-04,
     + -0.52594766E-03,-0.39612842E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      t_sp_den  =t_sp_den  
     9  +coeff(  9)                *x51
     1  +coeff( 10)            *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21        *x52
      t_sp_den  =t_sp_den  
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)*x11        *x43    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)*x11*x22    *x42    
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)            *x42*x51
      t_sp_den  =t_sp_den  
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)    *x22        *x51
     5  +coeff( 32)            *x41*x52
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)    *x22    *x41*x51
     8  +coeff( 35)    *x21    *x42*x51
      t_sp_den  =t_sp_den  
     9  +coeff( 36)    *x22        *x52
     1  +coeff( 37)*x11*x21    *x43    
     2  +coeff( 38)*x12*x22    *x41    
     3  +coeff( 39)*x11*x22    *x43    
     4  +coeff( 40)*x11*x23        *x52
     5  +coeff( 41)*x11*x22    *x41*x52
     6  +coeff( 42)*x11*x23    *x41*x52
c
      return
      end
      function y_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.4713811E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.31306199E-02, 0.67911580E-01, 0.16590914E-01,-0.53495886E-02,
     +  0.30774388E-02, 0.14186722E-01, 0.19191395E-02, 0.22967993E-02,
     + -0.11721007E-01,-0.11691705E-02,-0.55565028E-02, 0.31203208E-02,
     + -0.81726387E-02,-0.16448221E-02,-0.67748791E-02, 0.27288531E-03,
     + -0.80801884E-03,-0.15407763E-02,-0.68682706E-03,-0.18809831E-02,
     +  0.28078402E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_den  =y_sp_den  
     9  +coeff(  9)    *x22            
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x22    *x43    
     7  +coeff( 16)    *x21        *x51
     8  +coeff( 17)            *x41*x52
      y_sp_den  =y_sp_den  
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)    *x23    *x41    
c
      return
      end
      function p_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.7511566E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37622619E-02, 0.55136187E-02,-0.81334345E-01,-0.10808811E-01,
     +  0.97114835E-02,-0.21333363E-01,-0.63938368E-03,-0.44747558E-02,
     +  0.13848591E-01, 0.30591616E-02,-0.10145593E-01,-0.41133412E-02,
     +  0.64681313E-03, 0.81394817E-03, 0.41222563E-02, 0.59389290E-02,
     + -0.11977333E-02, 0.94783690E-03,-0.13442024E-02,-0.86945528E-03,
     +  0.11680829E-02,-0.11630149E-02,-0.15372341E-03, 0.37094473E-02,
     + -0.21034994E-02, 0.12040242E-03, 0.21127949E-03,-0.26864879E-03,
     +  0.25877010E-02, 0.40719294E-03,-0.10890982E-03,-0.28777556E-02,
     + -0.50197780E-03, 0.44659691E-03,-0.30730416E-02, 0.21112337E-04,
     + -0.15445873E-03, 0.14236987E-03,-0.18484148E-03,-0.95095653E-04,
     +  0.12353164E-03, 0.22504214E-03,-0.26668675E-03, 0.82171580E-03,
     + -0.34675666E-03,-0.11012919E-02,-0.73330721E-03, 0.10830719E-02,
     +  0.18840343E-02,-0.66871912E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_den  =p_sp_den  
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)            *x41*x52
      p_sp_den  =p_sp_den  
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)    *x23            
     4  +coeff( 22)            *x43    
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)    *x22    *x43    
     7  +coeff( 25)    *x22    *x42*x51
     8  +coeff( 26)*x11            *x51
      p_sp_den  =p_sp_den  
     9  +coeff( 27)    *x21        *x52
     1  +coeff( 28)                *x53
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x23        *x51
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)*x11*x23    *x41    
     7  +coeff( 34)*x11*x21    *x43    
     8  +coeff( 35)    *x23    *x43    
      p_sp_den  =p_sp_den  
     9  +coeff( 36)*x11                
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)*x11*x21        *x51
     4  +coeff( 40)*x11        *x41*x51
     5  +coeff( 41)    *x22        *x52
     6  +coeff( 42)    *x21    *x41*x52
     7  +coeff( 43)            *x42*x52
     8  +coeff( 44)*x11*x21    *x42    
      p_sp_den  =p_sp_den  
     9  +coeff( 45)*x11*x21    *x42*x51
     1  +coeff( 46)*x11*x23    *x42    
     2  +coeff( 47)*x11*x22    *x41*x52
     3  +coeff( 48)    *x22    *x43*x52
     4  +coeff( 49)*x11*x23    *x42*x51
     5  +coeff( 50)*x11*x23    *x41*x52
c
      return
      end
      function l_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.3990354E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.75977365E-02, 0.22288786E+00,-0.68230387E-02,-0.60518216E-02,
     + -0.25638392E-01, 0.13829921E-01,-0.12027704E-01, 0.87228343E-02,
     + -0.63989912E-02, 0.71439235E-02,-0.12918758E-01, 0.15910666E-01,
     + -0.15720094E-02, 0.16801985E-02, 0.70194103E-03, 0.16058641E-02,
     +  0.69957301E-02, 0.21148798E-03, 0.78422166E-02, 0.63787866E-03,
     + -0.12408050E-02,-0.13886839E-02, 0.11214203E-01, 0.25295198E-02,
     +  0.44935520E-03,-0.19019294E-03,-0.16055334E-02,-0.16262985E-02,
     +  0.83087228E-03,-0.59699104E-03,-0.75974601E-03, 0.38696618E-02,
     +  0.18567109E-01,-0.51528029E-03,-0.54747080E-02, 0.72725269E-03,
     +  0.10266728E-02, 0.54688437E-04,-0.71788288E-03,-0.25170757E-02,
     + -0.20024655E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      l_sp_den  =l_sp_den  
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x22    *x42    
      l_sp_den  =l_sp_den  
     9  +coeff( 18)*x11        *x43    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)                *x51
     8  +coeff( 26)                *x52
      l_sp_den  =l_sp_den  
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)            *x43    
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)*x11        *x42    
     4  +coeff( 31)*x11*x23            
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)    *x23    *x43    
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)    *x21    *x43    
      l_sp_den  =l_sp_den  
     9  +coeff( 36)    *x21    *x42*x51
     1  +coeff( 37)*x11*x21    *x43    
     2  +coeff( 38)*x11*x22        *x52
     3  +coeff( 39)*x12*x22    *x41    
     4  +coeff( 40)*x11*x22    *x41*x52
     5  +coeff( 41)*x11*x23    *x41*x52
c
      return
      end
      function x_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3488090E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23303491E-02, 0.31215775E+00,-0.59613464E-02, 0.13342120E+00,
     + -0.11632902E-01, 0.22775318E-01,-0.41492160E-01, 0.42424414E-01,
     +  0.16802428E-01,-0.35101741E-02, 0.13982259E-01,-0.22881454E-01,
     +  0.27753811E-01,-0.39335643E-02,-0.29291681E-02, 0.30644787E-02,
     + -0.68174661E-02,-0.21611566E-02,-0.12778412E-02, 0.82139932E-02,
     +  0.40167966E-02,-0.45936438E-02, 0.11439581E-01,-0.56546152E-03,
     +  0.10221186E-02, 0.31456666E-03, 0.22286361E-02, 0.30419743E-02,
     + -0.13548157E-02, 0.30098842E-01, 0.13489387E-02,-0.94701918E-02,
     + -0.73296425E-03,-0.96037716E-03,-0.22161534E-03,-0.87981287E-03,
     +  0.51340635E-03, 0.16049974E-03,-0.50500402E-03, 0.56331689E-02,
     +  0.91391301E-03,-0.30119850E-02,-0.35380959E-02, 0.35738740E-02,
     + -0.44490434E-02,-0.65174932E-03, 0.36800466E-02,-0.62780717E-03,
     + -0.13461944E-02,-0.14385364E-02, 0.13907839E-02,-0.18487206E-02,
     +  0.24629624E-02,-0.12257314E-01,-0.34452265E-02,-0.10710343E-02,
     +  0.69140331E-02, 0.49929465E-02,-0.17983671E-01, 0.68857255E-04,
     + -0.16176212E-03,-0.59077207E-04, 0.58521499E-03,-0.15023865E-02,
     + -0.22401227E-03, 0.11635725E-03,-0.20782598E-02,-0.98677329E-03,
     +  0.39993340E-03, 0.13302983E-02,-0.48975233E-03, 0.36663408E-03,
     +  0.81434811E-03, 0.60078236E-02,-0.84414461E-03, 0.13537899E-02,
     + -0.23993850E-03, 0.54128761E-02, 0.34988734E-02, 0.77348598E-03,
     +  0.18935562E-02,-0.42870457E-03,-0.26936755E-02,-0.20659904E-02,
     +  0.19428483E-02, 0.15971331E-02, 0.73637447E-03,-0.24400786E-02,
     + -0.24282765E-02, 0.20393042E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x24            
      x_sp_dex    =x_sp_dex    
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)    *x24    *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)            *x41*x51
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)*x11            *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)    *x22        *x51
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x24            
     7  +coeff( 34)            *x43    
     8  +coeff( 35)*x11*x21        *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 36)*x11        *x42    
     1  +coeff( 37)*x12*x21            
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x11        *x43    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)    *x22    *x42*x51
     6  +coeff( 42)*x11*x23    *x41    
     7  +coeff( 43)    *x21    *x43*x51
     8  +coeff( 44)*x11*x22    *x42    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 45)    *x24    *x42    
     1  +coeff( 46)    *x24        *x52
     2  +coeff( 47)    *x23    *x42*x51
     3  +coeff( 48)*x11*x21    *x41*x52
     4  +coeff( 49)    *x23        *x53
     5  +coeff( 50)*x11*x24    *x41    
     6  +coeff( 51)*x11*x23    *x42    
     7  +coeff( 52)*x11*x22    *x42*x51
     8  +coeff( 53)*x11*x22    *x41*x52
      x_sp_dex    =x_sp_dex    
     9  +coeff( 54)*x12*x24    *x41    
     1  +coeff( 55)*x12*x23    *x42    
     2  +coeff( 56)*x12*x22    *x43    
     3  +coeff( 57)*x11*x24    *x42*x51
     4  +coeff( 58)*x12*x24    *x41*x51
     5  +coeff( 59)*x11*x24    *x43*x52
     6  +coeff( 60)*x12                
     7  +coeff( 61)*x11        *x41*x51
     8  +coeff( 62)    *x21    *x41*x52
      x_sp_dex    =x_sp_dex    
     9  +coeff( 63)    *x21        *x53
     1  +coeff( 64)    *x24        *x51
     2  +coeff( 65)            *x41*x53
     3  +coeff( 66)    *x23    *x41*x51
     4  +coeff( 67)    *x22    *x41*x52
     5  +coeff( 68)    *x22        *x53
     6  +coeff( 69)*x12*x21    *x41    
     7  +coeff( 70)*x11*x23        *x51
     8  +coeff( 71)    *x21    *x42*x52
      x_sp_dex    =x_sp_dex    
     9  +coeff( 72)*x11*x21    *x42*x51
     1  +coeff( 73)*x11*x21        *x53
     2  +coeff( 74)*x12*x22    *x41    
     3  +coeff( 75)*x11*x24        *x51
     4  +coeff( 76)    *x21    *x43*x52
     5  +coeff( 77)*x12        *x43    
     6  +coeff( 78)*x11*x22    *x43    
     7  +coeff( 79)    *x24    *x41*x52
     8  +coeff( 80)*x11*x22        *x53
      x_sp_dex    =x_sp_dex    
     9  +coeff( 81)    *x24        *x53
     1  +coeff( 82)*x12*x23        *x51
     2  +coeff( 83)    *x23    *x41*x53
     3  +coeff( 84)*x12*x22    *x41*x51
     4  +coeff( 85)    *x22    *x43*x52
     5  +coeff( 86)*x11*x23    *x43    
     6  +coeff( 87)*x12*x21    *x42*x51
     7  +coeff( 88)*x11*x23    *x42*x51
     8  +coeff( 89)*x11*x23        *x53
      x_sp_dex    =x_sp_dex    
     9  +coeff( 90)    *x21    *x43*x53
c
      return
      end
      function t_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 38)
      data ncoeff/ 37/
      data avdat/  0.5345973E+00/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11696648E-02,-0.71709476E-01, 0.13199901E-02, 0.22726601E-01,
     +  0.19446103E-02,-0.48324461E-02, 0.87519428E-02,-0.40068193E-02,
     + -0.37231103E-02,-0.13033279E-02,-0.10911556E-02,-0.53762961E-02,
     +  0.37467571E-02, 0.54006244E-03,-0.54312369E-03,-0.16962084E-02,
     + -0.20340628E-03, 0.43654855E-03,-0.20740539E-03,-0.33188274E-03,
     +  0.56379591E-03,-0.81033073E-03,-0.25490138E-02,-0.34515184E-03,
     +  0.26866078E-03, 0.23300794E-03,-0.23079340E-02,-0.63299406E-02,
     + -0.64543994E-04, 0.36697075E-03, 0.85108746E-04,-0.76852979E-04,
     +  0.19325701E-02,-0.16491445E-03,-0.67384937E-03, 0.28186344E-03,
     +  0.56043704E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sp_dex    =t_sp_dex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)            *x41*x51
      t_sp_dex    =t_sp_dex    
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)    *x22        *x52
     8  +coeff( 26)*x11*x23            
      t_sp_dex    =t_sp_dex    
     9  +coeff( 27)    *x22    *x43    
     1  +coeff( 28)    *x23    *x43    
     2  +coeff( 29)*x11*x21            
     3  +coeff( 30)            *x43    
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x12*x21            
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)            *x42*x52
     8  +coeff( 35)*x11*x22    *x42    
      t_sp_dex    =t_sp_dex    
     9  +coeff( 36)*x11        *x42*x52
     1  +coeff( 37)*x11*x23    *x41*x52
c
      return
      end
      function y_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/  0.6287763E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25956267E-02,-0.44369400E-01, 0.72167832E-02, 0.15597047E-02,
     +  0.32638742E-02, 0.46117514E-01, 0.83063571E-02,-0.30842822E-01,
     + -0.83022602E-02,-0.40021658E-03,-0.19295836E-02, 0.94285433E-03,
     + -0.33900915E-02,-0.96895837E-03,-0.47304649E-02,-0.22008706E-01,
     + -0.91860779E-02, 0.68203835E-02, 0.28554781E-03, 0.39431741E-02,
     + -0.62550948E-03,-0.44132951E-02,-0.32736869E-02,-0.31998856E-02,
     + -0.47079223E-03, 0.15413454E-02,-0.45798030E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_dex    =y_sp_dex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)    *x22            
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x22        *x51
      y_sp_dex    =y_sp_dex    
     9  +coeff( 18)    *x23            
     1  +coeff( 19)            *x45    
     2  +coeff( 20)    *x21    *x44    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)    *x23        *x51
      y_sp_dex    =y_sp_dex    
     9  +coeff( 27)    *x22    *x42*x51
c
      return
      end
      function p_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 49)
      data ncoeff/ 48/
      data avdat/  0.2396113E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11878654E-02, 0.89041301E-03,-0.22729566E-01,-0.85213181E-03,
     +  0.11720856E-02,-0.64768973E-02, 0.29571913E-03,-0.20061631E-02,
     +  0.90583907E-02, 0.18317975E-02,-0.34000375E-02,-0.20597281E-02,
     + -0.11162928E-02, 0.22175914E-03, 0.10175956E-02, 0.85179030E-03,
     + -0.99242746E-03,-0.61458297E-03,-0.22095472E-03,-0.20744152E-03,
     +  0.62877854E-03,-0.97168417E-03,-0.13402371E-03, 0.10771247E-03,
     + -0.14105585E-03, 0.56256034E-03, 0.66345790E-04, 0.17504583E-04,
     +  0.66004670E-03,-0.19007255E-03, 0.16315168E-03,-0.52275584E-04,
     +  0.30153071E-04, 0.70721994E-03, 0.16373827E-03,-0.72780233E-04,
     + -0.78290715E-04,-0.53426083E-04, 0.40413226E-04,-0.78805089E-04,
     + -0.56880841E-03, 0.31686228E-03, 0.12117721E-03,-0.68188226E-03,
     + -0.65592874E-03,-0.17403055E-03,-0.12680501E-03, 0.16930210E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_dex    =p_sp_dex    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22    *x41*x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)    *x23        *x51
     4  +coeff( 22)    *x22    *x42*x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)    *x22    *x42    
      p_sp_dex    =p_sp_dex    
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x21    *x41*x52
     4  +coeff( 31)    *x23    *x41*x51
     5  +coeff( 32)            *x42*x51
     6  +coeff( 33)*x11        *x42    
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)            *x43*x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 36)            *x41*x53
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)*x11*x22        *x51
     4  +coeff( 40)*x11*x21    *x41*x51
     5  +coeff( 41)    *x23    *x42    
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x23        *x51
     8  +coeff( 44)    *x23    *x43    
      p_sp_dex    =p_sp_dex    
     9  +coeff( 45)    *x22    *x43*x51
     1  +coeff( 46)*x11*x22    *x43    
     2  +coeff( 47)*x12*x22    *x42    
     3  +coeff( 48)*x12*x22    *x41*x51
c
      return
      end
      function l_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.1278484E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.86215278E-02,-0.44609931E+00, 0.42620976E-02,-0.98771535E-01,
     +  0.14877189E-01,-0.39648384E-01, 0.58067765E-01,-0.42922355E-01,
     + -0.24095081E-01,-0.73344144E-02, 0.30813307E-01,-0.37948854E-01,
     + -0.13553023E-02,-0.17888561E-02, 0.17310387E-02,-0.33678804E-02,
     + -0.97059766E-02,-0.94245697E-04,-0.17350364E-01,-0.95062685E-03,
     + -0.22162714E-02, 0.33748795E-02, 0.22068464E-02,-0.56869015E-02,
     +  0.31397452E-02,-0.21350579E-02, 0.12002179E-02,-0.27093329E-02,
     + -0.71731564E-02,-0.12791391E-01,-0.42095818E-01, 0.18238588E-02,
     +  0.80831075E-03, 0.41730647E-03,-0.51050028E-03, 0.12310297E-01,
     + -0.14644387E-02, 0.56584319E-03,-0.20588655E-02, 0.20223309E-02,
     +  0.38399748E-02, 0.57051452E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_sp_dex    =l_sp_dex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x22    *x42    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 18)*x11        *x43    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)*x11*x21    *x41    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)    *x22    *x43    
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)            *x43    
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)*x12*x21            
      l_sp_dex    =l_sp_dex    
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)*x11*x21    *x42    
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x22    *x41*x52
     6  +coeff( 42)*x11*x23    *x42*x52
c
      return
      end
      function x_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3959992E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13935241E-02, 0.19682968E+00,-0.37743514E-02, 0.13865240E+00,
     + -0.81181712E-02, 0.24504457E-01,-0.27642675E-01, 0.33161711E-01,
     +  0.12318333E-01,-0.59591304E-02, 0.61197532E-02,-0.15822114E-01,
     +  0.19149899E-01,-0.30839450E-02, 0.37720893E-02,-0.28248799E-02,
     +  0.21135889E-02,-0.39799758E-02,-0.71442273E-03,-0.16017917E-02,
     +  0.26813538E-02, 0.48261047E-02, 0.78265378E-02, 0.48735049E-02,
     +  0.13488324E-02,-0.83945168E-03, 0.18647691E-01, 0.20207610E-03,
     +  0.11280699E-03,-0.69361209E-03,-0.67468570E-03, 0.39186308E-03,
     +  0.41671437E-02,-0.56815175E-02,-0.53446845E-04,-0.76153426E-03,
     + -0.26323011E-02,-0.51930879E-03,-0.14620865E-03, 0.29788172E-03,
     +  0.31681440E-03, 0.37306716E-03, 0.58369775E-03, 0.11258656E-02,
     + -0.21386335E-02,-0.15387961E-02,-0.19685966E-02, 0.17081758E-02,
     +  0.59422418E-04, 0.21843307E-02,-0.39714198E-02, 0.79105003E-03,
     +  0.41701025E-02,-0.15744879E-02,-0.12157859E-01,-0.35909098E-02,
     + -0.19542812E-02, 0.23458926E-02, 0.35940851E-02, 0.11416777E-01,
     +  0.61827619E-03,-0.17465638E-01,-0.71837567E-04,-0.15094574E-02,
     + -0.12891568E-03,-0.17193765E-03, 0.36426386E-03, 0.40341899E-03,
     + -0.65809063E-03, 0.60498207E-02,-0.92873978E-03, 0.44719168E-03,
     + -0.38191627E-03,-0.19272220E-03, 0.95619919E-03, 0.98510436E-03,
     +  0.49794052E-03,-0.12321502E-02,-0.52049011E-03, 0.99035256E-04,
     +  0.10276294E-02,-0.22743222E-03, 0.13904993E-02, 0.22304074E-02,
     + -0.14192832E-03, 0.42134495E-02, 0.56779460E-03, 0.51600654E-02,
     + -0.83575380E-03,-0.50000111E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)*x11*x22            
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 18)    *x24            
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x23        *x51
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)    *x22    *x43    
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x23            
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 27)    *x23    *x43    
     1  +coeff( 28)*x11*x21            
     2  +coeff( 29)*x11            *x51
     3  +coeff( 30)    *x21        *x52
     4  +coeff( 31)            *x43    
     5  +coeff( 32)                *x53
     6  +coeff( 33)    *x22    *x42    
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)    *x24    *x41    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 36)*x11*x24            
     1  +coeff( 37)    *x24    *x42    
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)*x11        *x41*x51
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)    *x21    *x42*x51
     6  +coeff( 42)*x11*x22        *x51
     7  +coeff( 43)    *x23        *x52
     8  +coeff( 44)    *x22    *x42*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 45)    *x22        *x53
     1  +coeff( 46)*x11*x23    *x41    
     2  +coeff( 47)    *x21    *x43*x51
     3  +coeff( 48)*x11*x22    *x42    
     4  +coeff( 49)*x11*x22    *x41*x51
     5  +coeff( 50)    *x23    *x42*x51
     6  +coeff( 51)*x11*x24    *x41    
     7  +coeff( 52)*x11*x23    *x42    
     8  +coeff( 53)    *x24        *x53
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 54)*x11*x23    *x41*x52
     1  +coeff( 55)*x12*x24    *x41    
     2  +coeff( 56)*x12*x23    *x42    
     3  +coeff( 57)*x12*x22    *x43    
     4  +coeff( 58)*x11*x24    *x41*x52
     5  +coeff( 59)*x11*x24    *x43*x51
     6  +coeff( 60)*x12*x24    *x41*x52
     7  +coeff( 61)*x12*x22    *x43*x52
     8  +coeff( 62)*x11*x24    *x43*x52
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 63)    *x21    *x41*x52
     1  +coeff( 64)    *x24        *x51
     2  +coeff( 65)            *x41*x53
     3  +coeff( 66)*x11        *x42*x51
     4  +coeff( 67)*x12*x21    *x41    
     5  +coeff( 68)*x12*x21        *x51
     6  +coeff( 69)    *x21    *x42*x52
     7  +coeff( 70)*x12*x22    *x41    
     8  +coeff( 71)*x11*x24        *x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 72)*x12*x21    *x42    
     1  +coeff( 73)*x11*x23    *x41*x51
     2  +coeff( 74)*x12        *x43    
     3  +coeff( 75)*x11*x22    *x43    
     4  +coeff( 76)    *x24    *x41*x52
     5  +coeff( 77)*x11*x22        *x53
     6  +coeff( 78)*x12*x23        *x51
     7  +coeff( 79)    *x23    *x42*x52
     8  +coeff( 80)*x11*x21    *x41*x53
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 81)*x12*x22    *x42    
     1  +coeff( 82)*x12*x22        *x52
     2  +coeff( 83)*x11*x23    *x43    
     3  +coeff( 84)*x11*x22    *x42*x52
     4  +coeff( 85)    *x24    *x42*x52
     5  +coeff( 86)    *x23    *x43*x52
     6  +coeff( 87)*x11*x21    *x42*x53
     7  +coeff( 88)*x11*x24    *x43    
     8  +coeff( 89)*x12*x22    *x42*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 90)*x12*x22    *x41*x52
c
      return
      end
      function t_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/ -0.2856553E-04/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55983564E-03,-0.72211131E-01, 0.11995035E-02, 0.22753982E-01,
     +  0.19327936E-02,-0.59924386E-02, 0.86752791E-02,-0.49293898E-02,
     + -0.34907842E-02,-0.14117664E-02,-0.12714704E-02, 0.43269387E-02,
     + -0.57709273E-02,-0.49943710E-03, 0.22828353E-05, 0.43803218E-03,
     + -0.16056329E-03,-0.54703763E-03, 0.66197134E-03,-0.90148638E-03,
     + -0.24083273E-02,-0.11517211E-02,-0.33489934E-04,-0.18689894E-04,
     +  0.45130774E-03,-0.23495934E-03,-0.29789531E-03, 0.18318095E-03,
     + -0.61640644E-03, 0.33100194E-03,-0.87184372E-03,-0.20052528E-02,
     +  0.25075296E-03,-0.60431822E-02, 0.27862965E-03, 0.20311866E-03,
     +  0.90977454E-04,-0.56905057E-04, 0.18214406E-02,-0.22003730E-03,
     + -0.18129279E-03, 0.16339077E-03,-0.41342969E-03, 0.73803263E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)*x11        *x43    
     7  +coeff( 16)            *x42    
     8  +coeff( 17)*x11            *x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)*x11*x22    *x42    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)            *x42*x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x22    *x43    
     6  +coeff( 33)    *x21    *x43*x51
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)            *x43    
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 36)    *x21    *x41*x51
     1  +coeff( 37)                *x53
     2  +coeff( 38)*x12*x21            
     3  +coeff( 39)    *x21    *x43    
     4  +coeff( 40)    *x21    *x42*x51
     5  +coeff( 41)    *x22        *x52
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)    *x23        *x52
     8  +coeff( 44)*x11*x22    *x42*x52
c
      return
      end
      function y_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/  0.8666045E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.38104616E-02,-0.66021673E-01, 0.66568307E-02, 0.19248382E-02,
     +  0.32329806E-02, 0.57547923E-01, 0.10354813E-01,-0.34951776E-01,
     + -0.10331575E-01,-0.16114744E-02,-0.43467153E-02,-0.60850051E-02,
     + -0.25498452E-01,-0.11806278E-01, 0.78326426E-02,-0.55625066E-02,
     + -0.39954046E-02,-0.76137303E-03, 0.10501109E-02, 0.35912378E-02,
     + -0.13028379E-02,-0.64889941E-03,-0.85165503E-03,-0.62517944E-03,
     +  0.17210176E-02,-0.47593275E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_q3en   =y_sp_q3en   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)            *x41*x52
     3  +coeff( 12)    *x21    *x41*x51
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22    *x41*x51
     8  +coeff( 17)            *x43    
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)                *x53
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x22    *x42*x51
c
      return
      end
      function p_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 48)
      data ncoeff/ 47/
      data avdat/  0.2475173E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10721337E-02, 0.10074286E-02,-0.21946477E-01, 0.10501429E-02,
     + -0.82416320E-02, 0.42930289E-03,-0.21151982E-02, 0.10229197E-01,
     +  0.20579838E-02,-0.45746421E-02,-0.23271940E-02,-0.36485129E-03,
     +  0.26858944E-03, 0.13222403E-02, 0.10135757E-02,-0.95118949E-03,
     + -0.19464494E-02,-0.65014733E-03, 0.14352775E-03, 0.51053468E-03,
     + -0.10088186E-02,-0.80337457E-04,-0.18541249E-03,-0.14245178E-03,
     +  0.11792550E-03,-0.15610707E-03, 0.28403458E-03, 0.69349841E-03,
     + -0.13158901E-03, 0.13825515E-04,-0.13309841E-03,-0.58232865E-03,
     + -0.13622182E-04,-0.52683885E-04, 0.37606525E-04,-0.82187653E-04,
     + -0.25493271E-04,-0.80674377E-04, 0.51009421E-04,-0.85659776E-04,
     +  0.39056331E-03,-0.78954741E-04, 0.14820675E-03,-0.43617978E-03,
     + -0.27765491E-03, 0.21226614E-03, 0.16606793E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x22        *x51
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x22    *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x23        *x51
     3  +coeff( 21)    *x22    *x42*x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)                *x53
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x21        *x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)    *x21    *x41*x52
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)            *x41*x53
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)*x11                
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)*x11        *x42    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 36)    *x22        *x52
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)*x11*x22        *x51
     4  +coeff( 40)*x11*x21    *x41*x51
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)    *x23        *x52
     7  +coeff( 43)*x11*x23        *x51
     8  +coeff( 44)    *x23    *x43    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 45)    *x23    *x41*x52
     1  +coeff( 46)    *x22    *x41*x53
     2  +coeff( 47)*x12*x22    *x41*x51
c
      return
      end
      function l_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.6459657E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41118697E-02,-0.28979233E+00,-0.32086760E-01, 0.90512084E-02,
     + -0.33979975E-01, 0.37512898E-01,-0.19616937E-01,-0.15785994E-01,
     + -0.37413556E-02,-0.21984349E-02, 0.18696314E-01,-0.24671704E-01,
     +  0.17987404E-02,-0.36226169E-02,-0.24795667E-02,-0.20488033E-02,
     + -0.47725909E-04,-0.77519094E-03, 0.11924255E-02, 0.17765121E-02,
     + -0.47917049E-02,-0.10521615E-01,-0.36958538E-02,-0.40869132E-03,
     +  0.16264539E-04, 0.19454315E-02,-0.13125715E-02, 0.77818893E-03,
     + -0.20696942E-02,-0.44281445E-02,-0.70291469E-02,-0.26801929E-01,
     +  0.99683041E-03, 0.11621285E-02, 0.44440458E-03,-0.25000836E-03,
     +  0.79625100E-02,-0.11254795E-02, 0.97133929E-03,-0.11205580E-02,
     +  0.19284280E-02, 0.21503777E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      l_sp_q3en   =l_sp_q3en   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x41    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11        *x43    
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)*x11*x22    *x42    
     6  +coeff( 24)                *x52
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)*x11        *x41    
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)    *x22    *x43    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)            *x43    
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)            *x41*x52
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 36)*x12*x21            
     1  +coeff( 37)    *x21    *x43    
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)*x11*x22    *x41*x52
     6  +coeff( 42)*x11*x23    *x41*x52
c
      return
      end
      function x_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4776838E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27855753E-02, 0.10534061E+00,-0.21628614E-02, 0.19110584E+00,
     + -0.60017030E-02, 0.17307634E-01,-0.16949685E-01, 0.28735202E-01,
     + -0.97267823E-02, 0.78054788E-02,-0.10862025E-01, 0.12079732E-01,
     + -0.26964755E-02, 0.50067399E-02, 0.37170381E-02,-0.31380621E-02,
     + -0.38278964E-02,-0.86426217E-03, 0.12075955E-02, 0.45390064E-02,
     + -0.50155970E-04,-0.23659965E-03, 0.65913098E-03, 0.13500294E-02,
     +  0.26313332E-02,-0.11318069E-02, 0.19785394E-02,-0.98222774E-03,
     + -0.49584516E-03, 0.82253345E-03,-0.32490387E-03,-0.16665293E-02,
     +  0.31870331E-02, 0.80650690E-03, 0.11285966E-01,-0.75379212E-04,
     + -0.18129972E-03, 0.18421540E-03,-0.33488940E-02, 0.97532851E-04,
     + -0.24580027E-02, 0.24924253E-03, 0.32747857E-03,-0.33474172E-03,
     + -0.71174750E-03, 0.17087056E-02,-0.21804508E-02, 0.28179178E-02,
     +  0.37151456E-02,-0.51682588E-03,-0.10974534E-03, 0.99468532E-04,
     + -0.28398781E-03,-0.98809185E-04,-0.19549560E-03, 0.22423337E-02,
     + -0.13751502E-03, 0.72926696E-03,-0.15794205E-02,-0.13419575E-02,
     + -0.73183427E-03,-0.87072415E-03, 0.14419514E-02, 0.32969078E-03,
     +  0.11522259E-02,-0.67548989E-03,-0.74919447E-03,-0.17007480E-02,
     + -0.81613875E-03,-0.12531925E-02, 0.29431358E-02, 0.23347524E-02,
     + -0.26855038E-02, 0.65493741E-03, 0.28474173E-02,-0.15254075E-02,
     +  0.25830902E-02, 0.11332198E-03,-0.10743142E-04,-0.44843258E-04,
     +  0.34069421E-04, 0.16355301E-03,-0.16015294E-03, 0.86135558E-04,
     + -0.32805986E-03, 0.53386355E-03, 0.59665676E-04, 0.20148318E-03,
     + -0.81156504E-04, 0.10572888E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x24            
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)*x11        *x43    
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)                *x53
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x24    *x41    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 27)*x11*x22    *x42    
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)            *x42*x51
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)    *x22        *x52
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x23    *x43    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 36)*x11            *x51
     1  +coeff( 37)    *x22    *x41*x51
     2  +coeff( 38)*x12*x21            
     3  +coeff( 39)    *x21    *x43    
     4  +coeff( 40)*x11*x22        *x51
     5  +coeff( 41)    *x24        *x51
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)*x11*x21    *x42    
     8  +coeff( 44)*x11*x24            
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 45)*x11*x23    *x41    
     1  +coeff( 46)    *x23    *x42*x51
     2  +coeff( 47)*x11*x24    *x41    
     3  +coeff( 48)    *x24        *x53
     4  +coeff( 49)*x11*x24    *x42*x53
     5  +coeff( 50)            *x43    
     6  +coeff( 51)*x11        *x41*x51
     7  +coeff( 52)*x11            *x52
     8  +coeff( 53)*x11*x23            
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 54)    *x21    *x41*x52
     1  +coeff( 55)    *x21        *x53
     2  +coeff( 56)    *x22    *x43    
     3  +coeff( 57)*x11        *x42*x51
     4  +coeff( 58)    *x22    *x42*x51
     5  +coeff( 59)    *x22        *x53
     6  +coeff( 60)    *x21    *x43*x51
     7  +coeff( 61)    *x21    *x42*x52
     8  +coeff( 62)    *x24    *x42    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 63)    *x24        *x52
     1  +coeff( 64)    *x23    *x41*x52
     2  +coeff( 65)    *x22    *x42*x52
     3  +coeff( 66)*x12*x24            
     4  +coeff( 67)*x11*x23    *x41*x52
     5  +coeff( 68)*x12*x24    *x41    
     6  +coeff( 69)*x11*x22    *x42*x52
     7  +coeff( 70)*x12*x23    *x42    
     8  +coeff( 71)    *x23    *x43*x52
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 72)*x12*x22    *x43    
     1  +coeff( 73)*x12*x24    *x42    
     2  +coeff( 74)*x11*x21    *x43*x53
     3  +coeff( 75)*x12*x24    *x41*x52
     4  +coeff( 76)*x12*x22    *x43*x52
     5  +coeff( 77)*x12*x23    *x43*x53
     6  +coeff( 78)*x11*x21            
     7  +coeff( 79)*x12                
     8  +coeff( 80)            *x41*x52
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 81)*x11*x21        *x51
     1  +coeff( 82)    *x21    *x42*x51
     2  +coeff( 83)*x11*x21    *x41*x51
     3  +coeff( 84)*x11*x21        *x52
     4  +coeff( 85)    *x23        *x52
     5  +coeff( 86)*x12*x22            
     6  +coeff( 87)*x11        *x41*x52
     7  +coeff( 88)*x12*x21    *x41    
     8  +coeff( 89)*x12*x21        *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 90)    *x24    *x41*x51
c
      return
      end
      function t_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/ -0.1741489E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18204733E-02,-0.40451668E-01, 0.66576287E-03, 0.64082094E-01,
     + -0.54657315E-02,-0.20862045E-02, 0.42654853E-02, 0.38657393E-03,
     + -0.14519593E-02, 0.55089436E-03,-0.29798574E-02,-0.24020491E-03,
     + -0.17410374E-03,-0.53411414E-03, 0.13721932E-02, 0.49520587E-03,
     +  0.24310186E-04, 0.20917975E-03,-0.86824410E-04,-0.32690156E-03,
     + -0.34254667E-03,-0.30477485E-03,-0.46638827E-03,-0.15524093E-02,
     + -0.22617336E-03,-0.12367812E-03,-0.84602209E-04,-0.21526971E-03,
     +  0.25969456E-03, 0.22349016E-03,-0.63099782E-03,-0.88958815E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11                
      t_sp_q3m    =t_sp_q3m    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)                *x53
     8  +coeff( 17)*x11*x21            
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)*x11*x21    *x41    
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)            *x42*x52
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)    *x22    *x43    
     5  +coeff( 32)    *x23    *x43    
c
      return
      end
      function y_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.1116186E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48525026E-02,-0.89280561E-01, 0.48446818E-02, 0.31660297E-02,
     +  0.35682882E-02, 0.64916521E-01, 0.12025579E-01,-0.42261567E-01,
     + -0.12091410E-01,-0.39790593E-04,-0.69742249E-02,-0.28307153E-01,
     + -0.14340403E-01, 0.86128609E-02,-0.21780176E-03,-0.94032753E-02,
     +  0.12769973E-02, 0.45393612E-02,-0.24886138E-02,-0.41886191E-02,
     + -0.38191043E-02,-0.57415460E-03,-0.10096730E-02,-0.75764128E-03,
     +  0.23102763E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_q3m    =y_sp_q3m    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21    *x41*x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x23            
     6  +coeff( 15)            *x43*x52
     7  +coeff( 16)    *x22    *x41*x51
     8  +coeff( 17)*x11        *x41    
      y_sp_q3m    =y_sp_q3m    
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)            *x42*x53
     2  +coeff( 20)            *x43    
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x23        *x51
c
      return
      end
      function p_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1248721E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12785489E-03, 0.31157941E-03,-0.97713876E-03,-0.16700301E-02,
     +  0.10850773E-02, 0.23834310E-02,-0.50336419E-03, 0.86497399E-03,
     + -0.68775476E-02,-0.80906111E-03,-0.92006259E-03, 0.26457398E-02,
     +  0.33855194E-03, 0.95952186E-03, 0.17280924E-02, 0.35816224E-03,
     +  0.18366260E-03, 0.13517165E-03, 0.20842091E-03, 0.66056149E-03,
     + -0.50081850E-04, 0.37990878E-04, 0.38063063E-03, 0.99575860E-04,
     + -0.28967013E-03,-0.31067134E-03,-0.36439401E-03,-0.11885707E-03,
     +  0.44476395E-03, 0.65696899E-04,-0.84541876E-04, 0.22643365E-03,
     + -0.36839873E-03, 0.52638643E-03,-0.15891124E-03,-0.50327484E-03,
     +  0.65409089E-03,-0.94929337E-05,-0.11022313E-03,-0.25924581E-04,
     + -0.12828478E-03, 0.61384108E-05, 0.18175309E-03,-0.84503336E-04,
     + -0.21444050E-03,-0.18898334E-03,-0.77861405E-04,-0.64249019E-04,
     +  0.96551987E-04, 0.22652777E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_q3m    =p_sp_q3m    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x23    *x41*x51
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x22        *x52
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 27)    *x21    *x41*x52
     1  +coeff( 28)            *x41*x53
     2  +coeff( 29)*x11*x23    *x41*x51
     3  +coeff( 30)*x11            *x52
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)*x11*x21    *x42    
     6  +coeff( 33)    *x23    *x42    
     7  +coeff( 34)    *x22    *x43    
     8  +coeff( 35)*x11*x22        *x52
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 36)    *x22    *x41*x53
     1  +coeff( 37)    *x22    *x42*x53
     2  +coeff( 38)*x11            *x51
     3  +coeff( 39)    *x21        *x52
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)            *x43*x51
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)    *x23        *x52
     8  +coeff( 44)            *x42*x53
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 45)*x11*x23        *x51
     1  +coeff( 46)*x11*x23    *x42    
     2  +coeff( 47)*x11*x23        *x52
     3  +coeff( 48)*x11*x22        *x53
     4  +coeff( 49)*x12*x21    *x43    
     5  +coeff( 50)*x12*x22    *x42*x51
c
      return
      end
      function l_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.7611644E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.52851294E-02,-0.28980955E+00,-0.32017931E-01, 0.91627752E-02,
     + -0.37113547E-01, 0.37602626E-01,-0.16731950E-01,-0.16205970E-01,
     + -0.39393990E-02,-0.22734436E-02, 0.18524259E-01, 0.99054212E-03,
     + -0.24857257E-01,-0.14775789E-02,-0.32907589E-02,-0.25374941E-02,
     + -0.22392417E-02,-0.48901713E-02,-0.99959504E-02, 0.18878675E-02,
     +  0.19012307E-02,-0.85854513E-03, 0.13802093E-02, 0.17535090E-02,
     + -0.40007504E-02, 0.16964159E-03,-0.12616679E-02,-0.64579765E-02,
     + -0.27350426E-01, 0.48559357E-03,-0.24851409E-03, 0.83623677E-02,
     + -0.15100685E-02, 0.89976337E-03, 0.73788437E-03,-0.11760864E-02,
     + -0.23080695E-02, 0.10134686E-02,-0.63496595E-02, 0.76566064E-02,
     +  0.22541501E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      l_sp_q3m    =l_sp_q3m    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x22            
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)            *x41    
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)*x11*x22    *x41    
     8  +coeff( 26)*x11*x21            
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)    *x23    *x43    
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)*x12*x21            
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)*x11*x23            
     8  +coeff( 35)    *x23    *x41*x51
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 36)    *x23        *x52
     1  +coeff( 37)*x11*x22    *x42    
     2  +coeff( 38)*x11        *x42*x52
     3  +coeff( 39)    *x23    *x42*x51
     4  +coeff( 40)    *x23    *x43*x51
     5  +coeff( 41)*x11*x23    *x41*x52
c
      return
      end
      function x_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.6104513E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15033337E-02, 0.69587328E-01,-0.17860141E-02, 0.32967678E+00,
     + -0.69292220E-02,-0.15321297E-01, 0.35004091E-01,-0.22200098E-01,
     + -0.53037945E-02, 0.19017056E-01, 0.75820182E-02,-0.39148391E-02,
     +  0.51983879E-02,-0.12029795E-01, 0.24975119E-02, 0.11442260E-01,
     + -0.15044493E-02, 0.36646000E-02,-0.47968421E-02,-0.91257878E-03,
     + -0.82785141E-03, 0.88061747E-03, 0.71689987E-03, 0.31088961E-02,
     +  0.45857094E-02, 0.21981404E-03,-0.54288562E-03,-0.95614826E-03,
     +  0.90360513E-03, 0.30587190E-02,-0.10197343E-02,-0.22336643E-02,
     +  0.88724296E-03, 0.29241936E-02, 0.68133208E-03, 0.99381097E-02,
     + -0.28438136E-03,-0.22265068E-02,-0.24492533E-02, 0.14213292E-02,
     + -0.98968926E-03, 0.28307054E-02, 0.12365510E-02, 0.10110728E-02,
     + -0.20910928E-02,-0.30166714E-02,-0.23447209E-03,-0.14408812E-03,
     +  0.10638648E-03, 0.35656671E-03, 0.31507816E-03, 0.24303785E-02,
     + -0.92370348E-03,-0.14307295E-02, 0.18090651E-02,-0.47361961E-03,
     +  0.19921486E-03, 0.57054325E-02,-0.18692749E-02,-0.40194974E-02,
     +  0.12287115E-02, 0.12843208E-03, 0.28639535E-04, 0.21926679E-03,
     + -0.15884852E-03, 0.19166624E-03, 0.19744407E-03, 0.88059918E-04,
     + -0.23840701E-02,-0.74783171E-03,-0.31277575E-03, 0.12555977E-03,
     +  0.15217790E-02, 0.16973107E-02,-0.17535689E-02,-0.11096186E-03,
     + -0.15126124E-02, 0.32728657E-03, 0.25338721E-02,-0.19615968E-02,
     +  0.23403121E-02,-0.20511544E-02, 0.14127219E-02, 0.97598246E-03,
     + -0.39628558E-02, 0.33924691E-03,-0.22013103E-03,-0.92138525E-03,
     +  0.80920878E-03, 0.63363486E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff(  9)    *x24            
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)            *x41*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)*x11*x21            
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 27)            *x43    
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)    *x24    *x41    
     5  +coeff( 32)    *x24        *x51
     6  +coeff( 33)            *x42*x52
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x24    *x41*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 36)    *x23    *x43    
     1  +coeff( 37)*x11            *x51
     2  +coeff( 38)    *x22        *x52
     3  +coeff( 39)    *x21    *x43    
     4  +coeff( 40)*x11*x22        *x51
     5  +coeff( 41)    *x23        *x52
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x22    *x42    
     8  +coeff( 44)    *x24        *x52
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 45)*x11*x24    *x41    
     1  +coeff( 46)*x11*x23    *x41*x52
     2  +coeff( 47)*x11        *x42    
     3  +coeff( 48)*x11        *x41*x51
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x21    *x42*x51
     6  +coeff( 51)*x11*x21    *x42    
     7  +coeff( 52)    *x22    *x42*x51
     8  +coeff( 53)    *x21    *x43*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 54)    *x24    *x42    
     1  +coeff( 55)    *x23    *x41*x52
     2  +coeff( 56)*x12*x24            
     3  +coeff( 57)*x11*x22    *x41*x52
     4  +coeff( 58)    *x24    *x41*x52
     5  +coeff( 59)*x11*x23    *x42*x51
     6  +coeff( 60)*x12*x24    *x41    
     7  +coeff( 61)    *x23    *x42*x53
     8  +coeff( 62)*x12*x22    *x43    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 63)            *x41*x52
     1  +coeff( 64)*x11*x21        *x51
     2  +coeff( 65)    *x21    *x41*x52
     3  +coeff( 66)            *x41*x53
     4  +coeff( 67)*x12*x22            
     5  +coeff( 68)*x11        *x42*x51
     6  +coeff( 69)    *x22    *x41*x52
     7  +coeff( 70)*x11*x23    *x41    
     8  +coeff( 71)*x11*x23        *x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 72)            *x43*x52
     1  +coeff( 73)    *x23    *x42*x51
     2  +coeff( 74)*x12*x22    *x41    
     3  +coeff( 75)*x11*x24        *x51
     4  +coeff( 76)*x12        *x43    
     5  +coeff( 77)*x11*x22    *x42*x51
     6  +coeff( 78)*x11*x21    *x43*x51
     7  +coeff( 79)*x11*x24    *x41*x51
     8  +coeff( 80)    *x22    *x42*x53
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 81)    *x24    *x42*x52
     1  +coeff( 82)*x11*x22    *x41*x53
     2  +coeff( 83)    *x24    *x41*x53
     3  +coeff( 84)*x11*x24    *x43    
     4  +coeff( 85)*x11*x24    *x41*x52
     5  +coeff( 86)*x11        *x43*x53
     6  +coeff( 87)*x12*x21    *x42*x52
     7  +coeff( 88)*x11*x23    *x41*x53
     8  +coeff( 89)*x12*x24    *x41*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 90)*x12*x24        *x52
c
      return
      end
      function t_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.2376385E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.68292627E-03,-0.18205313E-01, 0.14938215E-03, 0.11589102E+00,
     +  0.66292146E-02,-0.10701562E-01,-0.11094920E-02, 0.13013072E-02,
     + -0.10396538E-02,-0.31215913E-03, 0.42838923E-03,-0.13466791E-02,
     +  0.19936121E-03,-0.11208061E-02,-0.74825675E-03, 0.56628278E-03,
     +  0.23066039E-03, 0.16035393E-03, 0.31542778E-03,-0.51041675E-03,
     +  0.91528433E-03, 0.55753900E-03,-0.39736729E-03, 0.80345862E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x21    *x41    
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)            *x42*x52
     5  +coeff( 23)    *x22        *x52
     6  +coeff( 24)    *x23    *x41*x51
c
      return
      end
      function y_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.8561194E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.33331865E-02,-0.68076521E-01, 0.22302258E-02, 0.29629113E-02,
     +  0.20340167E-02, 0.36659073E-01, 0.71266075E-02,-0.28998056E-01,
     + -0.76669194E-02, 0.12411295E-02,-0.17667681E-01,-0.85820118E-02,
     +  0.51105400E-02, 0.41224929E-02, 0.12652412E-03,-0.25724825E-02,
     + -0.22744817E-03, 0.92713506E-03, 0.12397190E-02,-0.47538667E-02,
     + -0.57789646E-02,-0.75080566E-03, 0.17805205E-02,-0.44254242E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21    *x44    
     6  +coeff( 15)    *x21    *x43*x51
     7  +coeff( 16)            *x43    
     8  +coeff( 17)            *x42*x51
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)    *x22    *x42*x51
c
      return
      end
      function p_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2922950E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11159079E-02,-0.31140825E-03, 0.19449795E-01,-0.26210607E-02,
     +  0.12609024E-01,-0.13143973E-02, 0.33150772E-02,-0.21291101E-01,
     + -0.35282951E-02,-0.30129619E-02, 0.90805003E-02, 0.35833204E-02,
     +  0.41802675E-02, 0.95609704E-03,-0.94144698E-03, 0.12709551E-02,
     +  0.17023805E-02, 0.29340503E-03,-0.35375843E-03, 0.37934922E-03,
     +  0.41589065E-03, 0.16460194E-02,-0.11935177E-02, 0.95570815E-03,
     +  0.17858882E-02, 0.60606498E-03,-0.59745536E-03,-0.34271550E-03,
     + -0.25134030E-03,-0.91753223E-04, 0.72865526E-03,-0.47385282E-03,
     + -0.60237038E-04, 0.27819595E-04, 0.11903855E-03,-0.89977705E-03,
     +  0.13902069E-04,-0.23938883E-03, 0.39389881E-03,-0.22267488E-03,
     + -0.23031278E-03, 0.10209035E-02, 0.10988211E-02,-0.79265774E-04,
     + -0.34890862E-03, 0.13872725E-02, 0.56843256E-03,-0.59442478E-03,
     + -0.47485036E-03, 0.27541842E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41*x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x21    *x41*x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)            *x41*x53
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)    *x22    *x42*x51
     8  +coeff( 26)*x11*x21    *x42*x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)            *x43*x51
     2  +coeff( 29)            *x42*x52
     3  +coeff( 30)*x11        *x42*x51
     4  +coeff( 31)    *x22    *x43    
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)*x11            *x52
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 36)    *x23    *x41    
     1  +coeff( 37)    *x22        *x52
     2  +coeff( 38)    *x21    *x41*x52
     3  +coeff( 39)*x11*x21    *x42    
     4  +coeff( 40)*x11*x22        *x51
     5  +coeff( 41)*x11*x21        *x52
     6  +coeff( 42)    *x23        *x52
     7  +coeff( 43)    *x22        *x53
     8  +coeff( 44)*x11*x23    *x41    
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 45)*x11*x22        *x52
     1  +coeff( 46)    *x23    *x41*x52
     2  +coeff( 47)    *x23        *x53
     3  +coeff( 48)*x11*x23    *x42    
     4  +coeff( 49)*x12*x23        *x51
     5  +coeff( 50)*x12*x21        *x53
c
      return
      end
      function l_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.7373984E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39978E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34596906E-02,-0.29091799E+00,-0.31716004E-01, 0.89882789E-02,
     + -0.37982002E-01, 0.37658960E-01,-0.16665386E-01,-0.12713615E-01,
     + -0.89171696E-02,-0.41350937E-02,-0.22990014E-02, 0.18310837E-01,
     +  0.90592913E-03,-0.24182234E-01,-0.36522769E-02, 0.18438520E-02,
     + -0.20576222E-02,-0.52549317E-02, 0.59844268E-03,-0.11797555E-01,
     +  0.20646267E-02, 0.28429841E-03,-0.71470550E-03,-0.21573396E-02,
     +  0.15881804E-02,-0.34023961E-02, 0.17257000E-02,-0.10303803E-02,
     +  0.88521180E-03,-0.16002097E-02,-0.30950261E-02,-0.64910818E-02,
     + -0.27964009E-01, 0.39746313E-03, 0.12130234E-03, 0.85886801E-02,
     +  0.66835224E-03,-0.15501109E-02, 0.21502732E-02,-0.49924641E-02,
     +  0.26361865E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x21        *x51
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)            *x42    
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)*x11*x22            
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)*x11        *x43    
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)            *x41    
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)                *x53
     8  +coeff( 26)*x11*x22    *x42    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11        *x42    
     3  +coeff( 30)    *x23        *x51
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x22    *x43    
     6  +coeff( 33)    *x23    *x43    
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)*x11            *x52
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x21    *x43    
     3  +coeff( 39)*x11*x21    *x41*x52
     4  +coeff( 40)*x11*x22    *x43    
     5  +coeff( 41)*x11*x22    *x41*x52
c
      return
      end
      function x_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.1179949E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12986526E-03,-0.10567220E-04,-0.29915418E-01,-0.19968271E-04,
     + -0.39090617E-04, 0.52949516E-04, 0.28759205E-04, 0.82994311E-05,
     + -0.61949931E-05,-0.11855483E-05, 0.32094135E-05,-0.73929964E-05,
     + -0.10290519E-04, 0.38495041E-05,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      x_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21    *x41    
      x_sp_sen    =x_sp_sen    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)*x11                
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)*x11*x21    *x41    
c
      return
      end
      function t_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.9684317E-01/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13081873E-03,-0.79689024E-04,-0.24022397E-01,-0.24724574E-03,
     + -0.53513760E-03, 0.42357238E-03, 0.20327129E-03,-0.19933474E-03,
     + -0.23867195E-04,-0.12438762E-03, 0.18426233E-03, 0.32840300E-04,
     +  0.91896123E-04,-0.32620253E-04, 0.98783195E-04,-0.14937506E-03,
     +  0.23502225E-04,-0.55505610E-04,-0.86497683E-04,-0.83114588E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      t_sp_sen    =t_sp_sen    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)*x11*x21    *x41    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x22        *x51
      t_sp_sen    =t_sp_sen    
     9  +coeff( 18)*x11*x23            
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)*x11*x23    *x41    
c
      return
      end
      function y_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.1015245E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.94468438E-03,-0.17227860E-07, 0.15608674E-05, 0.66365875E-01,
     +  0.39879410E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
c
      return
      end
      function p_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.8532951E-03/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56719035E-03, 0.50793618E-01,-0.10544643E-02,-0.15845473E-03,
     + -0.98367287E-04, 0.18315678E-03, 0.36813351E-03, 0.65436921E-04,
     + -0.43348176E-04, 0.15503194E-03,-0.16648673E-03, 0.28580654E-03,
     +  0.40326311E-04, 0.68247435E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x21        *x51
      p_sp_sen    =p_sp_sen    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)*x11*x22            
c
      return
      end
      function l_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1274153E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73294551E-03, 0.49782002E-05, 0.27576557E-02,-0.17160119E-02,
     + -0.35103349E-03, 0.17931454E-05,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)                *x51
c
      return
      end
      function x_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1790629E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.41872545E-03, 0.24351390E-03,-0.41963819E-01,-0.84961049E-03,
     +  0.83617860E-03, 0.99541433E-03, 0.27145250E-03,-0.19102635E-03,
     + -0.36448083E-03, 0.91489125E-03,-0.44140685E-03,-0.15067335E-03,
     +  0.81644233E-04, 0.20704635E-03,-0.21710174E-03, 0.57749738E-04,
     +  0.44142373E-04,-0.63115323E-04,-0.39976181E-04, 0.14774222E-03,
     + -0.63655258E-04,-0.21875856E-03, 0.48749059E-03,-0.12648918E-03,
     + -0.57508011E-03,-0.68995547E-06,-0.90231326E-04,-0.37557405E-03,
     +  0.14258396E-04, 0.23789060E-04,-0.81685039E-05,-0.89066722E-04,
     +  0.46403752E-05, 0.10992764E-03,-0.80099941E-04, 0.13918013E-03,
     +  0.68650363E-04, 0.17421842E-04,-0.13010632E-03,-0.13742522E-03,
     + -0.72441195E-04,-0.65056025E-04, 0.60243216E-04, 0.42430303E-04,
     + -0.29937761E-04, 0.36353646E-04,-0.27270602E-04, 0.15549853E-04,
     +  0.35011262E-05, 0.70825146E-04, 0.37343871E-04,-0.28461072E-04,
     +  0.20622416E-04,-0.55551723E-05,-0.19808767E-04, 0.15952509E-03,
     + -0.12565264E-03, 0.92507864E-04, 0.10219603E-03, 0.12388194E-04,
     + -0.53156273E-05, 0.44022505E-06,-0.11882915E-04,-0.82512342E-05,
     +  0.75895787E-05,-0.84208214E-05,-0.16978971E-04, 0.51590978E-05,
     +  0.11362852E-04,-0.14006022E-04,-0.12590411E-04,-0.82198549E-05,
     + -0.14411634E-04,-0.12895235E-04,-0.28051147E-04,-0.11178549E-04,
     + -0.16649639E-04,-0.11981630E-04,-0.46411842E-04, 0.96180884E-05,
     +  0.13438520E-04, 0.73679832E-04,-0.17443532E-04, 0.24674118E-04,
     + -0.71967079E-05,-0.30635227E-04, 0.71567033E-05,-0.24767543E-04,
     +  0.35314926E-04,-0.11461236E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      x_sp_sm     =x_sp_sm     
     9  +coeff(  9)    *x24            
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x24    *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)                *x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x43    
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)*x11*x23            
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x23    *x41    
     7  +coeff( 25)    *x24    *x42    
     8  +coeff( 26)    *x21        *x51
      x_sp_sm     =x_sp_sm     
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)    *x24    *x43    
     2  +coeff( 29)*x11*x24    *x42    
     3  +coeff( 30)*x11                
     4  +coeff( 31)    *x21    *x41*x51
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)            *x41*x52
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)*x11*x22    *x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)*x11*x24            
     2  +coeff( 38)*x11        *x43    
     3  +coeff( 39)    *x23    *x43    
     4  +coeff( 40)*x11*x23    *x42    
     5  +coeff( 41)*x11*x22    *x43    
     6  +coeff( 42)*x11*x21    *x42*x52
     7  +coeff( 43)*x12*x23    *x41*x51
     8  +coeff( 44)*x11*x24    *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 45)    *x21    *x42*x51
     1  +coeff( 46)    *x24        *x51
     2  +coeff( 47)            *x43*x51
     3  +coeff( 48)            *x41*x53
     4  +coeff( 49)*x11*x21        *x52
     5  +coeff( 50)    *x24    *x41*x51
     6  +coeff( 51)    *x23    *x42*x51
     7  +coeff( 52)*x11*x21    *x41*x52
     8  +coeff( 53)*x11*x23        *x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 54)*x12*x23    *x41    
     1  +coeff( 55)*x12*x24        *x51
     2  +coeff( 56)*x12*x23    *x41*x52
     3  +coeff( 57)*x12*x21    *x43*x52
     4  +coeff( 58)*x11*x22    *x43*x53
     5  +coeff( 59)*x12*x24    *x42*x52
     6  +coeff( 60)*x11        *x41    
     7  +coeff( 61)    *x21        *x52
     8  +coeff( 62)*x12                
      x_sp_sm     =x_sp_sm     
     9  +coeff( 63)*x11*x21        *x51
     1  +coeff( 64)    *x23        *x51
     2  +coeff( 65)*x11        *x41*x51
     3  +coeff( 66)*x12        *x41    
     4  +coeff( 67)            *x42*x52
     5  +coeff( 68)*x11*x21    *x41*x51
     6  +coeff( 69)    *x23    *x41*x51
     7  +coeff( 70)*x12*x22            
     8  +coeff( 71)    *x22    *x41*x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 72)*x12*x21        *x51
     1  +coeff( 73)    *x21    *x43*x51
     2  +coeff( 74)    *x21    *x41*x53
     3  +coeff( 75)*x11*x22    *x42    
     4  +coeff( 76)*x12        *x41*x51
     5  +coeff( 77)*x11*x22    *x41*x51
     6  +coeff( 78)            *x43*x52
     7  +coeff( 79)    *x23    *x41*x52
     8  +coeff( 80)*x11*x21        *x53
      x_sp_sm     =x_sp_sm     
     9  +coeff( 81)*x12*x22    *x41    
     1  +coeff( 82)*x11*x24    *x41    
     2  +coeff( 83)*x11        *x43*x51
     3  +coeff( 84)    *x22    *x43*x51
     4  +coeff( 85)*x11        *x42*x52
     5  +coeff( 86)    *x22    *x41*x53
     6  +coeff( 87)*x12*x21    *x42    
     7  +coeff( 88)*x12*x21    *x41*x51
     8  +coeff( 89)    *x21    *x43*x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 90)    *x21    *x42*x53
c
      return
      end
      function t_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1610043E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.89768355E-03, 0.82306110E-03,-0.25079347E-01,-0.34176514E-02,
     +  0.21951371E-02, 0.16738209E-02, 0.57043589E-03, 0.19538625E-03,
     + -0.52324153E-03,-0.34736106E-03, 0.20386362E-03, 0.53171866E-03,
     +  0.51216112E-03, 0.12527423E-02,-0.17078209E-03, 0.53791959E-04,
     +  0.82939380E-03, 0.14974247E-03,-0.46633778E-03, 0.91580390E-04,
     + -0.12408308E-03, 0.10042501E-04,-0.57667703E-03,-0.24434691E-03,
     + -0.48921420E-03,-0.28621883E-03,-0.31527638E-03,-0.55807730E-03,
     + -0.38300434E-04, 0.95542135E-04, 0.77843819E-04, 0.46927686E-03,
     + -0.10234798E-03, 0.69582229E-03,-0.79268287E-03,-0.90125622E-03,
     + -0.22871388E-03, 0.30642218E-03, 0.55501075E-03,-0.56533434E-03,
     +  0.15975359E-03,-0.29728444E-04,-0.43041116E-04, 0.46962461E-04,
     +  0.43737451E-04,-0.56469769E-04, 0.55524026E-04, 0.78977435E-04,
     +  0.46753139E-04,-0.37449379E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      t_sp_sm     =t_sp_sm     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)            *x42*x52
     7  +coeff( 16)*x12        *x42    
     8  +coeff( 17)    *x22    *x43    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 18)    *x22    *x42*x52
     1  +coeff( 19)            *x42    
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)*x11*x23            
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)            *x43*x52
     7  +coeff( 25)*x11*x23    *x41    
     8  +coeff( 26)*x11*x21    *x41*x53
      t_sp_sm     =t_sp_sm     
     9  +coeff( 27)*x12*x22    *x42    
     1  +coeff( 28)*x12*x22    *x43    
     2  +coeff( 29)    *x21        *x51
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)*x12        *x41    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)*x11*x21    *x42    
     8  +coeff( 35)    *x23    *x43    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 36)*x11*x23    *x42    
     1  +coeff( 37)*x11*x21    *x42*x53
     2  +coeff( 38)*x11*x23    *x43*x51
     3  +coeff( 39)*x12*x23    *x43    
     4  +coeff( 40)*x12*x23    *x41*x52
     5  +coeff( 41)*x12*x23        *x53
     6  +coeff( 42)*x11        *x41    
     7  +coeff( 43)*x11*x22            
     8  +coeff( 44)*x11            *x52
      t_sp_sm     =t_sp_sm     
     9  +coeff( 45)*x12*x21            
     1  +coeff( 46)            *x43*x51
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)            *x41*x53
     4  +coeff( 49)*x11*x22        *x51
     5  +coeff( 50)*x11        *x42*x51
c
      return
      end
      function y_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1419671E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12375814E-02,-0.77891826E-04, 0.13288310E-04, 0.90704642E-01,
     +  0.39470121E-02,-0.10026846E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
c
      return
      end
      function p_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.8581974E-03/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.39444552E-03, 0.49655240E-01,-0.36798054E-02, 0.16759750E-02,
     +  0.61067479E-03, 0.23285388E-02,-0.45775610E-03, 0.71447727E-03,
     + -0.14209021E-02,-0.18950296E-03, 0.24895702E-03, 0.13766111E-02,
     +  0.83047911E-04, 0.12708911E-03,-0.19315984E-03, 0.15327320E-03,
     +  0.43144665E-03,-0.46314566E-04, 0.33039376E-03, 0.66502536E-04,
     +  0.15504547E-04,-0.16976280E-03, 0.15584986E-03,-0.12229964E-03,
     +  0.11543134E-03,-0.15929306E-03, 0.11665253E-03,-0.90040900E-04,
     + -0.21325753E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23    *x41    
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22    *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)*x11                
     2  +coeff( 11)*x11*x22            
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x22    *x42    
      p_sp_sm     =p_sp_sm     
     9  +coeff( 18)            *x42*x52
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)    *x22    *x42*x52
     3  +coeff( 21)                *x51
     4  +coeff( 22)            *x42    
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)*x11*x23    *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff( 27)*x11*x22    *x42    
     1  +coeff( 28)*x11*x21    *x41*x53
     2  +coeff( 29)*x12*x21    *x41*x53
c
      return
      end
      function l_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.1875851E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10341090E-02,-0.25060746E-04, 0.42738244E-02, 0.11528481E-03,
     + -0.24196862E-02,-0.47784348E-03,-0.25713744E-04,-0.38276707E-04,
     + -0.27809718E-04, 0.20855836E-04,-0.61618994E-05,-0.12735727E-04,
     + -0.16118365E-04, 0.21037529E-04,-0.30942763E-04,-0.69387997E-05,
     + -0.57234065E-05, 0.47039825E-05, 0.79077809E-05, 0.15574085E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_sm     =l_sp_sm     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)                *x52
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)            *x43    
      l_sp_sm     =l_sp_sm     
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x23            
     2  +coeff( 20)*x11*x23    *x41    
c
      return
      end
      function x_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2741545E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12635069E-02, 0.86288073E-03,-0.54637704E-01,-0.34126665E-02,
     +  0.34458342E-02, 0.25412063E-02, 0.68535120E-03, 0.31935587E-03,
     + -0.54658140E-03,-0.11071911E-02, 0.23424951E-02,-0.11730590E-02,
     + -0.47031883E-03, 0.18081391E-03, 0.51118806E-03,-0.55446109E-03,
     +  0.15560605E-03,-0.21818513E-03, 0.25770042E-03,-0.19842171E-03,
     + -0.54135820E-03, 0.11702139E-02,-0.17358171E-03,-0.14226430E-02,
     +  0.64434265E-04, 0.12966428E-05,-0.19748167E-03, 0.30999476E-03,
     + -0.29809409E-03,-0.78409488E-04, 0.11773873E-04,-0.21207395E-03,
     +  0.23383144E-03,-0.28562467E-03, 0.11104241E-03,-0.75191138E-05,
     + -0.30944229E-03,-0.85321342E-03, 0.56897665E-04,-0.12315289E-03,
     + -0.42347543E-04, 0.97294331E-04,-0.72634983E-04,-0.49219354E-04,
     +  0.36259520E-04, 0.10418652E-03, 0.13946752E-03,-0.41030467E-04,
     + -0.60417173E-04, 0.23928849E-03, 0.65184468E-04,-0.10890573E-03,
     +  0.26344109E-03,-0.29614122E-03,-0.34803768E-04,-0.11610422E-03,
     + -0.71887232E-04,-0.37776557E-04,-0.37598363E-05,-0.19663337E-04,
     + -0.10103000E-04,-0.12017749E-04,-0.57704565E-04,-0.25549358E-04,
     +  0.17356694E-03, 0.23927834E-04, 0.14575476E-03, 0.24226560E-04,
     + -0.93858929E-04, 0.75149146E-04, 0.47729541E-05, 0.13507674E-04,
     + -0.48986847E-04, 0.39613355E-04,-0.26205656E-04, 0.48596958E-04,
     +  0.13680174E-03, 0.54041044E-04,-0.66825091E-05, 0.67426256E-04,
     +  0.39561182E-04, 0.57620102E-04,-0.39179577E-04, 0.80743928E-04,
     + -0.26919673E-04, 0.34365901E-04,-0.22398149E-03, 0.15038677E-03,
     + -0.37434096E-04,-0.17108044E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      x_sp_sex    =x_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x24    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)            *x41*x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)*x11*x23            
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)*x11*x23    *x41    
     6  +coeff( 24)    *x24    *x42    
     7  +coeff( 25)*x11                
     8  +coeff( 26)    *x21        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)*x11*x21    *x42    
     2  +coeff( 29)*x11*x23    *x42    
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)            *x41*x52
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)*x11*x22    *x41    
     8  +coeff( 35)*x11*x24            
      x_sp_sex    =x_sp_sex    
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)    *x23    *x43    
     2  +coeff( 38)    *x24    *x43    
     3  +coeff( 39)*x11        *x41    
     4  +coeff( 40)            *x43    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x24        *x51
     7  +coeff( 43)            *x43*x51
     8  +coeff( 44)            *x42*x52
      x_sp_sex    =x_sp_sex    
     9  +coeff( 45)            *x41*x53
     1  +coeff( 46)    *x23    *x41*x51
     2  +coeff( 47)    *x24    *x41*x51
     3  +coeff( 48)            *x43*x52
     4  +coeff( 49)*x11*x21    *x41*x52
     5  +coeff( 50)*x11*x24    *x41    
     6  +coeff( 51)    *x22    *x42*x52
     7  +coeff( 52)*x11*x21    *x42*x52
     8  +coeff( 53)*x12*x23    *x41*x52
      x_sp_sex    =x_sp_sex    
     9  +coeff( 54)*x12*x21    *x43*x52
     1  +coeff( 55)*x11*x21        *x51
     2  +coeff( 56)    *x21    *x42*x51
     3  +coeff( 57)    *x21    *x41*x52
     4  +coeff( 58)*x12        *x41    
     5  +coeff( 59)*x12            *x51
     6  +coeff( 60)*x11*x21    *x41*x51
     7  +coeff( 61)*x11*x21        *x52
     8  +coeff( 62)    *x22        *x53
      x_sp_sex    =x_sp_sex    
     9  +coeff( 63)*x12*x21        *x51
     1  +coeff( 64)*x12        *x41*x51
     2  +coeff( 65)*x11*x21    *x43    
     3  +coeff( 66)*x11*x21    *x42*x51
     4  +coeff( 67)    *x23    *x42*x51
     5  +coeff( 68)*x11*x21        *x53
     6  +coeff( 69)*x12*x22    *x41    
     7  +coeff( 70)    *x22    *x43*x51
     8  +coeff( 71)*x11        *x42*x52
      x_sp_sex    =x_sp_sex    
     9  +coeff( 72)*x11        *x41*x53
     1  +coeff( 73)    *x22    *x41*x53
     2  +coeff( 74)*x12*x21    *x42    
     3  +coeff( 75)*x12*x21        *x52
     4  +coeff( 76)*x11*x23        *x52
     5  +coeff( 77)    *x21    *x43*x52
     6  +coeff( 78)*x12        *x43    
     7  +coeff( 79)*x11*x22    *x43    
     8  +coeff( 80)*x12*x23        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 81)    *x23    *x42*x52
     1  +coeff( 82)*x11*x21    *x41*x53
     2  +coeff( 83)    *x23    *x41*x53
     3  +coeff( 84)*x11*x24        *x52
     4  +coeff( 85)*x11        *x43*x52
     5  +coeff( 86)*x12*x21    *x43    
     6  +coeff( 87)*x11*x23    *x43    
     7  +coeff( 88)*x12*x24    *x41    
     8  +coeff( 89)*x12*x24        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 90)*x11*x22    *x42*x52
c
      return
      end
      function t_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2221833E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20802314E-02, 0.16480152E-02,-0.25442693E-01,-0.66116122E-02,
     +  0.56315339E-02, 0.17597127E-02, 0.79321227E-03, 0.45351972E-03,
     + -0.11071280E-02, 0.30017053E-03,-0.72930433E-03, 0.60231017E-03,
     +  0.10107230E-03, 0.64879103E-03,-0.71233249E-03, 0.15871259E-02,
     + -0.14391076E-03,-0.70647994E-03,-0.89045352E-04, 0.27374244E-04,
     + -0.14992090E-03, 0.10491637E-03, 0.92511630E-03,-0.63792884E-03,
     +  0.99534751E-03,-0.34998910E-03,-0.51878614E-03,-0.30703904E-04,
     + -0.11224294E-02,-0.17760963E-03, 0.10654036E-02,-0.10028256E-02,
     +  0.10225402E-03,-0.45865447E-04, 0.21174533E-03, 0.13155417E-03,
     +  0.10870973E-03,-0.17509308E-03,-0.62714542E-04, 0.64047344E-04,
     +  0.13239252E-03, 0.12385770E-03, 0.15951585E-03,-0.13421923E-03,
     +  0.90880909E-04, 0.34047487E-04,-0.16096741E-03,-0.51934249E-03,
     +  0.11291015E-03,-0.28115403E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      t_sp_sex    =t_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)*x11*x23            
      t_sp_sex    =t_sp_sex    
     9  +coeff( 18)*x11*x23    *x41    
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11            *x52
     5  +coeff( 23)*x11*x21    *x42    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)            *x43*x52
      t_sp_sex    =t_sp_sex    
     9  +coeff( 27)*x11*x22        *x52
     1  +coeff( 28)*x12*x22    *x41    
     2  +coeff( 29)*x11*x23    *x42    
     3  +coeff( 30)    *x22    *x42*x53
     4  +coeff( 31)*x11*x22    *x42*x52
     5  +coeff( 32)*x12*x22    *x43    
     6  +coeff( 33)*x11                
     7  +coeff( 34)*x11        *x41    
     8  +coeff( 35)            *x42*x51
      t_sp_sex    =t_sp_sex    
     9  +coeff( 36)            *x41*x52
     1  +coeff( 37)                *x53
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)*x11*x21        *x51
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)*x12        *x41    
     6  +coeff( 42)    *x21    *x43    
     7  +coeff( 43)    *x23        *x51
     8  +coeff( 44)    *x21    *x42*x51
      t_sp_sex    =t_sp_sex    
     9  +coeff( 45)            *x41*x53
     1  +coeff( 46)*x11            *x53
     2  +coeff( 47)*x12*x22            
     3  +coeff( 48)    *x22    *x42*x51
     4  +coeff( 49)*x12            *x52
     5  +coeff( 50)    *x22        *x53
c
      return
      end
      function y_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.1852641E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12859759E-02,-0.18188675E-03, 0.29027566E-04, 0.11363956E+00,
     +  0.38901742E-02,-0.35011182E-02, 0.68296766E-03, 0.18543696E-02,
     +  0.20657571E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      y_sp_sex    =y_sp_sex    
     9  +coeff(  9)    *x23    *x41    
c
      return
      end
      function p_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.1063462E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35744937E-03, 0.47568046E-01,-0.42276108E-02, 0.16149377E-02,
     +  0.25674829E-02,-0.53457444E-03, 0.61357929E-03, 0.88637817E-03,
     + -0.28721563E-03, 0.21459232E-03,-0.15085072E-02, 0.14201850E-02,
     +  0.28745527E-03, 0.33199602E-04, 0.99134988E-04,-0.22678013E-03,
     +  0.13268647E-03, 0.45717543E-03,-0.66486944E-03,-0.56059929E-04,
     +  0.41866335E-03, 0.81187766E-03, 0.53503038E-04,-0.17825339E-03,
     +  0.17351232E-03,-0.11077688E-03, 0.76486736E-04, 0.11270749E-03,
     + -0.16767849E-03,-0.12687256E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x23    *x41    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x22    *x41    
      p_sp_sex    =p_sp_sex    
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)                *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x21    *x41*x51
      p_sp_sex    =p_sp_sex    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)            *x42*x52
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)    *x23    *x43    
     5  +coeff( 23)    *x22    *x42*x52
     6  +coeff( 24)            *x42    
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x23            
      p_sp_sex    =p_sp_sex    
     9  +coeff( 27)*x11*x21    *x42    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)*x11*x23    *x41    
     3  +coeff( 30)*x11*x21    *x41*x53
c
      return
      end
      function l_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 37)
      data ncoeff/ 36/
      data avdat/ -0.2791778E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14539291E-02,-0.13385137E-03, 0.67096828E-02, 0.61908388E-03,
     + -0.34270466E-02,-0.58458175E-03,-0.10971865E-03,-0.10888673E-03,
     + -0.62348699E-04,-0.39499919E-04,-0.40136718E-04, 0.72229341E-04,
     + -0.14000913E-03,-0.56794266E-04, 0.36445806E-04, 0.38032995E-04,
     +  0.75118519E-05,-0.36612125E-05,-0.55042554E-04, 0.73827623E-05,
     +  0.51155395E-04,-0.80256264E-04, 0.18863284E-04, 0.53333068E-04,
     + -0.85455231E-05, 0.15348789E-04, 0.58718624E-05,-0.44900771E-04,
     +  0.12368737E-04,-0.73477313E-04,-0.22489183E-04, 0.51814623E-04,
     +  0.86266977E-04, 0.21623315E-04, 0.16649918E-04, 0.68361238E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_sex    =l_sp_sex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21        *x51
      l_sp_sex    =l_sp_sex    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)*x11*x23            
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)            *x43*x52
     6  +coeff( 24)*x11*x23    *x41    
     7  +coeff( 25)*x11                
     8  +coeff( 26)*x11*x22            
      l_sp_sex    =l_sp_sex    
     9  +coeff( 27)*x11*x21        *x51
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)*x11*x21    *x42    
     4  +coeff( 31)*x12*x22    *x41    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)*x11*x23    *x42    
     7  +coeff( 34)*x12*x22    *x42    
     8  +coeff( 35)*x12        *x43*x51
      l_sp_sex    =l_sp_sex    
     9  +coeff( 36)*x12*x22    *x43    
c
      return
      end
      function x_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1178958E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23279032E-02, 0.69028638E-01, 0.70177950E-01, 0.54798194E-03,
     +  0.23404376E-02,-0.30422898E-03,-0.31482245E-03,-0.14190312E-03,
     +  0.22051710E-04, 0.84733067E-04,
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
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21        *x52
      x_sp_cq1x   =x_sp_cq1x   
     9  +coeff(  9)                *x51
     1  +coeff( 10)    *x21*x31        
c
      return
      end
      function t_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.2296339E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23176309E-03, 0.34201540E-01,-0.35895897E-04,-0.44594191E-01,
     +  0.49841870E-04,-0.39394934E-04, 0.55368931E-03, 0.12382806E-04,
     +  0.20878282E-02,-0.15950741E-05,-0.11219324E-03,-0.12739027E-03,
     +  0.44465858E-04,-0.17636716E-03, 0.11605780E-02, 0.19918659E-03,
     + -0.41106545E-04, 0.20022196E-03,-0.42604448E-03,-0.19956286E-03,
     + -0.12827340E-02,
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
     3  +coeff(  3)        *x31        
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x51
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff(  9)*x11            *x51
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)            *x41    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x23*x32        
     7  +coeff( 16)*x12*x21    *x41    
     8  +coeff( 17)    *x21*x32*x42*x51
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x31*x41    
     2  +coeff( 20)    *x21*x33*x41    
     3  +coeff( 21)*x11*x22*x32        
c
      return
      end
      function y_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.9726789E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.97567430E-02, 0.71794532E-01, 0.52567977E-01,-0.14819612E-02,
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
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.4125901E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.36679697E-02, 0.36331885E-05, 0.21243397E-01, 0.32554068E-01,
     + -0.10782029E-02, 0.23123217E-03,-0.21351330E-03,-0.43444437E-03,
     +  0.75329859E-04, 0.15432463E-04, 0.72742303E-04,-0.27671913E-03,
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
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)                *x51
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)            *x41*x51
      p_sp_cq1x   =p_sp_cq1x   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)        *x31    *x52
     3  +coeff( 12)*x11*x23    *x41    
c
      return
      end
      function l_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.5100312E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.47756752E-03, 0.27182582E-03, 0.17973357E-03, 0.27314358E-04,
     +  0.33856782E-05,-0.15080121E-02,-0.18202079E-04,-0.72198833E-03,
     +  0.17119064E-04, 0.54792478E-03,-0.15098462E-02, 0.45797438E-04,
     +  0.80532445E-04,-0.34055291E-03, 0.20686258E-02,-0.11049407E-02,
     + -0.16211432E-04, 0.37829936E-03,-0.13746567E-04,-0.46627360E-04,
     +  0.13154068E-04,-0.25917001E-04, 0.43054777E-04,-0.13453490E-02,
     + -0.89057052E-04, 0.25867799E-04,-0.83080886E-06, 0.22340888E-03,
     + -0.28971193E-03, 0.25175787E-02, 0.57973234E-05,-0.10808886E-04,
     + -0.11027701E-02, 0.10658231E-04, 0.12417911E-04,-0.34658566E-04,
     + -0.10567222E-04, 0.14969837E-04,-0.90013682E-05,-0.12721206E-04,
     +  0.24600145E-04,
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
      l_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)        *x32    *x51
     4  +coeff( 13)        *x31*x41*x51
     5  +coeff( 14)        *x31    *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x12                
     8  +coeff( 17)    *x21            
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)        *x31*x42*x51
     2  +coeff( 20)                *x52
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)        *x33        
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)            *x41*x52
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 27)                *x53
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)        *x33*x41    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)*x12            *x51
     7  +coeff( 34)    *x22*x31    *x51
     8  +coeff( 35)        *x33    *x51
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22*x31        
     3  +coeff( 39)*x11*x21*x32        
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)*x11*x21*x31*x41    
c
      return
      end
      function x_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 54)
      data ncoeff/ 53/
      data avdat/ -0.5097397E+01/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.39183921E-02,-0.18073946E+00, 0.61359458E-01, 0.27038695E-02,
     +  0.77599194E-02, 0.21831453E-03, 0.86637934E-04,-0.43740524E-02,
     + -0.11152626E-02,-0.84723271E-02,-0.43670106E-03, 0.21716527E-03,
     +  0.10567388E-03,-0.19093946E-03,-0.59786576E-04,-0.15127295E-03,
     + -0.19170864E-02,-0.54244312E-04, 0.35888972E-03, 0.12251373E-02,
     + -0.29992862E-03,-0.15696074E-03, 0.19503052E-02,-0.35770365E-04,
     + -0.98387842E-04, 0.18608883E-02, 0.24756347E-03, 0.11278344E-02,
     +  0.70606061E-03,-0.47925790E-03, 0.59371168E-03, 0.82692434E-03,
     +  0.82841021E-03,-0.38167182E-04,-0.21521980E-02, 0.90726657E-03,
     + -0.27784379E-03,-0.19268371E-04,-0.30613523E-04, 0.33650638E-04,
     +  0.23207711E-04,-0.86123531E-04, 0.60124050E-04,-0.58779438E-04,
     +  0.24083100E-03,-0.55473234E-03,-0.57182286E-03, 0.52093860E-03,
     + -0.18188913E-02,-0.74383301E-04, 0.17298228E-02, 0.10759631E-02,
     +  0.53071615E-03,
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
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      x_sp_cden   =x_sp_cden   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)*x11            *x51
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)                *x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)    *x23*x32        
      x_sp_cden   =x_sp_cden   
     9  +coeff( 18)    *x21*x32*x41    
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)*x11*x22*x31        
     3  +coeff( 21)*x11*x21*x32        
     4  +coeff( 22)    *x21*x33*x41    
     5  +coeff( 23)*x12*x21*x31*x41    
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)*x11*x22            
      x_sp_cden   =x_sp_cden   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11    *x31*x41    
     2  +coeff( 29)*x11            *x52
     3  +coeff( 30)*x11        *x43    
     4  +coeff( 31)    *x21*x32*x42    
     5  +coeff( 32)*x11*x24*x32        
     6  +coeff( 33)    *x23    *x41    
     7  +coeff( 34)    *x21        *x53
     8  +coeff( 35)*x11*x22    *x41    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)*x11*x22*x31    *x51
     2  +coeff( 38)        *x33*x41*x51
     3  +coeff( 39)*x11    *x31        
     4  +coeff( 40)    *x23        *x51
     5  +coeff( 41)    *x22    *x42    
     6  +coeff( 42)    *x22    *x41*x51
     7  +coeff( 43)*x12            *x51
     8  +coeff( 44)*x11    *x32    *x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 45)*x12    *x31*x41    
     1  +coeff( 46)*x11*x24*x31        
     2  +coeff( 47)    *x21*x32*x43    
     3  +coeff( 48)*x11*x22*x33        
     4  +coeff( 49)*x11*x22*x32*x41    
     5  +coeff( 50)    *x24    *x43    
     6  +coeff( 51)    *x23*x33*x41    
     7  +coeff( 52)*x11*x24*x31*x41    
     8  +coeff( 53)    *x21*x33*x43    
c
      return
      end
      function t_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.1298158E+01/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50902099E-03, 0.95163032E-01,-0.41511271E-01, 0.58628125E-02,
     + -0.25261512E-02, 0.53420750E-03, 0.22276046E-02,-0.15574138E-03,
     +  0.10691692E-03,-0.27665327E-03,-0.13590235E-02,-0.79948426E-03,
     +  0.20808033E-02, 0.46002115E-02,-0.15419634E-02,-0.43556638E-03,
     + -0.15152665E-02,-0.76736452E-03, 0.61608857E-03,-0.44152202E-04,
     +  0.86782486E-04,-0.85760345E-03,-0.36641720E-03, 0.10844576E-03,
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
      t_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x51
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)        *x31        
      t_sp_cden   =t_sp_cden   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21        *x52
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)        *x32        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)*x12                
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x21*x31*x41    
      t_sp_cden   =t_sp_cden   
     9  +coeff( 18)*x11*x23            
     1  +coeff( 19)            *x43    
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)*x11*x21        *x51
     4  +coeff( 22)            *x42    
     5  +coeff( 23)        *x33        
     6  +coeff( 24)*x11*x21*x31        
c
      return
      end
      function y_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.4873689E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.43988763E-02, 0.90335887E-02, 0.64736232E-01,-0.73136250E-02,
     +  0.92114741E-02, 0.57317889E-02, 0.69538215E-02, 0.55424082E-02,
     + -0.45447858E-03, 0.42541500E-03,-0.11664922E-02,-0.11584302E-02,
     + -0.11598447E-01, 0.26625972E-01,-0.16638331E-01,-0.63187699E-03,
     + -0.23314485E-02,
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
      x51 = x5
c
c                  function
c
      y_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)        *x31*x41    
     5  +coeff(  5)            *x42    
     6  +coeff(  6)        *x31    *x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x21*x31        
      y_sp_cden   =y_sp_cden   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)                *x51
     3  +coeff( 12)    *x21            
     4  +coeff( 13)        *x33        
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)    *x21*x31    *x51
     8  +coeff( 17)    *x22*x31        
c
      return
      end
      function p_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 37)
      data ncoeff/ 36/
      data avdat/  0.7591116E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.63035609E-02, 0.62215123E-02,-0.70744306E-01,-0.82376301E-02,
     + -0.23706285E-02,-0.22563227E-01,-0.18910162E-01,-0.67023881E-01,
     +  0.41285954E-01, 0.27453618E-01,-0.25378990E-02, 0.77325692E-02,
     +  0.31135972E-02,-0.55113495E-02,-0.68593980E-02,-0.23811024E-04,
     + -0.46305363E-02, 0.41405432E-03,-0.49906963E-03, 0.12260364E-01,
     + -0.95521845E-03, 0.16066801E-01,-0.10586130E-01,-0.48332158E-03,
     +  0.10090739E-01, 0.38454297E-02,-0.15124473E-03,-0.20977719E-01,
     +  0.20336743E-01,-0.37415705E-02,-0.42043670E-03,-0.71379676E-03,
     +  0.10450153E-02,-0.10474495E-02, 0.22804421E-03, 0.49255299E-03,
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
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cden   =p_sp_cden   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x32        
     2  +coeff( 11)                *x51
     3  +coeff( 12)*x11    *x31        
     4  +coeff( 13)    *x22            
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)        *x33        
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11*x21*x31        
      p_sp_cden   =p_sp_cden   
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x22*x33        
     2  +coeff( 20)        *x31    *x51
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)        *x32*x41    
     5  +coeff( 23)        *x31*x42    
     6  +coeff( 24)        *x31    *x52
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)    *x21    *x41*x51
      p_sp_cden   =p_sp_cden   
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)*x11        *x41*x51
     4  +coeff( 31)    *x21*x33        
     5  +coeff( 32)    *x23    *x41    
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)    *x22*x31    *x51
     8  +coeff( 35)    *x21*x31    *x52
      p_sp_cden   =p_sp_cden   
     9  +coeff( 36)*x11*x21*x32        
c
      return
      end
      function l_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.2090388E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10209983E-01, 0.35107616E+00,-0.12013501E+00,-0.31077741E-01,
     + -0.68690251E-02, 0.12132378E-02, 0.53477462E-03,-0.50137914E-02,
     +  0.16759811E-01,-0.45382869E-02,-0.59656603E-02, 0.84605981E-02,
     +  0.17570401E-01, 0.61547436E-03,-0.98142039E-03,-0.11291183E-02,
     +  0.23422015E-03,-0.93087065E-03,
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
c
c                  function
c
      l_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)        *x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      l_sp_cden   =l_sp_cden   
     9  +coeff(  9)*x11            *x51
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)        *x31*x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)*x11*x21        *x51
     8  +coeff( 17)                *x51
      l_sp_cden   =l_sp_cden   
     9  +coeff( 18)    *x23            
c
      return
      end
      function x_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 39)
      data ncoeff/ 38/
      data avdat/ -0.3959992E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.79159858E-02, 0.38615480E+00, 0.20844664E-02, 0.13851066E+00,
     + -0.17830612E+00, 0.33795852E-01, 0.73881622E-03, 0.25929360E-01,
     + -0.59239967E-02,-0.12585418E-01, 0.11127192E-03, 0.22194730E-02,
     + -0.15353381E-01,-0.19418783E-02, 0.23572873E-01, 0.56645680E-04,
     + -0.18361636E-02, 0.56216274E-02, 0.35385997E-03, 0.81616491E-02,
     + -0.15432357E-01,-0.67647663E-02, 0.33978108E-03,-0.59704180E-03,
     +  0.20785849E-02, 0.17702263E-02,-0.31451698E-03, 0.68259111E-03,
     +  0.17474948E-03, 0.29350075E-03,-0.62615150E-02,-0.71242510E-03,
     + -0.73086791E-03,-0.30857923E-02,-0.13230129E-02, 0.45971436E-03,
     +  0.27509587E-03,-0.26753507E-03,
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)    *x24            
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x23            
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)*x11*x21*x31        
     6  +coeff( 24)    *x22*x32        
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)*x11*x22*x31        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 27)    *x24        *x51
     1  +coeff( 28)    *x23*x32        
     2  +coeff( 29)    *x22    *x43    
     3  +coeff( 30)    *x21*x33*x41    
     4  +coeff( 31)*x11*x22*x31*x41    
     5  +coeff( 32)*x12*x23            
     6  +coeff( 33)*x12*x21    *x42*x51
     7  +coeff( 34)    *x21*x32*x42    
     8  +coeff( 35)        *x31        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 36)*x11    *x31        
     1  +coeff( 37)        *x32    *x51
     2  +coeff( 38)        *x31*x41*x51
c
      return
      end
      function t_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 23)
      data ncoeff/ 22/
      data avdat/ -0.2856553E-04/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13903702E-02,-0.11428230E+00,-0.33677719E-04, 0.22788001E-01,
     +  0.38867537E-01,-0.11834206E-01, 0.24102013E-03,-0.14358137E-02,
     +  0.72355075E-02, 0.72635501E-03, 0.92582626E-03,-0.42802789E-02,
     +  0.60396781E-03, 0.39705049E-03,-0.16428344E-02,-0.31504143E-03,
     + -0.19488193E-03, 0.15976615E-03, 0.11088566E-03, 0.30986019E-03,
     +  0.83729369E-03,-0.16903678E-03,
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
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)        *x32    *x51
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 18)        *x32    *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)*x12*x21        *x51
c
      return
      end
      function y_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/  0.8666045E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72235470E-02,-0.15533943E+00, 0.10603420E+00,-0.88947769E-02,
     +  0.94082113E-02,-0.36086746E-01, 0.98467320E-01,-0.63358627E-01,
     + -0.38208885E-02,-0.28245430E-01,-0.31567771E-01, 0.13556974E+00,
     + -0.31707659E+00, 0.47082562E-01, 0.40227543E-02, 0.22135919E-01,
     + -0.27954642E-01,-0.56840125E-02, 0.11331489E-02,-0.37227995E-02,
     + -0.20903067E-02, 0.17171958E-01, 0.18080588E+00,-0.54754842E-01,
     +  0.18091239E-01, 0.45280188E-01,-0.20323485E-01,
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
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
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
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)        *x33        
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)    *x22            
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x42    
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21*x32*x41    
     2  +coeff( 20)*x11*x21*x31        
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)        *x31*x42    
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)            *x41*x51
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff( 27)    *x22    *x41    
c
      return
      end
      function p_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/  0.2475173E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19877907E-02, 0.21980733E-02,-0.33719823E-01, 0.14635437E-01,
     + -0.93178137E-03,-0.54958300E-02,-0.93317544E-02, 0.43759448E-02,
     + -0.20483476E-02, 0.33641937E-02, 0.85785933E-03, 0.29215729E-02,
     + -0.12963292E-02,-0.17793191E-02, 0.17078067E-02, 0.14193046E-03,
     + -0.19623952E-02, 0.20364618E-03, 0.70764851E-02,-0.33834818E-03,
     + -0.77989200E-04, 0.37870908E-02, 0.33326314E-02,-0.34065023E-02,
     +  0.57090411E-03,-0.14117222E-02,-0.20504419E-02,-0.21722792E-03,
     + -0.25472612E-03,-0.59460234E-02,-0.27730565E-02, 0.25300245E-03,
     +  0.25988277E-02, 0.22867505E-03, 0.28890485E-03,-0.12415042E-03,
     +  0.26914091E-02,-0.12935746E-03,-0.78780408E-03,-0.14001707E-02,
     +  0.27041425E-03,-0.55949675E-03,-0.86689153E-03,-0.16783226E-03,
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
      p_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)    *x22            
     3  +coeff( 12)*x11    *x31        
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)        *x33        
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)*x11*x21*x31        
     8  +coeff( 17)        *x32        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 18)    *x21        *x51
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x23            
     4  +coeff( 22)        *x32*x41    
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)                *x51
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)    *x22*x32        
     2  +coeff( 29)    *x21*x31    *x52
     3  +coeff( 30)    *x22    *x41    
     4  +coeff( 31)        *x31*x42    
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)    *x21*x33        
     8  +coeff( 35)    *x23    *x41    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 36)        *x31    *x53
     1  +coeff( 37)*x11        *x41    
     2  +coeff( 38)    *x21*x32        
     3  +coeff( 39)*x11    *x31    *x51
     4  +coeff( 40)    *x22*x31    *x51
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)    *x22*x31*x42    
     7  +coeff( 43)*x11*x23*x31        
     8  +coeff( 44)*x11*x22*x31    *x51
c
      return
      end
      function l_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.5902833E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.55175265E-02,-0.49472320E+00, 0.35973601E-02,-0.32495260E-01,
     +  0.19102311E+00,-0.53215947E-01, 0.38724551E-02,-0.14330918E-03,
     + -0.51951781E-02, 0.44521950E-02,-0.18386341E-02, 0.29015495E-01,
     + -0.20658571E-01,-0.47818097E-03,-0.15588992E-02, 0.59700860E-02,
     +  0.13756087E-02,-0.38725517E-02,-0.16099584E-02,-0.17491962E-03,
     +  0.96900435E-03,-0.62641175E-02, 0.69129886E-02,
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
      l_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x32        
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31*x41*x51
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff( 18)    *x21*x31    *x51
     1  +coeff( 19)        *x31        
     2  +coeff( 20)                *x52
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)    *x23    *x43    
     5  +coeff( 23)    *x23*x33*x41    
c
      return
      end
      function x_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 73)
      data ncoeff/ 72/
      data avdat/  0.6104513E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51965890E-02, 0.24591079E+00, 0.19926038E-02, 0.32966915E+00,
     + -0.16953263E+00,-0.13757939E-02, 0.34884531E-01,-0.22467298E-01,
     + -0.97348867E-03, 0.18434724E-01, 0.71323570E-02,-0.13346652E-01,
     +  0.26609234E-02, 0.41582491E-02, 0.26635178E-02,-0.18831210E-01,
     +  0.28579853E-01, 0.23545089E-03,-0.13038441E-01, 0.84748752E-02,
     + -0.34916896E-03,-0.22483527E-03,-0.98940905E-03,-0.30609425E-02,
     +  0.10763162E-01,-0.87083867E-02,-0.27355790E-03,-0.20150468E-02,
     +  0.23713452E-02,-0.24506981E-02, 0.31036632E-02,-0.56403279E-02,
     +  0.13069005E-03,-0.79253996E-02,-0.19142843E-02,-0.24232338E-02,
     + -0.99569246E-04,-0.33957392E-03, 0.60794451E-02, 0.97216800E-03,
     +  0.88964002E-02,-0.20760608E-03,-0.63831038E-02, 0.23454039E-03,
     + -0.20733180E-01, 0.69712894E-03,-0.34301991E-02,-0.17806733E-02,
     + -0.63351140E-03,-0.21245760E-03,-0.48365924E-03,-0.18703284E-02,
     + -0.15293077E-02, 0.11777836E-02,-0.14425388E-03, 0.21260760E-03,
     + -0.15228351E-02, 0.18057232E-02, 0.65339038E-04,-0.86497385E-02,
     +  0.18967079E-01,-0.10533996E-01,-0.11742421E-03, 0.89874183E-03,
     + -0.10990940E-01, 0.24725737E-01,-0.16458066E-01, 0.44479170E-02,
     +  0.66669549E-04, 0.36649529E-02, 0.11055162E-02, 0.53617242E-02,
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
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x22            
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)                *x53
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)    *x21*x31*x41    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)    *x24            
     2  +coeff( 20)            *x42    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)        *x31    *x51
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)    *x23            
     8  +coeff( 26)*x11*x22            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 27)*x11*x21        *x51
     1  +coeff( 28)    *x22*x31    *x51
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)*x11*x22        *x51
     6  +coeff( 33)            *x43*x51
     7  +coeff( 34)        *x32    *x52
     8  +coeff( 35)    *x23*x32        
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 36)    *x23        *x52
     1  +coeff( 37)*x11    *x33        
     2  +coeff( 38)    *x21*x33*x41    
     3  +coeff( 39)*x11*x23        *x51
     4  +coeff( 40)*x11*x24*x31        
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)    *x21*x31    *x52
     7  +coeff( 43)*x11*x24            
     8  +coeff( 44)    *x22*x32    *x51
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 45)    *x21    *x42    
     1  +coeff( 46)    *x21*x31    *x51
     2  +coeff( 47)    *x22*x31    *x52
     3  +coeff( 48)    *x21*x32*x42    
     4  +coeff( 49)    *x21*x32    *x52
     5  +coeff( 50)        *x32*x42*x51
     6  +coeff( 51)        *x32    *x53
     7  +coeff( 52)*x11*x21        *x53
     8  +coeff( 53)*x11*x24*x31*x41    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 54)            *x41*x51
     1  +coeff( 55)*x11    *x31        
     2  +coeff( 56)        *x33        
     3  +coeff( 57)    *x21*x33        
     4  +coeff( 58)    *x21    *x42*x51
     5  +coeff( 59)        *x33*x41    
     6  +coeff( 60)    *x24        *x51
     7  +coeff( 61)        *x31*x41*x52
     8  +coeff( 62)            *x42*x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 63)    *x23*x31*x41    
     1  +coeff( 64)*x11        *x43    
     2  +coeff( 65)*x11    *x32    *x51
     3  +coeff( 66)*x11    *x31*x41*x51
     4  +coeff( 67)*x11        *x42*x51
     5  +coeff( 68)    *x22    *x41*x52
     6  +coeff( 69)*x11            *x53
     7  +coeff( 70)*x12*x21        *x51
     8  +coeff( 71)*x12            *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 72)*x12*x23            
c
      return
      end
      function t_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.2376385E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55087812E-03, 0.14761161E-01,-0.11059318E-01, 0.11799233E+00,
     +  0.79269186E-02,-0.10805414E-01,-0.32949228E-01, 0.10036961E-01,
     + -0.25153169E-02, 0.95729474E-02, 0.57931994E-02,-0.49046648E-03,
     + -0.63093641E-03,-0.91572008E-04,-0.62157749E-02,-0.50528126E-03,
     +  0.57217578E-03, 0.54583512E-03, 0.19752239E-02,-0.53071808E-02,
     + -0.58806883E-02,-0.39456965E-03,-0.76640025E-03,
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
      t_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)        *x31        
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)        *x32    *x51
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)                *x53
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff( 18)        *x32    *x52
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)    *x22        *x52
     5  +coeff( 23)*x11*x22        *x51
c
      return
      end
      function y_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.8561194E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.76427991E-02,-0.13263530E+00, 0.76093599E-01,-0.63723177E-02,
     +  0.75302636E-02, 0.51221255E-01,-0.11182309E+00, 0.62370908E-01,
     + -0.32002407E-02, 0.23548878E-02,-0.26998140E-01,-0.20185735E-01,
     +  0.34902006E-01,-0.84845684E-01, 0.48397351E-01,-0.74186147E-03,
     +  0.13494731E-02,-0.17841479E-01, 0.35053276E-01, 0.37102064E-02,
     +  0.15527426E-01,-0.53391871E-02,-0.24528243E-02,
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
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff(  9)*x11                
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)        *x33        
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)    *x21*x31    *x51
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)        *x31    *x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)*x11    *x31        
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)    *x22        *x51
c
      return
      end
      function p_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2922950E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.28641182E-02,-0.24984756E-02, 0.53010944E-01,-0.38554695E-01,
     +  0.93132886E-03,-0.79233776E-03, 0.22670709E-01,-0.69326557E-01,
     +  0.42522773E-01, 0.49395589E-04, 0.36560197E-02, 0.30616019E-02,
     +  0.23355128E-02, 0.31836280E-02, 0.19526426E-02, 0.27605381E-01,
     + -0.23736434E-01,-0.69318863E-03, 0.14556956E-03, 0.74844022E-03,
     + -0.48170221E-03,-0.14600404E-03, 0.82390725E-04,-0.77719740E-02,
     +  0.50232909E-02, 0.10506097E-01, 0.54325006E-03,-0.17550233E-02,
     +  0.94023440E-03,-0.47738737E-03,-0.88283111E-03, 0.96869143E-03,
     +  0.19842750E-03, 0.19193867E-04,-0.37379310E-03,-0.78793784E-03,
     +  0.18728377E-03,-0.76154392E-03,-0.40770471E-03,-0.39862408E-03,
     +  0.44589108E-03,-0.73745274E-02, 0.30653688E-03,-0.11375468E-02,
     + -0.11238361E-01, 0.12369077E-01, 0.22924405E-02,-0.19806386E-02,
     + -0.80241077E-03,-0.18245943E-02,
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
      p_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x43    
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)*x11*x23*x31        
     5  +coeff( 14)                *x51
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)        *x32        
     8  +coeff( 17)            *x41*x51
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11    *x31        
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)        *x31    *x53
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)        *x33    *x51
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)        *x33        
     8  +coeff( 26)    *x22    *x41    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)    *x22*x32*x41    
     3  +coeff( 30)    *x22            
     4  +coeff( 31)    *x21        *x51
     5  +coeff( 32)                *x52
     6  +coeff( 33)    *x21*x32        
     7  +coeff( 34)                *x53
     8  +coeff( 35)    *x21*x33        
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 36)    *x23    *x41    
     1  +coeff( 37)    *x21*x31    *x52
     2  +coeff( 38)            *x41*x53
     3  +coeff( 39)    *x21*x33    *x51
     4  +coeff( 40)    *x21*x31    *x53
     5  +coeff( 41)    *x23            
     6  +coeff( 42)        *x32*x41    
     7  +coeff( 43)        *x32    *x51
     8  +coeff( 44)*x11*x21        *x51
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 45)    *x22*x32        
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)        *x33*x41    
     3  +coeff( 48)        *x32*x42    
     4  +coeff( 49)        *x31*x43    
     5  +coeff( 50)    *x22*x31    *x51
c
      return
      end
      function l_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.6849965E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.55033369E-02,-0.49491718E+00,-0.31197719E-01, 0.19046438E+00,
     + -0.57875227E-01, 0.79192361E-02,-0.90502016E-02, 0.76532853E-02,
     + -0.22234073E-01,-0.86048506E-02, 0.15038287E-02, 0.29290289E-01,
     +  0.15556106E-01,-0.16263634E-02, 0.91834627E-02,-0.17375793E-01,
     +  0.99678095E-02, 0.15049648E-02,-0.29565852E-01, 0.19928824E-01,
     + -0.10723788E-02,-0.90422778E-03, 0.16170931E-02, 0.28409036E-02,
     + -0.54314954E-03,
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
      x53 = x52*x5
c
c                  function
c
      l_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)        *x32        
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)        *x33        
     2  +coeff( 11)        *x31        
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)            *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)        *x32*x41    
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff( 18)                *x53
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)*x12*x23            
c
      return
      end
      function x_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 50)
      data ncoeff/ 49/
      data avdat/  0.1437417E-01/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.75659929E-02, 0.29664576E+00,-0.99351136E-02, 0.73512805E+00,
     + -0.28234470E+00, 0.62342767E-01,-0.59807289E-01, 0.10469047E-01,
     +  0.21411458E-01, 0.81554111E-02,-0.16637562E-01,-0.21527032E-02,
     +  0.67993985E-02,-0.13790219E-02, 0.68961559E-02,-0.48865704E-02,
     + -0.24211258E-03, 0.59486721E-02,-0.29341504E-02, 0.28576098E-01,
     + -0.11206361E+00,-0.34770563E-02, 0.17711909E+00, 0.52372282E-02,
     + -0.62450059E-02, 0.28099369E-02,-0.26940506E-01, 0.22662000E-02,
     +  0.12721433E-02,-0.43153581E-02,-0.52708732E-02, 0.88664079E-02,
     + -0.21398462E-01,-0.46919845E-02,-0.22388510E-03,-0.25672307E-02,
     + -0.33951629E-03,-0.20989177E-02,-0.19032592E-02, 0.21165700E-02,
     +  0.12347292E-01,-0.15281799E-04,-0.99266553E-02, 0.28276999E-03,
     +  0.21782299E-02,-0.20820084E-02,-0.73408730E-01,-0.20320870E-02,
     + -0.11572782E-02,
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
      x_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)        *x31        
      x_sp_cfp    =x_sp_cfp    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x32    *x51
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x21    *x41    
     8  +coeff( 17)            *x41*x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x22        *x52
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)*x11*x22        *x51
     7  +coeff( 25)    *x24        *x51
     8  +coeff( 26)        *x32    *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 27)        *x31*x41    
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)    *x21*x32    *x51
     3  +coeff( 30)    *x23        *x52
     4  +coeff( 31)*x11*x22*x32        
     5  +coeff( 32)    *x23            
     6  +coeff( 33)    *x21    *x42    
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)*x11    *x31*x41*x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 36)    *x21*x33*x41    
     1  +coeff( 37)        *x32*x43    
     2  +coeff( 38)        *x32    *x53
     3  +coeff( 39)*x11    *x32    *x52
     4  +coeff( 40)*x12*x23*x31        
     5  +coeff( 41)            *x42    
     6  +coeff( 42)    *x21*x31    *x51
     7  +coeff( 43)    *x23        *x51
     8  +coeff( 44)        *x33*x41    
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 45)        *x31    *x53
     1  +coeff( 46)            *x41*x53
     2  +coeff( 47)*x12*x22            
     3  +coeff( 48)*x11*x24            
     4  +coeff( 49)    *x22        *x53
c
      return
      end
      function t_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.2376542E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55064959E-03, 0.14760122E-01,-0.11060535E-01, 0.11799386E+00,
     +  0.79252627E-02,-0.10805252E-01,-0.32947969E-01, 0.10037976E-01,
     + -0.25155516E-02, 0.94905607E-02, 0.57867523E-02,-0.49290346E-03,
     + -0.63407241E-03,-0.91736867E-04,-0.62080733E-02,-0.50452171E-03,
     +  0.56866120E-03, 0.54915901E-03, 0.19752756E-02,-0.52728849E-02,
     + -0.58310791E-02,-0.39534803E-03,-0.76401839E-03,
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
      t_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)        *x31        
      t_sp_cfp    =t_sp_cfp    
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)        *x32    *x51
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)                *x53
      t_sp_cfp    =t_sp_cfp    
     9  +coeff( 18)        *x32    *x52
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)    *x22        *x52
     5  +coeff( 23)*x11*x22        *x51
c
      return
      end
      function y_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/ -0.1608319E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37252971E-02, 0.82986720E-01,-0.95140368E-01, 0.11377689E-01,
     + -0.45907199E-02,-0.32494411E-01, 0.90713762E-01,-0.60148131E-01,
     + -0.94035991E-01, 0.69736794E-01, 0.29942424E-02,-0.91109145E-02,
     +  0.19737794E-02, 0.69599282E-02,-0.21094589E+00, 0.12756954E+01,
     + -0.21662061E+01, 0.11279694E+01, 0.12508747E+00,-0.91826856E-01,
     + -0.74868329E-01, 0.34047756E-01,-0.73286588E-02,-0.12107715E-02,
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
      y_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cfp    =y_sp_cfp    
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)        *x33        
     7  +coeff( 16)        *x32*x41    
     8  +coeff( 17)        *x31*x42    
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)        *x32    *x51
     2  +coeff( 20)        *x31*x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)        *x31    *x52
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)                *x53
c
      return
      end
      function p_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2923254E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.29086184E-02,-0.23018096E-02, 0.56080144E-01,-0.42095210E-01,
     +  0.60300383E-03, 0.37218686E-02, 0.17325111E-01, 0.26151387E-01,
     + -0.13392090E-01,-0.15959289E-01, 0.91092632E-03, 0.37990932E-02,
     +  0.46929106E-03, 0.38557516E-02, 0.18205750E-02,-0.13115953E-01,
     + -0.49971906E-02,-0.26238773E-02,-0.17842014E-03, 0.13759147E-02,
     + -0.11898729E-03,-0.45785523E-03,-0.43721512E-03,-0.73291287E-02,
     +  0.49032184E-03, 0.98010115E-02,-0.15110447E-03, 0.38782130E-02,
     +  0.12492425E-02,-0.14019954E-03,-0.72595780E-02,-0.41656796E-03,
     + -0.30033197E-03, 0.15177023E-03, 0.20438091E-03,-0.12244742E-02,
     + -0.52009546E-03,-0.22125515E-03, 0.51900244E-03, 0.81058315E-04,
     +  0.79909328E-03,-0.22139408E-03,-0.58298273E-03, 0.10075546E-02,
     +  0.16077103E-03,-0.79292513E-03, 0.71040203E-03,-0.65681717E-03,
     + -0.14637415E-03,-0.11971498E-03,
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
      p_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x43    
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)*x11*x23*x31        
     5  +coeff( 14)                *x51
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)        *x32        
     8  +coeff( 17)            *x41*x51
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11    *x31        
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)        *x31    *x53
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)        *x33    *x51
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)        *x33        
     8  +coeff( 26)    *x22    *x41    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)    *x22*x33        
     3  +coeff( 30)*x12    *x31*x41    
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)    *x22*x32        
     6  +coeff( 33)    *x21*x33        
     7  +coeff( 34)    *x22*x31    *x51
     8  +coeff( 35)    *x21*x31    *x52
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 36)            *x41*x53
     1  +coeff( 37)    *x21*x33    *x51
     2  +coeff( 38)    *x21*x31    *x53
     3  +coeff( 39)            *x42*x53
     4  +coeff( 40)*x11    *x33*x41    
     5  +coeff( 41)*x12*x22    *x41    
     6  +coeff( 42)        *x33*x41*x52
     7  +coeff( 43)                *x52
     8  +coeff( 44)*x12                
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 45)    *x21*x32        
     1  +coeff( 46)    *x21        *x52
     2  +coeff( 47)*x11            *x52
     3  +coeff( 48)    *x23    *x41    
     4  +coeff( 49)    *x21        *x53
     5  +coeff( 50)*x11*x21        *x52
c
      return
      end
      function l_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.1386956E-01/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.46948E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22247E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15340016E-02,-0.49542394E+00,-0.33125766E-01, 0.19127467E+00,
     + -0.58178667E-01,-0.32401524E-01,-0.57595596E-02,-0.10988595E-03,
     +  0.58505004E-02, 0.24885547E-02, 0.29597878E-01, 0.42446773E-02,
     +  0.35595193E-02, 0.41206642E-02,-0.13331049E-02,-0.12665669E-02,
     + -0.35287256E-02,-0.17019850E-02,-0.12172788E-02, 0.14447203E-02,
     +  0.25467377E-02,-0.53337933E-02,-0.48807566E-02, 0.10686594E-02,
     +  0.26354443E-02,
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
      l_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x52
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31        
      l_sp_cfp    =l_sp_cfp    
     9  +coeff(  9)                *x53
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)        *x32    *x51
     5  +coeff( 14)        *x33    *x51
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)    *x21        *x52
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 18)*x11*x21        *x51
     1  +coeff( 19)    *x23        *x51
     2  +coeff( 20)    *x21    *x42*x51
     3  +coeff( 21)            *x41    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)        *x32*x41*x51
     6  +coeff( 24)    *x21        *x53
     7  +coeff( 25)    *x23*x31*x41    
c
      return
      end
      function x_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/ -0.3488090E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.76194233E-02, 0.57868373E+00, 0.94185881E-02, 0.13184604E+00,
     + -0.25044280E+00, 0.13301482E-03, 0.27759058E-01,-0.76482417E-02,
     +  0.18829500E-01,-0.33774499E-02,-0.66588179E-03,-0.49899689E-04,
     +  0.10366191E-02, 0.18293473E-02,-0.58230213E-02, 0.12104967E-01,
     + -0.10135247E-02,-0.39666747E-02, 0.48171406E-03,-0.54783807E-02,
     +  0.10548783E-01,-0.22648855E-02, 0.20652004E-02, 0.36637604E-03,
     + -0.39123595E-02,-0.14789660E-02, 0.41348725E-02, 0.12505273E-02,
     +  0.50100905E-03,-0.20070108E-03,-0.10259192E-01,-0.85533848E-02,
     + -0.28136346E-03, 0.32168915E-02, 0.13077750E-02,-0.16078033E-02,
     + -0.61591546E-03, 0.60284097E-03,-0.39684265E-02,
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
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x31        
      x_sp_cdex   =x_sp_cdex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)        *x33*x41    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)        *x32        
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21        *x52
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)        *x31    *x51
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x23*x32        
     4  +coeff( 22)    *x24            
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)*x12*x21            
     8  +coeff( 26)    *x21*x33        
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 27)    *x21*x32*x41    
     1  +coeff( 28)*x12            *x51
     2  +coeff( 29)    *x22*x31*x42    
     3  +coeff( 30)    *x21*x33*x41    
     4  +coeff( 31)*x11*x22*x32        
     5  +coeff( 32)*x11*x22*x31*x41    
     6  +coeff( 33)            *x42    
     7  +coeff( 34)    *x23            
     8  +coeff( 35)*x11    *x31        
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 36)    *x22        *x51
     1  +coeff( 37)        *x31    *x52
     2  +coeff( 38)            *x41*x52
     3  +coeff( 39)    *x21*x32*x42    
c
      return
      end
      function t_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.5345973E+00/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.79004880E-03,-0.11354718E+00, 0.17439484E-03, 0.22769039E-01,
     +  0.38684666E-01,-0.89330683E-02, 0.20987610E-02,-0.12068608E-02,
     + -0.38488553E-03, 0.56456742E-02,-0.53109284E-02, 0.48493169E-03,
     + -0.88993961E-03, 0.13725765E-02, 0.72147925E-02,-0.31551374E-02,
     +  0.23180233E-03,-0.24296036E-03,-0.38325670E-02,-0.26078790E-03,
     +  0.39648282E-03,-0.88095909E-03, 0.96793054E-04, 0.57204784E-03,
     +  0.53855404E-03,
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
      t_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sp_cdex   =t_sp_cdex   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)        *x32    *x51
      t_sp_cdex   =t_sp_cdex   
     9  +coeff( 18)    *x21*x31        
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)        *x31*x42*x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)        *x33    *x51
     7  +coeff( 25)    *x23*x31*x41    
c
      return
      end
      function y_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.6287763E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.52732937E-02,-0.12702617E+00, 0.97353756E-01,-0.80547109E-02,
     +  0.85519329E-02, 0.70087470E-01,-0.16107512E+00, 0.94463237E-01,
     +  0.43022301E-01,-0.37409947E-02, 0.18093650E-02,-0.29209642E-01,
     + -0.21495596E-01,-0.33521382E-02, 0.17223561E-01,-0.30255336E-02,
     + -0.39839894E-02,-0.18258097E-01, 0.63264486E-03, 0.36046740E-02,
     + -0.47639303E-03,-0.32048933E-02,-0.22564754E-02,
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
     1  +coeff( 10)*x11                
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x21*x31    *x51
      y_sp_cdex   =y_sp_cdex   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x23    *x41    
     5  +coeff( 23)*x11*x23*x31    *x51
c
      return
      end
      function p_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.2396113E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20830778E-02, 0.15493822E-02,-0.32189809E-01, 0.12069057E-01,
     + -0.59015531E-03,-0.31509225E-02,-0.86720595E-02, 0.12493253E-02,
     + -0.50118787E-03, 0.33684182E-02, 0.37861097E-04,-0.13325355E-02,
     +  0.15092780E-02, 0.36473793E-03,-0.74253511E-03,-0.46262547E-03,
     + -0.48261660E-03, 0.44947350E-03,-0.10408828E-02,-0.12474804E-02,
     +  0.13013890E-03, 0.54886244E-03, 0.59588440E-02,-0.12750132E-03,
     + -0.10044823E-02, 0.90804562E-04, 0.14727442E-02,-0.15785196E-02,
     + -0.18946707E-03,-0.76544886E-04,-0.45412118E-02, 0.84455934E-03,
     +  0.21590739E-03,-0.20259740E-03,-0.10024138E-03, 0.28985513E-02,
     + -0.10211968E-03, 0.13766075E-02, 0.80360107E-04, 0.16336254E-03,
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
      p_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)    *x21*x31    *x51
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11    *x31        
     5  +coeff( 14)        *x33        
     6  +coeff( 15)*x11*x23*x31        
     7  +coeff( 16)*x12    *x33        
     8  +coeff( 17)        *x32        
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)*x11*x23            
     4  +coeff( 22)    *x22            
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)        *x31    *x52
     7  +coeff( 25)*x11*x21*x31        
     8  +coeff( 26)    *x23*x31        
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 27)*x12        *x41    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x22*x31    *x51
     3  +coeff( 30)*x11*x21            
     4  +coeff( 31)    *x22    *x41    
     5  +coeff( 32)    *x23    *x41    
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)    *x21*x31    *x52
     8  +coeff( 35)        *x31    *x53
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 36)*x11        *x41    
     1  +coeff( 37)    *x21*x32        
     2  +coeff( 38)    *x22*x32        
     3  +coeff( 39)        *x31*x43    
     4  +coeff( 40)    *x23    *x41*x51
c
      return
      end
      function l_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.1220003E-01/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.56346604E-02,-0.78470778E+00, 0.11882194E-01,-0.10072772E+00,
     +  0.31689569E+00,-0.41343320E-01,-0.38003575E-02,-0.12601081E-01,
     +  0.16036464E-01, 0.18489134E-01, 0.10910170E-02,-0.47537383E-01,
     +  0.25458856E-01,-0.25199790E-01, 0.15410097E-01,-0.26256375E-01,
     +  0.22776879E-01,-0.90477802E-02, 0.11417731E-03, 0.89499465E-03,
     +  0.82955364E-03,
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
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_sp_cdex   =l_sp_cdex   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)        *x32        
     2  +coeff( 11)                *x52
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)*x12                
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x42    
      l_sp_cdex   =l_sp_cdex   
     9  +coeff( 18)        *x31        
     1  +coeff( 19)            *x43    
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)    *x21        *x52
c
      return
      end
