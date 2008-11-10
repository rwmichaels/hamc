C Forward transfer functions for right hrs with septum based on s5a_dir.dat
c HRS + PREX room temperature septum (right side)
c                     -JJL 7/30/08


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
      function x_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  8)
      data ncoeff/  7/
      data avdat/  0.1188560E+00/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37323244E-03,-0.18724430E-04,-0.31158429E-01,-0.61022944E-04,
     +  0.59076585E-04, 0.12867334E-04, 0.40358996E-05,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      x_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
c
      return
      end
      function t_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.9747868E-01/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53280435E-03,-0.17382552E-03,-0.25011219E-01,-0.70029759E-03,
     +  0.18654413E-03, 0.85888583E-04,-0.63997781E-04,-0.48294671E-04,
     + -0.50924013E-04, 0.67486879E-04, 0.25477413E-04,-0.11976348E-04,
     + -0.14086846E-04, 0.15457241E-04, 0.22179374E-04,-0.71859245E-05,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x23            
      t_sp_sen    =t_sp_sen    
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x22    *x43    
     6  +coeff( 15)    *x23    *x43    
     7  +coeff( 16)    *x21    *x43    
c
      return
      end
      function y_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.7581788E-03/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72351907E-03,-0.71881473E-05, 0.66434227E-01,-0.18467449E-03,
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
c          set up monomials   functions
      x21 = x2
      x41 = x4
c
c                  function
c
      y_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
c
      return
      end
      function p_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.7002255E-03/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54072577E-03, 0.51371388E-01,-0.76089817E-03,-0.10468567E-03,
     +  0.58859554E-04, 0.11959753E-03,-0.11247696E-03, 0.56211269E-04,
     + -0.23050099E-04, 0.10306026E-03, 0.27099337E-04, 0.21429563E-04,
     +  0.13945641E-04,-0.24231011E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x22    *x41    
      p_sp_sen    =p_sp_sen    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x23    *x43    
c
      return
      end
      function l_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1302797E-02/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.63202257E-03, 0.35450248E-05, 0.29043339E-02,-0.17179090E-02,
     + -0.38120802E-03, 0.28926179E-05,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
c
      return
      end
      function x_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.1799058E+00/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48627629E-03, 0.14867581E-04,-0.43281436E-01, 0.65161628E-04,
     +  0.35393369E-03, 0.92617345E-04,-0.55165208E-04, 0.36175734E-04,
     + -0.78375066E-04,-0.16599525E-03, 0.63480053E-04,-0.99215307E-04,
     +  0.24785969E-03,-0.14658405E-03, 0.36427322E-04,-0.65509528E-04,
     +  0.26551916E-04, 0.68284302E-04,-0.26512185E-04,-0.55615768E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)            *x42    
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x24    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x24    *x42    
     6  +coeff( 15)    *x22            
     7  +coeff( 16)    *x23    *x42    
     8  +coeff( 17)    *x21    *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)    *x23    *x43    
     2  +coeff( 20)    *x24    *x43    
c
      return
      end
      function t_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.1634534E+00/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.36099529E-04, 0.24439383E-03,-0.25618961E-01, 0.65506191E-03,
     +  0.50225691E-03, 0.39006362E-03, 0.21170665E-03, 0.75636170E-04,
     + -0.14843247E-03, 0.12617536E-03,-0.19066659E-03, 0.13665338E-03,
     +  0.51820236E-04,-0.97734228E-04, 0.12319337E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x22    *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      t_sp_sm     =t_sp_sm     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x22    *x43    
c
      return
      end
      function y_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1042889E-02/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.99375797E-03,-0.51507184E-04, 0.91575466E-01,-0.58780430E-03,
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
c          set up monomials   functions
      x21 = x2
      x41 = x4
c
c                  function
c
      y_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
c
      return
      end
      function p_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.5671039E-03/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51187351E-03, 0.52666925E-01,-0.13187571E-02, 0.34926637E-03,
     + -0.16581564E-03, 0.18185412E-03, 0.22464215E-03, 0.54439472E-03,
     + -0.27820893E-03, 0.32738902E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      p_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x23    *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x42    
c
      return
      end
      function l_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1900513E-02/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87309466E-03,-0.22809227E-05, 0.44363430E-02,-0.23805639E-02,
     + -0.53183892E-03,-0.58287992E-05, 0.53313192E-05,-0.74474833E-05,
     + -0.50335270E-05,-0.10539096E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      l_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      l_sp_sm     =l_sp_sm     
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x22    *x42    
c
      return
      end
      function x_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.2757092E+00/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23163721E-03, 0.26926174E-03,-0.55834062E-01, 0.12513703E-02,
     +  0.10012292E-03, 0.75123779E-03, 0.19759031E-03,-0.10312464E-03,
     + -0.16573553E-03, 0.91742011E-04,-0.17677790E-03, 0.60238136E-03,
     + -0.36650780E-03, 0.15163903E-03,-0.15477849E-03,-0.35440581E-03,
     +  0.54332366E-04, 0.14185028E-03,-0.10653086E-03,-0.45043529E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      x_sp_sex    =x_sp_sex    
     9  +coeff(  9)    *x24            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x24    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23    *x42    
     7  +coeff( 16)    *x24    *x42    
     8  +coeff( 17)    *x21    *x43    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)    *x24    *x43    
     2  +coeff( 20)    *x23    *x43    
c
      return
      end
      function t_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.2262530E+00/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.87436219E-03, 0.68296707E-03,-0.26512470E-01, 0.35347044E-02,
     +  0.10901342E-03, 0.42389025E-03, 0.16555993E-03, 0.46196557E-03,
     + -0.37697377E-04, 0.15646023E-03, 0.20399608E-03,-0.14984079E-03,
     + -0.21449316E-03, 0.61618040E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22    *x42    
      t_sp_sex    =t_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x21    *x43    
c
      return
      end
      function y_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1342142E-02/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12638627E-02,-0.97658231E-04, 0.11678908E+00,-0.10685000E-02,
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
c          set up monomials   functions
      x21 = x2
      x41 = x4
c
c                  function
c
      y_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
c
      return
      end
      function p_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.8102019E-03/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51823992E-03, 0.50633494E-01,-0.11281696E-02, 0.74019126E-03,
     + -0.13998860E-03, 0.25836486E-03, 0.93368042E-04, 0.35972826E-04,
     + -0.29974082E-03, 0.23461669E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      p_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23    *x41    
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      p_sp_sex    =p_sp_sex    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x42    
c
      return
      end
      function l_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.2783210E-02/
      data xmin/
     1  0.00000E+00,-0.52020E-01, 0.00000E+00,-0.29973E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.51993E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11814088E-02,-0.52880572E-04, 0.68873218E-02,-0.32744755E-02,
     + -0.69446064E-03,-0.20969106E-04, 0.13871773E-04,-0.17346734E-04,
     + -0.38450205E-04,-0.94078687E-05,-0.98461696E-05, 0.17482245E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      l_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      l_sp_sex    =l_sp_sex    
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
c
      return
      end
      function x_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.1644420E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14141826E-02, 0.14290738E+00,-0.18682294E-02,-0.37562812E-03,
     +  0.36283693E-03, 0.13183795E-02, 0.94675273E-03,-0.69899863E-03,
     + -0.19036471E-03,-0.10097099E-03, 0.40320511E-03,-0.53781469E-03,
     +  0.81232924E-03,-0.26807920E-04,-0.25591182E-04,-0.39324156E-03,
     + -0.17389143E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23    *x41    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21    *x42    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff(  9)    *x24            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x24    *x41    
     4  +coeff( 13)    *x23    *x43    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x24    *x42    
c
      return
      end
      function t_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.7995618E-04/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11260309E-03,-0.10176018E-01,-0.41335673E-03,-0.54452149E-03,
     +  0.69488684E-03,-0.56891862E-04,-0.31243995E-03, 0.11025835E-03,
     + -0.32744513E-03,-0.29446060E-05, 0.92882015E-04, 0.49884235E-04,
     + -0.12014406E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23    *x42    
     5  +coeff(  5)    *x23    *x41    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x22    *x43    
     4  +coeff( 13)    *x23    *x43    
c
      return
      end
      function y_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.1449794E-01/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.86393310E-02, 0.12606063E+00,-0.14492661E-02,-0.64986073E-02,
     + -0.21614390E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
c
c                  function
c
      y_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22    *x41    
c
      return
      end
      function p_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.6013809E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.35145064E-02,-0.72684936E-03, 0.52931421E-01,-0.30749957E-02,
     + -0.11238578E-02,-0.46315792E-04, 0.10781990E-03,-0.22515963E-03,
     + -0.33391680E-03, 0.22975454E-03, 0.26395466E-03,-0.10002103E-02,
     + -0.38937858E-03, 0.52035402E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x21    *x43    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x42    
c
      return
      end
      function l_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.1123046E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19522836E-02,-0.10944255E-04,-0.43136869E-02,-0.36344177E-02,
     + -0.21845947E-02, 0.27217460E-03, 0.86444583E-04, 0.62356557E-04,
     + -0.10715291E-04, 0.44129523E-04,-0.43925811E-04, 0.11620600E-03,
     + -0.12494361E-04,-0.53398231E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      l_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22    *x42    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x43    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x23    *x42    
c
      return
      end
      function x_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.2996161E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18998638E-02, 0.19673347E+00,-0.87234395E-04,-0.77609211E-03,
     +  0.16790485E-02, 0.69862482E-03,-0.44872835E-02,-0.29393071E-02,
     +  0.48638671E-02,-0.18942271E-02,-0.12574262E-02,-0.47267959E-03,
     + -0.19692369E-03, 0.63166418E-03,-0.53427526E-03,-0.16252226E-03,
     +  0.62747509E-03, 0.12268220E-02,-0.41185060E-04,-0.15605814E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x43    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21    *x42    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x23    *x42    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x24            
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x23    *x43    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)    *x24    *x42    
     2  +coeff( 20)    *x24    *x43    
c
      return
      end
      function t_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.8726257E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.44780114E-03, 0.47984928E-01,-0.13330900E-02,-0.13414321E-02,
     + -0.25340146E-03,-0.58712170E-03, 0.46034297E-03, 0.17855800E-02,
     +  0.10788308E-03,-0.73848356E-03, 0.10008808E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x23    *x41    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x23    *x42    
     2  +coeff( 11)    *x22    *x42    
c
      return
      end
      function y_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.2008815E-01/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11688099E-01, 0.17596337E+00,-0.23197662E-02,-0.10280775E-01,
     + -0.41991938E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
c
c                  function
c
      y_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22    *x41    
c
      return
      end
      function p_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.2312955E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14145254E-02, 0.16227552E-03,-0.19697649E-01, 0.79550367E-03,
     +  0.10701868E-03,-0.30088349E-03, 0.29447000E-03, 0.30780138E-04,
     + -0.30236835E-04,-0.12809162E-03, 0.10447978E-03,-0.82288418E-04,
     + -0.12087290E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x22    *x42    
     8  +coeff(  8)            *x42    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x23    *x42    
c
      return
      end
      function l_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/ -0.1887019E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26645211E-02, 0.11237672E-05,-0.37511154E-02,-0.48387502E-02,
     + -0.37165245E-02, 0.44610057E-03, 0.10808307E-03, 0.12985426E-03,
     + -0.16809847E-04,-0.14895649E-04, 0.23059626E-03, 0.56352186E-04,
     + -0.31479569E-04,-0.43057702E-04, 0.76894867E-04,-0.15557263E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      l_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22    *x42    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x22    *x43    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x23    *x43    
c
      return
      end
      function x_sp_den     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.5097396E+01/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18472280E-02,-0.11459999E+00, 0.25357101E-02, 0.29348705E-02,
     +  0.24809076E-02, 0.55129989E-03,-0.12760541E-02,-0.35467283E-02,
     +  0.10328604E-02, 0.28028595E-03,-0.53125399E-03, 0.15320389E-02,
     +  0.54097502E-03, 0.27788181E-04, 0.10917297E-03,-0.11539580E-03,
     + -0.83360530E-03,-0.63701125E-04,-0.17577615E-03, 0.10331895E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_den     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x23    *x41    
      x_sp_den     =x_sp_den     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x24    *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x22    *x43    
      x_sp_den     =x_sp_den     
     9  +coeff( 18)    *x24    *x42    
     1  +coeff( 19)    *x23    *x43    
     2  +coeff( 20)    *x24    *x43    
c
      return
      end
      function t_sp_den     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/  0.3612545E+01/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.42533036E-01, 0.68007201E+00, 0.17380103E+00,-0.10196894E-01,
     + -0.22128209E-01, 0.16447881E-01, 0.30989170E-01,-0.17225152E-01,
     +  0.26779411E-01,-0.10170405E-01, 0.30668424E-02,-0.10965238E-01,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_den     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x21    *x42    
      t_sp_den     =t_sp_den     
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x23    *x42    
c
      return
      end
      function y_sp_den     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1122525E-01/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.64281803E-02, 0.10118288E+00,-0.23863325E-02, 0.32201845E-02,
     + -0.67921230E-02,-0.61127385E-02, 0.78459480E-03,-0.19332048E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      y_sp_den     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x43    
c
      return
      end
      function p_sp_den     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.7253873E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.38309072E-02, 0.43382668E-02,-0.57711404E-01, 0.42972015E-02,
     + -0.21210851E-01,-0.11175048E-01, 0.55458804E-03, 0.84412267E-03,
     +  0.66136452E-03,-0.95222960E-03, 0.14163859E-02, 0.82810660E-03,
     + -0.73279371E-03, 0.79496740E-03,-0.79158484E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_den     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x23            
      p_sp_den     =p_sp_den     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x23    *x42    
c
      return
      end
      function l_sp_den     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.6188721E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.76292530E-02, 0.22150622E+00,-0.45719263E-02,-0.14579500E-01,
     + -0.53480482E-02,-0.44846619E-02,-0.46779639E-02, 0.22146762E-02,
     +  0.67160856E-02,-0.17581257E-02, 0.11278550E-02,-0.27935561E-02,
     +  0.10573665E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      l_sp_den     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x22    *x41    
      l_sp_den     =l_sp_den     
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x22    *x43    
c
      return
      end
      function x_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/ -0.3152489E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.86884771E-03, 0.30341056E+00, 0.12141213E-01,-0.81259683E-02,
     + -0.10284689E-01,-0.28954260E-02, 0.11644668E-01,-0.12059503E-02,
     +  0.31125546E-02,-0.24435660E-02,-0.39124527E-03,-0.64371125E-03,
     + -0.50229831E-02, 0.14888630E-02,-0.43595769E-03, 0.80584502E-03,
     +  0.28396263E-02,-0.40504080E-03,-0.33912286E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x23            
      x_sp_dex    =x_sp_dex    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x24            
     2  +coeff( 11)            *x43    
     3  +coeff( 12)            *x41    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x22    *x43    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 18)    *x24    *x42    
     1  +coeff( 19)    *x24    *x43    
c
      return
      end
      function t_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/  0.5936885E+00/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.89536194E-03,-0.96786618E-01, 0.22584465E-02, 0.61492715E-03,
     +  0.26310079E-02,-0.46748208E-03,-0.30770660E-02,-0.62307494E-03,
     +  0.11361797E-02, 0.11748408E-03,-0.34062951E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x22    *x41    
      t_sp_dex    =t_sp_dex    
     9  +coeff(  9)    *x23    *x42    
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x22    *x43    
c
      return
      end
      function y_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.4437830E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.31422672E-02, 0.53553045E-01, 0.50140088E-02,-0.40406458E-01,
     + -0.32793875E-02,-0.65109262E-03,-0.24444377E-01, 0.25080277E-02,
     +  0.96141739E-03, 0.29710245E-02,-0.46953447E-02, 0.35102065E-02,
     +  0.15711713E-02,-0.17931686E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
c
c                  function
c
      y_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x44    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x23            
      y_sp_dex    =y_sp_dex    
     9  +coeff(  9)            *x45    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x22    *x42    
c
      return
      end
      function p_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.6914783E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28940185E-03, 0.13907087E-02,-0.29157556E-02, 0.18085563E-03,
     + -0.91780331E-02, 0.43802231E-03, 0.43324466E-03,-0.51918677E-02,
     + -0.66993624E-03, 0.12221587E-02, 0.39363025E-04, 0.32456958E-03,
     + -0.24485184E-03,-0.15750188E-03, 0.22164613E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      p_sp_dex    =p_sp_dex    
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x23    *x43    
c
      return
      end
      function l_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.4811145E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40943502E-02,-0.43798047E+00,-0.28220754E-01, 0.11438667E-01,
     + -0.18524256E-02, 0.12620423E-01,-0.21328509E-02,-0.15592709E-01,
     + -0.28756440E-02, 0.80677652E-03, 0.57997955E-02,-0.18030908E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x23    *x41    
      l_sp_dex    =l_sp_dex    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23    *x42    
     3  +coeff( 12)    *x22    *x42    
c
      return
      end
      function x_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/  0.1586056E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18258899E-02, 0.18948151E+00, 0.17557595E-01,-0.51607294E-02,
     + -0.74611525E-02,-0.25196194E-02, 0.78830197E-02, 0.15778140E-02,
     + -0.20642972E-02,-0.26665503E-03,-0.19592071E-03,-0.33990911E-02,
     +  0.19687959E-02, 0.59046684E-03, 0.32210562E-03,-0.23339062E-02,
     +  0.83420268E-03,-0.61861932E-03, 0.30212808E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x22    *x41    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff(  9)    *x24            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)            *x41    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)    *x21    *x43    
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)    *x24    *x43    
     8  +coeff( 17)    *x22    *x42    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 18)    *x24    *x42    
     1  +coeff( 19)    *x23            
c
      return
      end
      function t_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.4090091E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27847470E-03,-0.72157845E-01,-0.39366912E-02, 0.18896004E-02,
     +  0.13720599E-02, 0.35478637E-03,-0.22582356E-02,-0.61976811E-03,
     +  0.11313686E-03, 0.13221963E-03, 0.90075948E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x22    *x41    
      t_sp_q3en   =t_sp_q3en   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23    *x42    
c
      return
      end
      function y_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.4018309E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30081538E-02, 0.53183392E-01, 0.63780746E-02,-0.50072625E-01,
     + -0.34226580E-02,-0.77511132E-03,-0.29168444E-01, 0.29444082E-02,
     +  0.11057055E-02, 0.54810750E-02, 0.35662607E-02,-0.56271171E-02,
     +  0.19932003E-02,-0.19783592E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
c
c                  function
c
      y_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x44    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x23            
      y_sp_q3en   =y_sp_q3en   
     9  +coeff(  9)            *x45    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x22    *x42    
c
      return
      end
      function p_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.3861708E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.89865716E-04, 0.14815822E-02, 0.11352463E-03, 0.89524183E-04,
     + -0.98956982E-02, 0.50888804E-03, 0.46336575E-03,-0.62064757E-02,
     + -0.77786791E-03, 0.96780626E-03, 0.33370062E-03,-0.68567271E-04,
     + -0.21563839E-03, 0.19736908E-03,-0.56062676E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x22    *x43    
     5  +coeff( 14)    *x23    *x43    
     6  +coeff( 15)    *x23    *x42    
c
      return
      end
      function l_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.5517609E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53385654E-02,-0.28619632E+00,-0.21836793E-02,-0.25882991E-01,
     +  0.76154219E-02,-0.37282687E-02, 0.75277048E-02,-0.99299867E-02,
     + -0.16774572E-02, 0.38646851E-03, 0.34903085E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23    *x41    
      l_sp_q3en   =l_sp_q3en   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23    *x42    
c
      return
      end
      function x_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.1130598E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15888128E-02, 0.97823545E-01, 0.13303694E-01,-0.62335874E-02,
     + -0.26718623E-02, 0.30147817E-03, 0.83534443E-03,-0.22285082E-02,
     +  0.52299147E-03,-0.21224902E-02,-0.15331164E-03, 0.52002757E-02,
     + -0.26766597E-02, 0.11239180E-02, 0.28409055E-03, 0.13923994E-03,
     + -0.16920757E-03,-0.48442822E-03, 0.13100064E-02,-0.15999530E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x24            
      x_sp_q3m    =x_sp_q3m    
     9  +coeff(  9)    *x21    *x43    
     1  +coeff( 10)    *x23    *x42    
     2  +coeff( 11)    *x23    *x43    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)    *x23            
     8  +coeff( 17)            *x43    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 18)    *x24    *x42    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)    *x24    *x43    
c
      return
      end
      function t_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.9397559E-04/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.83821593E-04,-0.41007005E-01, 0.11094816E-02, 0.38766669E-03,
     + -0.80115150E-03,-0.57764089E-03, 0.22479633E-03,-0.42797090E-03,
     + -0.24485285E-03,-0.80931064E-03, 0.55048364E-03, 0.33262861E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      t_sp_q3m    =t_sp_q3m    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x23    *x42    
c
      return
      end
      function y_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.3007710E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24804010E-02, 0.46947777E-01, 0.76933461E-02, 0.32980850E-02,
     + -0.57962764E-01,-0.51732217E-02,-0.29036268E-02,-0.34528010E-01,
     +  0.32092347E-02, 0.63090990E-02, 0.23054543E-02,-0.17730613E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      y_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)            *x42    
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x43    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x22    *x41    
      y_sp_q3m    =y_sp_q3m    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x22    *x42    
c
      return
      end
      function p_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.1275981E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78862353E-03,-0.49253798E-03,-0.12909058E-01, 0.94349548E-03,
     +  0.49069826E-02, 0.24074421E-02,-0.36154079E-03,-0.33448281E-03,
     +  0.49476733E-03,-0.49692666E-03, 0.28678917E-03,-0.24445224E-03,
     +  0.11161844E-03,-0.10552263E-03, 0.12293714E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x42    
      p_sp_q3m    =p_sp_q3m    
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x22    *x43    
c
      return
      end
      function l_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.6287594E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.61055254E-02,-0.28615424E+00,-0.21667702E-02,-0.28776851E-01,
     +  0.76519102E-02,-0.37416788E-02, 0.74087507E-02,-0.99571282E-02,
     + -0.15951457E-02, 0.36734561E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23    *x41    
      l_sp_q3m    =l_sp_q3m    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23    *x42    
c
      return
      end
      function x_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.1285031E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22065639E-02, 0.58846593E-01, 0.95934595E-03, 0.16194604E-01,
     + -0.42569889E-02,-0.83371811E-02,-0.35465818E-02, 0.53154971E-02,
     +  0.54974778E-03,-0.16431458E-02, 0.19850142E-02, 0.53354999E-03,
     +  0.39759369E-03,-0.20041624E-02,-0.40197806E-03, 0.74152095E-03,
     + -0.59176667E-03,-0.16275408E-03, 0.13527314E-02,-0.17152557E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x24            
     8  +coeff(  8)    *x23    *x41    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x24    *x41    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x23    *x43    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x24    *x42    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)    *x24    *x43    
c
      return
      end
      function t_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.1963517E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.46971260E-03,-0.20013655E-01, 0.51218580E-03, 0.21476198E-02,
     + -0.13886576E-02,-0.19714409E-02, 0.40934613E-03, 0.34020975E-03,
     +  0.11002687E-02,-0.29456685E-03, 0.75126166E-03, 0.38068162E-03,
     + -0.72550197E-05,-0.69373188E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x23    *x43    
c
      return
      end
      function y_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.4795043E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.79713296E-03, 0.17294960E-01, 0.51178010E-02, 0.24150936E-02,
     + -0.36616158E-01,-0.37918712E-02,-0.85385930E-03,-0.50439587E-03,
     + -0.22366066E-01, 0.18452428E-02,-0.63378851E-04, 0.40962840E-02,
     +  0.15326652E-02,-0.10104737E-02, 0.68027910E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
c
c                  function
c
      y_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)            *x42    
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x43    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)            *x44    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21    *x44    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)            *x45    
c
      return
      end
      function p_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.1862606E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13011559E-02,-0.22853021E-01, 0.15606731E-02, 0.18476494E-01,
     + -0.10935896E-02,-0.10926145E-02, 0.10444402E-01, 0.14637440E-03,
     +  0.16376949E-02,-0.23381212E-02,-0.18861920E-02, 0.56329521E-03,
     + -0.65934571E-03, 0.46029009E-03,-0.43684145E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21    *x42    
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x21            
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x22    *x43    
     6  +coeff( 15)    *x23    *x43    
c
      return
      end
      function l_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.6512252E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.63208840E-02,-0.28612861E+00,-0.20980169E-02,-0.29324466E-01,
     +  0.75665321E-02,-0.40096724E-02, 0.77267378E-02,-0.99984538E-02,
     + -0.14986288E-02, 0.35946539E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23    *x41    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23    *x42    
c
      return
      end
      function x_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.1849760E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38154849E-02, 0.12945750E-02, 0.23631311E-02, 0.25019847E-01,
     + -0.13134528E-03,-0.80703245E-02, 0.13780674E-02, 0.14191330E-03,
     + -0.14008746E-01,-0.67226868E-02, 0.67031505E-02, 0.40238951E-02,
     + -0.22311250E-02, 0.13002659E-02, 0.69584162E-03,-0.94572065E-03,
     + -0.93882199E-03,-0.18683435E-03, 0.17319495E-02,-0.23137880E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      x_sp_fp     =x_sp_fp     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x21    *x43    
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)    *x24    *x42    
     8  +coeff( 17)    *x23    *x43    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)    *x24    *x43    
c
      return
      end
      function t_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.1963577E-03/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.46972375E-03,-0.20013280E-01, 0.51219109E-03, 0.21476890E-02,
     + -0.13886783E-02,-0.19714755E-02, 0.40932014E-03, 0.34015204E-03,
     +  0.11003034E-02,-0.29459054E-03, 0.75132481E-03, 0.38070814E-03,
     + -0.74276236E-05,-0.69369777E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      t_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      t_sp_fp     =t_sp_fp     
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x23    *x43    
c
      return
      end
      function y_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.4882433E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29908926E-02,-0.48695300E-01, 0.16833602E-01, 0.35506459E-02,
     +  0.29179602E-03, 0.81135603E-02,-0.92455966E-03,-0.13295314E-02,
     + -0.15762027E-02,-0.11751843E-02, 0.17732368E-02, 0.12654557E-02,
     + -0.18406743E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      y_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x43    
     8  +coeff(  8)    *x23            
      y_sp_fp     =y_sp_fp     
     9  +coeff(  9)    *x21            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x23    *x41    
c
      return
      end
      function p_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.1862606E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13011607E-02,-0.22853123E-01, 0.15606756E-02, 0.18476726E-01,
     + -0.10935984E-02,-0.10925996E-02, 0.10444529E-01, 0.14639380E-03,
     +  0.16377149E-02,-0.23381594E-02,-0.18863326E-02, 0.56327024E-03,
     + -0.65939379E-03, 0.46032626E-03,-0.43679305E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21    *x42    
      p_sp_fp     =p_sp_fp     
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x21            
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x22    *x43    
     6  +coeff( 15)    *x23    *x43    
c
      return
      end
      function l_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.6863522E-02/
      data xmin/
     1  0.00000E+00,-0.50132E-01, 0.00000E+00,-0.26299E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.50362E-01, 0.00000E+00, 0.18431E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.66504558E-02,-0.28618294E+00,-0.19392535E-02,-0.29915592E-01,
     +  0.73159421E-02,-0.46674726E-02, 0.88622654E-02,-0.10098234E-01,
     +  0.26819581E-03,-0.12696118E-02, 0.28544185E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23    *x41    
      l_sp_fp     =l_sp_fp     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23    *x42    
c
      return
      end
      function x_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.1069893E-03/
      data xmin/
     1  0.00000E+00,-0.49437E-01, 0.00000E+00,-0.24265E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.49466E-01, 0.00000E+00, 0.17427E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21655473E-04, 0.11964185E+00,-0.77875261E-03,-0.18507478E-03,
     +  0.25391270E-03, 0.31848933E-03, 0.37851796E-03, 0.49670209E-03,
     + -0.25068913E-03,-0.11876808E-03, 0.19527243E-03,-0.20089356E-03,
     + -0.12231429E-03,-0.53779182E-04, 0.21252428E-03,-0.30504458E-04,
     +  0.24141744E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      x_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x23    *x41    
      x_sp_col    =x_sp_col    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x23    *x42    
     3  +coeff( 12)    *x24    *x41    
     4  +coeff( 13)    *x24    *x42    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x22    *x43    
c
      return
      end
      function t_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.2877919E-03/
      data xmin/
     1  0.00000E+00,-0.49437E-01, 0.00000E+00,-0.24265E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.49466E-01, 0.00000E+00, 0.17427E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10450845E-04, 0.46518188E-01,-0.97771781E-03,-0.14618853E-03,
     +  0.10615338E-03, 0.22706302E-03, 0.60993835E-03, 0.12919022E-03,
     + -0.65926601E-04,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x23            
      t_sp_col    =t_sp_col    
     9  +coeff(  9)    *x21    *x42    
c
      return
      end
      function y_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.9406090E-02/
      data xmin/
     1  0.00000E+00,-0.49437E-01, 0.00000E+00,-0.24265E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.49466E-01, 0.00000E+00, 0.17427E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.64322422E-03, 0.50915066E-01,-0.15523047E-02,-0.31810816E-03,
     + -0.50937588E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
c
c                  function
c
      y_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21            
     5  +coeff(  5)    *x22    *x41    
c
      return
      end
      function p_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.3460658E-02/
      data xmin/
     1  0.00000E+00,-0.49437E-01, 0.00000E+00,-0.24265E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.49466E-01, 0.00000E+00, 0.17427E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.71326958E-03,-0.43912689E-03, 0.21617435E-01,-0.19616447E-02,
     + -0.34446325E-03,-0.19145336E-03,-0.33174592E-03, 0.67393405E-04,
     +  0.70421083E-04,-0.15825406E-03,-0.11708584E-03, 0.15585222E-03,
     +  0.12786932E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      p_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x22    *x42    
     8  +coeff(  8)            *x42    
      p_sp_col    =p_sp_col    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x23    *x42    
c
      return
      end
      function l_sp_col     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.3519796E-03/
      data xmin/
     1  0.00000E+00,-0.49437E-01, 0.00000E+00,-0.24265E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.00000E+00, 0.49466E-01, 0.00000E+00, 0.17427E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11204100E-02, 0.11199602E-04,-0.45803515E-02,-0.29398352E-02,
     + -0.53614093E-03, 0.81046754E-04, 0.27625254E-04,-0.18395278E-04,
     + -0.20874988E-04, 0.22561895E-04, 0.56879903E-05,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_col     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      l_sp_col     =l_sp_col     
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x21    *x42    
c
      return
      end
      function x_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/  0.8805449E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.26676455E-02, 0.69660664E-01, 0.69676369E-01,-0.27584154E-03,
     + -0.29681003E-03, 0.75380725E-04, 0.11498701E-04,-0.17257515E-03,
     +  0.39927859E-05, 0.96072399E-04, 0.64521264E-04,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      x_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23*x32        
      x_sp_cq1x   =x_sp_cq1x   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x23*x31        
c
      return
      end
      function t_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.2693314E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17809260E-03, 0.54506164E-01,-0.64578220E-01,-0.57909335E-03,
     + -0.19178996E-02, 0.13182421E-01,-0.10892597E-01, 0.41205592E-04,
     + -0.14258438E-01, 0.30719515E-01,-0.65919280E-01, 0.35225831E-01,
     +  0.13537859E-01, 0.62423467E-03,-0.34857620E-03,
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
c
c                  function
c
      t_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x23*x32        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)        *x31        
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x12                
     6  +coeff( 15)    *x23            
c
      return
      end
      function y_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.9509774E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.94092898E-02, 0.72639905E-01, 0.42631917E-01,
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
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
c
      return
      end
      function p_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.3708671E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37825757E-02, 0.22040663E-01, 0.26241813E-01,-0.33924304E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22*x31        
c
      return
      end
      function l_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.4928363E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.46331095E-03,-0.44748598E-04, 0.19777422E-03, 0.20204423E-03,
     +  0.59997019E-04,-0.52936526E-03,-0.11712846E-05, 0.31709464E-03,
     + -0.16305980E-02,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      l_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff(  9)        *x31*x41    
c
      return
      end
      function x_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.5099465E+01/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14918935E-02,-0.17886765E+00, 0.23501199E-02, 0.67369208E-01,
     +  0.56835240E+00, 0.92030159E-03, 0.73312074E-02,-0.19184297E-02,
     + -0.74613301E-04, 0.22445959E-02,-0.53937358E+00,-0.24904206E-03,
     + -0.23052975E-03, 0.10251175E-02,-0.72626104E-02, 0.15576654E-02,
     + -0.57882011E+00, 0.55055332E+00,
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
c
c                  function
c
      x_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)    *x23            
     7  +coeff(  7)        *x31        
     8  +coeff(  8)    *x21*x31        
      x_sp_cden   =x_sp_cden   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x23*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21*x33*x41    
     4  +coeff( 13)*x12    *x31*x41    
     5  +coeff( 14)*x11    *x33*x41    
     6  +coeff( 15)            *x41    
     7  +coeff( 16)    *x21    *x41    
     8  +coeff( 17)*x11    *x32        
      x_sp_cden   =x_sp_cden   
     9  +coeff( 18)*x11    *x31*x41    
c
      return
      end
      function t_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/  0.3619207E+01/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.57060584E-01, 0.11115807E+01, 0.16232762E+00, 0.13100660E-01,
     +  0.22418457E+03,-0.48660007E-02,-0.21151200E-01,-0.45349306E+00,
     + -0.44992459E+03, 0.22576875E+03, 0.33721852E-02,
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
c
c                  function
c
      t_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)            *x42    
     5  +coeff(  5)    *x23            
     6  +coeff(  6)        *x31        
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)*x11                
      t_sp_cden   =t_sp_cden   
     9  +coeff(  9)*x11*x22            
     1  +coeff( 10)*x12*x21            
     2  +coeff( 11)*x12    *x31        
c
      return
      end
      function y_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.6402886E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.68252799E-02, 0.26214095E-01, 0.65260425E-01, 0.34051086E-02,
     + -0.67562505E-03,-0.25416471E-02,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21            
     6  +coeff(  6)    *x22*x31        
c
      return
      end
      function p_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.5057726E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.44935271E-02,-0.41356571E-01,-0.22345006E-01,-0.49664001E-02,
     +  0.90289470E-02,-0.82189785E-02, 0.27884413E-02,-0.10955267E-01,
     +  0.25642231E+00, 0.51471783E-03,-0.91240596E-03,-0.23734906E+00,
     + -0.18472020E-01,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      p_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)*x11    *x31        
     6  +coeff(  6)    *x22*x31        
     7  +coeff(  7)    *x21            
     8  +coeff(  8)            *x41    
      p_sp_cden   =p_sp_cden   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)        *x33        
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x12                
c
      return
      end
      function l_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.1274115E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.66557672E-03, 0.37669790E+00,-0.11456131E-01,-0.28743206E-02,
     +  0.49938550E-02,-0.16124129E+00,-0.16903842E-02,-0.13751920E-02,
     +  0.58654465E-01,-0.39555635E-02,-0.57096299E-01,-0.30710944E-02,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      l_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)        *x32        
     5  +coeff(  5)        *x31        
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x23            
      l_sp_cden   =l_sp_cden   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)            *x41    
     2  +coeff( 11)*x11    *x31        
     3  +coeff( 12)    *x23*x31*x41    
c
      return
      end
      function x_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.2041022E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.78784917E-02, 0.11219431E+01,-0.20052146E-02,-0.84581785E-01,
     +  0.98890841E-01, 0.23659982E+00,-0.83049655E+00,-0.33730680E+00,
     +  0.77717349E-01, 0.17029047E+00,-0.16828825E+00,-0.69915196E-02,
     +  0.32376882E-01, 0.62726173E-02,-0.13161300E-01,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      x_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21*x32        
     5  +coeff(  5)        *x32        
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11                
     8  +coeff(  8)        *x31*x41    
      x_sp_cdex   =x_sp_cdex   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)        *x31        
     2  +coeff( 11)            *x41    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)*x11        *x41    
c
      return
      end
      function t_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/  0.5920142E+00/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18669374E-02,-0.94290026E-01,-0.25073748E-01, 0.27339888E-03,
     + -0.46703854E-03, 0.17130730E-03, 0.73219561E-02, 0.85299071E-04,
     +  0.19553013E-01, 0.64846448E-03,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      t_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21*x32        
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)*x11*x22            
      t_sp_cdex   =t_sp_cdex   
     9  +coeff(  9)*x11    *x32        
     1  +coeff( 10)    *x23*x31*x41    
c
      return
      end
      function y_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.7661504E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23389887E-02,-0.40522560E+00, 0.45254165E+00, 0.73188144E+00,
     + -0.72275406E+00,-0.16473278E-02,-0.23491209E-01,-0.23737223E+01,
     +  0.24025190E+01,-0.45110440E+00, 0.41174459E+00,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22            
      y_sp_cdex   =y_sp_cdex   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)*x11*x21*x31        
c
      return
      end
      function p_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/  0.8586791E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48617666E-03, 0.66417925E-01, 0.89432016E-01,-0.91495745E-01,
     + -0.67040078E-01,-0.78940988E+00,-0.68825674E-02, 0.23544803E+01,
     + -0.15727535E+01,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)*x11*x21            
      p_sp_cdex   =p_sp_cdex   
     9  +coeff(  9)*x12                
c
      return
      end
      function l_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.1112376E-01/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14712544E-01,-0.87934977E+00,-0.24040615E-01, 0.71987677E-02,
     +  0.42476697E-03, 0.75777536E-02, 0.45474282E+00,-0.48534572E-02,
     +  0.53284764E-02,-0.94849681E-02, 0.10285111E-02,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      l_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x21*x32        
     5  +coeff(  5)        *x31        
     6  +coeff(  6)        *x32        
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x21*x31        
      l_sp_cdex   =l_sp_cdex   
     9  +coeff(  9)    *x23*x32        
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)*x11*x23            
c
      return
      end
      function x_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/  0.4626027E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.70136744E-02, 0.48075646E+00, 0.14882736E-01,-0.21070852E-02,
     + -0.80543188E-02, 0.16786253E-01,-0.16027903E-01,-0.29755646E+00,
     + -0.20920938E+00,-0.13658800E-02, 0.43405764E-03, 0.21282898E+00,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      x_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)        *x32        
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)            *x41    
     8  +coeff(  8)*x11                
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)*x11    *x31        
c
      return
      end
      function t_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.1573276E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22588035E-02,-0.12591366E+00,-0.34892571E-02, 0.55900678E-01,
     +  0.67219522E-03,-0.43787709E-05,-0.34852729E-02, 0.29396103E-02,
     +  0.24240377E-03, 0.25353124E-03,-0.15364804E-03, 0.85540197E-03,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      t_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22*x32        
     2  +coeff( 11)*x11*x23            
     3  +coeff( 12)    *x23*x32        
c
      return
      end
      function y_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.3687025E-04/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20734449E-02,-0.53917319E+00, 0.58607036E+00, 0.10794742E+01,
     + -0.10689353E+01,-0.17365105E+00,-0.41910782E-01,-0.29791226E+01,
     +  0.18647638E+00, 0.30167410E+01,-0.54129869E+00, 0.49286652E+00,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22            
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)*x11*x21*x31        
c
      return
      end
      function p_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/  0.7061802E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13760751E-03, 0.62320329E-01, 0.10412099E+00,-0.10334565E+00,
     + -0.63082971E-01,-0.92696017E+00,-0.75714900E-02, 0.44478543E-03,
     +  0.27905939E+01,-0.18722508E+01,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      p_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)            *x42    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)*x12                
c
      return
      end
      function l_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.9169514E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11414282E-01,-0.56928802E+00,-0.20725615E-01,-0.34676322E-02,
     + -0.32953715E-02, 0.83824329E-03, 0.29221955E+00,-0.42465284E+00,
     +  0.10663365E-01, 0.41172245E+00, 0.17459477E+00,-0.16082302E+00,
     +  0.15702934E-02,-0.19296344E-02,
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
c
c                  function
c
      l_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)            *x42    
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x21*x31        
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11    *x31        
     2  +coeff( 11)    *x23*x32        
     3  +coeff( 12)*x11*x22*x32        
     4  +coeff( 13)*x11*x21*x33        
     5  +coeff( 14)*x11*x22*x31*x41    
c
      return
      end
      function x_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.2200538E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37452008E-02, 0.56259919E-01, 0.50085008E-01,-0.49172521E-01,
     +  0.11443295E-01,-0.34791425E-01, 0.18718439E+00, 0.43332605E+01,
     + -0.33167205E-02, 0.31420749E-01,-0.19360600E+00,-0.74380279E+01,
     +  0.31035957E+01,-0.45802845E-02,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
c
c                  function
c
      x_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x23            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff(  9)    *x24            
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)*x12*x21            
     5  +coeff( 14)*x11*x24*x32        
c
      return
      end
      function t_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.1370570E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.65654840E-05,-0.13505100E-01, 0.30637514E-02,-0.27517893E-02,
     +  0.27893225E-02,-0.10957087E-02,-0.61974555E-01, 0.70503249E-03,
     + -0.93986094E-03, 0.12884454E+00,-0.61226422E-02, 0.49043581E-03,
     + -0.68510234E-01,-0.51088475E-04,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x23            
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff(  9)*x11*x23            
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x22*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22*x31        
c
      return
      end
      function y_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1453297E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.90209731E-04,-0.60699306E-01, 0.77294886E-01, 0.23544034E-01,
     +  0.24939182E+00,-0.83278853E+00, 0.58342302E+00,-0.19092897E-01,
     + -0.26427830E-01, 0.51895612E+00,-0.89159781E+00, 0.35925746E+00,
     +  0.10673750E-01,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      y_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)*x11                
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)        *x33        
     2  +coeff( 11)        *x32*x41    
     3  +coeff( 12)        *x31*x42    
     4  +coeff( 13)            *x43    
c
      return
      end
      function p_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.2409164E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72516123E-03,-0.19030863E+00, 0.16849877E+00, 0.53301327E+02,
     +  0.35434432E-01,-0.19444218E-01,-0.10813697E+03, 0.54850552E+02,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x12                
c
      return
      end
      function l_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1004987E-01/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12359481E-01,-0.52740657E+00,-0.24916641E-01,-0.28650069E-02,
     +  0.77140965E-02, 0.92173449E-03, 0.24989888E+00, 0.18292344E+00,
     + -0.18577358E+00, 0.72899141E-03,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
c
c                  function
c
      l_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)        *x32        
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x21*x31        
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)*x11*x23            
c
      return
      end
      function x_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.1805452E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38226414E-02, 0.94131706E-02,-0.79365820E-02, 0.75835944E-03,
     +  0.17635282E+01, 0.16149877E+00, 0.42934712E-01,-0.60066748E-02,
     + -0.36104145E+01,-0.35013416E+00,-0.20586331E+00,-0.88766376E-02,
     +  0.18721707E+01, 0.35551256E+00, 0.16490655E-02,-0.89169452E-02,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      x_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x32        
      x_sp_cfp    =x_sp_cfp    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)*x11    *x31        
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x12                
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)    *x24            
     7  +coeff( 16)*x11*x23            
c
      return
      end
      function t_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.1370436E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.65458180E-05,-0.13504932E-01, 0.30623144E-02,-0.27503453E-02,
     +  0.27895388E-02,-0.10957246E-02,-0.61965559E-01, 0.70495351E-03,
     + -0.93996112E-03, 0.12882711E+00,-0.61224401E-02, 0.49039250E-03,
     + -0.68501852E-01,-0.51096893E-04,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      t_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x23            
      t_sp_cfp    =t_sp_cfp    
     9  +coeff(  9)*x11*x23            
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x22*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22*x31        
c
      return
      end
      function y_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/  0.2146685E-02/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29043653E-02,-0.11591473E+00, 0.70562474E-01,-0.76317797E+01,
     +  0.75994248E+01,-0.10054087E+02, 0.76521115E+01,-0.76080894E+01,
     +  0.18942253E+02,-0.88784399E+01,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      y_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)*x11        *x41    
      y_sp_cfp    =y_sp_cfp    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)*x12                
c
      return
      end
      function p_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.2409057E-03/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72515634E-03,-0.19031121E+00, 0.16850124E+00, 0.53302105E+02,
     +  0.35434946E-01,-0.19444562E-01,-0.10813855E+03, 0.54851353E+02,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x12                
c
      return
      end
      function l_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.1033247E-01/
      data xmin/
     1 -0.11988E+00,-0.46782E-01,-0.59070E-01,-0.23971E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11750E+00, 0.45342E-01, 0.40924E-01, 0.18288E-01, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12567790E-01,-0.27802187E+00,-0.25549928E-01,-0.36515025E-02,
     +  0.85841388E+00, 0.14103093E-02,-0.27307551E-02, 0.13230976E-02,
     + -0.17294620E+01, 0.87953430E+00,-0.54787681E-01, 0.52984070E-01,
     +  0.21854371E-02,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
     4  +coeff(  4)        *x32        
     5  +coeff(  5)    *x21*x32        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)*x11    *x32        
     8  +coeff(  8)    *x21*x31        
      l_sp_cfp    =l_sp_cfp    
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23*x31        
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x23*x32        
c
      return
      end
