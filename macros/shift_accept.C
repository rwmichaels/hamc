{

  // Shift and reduce acceptance in a smooth way.

    Int_t toshift = 0;  // = 1 if shifting to new_central
 
    Float_t theta_central = 5;
    Float_t new_central = 4;

    Int_t MAXPOINT=291;
    FILE *fd1, *fd2;
    char schar[200],scout[200];
    fd1=fopen("accept.dat", "r");
    fd2=fopen("accout.dat", "w");

    if (fd1 == NULL) {
      cout << "Cannot find accept.dat";
      exit(0);
    }

    TCanvas c1;
    TH2F *h2 = new TH2F("hacc","Acceptance Function",200,2.5,8,100,0,1.0);
    
    Float_t theta_degr, xacc; 
    Float_t th[2*MAXPOINT],acc[2*MAXPOINT];
    Float_t thshift[2*MAXPOINT], accshift[2*MAXPOINT];
    Int_t ipt=0;

    h2->Draw();

    while(fgets(schar,100,fd1)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
	//        cout << "theta "<<theta_degr<<"   acc "<<xacc<<endl;
        th[ipt] = theta_degr;
        acc[ipt] = 100.*xacc;
        ipt++;
        if (ipt > MAXPOINT) {
	  cout << "ERROR in ipt"<<endl;
          exit(0);
        }
    }

    // Analysis: find the peak, and the 1/2 peak.  This determines the
    // width.  Then reduce this width by tan(5)/tan(theta_central), and
    // take that fraction of bins from about the central max.

    Float_t accmax = -999;
    Int_t imax=0;
    for (Int_t i=0; i < MAXPOINT; i++) {
      if (acc[i] > accmax) {
        accmax = acc[i];
        imax = i;
      }
    }
    cout << "Max accept "<<accmax<<"   "<<imax<<endl;

    Float_t thetaL = 0;
    Float_t thetaR = 0;
    Float_t minleft = 999;
    Float_t minright = 999;
    Float_t diff;
    for (Int_t i=0; i < MAXPOINT; i++) {
      diff = acc[i]-0.5*accmax;
      if (diff < 0)  diff = -1.0*diff;
      if (th[i] < theta_central && diff < minleft) {
        thetaL = th[i];
        minleft = diff;
      }
      if (th[i] > theta_central && diff < minright) {
        thetaR = th[i];
        minright = diff;
      }
    }
    cout << "thetaL = "<<thetaL<<"   thetaR = "<<thetaR<<endl;

    Int_t ncentral=0;
    for (Int_t i=0; i < MAXPOINT; i++) {
      if (th[i] > thetaL && th[i] < thetaR) ncentral++;
    }
    cout << "number in central region  "<<ncentral<<endl;
    
    Float_t pi = 3.1415926;
    Float_t xrfact = 1.-(TMath::Tan(new_central*pi/180.) / TMath::Tan(5.0*pi/180.));
    cout << "xrfact = "<<xrfact<<endl;
    Float_t xreduce = 0.5*xrfact*ncentral;
    Int_t nreduce = (Int_t) xreduce;
    cout << "number to reduce by "<<xrfact<< "  "<<xreduce<<"  "<<nreduce<<endl;

    Int_t ipt=0;
    Int_t icenter=0;
    for (Int_t i=0; i < MAXPOINT; i++) {
      if ((i > imax-nreduce) && (i < imax+nreduce)) {
        icenter = ipt;
        continue;
      }
      thshift[ipt] = th[i];
      accshift[ipt] = acc[i];
      ipt++;      
    }
    Int_t nshpts = ipt;
    cout << "Num points in new acc fcn "<<nshpts;
    cout << "center "<< icenter<<endl;

    for (Int_t i=0; i < nshpts; i++) {
      thshift[i] = thshift[i] - (theta_central-new_central);
      if (i > icenter) {
        thshift[i] = thshift[i] - (th[imax+nreduce]-th[imax-nreduce]);
      } 
    }
 
    for (Int_t i=0; i < nshpts; i++) {
      theta_degr = thshift[i];
      xacc = accshift[i];
      sprintf(scout, "%f %f \n", theta_degr, xacc);
      fputs(scout,fd2);
    }

    gr1 = new TGraph(MAXPOINT, th, acc);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.2);
    gr1->SetMarkerStyle(3);
    gr1->Draw("P");

    gr2 = new TGraph(nshpts, thshift, accshift);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerSize(1.2);
    gr2->SetMarkerStyle(3);
    gr2->Draw("P");

    c1->Update();


    fclose(fd1);
    fclose(fd2);

}

