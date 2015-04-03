{

   gROOT->Reset();

   gStyle->SetOptStat(0);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.07);
   gStyle->SetTitleW(0.6);
   gStyle->SetLabelSize(0.052,"x");
   gStyle->SetLabelSize(0.052,"y");
   gROOT->ForceStyle();


    Float_t xscale = 1.0;
    FILE *fd1, *fd2;
    char schar[200],scout[200];
    fd1=fopen("accept4.0.dat", "r");
    fd2=fopen("accept4.5.dat", "r");
    fd3=fopen("accept5.0.dat", "r");

    if (fd1 == NULL) {
      cout << "Cannot find accept4.0.dat";
      exit(0);
    }
    if (fd2 == NULL) {
      cout << "Cannot find accept4.5.dat";
      exit(0);
    }
    if (fd3 == NULL) {
      cout << "Cannot find accept5.0.dat";
      exit(0);
    }

    TCanvas c1;
    TH2F *h2 = new TH2F("hacc","Acceptance Function",200,2.5,8,100,0,1.0);
    
    Float_t theta_degr, xacc;
    Int_t npts; 
    Int_t ipt=0;

    h2->Draw();

    Int_t npts1=0;
    while(fgets(schar,100,fd1)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
        npts1++;
    }
    Int_t npts2=0;
    while(fgets(schar,100,fd2)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
        npts2++;
    }
    Int_t npts3=0;
    while(fgets(schar,100,fd3)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
        npts3++;
    }
    fclose(fd1);
    fclose(fd2);
    fclose(fd3);

    npts = npts1;
    if (npts2 > npts) npts = npts2;
    if (npts3 > npts) npts = npts3;

    cout << "num points"<<npts<<endl;

    Float_t th[npts],acc[npts];

    fd1=fopen("accept4.0.dat", "r");
    while(fgets(schar,100,fd1)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
	//        cout << "theta "<<theta_degr<<"   acc "<<xacc<<endl;

        th[ipt] = theta_degr;
        acc[ipt] = xscale*xacc;

	//        if (th[ipt]>4.5&&th[ipt]<4.8) cout << "big ? "<<theta_degr<<"   acc "<<xacc<<endl;
        ipt++;

    }
    fclose(fd1);

    cout << "num points = "<<npts<<endl;

    gr1 = new TGraph(npts, th, acc);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.2);
    gr1->SetMarkerStyle(3);
    gr1->Draw("P");

    // ------------------

    fd2=fopen("accept4.5.dat", "r");
    ipt=0;
    while(fgets(schar,100,fd2)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
	//        cout << "theta "<<theta_degr<<"   acc "<<xacc<<endl;

        th[ipt] = theta_degr;
        acc[ipt] = xscale*xacc;

	//        if (th[ipt]>4.5&&th[ipt]<4.8) cout << "big ? "<<theta_degr<<"   acc "<<xacc<<endl;
        ipt++;

    }
    fclose(fd2);

    cout << "num points = "<<npts<<endl;

    gr2 = new TGraph(npts, th, acc);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerSize(1.2);
    gr2->SetMarkerStyle(3);
    gr2->Draw("P");

    // ------------------

    fd3=fopen("accept5.0.dat", "r");
    ipt=0;
    while(fgets(schar,100,fd3)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
	//        cout << "theta "<<theta_degr<<"   acc "<<xacc<<endl;

        th[ipt] = theta_degr;
        acc[ipt] = xscale*xacc;

	//        if (th[ipt]>4.5&&th[ipt]<4.8) cout << "big ? "<<theta_degr<<"   acc "<<xacc<<endl;
        ipt++;

    }
    fclose(fd3);

    cout << "num points = "<<npts<<endl;

    gr3 = new TGraph(npts, th, acc);
    gr3->SetMarkerColor(8);
    gr3->SetMarkerSize(1.2);
    gr3->SetMarkerStyle(3);
    gr3->Draw("P");

    TLegend *leg = new TLegend(0.55,0.70, 0.89, 0.89);
    leg->AddEntry(gr1,"4.0 deg, FWHM = 1.40 deg","p");
    leg->AddEntry(gr2,"4.5 deg, FWHM = 1.56 deg","p");
    leg->AddEntry(gr3,"5.0 deg, FWHM = 1.70 deg","p");
    leg->Draw();



    c1->Update();


}

