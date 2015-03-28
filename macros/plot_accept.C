{

    Float_t xscale = 1.0;
    FILE *fd1, *fd2;
    char schar[200],scout[200];
    fd1=fopen("accept.dat", "r");

    if (fd1 == NULL) {
      cout << "Cannot find accept.dat";
      exit(0);
    }

    TCanvas c1;
    TH2F *h2 = new TH2F("hacc","Acceptance Function",200,2.5,8,100,0,1.0);
    
    Float_t theta_degr, xacc; 
    Int_t ipt=0;

    h2->Draw();

    Int_t npts=0;
    while(fgets(schar,100,fd1)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
        npts++;
    }

    Float_t th[npts],acc[npts];

    fclose(fd1);
    fd1=fopen("accept.dat", "r");
    while(fgets(schar,100,fd1)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
	//        cout << "theta "<<theta_degr<<"   acc "<<xacc<<endl;

        th[ipt] = theta_degr;
        acc[ipt] = xscale*xacc;

	//        if (th[ipt]>4.5&&th[ipt]<4.8) cout << "big ? "<<theta_degr<<"   acc "<<xacc<<endl;
        ipt++;

    }

    cout << "num points = "<<npts<<endl;

    gr1 = new TGraph(npts, th, acc);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.2);
    gr1->SetMarkerStyle(3);
    gr1->Draw("P");

    c1->Update();

    fclose(fd1);

}

