{

    Int_t MAXPOINT=291;
    FILE *fd1, *fd2;
    char schar[200],scout[200];
    fd1=fopen("accept.dat", "r");

    if (fd1 == NULL) {
      cout << "Cannot find accept.dat";
      exit(0);
    }

    TCanvas c1;
    TH2F *h2 = new TH2F("hacc","Acceptance Fcn",200,2.5,8,100,0,1.0);
    
    Float_t theta_degr, xacc; 
    Float_t th[MAXPOINT],acc[MAXPOINT];
    Int_t ipt=0;

    h2->Draw();

    while(fgets(schar,100,fd1)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);
        cout << "theta "<<theta_degr<<"   acc "<<xacc<<endl;

        th[ipt] = theta_degr;
        acc[ipt] = 100.*xacc;
        ipt++;
        if (ipt > MAXPOINT) {
	  cout << "ERROR in ipt"<<endl;
          exit(0);
        }


    }

    cout << "num points = "<<ipt<<endl;

    gr1 = new TGraph(MAXPOINT, th, acc);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.2);
    gr1->SetMarkerStyle(3);
    gr1->Draw("P");

    c1->Update();

    fclose(fd1);

}

