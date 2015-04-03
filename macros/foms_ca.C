{

  // Plot of figure of merit info for quick_fom study of hamc.
  // need to control PREX_model in prex.dat
  // Run hamc with quick_fom=1 in hamcPhyPREX.h, then do this:
  // ./prex 20 > x.x ;  grep FOMS > foms_pb.tx ;  root .x foms_pb.C

   gROOT->Reset();

   gStyle->SetOptStat(0);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.07);
   gStyle->SetTitleW(0.7);
   gStyle->SetLabelSize(0.04,"x");
   gStyle->SetLabelSize(0.04,"y");
   gROOT->ForceStyle();

   Float_t ee,qq,aa,rr,ss,ddaa,ddrr,ddrrtt;
   char strin[100];

   Int_t npt = 61;

   Int_t i;

   Int_t which = 6; // 1 = Rate vs E;  2 = A vs E;  3 = daa vs E;  4 = drr vs E;  5 = drrtot vs E;  6 = 4&5; 7 = sensi

   Int_t logy   = 0;   // 0 or 1  (to use Log scale for Y)

   Float_t energy40[npt], qsq40[npt], asy40[npt], rate40[npt], sensi40[npt];
   Float_t daa40[npt], drr40[npt], drrtot40[npt];

   Float_t energy48[npt], qsq48[npt], asy48[npt], rate48[npt], sensi48[npt];
   Float_t daa48[npt], drr48[npt], drrtot48[npt];

    c1 = new TCanvas("c1","",200,10,700,500);

//    c1->SetFillColor(42);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    if (logy) c1->SetLogy(1);

    TH2F *h1b = new TH2F("h1b","Ca48: Rate vs Energy ",200,0.5,3.5,200,1e5,1e10);
    TH2F *h2b = new TH2F("h2b","Ca48: Asy vs Energy ",200,0.5,3.5,200,0,12);
    TH2F *h3b = new TH2F("h3b","Ca48: Stat Err in A vs Energy ",200,0.5,3.5,200,0.002,0.08);
    TH2F *h4b = new TH2F("h4b","Ca48: Stat Error in R_n vs Energy ",200,0.5,3.5,200,0.002,0.05);
    TH2F *h5b = new TH2F("h5b","Ca48: Tot Error in R_n vs Energy ",200,0.5,3.5,200,0.004,0.025);
    TH2F *h7b = new TH2F("h7b","Ca48: Sensitivity to R_n vs Energy ",200,0.5,3.5,200,-0.02,0.085);


    if (which==1) h1b->Draw();
    if (which==2) h2b->Draw();
    if (which==3) h3b->Draw();
    if (which==4) h4b->Draw();
    if (which==5 || which ==6) h5b->Draw();
    if (which==7) h7b->Draw();

    FILE *fd = fopen("foms_ca.txt","r");
    if (fd == NULL) {
      cout << "No input data !.  Quitting. "<<endl;
      exit(0);
   }

   i = 0;
   while (fgets(strin,100,fd) != NULL)  {

     sscanf(strin,"FOMS %f %f %f %f %f %f %f %f",&ee,&qq,&aa,&rr,&ss,&ddaa,&ddrr,&ddrrtt);

     energy48[i] = ee;
     qsq48[i] = qq;
     asy48[i] = aa;
     rate48[i] = rr;
     sensi48[i] = -1.0*ss;
     daa48[i] = ddaa;
     drr48[i] = ddrr;
     drrtot48[i] = ddrrtt;

     cout << "\n\n ----------- \n"<<endl;

     cout << "Ca48  "<<i+1<<"  E = "<<energy48[i]<<"   qsq = "<<qsq48[i]<<"   asy = "<<asy48[i]<<" ppm";
     cout << "  rate = "<<rate48[i]<<endl;
     cout << "  Sensi = "<<sensi48[i]<<"   daa = "<<daa48[i]<<"   drr = "<<drr48[i]<<"   drrtot = "<<drrtot48[i]<<endl;

     i++;
     if (i > npt) break;
   
   }

   if (which == 1) gr3 = new TGraph(npt, energy48, rate48);
   if (which == 2) gr3 = new TGraph(npt, energy48, asy48);
   if (which == 3) gr3 = new TGraph(npt, energy48, daa48);
   if (which == 4) gr3 = new TGraph(npt, energy48, drr48);
   if (which == 5) gr3 = new TGraph(npt, energy48, drrtot48);
   if (which == 6) gr3 = new TGraph(npt, energy48, drr48);
   if (which == 7) gr3 = new TGraph(npt, energy48, sensi48);

    gr3->SetMarkerColor(4);
    gr3->SetMarkerSize(1.2);
    gr3->SetMarkerStyle(3);
    gr3->Draw("P");

    if (which == 6) {
       gr4 = new TGraph(npt, energy48, drrtot48);
       gr4->SetMarkerColor(2);
       gr4->SetMarkerSize(1.2);
       gr4->SetMarkerStyle(3);
       gr4->Draw("P");
       TLegend *leg = new TLegend(0.55,0.7,0.89,0.89);
       leg->AddEntry(gr3,"Stat Error in R_n","p");
       leg->AddEntry(gr4,"Total Err inclu P_beam (1%)","p");
       leg->Draw();
    }

   fclose(fd);


}













}
