{

  // comparing angles for Ca and for Pb

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

   Int_t which = 5; // 1 = Rate vs E;  2 = A vs E;  3 = daa vs E;  4 = drr vs E;  5 = drrtot vs E;  6 = 4&5; 7 = sensi

   Int_t logy   = 0;   // 0 or 1  (to use Log scale for Y)

   Float_t energypb1[npt], qsqpb1[npt], asypb1[npt], ratepb1[npt], sensipb1[npt];
   Float_t daapb1[npt], drrpb1[npt], drrtotpb1[npt];
   Float_t energypb2[npt], qsqpb2[npt], asypb2[npt], ratepb2[npt], sensipb2[npt];
   Float_t daapb2[npt], drrpb2[npt], drrtotpb2[npt];
   Float_t energypb3[npt], qsqpb3[npt], asypb3[npt], ratepb3[npt], sensipb3[npt];
   Float_t daapb3[npt], drrpb3[npt], drrtotpb3[npt];


    c1 = new TCanvas("c1","",200,10,700,500);

//    c1->SetFillColor(42);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    if (logy) c1->SetLogy(1);

    Float_t elo = 0.7;
    Float_t ehi = 1.3;

    TH2F *h1a = new TH2F("h1a"," Rate vs Energy ",200,elo,ehi,200,1e8,5e10);
    TH2F *h2a = new TH2F("h2a"," Asy vs Energy ",200,elo,ehi,200,0,1.0);
    TH2F *h3a = new TH2F("h3a"," Stat Err in A vs Energy ",200,elo,ehi,200,0.01,0.05);
    TH2F *h4a = new TH2F("h4a"," Stat Error in R_n vs Energy ",200,elo,ehi,200,0.002,0.05);
    TH2F *h5a = new TH2F("h5a"," Tot Error in R_n vs Energy ",200,elo,ehi,200,0.004,0.027);
    TH2F *h7a = new TH2F("h7a"," Sensitivity to R_n vs Energy ",200,elo,ehi,200,0.0,0.04);

    if (which==1) h1a->Draw();
    if (which==2) h2a->Draw();
    if (which==3) h3a->Draw();
    if (which==4) h4a->Draw();
    if (which==5 || which ==6) h5a->Draw();
    if (which==7) h7a->Draw();

    // 1st of 3 files

    FILE *fd = fopen("foms1.txt","r");
    if (fd == NULL) {
      cout << "No input data !.  Quitting. "<<endl;
      exit(0);
   }

   i = 0;
   while (fgets(strin,100,fd) != NULL)  {

     sscanf(strin,"FOMS %f %f %f %f %f %f %f %f",&ee,&qq,&aa,&rr,&ss,&ddaa,&ddrr,&ddrrtt);

     energypb1[i] = ee;
     qsqpb1[i] = qq;
     asypb1[i] = aa;
     ratepb1[i] = rr;
     sensipb1[i] = -1.0*ss;
     daapb1[i] = ddaa;
     drrpb1[i] = ddrr;
     drrtotpb1[i] = ddrrtt;

     i++;
     if (i > npt) break;
   
   }

   if (which == 1) gr1 = new TGraph(npt, energypb1, ratepb1);
   if (which == 2) gr1 = new TGraph(npt, energypb1, asypb1);
   if (which == 3) gr1 = new TGraph(npt, energypb1, daapb1);
   if (which == 4) gr1 = new TGraph(npt, energypb1, drrpb1);
   if (which == 5) gr1 = new TGraph(npt, energypb1, drrtotpb1);
   if (which == 6) gr1 = new TGraph(npt, energypb1, drrpb1);
   if (which == 7) gr1 = new TGraph(npt, energypb1, sensipb1);

    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.2);
    gr1->SetMarkerStyle(3);
    gr1->Draw("P");

    if (which == 6) {
       gr1a = new TGraph(npt, energypb1, drrtotpb1);
       gr1a->SetMarkerColor(2);
       gr1a->SetMarkerSize(1.2);
       gr1a->SetMarkerStyle(3);
       gr1a->Draw("P");
       TLegend *leg = new TLegend(0.55,0.7,0.89,0.89);
       leg->AddEntry(gr1,"Stat Error in R_n","p");
       leg->AddEntry(gr1a,"Total Err inclu P_beam (1%)","p");
       leg->Draw();
    }

    // 2nd file

    fclose(fd);
    FILE *fd = fopen("foms2.txt","r");
    if (fd == NULL) {
      cout << "No input data !.  Quitting. "<<endl;
      exit(0);
   }

   i = 0;
   while (fgets(strin,100,fd) != NULL)  {

     sscanf(strin,"FOMS %f %f %f %f %f %f %f %f",&ee,&qq,&aa,&rr,&ss,&ddaa,&ddrr,&ddrrtt);

     energypb2[i] = ee;
     qsqpb2[i] = qq;
     asypb2[i] = aa;
     ratepb2[i] = rr;
     sensipb2[i] = -1.0*ss;
     daapb2[i] = ddaa;
     drrpb2[i] = ddrr;
     drrtotpb2[i] = ddrrtt;

     i++;
     if (i > npt) break;
   
   }

   if (which == 1) gr2 = new TGraph(npt, energypb2, ratepb2);
   if (which == 2) gr2 = new TGraph(npt, energypb2, asypb2);
   if (which == 3) gr2 = new TGraph(npt, energypb2, daapb2);
   if (which == 4) gr2 = new TGraph(npt, energypb2, drrpb2);
   if (which == 5) gr2 = new TGraph(npt, energypb2, drrtotpb2);
   if (which == 6) gr2 = new TGraph(npt, energypb2, drrpb2);
   if (which == 7) gr2 = new TGraph(npt, energypb2, sensipb2);

    gr2->SetMarkerColor(2);
    gr2->SetMarkerSize(1.2);
    gr2->SetMarkerStyle(3);
    gr2->Draw("P");

    if (which == 6) {
       gr2a = new TGraph(npt, energypb2, drrtotpb2);
       gr2a->SetMarkerColor(4);
       gr2a->SetMarkerSize(1.2);
       gr2a->SetMarkerStyle(3);
       gr2a->Draw("P");
       TLegend *leg = new TLegend(0.55,0.7,0.89,0.89);
       leg->AddEntry(gr2,"Stat Error in R_n","p");
       leg->AddEntry(gr2a,"Total Err inclu P_beam (1%)","p");
       leg->Draw();
    }


    // 3rd file

    FILE *fd = fopen("foms3.txt","r");
    if (fd == NULL) {
      cout << "No input data !.  Quitting. "<<endl;
      exit(0);
   }

   i = 0;
   while (fgets(strin,100,fd) != NULL)  {

     sscanf(strin,"FOMS %f %f %f %f %f %f %f %f",&ee,&qq,&aa,&rr,&ss,&ddaa,&ddrr,&ddrrtt);

     energypb3[i] = ee;
     qsqpb3[i] = qq;
     asypb3[i] = aa;
     ratepb3[i] = rr;
     sensipb3[i] = -1.0*ss;
     daapb3[i] = ddaa;
     drrpb3[i] = ddrr;
     drrtotpb3[i] = ddrrtt;

     i++;
     if (i > npt) break;
   
   }

   if (which == 1) gr3 = new TGraph(npt, energypb3, ratepb3);
   if (which == 2) gr3 = new TGraph(npt, energypb3, asypb3);
   if (which == 3) gr3 = new TGraph(npt, energypb3, daapb3);
   if (which == 4) gr3 = new TGraph(npt, energypb3, drrpb3);
   if (which == 5) gr3 = new TGraph(npt, energypb3, drrtotpb3);
   if (which == 6) gr3 = new TGraph(npt, energypb3, drrpb3);
   if (which == 7) gr3 = new TGraph(npt, energypb3, sensipb3);

    gr3->SetMarkerColor(8);
    gr3->SetMarkerSize(1.2);
    gr3->SetMarkerStyle(3);
    gr3->Draw("P");

    if (which == 6) {
       gr3a = new TGraph(npt, energypb3, drrtotpb3);
       gr3a->SetMarkerColor(8);
       gr3a->SetMarkerSize(1.2);
       gr3a->SetMarkerStyle(3);
       gr3a->Draw("P");
       TLegend *leg = new TLegend(0.55,0.7,0.89,0.89);
       leg->AddEntry(gr3,"Stat Error in R_n","p");
       leg->AddEntry(gr3a,"Total Err inclu P_beam (1%)","p");
       leg->Draw();
    }


    if (which == 5) {
       TLegend *leg = new TLegend(0.55,0.7,0.89,0.89);
       leg->AddEntry(gr1,"theta = 4.0","p");
       leg->AddEntry(gr2,"theta = 4.5","p");
       leg->AddEntry(gr3,"theta = 5.0","p");
       leg->Draw();
    }        

}
