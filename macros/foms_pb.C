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

   Int_t which = 1; // 1 = Rate vs E;  2 = A vs E;  3 = daa vs E;  4 = drr vs E;  5 = drrtot vs E;  6 = 4&5; 7 = sensi

   Int_t logy   = 1;   // 0 or 1  (to use Log scale for Y)

   Float_t energypb[npt], qsqpb[npt], asypb[npt], ratepb[npt], sensipb[npt];
   Float_t daapb[npt], drrpb[npt], drrtotpb[npt];

    c1 = new TCanvas("c1","",200,10,700,500);

//    c1->SetFillColor(42);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    if (logy) c1->SetLogy(1);

    TH2F *h1a = new TH2F("h1a","Pb208: Rate vs Energy ",200,0.5,1.3,200,1e8,5e10);
    TH2F *h2a = new TH2F("h2a","Pb208: Asy vs Energy ",200,0.5,1.3,200,0,1.0);
    TH2F *h3a = new TH2F("h3a","Pb208: Stat Err in A vs Energy ",200,0.5,1.3,200,0.01,0.05);
    TH2F *h4a = new TH2F("h4a","Pb208: Stat Error in R_n vs Energy ",200,0.5,1.3,200,0.002,0.05);
    TH2F *h5a = new TH2F("h5a","Pb208: Tot Error in R_n vs Energy ",200,0.5,1.3,200,0.006,0.027);
    TH2F *h7a = new TH2F("h7a","Pb208: Sensitivity to R_n vs Energy ",200,0.5,1.3,200,0.0,0.04);

    if (which==1) h1a->Draw();
    if (which==2) h2a->Draw();
    if (which==3) h3a->Draw();
    if (which==4) h4a->Draw();
    if (which==5 || which ==6) h5a->Draw();
    if (which==7) h7a->Draw();

    FILE *fd = fopen("foms_pb.txt","r");
    if (fd == NULL) {
      cout << "No input data !.  Quitting. "<<endl;
      exit(0);
   }

   i = 0;
   while (fgets(strin,100,fd) != NULL)  {

     sscanf(strin,"FOMS %f %f %f %f %f %f %f %f",&ee,&qq,&aa,&rr,&ss,&ddaa,&ddrr,&ddrrtt);

     energypb[i] = ee;
     qsqpb[i] = qq;
     asypb[i] = aa;
     ratepb[i] = rr;
     sensipb[i] = -1.0*ss;
     daapb[i] = ddaa;
     drrpb[i] = -1*ddrr;
     drrtotpb[i] = -1*ddrrtt;

     cout << "\n\n ----------- \n"<<endl;

     cout << "Pb208  "<<i+1<<"  E = "<<energypb[i]<<"   qsq = "<<qsqpb[i]<<"   asy = "<<asypb[i]<<" ppm";
     cout << "  rate = "<<ratepb[i]<<endl;
     cout << "  Sensi = "<<sensipb[i]<<"   daa = "<<daapb[i]<<"   drr = "<<drrpb[i]<<"   drrtot = "<<drrtotpb[i]<<endl;

     i++;
     if (i > npt) break;
   
   }

   if (which == 1) gr1 = new TGraph(npt, energypb, ratepb);
   if (which == 2) gr1 = new TGraph(npt, energypb, asypb);
   if (which == 3) gr1 = new TGraph(npt, energypb, daapb);
   if (which == 4) gr1 = new TGraph(npt, energypb, drrpb);
   if (which == 5) gr1 = new TGraph(npt, energypb, drrtotpb);
   if (which == 6) gr1 = new TGraph(npt, energypb, drrpb);
   if (which == 7) gr1 = new TGraph(npt, energypb, sensipb);

    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.2);
    gr1->SetMarkerStyle(3);
    gr1->Draw("P");

    if (which == 6) {
       gr2 = new TGraph(npt, energypb, drrtotpb);
       gr2->SetMarkerColor(3);
       gr2->SetMarkerSize(1.2);
       gr2->SetMarkerStyle(3);
       gr2->Draw("P");
       TLegend *leg = new TLegend(0.55,0.7,0.89,0.89);
       leg->AddEntry(gr1,"Stat Error in R_n","p");
       leg->AddEntry(gr2,"Total Err inclu P_beam (1%)","p");
       leg->Draw();
    }


}
