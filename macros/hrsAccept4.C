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

   TCanvas *c4a = new TCanvas();
   c4a->Divide(2,2);

   thpha2->SetMarkerColor(2);

   c4a->cd(1);

   thph2->Draw();
   thpha2->Draw("same");

   c4a->cd(2);

   thpha2->Draw("");

   c4a->cd(3);


}
