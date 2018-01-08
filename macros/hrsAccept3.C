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

   TCanvas *c1 = new TCanvas();
   c1->Divide(2,2);

   xyfoc1a->SetMarkerColor(2);

   c1->cd(1);

   xyfoc1->Draw();
   xyfoc1a->Draw("same");

   c1->cd(2);

   xyfoc2->Draw();

   c1->cd(3);

   xyfoc2->Draw("LEGO");


}
