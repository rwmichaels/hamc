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

   xysepia->SetMarkerColor(2);
   xysepoa->SetMarkerColor(2);

   xyq1a->SetMarkerColor(2);
   xycoll2a->SetMarkerColor(2);

   c1->cd(1);

   xysepi->Draw();
   xysepia->Draw("same");

   c1->cd(2);

   xysepo->Draw();
   xysepoa->Draw("same");

   c1->cd(3);

   xyq1->Draw();
   xyq1a->Draw("same");

   c1->cd(4);

   xycoll2->Draw();
   xycoll2a->Draw("same");

}
