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

   xydipia->SetMarkerColor(2);
   xydipoa->SetMarkerColor(2);

   xyq3ia->SetMarkerColor(2);
   xyq3oa->SetMarkerColor(2);

   c1->cd(1);

   xydipi->Draw();
   xydipia->Draw("same");

   c1->cd(2);

   xydipo->Draw();
   xydipoa->Draw("same");

   c1->cd(3);

   xyq3i->Draw();
   xyq3ia->Draw("same");

   c1->cd(4);

   xyq3o->Draw();
   xyq3oa->Draw("same");

}
