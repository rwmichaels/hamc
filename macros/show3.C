{

// Show some ntuple data
// Ntuple:
hamc->Print();

TCanvas *c3 = new TCanvas;

c3->Divide(2,2);

c3->cd(1);
hamc->Draw("pmom");

c3->cd(2);
hamc->Draw("yfoc:xfoc","inaccept==1");

c3->cd(3);
hamc->Draw("crsec:theta","abs((180*theta/3.14)-5)<3.5");

c3->cd(4);
hamc->Draw("zscat");


}
