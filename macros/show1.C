{

// Show some histograms
// Most of these are defined in hamcTrackOut::Init at the moment.

TCanvas *c1 = new TCanvas;

c1->Divide(2,2);

c1->cd(1);
xyfoc1->Draw();

c1->cd(2);
xyfoc2->Draw();

c1->cd(3);
xyfoc2->Draw("LEGO");


}
