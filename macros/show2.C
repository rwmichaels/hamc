{

// Show some histograms

TCanvas *c2 = new TCanvas;

c2->Divide(2,2);

c2->cd(1);
th1->Draw();

c2->cd(2);
th2->Draw();

c2->cd(3);
th3->Draw();

c2->cd(4);
th4->Draw();


}
