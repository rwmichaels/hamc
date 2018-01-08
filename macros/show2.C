{

// Show some histograms

TCanvas *c1 = new TCanvas;

c1->Divide(2,2);

c1->cd(1);
th1->Draw();

c1->cd(2);
th2->Draw();

c1->cd(3);
th3->Draw();

c1->cd(4);
th4->Draw();


}
