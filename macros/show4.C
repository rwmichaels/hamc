{

// Show the generated angles

TCanvas *c1 = new TCanvas();

c1->Divide(2,2);

c1->cd(1);
thph->Draw();

c1->cd(2);
thphc->Draw("");

c1->cd(3);
thpha->Draw("");

c1->cd(4);
thphaw->Draw("");


}
