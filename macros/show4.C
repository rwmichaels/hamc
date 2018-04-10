{

// Show the generated angles

TCanvas *c4 = new TCanvas();

c4->Divide(2,2);

c4->cd(1);
thph->Draw();

c4->cd(2);
thphc->Draw("");

c4->cd(3);
thpha->Draw("");

c4->cd(4);
thphaw->Draw("");


}
