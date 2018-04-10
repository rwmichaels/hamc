{

// Show the dispersion works if you set in prex.dat:
// iterate
// kick:track P0 0.003

TCanvas *c5 = new TCanvas;

c5->Divide(2,1);

c5->cd(1);
xyfoc5->Draw();

c5->cd(2);
// Note, the "x" after xyfoc5 indicates the 2nd iteration histogram
// that automatically appears if you set "iterate" 
// in the control file prex.dat
xyfoc5x->Draw("");


}
