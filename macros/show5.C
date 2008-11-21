{

// Show the dispersion works if you set in prex.dat:
// iterate
// kick:track P0 0.003

TCanvas c1;

c1->Divide(2,1);

c1->cd(1);
xyfoc5->Draw();

c1->cd(2);
// Note, the "x" after xyfoc5 indicates the 2nd iteration histogram
// that automatically appears if you set "iterate" 
// in the control file prex.dat
xyfoc5x->Draw("");


}
