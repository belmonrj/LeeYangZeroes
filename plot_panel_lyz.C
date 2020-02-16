const int multbins = 10;
const int thetabins = 10;


void plot_panel_lyz()
{

  TCanvas* c1 = new TCanvas("c1","",1000,500);
  c1->Divide(5,2);

  TCanvas* c2 = new TCanvas("c2","",1000,500);
  c2->Divide(5,2);

  TFile* file = TFile::Open("Output_LYZ.root","read");

  TProfile* lyz_summ[multbins][thetabins];
  TProfile* lyz_prod[multbins][thetabins];

  for ( int imult = 0; imult < multbins; ++imult )
    {
      for ( int itheta = 0; itheta < thetabins; ++itheta )
        {
          lyz_summ[imult][itheta] = (TProfile*)file->Get(Form("LeeYangS_mult%d_theta%d",imult,itheta));
          lyz_prod[imult][itheta] = (TProfile*)file->Get(Form("LeeYangP_mult%d_theta%d",imult,itheta));
          c1->cd(itheta+1);
          lyz_summ[imult][itheta]->Draw();
          c2->cd(itheta+1);
          lyz_prod[imult][itheta]->Draw();
        }
      c1->cd(0);
      c1->Print(Form("Figures/test_panels_summ_mult%d.png",imult));
      c2->cd(0);
      c2->Print(Form("Figures/test_panels_prod_mult%d.png",imult));
    }


}
