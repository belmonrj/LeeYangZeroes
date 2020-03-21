// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

// this is a copy with very minor modification by Ron Belmont
// Also edited to by Aidan Lytle to produce an eta distribution

#include "TROOT.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TComplex.h"

#include "Pythia8/Pythia.h"

using namespace Pythia8;

const double pi = 3.1415926;

// Stuff for cumulants greater than 2
static const int maxHarmonic = 10; // Need to assess on case-by-case basis, but this gets you v2{8} and v3{6}
static const int maxPower = 9;
TComplex Qvector[maxHarmonic][maxPower]; // All needed Q-vector components

static const int rbins = 100;
static const int thetabins = 10;
TComplex Qtheta[maxHarmonic][thetabins];
TComplex genfunS[rbins][thetabins]; // sum
TComplex genfunP[rbins][thetabins]; // product

int main()
{

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.init();

  // --- histograms
  TH1D* hmult = new TH1D("hmult","charged multiplicity", 100, -0.5, 499.5);
  TH1D* hmult_selected = new TH1D("hmult_selected","charged multiplicity for selected particles", 100, -0.5, 499.5);

  TH1D* heta = new TH1D("heta", "eta distribution", 100, -3, 3); // this is correct
  TH1D* hPhi = new TH1D("hPhi", "phi distribution", 100, -pi, pi);//This should work but idk

  TH1D* heta_vec = new TH1D("heta_vec", "eta distribution", 100, -3, 3);
  TH1D* hPhi_vec = new TH1D("hPhi_vec", "phi distribution", 100, -pi, pi);


  TProfile* tp1f_c22mult = new TProfile("tp1f_c22mult","c22 vs mult",100,-0.5,499.5,-1e10,1e10);

  //LeeYang Histograms
  TProfile* LeeYang = new TProfile("LeeYang","",rbins,0,1,-1e10,1e10);
  const int multbins = 10; // we'll do just a few for now, may do more later
  TProfile* LeeYangHistosS[multbins][thetabins];
  TProfile* LeeYangHistosP[multbins][thetabins];



  for (Int_t i = 0; i < multbins; i++){
    for (Int_t j = 0; j < thetabins; j++){
      LeeYangHistosS[i][j] = new TProfile(Form("LeeYangS_mult%d_theta%d",i,j),"",rbins,0,1,-1e10,1e10);
      LeeYangHistosP[i][j] = new TProfile(Form("LeeYangP_mult%d_theta%d",i,j),"",rbins,0,1,-1e10,1e10);
    }
  }




  for (int iEvent = 0; iEvent < 5000; ++iEvent)
    {
      // --- for the generic formulas ---------
      for(int h=0;h<maxHarmonic;h++)
        {
          for(int w=0;w<maxPower;w++)
            {
              Qvector[h][w] = TComplex(0.,0.);
            } //  for(int p=0;p<maxPower;p++)
          for ( int i = 0; i < thetabins; ++i )
            {
              Qtheta[h][i] = TComplex(0.0,0.0);
            }
        } // for(int h=0;h<maxHarmonic;h++)
      // --------------------------------------
      for ( int i = 0; i < rbins; ++i )
        {
          for ( int j = 0; j < thetabins; ++j )
            {
              genfunS[i][j] = TComplex(0.0,0.0); // initialize to 0, calculate directly
              genfunP[i][j] = TComplex(1.0,1.0); // initialize to 1, calcualte via product
            }
        }

      if (!pythia.next()) continue;
      // Find number of all final charged particles and fill histogram.

      Event& event = pythia.event;
      // First we find the particles we need for the analysis,
      // as well as event multiplicity.
      vector<Particle*> parts;
      int mult = 0;
      int nCharged = 0;
      double Q2x = 0;
      double Q2y = 0;
      //int heta = 0; // this is extremely bad

      for (int i = 0; i < event.size(); ++i)
        {
          // Particle short notation
          Particle& p = event[i];

          // count all charged particles
          if (p.isFinal() && p.isCharged()) ++nCharged;

          // Apply simple, particle level, cuts.

          if ( p.isFinal() && p.isCharged() ) heta->Fill(p.eta()); // this is correct
          //eta filler

          if ( p.isFinal() && p.isCharged() ) hPhi->Fill(p.phi());
          // fill hPhi

          if(p.isFinal() && p.isCharged() && abs(p.eta()) < 2.5 && p.pT() > 0.5)
            {
              ++mult;
              //if ( p.isFinal() && p.isCharged() ) heta_vec->Fill(p.eta()); // this is redundant
              heta_vec->Fill(p.eta()); // this is correct
              //eta filler

              //if ( p.isFinal() && p.isCharged() ) hPhi_vec->Fill(p.phi());
              hPhi_vec->Fill(p.phi());
              // fill hPhi

              if(p.pT() > 1.0 && p.pT() < 3.0) parts.push_back(&p);
            }

          double phi = p.phi();
          Q2x += cos(2*phi);
          // Qx
          Q2y += sin(2*phi);
          // Qy

          // ---------------------------
          // --- calculate the Q-vectors
          // ---
          for(int h=0;h<maxHarmonic;h++)
            {
              for(int w=0;w<maxPower;w++)
                {
                  //if(bUseWeights){wPhiToPowerP = pow(wPhi,p);} // no weights for us...
                  Qvector[h][w] += TComplex(cos(h*phi),sin(h*phi));
                } //  for(int w=0;w<maxPower;w++)
              for ( int thetabin = 0; thetabin < thetabins; ++thetabin )
                {
                  float theta_width = (2*pi)/thetabins;
                  float theta = thetabin*theta_width;
                  float angle = phi-theta;
                  Qtheta[h][thetabin] += TComplex(cos(h*angle),sin(h*angle));
                }
            } // loop over harmonics

          // ---------------------------------------------
          // --- calculate the product generating function
          // ---
          for ( int rbin = 0; rbin < rbins; ++rbin )
            {
              double r = (double(rbin)/rbins);
              for ( int thetabin = 0; thetabin < thetabins; ++thetabin )
                {
                  float theta_width = (2*pi)/thetabins;
                  float theta = thetabin*theta_width;
                  float angle = phi-theta;
                  genfunP[rbin][thetabin] *= TComplex(1.0,r*cos(2*(angle)));
                }
            }

        } // end loop over particles


      // -------------------------------------------
      // --- fill the generating function histograms
      // ---
      int multbin = mult/10;
      if ( multbin >= multbins )
        {
          cout << "Skipping event with multiplicity " << mult << ", event index " << iEvent << endl;
          continue;
        }
      for ( int rbin = 0; rbin < rbins; ++rbin )
        {
          double r = (double(rbin)/rbins);
          for ( int thetabin = 0; thetabin < thetabins; ++thetabin )
            {
              genfunS[rbin][thetabin] = TComplex::Exp(r*Qtheta[2][thetabin]); // generating function from Q-vectors
              LeeYangHistosS[multbin][thetabin]->Fill(r,genfunS[rbin][thetabin].Rho());
              LeeYangHistosP[multbin][thetabin]->Fill(r,genfunP[rbin][thetabin].Rho());
            }
        }



      // calculate the 2 cumulant
      double c22 = (((Q2x*Q2x) + (Q2y*Q2y)) - mult)/(mult*(mult-1));
      // fill TProfile to compute average cumulant vs multiplicity
      tp1f_c22mult->Fill(mult,c22);

      // multiplicity distributions
      hmult->Fill( nCharged );
      hmult_selected->Fill( mult );
    } // end of loop over events



  pythia.stat(); // tell about some statistics for this run



  //Tfile for I/O stuff
  TFile* HistFile = new TFile("Output_LYZ.root","recreate");
  HistFile->cd();
  heta->Write();
  heta_vec->Write();
  hPhi->Write();
  hPhi_vec->Write();
  tp1f_c22mult->Write();
  LeeYang->Write();
  for ( int i = 0; i < multbins; ++i )
    {
      for ( int j = 0; j < thetabins; ++j )
        {
          LeeYangHistosS[i][j]->Write();
          LeeYangHistosP[i][j]->Write();
        }
    }
  HistFile->Close();

  // we'll make a file to output some histograms...

  return 0;

}
