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
static const int maxCorrelator = 12; // Somewhat abusing the setup as it is...
static const int maxHarmonic = 10; // Need to assess on case-by-case basis, but this gets you v2{8} and v3{6}
static const int maxPower = 9;
TComplex Qvector[maxHarmonic][maxPower]; // All needed Q-vector components
TComplex Q(int, int);
TComplex Recursion(int, int*);
TComplex Recursion(int, int*, int, int);


int main()
{

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();

  // --- histograms
  TH1D* hmult = new TH1D("hmult","charged multiplicity", 100, -0.5, 499.5);
  TH1D* hmult_selected = new TH1D("hmult_selected","charged multiplicity for selected particles", 100, -0.5, 499.5);

  TH1D* heta = new TH1D("heta", "eta distribution", 100, -3, 3); // this is correct
  TH1D* hphi = new TH1D("hphi", "phi distribution", 100, -pi, pi);//This should work but idk



  TH1D* heta_vec = new TH1D("heta_vec", "eta distribution", 100, -3, 3);
  TH1D* hphi_vec = new TH1D("hphi_vec", "phi distribution", 100, -pi, pi);


  TProfile* hmult_recursion[2][maxCorrelator];
  for ( int cs = 0; cs < 2; ++cs )
    {
      for(int c = 0; c < maxCorrelator; ++c )
        {
          hmult_recursion[cs][c] = new TProfile(Form("hmult_recursion_%d_%d",cs,c),"",700,-0.5,699.5,-1.1,1.1);
        }
    }


  for (int iEvent = 0; iEvent < 5; ++iEvent)
    {
      // --- for the generic formulas ---------
      for(int h=0;h<maxHarmonic;h++)
        {
          for(int w=0;w<maxPower;w++)
            {
              Qvector[h][w] = TComplex(0.,0.);
            } //  for(int p=0;p<maxPower;p++)
        } // for(int h=0;h<maxHarmonic;h++)
      // --------------------------------------

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

          if ( p.isFinal() && p.isCharged() ) hphi->Fill(p.phi());
          // fill hphi

          if(p.isFinal() && p.isCharged() && abs(p.eta()) < 2.5 && p.pT() > 0.5)
            {
              ++mult;
              //if ( p.isFinal() && p.isCharged() ) heta_vec->Fill(p.eta()); // this is redundant
              heta_vec->Fill(p.eta()); // this is correct
              //eta filler

              //if ( p.isFinal() && p.isCharged() ) hphi_vec->Fill(p.phi());
              hphi_vec->Fill(p.phi());
              // fill hphi

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
            } // loop over harmonics

        } // end loop over particles


      ////Defining usage and objects necessary for Recursion and
      ////higher cumulants after particle loop, but over
      ////event loop

      int harmonics_Two_Num[2] = {2,-2}; // 2, -2
      int harmonics_Two_Den[2] = {0,0}; // recursion gives right combinatorics
      TComplex twoRecursion = Recursion(2,harmonics_Two_Num)/Recursion(2,harmonics_Two_Den).Re();
      //double spwTwoRecursion = Recursion(2,harmonics_Two_Den).Re();
      double wTwoRecursion = 1.0;
      if(mult >= 2){
      hmult_recursion[0][0]->Fill(mult,twoRecursion.Re(),wTwoRecursion); // <<cos(h1*phi1+h2*phi2)>>
      hmult_recursion[1][0]->Fill(mult,twoRecursion.Im(),wTwoRecursion); // <<sin(h1*phi1+h2*phi2)>>
      }
      //  4-p correlations:
      //cout<<” => Calculating 4-p correlations (using recursion)...       \r”<<flush;
      int harmonics_Four_Num[4] = {2,2,-2,-2};
      int harmonics_Four_Den[4] = {0,0,0,0}; // recursion gives right combinatorics
      TComplex fourRecursion = Recursion(4,harmonics_Four_Num)/Recursion(4,harmonics_Four_Den).Re();
      //double spwFourRecursion = Recursion(4,harmonics_Four_Den).Re();
      double wFourRecursion = 1.0;
      if(mult >=4){
      hmult_recursion[0][2]->Fill(mult,fourRecursion.Re(),wFourRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
      hmult_recursion[1][2]->Fill(mult,fourRecursion.Im(),wFourRecursion); // <<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
      }
      //  6-p correlations:
      //cout<<” => Calculating 6-p correlations (using recursion)...       \r"<<flush;
      int harmonics_Six_Num[6] = {2,2,2,-2,-2,-2};
      int harmonics_Six_Den[6] = {0,0,0,0,0,0};
      TComplex sixRecursion = Recursion(6,harmonics_Six_Num)/Recursion(6,harmonics_Six_Den).Re();
      //double spwSixRecursion = Recursion(6,harmonics_Six_Den).Re();
      double wSixRecursion = 1.0;
      if(mult >= 6){
      hmult_recursion[0][4]->Fill(mult,sixRecursion.Re(),wSixRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
      hmult_recursion[1][4]->Fill(mult,sixRecursion.Im(),wSixRecursion); // <<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
      }
      //  8-p correlations:
      //cout<<” => Calculating 8-p correlations (using recursion)...       \r"<<flush;
      int harmonics_Eight_Num[8] = {2,2,2,2,-2,-2,-2,-2};
      int harmonics_Eight_Den[8] = {0,0,0,0,0,0,0,0};
      TComplex eightRecursion = Recursion(8,harmonics_Eight_Num)/Recursion(8,harmonics_Eight_Den).Re();
      //double spwEightRecursion = Recursion(8,harmonics_Eight_Den).Re();
      double wEightRecursion = 1.0;
      if(mult >= 8){
      hmult_recursion[0][6]->Fill(mult,eightRecursion.Re(),wEightRecursion);
      hmult_recursion[1][6]->Fill(mult,eightRecursion.Im(),wEightRecursion);
      }
      ////End of ripped codes from 'Boulder' codes



      // multiplicity distributions
      hmult->Fill( nCharged );
      hmult_selected->Fill( mult );
    } // end of loop over events



  pythia.stat(); // tell about some statistics for this run



  //Tfile for I/O stuff
  TFile* HistFile = new TFile("Output_Recursion.root","recreate");
  HistFile->cd();
  heta->Write();
  heta_vec->Write();
  hphi->Write();
  hphi_vec->Write();
  // --- write the recursion histograms
  for ( int cs = 0; cs < 2; ++cs )
    {
      for(int c = 0; c < maxCorrelator; ++c )
        {
          hmult_recursion[cs][c]->Write();
        }
    }
  HistFile->Close();

  // we'll make a file to output some histograms...

  return 0;

}
// --- from generic forumulas ----------------------------------------------------
TComplex Recursion(int n, int* harmonic)
{
  return Recursion(n,harmonic,1,0); // 1 and 0 are defaults from above
}

TComplex Recursion(int n, int* harmonic, int mult, int skip)
{
  // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
  // Kristjan Gulbrandsen (gulbrand@nbi.dk).

  int nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  int multp1 = mult+1;
  int nm2 = n-2;
  int counter1 = 0;
  int hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  int counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-double(mult)*c2;

}
TComplex Q(int n, int p)
{
  // --- for the generic formulas ---------

  // --------------------------------------
  // Using the fact that Q{-n,p} = Q{n,p}^*.
  if(n>=0){return Qvector[n][p];}
  return TComplex::Conjugate(Qvector[-n][p]);
} // TComplex Q(int n, int p)
// ------------------------
