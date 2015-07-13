// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/tools/Filter.hh"

using namespace fastjet;

namespace Rivet {

    double D2(const Jet& jet) {

        const PseudoJet& pj = jet.pseudojet();

        // beta always = 1.0
        double ECF1 = 
            contrib::EnergyCorrelator(1, 1.0).result(pj);

        double ECF2 = 
            contrib::EnergyCorrelator(2, 1.0).result(pj);

        double ECF3 = 
            contrib::EnergyCorrelator(3, 1.0).result(pj);

        if (ECF1 == 0 || ECF2 == 0 || ECF3 == 0)
            return -999;

        // E correlation ratios
        double e2 = ECF2/(ECF1*ECF1);

        double e3 = ECF3/(ECF1*ECF1*ECF1);

        return e3/(e2*e2*e2);
    }




    class MC_BOOSTEDBOSON : public Analysis {
        public:

            MC_BOOSTEDBOSON()
                : Analysis("MC_BOOSTEDBOSON")
            {    }

            string toJSON(const Jet& jet, double weight) {
                const Jet& tjet = trimmer(jet);

                ostringstream os;
                os << "{" << endl <<
                    "\t\"weight\" : " << weight << "," << endl <<

                    "\t\"momentum\" : [" <<
                    jet.px() << ", " <<
                    jet.py() << ", " <<
                    jet.pz() << ", " <<
                    jet.energy() << "]," <<
                    endl <<

                    "\t\"D_2\" : " << D2(jet) << "," << endl <<

                    "\t\"trimmed momentum\" : [" <<
                    tjet.px() << ", " <<
                    tjet.py() << ", " <<
                    tjet.pz() << ", " <<
                    tjet.energy() << "]," <<
                    endl <<

                    "\t\"trimmed D_2\" : " << D2(tjet) << endl <<
                    "}";

                return os.str();
            }

            void init() {

                /*
                hPt = bookHisto1D("hPt", 50, 0, 1000*GeV,
                        "anti-$k_t$ $R=1.0$ jets",
                        "$p_T$ / GeV",
                        "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                hPtVsTrimmedPt = bookHisto2D("hPtVsTrimmedPt",
                        50, 0, 1000*GeV,
                        50, 0, 1000*GeV,
                        "anti-$k_t$ $R=1.0$ jets",
                        "untrimmed $p_T$ / GeV", "trimmed $p_T$ / GeV",
                        "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                hM = bookHisto1D("hM", 50, 0, 500*GeV,
                        "anti-$k_t$ $R=1.0$ jets",
                        "mass / GeV",
                        "$\\frac{d\\sigma}{d{\\rm mass}} / \\frac{\\rm pb}{\\rm GeV}$");

                hD2 = bookHisto1D("hD2", 40, 0, 4,
                        "anti-$k_t$ $R=1.0$ jets",
                        "$D_2$",
                        "$\\frac{d\\sigma}{dD_2} / {\\rm pb}$");
                */

                FastJets jetFS(
                        FinalState(Cuts::pT > 0*GeV && Cuts::abseta < 2.0),
                        FastJets::ANTIKT, 1.0);
                jetFS.useInvisibles(true);
                addProjection(jetFS, "jetFS");

                trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.05));

                return;
            }


            void analyze(const Event& event) {
                const double weight = event.weight();

                const Jets& jets =
                    applyProjection<FastJets>(event, "jetFS").jetsByPt(
                            Cuts::pT > 250*GeV && Cuts::abseta < 2.0);


                foreach(const Jet& jet, jets)
                    cout << toJSON(jet, weight) << endl;

                return;
            }


            void finalize() {
                double norm = crossSection()/sumOfWeights();
                /*
                scale(hPt, norm);
                scale(hM, norm);
                scale(hD2, norm);
                */

                return;
            }


        private:

            Filter trimmer;

            /*
            Histo1DPtr hD2;

            Histo1DPtr hM;
            Histo1DPtr hTrimmedM;

            Profile1DPtr hMeanMVsTrimmedM;
            Profile1DPtr hMeanTrimmedMVsM;
            Profile1DPtr hMeanMVsPt;
            Profile1DPtr hMeanTrimmedMVsPt;
            Profile1DPtr hMeanTrimmedMVsTrimmedPt;

            Histo2DPtr hMVsTrimmedM;
            Histo2DPtr hTrimmedMVsM;
            Histo2DPtr hMVsTrimmedPt;
            Histo2DPtr hMVsPt;
            Histo2DPtr hTrimmedMVsPt;
            Histo2DPtr hTrimmedMVsTrimmedPt;

            Histo1DPtr hPt;
            Histo1DPtr hTrimmedPt;
            Histo1DPtr hMeanPtVsTrimmedPt;
            Histo1DPtr hMeanTrimmedPtPt;
            Histo1DPtr hPtVsTrimmedPt;
            */
    };

    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDBOSON);
}
