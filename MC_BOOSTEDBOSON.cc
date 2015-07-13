// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/tools/Filter.hh"

namespace Rivet {

    double D2(const Jet& jet) {

        const fastjet::PseudoJet& pj = jet.pseudojet();

        // beta always = 1.0
        double ECF1 = 
            fastjet::contrib::EnergyCorrelator(1, 1.0).result(pj);

        double ECF2 = 
            fastjet::contrib::EnergyCorrelator(2, 1.0).result(pj);

        double ECF3 = 
            fastjet::contrib::EnergyCorrelator(3, 1.0).result(pj);

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


            void init() {

                hPt = bookHisto1D("hPt", 50, 0, 1000*GeV,
                            "anti-$k_t$ $R=1.0$ jets",
                            "$p_T$ / GeV",
                            "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                hM = bookHisto1D("hM", 50, 0, 500*GeV,
                            "anti-$k_t$ $R=1.0$ jets",
                            "mass / GeV",
                            "$\\frac{d\\sigma}{d{\\rm mass}} / \\frac{\\rm pb}{\\rm GeV}$");

                hD2 = bookHisto1D("hD2", 40, 0, 4,
                        "anti-$k_t$ $R=1.0$ jets",
                        "$D_2$",
                        "$\\frac{d\\sigma}{dD_2} / {\\rm pb}$");

                FastJets jetFS(
                        FinalState(Cuts::pT > 0*GeV && Cuts::abseta < 2.0),
                        FastJets::ANTIKT, 1.0);
                jetFS.useInvisibles(true);
                addProjection(jetFS, "jetFS");

                return;
            }


            void analyze(const Event& event) {
                const double weight = event.weight();

                const Jets& jets =
                    applyProjection<FastJets>(event, "jetFS").jetsByPt(
                            Cuts::pT > 250*GeV && Cuts::abseta < 2.0);


                foreach(const Jet& jet, jets) {
                    hPt->fill(jet.pT(), weight);
                    hM->fill(jet.mass(), weight);
                    hD2->fill(D2(jet), weight);
                }

                return;
            }


            void finalize() {
                double norm = crossSection()/sumOfWeights();
                scale(hPt, norm);
                scale(hM, norm);
                scale(hD2, norm);

                return;
            }


        private:

            Histo1DPtr hD2;
            Histo1DPtr hM;
            Histo1DPtr hPt;
    };

    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDBOSON);
}
