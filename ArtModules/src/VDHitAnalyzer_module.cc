// Andrei Gaponenko, 2024

#include <string>
#include <vector>
#include <cmath>
#include <optional>
#include <limits>

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/TupleAs.h"

#include "art/Framework/Core/SharedAnalyzer.h"
#include "art/Framework/Principal/Event.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TH1.h"
#include "TH2.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/CompressedPDGCode.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/compressPdgId.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

namespace mu2e {

  // Main beam window drawing pos = (-3118.0, 122.8, -9398.0)
  double constexpr dumpwin_x0 = +3118.0; // Matt Slabaugh's X sign is flipped
  double constexpr dumpwin_y0 =   122.8;

  double constexpr hwin_x1 = 2500.;
  double constexpr hwin_x2 = 4500.;
  double constexpr hwin_y1 = -750;
  double constexpr hwin_y2 = +750;

  //================================================================
  namespace {
    class WinCut {
      double x0_;
      double y0_;
      double rCut_;
    public:

      WinCut(double x, double y, double r) : x0_{x}, y0_{y}, rCut_{r} {}

      bool pass(const StepPointMC& hit) const {
        return sqrt(std::pow(hit.position().x() - x0_, 2)
                    +std::pow(hit.position().y() - y0_, 2))
          < rCut_;
      }
    };
  }

  //================================================================
  class VDHitAnalyzer : public art::SharedAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> input {
        Name("input"),
        Comment("The StepPointMCCollection to process")
      };

      fhicl::Atom<std::string> vdName {
        Name("vdName"),
        Comment("The name of the VD we want to look at")
      };

      fhicl::OptionalAtom<int> selectCharge {
        Name("selectCharge"),
        Comment("Histogram positive, negative, or neutral particles per the sign of this integer.  No cut if not set.")
      };
    };

    using Parameters = art::SharedAnalyzer::Table<Config>;
    explicit VDHitAnalyzer(const Parameters& conf, const art::ProcessingFrame&);

    void analyze(const art::Event& event, const art::ProcessingFrame&) override;

  private:
    art::InputTag input_;
    VirtualDetectorId vdId_;

    std::optional<int> selectCharge_;

    WinCut cutdump_;
    WinCut cutq_;

    art::ServiceHandle<art::TFileService> tfs_;
    TH1D *pdgcounts_;

    TH1D *momentum1_;
    TH1D *ek1_;

    TH2D *hitxyCount_;
    TH2D *hitxyESum_;


    TH1D *cutdump_pdgcounts_;
    TH1D *cutdump_momentum_;

    TH1D *cutq_pdgcounts_;
    TH1D *cutq_momentum_;

    GlobalConstantsHandle<ParticleDataList> pdgList_;
  };

  //================================================================
  VDHitAnalyzer::VDHitAnalyzer(const Parameters& conf, const art::ProcessingFrame&)
    : SharedAnalyzer{conf}
    , input_{conf().input()}
    , vdId_{conf().vdName()}
    , selectCharge_{conf().selectCharge()}
    , cutdump_{dumpwin_x0, dumpwin_y0, 150. /* try to get the hot spot by eye, smaller than the window*/}
    , cutq_{3340., 85., 60. /* all by eye*/}

    , tfs_{art::ServiceHandle<art::TFileService>()}

    , pdgcounts_{compressPDGCodeHisto(tfs_)}

    , momentum1_{tfs_->make<TH1D>("momentum1", "Particle momentum",200,0.,10000.)}
    , ek1_{tfs_->make<TH1D>("ek1", "Particle kinetic energy",170,0.,8500.)}

    //, hitxyCount_{tfs_->make<TH2D>("hitxyCount", "Y vs X of hits, count",200,0.,0., 200, 0.,0.)}
    //, hitxyESum_{tfs_->make<TH2D>("hitxyESum", "Y vs X of hits, kinetic energy sum",200,0.,0., 200, 0.,0.)}

      // Main beam window drawing pos = (-3118.0, 122.8, -9398.0)
    , hitxyCount_{tfs_->make<TH2D>("hitxyCount", "Y vs X of hits, count",200,hwin_x1, hwin_x2, 150, hwin_y1, hwin_y2)}
    , hitxyESum_{tfs_->make<TH2D>("hitxyESum", "Y vs X of hits, kinetic energy sum",200,hwin_x1, hwin_x2, 150, hwin_y1, hwin_y2)}

    , cutdump_pdgcounts_{compressPDGCodeHisto(tfs_, "cutdump_pdgcounts")}
    , cutdump_momentum_{tfs_->make<TH1D>("cutdump_momentum", "Momentum for cutdump",200,0.,10000.)}

    , cutq_pdgcounts_{compressPDGCodeHisto(tfs_, "cutq_pdgcounts")}
    , cutq_momentum_{tfs_->make<TH1D>("cutq_momentum", "Momentum for cutq",200,0.,10000.)}

  {
    serialize(art::SharedResource<art::TFileService>);

    // normalize the charge value
    if(selectCharge_) {
      if(selectCharge_.value() < 0) selectCharge_ = -1;
      if(selectCharge_.value() + 0) selectCharge_ = +1;
    }
  }

  //================================================================
  void VDHitAnalyzer::analyze(const art::Event& evt, const art::ProcessingFrame&) {
    auto const steps = evt.getValidHandle<StepPointMCCollection>(input_);
    for(const auto& step: *steps) {
      pdgcounts_->Fill(compressPDGCode(step.simParticle()->pdgId()));

      if(!selectCharge_ || (selectCharge_ == pdgList_->particle(step.simParticle()->pdgId()).charge() )) {

        const double mom = step.momentum().mag();
        momentum1_->Fill(mom);

        const double mass = pdgList_->particle(step.simParticle()->pdgId()).mass();
        const double ek = sqrt(std::pow(mass,2)+std::pow(mom,2)) - mass;
        ek1_->Fill(ek);


        hitxyCount_->Fill(step.position().x(), step.position().y());
        hitxyESum_->Fill(step.position().x(), step.position().y(), ek);

        if(cutdump_.pass(step)) {
          cutdump_pdgcounts_->Fill(compressPDGCode(step.simParticle()->pdgId()));
          cutdump_momentum_->Fill(step.momentum().mag());
        }

        if(cutq_.pass(step)) {
          cutq_pdgcounts_->Fill(compressPDGCode(step.simParticle()->pdgId()));
          cutq_momentum_->Fill(step.momentum().mag());
        }

      }
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VDHitAnalyzer)
