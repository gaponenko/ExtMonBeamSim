// Andrei Gaponenko, 2024

#include <string>
#include <vector>
#include <cmath>

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

  //================================================================
  class PSEHitAnalyzer : public art::SharedAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> input {
        Name("input"),
        Comment("The StepPointMCCollection to process")
      };

    };

    using Parameters = art::SharedAnalyzer::Table<Config>;
    explicit PSEHitAnalyzer(const Parameters& conf, const art::ProcessingFrame&);

    void analyze(const art::Event& event, const art::ProcessingFrame&) override;

  private:
    art::InputTag input_;
    const double proton_mass_;
    art::ServiceHandle<art::TFileService> tfs_;
    TH1D *pdgcounts_;
    TH1D *ek1_;
    TH1D *momentum1_;
    TH2D *hitxy_;
    TH2D *hitzy_;
    TH2D *hitzx_;
  };

  //================================================================
  PSEHitAnalyzer::PSEHitAnalyzer(const Parameters& conf, const art::ProcessingFrame&)
    : SharedAnalyzer{conf}
    , input_{conf().input()}
    , proton_mass_(GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::proton).mass())
    , tfs_{art::ServiceHandle<art::TFileService>()}
    , pdgcounts_{compressPDGCodeHisto(tfs_)}
    , ek1_{tfs_->make<TH1D>("ek1", "Proton kinetic energy",170,0.,8500.)}
    , momentum1_{tfs_->make<TH1D>("momentum1", "Particle momentum",200,0.,10000.)}
    , hitxy_{tfs_->make<TH2D>("hitxy", "Y vs X of StepPointMC",200,0.,0., 200, 0.,0.)}
    , hitzy_{tfs_->make<TH2D>("hitzy", "Y vs Z of StepPointMC",200,0.,0., 200, 0.,0.)}
    , hitzx_{tfs_->make<TH2D>("hitzx", "X vs Z of StepPointMC",200,0.,0., 200, 0.,0.)}
  {
    serialize(art::SharedResource<art::TFileService>);
  }

  //================================================================
  void PSEHitAnalyzer::analyze(const art::Event& evt, const art::ProcessingFrame&) {
    auto const steps = evt.getValidHandle<StepPointMCCollection>(input_);
    for(const auto& step: *steps) {
      pdgcounts_->Fill(compressPDGCode(step.simParticle()->pdgId()));
      hitxy_->Fill(step.position().x(), step.position().y());
      hitzy_->Fill(step.position().z(), step.position().y());
      hitzx_->Fill(step.position().z(), step.position().x());

      const double p = step.momentum().mag();
      momentum1_->Fill(p);
      if(PDGCode::proton == step.simParticle()->pdgId()) {
        ek1_->Fill(sqrt(std::pow(proton_mass_,2)+std::pow(p,2)) - proton_mass_);
      }
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PSEHitAnalyzer)
