// Andrei Gaponenko, 2024

#include <string>
#include <vector>
#include <cmath>
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
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/DataProducts/inc/CompressedPDGCode.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/compressPdgId.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

namespace mu2e {

  //================================================================
  class ExtMonSignalAnalyzer : public art::SharedAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> inputHits {
        Name("inputHits"),
        Comment("The StepPointMCCollection to process")
      };

    };

    using Parameters = art::SharedAnalyzer::Table<Config>;
    explicit ExtMonSignalAnalyzer(const Parameters& conf, const art::ProcessingFrame&);

    void analyze(const art::Event& event, const art::ProcessingFrame&) override;

  private:
    art::InputTag inputHits_;
    const double proton_mass_;
    art::ServiceHandle<art::TFileService> tfs_;
    TH1D *vdnames_;
  };

  //================================================================
  ExtMonSignalAnalyzer::ExtMonSignalAnalyzer(const Parameters& conf, const art::ProcessingFrame&)
    : SharedAnalyzer{conf}
    , inputHits_{conf().inputHits()}
    , proton_mass_(GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::proton).mass())
    , tfs_{art::ServiceHandle<art::TFileService>()}
    , vdnames_{
        (TH1::SetDefaultBufferSize(10000), // the comma operator
         tfs_->make<TH1D>("vdnames", "vdnames", 100, 0., -1.)
         )
      }
  {
    serialize(art::SharedResource<art::TFileService>);
  }

  //================================================================
  void ExtMonSignalAnalyzer::analyze(const art::Event& evt, const art::ProcessingFrame&) {
    auto const steps = evt.getValidHandle<StepPointMCCollection>(inputHits_);
    for(const auto& step: *steps) {
      const std::string nn{VirtualDetectorId::name(VirtualDetectorId::enum_type(step.volumeId()))};
      vdnames_->Fill(nn.c_str(), 1.);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonSignalAnalyzer)
