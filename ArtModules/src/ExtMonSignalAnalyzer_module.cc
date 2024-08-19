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

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"


namespace mu2e {

  //================================================================
  class ExtMonSignalAnalyzer : public art::SharedAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<double> nominalMomentum {
        Name("nominalMomentum"),
        Comment("Central momentum value to be used for some cuts"),
        4200.
      };

      fhicl::Atom<double> momentumTolerance {
        Name("momentumTolerance"),
        Comment("delta p value to be used for some cuts, along with nominalMomentum"),
        250.
      };

      fhicl::Atom<art::InputTag> inputHits {
        Name("inputHits"),
        Comment("The StepPointMCCollection to process")
      };

    };

    using Parameters = art::SharedAnalyzer::Table<Config>;
    explicit ExtMonSignalAnalyzer(const Parameters& conf, const art::ProcessingFrame&);

    void beginRun(const art::Run&, const art::ProcessingFrame&) override;
    void analyze(const art::Event& event, const art::ProcessingFrame&) override;

  private:
    art::InputTag inputHits_;
    const double proton_mass_;
    const double nominalMomentum_;
    const double momentumTolerance_;
    art::ServiceHandle<art::TFileService> tfs_;
    TH1D *vdnames_;
    TH1D *entrance_mom_;
    TH2D *entrance_mom_dir_;
    TH1D *coll1out_mom_;
    TH1D *exit_mom_;
    TH2D *exit_yvsxall_;
    TH2D *exit_yvsxzoom_;
  };

  //================================================================
  ExtMonSignalAnalyzer::ExtMonSignalAnalyzer(const Parameters& conf, const art::ProcessingFrame&)
    : SharedAnalyzer{conf}
    , inputHits_{conf().inputHits()}
    , proton_mass_(GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::proton).mass())
    , nominalMomentum_{conf().nominalMomentum()}
    , momentumTolerance_{conf().momentumTolerance()}
    , tfs_{art::ServiceHandle<art::TFileService>()}
    , vdnames_{
        (TH1::SetDefaultBufferSize(10000), // the comma operator
         tfs_->make<TH1D>("vdnames", "vdnames", 100, 0., -1.)
         )
      }
    , entrance_mom_{tfs_->make<TH1D>("entrance_mom", "entrance_mom", 200, 2500., 6000.)}
    , entrance_mom_dir_{tfs_->make<TH2D>("entrance_mom_dir", "entrance_mom_dir", 200, 0, -1, 200, 0, -1)}
    , coll1out_mom_{tfs_->make<TH1D>("coll1out_mom", "coll1out_mom", 200, 2500., 6000.)}
    , exit_mom_{tfs_->make<TH1D>("exit_mom", "exit_mom", 200, 2500., 6000.)}
    , exit_yvsxall_{tfs_->make<TH2D>("exit_yvsxall", "exit_yvsx all", 200, 0., -1., 200, 0., -1)}
    , exit_yvsxzoom_{nullptr}
  {
    serialize(art::SharedResource<art::TFileService>);
    vdnames_->SetOption("HIST,TEXT");
  }

  //================================================================
  //void ExtMonSignalAnalyzer::beginJob(const art::ProcessingFrame&) {
  void ExtMonSignalAnalyzer::beginRun(const art::Run&, const art::ProcessingFrame&) {
    if(!exit_yvsxzoom_) {
      GeomHandle<ExtMonFNALBuilding> emfb;
      const auto exitpoint = emfb->filter().exitInMu2e();
      const double detsize = 40; // mm
      exit_yvsxzoom_ = tfs_->make<TH2D>("exit_yvsxzoom",
                                        "exit_yvsx 4 cm x 4 cm",
                                        40, exitpoint.x()-detsize/2, exitpoint.x()+detsize/2,
                                        40, exitpoint.y()-detsize/2, exitpoint.y()+detsize/2);

      std::cout<<"ExtMonSignalAnalyzer: filter().exitInMu2e() = "<<exitpoint<<std::endl;
    }
  }

  //================================================================
  void ExtMonSignalAnalyzer::analyze(const art::Event& evt, const art::ProcessingFrame&) {
    auto const steps = evt.getValidHandle<StepPointMCCollection>(inputHits_);
    for(const auto& step: *steps) {
      const std::string nn{VirtualDetectorId::name(VirtualDetectorId::enum_type(step.volumeId()))};
      vdnames_->Fill(nn.c_str(), 1.);

      if(step.volumeId() == VirtualDetectorId::EMFC1Entrance) {
        const double pmag = step.momentum().mag();
        entrance_mom_->Fill(pmag);
        if(std::abs(pmag - nominalMomentum_)<momentumTolerance_) {
          const double dxdz = step.momentum().x() / step.momentum().z();
          const double dydz = step.momentum().y() / step.momentum().z();
          entrance_mom_dir_->Fill(dxdz, dydz);
        }
      }

      if(step.volumeId() == VirtualDetectorId::EMFC1Exit) {
        coll1out_mom_->Fill(step.momentum().mag());
      }

      if(step.volumeId() == VirtualDetectorId::EMFC2Exit) {
        exit_mom_->Fill(step.momentum().mag());
        exit_yvsxall_->Fill(step.position().x(), step.position().y());
        exit_yvsxzoom_->Fill(step.position().x(), step.position().y());
      }

    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonSignalAnalyzer)
