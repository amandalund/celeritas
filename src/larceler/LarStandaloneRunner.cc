//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file larceler/LarStandaloneRunner.cc
//---------------------------------------------------------------------------//
#include "LarStandaloneRunner.hh"

#include <memory>
#include <utility>
#include <lardataobj/Simulation/OpDetBacktrackerRecord.h>
#include <lardataobj/Simulation/SimEnergyDeposit.h>

#include "corecel/Assert.hh"
#include "corecel/Macros.hh"
#include "corecel/io/Logger.hh"
#include "geocel/DetectorParams.hh"
#include "geocel/GeantGeoParams.hh"
#include "geocel/Types.hh"
#include "geocel/VolumeParams.hh"
#include "geocel/detail/LengthUnits.hh"
#include "celeritas/Quantities.hh"
#include "celeritas/Types.hh"
#include "celeritas/inp/StandaloneInput.hh"
#include "celeritas/optical/CoreParams.hh"
#include "celeritas/optical/Runner.hh"

#include "Convert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with problem setup and detector ID coordinates.
 *
 * The detector "channels" (coordinates) should be input as a vector.
 */
LarStandaloneRunner::LarStandaloneRunner(Input&& i, VecReal3 const& det_coords)
{
    CELER_EXPECT(!det_coords.empty());

    i.problem.detectors.callback
        = [this](SpanCelerHits h) { return this->hit(h); };
    runner_ = std::make_shared<optical::Runner>(std::move(i));

    // Map detector coordinates
    // TODO: probably better to use native geometry than G4
    // auto geo = runner_->params()->geometry();
    auto geo = global_geant_geo().lock();
    CELER_ASSERT(geo);
    // TODO: volumes should be in the params
    auto vols = global_volumes().lock();
    CELER_ASSERT(vols);
    auto dets = runner_->params()->detectors();
    CELER_ASSERT(dets);

    detids_.resize(det_coords.size());
    for (auto i : range(det_coords.size()))
    {
        auto inst_id = geo->find_volume_instance_at(det_coords[i]);
        CELER_VALIDATE(inst_id,
                       << "could not find a volume at " << det_coords[i]
                       << " [" << lengthunits::native_label << "]");
        VolumeId vol_id = vols->volume(inst_id);
        CELER_ASSERT(vol_id);
        auto det_id = dets->detector_id(vol_id);
        CELER_VALIDATE(det_id,
                       << "detector " << i << " point " << det_coords[i]
                       << " [" << lengthunits::native_label
                       << "] is not associated with a detector");
        // Save lar-to-celer mapping
        detids_[i] = det_id;

        // Save celer-to-lar helpers: this point shouldn't be associated with
        // any existing detector channel
        if (det_id < btr_helpers_.size())
        {
            btr_helpers_.resize(det_id.get() + 1);
        }
        CELER_VALIDATE(!btr_helpers_[det_id.get()],
                       << "multiple points map to detector " << i << " = '"
                       << dets->detector_labels().at(det_id) << "'");
        btr_helpers_[det_id.get()] = std::make_unique<sim::OBTRHelper>(i);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Run scintillation optical photons from a single set of energy steps.
 *
 * \todo With Cherenkov enabled we would need to determine the incident
 * particle's charge and the pre- and post-step speed.
 */
auto LarStandaloneRunner::operator()(VecSED const& sed) -> VecBTR
{
    CELER_EXPECT(!sed.empty());

    std::vector<celeritas::optical::GeneratorDistributionData> gdd;
    gdd.reserve(sed.size());

    for (auto const& edep : sed)
    {
        // Convert LArSoft sim edeps to Celeritas generator distribution
        // data
        celeritas::optical::GeneratorDistributionData data;
        data.type = GeneratorType::scintillation;
        data.num_photons = edep.NumPhotons();
        data.primary = id_cast<PrimaryId>(edep.TrackID());
        data.step_length = convert_from_larsoft<LarsoftLen>(edep.StepLength());
        //! XXX Given post-step point find optical material
        data.material = OptMatId{0};
        // Assume continuous energy loss along the step
        //! \todo For neutral particles, set this to 0 (LED at post-step
        //! point)
        data.continuous_edep_fraction = 1;
        data.points[StepPoint::pre].time
            = convert_from_larsoft<LarsoftTime>(edep.StartT());
        data.points[StepPoint::pre].pos
            = convert_from_larsoft<LarsoftLen>(edep.Start());
        data.points[StepPoint::post].time
            = convert_from_larsoft<LarsoftTime>(edep.EndT());
        data.points[StepPoint::post].pos
            = convert_from_larsoft<LarsoftLen>(edep.End());
        CELER_ASSERT(data);
        gdd.push_back(data);
    }

    // Execute
    auto result = (*runner_)(make_span(std::as_const(gdd)));

    CELER_ASSERT(result.counters.generators.size() == 1);
    auto const& gen = result.counters.generators.front();
    CELER_LOG(debug) << "Transported " << gen.num_generated
                     << " optical photons from " << gen.buffer_size
                     << " sim energy deposits a total of "
                     << result.counters.steps << " steps over "
                     << result.counters.step_iters << " step iterations";

    // Convert BTR helpers to BTRs in the LarSoft order
    VecBTR btrs;
    btrs.reserve(detids_.size());
    for (auto i : range(detids_.size()))
    {
        // "move" with lvalue reference and then recreate
        btrs.emplace_back(*btr_helpers_[i]);
        btr_helpers_[i] = std::make_unique<sim::OBTRHelper>(i);
    }

    return btrs;
}

//---------------------------------------------------------------------------//
/*!
 * Convert Celeritas hits to optical backtracker records.
 */
void LarStandaloneRunner::hit(Span<optical::DetectorHit const> hits)
{
    CELER_LOG(info) << "Processing " << hits.size() << "hits";
    for (auto& h : hits)
    {
        CELER_ASSERT(h.detector < btr_helpers_.size());
        CELER_ASSERT(btr_helpers_[h.detector.get()]);

        Real3 larpos{convert_to_larsoft<LarsoftLen>(h.position[0]),
                     convert_to_larsoft<LarsoftLen>(h.position[1]),
                     convert_to_larsoft<LarsoftLen>(h.position[2])};
        btr_helpers_[h.detector.get()]->AddScintillationPhotonsToMap(
            h.primary.get(),
            convert_to_larsoft<LarsoftTime>(h.time),
            /* num photons = */ 1,
            larpos.data(),
            value_as<units::MevEnergy>(h.energy));
    }
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
