//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file EventReader.test.cc
//---------------------------------------------------------------------------//
#include "io/EventReader.hh"

#include "gtest/Main.hh"
#include "gtest/Test.hh"
#include "base/Span.hh"
#include "physics/base/ParticleParams.hh"

using celeritas::EventReader;
using celeritas::ParticleParams;
using celeritas::Primary;
namespace pdg = celeritas::pdg;

//---------------------------------------------------------------------------//
// TEST HARNESS
//---------------------------------------------------------------------------//

class EventReaderTest : public celeritas::Test,
                        public testing::WithParamInterface<const char*>
{
  protected:
    void SetUp() override
    {
        using celeritas::ParticleDef;
        using celeritas::PDGNumber;

        // Create shared standard model particle data
        ParticleParams::VecAnnotatedDefs defs
            = {{{"proton", pdg::proton()},
                {938.27208816, 1, ParticleDef::stable_decay_constant()}},
               {{"d_quark", PDGNumber(1)},
                {4.7, -1 / 3, ParticleDef::stable_decay_constant()}},
               {{"anti_u_quark", PDGNumber(-2)},
                {2.2, -2 / 3, ParticleDef::stable_decay_constant()}},
               {{"w_minus", PDGNumber(-24)}, {8.0379e4, 0, 3.168e24}},
               {{"gamma", pdg::gamma()},
                {0, 0, ParticleDef::stable_decay_constant()}}};
        particle_params_ = std::make_shared<ParticleParams>(std::move(defs));
    }

    std::string                     filename_;
    std::shared_ptr<ParticleParams> particle_params_;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_P(EventReaderTest, all)
{
    // Determine the event record format and open the file
    filename_ = this->test_data_path("io", GetParam());

    // Read event and primary particle information from event record
    EventReader read_event(filename_.c_str(), particle_params_);
    auto        event_record = read_event();

    EXPECT_EQ(8, event_record.primaries.size());

    // Expected PDG: 2212, 1, 2212, -2, 22, -24, 1, -2
    const int expected_def_id[] = {0, 1, 0, 2, 4, 3, 1, 2};

    const double expected_energy[] = {
        7.e6, 3.2238e4, 7.e6, 5.7920e4, 4.233e3, 8.5925e4, 2.9552e4, 5.6373e4};

    const double expected_direction[][3] = {
        {0, 0, 1},
        {2.326451417389850e-2, -4.866936365179566e-2, 9.985439676959555e-1},
        {0, 0, -1},
        {-5.260794237813896e-2, -3.280442747570201e-1, -9.431963518790131e-1},
        {-9.009470900796461e-1, 2.669997932835038e-2, -4.331067443262500e-1},
        {5.189457940206315e-2, -7.074356638330033e-1, -7.048700122475354e-1},
        {-8.273504806466310e-2, 9.750892208717103e-1, 2.058055469649411e-1},
        {7.028153760960004e-2, -8.780402697122620e-1, -4.733981307893478e-1}};

    for (std::size_t i = 0; i < event_record.primaries.size(); ++i)
    {
        const auto& primary = event_record.primaries[i];

        // Check that the particle types were read correctly
        EXPECT_EQ(expected_def_id[i], primary.def_id.get());

        // Check that the event IDs match
        EXPECT_EQ(0, primary.event_id.get());

        // Check that the position, direction, and energy  was read correctly
        const double expected_position[] = {0, 0, 0};
        EXPECT_VEC_SOFT_EQ(expected_position, primary.position);
        EXPECT_VEC_SOFT_EQ(expected_direction[i], primary.direction);
        EXPECT_DOUBLE_EQ(expected_energy[i], primary.energy);
    }
}

INSTANTIATE_TEST_SUITE_P(EventReaderTests,
                         EventReaderTest,
                         testing::Values("event-record.hepmc3",
                                         "event-record.hepmc2",
                                         "event-record.hepevt"));
