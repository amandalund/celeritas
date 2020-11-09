//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file InitializeTracks.test.hh
//---------------------------------------------------------------------------//
#include "physics/base/Interaction.hh"
#include "physics/base/SecondaryAllocatorPointers.hh"
#include "physics/base/SecondaryAllocatorView.hh"
#include "physics/base/ParticleStateStore.hh"
#include "physics/base/ParticleTrackView.hh"
#include "sim/TrackInitializerStore.hh"
#include <vector>

namespace celeritas_test
{
using namespace celeritas;

//---------------------------------------------------------------------------//
// TESTING INTERFACE
//---------------------------------------------------------------------------//
//! Interactor
struct Interactor
{
    CELER_FUNCTION Interactor(const ParticleTrackView& particle,
                              SecondaryAllocatorView&  allocate_secondaries)
        : inc_energy(particle.energy())
        , allocate_secondaries(allocate_secondaries)
    {
    }

    CELER_FUNCTION Interaction operator()()
    {
        // Create secondary particles for some of the tracks
        unsigned long int alloc_size
            = static_cast<int>(inc_energy.value() / 100 - 1) % 4;
        Secondary* allocated = this->allocate_secondaries(alloc_size);
        if (!allocated)
        {
            return Interaction::from_failure();
        }
        Interaction result;

        // Kill some tracks to create vacancies in the track vector
        if (static_cast<int>(inc_energy.value()) % 3)
        {
            result.action = Action::scattered;
        }
        else
        {
            result.action = Action::absorbed;
        }

        // Initialize secondaries, including one with zero energy to account
        // for secondaries being pruned when applying cutoffs
        result.secondaries = {allocated, alloc_size};
        for (unsigned long int i = 0; i < alloc_size; ++i)
        {
            units::MevEnergy energy{i * inc_energy.value() / 100};
            if (energy.value() > 0)
            {
                result.secondaries[i].def_id = ParticleDefId(0);
                result.secondaries[i].energy = energy;
            }
            else
            {
                result.secondaries[i] = Secondary::from_cutoff();
            }
        }
        return result;
    }

    units::MevEnergy        inc_energy;
    SecondaryAllocatorView& allocate_secondaries;
};

//---------------------------------------------------------------------------//
//! Launch a kernel to process interactions
void interact(StatePointers              states,
              ParamPointers              params,
              SecondaryAllocatorPointers secondaries);

//---------------------------------------------------------------------------//
//! Launch a kernel to get the energies of the initialized tracks
std::vector<double> tracks_test(StatePointers states, ParamPointers params);

//---------------------------------------------------------------------------//
//! Launch a kernel to get the energies of the track initializers created from
//! primaries or secondaries
std::vector<double> initializers_test(TrackInitializerPointers initializers);

//---------------------------------------------------------------------------//
//! Launch a kernel to get the indices of the vacant slots in the track vector
std::vector<size_type> vacancies_test(TrackInitializerPointers initializers);

//---------------------------------------------------------------------------//
} // namespace celeritas_test
