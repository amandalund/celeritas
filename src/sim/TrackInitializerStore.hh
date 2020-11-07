//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file TrackInitializerStore.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/DeviceVector.hh"
#include "physics/base/ParticleStateStore.hh"
#include "physics/base/SecondaryAllocatorStore.hh"
#include "ParamPointers.hh"
#include "StatePointers.hh"
#include "TrackInitializerPointers.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Manage device data for track initializers.
 */
class TrackInitializerStore
{
  public:
    // Construct with the maximum number of track initializers to store on
    // device
    explicit TrackInitializerStore(ParticleStateStore&      particles,
                                   SecondaryAllocatorStore& secondaries);

    // Get a view to the managed data
    TrackInitializerPointers device_pointers();

    // Initialize track states on device.
    void initialize_tracks(StatePointers states, ParamPointers params);

    // Find vacant slots and count surviving secondaries per track
    void find_vacancies(StatePointers states);

    // Create track initializers on device from primary particles
    void create_from_primaries(span<const Primary> primaries);

    // Create track initializers on device from secondary particles.
    void create_from_secondaries(StatePointers states, ParamPointers params);

  private:
    DeviceVector<TrackInitializer> initializers_;
    DeviceVector<size_type>        vacancies_;
    DeviceVector<size_type>        secondary_counts_;
    size_type                      track_count_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
