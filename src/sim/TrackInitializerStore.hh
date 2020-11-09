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

    // Create track initializers on device from primary particles
    void create_from_primaries(span<const Primary> primaries);

    // Create track initializers on device from secondary particles.
    void create_from_secondaries(StatePointers states, ParamPointers params);

    // Initialize track states on device.
    void initialize_tracks(StatePointers states, ParamPointers params);

  private:
    // Track initializers created from primaries or secondaries
    DeviceVector<TrackInitializer> initializers_;

    // Thread ID of the secondary's parent
    DeviceVector<size_type> parent_;

    // Index of empty slots in track vector
    DeviceVector<size_type> vacancies_;

    // Numbe of surviving secondaries produced in each interaction
    DeviceVector<size_type> secondary_counts_;

    // Total number of tracks initialized in the simulation
    size_type track_count_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
