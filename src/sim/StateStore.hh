//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file StateStore.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/DeviceVector.hh"
#include "physics/base/ParticleStateStore.hh"
#include "geometry/GeoStateStore.hh"
#include "random/cuda/RngStateStore.hh"
#include "SimStateStore.hh"
#include "StatePointers.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Manage device data for tracks.
 */
class StateStore
{
  public:
    // Construct with the track state storage objects
    StateStore(ParticleStateStore& particle_states,
               GeoStateStore&      geo_states,
               SimStateStore&      sim_states,
               RngStateStore&      rng_states);

    // Get the total number of tracks
    size_type size() const { return particle_states_.size(); }

    // Get a view to the managed data
    StatePointers device_pointers();

  private:
    ParticleStateStore& particle_states_;
    GeoStateStore&            geo_states_;
    SimStateStore&            sim_states_;
    RngStateStore&            rng_states_;
    DeviceVector<Interaction> interactions_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
