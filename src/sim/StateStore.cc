//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file StateStore.cc
//---------------------------------------------------------------------------//
#include "StateStore.hh"

#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with the track state storage objects.
 */
StateStore::StateStore(ParticleStateStore& particle_states,
                       GeoStateStore&      geo_states,
                       SimStateStore&      sim_states,
                       RngStateStore&      rng_states)
    : particle_states_(particle_states)
    , geo_states_(geo_states)
    , sim_states_(sim_states)
    , rng_states_(rng_states)
    , interactions_(this->size())
{
    REQUIRE(geo_states_.size() == this->size());
    REQUIRE(sim_states_.size() == this->size());
    REQUIRE(rng_states_.size() == this->size());
}

//---------------------------------------------------------------------------//
/*!
 * Get a view to the managed data.
 */
StatePointers StateStore::device_pointers()
{
    StatePointers result;
    result.particle = particle_states_.device_pointers();
    result.geo          = geo_states_.device_pointers();
    result.sim          = sim_states_.device_pointers();
    result.rng          = rng_states_.device_pointers();
    result.interactions = interactions_.device_pointers();
    ENSURE(result);
    return result;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
