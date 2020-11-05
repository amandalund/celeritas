//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file TrackInitializer.hh
//---------------------------------------------------------------------------//
#pragma once

#include "Types.hh"
#include "physics/base/Interaction.hh"
#include "physics/base/Primary.hh"
#include "physics/base/Secondary.hh"
#include "physics/base/ParticleTrackView.hh"
#include "geometry/GeoTrackView.hh"
#include "SimTrackView.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Lightweight version of a track used to initialize new tracks from primaries
 * or secondaries
 */
struct TrackInitializer
{
    // Simulation state
    SimTrackView::Initializer_t      sim;
    GeoTrackView::Initializer_t      geo;
    ParticleTrackView::Initializer_t particle;

    // Initialize from a primary particle
    CELER_FUNCTION TrackInitializer& operator=(const Primary& primary)
    {
        particle.def_id = primary.def_id;
        particle.energy = primary.energy;
        geo.dir         = primary.direction;
        geo.pos         = primary.position;
        sim.event_id    = primary.event_id;
        sim.track_id    = primary.track_id;
        sim.parent_id   = TrackId{0};
        sim.alive       = true;
        return *this;
    }
};

//---------------------------------------------------------------------------//
} // namespace celeritas
