//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file InitializeTracks.hh
//---------------------------------------------------------------------------//
#pragma once

#include "TrackInitializerStore.hh"
#include "base/KernelParamCalculator.cuda.hh"
#include "base/Span.hh"
#include "physics/base/Primary.hh"
#include "ParamPointers.hh"
#include "StatePointers.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// Predicate used to check whether the track at a given index in the track
// vector is alive
struct occupied
{
    size_type flag;

    occupied(size_type flag) : flag(flag){};

    CELER_FUNCTION bool operator()(const size_type x) { return x == flag; }
};

//---------------------------------------------------------------------------//
// Mark a track state as occupied, i.e., the particle in that slot is still
// alive, so a new track can't be initialized there.
CELER_CONSTEXPR_FUNCTION size_type occupied_flag()
{
    return numeric_limits<size_type>::max();
}

//---------------------------------------------------------------------------//
// Initialize the track states on device.
void initialize_tracks(StatePointers          states,
                       ParamPointers          params,
                       TrackInitializerStore& initializers);

//---------------------------------------------------------------------------//
// Find empty slots in the vector of track states and count the number of
// secondaries that survived cutoffs for each interaction.
void find_vacancies(StatePointers states, TrackInitializerStore& initializers);

//---------------------------------------------------------------------------//
// Create track initializers on device from primary particle.s
void process_primaries(span<const Primary>    primaries,
                       TrackInitializerStore& initializers);

//---------------------------------------------------------------------------//
// Create track initializers on device from secondary particles.
void process_secondaries(StatePointers          states,
                         ParamPointers          params,
                         TrackInitializerStore& initializers);

//---------------------------------------------------------------------------//
} // namespace celeritas
