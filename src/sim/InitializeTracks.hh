//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file InitializeTracks.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/KernelParamCalculator.cuda.hh"
#include "base/Span.hh"
#include "physics/base/Primary.hh"
#include "ParamPointers.hh"
#include "StatePointers.hh"
#include "TrackInitializerPointers.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// Predicate used to check whether the track at a given index in the track
// vector is alive
struct alive
{
    size_type flag;

    alive(size_type flag) : flag(flag){};

    CELER_FUNCTION bool operator()(const size_type x) { return x == flag; }
};

//---------------------------------------------------------------------------//
// Mark a track state as alive to indicate a new track can't be initialized
// there.
CELER_CONSTEXPR_FUNCTION size_type flag_alive()
{
    return numeric_limits<size_type>::max();
}

//---------------------------------------------------------------------------//
// Initialize the track states on device.
void init_tracks(StatePointers            states,
                 ParamPointers            params,
                 TrackInitializerPointers inits);

//---------------------------------------------------------------------------//
// Identify which tracks are still alive and count the number of secondaries
// that survived cutoffs for each interaction.
void locate_alive(StatePointers            states,
                  ParamPointers            params,
                  TrackInitializerPointers inits);

//---------------------------------------------------------------------------//
// Create track initializers on device from primary particles
void process_primaries(span<const Primary>      primaries,
                       TrackInitializerPointers inits);

//---------------------------------------------------------------------------//
// Create track initializers on device from secondary particles.
void process_secondaries(StatePointers            states,
                         ParamPointers            params,
                         TrackInitializerPointers inits);

//---------------------------------------------------------------------------//
// Remove all elements in the vacancy vector that were flagged as alive
size_type remove_if_alive(span<size_type> vacancies);

//---------------------------------------------------------------------------//
// Sum the total number of surviving secondaries.
size_type reduce_counts(span<size_type> counts);

//---------------------------------------------------------------------------//
// Calculate the exclusive prefix sum of the number of surviving secondaries
void exclusive_scan_counts(span<size_type> counts);

//---------------------------------------------------------------------------//
} // namespace celeritas
