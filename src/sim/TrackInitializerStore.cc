//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file TrackInitializerStore.cc
//---------------------------------------------------------------------------//
#include "TrackInitializerStore.hh"

#include <numeric>
#include "InitializeTracks.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with the maximum number of elements to allocate on device.
 */
TrackInitializerStore::TrackInitializerStore(
    ParticleStateStore& particles, SecondaryAllocatorStore& secondaries)
    : initializers_(secondaries.capacity())
    , vacancies_(particles.size())
    , secondary_counts_(particles.size())
    , track_count_(0)
{
    // Initialize vacancies to mark all track slots as initially empty
    std::vector<size_type> host_vacancies(vacancies_.size());
    std::iota(host_vacancies.begin(), host_vacancies.end(), 0);
    vacancies_.copy_to_device(make_span(host_vacancies));

    // Start with an empty vector of track initializers
    initializers_.resize(0);
}

//---------------------------------------------------------------------------//
/*!
 * Get a view to the managed data.
 */
TrackInitializerPointers TrackInitializerStore::device_pointers()
{
    TrackInitializerPointers result;
    result.initializers     = initializers_.device_pointers();
    result.vacancies        = vacancies_.device_pointers();
    result.secondary_counts = secondary_counts_.device_pointers();
    result.track_count      = track_count_;

    ENSURE(result);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Initialize track states on device.
 */
void TrackInitializerStore::initialize_tracks(StatePointers states,
                                              ParamPointers params)
{
    // The number of new tracks to initialize is the smaller of the number of
    // empty slots in the track vector and the number of track initializers
    size_type num_tracks = std::min(vacancies_.size(), initializers_.size());
    vacancies_.resize(num_tracks);

    // Launch a kernel to initialize tracks on device
    process_tracks(states, params, this->device_pointers());

    // Resize the vector of track initializers
    initializers_.resize(initializers_.size() - num_tracks);

    // Update the total number of tracks initialized in the simulation
    track_count_ += num_tracks;
}

//---------------------------------------------------------------------------//
/*!
 * Find empty slots in the vector of track states and count the number of
 * secondaries that survived cutoffs for each interaction.
 */
void TrackInitializerStore::find_vacancies(StatePointers states)
{
    // Resize the vector of vacancies to be equal to the number of tracks
    vacancies_.resize(states.size());

    // Launch a kernel to find the indices of the empty track slots and the
    // number of surviving secondaries per track
    process_post_interaction(states, this->device_pointers());

    // Remove all elements in the vacancy vector that were flagged as active
    // tracks, leaving the (sorted) indices of the empty slots
    size_type num_vac = remove_occupied(vacancies_.device_pointers());

    // Resize the vector of vacancies to be equal to the number of empty slots
    vacancies_.resize(num_vac);
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from primary particles
 */
void TrackInitializerStore::create_from_primaries(span<const Primary> primaries)
{
    REQUIRE(primaries.size()
            <= initializers_.capacity() - initializers_.size());

    // Launch a kernel to create track initializers from primaries
    process_primaries(primaries, this->device_pointers());

    // Resize the vector of track initializers
    initializers_.resize(initializers_.size() + primaries.size());
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from secondary particles
 */
void TrackInitializerStore::create_from_secondaries(StatePointers states,
                                                    ParamPointers params)
{
    // Sum the total number secondaries produced in all interactions
    size_type num_sec = reduce_counts(secondary_counts_.device_pointers());

    // TODO: if we don't have space for all the secondaries, we will need to
    // buffer the current track initializers to create room
    REQUIRE(num_sec <= initializers_.capacity() - initializers_.size());

    // The exclusive prefix sum of the number of secondaries produced in each
    // interaction is used to get the starting index for each track in the
    // vector of track initializers. The tracks will create track initializers
    // for all of the secondaries produced in the interaction.
    exclusive_scan_counts(secondary_counts_.device_pointers());

    // Launch a kernel to create track initializers from secondaries
    process_secondaries(states, params, this->device_pointers());

    // Resize the vector of track initializers
    initializers_.resize(initializers_.size() + num_sec);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
