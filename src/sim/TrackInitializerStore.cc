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
    , parent_(secondaries.capacity())
    , vacancies_(particles.size())
    , secondary_counts_(particles.size())
    , track_count_(0)
{
    // Initialize vacancies to mark all track slots as initially empty
    std::vector<size_type> host_vacancies(vacancies_.size());
    std::iota(host_vacancies.begin(), host_vacancies.end(), 0);
    vacancies_.copy_to_device(make_span(host_vacancies));

    // Start with an empty vector of track initializers and parent thread IDs
    initializers_.resize(0);
    parent_.resize(0);
}

//---------------------------------------------------------------------------//
/*!
 * Get a view to the managed data.
 */
TrackInitializerPointers TrackInitializerStore::device_pointers()
{
    TrackInitializerPointers result;
    result.initializers     = initializers_.device_pointers();
    result.parent           = parent_.device_pointers();
    result.vacancies        = vacancies_.device_pointers();
    result.secondary_counts = secondary_counts_.device_pointers();
    result.track_count      = track_count_;

    ENSURE(result);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from primary particles.
 */
void TrackInitializerStore::create_from_primaries(span<const Primary> primaries)
{
    REQUIRE(primaries.size()
            <= initializers_.capacity() - initializers_.size());

    // Launch a kernel to create track initializers from primaries
    process_primaries(primaries, this->device_pointers());
    initializers_.resize(initializers_.size() + primaries.size());
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from secondary particles.
 *
 * Secondaries produced by each track are ordered arbitrarily in memory, and
 * the memory may be fragmented if not all secondaries survived cutoffs. For
 * example, after the interactions have been processed and cutoffs applied, the
 * track states and their secondaries might look like the following (where 'X'
 * indicates a track or secondary that did not survive):
 * \verbatim

   thread ID   | 0   1 2           3       4   5 6           7       8   9
   track ID    | 10  X 8           7       X   5 4           X       2   1
   secondaries | [X]   [X, 11, 12] [13, X] [X]   [14, X, 15] [X, 16] [X]

   \endverbatim
 *
 * Because the order in which threads receive a chunk of memory from the
 * secondary allocator is nondeterministic, the actual ordering of the
 * secondaries in memory is unpredictable; for instance:
 * \verbatim

  secondary storage | [X, 13, X, X, 11, 12, X, X, 16, 14, X, 15, X]

  \endverbatim
 *
 * When track initializers are created from secondaries, they are ordered by
 * thread ID to ensure reproducibility. If a track that produced secondaries
 * has died (e.g., thread ID 7 in the example above), one of its secondaries is
 * immediately used to fill that track slot:
 * \verbatim

   thread ID   | 0   1 2           3       4   5 6           7       8   9
   track ID    | 10  X 8           7       X   5 4           16      2   1
   secondaries | [X]   [X, 11, 12] [13, X] [X]   [14, X, 15] [X, X]  [X]

   \endverbatim
 *
 * This way, the geometry state is reused rather than initialized from the
 * position (which is expensive). This also prevents the geometry state from
 * being overwritten by another track's secondary, so if the track produced
 * multiple secondaries, the rest are still able to copy the parent's state.
 *
 * Track initializers are created from the remaining secondaries and are added
 * to the back of the vector. The thread ID of each secondary's parent is also
 * stored, so any new tracks initialized from secondaries produced in this
 * step can copy the geometry state from the parent. The indices of the empty
 * slots in the track vector are identified and stored as a sorted vector of
 * vacancies.
 * \verbatim

   track initializers | 11 12 13 14 15
   parent             | 2  2  3  6  6
   vacancies          | 1  4

   \endverbatim
 */
void TrackInitializerStore::create_from_secondaries(StatePointers states,
                                                    ParamPointers params)
{
    // Resize the vector of vacancies to be equal to the number of tracks
    vacancies_.resize(states.size());

    // Launch a kernel to identify which track slots are still alive and count
    // the number of surviving secondaries per track
    locate_alive(states, params, this->device_pointers());

    // Remove all elements in the vacancy vector that were flagged as active
    // tracks, leaving the (sorted) indices of the empty slots
    size_type num_vac = remove_if_alive(vacancies_.device_pointers());
    vacancies_.resize(num_vac);

    // Sum the total number secondaries produced in all interactions
    // TODO: if we don't have space for all the secondaries, we will need to
    // buffer the current track initializers to create room
    size_type num_sec = reduce_counts(secondary_counts_.device_pointers());
    REQUIRE(num_sec <= initializers_.capacity() - initializers_.size());

    // The exclusive prefix sum of the number of secondaries produced by each
    // track is used to get the start index in the vector of track initializers
    // for each thread. Starting at that index, each thread creates track
    // initializers from all surviving secondaries produced in its
    // interaction.
    exclusive_scan_counts(secondary_counts_.device_pointers());

    // Launch a kernel to create track initializers from secondaries
    process_secondaries(states, params, this->device_pointers());
    initializers_.resize(initializers_.size() + num_sec);
    parent_.resize(num_sec);
}

//---------------------------------------------------------------------------//
/*!
 * Initialize track states on device.
 *
 * Tracks created from secondaries produced in this step will have the geometry
 * state copied over from the parent instead of initialized from the position.
 * If there are more empty slots than new secondaries, they will be filled by
 * any track initializers remaining from previous steps using the position.
 */
void TrackInitializerStore::initialize_tracks(StatePointers states,
                                              ParamPointers params)
{
    // The number of new tracks to initialize is the smaller of the number of
    // empty slots in the track vector and the number of track initializers
    size_type num_tracks = std::min(vacancies_.size(), initializers_.size());
    vacancies_.resize(num_tracks);

    // Launch a kernel to initialize tracks on device
    init_tracks(states, params, this->device_pointers());
    initializers_.resize(initializers_.size() - num_tracks);

    // Update the total number of tracks initialized in the simulation
    track_count_ += num_tracks;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
