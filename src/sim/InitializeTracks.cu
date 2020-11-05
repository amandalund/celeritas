//---------------------------------*-CUDA-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file InitializeTracks.cu
//---------------------------------------------------------------------------//
#include "InitializeTracks.hh"

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/remove.h>
#include <thrust/scan.h>
#include <vector>
#include "base/Atomics.hh"
#include "base/DeviceVector.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
/*!
 * Initialize the track states on device. The track initializers are created
 * from either primary particles or secondaries. The new tracks are inserted
 * into empty slots (vacancies) in the track vector.
 */
__global__ void
initialize_tracks_kernel(const StatePointers            states,
                         const ParamPointers            params,
                         const TrackInitializerPointers initializers)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < initializers.vacancies.size())
    {
        // Get the track initializer, starting from the back of the vector
        const TrackInitializer& init
            = initializers.tracks[initializers.tracks.size() - thread_id - 1];

        // Index of the empty slot to create the new track in
        size_type empty_id = initializers.vacancies[thread_id];

        // Initialize particle physics data
        ParticleTrackView particle(
            params.particle, states.particle, ThreadId(empty_id));
        particle = init.particle;

        // Initialize geometry state
        // GeoTrackView geo(params.geo, states.geo, ThreadId(empty_id));
        // geo = init.geo;

        // Initialize simulation state
        SimTrackView sim(states.sim, ThreadId(empty_id));
        sim = init.sim;
    }
}

//---------------------------------------------------------------------------//
/*!
 * Find empty slots in the track vector and count the number of secondaries
 * that survived cutoffs for each interaction.
 */
__global__ void
find_vacancies_kernel(const StatePointers            states,
                      const TrackInitializerPointers initializers)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < states.size())
    {
        // Determine which indices in the track states are available for
        // initializing new particles
        SimTrackView sim(states.sim, ThreadId(thread_id));
        if (sim.alive())
        {
            initializers.vacancies[thread_id] = occupied_flag();
        }
        else
        {
            initializers.vacancies[thread_id] = thread_id;
        }

        // Count how many secondaries survived cutoffs for each track
        initializers.secondary_counts[thread_id] = 0;
        const Interaction& result = states.interactions[thread_id];
        for (const auto& secondary : result.secondaries)
        {
            if (secondary)
            {
                ++initializers.secondary_counts[thread_id];
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from primary particles.
 */
__global__ void
process_primaries_kernel(const span<const Primary>      primaries,
                         const TrackInitializerPointers initializers)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < primaries.size())
    {
        size_type         offset_id = initializers.tracks.size();
        TrackInitializer& init = initializers.tracks[offset_id + thread_id];

        // Construct a track initializer from a primary particle
        init = primaries[thread_id];
    }
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from secondary particles.
 */
__global__ void
process_secondaries_kernel(const StatePointers            states,
                           const ParamPointers            params,
                           const TrackInitializerPointers initializers)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < states.size())
    {
        // Construct the state accessors
        // GeoTrackView geo(params.geo, states.geo, ThreadId(thread_id));
        SimTrackView sim(states.sim, ThreadId(thread_id));

        // Starting index in the vector of track initializers
        size_type index = initializers.secondary_counts[thread_id];

        const Interaction& result = states.interactions[thread_id];
        for (const auto& secondary : result.secondaries)
        {
            // If the secondary survived cutoffs
            if (secondary)
            {
                // Construct a track initializer from a secondary
                TrackInitializer& init = initializers.tracks[index];
                init.particle.def_id   = secondary.def_id;
                init.particle.energy   = secondary.energy;
                init.geo.dir           = secondary.direction;
                // init.geo.pos           = geo.pos();
                unsigned int track_id = initializers.track_count + 1 + index++;
                init.sim.track_id     = TrackId{track_id};
                init.sim.parent_id    = sim.track_id();
                init.sim.event_id     = sim.event_id();
                init.sim.alive        = true;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// KERNEL INTERFACE
//---------------------------------------------------------------------------//
/*!
 * Initialize the track states on device.
 */
void initialize_tracks(StatePointers          states,
                       ParamPointers          params,
                       TrackInitializerStore& initializers)
{
    // The number of new tracks to initialize is the smaller of the number of
    // empty slots in the track vector and the number of track initializers
    size_type num_new_tracks
        = std::min(initializers.num_vacancies(), initializers.size());
    initializers.num_vacancies() = num_new_tracks;

    // Initialize tracks on device
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(num_new_tracks);
    initialize_tracks_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, initializers.device_pointers());

    CELER_CUDA_CALL(cudaDeviceSynchronize());

    // Resize the vector of track initializers
    initializers.resize(initializers.size() - num_new_tracks);
    initializers.track_count() += num_new_tracks;
}

//---------------------------------------------------------------------------//
/*!
 * Find empty slots in the vector of tracks and count the number of secondaries
 * that survived cutoffs for each interaction.
 */
void find_vacancies(StatePointers states, TrackInitializerStore& initializers)
{
    // Resize the vector of vacancies to be equal to the number of tracks
    initializers.num_vacancies() = states.size();

    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(states.size());
    find_vacancies_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, initializers.device_pointers());

    CELER_CUDA_CALL(cudaDeviceSynchronize());

    // Remove all the elements in the vacancy vector that were flagged as
    // active tracks, so we are left with a vector containing the (sorted)
    // indices of the empty slots
    span<size_type> vacancies = initializers.device_pointers().vacancies;
    thrust::device_ptr<size_type> end = thrust::remove_if(
        thrust::device_pointer_cast(vacancies.data()),
        thrust::device_pointer_cast(vacancies.data() + vacancies.size()),
        occupied(occupied_flag()));

    // Resize the vector of vacancies to be equal to the number of empty slots
    initializers.num_vacancies() = thrust::raw_pointer_cast(end)
                                   - vacancies.data();
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers from primary particles.
 */
void process_primaries(span<const Primary>    primaries,
                       TrackInitializerStore& initializers)
{
    REQUIRE(primaries.size() <= initializers.capacity() - initializers.size());

    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(primaries.size());
    process_primaries_kernel<<<lparams.grid_size, lparams.block_size>>>(
        primaries, initializers.device_pointers());

    CELER_CUDA_CALL(cudaDeviceSynchronize());

    // Resize the vector of track initializers
    initializers.resize(initializers.size() + primaries.size());
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers from secondary particles.
 */
void process_secondaries(StatePointers          states,
                         ParamPointers          params,
                         TrackInitializerStore& initializers)
{
    // Sum the total number secondaries produced in all interactions
    span<size_type> counts = initializers.device_pointers().secondary_counts;
    size_type       num_secondaries = thrust::reduce(
        thrust::device_pointer_cast(counts.data()),
        thrust::device_pointer_cast(counts.data()) + counts.size(),
        0,
        thrust::plus<size_type>());

    // TODO: if we don't have space for all the secondaries, we will need to
    // buffer the current track initializers
    REQUIRE(num_secondaries <= initializers.capacity() - initializers.size());

    // The exclusive prefix sum of the number of secondaries produced in each
    // interaction is used to get the starting index in the vector of track
    // initializers for creating initializers from secondaries from an
    // interaction
    thrust::exclusive_scan(
        thrust::device_pointer_cast(counts.data()),
        thrust::device_pointer_cast(counts.data()) + counts.size(),
        counts.data(),
        0);

    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(states.size());
    process_secondaries_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, initializers.device_pointers());

    CELER_CUDA_CALL(cudaDeviceSynchronize());

    // Resize the vector of track initializers
    initializers.resize(initializers.size() + num_secondaries);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
