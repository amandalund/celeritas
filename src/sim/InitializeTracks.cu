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
__global__ void process_tracks_kernel(const StatePointers            states,
                                      const ParamPointers            params,
                                      const TrackInitializerPointers inits)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < inits.vacancies.size())
    {
        // Get the track initializer, starting from the back of the vector
        const TrackInitializer& track
            = inits.initializers[inits.initializers.size() - thread_id - 1];

        // Index of the empty slot to create the new track in
        size_type empty_id = inits.vacancies[thread_id];

        // Initialize particle physics data
        ParticleTrackView particle(
            params.particle, states.particle, ThreadId(empty_id));
        particle = track.particle;

        // Copy the geometry state from the parent if possible and initialize
        // it otherwise
        GeoTrackView geo(params.geo, states.geo, ThreadId(empty_id));
        if (thread_id < inits.parent.size())
        {
            unsigned int parent_id
                = inits.parent[inits.parent.size() - thread_id - 1];
            GeoTrackView parent(params.geo, states.geo, ThreadId(parent_id));
            geo.copy_state(parent, track.geo.dir);
        }
        else
        {
            geo = track.geo;
        }

        // Initialize simulation state
        SimTrackView sim(states.sim, ThreadId(empty_id));
        sim = track.sim;
    }
}

//---------------------------------------------------------------------------//
/*!
 * Find empty slots in the track vector and count the number of secondaries
 * that survived cutoffs for each interaction. If the track is dead and
 * produced secondaries, fill the empty track slot with one of the secondaries.
 */
__global__ void
process_interaction_change_kernel(const StatePointers            states,
                                  const ParamPointers            params,
                                  const TrackInitializerPointers inits)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < states.size())
    {
        const Interaction& result = states.interactions[thread_id];

        // Index of the secondary to copy to the parent's track slot if the
        // parent has died
        long int secondary_id = -1;

        // Count how many secondaries survived cutoffs for each track
        inits.secondary_counts[thread_id] = 0;
        for (auto i = 0; i < result.secondaries.size(); ++i)
        {
            const auto& secondary = result.secondaries[i];
            if (secondary)
            {
                if (secondary_id < 0)
                {
                    secondary_id = i;
                }
                ++inits.secondary_counts[thread_id];
            }
        }

        // Determine which indices in the track states are available for
        // initializing new particles
        SimTrackView sim(states.sim, ThreadId(thread_id));
        if (sim.alive())
        {
            // The track is alive
            inits.vacancies[thread_id] = flag_alive();
        }
        else if (secondary_id >= 0)
        {
            // The track is dead and produced secondaries -- fill the empty
            // track slot with the first secondary
            Secondary& secondary = result.secondaries[secondary_id];

            // Keep the parent's geometry state
            GeoTrackView geo(params.geo, states.geo, ThreadId(thread_id));
            geo.copy_state(geo, secondary.direction);

            // Initialize the particle state from the secondary
            ParticleTrackView particle(
                params.particle, states.particle, ThreadId(thread_id));
            particle = {secondary.def_id, secondary.energy};

            // TODO: correct track ID
            // Initialize the simulation state
            sim = {TrackId{}, sim.track_id(), sim.event_id(), true};

            // Mark this track as alive and the secondary as processed
            inits.vacancies[thread_id] = flag_alive();
            secondary                  = Secondary::from_cutoff();
            --inits.secondary_counts[thread_id];
        }
        else
        {
            // The track is dead and did not produce secondaries
            inits.vacancies[thread_id] = thread_id;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from primary particles.
 */
__global__ void process_primaries_kernel(const span<const Primary> primaries,
                                         const TrackInitializerPointers inits)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < primaries.size())
    {
        size_type         offset_id = inits.initializers.size();
        TrackInitializer& track = inits.initializers[offset_id + thread_id];

        // Construct a track initializer from a primary particle
        track = primaries[thread_id];
    }
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers on device from secondary particles.
 */
__global__ void process_secondaries_kernel(const StatePointers states,
                                           const ParamPointers params,
                                           const TrackInitializerPointers inits)
{
    auto thread_id = KernelParamCalculator::thread_id().get();
    if (thread_id < states.size())
    {
        // Construct the state accessors
        GeoTrackView geo(params.geo, states.geo, ThreadId(thread_id));
        SimTrackView sim(states.sim, ThreadId(thread_id));

        // Offset in the vector of track initializers
        size_type index = inits.secondary_counts[thread_id];

        const Interaction& result = states.interactions[thread_id];
        for (const auto& secondary : result.secondaries)
        {
            // If the secondary survived cutoffs
            if (secondary)
            {
                // Thread ID of the secondary's parent
                inits.parent[index] = thread_id;

                // Construct a track initializer from a secondary
                TrackInitializer& track
                    = inits.initializers[index + inits.initializers.size()];
                track.particle.def_id = secondary.def_id;
                track.particle.energy = secondary.energy;
                track.geo.dir         = secondary.direction;
                track.geo.pos         = geo.pos();
                unsigned int track_id = inits.track_count + 1 + index++;
                track.sim.track_id    = TrackId{track_id};
                track.sim.parent_id   = sim.track_id();
                track.sim.event_id    = sim.event_id();
                track.sim.alive       = true;
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
void process_tracks(StatePointers            states,
                    ParamPointers            params,
                    TrackInitializerPointers inits)
{
    // Initialize tracks on device
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(inits.vacancies.size());
    process_tracks_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, inits);

    CELER_CUDA_CALL(cudaDeviceSynchronize());
}

//---------------------------------------------------------------------------//
/*!
 * Find empty slots in the vector of tracks and count the number of secondaries
 * that survived cutoffs for each interaction.
 */
void process_interaction_change(StatePointers            states,
                                ParamPointers            params,
                                TrackInitializerPointers inits)
{
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(states.size());
    process_interaction_change_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, inits);

    CELER_CUDA_CALL(cudaDeviceSynchronize());
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers from primary particles.
 */
void process_primaries(span<const Primary>      primaries,
                       TrackInitializerPointers inits)
{
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(primaries.size());
    process_primaries_kernel<<<lparams.grid_size, lparams.block_size>>>(
        primaries, inits);

    CELER_CUDA_CALL(cudaDeviceSynchronize());
}

//---------------------------------------------------------------------------//
/*!
 * Create track initializers from secondary particles.
 */
void process_secondaries(StatePointers            states,
                         ParamPointers            params,
                         TrackInitializerPointers inits)
{
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(states.size());
    process_secondaries_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, inits);

    CELER_CUDA_CALL(cudaDeviceSynchronize());
}

//---------------------------------------------------------------------------//
/*!
 * Remove all elements in the vacancy vector that were flagged as active
 * tracks.
 */
size_type remove_if_occupied(span<size_type> vacancies)
{
    thrust::device_ptr<size_type> end = thrust::remove_if(
        thrust::device_pointer_cast(vacancies.data()),
        thrust::device_pointer_cast(vacancies.data() + vacancies.size()),
        occupied(flag_alive()));

    CELER_CUDA_CALL(cudaDeviceSynchronize());

    // New size of the vacancy vector
    size_type result = thrust::raw_pointer_cast(end) - vacancies.data();
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Sum the total number of surviving secondaries.
 */
size_type reduce_counts(span<size_type> counts)
{
    size_type result = thrust::reduce(
        thrust::device_pointer_cast(counts.data()),
        thrust::device_pointer_cast(counts.data()) + counts.size(),
        size_type(0),
        thrust::plus<size_type>());

    CELER_CUDA_CALL(cudaDeviceSynchronize());
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the exclusive prefix sum of the number of surviving secondaries
 * from each interaction.
 */
void exclusive_scan_counts(span<size_type> counts)
{
    thrust::exclusive_scan(
        thrust::device_pointer_cast(counts.data()),
        thrust::device_pointer_cast(counts.data()) + counts.size(),
        counts.data(),
        size_type(0));

    CELER_CUDA_CALL(cudaDeviceSynchronize());
}

//---------------------------------------------------------------------------//
} // namespace celeritas
