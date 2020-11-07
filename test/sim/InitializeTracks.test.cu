//---------------------------------*-CUDA-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file InitializeTracks.test.cu
//---------------------------------------------------------------------------//
#include "sim/InitializeTracks.hh"
#include "InitializeTracks.test.hh"

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include "base/KernelParamCalculator.cuda.hh"

namespace celeritas_test
{
using namespace celeritas;

//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

__global__ void interact_kernel(StatePointers              states,
                                ParamPointers              params,
                                SecondaryAllocatorPointers secondaries)
{
    auto thread_id = celeritas::KernelParamCalculator::thread_id();
    if (thread_id.get() < states.size())
    {
        ParticleTrackView particle(params.particle, states.particle, thread_id);
        SimTrackView      sim(states.sim, thread_id);
        SecondaryAllocatorView allocate_secondaries(secondaries);

        // Allow the particle to interact and create secondaries
        Interactor interact(particle, allocate_secondaries);
        states.interactions[thread_id.get()] = interact();

        // Mark the track as dead if the particle was killed
        if (action_killed(states.interactions[thread_id.get()].action))
        {
            sim.alive() = false;
        }
    }
}

__global__ void
tracks_test_kernel(StatePointers states, ParamPointers params, double* output)
{
    auto thread_id = celeritas::KernelParamCalculator::thread_id();
    if (thread_id.get() < states.size())
    {
        ParticleTrackView particle(params.particle, states.particle, thread_id);
        output[thread_id.get()] = particle.energy().value();
    }
}

__global__ void
initializers_test_kernel(TrackInitializerPointers inits, double* output)
{
    auto thread_id = celeritas::KernelParamCalculator::thread_id().get();
    if (thread_id < inits.initializers.size())
    {
        TrackInitializer& track = inits.initializers[thread_id];
        output[thread_id]       = track.particle.energy.value();
    }
}

__global__ void
vacancies_test_kernel(TrackInitializerPointers inits, size_type* output)
{
    auto thread_id = celeritas::KernelParamCalculator::thread_id().get();
    if (thread_id < inits.vacancies.size())
    {
        output[thread_id] = inits.vacancies[thread_id];
    }
}

//---------------------------------------------------------------------------//
// TESTING INTERFACE
//---------------------------------------------------------------------------//

void interact(StatePointers              states,
              ParamPointers              params,
              SecondaryAllocatorPointers secondaries)
{
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(states.size());
    interact_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, secondaries);

    CELER_CUDA_CHECK_ERROR();
}

std::vector<double> tracks_test(StatePointers states, ParamPointers params)
{
    // Allocate memory for results
    thrust::device_vector<double> output(states.size());

    // Launch a kernel to check the properties of the initialized tracks
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(states.size());
    tracks_test_kernel<<<lparams.grid_size, lparams.block_size>>>(
        states, params, thrust::raw_pointer_cast(output.data()));

    CELER_CUDA_CHECK_ERROR();

    // Copy data back to host
    std::vector<double> host_output(states.size());
    thrust::copy(output.begin(), output.end(), host_output.begin());

    return host_output;
}

std::vector<double> initializers_test(TrackInitializerPointers inits)
{
    // Allocate memory for results
    thrust::device_vector<double> output(inits.initializers.size());

    // Launch a kernel to check the properties of the track initializers
    KernelParamCalculator calc_launch_params;
    auto lparams = calc_launch_params(inits.initializers.size());
    initializers_test_kernel<<<lparams.grid_size, lparams.block_size>>>(
        inits, thrust::raw_pointer_cast(output.data()));

    CELER_CUDA_CHECK_ERROR();

    // Copy data back to host
    std::vector<double> host_output(inits.initializers.size());
    thrust::copy(output.begin(), output.end(), host_output.begin());

    return host_output;
}

std::vector<size_type> vacancies_test(TrackInitializerPointers inits)
{
    // Allocate memory for results
    thrust::device_vector<size_type> output(inits.vacancies.size());

    // Launch a kernel to check the indices of the empty slots
    KernelParamCalculator calc_launch_params;
    auto                  lparams = calc_launch_params(inits.vacancies.size());
    vacancies_test_kernel<<<lparams.grid_size, lparams.block_size>>>(
        inits, thrust::raw_pointer_cast(output.data()));

    CELER_CUDA_CHECK_ERROR();

    // Copy data back to host
    std::vector<size_type> host_output(inits.vacancies.size());
    thrust::copy(output.begin(), output.end(), host_output.begin());

    return host_output;
}

//---------------------------------------------------------------------------//
} // namespace celeritas_test
