//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file TrackInitializerStore.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/DeviceVector.hh"
#include "TrackInitializerPointers.hh"
#include "physics/base/ParticleStateStore.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Manage device data for track initializers.
 */
class TrackInitializerStore
{
  public:
    // Construct with the maximum number of track initializers to store on
    // device
    explicit TrackInitializerStore(ParticleStateStore& states,
                                   size_type           capacity);

    // Get the number of elements
    size_type capacity() const { return initializers_.size(); }

    // Get the number of elements
    size_type size() const { return size_; }

    // Change the size (without changing capacity)
    void resize(size_type size);

    // Modifier for the number of empty track slots
    size_type& num_vacancies() { return num_vacancies_; }

    // Modifier for the total number of initialized tracks
    size_type& track_count() { return track_count_; }

    // Get a view to the managed data
    TrackInitializerPointers device_pointers();

  private:
    // TODO: shrinkable device vector
    // Track initializers
    DeviceVector<TrackInitializer> initializers_;
    size_type                      size_;

    // Indices of empty slots in the track vector
    DeviceVector<size_type> vacancies_;
    size_type               num_vacancies_;

    // Number of secondaries produced in each interaction
    DeviceVector<size_type> secondary_counts_;

    // Total number of tracks that have been initialized in the simulation.
    // This is used to set the track ID of secondaries
    size_type track_count_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
