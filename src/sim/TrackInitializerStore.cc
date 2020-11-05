//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file TrackInitializerStore.cc
//---------------------------------------------------------------------------//
#include "TrackInitializerStore.hh"

#include "base/Assert.hh"
#include <numeric>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with the maximum number of elements to allocate on device.
 */
TrackInitializerStore::TrackInitializerStore(ParticleStateStore& states,
                                             size_type           capacity)
    : initializers_(capacity)
    , size_(0)
    , vacancies_(states.size())
    , num_vacancies_(0)
    , secondary_counts_(states.size())
    , track_count_(0)
{
    REQUIRE(capacity > 0);

    // Initialize vacancies to mark all track slots as empty
    std::vector<size_type> host_vacancies(vacancies_.size());
    std::iota(host_vacancies.begin(), host_vacancies.end(), 0);
    vacancies_.copy_to_device(make_span(host_vacancies));
    num_vacancies_ = vacancies_.size();
}

//---------------------------------------------------------------------------//
/*!
 * Change the size (without changing capacity).
 */
void TrackInitializerStore::resize(size_type size)
{
    REQUIRE(size <= this->capacity());
    size_ = size;
}

//---------------------------------------------------------------------------//
/*!
 * Get a view to the managed data.
 */
TrackInitializerPointers TrackInitializerStore::device_pointers()
{
    TrackInitializerPointers result;
    result.tracks = {initializers_.device_pointers().data(), this->size()};
    result.vacancies
        = {vacancies_.device_pointers().data(), this->num_vacancies()};
    result.secondary_counts = secondary_counts_.device_pointers();
    result.track_count      = this->track_count();

    ENSURE(result);
    return result;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
