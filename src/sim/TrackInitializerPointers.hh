//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file TrackInitializerPointers.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/Types.hh"
#include "TrackInitializer.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * View to the data used to initialize new tracks.
 */
struct TrackInitializerPointers
{
    span<TrackInitializer> initializers;
    span<size_type>        parent;
    span<size_type>        vacancies;
    span<size_type>        secondary_counts;
    size_type              track_count;

    explicit CELER_FUNCTION operator bool() const
    {
        return !secondary_counts.empty();
    }
};

//---------------------------------------------------------------------------//
} // namespace celeritas
