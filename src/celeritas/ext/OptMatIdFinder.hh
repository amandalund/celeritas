//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file geocel/OptMatIdFinder.hh
//---------------------------------------------------------------------------//
#pragma once

#include <memory>

#include "geocel/GeoOpticalIdMap.hh"
#include "geocel/Types.hh"
#include "celeritas/geo/GeoFwd.hh"

namespace celeritas
{
class GeantGeoParams;
class VolumeParams;

//---------------------------------------------------------------------------//
/*!
 * Find an optical material ID from a global point.
 */
class OptMatIdFinder
{
  public:
    // Construct from geometry and volume data
    explicit OptMatIdFinder(CoreGeoParams const&);

    // Return the optical ID corresponding to a global point
    OptMatId operator()(Real3 const&) const;

  private:
    CoreGeoParams const& geometry_;
    std::shared_ptr<VolumeParams const> volumes_;
    std::shared_ptr<GeoOpticalIdMap const> geo_to_opt_;
};

//---------------------------------------------------------------------------//
}  // namespace celeritas
