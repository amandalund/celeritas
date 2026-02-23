//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file geocel/OptMatIdFinder.cc
//---------------------------------------------------------------------------//
#include "OptMatIdFinder.hh"

#include "corecel/Assert.hh"
#include "geocel/GeantGeoParams.hh"
#include "geocel/VolumeParams.hh"
#include "celeritas/geo/CoreGeoParams.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from geometry and volume data.
 */
OptMatIdFinder::OptMatIdFinder(CoreGeoParams const& geometry)
    : geometry_(geometry), volumes_(celeritas::global_volumes().lock())
{
    CELER_VALIDATE(volumes_, << "global canonical volumes are not loaded");

    auto geant_geo = celeritas::global_geant_geo().lock();
    CELER_VALIDATE(geant_geo, << "global Geant4 geometry is not loaded");

    geo_to_opt_ = geant_geo->geo_opt_id_map();
}

//---------------------------------------------------------------------------//
/*!
 * Return the optical ID corresponding to a global point.
 */
OptMatId OptMatIdFinder::operator()(Real3 const& global_point) const
{
    auto vi_id = geometry_.find_volume_instance_at(global_point);
    auto mat_id = volumes_->material(volumes_->volume(vi_id));
    return (*geo_to_opt_)[mat_id];
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
