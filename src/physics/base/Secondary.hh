//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file Secondary.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/Types.hh"
#include "ParticleDef.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * New particle created via an Interaction.
 *
 * It will be converted into a "track initializer" using the parent track's
 * information.
 */
struct Secondary
{
    ParticleDefId    def_id;    //!< New particle type
    units::MevEnergy energy;    //!< New kinetic energy
    Real3            direction; //!< New direction

    // Secondary did not survive cutoffs and has been processed
    static inline CELER_FUNCTION Secondary from_cutoff();

    // Whether the secondary survived cutoffs
    explicit inline CELER_FUNCTION operator bool() const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct a secondary that has been killed after having cutoffs applied and
 * energy deposited locally
 */
CELER_FUNCTION Secondary Secondary::from_cutoff()
{
    Secondary result;
    result.def_id = {};
    result.energy = zero_quantity();
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Whether the Secondary succeeded.
 */
CELER_FUNCTION Secondary::operator bool() const
{
    return static_cast<bool>(this->def_id);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
