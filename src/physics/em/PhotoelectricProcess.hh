//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file PhotoelectricProcess.hh
//---------------------------------------------------------------------------//
#pragma once

#include "physics/base/Process.hh"

#include "physics/base/ImportedProcessAdapter.hh"
#include "physics/base/ParticleParams.hh"
#include "physics/em/AtomicRelaxationParams.hh"
#include "physics/material/MaterialParams.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Photoelectric effect process for gammas.
 */
class PhotoelectricProcess : public Process
{
  public:
    //!@{
    //! Type aliases
    using SPConstParticles   = std::shared_ptr<const ParticleParams>;
    using SPConstMaterials   = std::shared_ptr<const MaterialParams>;
    using SPConstImported    = std::shared_ptr<const ImportedProcesses>;
    using SPConstAtomicRelax = std::shared_ptr<const AtomicRelaxationParams>;
    //!@}

  public:
    // Construct from Livermore photoelectric data
    PhotoelectricProcess(SPConstParticles particles,
                         SPConstMaterials materials,
                         SPConstImported  process_data);

    // Construct from Livermore data and EADL atomic relaxation data
    PhotoelectricProcess(SPConstParticles   particles,
                         SPConstMaterials   materials,
                         SPConstImported    process_data,
                         SPConstAtomicRelax atomic_relaxation,
                         size_type          vacancy_stack_size);

    // Construct the models associated with this process
    VecModel build_models(ModelIdGenerator next_id) const final;

    // Get the interaction cross sections for the given energy range
    StepLimitBuilders step_limits(Applicability range) const final;

    // Name of the process
    std::string label() const final;

  private:
    SPConstParticles       particles_;
    SPConstMaterials       materials_;
    ImportedProcessAdapter imported_;
    SPConstAtomicRelax     atomic_relaxation_;
    size_type              vacancy_stack_size_{};
};

//---------------------------------------------------------------------------//
} // namespace celeritas
