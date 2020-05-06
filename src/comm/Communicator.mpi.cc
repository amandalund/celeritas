//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file Communicator.mpi.cc
//---------------------------------------------------------------------------//
#include "Communicator.hh"

#include "base/Assert.hh"
#include "ScopedMpiInit.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with a native MPI communicator
 */
Communicator::Communicator(MpiComm comm) : comm_(comm), rank_(-1), size_(-1)
{
    REQUIRE(ScopedMpiInit::initialized());

    // Save rank and size
    int err;
    err = MPI_Comm_rank(comm_, &rank_);
    CHECK(err == MPI_SUCCESS);
    err = MPI_Comm_size(comm_, &size_);
    CHECK(err == MPI_SUCCESS);

    ENSURE(this->rank() >= 0 && this->rank() < this->size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait for all processes in this communicator to reach the barrier
 */
void Communicator::barrier() const
{
    int err = MPI_Barrier(comm_);
    CHECK(err == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
