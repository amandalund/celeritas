//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file testCommunicator.cc
//---------------------------------------------------------------------------//
#include "comm/Communicator.hh"

#include "../Main.hh"
#include "../Test.hh"

using celeritas::Communicator;

//---------------------------------------------------------------------------//
// TEST HARNESS
//---------------------------------------------------------------------------//

class CommunicatorTest : public celeritas::Test
{
  protected:
    void SetUp() { cout << "Hello!" << endl; }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CommunicatorTest, rank)
{
    Communicator comm;

    EXPECT_EQ(MPI_COMM_WORLD, comm.mpi_comm());
#ifdef CELERITAS_USE_MPI
    int expected_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &expected_rank);
    EXPECT_EQ(expected_rank, comm.rank());

    int expected_size;
    MPI_Comm_size(MPI_COMM_WORLD, &expected_size);
    EXPECT_EQ(expected_size, comm.size());
#endif

    comm.barrier();
}
