//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file NonMasterResultPrinter.hh
//---------------------------------------------------------------------------//
#ifndef test_detail_NonMasterResultPrinter_hh
#define test_detail_NonMasterResultPrinter_hh

#include <gtest/gtest.h>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Print test results on non-rank-zero processes.
 */
class NonMasterResultPrinter : public ::testing::EmptyTestEventListener
{
  public:
    // Construct with MPI rank
    explicit NonMasterResultPrinter(int rank);

    void OnTestPartResult(const ::testing::TestPartResult& result) override;

  private:
    // >>> DATA
    int rank_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

#endif // test_detail_NonMasterResultPrinter_hh
