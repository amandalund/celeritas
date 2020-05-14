//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file Main.hh
//---------------------------------------------------------------------------//
#ifndef test_Main_hh
#define test_Main_hh

#include "detail/TestMain.hh"

#include <string>
#include <vector>

using std::cout;
using std::endl;

//! Define main
int main(int argc, char** argv)
{
    return celeritas::detail::test_main(argc, argv);
}

#endif // test_Main_hh
