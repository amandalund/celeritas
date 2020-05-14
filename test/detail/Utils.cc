//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file Utils.cc
//---------------------------------------------------------------------------//
#include "Utils.hh"

#include <cstdlib>
#include <string>

namespace
{
bool use_color()
{
    static bool result = (std::getenv("GTEST_COLOR") != nullptr);
    return result;
}
} // namespace

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Get an ANSI color code for [y]ellow / [r]ed / [d]efault.
 */
const char* color_code(char abbrev)
{
    if (!use_color())
        return "";

    switch (abbrev)
    {
        case 'y':
            return "\033[33m";
        case 'r':
            return "\033[31m";
        case ' ':
            // Reset color
            return "\033[0m";
    }

    // Unknown color code: ignore
    return "";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the "skip" message for the skip macro
 */
const char* skip_cstring()
{
    static const std::string str = std::string(color_code('y'))
                                   + std::string("[   SKIP   ]")
                                   + std::string(color_code('d'));
    return str.c_str();
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
