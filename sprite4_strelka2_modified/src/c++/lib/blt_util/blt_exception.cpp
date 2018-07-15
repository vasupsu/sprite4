//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///
#include "blt_util/blt_exception.hh"

#ifdef KILL_EXCEPTIONS
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>
#endif



blt_exception::
blt_exception(const char* s)
    : message(s)
{
#ifdef KILL_EXCEPTIONS
    log_os << "ERROR:: " << s << std::endl;
    abort();
#endif
}

