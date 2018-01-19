/**
    degLPF: Computing Longest Previous Factor (LPF) Array in a Degenerate String
    Copyright (C) 2018 Ritu Kundu, Fatima Vayani, and Steven Watts
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

/** @file utilDefs.hpp
 * @brief Defines the utility functions related to the usage of the tool
 * It provides methods for processing the input flags and print the help for the tool.
 */


#ifndef UTIL_DEFS
#define UTIL_DEFS

#include <sys/time.h>
#include <getopt.h>
#include <cctype>
#include <cassert>

#include "globalDefs.hpp"

namespace deglpf{
struct InputFlags{
  std::string input_filename;
  std::string output_filename;
  AlphabetType alphabet_type;
};

void usage (void);
ReturnStatus decodeFlags(int argc, char* argv [], struct InputFlags& flags);

} // end namespace

#endif
