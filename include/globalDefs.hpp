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

/** Declaration used by each of the other modules
 */

#ifndef GLOBAL_DEFS
#define GLOBAL_DEFS

#include <cstdint>
#include <fstream>
#include <iostream>
//#include <iterator>
//#include <stdio.h>
#include <string>
//#include <tuple>
#include <cassert>
#include <vector>

namespace deglpf {
//#define DEBUG

using UINT = uint32_t;
using INT = int64_t;
using ENCODED_CHAR = u_int8_t;

const std::string cGENAlphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const std::string cPROTAlphabet = "ACDEFGHIKLMNPQRSTUVWY";
const std::string cDegenerate_PROTAlphabet = "ABCDEFGHIJKLMNPQRSTUVWXYZ";
const std::string cDNAAlphabet = "ACGTU";
const std::string cDegenerate_DNAAlphabet = "ACGTUNRDHKMSWYVB";
const char cDegenerate_symbol_start = '{';
const char cDegenerate_symbol_stop = '}';
const ENCODED_CHAR cMAxUniqueSymbol = 255;

enum class ReturnStatus {
  SUCCESS,
  ERR_ARGS,
  ERR_FILE_OPEN,
  ERR_INVALID_INPUT,
  ERR_INVALID_INDEX,
  ERR_LIMIT_EXCEEDS,
  HELP
};

enum class AlphabetType { DNA, PROT, GEN };

using SEED =
    std::vector<ENCODED_CHAR>; //< A seed is the vector of the encoded character
using SEEDS = std::vector<SEED>;
using DEGENERATE_SYMBOL =
    std::vector<bool>; // 0 or 1 corresponding to each letter of the alphabet in
                       // accordance with its presence (0=> absent, 1=> present)
using DEGENERATE_SYMBOLS = std::vector<DEGENERATE_SYMBOL>;
/** An index structure provides the index to a symbol in a degenerate string (collection of
 * seeds interleaved by the degenerate symbols)
 * **/
using INDEX = struct Index {
  bool is_seed;      // true if index is seed, false for the degenerate symbol
  UINT index;        // index of the seed/degenerate-symbol
  UINT inseed_index; // index within seed (valid only if the type is seed)
};

} // end namespace

#endif
