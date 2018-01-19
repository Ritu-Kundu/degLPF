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

/** Implements class Degenerate_string
 */
#include "../include/Degenerate_string.hpp"

namespace deglpf {

Degenerate_string::Degenerate_string(const UINT as)
    : _cAlphabet_size(as), _length(0) {}

void Degenerate_string::add_seed(SEED const &seed) {
  _seeds.push_back(std::move(seed));
  _length += seed.size();
}

void Degenerate_string::add_degenerate_symbol(
    std::vector<ENCODED_CHAR> const &deg) {
  DEGENERATE_SYMBOL ds(_cAlphabet_size + 1, false);
  for (auto l : deg) {
    ds[l] = true;
  }
  _degenerate_symbols.push_back(std::move(ds));
  _degenerate_indices.push_back(_length);
  ++_length;
}

UINT Degenerate_string::get_numberof_seeds() const { return _seeds.size(); }

UINT Degenerate_string::get_size() const { return _length; }

UINT Degenerate_string::get_alphabet_size() const { return _cAlphabet_size; }

const SEEDS &Degenerate_string::get_seeds() const { return _seeds; }

const DEGENERATE_SYMBOLS &Degenerate_string::get_degenerate_symbols() const { return _degenerate_symbols; }

const std::vector<UINT> &Degenerate_string::get_degenerate_indices() const {
  return _degenerate_indices;
}

// Assumes valid Index and non-empty seed
ENCODED_CHAR Degenerate_string::get_seed_lastletter(UINT ind) const {
  assert(ind < _length);
  assert(!_seeds[ind].empty());
  auto last = _seeds[ind].size() - 1;
  return _seeds[ind][last];
}

// Assumes valid Indexes
bool Degenerate_string::is_match(INDEX ind1, INDEX ind2) const {
  assert(ind1.index < _length);
  assert(ind2.index < _length);
  bool result = false;
  if (ind1.is_seed && ind2.is_seed) {
    auto letter1 = _seeds[ind1.index][ind1.inseed_index];
    auto letter2 = _seeds[ind2.index][ind2.inseed_index];
    result = letter1 == letter2;
  } else if (ind1.is_seed && !ind2.is_seed) {
    auto letter1 = _seeds[ind1.index][ind1.inseed_index];
    result = _degenerate_symbols[ind2.index][letter1];
  } else if (!ind1.is_seed && ind2.is_seed) {
    auto letter2 = _seeds[ind2.index][ind2.inseed_index];
    result = _degenerate_symbols[ind1.index][letter2];
  } else if (!ind1.is_seed && !ind2.is_seed) {
    for (auto i = 0; i < _cAlphabet_size; ++i) {
      if (_degenerate_symbols[ind1.index][i] &&
          _degenerate_symbols[ind2.index][i]) {
        result = true;
        break;
      }
    }
  }
  return result;
}

//////////////////////// private ////////////////////////

} // end namespace
