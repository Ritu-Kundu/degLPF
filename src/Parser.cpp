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

/** Implements class Parser
 */
#include "../include/Parser.hpp"

namespace deglpf {
Parser::Parser(const AlphabetType alphabetType, const std::string &alphabet)
    : _cAlphabetType(alphabetType), _cAlphabet(alphabet) {}

ReturnStatus Parser::parse_sequence(std::ifstream &infile,
                                    Degenerate_string &dgs) const {
  FCheckValidity fCheckValidity = &Parser::is_valid_char_general;
  FMapChar fMapChar = &Parser::map_char_general;

  if (_cAlphabetType == AlphabetType::DNA) {
    fCheckValidity = &Parser::is_valid_char_dna;
    fMapChar = &Parser::map_char_dna;
  } else if (_cAlphabetType == AlphabetType::PROT) {
    fCheckValidity = &Parser::is_valid_char_prot;
    fMapChar = &Parser::map_char_prot;
  }
  SEED seed;
  std::vector<ENCODED_CHAR> symbol;
  bool is_seed_mode = 1; // 0 for seed, 1 for symbol
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) {
      break; // end of this sequece
    }
    for (char c : line) {
      if (isspace(c)) {
        // Ignore
      } else if (c == cDegenerate_symbol_start) {
        seed.shrink_to_fit();
        dgs.add_seed(seed);
        seed.clear();
        is_seed_mode = false;
      } else if (c == cDegenerate_symbol_stop) {
        if (symbol.size() < 2) {
          std::cerr
              << "Invalid Input: Degenerate symbol has less than two letters."
              << std::endl;
          return ReturnStatus::ERR_INVALID_INPUT;
        }
        dgs.add_degenerate_symbol(symbol);
        symbol.clear();
        is_seed_mode = true;
      } else if ((this->*fCheckValidity)(c)) {
        if (is_seed_mode) { // currently collecting seed
          seed.push_back((this->*fMapChar)(c));
        } else { // currently collecting symbol
          symbol.push_back((this->*fMapChar)(c));
        }

      } else {
        std::cerr << "Invalid Input: Not a FASTA format: Invalid character: "
                  << c << std::endl;
        return ReturnStatus::ERR_INVALID_INPUT;
      }
    }
  } // sequence ends
  // Adding the last seed
  seed.shrink_to_fit();
  dgs.add_seed(seed);
  if (dgs.get_size() == 0) {
    std::cerr << "Invalid Input: Empty Sequence." << std::endl;
    return ReturnStatus::ERR_INVALID_INPUT;
  } else if (!is_seed_mode) {
    std::cerr << "Invalid Input: Sequence ended within a degenerate symbol (no "
              << cDegenerate_symbol_stop << " found)." << std::endl;
    return ReturnStatus::ERR_INVALID_INPUT;
  }

  return ReturnStatus::SUCCESS;
}

//////////////////////// private ////////////////////////

bool Parser::is_valid_char_general(const char c) const {
  auto pos = _cAlphabet.find(c);
  if (pos != std::string::npos) {
    return true;
  }
  return false;
}

// Assumes will always be a valid character
ENCODED_CHAR Parser::map_char_general(const char c) const {
  int pos = _cAlphabet.find(c);
  ENCODED_CHAR coded_char = pos + 1;
  return coded_char;
}

bool Parser::is_valid_char_dna(const char c) const {
  bool result = false;
  switch (std::toupper(c)) {
  case 'A':
  case 'C':
  case 'G':
  case 'T':
  case 'U':
    result = true;
  }
  return result;
}

// Assumes will always be a valid character
ENCODED_CHAR Parser::map_char_dna(const char c) const {
  ENCODED_CHAR result;
  switch (std::toupper(c)) {
  case 'A':
    result = 0;
    break;
  case 'C':
    result = 1;
    break;
  case 'G':
    result = 2;
    break;
  case 'T':
    result = 3;
    break;
  case 'U':
    result = 3;
    break;
  }
  return result+1;
}

bool Parser::is_valid_char_prot(const char c) const {
  bool result = false;
  switch (std::toupper(c)) {
  case 'A':
  case 'C':
  case 'D':
  case 'E':
  case 'F':
  case 'G':
  case 'H':
  case 'I':
  case 'K':
  case 'L':
  case 'M':
  case 'N':
  case 'O':
  case 'P':
  case 'Q':
  case 'R':
  case 'S':
  case 'T':
  case 'U':
  case 'V':
  case 'W':
  case 'Y':
    result = true;
    break;
  }
  return result;
}

ENCODED_CHAR Parser::map_char_prot(const char c) const {
  ENCODED_CHAR result;
  switch (std::toupper(c)) {
  case 'A':
    result = 0;
    break;
  case 'C':
    result = 1;
    break;
  case 'D':
    result = 2;
    break;
  case 'E':
    result = 3;
    break;
  case 'F':
    result = 4;
    break;
  case 'G':
    result = 5;
    break;
  case 'H':
    result = 6;
    break;
  case 'I':
    result = 7;
    break;
  case 'K':
    result = 8;
    break;
  case 'L':
    result = 9;
    break;
  case 'M':
    result = 10;
    break;
  case 'N':
    result = 11;
    break;
  case 'O':
    result = 12;
    break;
  case 'P':
    result = 13;
    break;
  case 'Q':
    result = 14;
    break;
  case 'R':
    result = 15;
    break;
  case 'S':
    result = 16;
    break;
  case 'T':
    result = 17;
    break;
  case 'U':
    result = 18;
    break;
  case 'V':
    result = 19;
    break;
  case 'W':
    result = 20;
    break;
  case 'Y':
    result = 21;
    break;
  }
  return result+1;
}

} // end namespace
