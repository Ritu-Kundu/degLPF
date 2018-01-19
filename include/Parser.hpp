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

/** @file Parser.hpp
 * @brief Defines the class Subgraph.
 * It provides methods for parsing Fasta files.
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include "Degenerate_string.hpp"
#include "globalDefs.hpp"

namespace deglpf {
/** Class Parser
 * A Parser contains the method for parsing a sequence into the encoded sequence of integers.
 *
 */
class Parser {
  // Type of function for checking validity of a character
  using FCheckValidity = bool (Parser::*)(const char) const;
  // Type of function for mapping a valid character to encoded character
  using FMapChar = ENCODED_CHAR (Parser::*)(const char) const;

public:
  /** @brief Constructor for Class Parser
   * @param at Alphabet type //DNA for genomic sequences, PROT for proteins
   * @param alphabet string containing all the valid letters of the alphabet
   * @see AlphabetType
   * @see cDegenerate_PROTAlphabet
   * @see cDegenerate_DNAAlphabet
   */
  Parser(const AlphabetType alphabetType, const std::string &alphabet);

  /** @brief parses sequence into encoded sequence of integers.
 * Length of one integer defined by ENCODED_CHAR: Currently 1 byte long
 * Input sequence should be in Fasta format: It returns only one sequence
 * Ignores spaces between charcters.
 * Degenerate symbol:
 *  - It starts from cDegenerate_symbol_start (currently '{')
 *  - It stops at cDegenerate_symbol_stop (currently '}')
 *  - Letters within a symbol may be separated by spaces or may not be separated
 at all.
 *  - Letters may be repeated but counted as one.
 *  - At least two letters should be present.
 * Returns on encountering an invalid character or if a symbol did not close or
 did not have at least two letters.
 * @param infile handle of the file (currently pointing at the beginning of the
 sequence)
 * @param dgs reference to the degenerate string to be set up from the read
 sequence
 * @return execution status // SUCCESS if input is valid, otherwise
 corresponding error code after logging the error
 **/
  ReturnStatus parse_sequence(std::ifstream &infile,
                              Degenerate_string &dgs) const;

private:
  const AlphabetType _cAlphabetType; ///< Type of alphabet: DNA, PROT or GEN
  const std::string _cAlphabet;      ///< Original alphabet

  /** @brief checks whether the given character is valid in the alphabet for
  which the parser is set
   *
   **/
  bool is_valid_char_general(const char c) const;

  /** @brief encodes the given character
   * Maps the character to its index in the alphabet for which the parser is set
   * Assumes the character to be valid always
   *
   **/
  ENCODED_CHAR map_char_general(const char c) const;

  /** @brief checks whether the given character is valid in the DNA alphabet
   *
   * @see cDegenerate_DNAAlphabet
   **/
  bool is_valid_char_dna(const char c) const;

  /** @brief encodes the given DNA character
 * Maps the character to its index in the alphabet for which the parser is set
 * Assumes the character to be valid always
 * 
 * @see cDegenerate_DNAAlphabet
 **/
  ENCODED_CHAR map_char_dna(const char c) const;

  /** @brief checks whether the given character is valid in the Protein alphabet
   *
   * @see cDegenerate_PROTAlphabet
   **/
  bool is_valid_char_prot(const char c) const;

  /** @brief encodes the given Protein character
   * Maps the character to its index in the alphabet for which the parser is set
   * Assumes the character to be valid always
   * 
   * @see cDegenerate_PROTAlphabet
   **/
  ENCODED_CHAR map_char_prot(const char c) const;
};

} // end namespace
#endif
