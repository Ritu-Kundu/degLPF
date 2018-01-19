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

/** @file Degenerate_string.hpp
 * @brief Defines the class Degenerate String.
 * It represents a degenerate string: Collection of seeds interleaved by the
 * degenerate symbols.
 * It must begin and end with a seed; the seed itself can be empty.
 * - the nuber of degenerate symbols = number of seeds - 1
 * The degenerate symbol is represented as a boolean vector of size equal to
 * alphabet. 0 if corresponding letter is absent else 1.
 * All the degenefrate symbols are collected together as a vector; all the seeds
 * are collected together as a vector.
 * It also contains the indices at which the degenerate symbol appears in the
 * string when seen as a sequence.
 *
 * Provides methods for the following:
 * - Adding seed.
 * - Adding a degenerate symbol.
 * - Getting a reference to the seeds, the total number of the seeds, the
 * total size of the string seen as a sequence, and the size of the alphabet
 * used.
 * - Checks whether the two symbols at the given indices match (degeneate
 * match).
 */

#ifndef DEGENERATESTRING_HPP
#define DEGENERATESTRING_HPP

#include "globalDefs.hpp"
#include "utilDefs.hpp"
// TODO: Add SNP support in pattern; map char to int
namespace deglpf {

class Degenerate_string {

public:
  /** @brief Constructor for Class Degenerate_string
   * @param as Alphabet-size. It implies that characters in the seeds or in
   * symbols will be from 1 to as.
   * @param as size of the alphabet
   *
   */
  Degenerate_string(const UINT as);

  /** @brief adds the given seed in the collection and changes the length of the
   *sequence correspondingly
   * @see SEED
   *
   **/
  void add_seed(SEED const &seed);

  /** @brief adds the given degenerate symbol in the collection and increases
   *the length of the sequence by 1.
   * It also updates the index of the symbol when the string is seens as the
   *sequence.
   * @param deg reference to the vector containing the letters present in the
   *symbol
   * @see DEGENERATE_SYMBOL
   *
   **/
  void add_degenerate_symbol(std::vector<ENCODED_CHAR> const &deg);

  /** @brief returns the number of the seeds in the collection
   *
   **/
  UINT get_numberof_seeds() const;

  /** @brief returns the length (size) of the string when seen as a sequence
   *
   **/
  UINT get_size() const;

  /** @brief returns the size of the alphabet used
   *
   **/
  UINT get_alphabet_size() const;

  /** @brief returns a reference to the collection of the seeds
   *
   **/
  const SEEDS &get_seeds() const;

  /** @brief returns a reference to the collection of the degenerate symbols
   *
   **/
  const DEGENERATE_SYMBOLS &get_degenerate_symbols() const;

  /** @brief returns a reference to the vector of indices of the degenerate
    *symbols in the string when seen as a sequence
   *
   **/
  const std::vector<UINT> &get_degenerate_indices() const;

  /** @brief returns the last symbol in the seed at the given index
   * Assumes the index is valid(i.e. it is in a seed) and that seed is non-empty
   * @see Index
   *
   **/
  ENCODED_CHAR get_seed_lastletter(UINT ind) const;

  /** @brief returns whether the symbol at the given indices match (degenerate)
   * The match is degnerate.
   * Assumes Indexes to be valid.
   * @see Index
   *
   **/
  bool is_match(INDEX ind1, INDEX ind2) const;

  //////////////////////// private ////////////////////////
private:
  const UINT _cAlphabet_size; //< Size of the alphabet
  SEEDS _seeds;               //< collection of all the seeds
                              /** Vector of the degenerate symbols
                               * A degenerate symbol has length 1 more than the alphabet size as the letters
                               * are being mapped from 1 to alphabet-size.
                              */
  DEGENERATE_SYMBOLS _degenerate_symbols;
  std::vector<UINT> _degenerate_indices; //< Vector of indices of the degenerate
                                         // symbols in the string
  UINT _length;                          //< Total length of the string
};

} // end namespace
#endif
