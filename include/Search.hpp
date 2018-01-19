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

/** @file Search.hpp
 * @brief Implements the algorithm.
 * For a given degnerate string (to which it set), it provides the methods to
 * calculate LPF-array using our algorithm and using the naive method as well.
 * Note that here, a solid sequence means the degnerate string where the
 * degenerate symbols have been replaced by unique letters.
 * ASSUMPTION: TOTAL NUMBER OF DEGENERATE SYMBOLS ARE SMALLER THAN 255 - ALPHABET SIZE
 */

#ifndef SEARCH_HPP
#define SEARCH_HPP

#include <list>
#include <algorithm>
#include <sdsl/lcp.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/suffix_arrays.hpp>

#include "Degenerate_string.hpp"
#include "globalDefs.hpp"

namespace deglpf {

class Search {
  /** A SearchDS structure provides the data-structures to make the LCP queries
   * in constant time
  * **/
  struct SearchDS {
    sdsl::csa_bitcompressed<> csa; //<Compreseed suffix array and its inverse
    sdsl::lcp_bitcompressed<> lcp; // < lcp array
    sdsl::rmq_succinct_sct<>
        rmq; // data-structure to answer rmq in constant time
  };

public:
  /** @brief Constructor for Class Search
     * @param dgs reference to the degenerate string for which it will be set
     *
     */
  Search(const Degenerate_string &dgs);

  /** @brief calculates the LPF-array using our algorithm
   * @param lpf reference to the vector in which result will be stored
   *
   **/
  ReturnStatus calculate_lpf(std::vector<UINT> &lpf);

  /** @brief checks whether the given LPF-array is same as would be calculated
   *using the naive approach
   * @param lpf reference to the LPF-array which is to be tested
   *
   **/
  bool naive_test(std::vector<UINT> &lpf) const;
  //////////////////////// private ////////////////////////
private:
  const Degenerate_string &_dgs; //< reference to the degenerate string
  const std::vector<UINT> &_degenerate_indices; //< reference to the positions
                                                // of the degenerate symbols
  const UINT _seq_size;                         // size of the string
  const UINT _k; //< number of the degenerate symbols

  Search::SearchDS
      _fwd_search_ds; //< Search Data-structures for the forward LCP queries
  Search::SearchDS
      _rev_search_ds; //< Search Data-structures for the reverse LCP queries
  /** For each letter of the alphabet, maintain the list of the indices of its
   * occurrence in the reverse sequence. Sorted in descending order wrt
   * reverse..*/
  std::vector<std::list<UINT>> _letter_ind_in_rev;
  /** vector of lpf in the solid sequence (obtained after substituting
   * degenerate symbol with unique letters) */
  std::vector<UINT> _solid_lpf;

  /** Table containing the longest k-lcp (degenerate lcp) beginning at each
   * symbol for each position
   * Number of rows = k; Number of columns = n
   * _longest_degenerate_prefix[i][j] = l => k-lcp of ith deg-symbol and jth
   * position (in solid-sequence) is l
   * The table is initialised to -1 in each cell
   */
  std::vector<std::vector<INT>> _longest_degenerate_prefix;

  /** @brief does the preprocessing:
   *  - Computes the data-structures to answer lcp queries (in constant time) in
   *forward as well as reverse of the solid sequence
   *  - Remembers the occurrences of each letter in the reverse solid sequence
   *  - Calculates the LPF-array of the solid sequence
   *  - Fills the table of longest degenerate match between each symbol and each
   *position
   *
   **/
  ReturnStatus preprocess(std::vector<UINT> &lpf);

  /** @brief Fills the longest degenerate match between given indices in the
   *table
   * - Recursively fills all the cells of the table which are made use of to
   *answer this query
   * - Note that a cell is filled using a single LCP call and making use of
   *results of other cells of the table
   * - Note that call to this function is made only when the cell was set to -1
   * - Assumes the first index to be that of some degenerate symbol
   * @param index_dg reference to the index of some degenerate symbol
   * @param index2 reference to the index of some position (may be a symbol or a
   *solid position in some seed)
   *
   * @see _longest_degenerate_prefix
   *
   **/
  void fill_longest_degenerate_match(const INDEX &index_dg,
                                     const INDEX &index2);

  /** @brief Compute the data-structures to answer lcp queries (in constant
   *time) in forward as well as reverse of the solid sequence
   * - Make the solid sequence and its reverse to set up the corresponding data
   *structures
   * - Also remembers the positions of each letter in the reverse solid sequence
   *to be used by TYPE 2 Search
   * @see _letter_ind_in_rev
   * @see _fwd_search_ds
   * @see _rev_search_ds
   * @see SearchDS
   *
   **/
  ReturnStatus setup_ds();

  /** @brief Calculate the LPF-array of the solid sequence
   * The LPF-Simple algorithm is used (as given in "Computing the Longest
Previous Factor" by
Maxime Crochemore, Lucian Ilie, Costas Iliopoulos, Marcin Kubica, Wojciech
Rytter, Tomasz Waleń)
   * @see _solid_lpf
   *
   **/
  ReturnStatus find_solid_lpf();

  /** @brief Helper to compute the data-structures to answer lcp queries (in
   *constant
   *time) in given sequence
   * Uses SDSL library
   * @param seq reference to the sequence for the data-structures are to be
   *computed
   * @param searchds reference to the structure where result will be stored
   * @see _fwd_search_ds
   * @see _rev_search_ds
   * @see SearchDS
   *
   **/
  void ds_helper(const std::string &seq, Search::SearchDS &searchds);

  /** @brief Answers the k-lcp (longest degenerate match) queries at the given
   *indices in the forward solid sequence
   * @param index1 reference to the first index
   * @param index2 reference to the second index
   * @return the k-lcp value;
   * @see _INDEX
   *
   **/
  UINT find_longest_degenerate_match(const INDEX &index1,
                                     const INDEX &index2) const;

  /** @brief Answers the lcp (exact) queries (in constant
   * time) at the given positions using the given data-structures
   * - Assumption: csa, lcp, rmq are valid (only size is being checked)
   * @param suff the first position
   * @param suff the second position
   * @param searchds reference to the structure containing all the required data
   *structures to be used for answering the queries
   * @return the lcp value; If any index > length of the sequence, returns 0
   * @see _INDEX
   *
   **/
  INT getLCP(const INT suff1, const INT suff2,
             const Search::SearchDS &searchds) const;
};

} // end namespace
#endif
