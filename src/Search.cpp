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

/** Implements class Search
 */
#include "../include/Search.hpp"

namespace deglpf {

Search::Search(const Degenerate_string &dgs)
    : _dgs(dgs),
      _letter_ind_in_rev(dgs.get_alphabet_size() + 1, std::list<UINT>{}),
      _degenerate_indices(dgs.get_degenerate_indices()),
      _seq_size(dgs.get_size()), _k(dgs.get_numberof_seeds() - 1),
      _longest_degenerate_prefix(dgs.get_numberof_seeds() - 1,
                                 std::vector<INT>(dgs.get_size(), -1)),
      _solid_lpf(_seq_size, 0) {}

ReturnStatus Search::calculate_lpf(std::vector<UINT> &lpf) {
  // std::cout << "Calculation started. " << std::endl;
  /* Preprocess */
  preprocess(lpf);

  /* Calculate */
  UINT block = 0;
  bool type2 = false;
  std::vector<UINT> type2_result{};
  for (auto i = 0; i < _seq_size; ++i) {
    if (block < _k && i == _degenerate_indices[block]) { // at degenerate symbol
#ifdef DEBUG
      // PRINTING FOR DEBUGGING
      std::cout << "At Symbol: " << block << std::endl;
#endif
      ++block;
      type2 = false;
      type2_result.clear();
    } else { // in seed
             // Note that we are here as seed is not empty
      auto solid_l = _solid_lpf[i];
      lpf[i] = solid_l;
/* Type 1 Search */
#ifdef DEBUG
      // PRINTING FOR DEBUGGING
      std::cout << "Type 1 at : i L: " << i << " " << solid_l << std::endl;
#endif
      for (auto j = 0; j < block; ++j) { // for each previous block (or seed)
        auto stop_pos = _degenerate_indices[j];
        // Check L-region in jth block (seed)
        auto first_pos = (j == 0) ? (0) : (_degenerate_indices[j - 1] + 1);
        UINT temp=  stop_pos - solid_l;
        if (static_cast<INT>(stop_pos) - static_cast<INT>(solid_l) < 0) {
          temp = 0;
        }
        INT start_pos = std::max(temp, first_pos);
#ifdef DEBUG
        // PRINTING FOR DEBUGGING
        std::cout << "Checking in block : block start_pos stop_pos: " << j << " " << start_pos
                  << " " << stop_pos << std::endl;
#endif
        for (auto pos = start_pos; pos < stop_pos; ++pos) {
          auto lcp = getLCP(pos, i, _fwd_search_ds);
          auto match_upto_pos = i + lcp;
          if ((pos + lcp) == stop_pos &&
              (match_upto_pos !=
               _seq_size)) { // Prefix of L is suffix of this seed
            UINT possible_lpf =
                lcp + _longest_degenerate_prefix[j][match_upto_pos];
#ifdef DEBUG
            // PRINTING FOR DEBUGGING
            std::cout << "Candiadte : pos possible_lpf: " << pos << " "
                      << possible_lpf << std::endl;
#endif
            lpf[i] = std::max(lpf[i], possible_lpf);
          }
        } // Checked L-region
        // Check jth symbol
        lpf[i] = std::max(lpf[i],
                          static_cast<UINT>(_longest_degenerate_prefix[j][i]));
      } // Checked each block(seed)

      /* Type 2 Search, if needed */
      if (block < _k && (i + solid_l == _degenerate_indices[block])) {
#ifdef DEBUG
        // PRINTING FOR DEBUGGING
        std::cout << "Type 2 at : i  " << i << std::endl;
#endif
        // Note that we will not be here for i=0 (as solid_l will be 0)
        if (!type2) { // enter into type 2 mode for the first time
          std::vector<UINT> type2_local(solid_l + 1, 0);

          ENCODED_CHAR c = _dgs.get_seed_lastletter(block);
          auto following_symb_pos = _degenerate_indices[block];
          auto rev_i = _seq_size - i - 1;
          auto rev_last_pos = _seq_size - following_symb_pos;
          auto rev_stop_pos =
              (block == 0) ? (_seq_size)
                           : (_seq_size - 1 - _degenerate_indices[block - 1]);
          // for each occurrence (succeeding) of letter in reverse, find
          // potential longer lpf
          for (auto p : _letter_ind_in_rev[c]) {
            if (p > rev_last_pos) {
              auto rev_lcp = 1;
              if ((rev_last_pos + 1) < rev_stop_pos &&
                  (p + 1 < _seq_size)) { // Take rev-lpf if  there are solid
                // letters preceeding it in the seed
                rev_lcp += getLCP(p + 1, rev_last_pos + 1, _rev_search_ds);
              }
              // find tail of the match
              auto tail_match = 0;
              auto reverse_next_p = _seq_size - p;
              if (reverse_next_p < _seq_size) {
                tail_match = _longest_degenerate_prefix[block][reverse_next_p];
              }
              UINT potential_lpf = rev_lcp + tail_match;

              // update the lpf-value of the corresponding length
              type2_local[rev_lcp] =
                  std::max(type2_local[rev_lcp], potential_lpf);
            } else {
              break;
            }
          }
          type2 = true;
          type2_result = std::move(type2_local);
#ifdef DEBUG
          // PRINTING FOR DEBUGGING
          std::cout << "L-Table : " << std::endl;
          for (auto c : type2_result) {
            std::cout << c << " ";
          }
          std::cout << std::endl;
#endif
        }
        // Use the stored result to find answer
        UINT potential_type2_lpf = (lpf[i - 1] == 0) ? (0) : (lpf[i - 1] - 1);
        potential_type2_lpf =
            std::max(type2_result[solid_l], potential_type2_lpf);
#ifdef DEBUG
        // PRINTING FOR DEBUGGING
        std::cout << "Candiadte : Type2 possible_lpf: " << potential_type2_lpf
                  << std::endl;
#endif
        lpf[i] = std::max(lpf[i], potential_type2_lpf);
      }
    }
#ifdef DEBUG
    // PRINTING FOR DEBUGGING
    std::cout << "Final ans: " << lpf[i] << std::endl;
#endif
  } // Filled each position

  // std::cout << "Search completed. " << std::endl;
  return ReturnStatus::SUCCESS;
}

bool Search::naive_test(std::vector<UINT> &lpf) const {
  std::cout << "NAIVE TESTING: ";
  bool result = true;
  UINT block_i = 0;
  UINT base_i = 0;
  for (auto i = 0; i < _seq_size; ++i) {
    INDEX index1{true, block_i, i - base_i}; // Assuming in seed
    if (block_i < _k &&
        i == _degenerate_indices[block_i]) { // at degenerate symbol
      index1 = INDEX{false, block_i, 0};     // Index of the seed position
      ++block_i;
      base_i = i + 1;
    }
    UINT longest_match = 0;
    UINT block_j = 0;
    UINT base_j = 0;
    for (auto j = 0; j < i; ++j) {
      INDEX index2{true, block_j, j - base_j}; // Index of the seed position
      if (block_j < _k &&
          j == _degenerate_indices[block_j]) { // at degenerate symbol
        index2 = INDEX{false, block_j, 0};     // Index of the symbol
        ++block_j;
        base_j = j + 1;
      }
      auto l = find_longest_degenerate_match(index1, index2);
      longest_match = std::max(longest_match, l);
    }
    if (lpf[i] != longest_match) {
      std::cerr << "ERROR: WRONG ANSWER (" << lpf[i] << ") at position: " << i
                << " : RIGHT ANS = " << longest_match << std::endl;
      result = false;
      break;
    }
  }
  if (result) {
    std::cout << " PASS\n";
  }
  else {
    std::cout << " FAIL\n";
  }
  return result;
}

//////////////////////// private ////////////////////////
ReturnStatus Search::preprocess(std::vector<UINT> &lpf) {
  /* Set-up the data-structures */
  setup_ds();

  /* Find the solid-lpf for each position */
  find_solid_lpf();

  /* Fill the table of the longest k-lcp at each symbol and each position */
  for (UINT symb = 0; symb < _k; ++symb) {
    auto symb_pos = _degenerate_indices[symb];
    INDEX index1{false, symb, 0}; // Index of the symbol
    _longest_degenerate_prefix[symb][symb_pos] =
        0; // The lcp of a symbol with itself is set to 0
    UINT block = 0;
    for (auto i = 0; i < _seq_size; ++i) {
      bool is_deg = false;
      if (block < _k &&
          i == _degenerate_indices[block]) { // at degenerate symbol
        ++block;
        is_deg = true;
      }
      if (_longest_degenerate_prefix[symb][i] ==
          -1) { // The cell is uninitalised
        INDEX index2{};
        if (is_deg) {                          // at degenerate symbol
          index2 = INDEX{false, block - 1, 0}; // Index of the symbol
          fill_longest_degenerate_match(index1, index2);
        } else { // in seed
          auto base = (block == 0) ? (0) : (_degenerate_indices[block - 1] + 1);
          index2 = INDEX{true, block, i - base}; // Index of the symbol
        }
        fill_longest_degenerate_match(index1, index2);
      } // This cell filled

      auto k_lcp = _longest_degenerate_prefix[symb][i];
      if (i < symb_pos &&
          k_lcp >
              lpf[symb_pos]) { // It influences the final LPF for this symbol
        lpf[symb_pos] = k_lcp;
      }
    } // each position done
  }   // each symbol done
#ifdef DEBUG
  // PRINTING FOR DEBUGGING
  for (int j = 0; j < _k; ++j) {
    std::cout << "TABLE: " << j << std::endl;
    for (auto i = 0; i < _seq_size; ++i) {
      std::cout << _longest_degenerate_prefix[j][i] << " ";
    }
    std::cout << std::endl;
  }
#endif
  return ReturnStatus::SUCCESS;
}

// Assumes index_dg is always for a degenerate symbol
void Search::fill_longest_degenerate_match(const INDEX &index_dg,
                                           const INDEX &index2) {
  // Map indices into sequence position
  auto symb_ind = index_dg.index;
  auto symb_pos = _degenerate_indices[symb_ind];
  UINT pos2 = 0;
  UINT pos2_next_symb_ind = index2.index;
  if (index2.is_seed) { // solid position
    auto block = index2.index;
    auto base = (block == 0) ? (0) : (_degenerate_indices[block - 1] + 1);
    pos2 = base + index2.inseed_index;
  } else { // degenerate symbol
    pos2 = _degenerate_indices[index2.index];
    ++pos2_next_symb_ind;
  }
  /* Find match */
  INT longest_match = 0;
  // As at least one symbol is degenerate, ask for approx match at this position
  if (_dgs.is_match(index_dg, index2)) { // these positions match; extend match
    longest_match = 1;                   // match is at least 1
    // If any of the positions exceeds the size, lcp will be returned as 0
    auto lcp = getLCP(symb_pos + 1, pos2 + 1, _fwd_search_ds);
    longest_match += lcp;
    UINT new_pos1 = symb_pos + longest_match;
    UINT new_pos2 = pos2 + longest_match;
    UINT next_stop_pos1 =
        (symb_ind < _k - 1) ? (_degenerate_indices[symb_ind + 1]) : (_seq_size);
    UINT next_stop_pos2 = (pos2_next_symb_ind < _k)
                              ? (_degenerate_indices[pos2_next_symb_ind])
                              : (_seq_size);
    // any new position goes outside string or both are in seed, we are done
    // otherwise, add the result in the cell of the new positions.
    bool new_pos1_deg = (new_pos1 == next_stop_pos1);
    bool new_pos2_deg = (new_pos2 == next_stop_pos2);
    if ((new_pos1 < _seq_size) && (new_pos2 < _seq_size) &&
        (new_pos1_deg || new_pos2_deg)) {
      UINT new_symb_ind =
          (new_pos1_deg) ? (symb_ind + 1) : (pos2_next_symb_ind);
      UINT pos = (new_pos1_deg) ? (new_pos2) : (new_pos1);
      if (_longest_degenerate_prefix[new_symb_ind][pos] ==
          -1) {                                   // check the cell
        INDEX new_index1{false, new_symb_ind, 0}; // Index of the symbol
        INDEX new_index2;
        if (new_pos1_deg) {   // Index1 created from pos1
          if (new_pos2_deg) { // Index2 is degenerate
            new_index2 = INDEX{false, pos2_next_symb_ind, 0};
          } else { // Index2 is in seeed
            auto base = (pos2_next_symb_ind == 0)
                            ? (0)
                            : (_degenerate_indices[pos2_next_symb_ind - 1] + 1);
            new_index2 = INDEX{true, pos2_next_symb_ind, new_pos2 - base};
          }
        } else {              // Index1 created from pos2
          if (new_pos1_deg) { // Index2 is degenerate
            new_index2 = INDEX{false, symb_ind + 1, 0};
          } else { // Index2 is in seeed
            auto base = (_degenerate_indices[symb_ind] + 1);
            new_index2 = INDEX{true, symb_ind + 1, new_pos1 - base};
          }
        }
        fill_longest_degenerate_match(new_index1, new_index2);
      }
      longest_match += _longest_degenerate_prefix[new_symb_ind][pos];
    }
  }
  /* Fill the cell/s */
  _longest_degenerate_prefix[symb_ind][pos2] = longest_match;
  // If the second position is also degenerate, fill the corresponding cell
  // It will definitely be -1; otherwise the second symbol would have already
  // filled this cell
  // And we wouldn't have been in this call.
  if (!index2.is_seed) {
    _longest_degenerate_prefix[index2.index][symb_pos] = longest_match;
  }
#ifdef DEBUG
  // PRINTING FOR DEBUGGING
  std::cout << "FILLED CELL: symb_ind pos: " << symb_ind << "  " << pos2
            << std::endl;
#endif
}

ReturnStatus Search::setup_ds() {
  const SEEDS &seeds = _dgs.get_seeds();
  UINT num_seeds = seeds.size();
  UINT delimiter = _dgs.get_alphabet_size() + 1;
  // Combine all seeds, replacing degenerate symbols with distinct unique
  // symbols (not in alphabet)
  std::string seq(_seq_size, 0);
  // Take reverse of the combined sequence for the reverse LCP queries
  std::string rev_seq(_seq_size, 0);
  UINT fwd = 0;
  UINT rev = _seq_size - 1;
  for (auto ind = 0; ind < num_seeds; ++ind) {
    for (auto c : seeds[ind]) {
      seq[fwd++] = static_cast<unsigned char>(c);
      _letter_ind_in_rev[c].push_back(
          rev); // Remember the indices of occurrences of each character
      rev_seq[rev--] = static_cast<unsigned char>(c);
    }
    if (ind < num_seeds - 1) { // Seed followed by lambda_i except the last seed
      if (delimiter > cMAxUniqueSymbol) {
        std::cerr << "Invalid Input: Sequence has more than allowed degenerate "
                     "symbols.: (allowed: "
                  << cMAxUniqueSymbol << " found: " << delimiter << std::endl;
      }

      seq[fwd++] = static_cast<unsigned char>(delimiter);
      rev_seq[rev--] = static_cast<unsigned char>(delimiter);
      ++delimiter;
    }
  }
#ifdef DEBUG
  // PRINTING FOR DEBUGGING
  std::cout << "SOLID SEQUENCE: \n";
  for (auto c : seq) {
    std::cout << (int)c << " ";
  }
  std::cout << std::endl;
#endif
  ds_helper(seq, _fwd_search_ds);
  ds_helper(rev_seq, _rev_search_ds);
  return ReturnStatus::SUCCESS;
}

ReturnStatus Search::find_solid_lpf() {
  std::vector<INT> prev(_seq_size + 1, 0);
  std::vector<INT> next(_seq_size + 1, 0);
  std::vector<UINT> lcp(_seq_size + 1, 0);
  for (auto i = 0; i < _seq_size; ++i) {
    lcp[i] = _fwd_search_ds.lcp[i];
    prev[i] = i - 1;
    next[i] = i + 1;
  }
  auto r = 0;
  for (int i = _seq_size - 1; i >= 0; --i) {
    r = _fwd_search_ds.csa.isa[i];
    _solid_lpf[i] = std::max(lcp[r], lcp[next[r]]);
    lcp[next[r]] = std::min(lcp[r], lcp[next[r]]);
    if (prev[r] >= 0) {
      next[prev[r]] = next[r];
    }
    if (next[r] < _seq_size) {
      prev[next[r]] = prev[r];
    }
  }
  return ReturnStatus::SUCCESS;
}

void Search::ds_helper(const std::string &seq, Search::SearchDS &searchds) {
  // sdsl::construct_im(searchds.csa, patternstr, 1); // 1 for alphabet type
  // std::cout << " i SA ISA T[SA[i]..SA[i]-1]" << std::endl;
  // sdsl::csXprintf(std::cout, "%2I %2S %3s %:3T", csa);
  // qsufsort::construct_sa(csa, );
  sdsl::construct_im(searchds.csa, seq, 1); // 1 for alpahabet type
  sdsl::construct_im(searchds.lcp, seq, 1); // 1 for alphabet type
  searchds.rmq = std::move(sdsl::rmq_succinct_sct<>(&(searchds.lcp)));
  // RMQ rmq(&(searchds.lcp));
  // rmq does not need its arg to answer the queries
  // sdsl::util::clear(lcp); // so we can free the space for v
}

UINT Search::find_longest_degenerate_match(const INDEX &index1,
                                           const INDEX &index2) const {
  // Map indices into sequence positions
  UINT pos1 = 0;
  UINT pos1_next_symb_ind = index1.index;
  if (index1.is_seed) { // solid position
    auto block = index1.index;
    auto base = (block == 0) ? (0) : (_degenerate_indices[block - 1] + 1);
    pos1 = base + index1.inseed_index;
  } else { // degenerate symbol
    pos1 = _degenerate_indices[index1.index];
    ++pos1_next_symb_ind;
  }
  UINT pos2 = 0;
  UINT pos2_next_symb_ind = index2.index;
  if (index2.is_seed) { // solid position
    auto block = index2.index;
    auto base = (block == 0) ? (0) : (_degenerate_indices[block - 1] + 1);
    pos2 = base + index2.inseed_index;
  } else { // degenerate symbol
    pos2 = _degenerate_indices[index2.index];
    ++pos2_next_symb_ind;
  }
  /* Find match */
  UINT longest_match = 0;
  // if any of the symbol is degenerate, first test the letters
  if (!index1.is_seed || !index2.is_seed) {
    if (_dgs.is_match(index1, index2)) { // these positions match; extend match
      longest_match += 1;                // match is at least 1
      ++pos1;
      ++pos2;
    } else { // no match at this position
      return longest_match;
    }
  }

  while (pos1 < _seq_size && pos2 < _seq_size) {
    // If any of the positions exceeds the size, lcp will be returned as 0
    auto lcp = getLCP(pos1, pos2, _fwd_search_ds);
    longest_match += lcp;
    pos1 += lcp;
    pos2 += lcp;
    auto next_stop_pos1 = (pos1_next_symb_ind < _k)
                              ? (_degenerate_indices[pos1_next_symb_ind])
                              : (_seq_size);
    auto next_stop_pos2 = (pos2_next_symb_ind < _k)
                              ? (_degenerate_indices[pos2_next_symb_ind])
                              : (_seq_size);
    // any of the new positions goes outside string or both are in seed, we are
    // done; otherwise, continue extending
    bool is_pos1_deg = (pos1 == next_stop_pos1);
    bool is_pos2_deg = (pos2 == next_stop_pos2);
    if ((pos1 < _seq_size) && (pos2 < _seq_size) &&
        (is_pos1_deg || is_pos2_deg)) {
      // Get the INDEX structutres from the pos (for checking match)
      INDEX new_index1;
      if (is_pos1_deg) { // Index1 is degenerate
        new_index1 = INDEX{false, pos1_next_symb_ind, 0};
        ++pos1_next_symb_ind;
      } else { // Index1 is in seeed
        auto base = (pos1_next_symb_ind == 0)
                        ? (0)
                        : (_degenerate_indices[pos1_next_symb_ind - 1] + 1);
        new_index1 = INDEX{true, pos1_next_symb_ind, pos1 - base};
      }
      INDEX new_index2;
      if (is_pos2_deg) { // Index2 is degenerate
        new_index2 = INDEX{false, pos2_next_symb_ind, 0};
        ++pos2_next_symb_ind;
      } else { // Index2 is in seeed
        auto base = (pos2_next_symb_ind == 0)
                        ? (0)
                        : (_degenerate_indices[pos2_next_symb_ind - 1] + 1);
        new_index2 = INDEX{true, pos2_next_symb_ind, pos2 - base};
      }
      // As at least one symbol is degenerate, ask for approx match at this
      // position
      if (_dgs.is_match(new_index1,
                        new_index2)) { // these positions match; extend match
        longest_match += 1;            // match is at least 1
        ++pos1;
        ++pos2;
      } else { // no match at this position
        break;
      }
    } else { // both real mismatches or hit the end of the string
      break;
    }
  }
  return longest_match;
}

// Assumption: csa, lcp, rmq are valid (only size is being checked)
// If any index > length of the sequence, returns 0
INT Search::getLCP(const INT suff1, const INT suff2,
                   const Search::SearchDS &searchds) const {
  assert(!searchds.csa.empty());
  assert(!searchds.lcp.empty());
  assert(searchds.rmq.size() > 0);
  if (suff1 >= _seq_size || suff2 >= _seq_size) {
    return 0;
  }
  INT l_rmq, r_rmq;
  INT r1 = searchds.csa.isa[suff1];
  INT r2 = searchds.csa.isa[suff2];
  if (r1 < r2) {
    l_rmq = r1;
    r_rmq = r2;
  } else {
    l_rmq = r2;
    r_rmq = r1;
  }
  INT minIndex = searchds.rmq(l_rmq + 1, r_rmq);
  INT lcp = searchds.lcp[minIndex];
  return lcp;
  // std::cout << "LCP: "<< lcp << std::endl;
}

} // end namespace
