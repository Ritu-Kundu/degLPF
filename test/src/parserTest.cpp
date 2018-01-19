#include "Degenerate_string.hpp"
#include "Parser.hpp"
#include "globalDefs.hpp"
#include "gtest/gtest.h"
#include <iterator>
#include <limits.h>
#include <string>
#include <vector>
#include <iostream>

using namespace deglpf;

TEST(parserTest, VariousDegeneratePatterns) {
  
  std::ifstream infile("test_files/testParser.txt");

  std::vector<SEEDS> s = {{{1, 1, 2}, {3, 3},  {},  {4, 1},},
  {{},{},{},{}},
   {{1, 1, 2, 3}},
   {{},{},{1, 1, 2}},
  {{1, 1, 2}, {}, {}},
   {{1, 1, 2}, {}, {2, 2}},
    {{},{},{1, 1, 2},{},{}}
    };

  std::vector<DEGENERATE_SYMBOLS> ds = {
      {{false, true, false, true, false, false},
       {false,true, true, false, false, false},
       {false,true, true, true, true, false}},
      {{false,true, false, true, false, false},
       {false,true, false, true, false, false},
       {false,true, true, true, false, false}},
      {},
      {{false,true, false, true, false, false}, {false,true, false, true, false, false}},
      {{false,true, false, true, false, false}, {false,true, false, true, false, false}},
      {{false,true, false, true, false, false}, {false,true, false, true, false, false}},
      {{false,true, false, true, false, false},
       {false,true, false, true, false, false},
       {false,true, false, true, false, false},
       {false,true, false, true, false, false}}

  };

  std::vector<std::vector<UINT>> ind = {
    {3,6,7},
    {0,1,2},
    {},
    {0,1},
    {3,4},
    {3,4},
    {0,1,5,6}
  };
std::string alphabet = "ACGTU";
  Parser parser(AlphabetType::DNA, alphabet);
  std::vector<Degenerate_string> dgs(s.size(),
                                     Degenerate_string(alphabet.size()));
  std::string line;
  // Get the first sequence
  std::getline(infile, line);
  int i=0;
  do {
    if (!line.empty()) {
      ;
      parser.parse_sequence(infile, dgs[i++]);
    }
  } while (std::getline(infile, line)); // file ends

  for (int i = 0; i < dgs.size(); ++i) {
    Degenerate_string &dstr = dgs[i];
    const SEEDS &seeds = dstr.get_seeds();
    assert(s[i].size() == seeds.size());
    for (int j = 0; j < seeds.size(); ++j) {
      for (int k = 0; k < seeds[j].size(); ++k) {
        EXPECT_EQ(s[i][j][k], seeds[j][k]);
      }
    }
    const DEGENERATE_SYMBOLS &symb = dstr.get_degenerate_symbols();
    const std::vector<UINT> &indices = dstr.get_degenerate_indices();
    for (int j = 0; j < seeds.size() - 1; ++j) {
      for (int k = 0; k < alphabet.size()+1; ++k) {
        EXPECT_EQ(ds[i][j][k], symb[j][k]);
      }
      EXPECT_EQ(ind[i][j], indices[j]);
    }
  }
}

