#include <limits.h>
#include <vector>
#include <string>
#include <iterator>
#include "Parser.hpp"
#include "globalDefs.hpp"
#include "Degenerate_string.hpp"
#include "Search.hpp"
#include "gtest/gtest.h"


using namespace deglpf;

// TODO: Add test invoving testing intermediate steps

const std::vector<int> num_occ = {5,8,5,7,15};
const std::vector<std::vector<UINT>> lpf = {
  {0, 0, 4, 3, 2, 1, 0, 2, 3, 2, 1},
  {0, 0, 4, 3, 2, 1},
  {0, 1, 1, 3, 2, 1},
  {0, 0, 4, 3, 2, 1},
  {0, 1, 1, 4, 3, 2, 1},
  {0,1,1,5,4,3,2,3,2,1},
  {0, 0, 5, 6, 5, 4, 3, 2, 2, 1 }
};
TEST(alsoTest, MultipleSimpleSeq) {
  std::vector<std::vector<UINT>> result;
  std::string alphabet = "ACGTU";
  Parser parser(AlphabetType::DNA, alphabet);
  std::ifstream infile("test_files/testAlgo.txt");
  std::string line;
  // Get the first sequence
  std::getline(infile, line);
  do {
    if (!line.empty()) {
      Degenerate_string dgs(alphabet.size());
      parser.parse_sequence(infile, dgs);
      UINT seq_size = dgs.get_size();
      std::vector<UINT> lpf(seq_size, 0);
      Search search(dgs);
      search.calculate_lpf(lpf);
      result.push_back(lpf);
    }
  } while (std::getline(infile, line)); // file ends


  for (int i=0; i < lpf.size(); ++i) {
    EXPECT_EQ(lpf[i].size(), result[i].size());
    for (auto j = 0; j < lpf[i].size(); ++j) {
      EXPECT_EQ(lpf[i][j], result[i][j]);
    }
   } 
}


