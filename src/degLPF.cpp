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

/** Module containing main() method.
 */

#include <cstdlib>

#include "../include/Degenerate_string.hpp"
#include "../include/Parser.hpp"
#include "../include/Search.hpp"
#include "../include/globalDefs.hpp"
#include "../include/utilDefs.hpp"

using namespace deglpf;
ReturnStatus calculate_lpf(const Parser &parser, const UINT alphabet_size,
                           std::ifstream &infile, std::ofstream &outfile);

int main(int argc, char **argv) {

  /* Decode arguments */
  struct InputFlags flags;
  if (decodeFlags(argc, argv, flags) != ReturnStatus::SUCCESS) {
    usage();
    return 1;
  }
  /* Input file */
  std::string filename = flags.input_filename;
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Cannot open input file \n";
    return static_cast<int>(ReturnStatus::ERR_FILE_OPEN);
  }
  /* Output file */
  filename = flags.output_filename;
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Cannot create output file \n";
    return static_cast<int>(ReturnStatus::ERR_FILE_OPEN);
  }
  /* Create Parser */
  std::string alphabet = cGENAlphabet;
  if (flags.alphabet_type == AlphabetType::DNA) {
    alphabet = cDNAAlphabet;
  } else if (flags.alphabet_type == AlphabetType::PROT) {
    alphabet = cPROTAlphabet;
  }
  Parser parser(flags.alphabet_type, alphabet);

  /* Calculate and test result */
  calculate_lpf(parser, alphabet.size(), infile, outfile);
}

ReturnStatus calculate_lpf(const Parser &parser, const UINT alphabet_size,
                           std::ifstream &infile, std::ofstream &outfile) {
  ReturnStatus status;
  std::string line;
  // Get the first sequence
  std::getline(infile, line);
  if (line.empty()) {
    std::cerr << "No Input: Empty File: " << std::endl;
    return ReturnStatus::ERR_INVALID_INPUT;
  }
  do {
    std::string seq_name;
    std::string seq_value;
    if (!line.empty()) {
      if (line[0] != '>') {
        std::cerr
            << "Invalid Input: Not a FASTA format: Expected '>' at line number"
            << std::endl;
        return ReturnStatus::ERR_INVALID_INPUT;
      }
      seq_name = line.substr(1);
      std::cout << "Processing Sequence: " << seq_name << std::endl;
      /* Encode the sequence */
      Degenerate_string dgs(alphabet_size);
      auto status = parser.parse_sequence(infile, dgs);
      if (status != ReturnStatus::SUCCESS) {
        std::cerr << "Invalid Input: Invalid sequence: " << seq_name
                  << std::endl;
        return status;
      }
      /* Calculate the LPF array and LPF-loc arrays for the sequence */
      auto seq_size = dgs.get_size();
      std::vector<UINT> lpf(seq_size, 0);
      Search search(dgs);

      std::clock_t startTime = clock();
      search.calculate_lpf(lpf);
      std::clock_t stopTime = clock();
      double exec_time =
          static_cast<double>(stopTime - startTime) / CLOCKS_PER_SEC;
#ifdef DEBUG
      // PRINTING FOR DEBUGGING
      std::cout << "LPF ARRAY: \n";
      for (UINT c : lpf) {
        std::cout << c << " ";
      }
      std::cout << "\n";
#endif
/**
      /* Test result 
      // The function assumes that the pre-processing has been done.
      if (!search.naive_test(lpf)) {
        std::cerr << "INCORRECT RESULT FOR THE SEQUENCE: " << seq_name
                  << std::endl;
        break;
      }

**/
      /* Print result */
      // First line of a block: > followed by the sequence name
      outfile << ">" << seq_name << std::endl;
      // Next line of the block: Execution time (in sec)
      outfile << exec_time << std::endl;
      // Next line: values of sequence size and number of degenerate symbols
      // (deleimited by a space)
      outfile << seq_size << " " << dgs.get_numberof_seeds() - 1 << std::endl;
      // Next line: lpf array : each cell deleimited by a space
      for (auto l : lpf) {
        outfile << l << " ";
      }
      outfile << std::endl;
      // The block ends with an empty line to delimit it from the following
      // block
      outfile << std::endl;
    }
  } while (std::getline(infile, line)); // sequence ends

  std::cout << "LPF calculated successfully: " << std::endl;
  return ReturnStatus::SUCCESS;
}