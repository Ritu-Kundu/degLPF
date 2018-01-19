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

/** Implemets the utility functions defined in utilDefs.hpp */

#include "../include/utilDefs.hpp"

namespace deglpf {

static struct option long_options[] = {
    {"alphabet", required_argument, NULL, 'a'},
    {"input-file", required_argument, NULL, 'i'},
    {"output-file", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

/** Decode the input flags
 */
ReturnStatus decodeFlags(int argc, char *argv[], struct InputFlags &flags) {
  int args = 0;
  int opt;
  std::string alph;
  /* initialisation */
  while ((opt = getopt_long(argc, argv, "a:i:o:h", long_options, nullptr)) !=
         -1) {
    switch (opt) {
    case 'a':
      alph = std::string(optarg);
      if (alph == "DNA") {
        flags.alphabet_type = AlphabetType::DNA;
      } else if (alph == "PROT") {
        flags.alphabet_type = AlphabetType::PROT;
      } else if (alph == "GEN") {
        flags.alphabet_type = AlphabetType::PROT;
      } else {
        std::cerr << "Invalid command: wrong alphabet type: " << std::endl;
        return (ReturnStatus::ERR_ARGS);
      }
      args++;
      break;

    case 'i':
      flags.input_filename = std::string(optarg);
      args++;
      break;

    case 'o':
      flags.output_filename = std::string(optarg);
      args++;
      break;

    case 'h':
      return (ReturnStatus::HELP);
    }
  }
  if (args < 3) {
    std::cerr << "Invalid command: Too few arguments: " << std::endl;
    return (ReturnStatus::ERR_ARGS);
  } else {
    return (ReturnStatus::SUCCESS);
  }
}

/*
 * Usage of the tool
 */
void usage(void) {
  std::cout << " Usage: eldes <options>\n";
  std::cout << " Standard (Mandatory):\n";
  std::cout << "  -a, --alphabet \t <str> \t \t `DNA' for nucleotide  "
               "sequences or `PROT' for protein  sequences or `GEN' for "
               "sequences containing A-Z. \n";
  std::cout << "  -i, --input-file \t <str> \t \t Input file  name for "
               "sequences (FASTA format currently).\n";
  std::cout << "  -o, --output-file \t <str> \t \t Output filename.\n";
}

} // end namespace