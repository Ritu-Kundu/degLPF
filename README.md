# degLPF: Computing Longest Previous Factor (LPF) Array in a Degenerate String


*degLPF* is a tool that computes the LPF (Longest Previous Factor) Array for a given degenerate sequence.
Currently it takes only the files mimicking FASTA format.

>A *degenerate string* is defined by the existence of one or more positions that are represented by sets of symbols from
>an alphabet Σ, unlike an accurate or certain (standard) string characterised by a single symbol at each position.
>For instance, {a,b} a c {b,c} a {b,c} is a degenerate string of length 6 over Σ = {a,b,c}.

>Thelongest previous factor array (LPF) of a degenerate string *T* of length *n* is an array storing the length of the longest factor at each position *i* in *T* which matches a factor with an occurrence >starting from the left of *i*. Formally stated:
>LPF[i]= max{l:T[i..i+l−1]≈T[j ..j+l−1],0≤ j < i}

The tool is based on the terminology and the algorithm described in the paper titled
[**Longest Previous Factor Array for Degenerate Strings**](https://link_to_the_paper)
by *Costas S. Iliopoulos, Ritu Kundu, Fatima Vayani, and Steven Watts*.

To compile degLPF, please follow the instructions given in file `INSTALL.md`


## Usage of the tool: 
```
Usage: degLPF <options>
 Standard (Mandatory):
  -a, --alphabet 	 	<str> 	 	 'DNA' for nucleotide  sequences or 'PROT' for protein  sequences or 'GEN' for general (A-Z)  sequences. 
  -i, --input-file 	 	<str> 	 	 Input file  name (Mimicing FASTA format currently).
  -o, --output-file		<str> 	 	 Output filename.
```

 **Example:** 

You can try the tool for compression with the given sample files (`sample` folder) via
```sh
./bin/degLPF -a DNA -i sample/input.txt -o sample/output.txt
```
Here, the degenerate sequences are given in the file "input.txt" which is in subfolder "sample" of current folder. 
LPF Array for each sequence will be calculated and result will be written in the file "output.txt" in subfolder "sample" of current folder.


## Notes
- Alphabet can be single letter codes representing bases in genomic sequences (for DNA) or amino acids in protiens (for PROT) or general A-Z (for GEN).
  * Valid DNA letters: ACGT (irrespective of case)
  * Valid Prot letters: ACDEFGHIKLMNPQRSTUVWY (irrespective of case)


- The current implementation does not have a provision to accommodate more than 255 symbols (maximum value of an unsigned character) due to the limitation of the external library used (SDSL). Therefore, the number of degnerate symbols in a sequence should not be more than (255-alphabet size).

- Input file is expected to be in format resembling a valid [FASTA format] (https://en.wikipedia.org/wiki/FASTA_format).
 * From the first line, a block representing a sequence starts. It ends with either an empty line or end of the file.
   - It starts with a '>' followed by an identifier for the text.
   - From the next line (until an empty line or end of file is hit), the sequence of letters from the specified alphabet and degenerate symbols starts.
   - A Degenerate symbol is represented as follows:
     -- It starts with symbol '{' (open curly braces).
     -- It ends with symbol '}' (close curly braces).
     -- Letters within a symbol may be separated by spaces or may not be delimited at all.
     -- Same letter may be repeated but will be seen as one.
    -- At least two letters should be present.
   - White-space characters and new lines within a sequence block are allowed (as these are being ignored).


- Broadly, working of the tool is as follows:
  * The tool parses and encodes the sequence (read from the input file) into numeric alphabet (1 to alphabet-size). 
  * It then preprocesses the sequence and computes the LPF table.
  * A function to test the resulting array (using the naive approach) has also been provided. However, currently the function is not being called (its call has been commented out).
  * The result is written in the output file.

- Output file is in the following format:
 * Corresponding to each sequence, there is a block (two blocks are separated by an empty line): 
  * The first line in the block begins with a '>' followed by the identifier (FASTA format) of the sequence.
  * The next line gives the time (in seconds) used for calculation (after input file has been read in memory up to calculating the array).
  * The next line contains the following pieces of information (separated by a blank space):
    - The length of the sequence
    - The number of the degenerate symbols in the sequence
  * The following line contains the LPF Array (each element delimited by a blank space).

## Running Experiments
To run the experiments, use the following command:
```sh
python3 scripts/experiments.py 
```
It generates the specified number of files; each containing the specified number of random sequences from DNA alphabet of specified length and the specified number of degenerate symbols. 
* Each file contains the (randomly generated) sequences of a given length (for each length specified in the array `text_size`); The number of sequences decided by `num_seq_in_file`  .
* Each sequence of the specified length will have copies equal to size of the array `k`.
* Each copy will have randomly distributed degenerate symbols equal to the chosen value from the array `k`.
* Each symbol will have a randomly chosen collection from the alphabet.


Setting parameters of the experiments:
=======================================
Change the following settings in the file `experiments.py` (in `scripts` folder) -
- `alphabet`: Alphabet used. Currently -> ['A','C', 'G', 'T']
- `text_size`: Array containing the lengths of the sequences. One file will begenerated for each length.
- `k`: Array containing the number of degenerate symbols. For each sequence, each value in this array is being tested. The positions of the symbols are randomly generted and the set of letters at that position is randomly generated as well. 
- `num_seq_in_file`: Number of sequences in a file. Currently -> 1
- `param_separator`: Separator of various stats parameters collected for a file. Currently -> '\t'

Stats collected from the experiments:
======================================
The files are generated in the `data` subfolder of the `experiments` folder. The `stats.txt` file is created in the `experiments` folder.

**Names of the files:**
- `input$i.txt`: Name of the sequence file (Mimics Fasta format). Here, `$i` is the ith file generated.
- `output$i.txt`: Name of the corresponding output file.

**Stats Info**

The following tab-separated parametres are recorded in the `stats.txt` for each file (One record corresponding to each sequence of each file, delimited by newline).
- `n`: The length of the sequence.
- `k`: The number of degenerate symbols in that sequence
- `time`: Time taken for calculation of the LPF Array (in sec).
 


## External Libraries

 * RMQ on LCP array is used to answer longest common prefix queries. For answering these queries, following libraries have been used:
   + [sdsl](https://github.com/simongog/sdsl-lite)
 * For testing [googletest](https://github.com/google/googletest) framework has been used.

