import string
import random
import re
import math
import os
import subprocess
import time
import sys
######### Settings of the experiments #############
alphabet = ['A', 'C', 'G', 'T']
# Each file contains the (randomly generated) sequences of a given length (from  'text_size'); The number of sequences decided by 'num_seq_in_file'  .
# Each sequence of the specified length will have copies equal to size of 'k' (Here 5).
# Each copy will have randomly distributed degenerate symbols equal to the given value from k.
# Each sybol will have a randomly chosen collection from the alphabet
text_size = [1000, 2000, 4000, 8000, 16000, 32000, 64000]
k = [5, 10, 20, 40, 80]
num_seq_in_file = 1
param_separator = '\t'
###################################################

FOLDER = './experiments/'
DATA_FOLDER = 'data/'
INPUT_FILE_NAME = 'input'
OUTPUT_FILE_NAME = 'output'
STATS_FILE_NAME = 'stats.txt'

stats_param = ['n', 'k', 'time']


def memory_usage_resource():
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # OSX produces the output in different units
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def collect_stats(o_file, stats_file):
    n = []
    d = []
    tm = []
    block_size = 5
    with open(o_file, "r") as f:
        mode = 0
        # 1st line seq name, next time, next n & k , next LPF array, next blank
        for line in f:
            if (mode == 0 or mode == 3 or mode == 4):
                # Ignore line                  
                dummy = line
            elif (mode ==1):
                tmf = re.findall("\d+\.\d+", line)
                tm.append(str(tmf[0]))
            else:
                seq_len, sd = line.split(' ')
                n.append(seq_len)
                d.append(sd)
            mode = (mode + 1) % block_size
    num = num_seq_in_file * len(k)
    for ind in range(num):
        stats_file.write(str(n[ind]) + param_separator + str(k[ind]) +
                         param_separator + tm[ind] + '\n')


def write_file(seq_file, seq_size, seq, d):
    alph_size = len(alphabet)
    symb_pos = random.sample(range(seq_size), d)
    symb_pos_sorted = sorted(symb_pos)
    symb_pos_sorted.append(seq_size)
    nxt = 0
    txt = ''
    for ind in range(seq_size):
        if ind != symb_pos_sorted[nxt]:
            txt += seq[ind]
        else:
            sym_size = random.randint(2, alph_size)
            symb_letters = random.sample(alphabet, sym_size)
            symb = ''.join(symb_letters)
            txt = txt + '{' + symb + '}'
            nxt += 1
    seq_file.write('>seq ' + str(seq_size) + '_' + str(d) + '\n')
    seq_file.write(txt + '\n\n')


def main():
    i_filename = FOLDER + DATA_FOLDER + INPUT_FILE_NAME
    o_filename = FOLDER + DATA_FOLDER + OUTPUT_FILE_NAME
    sf = open(FOLDER + STATS_FILE_NAME, 'w')
    sf.write(param_separator.join(stats_param))
    sf.write('\n')

    num_files = len(text_size)

    # Generate Data Files
    for i in range(num_files):
        seq_size = text_size[i]
        suff = str(i) + '.txt'
        seq_filename = i_filename + suff
        seq_file = open(seq_filename, 'w')
        for j in range(num_seq_in_file):
            seq = ''.join(random.choice(alphabet) for i in range(seq_size))
            for d in k:
                write_file(seq_file, seq_size, seq, d)
        seq_file.close()
    
    print("$$$$$$$$$$$$$$$$$$$$ FILE GENERATION COMPLETE $$$$$$$$$$$$$$$$$$$$$$")

    # Run the tool on the files
    for i in range(num_files):
        suff = str(i) + '.txt'

        # Call the tool
        tool = './bin/degLPF'
        cmd = tool + ' -a DNA -i ' + i_filename + suff + ' -o ' + o_filename + suff
        print('COMMAND: ' + cmd)
        comp = subprocess.Popen(cmd, shell=True)
        comp.wait()

    print("$$$$$$$$$$$$$$$$$$$$ FILE PROCESSING COMPLETE $$$$$$$$$$$$$$$$$$$$$$")

    # Analyse the output files
    for i in range(num_files):
        suff = str(i) + '.txt'
        collect_stats(o_filename + suff, sf)

    sf.close()


main()
