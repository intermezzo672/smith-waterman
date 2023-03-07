#!/usr/bin/python
__author__ = "Kelly Chen"
__email__ = "kelly.chen@yale.edu"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm
### EXTRA CREDIT command: pip install git+https://github.com/intermezzo672/smith-waterman.git
### Link to repository: https://github.com/intermezzo672/smith-waterman

import argparse
import csv

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

template_output = ['-----------', '|Sequences|', '-----------', 
                   'sequence1', 'insertseq1', 'sequence2', 'insertseq2', 
                   '--------------', '|Score Matrix|', '--------------']

### Preprocesses the score file into a nested dictionary 
def preprocessScoreFile(scoreFile):

    def removeSpaces(input):
        return input.split()
    
    score_dict = {}
    # reads in scores and removes white space
    with open(scoreFile, 'r') as csvfile:
        raw_txt = csv.reader(csvfile)
        headers = removeSpaces(next(raw_txt)[0])
        score_dict = {}
        for char in headers:
            score_dict[char] = {}
        for line in raw_txt:
            if line:
                scores = removeSpaces(line[0])
                # stores each score in score in corresponding nested dictionary
                for i, score in enumerate(scores[1:]):
                    score_dict[headers[i]][scores[0]] = int(score)
    return score_dict

### Sets up a matrix of 0s corresponding to input sequence lengths
def setUpMatrix(input_sequences):
    y = len(input_sequences[0])
    x = len(input_sequences[1])
    board = [[0 for j in range(y + 2)] for i in range(x + 2)]

    # empties top left corner of matrix
    clear = [(0, 0), (0, 1), (1, 0)]
    for z in clear:
        board[z[0]][z[1]] = ''
    return board

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap=-2, extGap=-1):
    # preprocess the input scores into a dictionary
    score_dict = preprocessScoreFile(scoreFile)

    # read in the input sequences
    with open(inputFile, 'r') as f_in:
        input_seq = f_in.readlines()
        input_seq = [x.strip() for x in input_seq]

    # initialize structures and variables
    # set up matrices with correct dimensions for score matrix and to keep track of states
    score_matrix = setUpMatrix(input_seq)
    state_matrix = setUpMatrix(input_seq)
    highestScore = 0
    highestLocation = (0, 0)

    ### calculation
    # iterates through each cell to assign a score
    for i, x in enumerate(input_seq[0]):
        # assigns sequence to headers of matrices
        score_matrix[0][i+2] = x
        state_matrix[0][i+2] = x

        # initialize default gap for new column
        ver_gap = openGap
        hor_gap = openGap

        for j, y in enumerate(input_seq[1]):
            # gets substitution from the dictionary
            score = score_dict[x][y]
            
            # assigns sequence to headers of matrices
            score_matrix[j+2][0] = y
            state_matrix[j+2][0] = y

            # determines horizontal and vertical gap
            if state_matrix[j+2][i+1] != 0:
                # continues extending horizontal gap if it already exists
                if 'he' in state_matrix[j+2][i+1] or 'hg' in state_matrix[j+2][i+1]:
                    hor_gap = extGap
                else:
                    hor_gap = openGap
            if state_matrix[j+1][i+2] != 0:
                # continues extending vertical gap if it already exists
                if 've' in state_matrix[j+1][i+2] or 'vg' in state_matrix[j+1][i+2]:
                    ver_gap = extGap
                else:
                    ver_gap = openGap

            # calculates diagonal match/mismatch, horizontal and vertical gap scores 
            diag = score_matrix[j+1][i+1] + score
            ver = score_matrix[j+1][i+2] + ver_gap
            hor = score_matrix[j+2][i+1] + hor_gap
            update = [hor, ver, diag]

            # finds the max value for the cell and where it came from
            max_num = max(update)
            indices = []
            for index, num in enumerate(update):
                if num == max_num:
                    indices.append(index)

            # stores trace of cell in the state matrix, including ties to properly record gaps
            if 0 in indices and hor_gap == extGap:
                state_matrix[j+2][i+2] = str(state_matrix[j+2][i+2]) + 'he'
            if 0 in indices and hor_gap == openGap:
                state_matrix[j+2][i+2] = str(state_matrix[j+2][i+2]) + 'hg'
            if 1 in indices and ver_gap == extGap:
                state_matrix[j+2][i+2] = str(state_matrix[j+2][i+2]) + 've'
            if 1 in indices and ver_gap == openGap:
                state_matrix[j+2][i+2] = str(state_matrix[j+2][i+2]) + 'vg'
            if 2 in indices:
                state_matrix[j+2][i+2] = str(state_matrix[j+2][i+2]) + 'm'

            # gets max positive value (0 if none) and updates cell
            update.append(0)
            update_value = max(update)
            score_matrix[j+2][i+2] = update_value
    
            # updates position of cell with highest score
            if update_value > highestScore:
                highestScore = update_value
                highestLocation = (j+2, i+2)
    
    ### traceback and alignment    
    # initialize output sequences
    seq1 = ''
    seq2 = ''
    pipes = ''
    
    # finds trailing sequence after alignment portion for each sequence and adds on to corresponding output sequence
    trailing_seq1 = []
    trailing_seq2 = []
    if highestLocation[1] != len(score_matrix[0]) - 1:
        trailing_seq1 = score_matrix[0][highestLocation[1]+1:]
    for row in range(highestLocation[0] + 1, len(score_matrix)):
        trailing_seq2.append(score_matrix[row][0])
    for char in trailing_seq2:
        seq2 = seq2 + char
    for char in trailing_seq1:
        seq1 = seq1 + char
    # adds on matching whitespace for trailing alignment
    if len(trailing_seq1) > len(trailing_seq2):
        seq2 = seq2 + (len(trailing_seq1) - len(trailing_seq2)) * ' '
        pipes = ' ' * len(trailing_seq1)
    else:
        seq1 = seq1 + (len(trailing_seq2) - len(trailing_seq1)) * ' '
        pipes = ' ' * len(trailing_seq2)
    
    seq1 = ')' + seq1
    seq2 = ')' + seq2
    pipes = ' ' + pipes

    # starts traceback from cell with highest score
    traceback = highestLocation
    # traces back until cell with score of 0 is reached
    while score_matrix[traceback[0]][traceback[1]] != 0:
        # recalls root states/cells from which the current cell stems from
        state = state_matrix[traceback[0]][traceback[1]]
        if 'm' in state:
            # adds matching/aligned cells to sequences
            char_seq1 = score_matrix[0][traceback[1]]
            char_seq2 = score_matrix[traceback[0]][0]
            seq1 = char_seq1 + seq1
            seq2 = char_seq2 + seq2
            if char_seq1 == char_seq2:
                pipes = '|' + pipes
            else:
                pipes = ' ' + pipes
            traceback = (traceback[0] - 1, traceback[1] - 1)
        elif ('ve' in state and score_matrix[traceback[0]][traceback[1]] + 1 == score_matrix[traceback[0] - 1][traceback[1]]) or ('vg' in state and score_matrix[traceback[0]][traceback[1]] + 2 == score_matrix[traceback[0] - 1][traceback[1]]):
            # adds gap to sequence 1 due to vertical gaps
            char_seq2 = score_matrix[traceback[0]][0]
            seq1 = '-' + seq1
            seq2 = char_seq2 + seq2
            pipes = ' ' + pipes
            traceback = (traceback[0] - 1, traceback[1])
        elif ('he' in state and score_matrix[traceback[0]][traceback[1]] + 1 == score_matrix[traceback[0]][traceback[1] - 1]) or ('hg' in state and score_matrix[traceback[0]][traceback[1]] + 2 == score_matrix[traceback[0]][traceback[1] - 1]):
            # adds gap to sequence 2 due to horizontal gaps
            char_seq1 = score_matrix[0][traceback[1]]
            seq1 = char_seq1 + seq1
            seq2 = '-' + seq2
            pipes = ' ' + pipes
            traceback = (traceback[0], traceback[1] - 1)

    seq1 = '(' + seq1
    seq2 = '(' + seq2
    pipes = ' ' + pipes
    
    # finds leading sequence before alignment portion for each sequence
    lead_str1 = ''
    lead_str2 = ''
    lead_seq1 = score_matrix[0][2:traceback[1]+1]
    lead_seq2 = []
    for char in range(2, traceback[0]+1):
        lead_seq2.append(score_matrix[char][0])
    
    #adds on whitespace for proper alignment depending on length of leading sequence
    if len(lead_seq1) > len(lead_seq2):
        lead_str2 = (len(lead_seq1) - len(lead_seq2)) * ' '
        pipes = len(lead_seq1) * ' ' + pipes
    else:
        pipes = len(lead_seq2) * ' ' + pipes
        lead_str1 = (len(lead_seq2) - len(lead_seq1)) * ' '

    # adds on leading sequence to corresponding output sequence
    for char in lead_seq1:
        lead_str1 = lead_str1 + char
    seq1 = lead_str1 + seq1
    for char in lead_seq2:
        lead_str2 = lead_str2 + char
    seq2 = lead_str2 + seq2

    # array for easy writing of alignment output
    align_output = [f'Alignment Score:{str(highestScore)}', 'Alignment Results:', seq1, pipes, seq2]

    ### write output
    with open('output.txt', 'w') as f_out:
        # outputs initial sequences
        for item in template_output:
            if item == 'insertseq1':
                item = input_seq[0]
            if item == 'insertseq2':
                item = input_seq[1]
            f_out.write(item + '\n')

        # outputs score matrix
        for row in score_matrix:
            for z in row:
                f_out.write(str(z) + '\t')
            f_out.write('\n')
        f_out.write('----------------------\n|Best Local Alignment|\n----------------------\n')

        # outputs alignment results
        for row in align_output:
            f_out.write(row + '\n')

    return