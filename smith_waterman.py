#########################################################################
## Name: Trevor Gahl, David Schwer, Procassius Loftin                  ##
## Course: Computational Biology                                       ##
## Date: 9/25/17                                                       ##
## About: Python implementation of smith-waterman algorithm            ##
##        based on work from https://gist.github.com/radaniba/11019717 ##
#########################################################################

import sys

seq1 = "ATAGACGACATGGGGACAGCAT"
seq2 = "TTTAGCATGCGCATATCAGCAATA"

#seq1 = 'AGCACACA'
#seq2 = 'ACACACTA'

match = 2
other = -1
maxScore = 0
maxPosition = (0, 0)

########################################
## Main method to run local alignment ##
########################################

def main():
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix, start_pos = createScoreMatrix(rows, cols)
    seq1_aligned, seq2_aligned = traceback(score_matrix, start_pos)
    assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'

# Pretty print the results. The printing follows the format of BLAST results
# as closely as possible.
    alignment_str, idents, gaps, mismatches = alignment_string(seq1_aligned, seq2_aligned)
    alength = len(seq1_aligned)
    print('\n')
    print('Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
        alength, idents / alength, gaps, alength, gaps / alength))
    print('\n')
    for i in range(0, alength, 60):
        seq1_slice = seq1_aligned[i:i+60]
        print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
        print('             {0}'.format(alignment_str[i:i+60]))
        seq2_slice = seq2_aligned[i:i+60]
        print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
        print('\n')

###############################################
## Creates scoring matrix from input strings ##
###############################################

def createScoreMatrix(rows,cols):
    global maxScore
    score_matrix = [[0 for col in range(cols)]for row in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            similarity = match if seq1[i - 1] == seq2[j - 1] else other
            diag_score = score_matrix[i - 1][j - 1] + similarity
            up_score = score_matrix[i - 1][j] + other
            left_score = score_matrix[i][j - 1] + other
            curMax = max(0, diag_score, up_score, left_score)
            if curMax > maxScore:
                maxScore = curMax
                maxPosition = (i, j)
            score_matrix[i][j] = curMax
            # print(similarity)
    '''
    print(score_matrix)
    print(maxScore)
    print(maxPosition)
    print(rows, cols)
    '''
    print('\nSCORE MATRIX: \n')
    print_matrix(score_matrix)
    return score_matrix, maxPosition

#######################################
## Creates best fit alignment string ##
## based on scoring matrix           ##
#######################################

def traceback(score_matrix, start_pos):

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos
    move         = nextMove(score_matrix, x, y)
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1

        move = nextMove(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq1[y - 1])

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))

################################################
## Determines the next move for the traceback ##
################################################

def nextMove(score_matrix, x, y):

    #Assign the diagonal score
    diag = score_matrix[x - 1][y - 1]

    #Assign insertion/deletion scores
    up   = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]

    #Check all three cases to find next character/insertion/deletion
    if diag >= up and diag >= left:
        return 1 if diag != 0 else 0
    elif up > diag and up >= left:
        return 2 if up != 0 else 0
    elif left > diag and left > up:
        return 3 if left != 0 else 0

    #Error detection
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')

###############################################
## Creates the alignment string for printing ##
###############################################

def alignment_string(aligned_seq1, aligned_seq2):

    #Sets initial values
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []

    #Runs through both strings
    for base1, base2 in zip(aligned_seq1, aligned_seq2):

        #Checks for match
        if base1 == base2:
            alignment_string.append('|')
            idents += 1

        #Checks for insertion/deletion
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1

        #If neither of the above, it's mismatch
        else:
            alignment_string.append(':')
            mismatches += 1

    #Returns the "alignment" string and the alignment characteristics
    return ''.join(alignment_string), idents, gaps, mismatches

def print_matrix(matrix):
    for row in matrix:
        for col in row:
            print '{:4}'.format(col),
        print

if __name__ == '__main__':
    sys.exit(main())
