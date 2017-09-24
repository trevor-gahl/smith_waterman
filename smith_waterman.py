seq1 = "ATAGACGACATACAGACAGCATACAGACAGCATACAGA"
seq2 = "TTTAGCATGCGCATATCAGCAATACAGACAGATACG"

rows = len(seq1) + 1
cols = len(seq2) + 1
match = 2
other = -1
maxScore = 0
maxPosition = (0, 0)


def createScoreMatrix(rows,cols):
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
    print(score_matrix)
    print(maxScore)
    print(maxPosition)
    print(rows, cols)


def traceback(score_matrix, start_pos):

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos
    move         = next_move(score_matrix, x, y)
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

        move = next_move(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq1[y - 1])

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def nextMove(score_matrix, x, y):

    diag = score_matrix[x - 1][y - 1]
    up   = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')
