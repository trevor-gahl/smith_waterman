seq1 = "ATAGACGACATACAGACAGCATACAGACAGCATACAGA"
seq2 = "TTTAGCATGCGCATATCAGCAATACAGACAGATACG"

rows = len(seq1) + 1
cols = len(seq2) + 1
match = 2
other = -1
maxScore = 0
maxPosition = (0, 0)

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
