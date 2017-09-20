seq1 = "ATAGACGACATACAGACAGCATACAGACAGCATACAGA"
seq2 = "TTTAGCATGCGCATATCAGCAATACAGACAGATACG"

rows = len(seq1) + 1
cols = len(seq2) + 1
match = 2
other = -1

score_matrix = [[0 for col in range(cols)]for row in range(rows)]

for i in range(1, rows):
    for j in range(1, cols):
        similarity = match if seq1[i - 1] == seq2[j - 1] else other
        print(similarity)
