import numpy as np

def print_result(seq1, seq2, score_matrix):
    print("score_matrix: ")
    print(score_matrix)
    print()
    print("aligned seqeunce 1: ", seq1)
    print("aligned seqeunce 2: ", seq2)

def matching(x, y, match_award, mismatch_penalty):
    if x == y:
        return match_award
    else:
        return mismatch_penalty

def water(seq1, seq2, match_award=3, mismatch_penalty=-3, gap_penalty=0):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    score_matrix = np.zeros((len_seq1+1, len_seq2+1))
    pointer_matrix = np.zeros((len_seq1+1, len_seq2+1))
    
    max_score = 0
    
    for i in range(1, len_seq1+1):
        for j in range(1, len_seq2+1):
            score_diag = score_matrix[i-1][j-1] + matching(seq1[i-1], seq2[j-1], match_award, mismatch_penalty)
            score_up   = score_matrix[i][j-1] + gap_penalty
            score_left = score_matrix[i-1][j] + gap_penalty
            score_matrix[i][j] = max(0, score_diag, score_up, score_left)
            if score_matrix[i][j] == 0:
                pointer_matrix[i][j] = 0
            if score_matrix[i][j] == score_diag:
                pointer_matrix[i][j] = 1
            if score_matrix[i][j] == score_up:
                pointer_matrix[i][j] = 2
            if score_matrix[i][j] == score_left:
                pointer_matrix[i][j] = 3
            if score_matrix[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score_matrix[i][j]
    align1 = []
    align2 = []
    i = max_i
    j = max_j
    del max_i, max_j
    
    while pointer_matrix[i][j] != 0:
        if pointer_matrix[i][j] == 1:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif pointer_matrix[i][j] == 2:
            align1.append(None)
            align2.append(seq2[j-1])
            j -= 1
        elif pointer_matrix[i][j] == 3:
            align1.append(seq1[i-1])
            align2.append(None)
            i -= 1
            
    align1.reverse()
    align2.reverse()
    return align1, align2, score_matrix

if __name__ == '__main__':
    l1 = "TGTTACGG"
    l2 = "GGTTGACTA"
    match_award=3
    mismatch_penalty=-3
    gap_penalty=0

    res1, res2, score_matrix = water(l1,l2, match_award, mismatch_penalty, gap_penalty)
    print_result(res1, res2, score_matrix)