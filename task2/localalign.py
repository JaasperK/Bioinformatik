import sys
import numpy as np

def parse_fasta(path: str) -> str:
    """
    Reads a fasta file with 2 sequences.

    Parameters
    ----------
    path : str
        The path to the fasta file.

    Returns
    -------
    List[str]
        List of length 2 containing both sequences.
    """
    with open(path, "r") as file:
        sequences = ["", ""]
        i = -1
        for line in file:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                i += 1
            else:
                sequences[i] += line
        return sequences


GAP = -1
MISMATCH = -1
MATCH = 1

def t(a: str, b: str):
    """
    Computes score updates depending on if the two input characters are equal.

    Parameters
    ----------
    a : str
        A character.
    b : str
        A character.

    Returns
    -------
    int
        MATCH if a and b are equal and MISMATCH else.
    """
    if a == b:  return MATCH
    else:       return MISMATCH

def local_align(seq1: str, seq2: str):
    """
    Computes the score and traceback matrices to align the input sequences.

    Parameters
    ----------
    seq1, seq2 : str
        Sequences to be aligned.
    
    Returns
    -------
    traceback_mat : np.ndarray
        Matrix used to trace the alignment cores.
    max_score : int
        Highest score in the score matrix.
    cores : list[tuples[int, int]]
        Indices of all alignment cores in the score and traceback matrices.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    m, n = len(seq1), len(seq2)
    score_mat = np.zeros((m + 1, n + 1))
    traceback_mat = np.zeros((m + 1, n + 1))
    max_score = 0
    cores: list[tuple[int, int]] = list()

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_up   = score_mat[i-1][j] + GAP
            score_left = score_mat[i][j-1] + GAP
            score_diag = score_mat[i-1][j-1] + t(seq1[i-1], seq2[j-1])
            score_mat[i][j] = max(0, score_up, score_left, score_diag)

            if score_mat[i][j] == score_up:
                traceback_mat[i][j] = 1  # 1 for trace up

            if score_mat[i][j] == score_left:
                traceback_mat[i][j] = 2  # 2 for trace left

            if score_mat[i][j] == score_diag:
                traceback_mat[i][j] = 3  # 3 for trace diagonal
                
            if score_mat[i][j] == 0:
                traceback_mat[i][j] = 0  # 0 for end of path 

            if score_mat[i][j] >= max_score:
                if score_mat[i][j] == max_score:  # add core to existing list
                    cores.append((i, j))
                else:                             # empty list and add core
                    cores = list()
                    cores.append((i, j))
                max_score = score_mat[i][j]

    if DEBUG: print(score_mat)
    return traceback_mat, max_score, cores

def trace(seq1: str,
          seq2: str,
          traceback_mat: np.ndarray,
          max_i: int,
          max_j: int):
    """
    Traces an alignment core from the traceback matrix.

    Parameters
    ----------
    seq1 : str
        Original sequence one.
    seq2 : str
        Original sequence two.
    traceback_mat : np.ndarray
        Matrix containing the traceback information.
    max_i : int
        First index of an alignment core.
    max_j : int
        Second index of an alignment core.

    Returns
    -------
    align1 : str
        The traceback of the alignment core in seq1.
    align2 : str
        The traceback of the alignment core in seq2.
    i : int
        The offset into seq1 at which the traceback begins.
    j : int
        The offset into seq2 at which the traceback begins.
    """
    i, j = max_i, max_j
    align1, align2 = str(), str()
    
    while traceback_mat[i][j] != 0:
        if traceback_mat[i][j] == 3:    # diag
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif traceback_mat[i][j] == 2:  # left
            align1 += "_"
            align2 += seq2[j-1]
            j -= 1
        elif traceback_mat[i][j] == 1:  # up
            align1 += seq1[i-1]
            align2 += "_"
            i -= 1
    
    # reverse strings since we constructed them backwards
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2, i, j

def print_alignment(seq1_len: int,
                    seq2_len: int,
                    align1: str,
                    align2: str,
                    i: int,
                    j: int):
    n = max(i, j)
    m = max(seq1_len + align1.count("_"), seq2_len + align2.count("_"))

    s1 = "*" * n + align1
    s2 = "*" * n + align2

    # pad strings with "*" until desired length is reached
    s1 = s1 + "*" * (m - len(s1))
    s2 = s2 + "*" * (m - len(s2))

    print(s1)
    print(s2)
    # print()

if __name__ == "__main__":
    DEBUG = False
    if DEBUG:
        import shutil
        width = shutil.get_terminal_size().columns
        np.set_printoptions(linewidth=width)

    s1, s2 = parse_fasta(sys.argv[1])
    traceback_mat, max_score, cores = local_align(s1, s2)
    
    if DEBUG: print(cores)
    print(int(max_score))
    for max_i, max_j in cores:
        align1, align2, i, j = trace(s1, s2, traceback_mat, max_i, max_j)
        print_alignment(len(s1), len(s2), align1, align2, i, j)
