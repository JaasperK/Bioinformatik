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

C_MAP = {
    "A" : 0,
    "C" : 1,
    "G" : 2,
    "T" : 3,
    "_" : 4
}

SCORE_MATRIX = np.array([
#     A   C   G   T    _
    [ 1, -1, -1 ,-1,  -1],  # A
    [-1,  1, -1, -1,  -1],  # C
    [-1, -1,  1, -1,  -1],  # G
    [-1, -1, -1,  1,  -1],  # T
    [-1, -1, -1, -1, -10]   # _
])

def la_score(a: str, b: str) -> int:
    """
    Calculates the local alignment score of the input sequences using score 1
    for a match and -1 for everything else.

    Parameters
    ----------
    a, b : str
        Sequences of equal length.
    
    Returns
    -------
    int
        The local alignment score of a and b.
    """
    assert len(a) == len(b)
    
    score = 0
    for i in range(len(a)):
        score += SCORE_MATRIX[C_MAP[a[i]]][C_MAP[b[i]]]
    return score

def local_align(seq1: str, seq2: str):
    pass

if __name__ == "__main__":
    s1, s2 = parse_fasta(sys.argv[1])
    local_align(s1, s2)