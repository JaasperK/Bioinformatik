import sys
import numpy as np
from typing import List

def parse_fasta(path: str) -> List[str]:
    """
    Parses the contents of a fasta file to a list of sequences.
    
    Parameters
    ----------
    path : str
        Path to the fasta file.

    Returns
    -------
    List[str]
        Contents of the fasta file.
    """
    with open(path, "r") as file:
        line = file.readline()  # skip first line
        patterns: List[str] = []
        while line != "":   # for patterns file
            line = file.readline()
            if line.startswith(">") or line == "":
                continue
            else:
                patterns.append(line.replace("\n", "", 1))
        return patterns

def to_matrix(content: List[str], context_size=15) -> List[str]:
    """
    Truncate each string in thhe input list to be of length context_size.

    Parameters
    ----------
    content : List[str]
        List of sequences.
    context_size : int (default=15)
        The target size that each string will be reduced to.

    Returns
    -------
    List[str]
        List of str where each str has length context_size.
    """
    mat = []
    for s in content:
        mat.append(s[:context_size])
    return mat

def pswm(mat: List[str]) -> np.ndarray:
    """
    Computes the Position-Specific-Weight-Matrix (PSWM) for the input mat.

    Parameters
    ----------
    mat : List[str]
        List of strings of equal length.

    Returns
    -------
    np.ndarray
        The PSWM of the input list.
    """
    # TODO
    return 1

if __name__ == "__main__":
    content = parse_fasta(sys.argv[1])
    mat = to_matrix(content)
    print(mat)
    print(pswm(mat))
