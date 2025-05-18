import sys
import numpy as np
from typing import List, Dict


def print_data(data: List[List[int]]) -> None:
    """
    Prints the data collected during the Boyer-Moore algorithm to the
    console using the required format.

    Parameters
    ----------
    data : List[List[int]]
        List containing a List[int] for each pattern, that contains the
        number BCR won over GSR, the total number of occurences of the
        pattern in the sequence and the indices of the first 10 occurences.
    """
    for pat_data in data:
        print(pat_data[0], end="")
        for d in pat_data[1:]:
            print(f" / {d}", end="")
        print()

class BoyerMoore:
    """
    Implementation of the Boyer-Moore algorithm for substring search.

    Parameters
    ----------
    sequence_path : str
        Path to the sequence file in fasta format.
    patterns_path : str
        Path to the patterns file in fasta format.
    
    Methods
    -------
    BCR(t, p, P_idx) -> int
        Calculates the shift according to the Bad-Character-Rule for indices t
        and p for the pattern that corresponds to P_idx.
    GSR(t, p, P_idx) -> int
        Calculates the shift according to the Good-Suffix-Rule for indices t
        and p for the pattern that corresponds to P_idx.
    boyer_moore() -> None
        Implementation of the Boyer-Moore algorithm. Prints the first 10
        occurences of each pattern to the console.
    """
    def __init__(self, sequence_path: str, patterns_path: str):
        self.sequence:       str = self._parse_fasta(sequence_path)
        self.patterns: List[str] = self._parse_fasta(patterns_path, multiple=True)
        self.c_map: Dict[str, int] = self._character_map()
        self.pat_map: Dict[str, int] = self._patterns_map()
        self.preprocessing = self._preprocessing()

    def _patterns_map(self) -> Dict[str, int]:
        """
        Creates the patterns map that maps each pattern to its index in
        self.patterns.

        Returns
        -------
        Dict[str, int]
            The patterns map.
        """
        return { p: idx for idx, p in enumerate(self.patterns)}

    def _character_map(self) -> Dict[str, int]:
        """
        Creates the character map that maps each character to its row in the
        precomputed lookup-tables.

        Returns
        -------
        Dict[str, int]
            The character map.
        """
        unique_seq = set(self.sequence)
        unique_pat = set()
        for pat in self.patterns:
            for c in set(pat):
                unique_pat.add(c)
        
        for c in unique_pat:
            unique_seq.add(c)
        
        return { c: idx for idx, c in enumerate(unique_seq) }
        
    def _preprocessing(self) -> List[np.ndarray]:
        """
        Precomputes the BCR lookup-tables for each pattern in self.patterns.

        Returns
        -------
        np.ndarray
            A 3-dimensional array that contains the 2-dimensional lookup tables
            for each pattern. Following [pattern][character][index].
        """
        prep = []
        for pat in self.patterns:
            arr = np.full((len(self.c_map), len(pat)), fill_value=-1)
            for c in self.c_map.keys():
                for i in range(len(pat)):
                    if pat[i] == c:
                        for j in range(i + 1, len(pat)):
                            arr[self.c_map[c], j] = i
            prep.append(arr)
        return prep
    
    def _parse_fasta(self, path: str, multiple=False):
        """
        Turns the contents of a fasta file to a str or to a List[str] if
        multiple=True.
        
        Parameters
        ----------
        path : str
            Path to the fasta file.
        multiple : bool
            Flag for when a sequence file or a patterns file is read.

        Returns
        -------
        str | List[str]
            Contents of the fasta file as str for multiple=False and list of
            patterns as List[str] for multiple=True.
        """
        with open(path, "r") as file:
            line = file.readline()  # skip first line
            if not multiple:        # for sequence file
                return file.read().replace("\n", "")
            else:
                patterns: List[str] = []
                while line != "":   # for patterns file
                    line = file.readline()
                    if line.startswith(">") or line == "":
                        continue
                    else:
                        patterns.append(line.replace("\n", "", 1))
                return patterns

    def BCR(self, t: int, p: int, P_idx: int) -> int:
        """
        Implementation of the Bad-Character-Rule.
        
        Parameters
        ----------
        t : int
            The index into the sequence T.
        p : int
            The index into the pattern P.
        P_idx : int
            The index of the pattern in the preprocessed lookup table.

        Returns
        -------
        int
            The number of characters, that can be skipped according to this
            rule.
        """
        x = self.sequence[t + p]
        prep: np.ndarray = self.preprocessing[P_idx]
        shift = prep[self.c_map[x]][p]

        if shift == -1:  # Case 1 & 3
            return p + 1
        else:            # Case 2
            return p - shift

    def GSR(self, t: int, p: int, P_idx: int) -> int:
        """
        Implementation of the Good-Suffix-Rule.
        
        Parameters
        ----------
        t : int
            The index into the sequence T.
        p : int
            The index into the pattern P.
        P_idx : int
            The index of the pattern in the preprocessed lookup table.

        Returns
        -------
        int
            The number of characters, that can be skipped according to this
            rule.
        """
        P = self.patterns[P_idx]
        suffix = self.sequence[t:t+p+1]
        suf_len = len(suffix)
        positions: List[int] = []
        for i in range(len(P) - suf_len):
            if P[i:i+suf_len] == suffix:
                if i < p:
                    positions.append(i)
        if not positions:
            return 1
        else:
            return p - positions[-1]

    def boyer_moore(self) -> List[List[int]]:
        """
        Implementation of the Boyer-Moore algorithm for substring search. 
        Prints the indices of the first 10 occurences of each pattern in
        self.patterns to the console.

        Returns
        -------
        List[List[int]]
            A table. Each row corresponds to a pattern. The first column
            contains number of times the Bad-Character-Rule suggested the
            larger shift, the second column contains the number each pattern
            occurs, the next up to ten columns contain the indices of those
            occurences.
        """
        T = self.sequence
        data: List[List[int]] = []
        for P in self.patterns:
            t = 0
            bcr_count = 0
            match_count = 0
            positions: List[int] = []
            while t <= len(T) - len(P):
                p = len(P) - 1
                match = True
                while match and p >= 0:
                    if T[t + p] == P[p]:
                        p -= 1
                    else:
                        match = False
                
                if match:
                    if match_count < 10: positions.append(t)
                    match_count += 1
                    t += 1
                else:
                    s1 = self.BCR(t, p, self.pat_map[P])
                    if p < len(P) / 2:  # No other occurence of suffix if we are len(P)/2 into the pattern.
                        bcr_count += 1
                        t += s1
                    else:
                        s2 = self.GSR(t, p, self.pat_map[P])
                        if (s1 >= s2): bcr_count += 1
                        t += max(s1, s2)
            
            data.append([bcr_count, match_count] + positions)
        return data

if __name__ == "__main__":
    bm = BoyerMoore(sys.argv[1], sys.argv[2])  # (sequence, pattern)
    data = bm.boyer_moore()
    print_data(data)
