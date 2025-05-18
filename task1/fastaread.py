import sys
from typing import List


def get_char_count(lines: List[str]) -> int:
    sum: int = 0
    for l in lines:
        l = l.replace("\n", "")
        sum += len(l)
    return sum

def print_fasta_seq_len(path: str) -> None:
    with open(path, "r", encoding="utf-8") as file:
        line: str = file.readline()  # we can always skip first line due to fasta format
        lengths: List[int] = []
        lines: List[str] = []
        while line != "":
            line = file.readline()
            if line.startswith(">") or line == "":
                lengths.append(get_char_count(lines))
                lines = []
            else:
                lines.append(line)
        
        for l in lengths:
            print(l)

if __name__ == "__main__":
    path = sys.argv[1]
    print_fasta_seq_len(path)
