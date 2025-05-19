import numpy as np

if __name__ == "__main__":
    freq_matrix = np.array([     # Frequency matrix from https://jaspar2020.genereg.net/matrix/MA0036.1/.
        [ 13,  0, 52,  0, 25 ],  # A
        [ 13,  5,  0,  0,  7 ],  # C
        [ 18, 48,  1,  0, 15 ],  # G
        [  9,  0,  0, 53,  6 ]   # T
    ])

    print(freq_matrix / 53)
