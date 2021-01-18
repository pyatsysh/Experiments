""" Daily Coding 84: Return number of islands of 1's in a matrix

    This uses dynamic programming by storing visited entries, and checking against them.
    The cool part is the recursive call inside the walker function, which explores the island.
    Algorithm:
        1. For {walk over cells}:
                move to New Cell if it is not in Visited table
                record index of New Cell.
           If {New Cell contains 1}: count++, {Step} in
        2. Function {Step}:
            loop over possible moves:
                Permit move if target cell contains 1
                                           not in Visited
                                           is on Board
                If {move Permitted}: update Visited
                                     call {Step}

TODO: visualize exploration process
"""


def count_islands(board):
    dim2 = len(board[0])
    dim1 = len(board)
    count = 0
    visited = [[False for _ in range(dim2)] for _ in range(dim1)]

    def step_in(i, j):
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for di, dj in moves:
            i1 = i+di; j1 = j+dj
            if (0 <= i1 < dim1 and 0 <= j1 < dim2 and
                    not visited[i1][j1] and
                    0 != board[i1][j1]):
                visited[i1][j1] = True
                step_in(i1, j1)

    for i in range(dim1):
        for j in range(dim2):
            if visited[i][j]: continue
            visited[i][j] = True
            if board[i][j]:
                count += 1
                step_in(i, j)

    return count


def main():
    a = [[1, 0, 0, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 1, 1, 0, 0],
         [0, 0, 1, 0, 0],
         [1, 1, 0, 0, 1],
         [1, 1, 0, 0, 1]]  # 4 islands

    print(count_islands(a))


main()