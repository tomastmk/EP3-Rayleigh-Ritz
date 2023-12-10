from solve_SDP.main import system_solve

A,bs = [
    [8, 20, 11, 17, 8],
    [7,  6,  3,  7, 0],
    [1,  7,  1,  0, 0]],2

B=[-13,-40,-4,-15,-38]

print("\nX =",system_solve(A,B,bs),"\n")