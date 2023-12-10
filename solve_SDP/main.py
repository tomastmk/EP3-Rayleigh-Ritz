from solve_SDP.decomposition import LDL
from solve_SDP.linear_systems import LyB,LtXDY

def system_solve(matriz,B,b):
    
    decomposition = LDL(matriz,b)
    Y = LyB(decomposition[1:],b,B)
    X = LtXDY(decomposition[1:],b,decomposition[0],Y)
    
    
    return X


    