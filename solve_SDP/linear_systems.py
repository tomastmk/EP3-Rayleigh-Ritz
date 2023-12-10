from solve_SDP import aij

class bcolors:
    BOLD = "\033[1m"
    ENDC = "\033[0m"

# Resolve o sistema Ly = B
def LyB(Lower,b,B):
    n = len(Lower[0])
    Y = []
    
    for i in range(n):
        soma = 0
        
        for k in range(0,i,1):
            soma+=aij(Lower,i,k,b)*Y[k]
        
        yi = B[i]-soma
        
        Y.append(yi)
        
    return Y    

# Resolve o sistema LtX = D^-1Y
def LtXDY(Lower,b,Diagonal,Y):

    n = len(Lower[0])
    X = ["a" for i in range(n)]
    for i in range(n-1,-1,-1):
        
        
        soma = 0
        for k in range(i+1,n,1):
            soma+=aij(Lower,k,i,b)*X[k]
        
        
        if Diagonal[i]==0:
            print()
            print(bcolors.BOLD+"O sistema possui uma linha nula"+bcolors.ENDC)
            print()
            return 0
        
        xi = Y[i]/Diagonal[i]-soma
        
        X[i] = xi
        
    return X

        