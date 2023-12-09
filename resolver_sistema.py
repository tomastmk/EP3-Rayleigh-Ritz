import random

# O elemento aij tem índices ij na matriz "original", mas índices i_mod e i na matriz "modificada". Essa função faz a conversão de índices e retorn o elemento aij.
def aij(matriz,i,j,b): 
    
    # Limite da quantidade de linhas da matriz modificada
    if j-i>=len(matriz):
        return 0
    
    # Simetria
    if i>j:
        i,j=j,i
        
    # Estrutura de banda
    if j>b+i:
        return 0
    
    i_mod = j-i
    
    return matriz[i_mod][i]

# Adaptações das fórmulas demonstradas no relatório
def lji(matriz,b,Lower,Diagonal,V,j,i):
    
    # Evita cálculos repetidos
    if aij(Lower,i,j,b)!="a":
        return aij(Lower,i,j,b)
    
    # Somatório
    soma = 0
    for k in range(0,i,1):
        soma += lji(matriz,b,Lower,Diagonal,V,j,k)*vj(matriz,b,Lower,Diagonal,V,i,k)
    
    d = di(matriz,b,Lower,Diagonal,V,i)
    
    Diagonal[i]=d # Salva os valores de d
    
    return (aij(matriz,j,i,b)-soma)/d

def vj(matriz,b,Lower,Diagonal,V,i,j):

    if V[i][j] != "a":
        return V[i][j]
    
    else:
        vj = lji(matriz,b,Lower,Diagonal,V,i,j)*di(matriz,b,Lower,Diagonal,V,j)
    	
    V[i][j]=vj
     
    return vj    

def di(matriz,b,Lower,Diagonal,V,i):
    
    # Evita cálculos repetidos
    if Diagonal[i]!="a":
        return Diagonal[i]
    
    # Somatório
    soma = 0
    for j in range(0,i,1):
        soma += lji(matriz,b,Lower,Diagonal,V,i,j)*vj(matriz,b,Lower,Diagonal,V,i,j)
    
    return aij(matriz,i,i,b)-soma


# A função abaixo retorna as matrizes Lower e Diagonal 
def LDL(matriz,b):
    
    n = len(matriz[0])

    decomposition = [["a" for i in range(n)] for i in range(len(matriz)+1)] # o "a" indica que aquele elemento ainda não foi calculado
    decomposition[1] = [1 for i in range(n)] # Diagonal principal
    
    V = [["a" for k in range(n)]for k in range(n)]
    
    for i in range(1,n): 
        for j in range(i):
            
             # Estrutura de banda
            if i>b+j:
                continue            
            
            j_mod= i-j+1
            decomposition[j_mod][j]= lji(matriz,b,decomposition[1:],decomposition[0],V,i,j)

    decomposition[0][-1]=di(matriz,b,decomposition[1:],decomposition[0],V,n-1) # O último elemento da diagonal precisa ser calculado separadamente pois o índice j do for acima vai até n-2 

    return decomposition
    

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


def solve_for_X(A,B,bs):
    

    decomposition = LDL(A,bs)
    n = len(decomposition[0])

    Y = LyB(decomposition[1:],bs,B)
    X = LtXDY(decomposition[1:],bs,decomposition[0],Y)
     
    if X == 0:
        return
    
    return X
