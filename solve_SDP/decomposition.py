from solve_SDP.coefficients import lji,di

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
    
