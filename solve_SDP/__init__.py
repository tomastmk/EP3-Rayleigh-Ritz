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
