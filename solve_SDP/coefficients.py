from solve_SDP import aij

# Adaptações das fórmulas demonstradas no relatório

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
