import splines as var
from splines.base_cubica import g

def v_barra(x,c):
    
    n = var.N
    
    soma = 0
    for i in range(n):
        soma += c[i]*g(i,x)
    
    return soma
    