import splines as var
from splines.gauss_quadrature import gauss_3_pontos, m_ap, m_aq, m_d

def aij(i,j):
   
    if i-2 < 0:
        P1 = 0
        Q1 = 0
    else:
        P1 = gauss_3_pontos(m_ap,i-2,i,j)
        Q1 = gauss_3_pontos(m_aq,i-2,i,j)

    
    if i-1 < 0:
        P2 = 0
        Q2 = 0
    else:
        P2 = gauss_3_pontos(m_ap,i-1,i,j)
        Q2 = gauss_3_pontos(m_aq,i-1,i,j)
    
    Q3 = gauss_3_pontos(m_aq,i,i,j)
    P3 = gauss_3_pontos(m_ap,i,i,j)
    
    if i >= var.N:
        P4 = 0
        Q4 = 0
    
    else: 
        P4 = gauss_3_pontos(m_ap,i+1,i,j)
        Q4 = gauss_3_pontos(m_aq,i+1,i,j)
    
    return var.H/2 * (P1+P2+P3+P4 + Q1+Q2+Q3+Q4)


def matrix_A():

    n = var.N
    
    A = [ [],[],[],[] ]
    
    for i in range(n):
        A[0].append(aij(i,i))
    for i in range(n-1):
        A[1].append(aij(i,i+1))
    for i in range(n-2):
        A[2].append(aij(i,i+2))
    for i in range(n-3):
        A[3].append(aij(i,i+3))
    
    return A


def di(i):
    
    h = var.H
    n = var.N
    
    if i == 0:
        
        return h/2 * (gauss_3_pontos(m_d, 0, 0,0) + gauss_3_pontos(m_d, 1 , 0,0)) - 2*h*gauss_3_pontos(m_d,0,-1,0)
    
    elif i == 1:
        
        return h/2 * (gauss_3_pontos(m_d,0,1,0)+gauss_3_pontos(m_d,1,1,0)+gauss_3_pontos(m_d,2,1,0)-gauss_3_pontos(m_d,0,-1,0))
    
    elif 2 <= i <= n-1:
        
        return h/2 * (gauss_3_pontos(m_d,i-2,i,0)+gauss_3_pontos(m_d,i-1,i,0)+gauss_3_pontos(m_d,i,i,0)+gauss_3_pontos(m_d,i+1,i,0))
    
    elif i == n:
        
        return h/2 * (gauss_3_pontos(m_d,n-2,n,0)+gauss_3_pontos(m_d,n-1,n,0)+gauss_3_pontos(m_d,n,n,0)+gauss_3_pontos(m_d,n,n+2,0))

    elif i == n+1:
        
        return h/2 * (gauss_3_pontos(m_d,n-1,n,0)+gauss_3_pontos(m_d,n,n,0)) - 2*h*gauss_3_pontos(m_d,n,n+2,0)

def matrix_D():
    
    n = var.N
    
    D = []
    
    
    for i in range(1,n+1):
        
        D.append(di(i))
    
    return D
  