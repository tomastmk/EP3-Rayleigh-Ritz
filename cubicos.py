def gauss_3_pontos(m,xi,i,j):
    
    h = H
     
    P1 = 5/9 * m( -0.6**(1/2), xi, i ,j)
    P2 = 8/9 * m( 0, xi, i , j)
    P3 = 5/9 * m( 0.6**(1/2), xi, i , j)
    
    return P1 + P2 + P3


def m_d(u , xi, bi, j):
    
    h = H
    x = X    
    
    
    return f( u*h/2 + x[xi]/2 + x[xi+1]/2 )*B( bi , u*h/2 + x[xi]/2 + x[xi+1]/2 )
    
    
def di(i):
    
    h = H
    n = N
    
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
    
    

def B(i,x):
    
    h = H
    a = 0

    t = (x-a-i*h)/h
    
    if t<=-2:   return 0
    
    elif -2<t and t<=-1:   return (1/4) * (2+t)**3
        
    elif -1<t and t<=0:   return (1/4) * ( (2+t)**3 - 4*(1+t)**3 )
        
    elif 0<t and t<=1:    return (1/4) * ( (2-t)**3 - 4*(1-t)**3 )
        
    elif 1<t and t<=2:   return (1/4) * (2-t)**3
        
    else:   return 0
    


def m_ap(u, xi, i, j):
    
    h = H
    x = X
    
    return p( u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g_derivada(i, u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g_derivada(j, u*h/2 + x[xi]/2 + x[xi+1]/2 )

def m_aq(u, xi, gi, gj):
    
    h = H
    x = X
    
    return q( u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g(gi, u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g(gj, u*h/2 + x[xi]/2 + x[xi+1]/2 )


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
    
    if i >= N:
        P4 = 0
        Q4 = 0
    
    else: 
        P4 = gauss_3_pontos(m_ap,i+1,i,j)
        Q4 = gauss_3_pontos(m_aq,i+1,i,j)
    
    return H/2 * (P1+P2+P3+P4 + Q1+Q2+Q3+Q4)



def g(i,x):
    
    n = N
    
    if i == 0:
        return B(0,x) - 4*B(-1,x)
        
    elif i == 1:
        return B(1,x) - B(-1,x)
    
    elif 2<=i and i<=n-1:
        return B(i,x)
        
    elif i == n:
        return B(n,x) - B(n+2,x)

    elif i ==n+1:
        return B(n+1,x) - 4*B(n+2,x)
    
  
def g_derivada(i,x):
    
    n = N
    
    if i == 0:
        return derivada_B(0,x) - 4*derivada_B(-1,x)
        
    elif i == 1:
        return derivada_B(1,x) - derivada_B(-1,x)
    
    elif 2<=i and i<=n-1:
        return derivada_B(i,x)
        
    elif i == n:
        return derivada_B(n,x) - derivada_B(n+2,x)

    elif i ==n+1:
        return derivada_B(n+1,x) - 4*derivada_B(n+2,x)    
           
           
def derivada_B(i,x):
    a = 0
    h = H
    t = (x-a-i*h)/h
    
    if t<=-2:   return 0
    
    elif -2<t and t<=-1:   return 3/4 * (2+t)**2
        
    elif -1<t and t<=0:   return 3/4 * (2+t)**2 - 3*(1+t)**2
        
    elif 0<t and t<=1:    return -3/4 * (2-t)**2 + 3*(1-t)**2
        
    elif 1<t and t<=2:   return -3/4 * (2-t)**2
        
    else:   return 0


from solve_SDP.main import system_solve  #EP2
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

class bcolors:
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    RED = "\033[0;31m"
    BROWN = "\033[0;33m"
    PURPLE = "\033[0;35m"
    CYAN = "\033[0;36m"
    BOLD = "\033[1m"
    ENDC = "\033[0m"

def matrix_D():
    
    x = X
    n = N
    h = H
    
    D = []
    
    
    for i in range(1,n+1):
        
        D.append(di(i))
    
    return D
  

def matrix_A():

    n = N
    h = H
    x = X
    
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


def f(x):
    return x+(2-x)*np.exp(x) #2*np.pi**2*np.sin(np.pi*x)

def p(x):
    return np.exp(x)

def q(x):
    return np.exp(x) #np.pi**2


def v_barra(x,i,c1,c2):
    return (X[i+1]-x)/H*c1 + (x-X[i])/H*c2


def erro(u,v):

    erros = []
    
    for ui,vi in zip(u,v):
        erro = abs(ui-vi)
        erros.append(erro)
        
    max = 0
    
    for erro in erros:
        if erro>max:
            max = erro
        
    return max

      
def main():
    
    global N
    global X
    global H
    
    entrada = [0,1,2]
    #list(map(int,input('Digite a, b e n no formato "a,b,n" \n\n').split(",")))
    
    a = entrada[0]
    b = entrada[1]
    N = entrada[2]
    
    H = (b-a)/(N+1)
    X = [ a+i*H for i in range(N+2) ]

    
    A = matrix_A()
    D = matrix_D()

    c = system_solve(A,D,1)
    
    # Adiciona c_0 e c_(n+1)
    c.insert(0,0)
    c.append(0)
    
    print("\nA = ",A)
    print("\nD = ",D)
    print("\nc = ",c)

    # u(x)
    eixo_x = np.linspace(0,1,N+1)
    eixo_y = list(map(lambda a: (a-1)*(np.exp(-a)-1),eixo_x))#list(map(lambda a: np.sin(np.pi*a),eixo_x))
    plt.plot(eixo_x,eixo_y)

    # Plot v_barra
    plt.scatter(X,c,color="r")
    plt.legend("uv")
    
    plt.xlabel('x')
    plt.ylabel('y')

    print(bcolors.BOLD+"\nErro máximo: "+bcolors.ENDC,erro(eixo_y,c),"\n")

    plt.show()
    
    
    
    
main()