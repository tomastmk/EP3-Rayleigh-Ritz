import sympy as sp



def matrix_D():
    
    x = X
    n = N
    h = H
    
    D = []
    
    
    for i in range(1,n+1):
        
        soma = 0
        for base in range(0,10000):
            soma += f(base/10000)*g(i,base/10000)*0.0001

        D.append(soma)
        
    return D


def matrix_A():

    n = N
    h = H
    x = X
    
    A = [ [],[],[],[] ]
    
    for i in range(n):
        soma = 0
        for base in range(0,10000):
           soma += ( p(base/10000)*g_derivada(i,base/10000)**2 + q(base/10000)*g(i,base/10000)**2)*0.0001
        A[0].append(soma)
        
    for i in range(n-1):
        soma = 0
        for base in range(0,10000):
           soma += ( p(base/10000)*g_derivada(i,base/10000)*g_derivada(i+1,base/10000) + q(base/10000)*g(i,base/10000)*g(i+1,base/10000))*0.0001
        A[1].append(soma)
        
    for i in range(n-2):
        soma = 0
        for base in range(0,10000):
           soma += ( p(base/10000)*g_derivada(i,base/10000)*g_derivada(i+2,base/10000) + q(base/10000)*g(i,base/10000)*g(i+2,base/10000))*0.0001
        A[2].append(soma)
        
    for i in range(n-3):
        soma = 0
        for base in range(0,10000):
           soma += ( p(base/10000)*g_derivada(i,base/10000)*g_derivada(i+3,base/10000) + q(base/10000)*g(i,base/10000)*g(i+3,base/10000))*0.0001
        A[3].append(soma)
    
    return A



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



def f(x):
    return x+(2-x)*np.exp(x) #2*np.pi**2*np.sin(np.pi*x)

def p(x):
    return np.exp(x)

def q(x):
    return np.exp(x) #np.pi**2


def v_barra(x,c):
    
    n = N
    
    soma = 0
    for i in range(n):
        soma += c[i]*g(i,x)
    
    return soma
    

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
    
    entrada = [0,1,50]
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

    linha = []
    for i in range(2001):
        linha.append(i*0.0005)
        
    aprox = []
    for x in linha:
        aprox.append(v_barra(x,c))

    # Plot v_barra
    plt.scatter(linha,aprox,color="r")
    plt.legend("uv")
    
    plt.xlabel('x')
    plt.ylabel('y')

    print(bcolors.BOLD+"\nErro máximo: "+bcolors.ENDC,erro(eixo_y,c),"\n")

    plt.show()
    
    
    
    
main()