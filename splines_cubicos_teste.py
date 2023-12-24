from solve_SDP.main import system_solve  #EP2
import matplotlib.pyplot as plt
import numpy as np

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
    
def derivada_g(i,x):
    
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
    
def m(i,u):
    
    a = 0
    b = 1
    
    x = (u * (b-a) + b + a)/2        
    
    Q1 = f(x)
    
    Q2 = g(i,x)
    
    return Q1*Q2

def di(i):
    
    integral = ( (5/9) * m(i,-(0.6)**(1/2)) + (8/9)*m(i,0) + (5/9) * m(i, (0.6)**(1/2)) )
    
    return 1/2*integral


def matrix_D():
    
    n = N
    
    D = []
    
    for i in range(n):
        D.append(di(i))
        
    return D
  

  
def l(i,j,u):
    
    a = 0
    b = 1
    
    x = (u * (b-a) + b + a)/2        
    
    Q1 = p(x)
    
    Q2 = derivada_g(i,x)
    
    Q3 = derivada_g(j,x)
    
    Q4 = q(x)
    
    Q5 = g(i,x)
    
    Q6 = g(j,x)
    
    return 1/2*( Q1*Q2*Q3 + Q4*Q5*Q6)



def matrix_A():

    n = N
    h = H
    x = X
    
    A = [[0 for i in range(n)] for i in range(n)]
    
    for i in range(n):
        for j in range(n):
            integral = ( (5/9) * l(i,j,-(0.6)**(1/2)) + (8/9)*l(i,j,0) + (5/9) * l(i,j, (0.6)**(1/2)) )

            A[i][j] = 1/2*integral
    
    return A


def f(x):
    return (x +(2-x)*np.exp(x))

def p(x):
    return np.exp(x)

def q(x):
    return np.exp(x)


def v_barra(x,i,c1,c2):
    return (X[i+1]-x)/H*c1 + (x-X[i])/H*c2


def erro(u,v):
    x = X
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
    
    entrada = [0,1,100]
    #list(map(int,input('Digite a, b e n no formato "a,b,n" \n\n').split(",")))
    
    a = entrada[0]
    b = entrada[1]
    N = entrada[2]
    
    H = (b-a)/(N)
    X = [ a+i*H for i in range(N+1) ]

    
    A = matrix_A()
    D = matrix_D()

    c = system_solve(A,D,1)
    
    # Adiciona c_0 e c_(n+1)
    c.insert(0,0)
    c.append(0)
    
    print(D)
    
    '''
    print("\nA = ",A)
    print("\nD = ",D)
    print("\nc = ",c)
    '''
    '''
    # u(x)
    eixo_x = np.linspace(0,1,N+1)
    eixo_y = list(map(lambda a: 0.5*a*(1-a),eixo_x))#list(map(lambda a: np.sin(np.pi*a),eixo_x))
    plt.plot(eixo_x,eixo_y)

    # Plot v_barra
    plt.scatter(X,c,color="r")
    plt.legend("uv")
    
    plt.xlabel('x')
    plt.ylabel('y')

    print(bcolors.BOLD+"\nErro máximo: "+bcolors.ENDC,erro(eixo_y,c),"\n")
    
    plt.show()
    '''
    
    
main()