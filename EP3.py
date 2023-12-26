from solve_SDP.main import system_solve  #EP2
import matplotlib.pyplot as plt
import numpy as np

class bcolors:
    BOLD = "\033[1m"
    ENDC = "\033[0m"


# Definindo funções do problema
def f(x):
    return x+(2-x)*np.exp(x) #2*np.pi**2*np.sin(np.pi*x)

def p(x):
    return np.exp(x)

def q(x):
    return np.exp(x) #np.pi**2

def u(x):
    return (x-1)*(np.exp(-x)-1)


'''

def f(x):
    return 1 

def p(x):
    return 1

def q(x):
    return 0 
    
def u(x):
    return 0.5*x*(1 - x)

'''

# Matrizes 
def matrix_D(linear = True):
    
    x = X
    n = N
    h = H
    
    D = []
    
    if linear is True:
        for i in range(1,n+1):
            
            Q1 = (h/2)*(  f( (x[i]+x[i-1]) /2 )  )
            Q2 = (h/2)*(  f( (x[i]+x[i+1]) /2 )  )
        
            D.append( Q1+Q2 )
    else:
        for i in range(1,n+1):
            D.append(di(i))
            
    return D
  
def matrix_A(linear = True):

    n = N
    h = H
    x = X
    
    if linear is True:
        
        A = [ [],[] ]
        
        # Primeira Diagonal Principal
        for i in range(1,n+1):
            
            Q3 = p( (x[i]+x[i-1])/2 ) / h
            Q4 = h/4 * q( (x[i]+x[i-1])/2 )
            
            Q5 = p( (x[i]+x[i+1])/2 ) / h
            Q6 = h/4 * q( (x[i]+x[i+1])/2 )

            A[0].append( Q3+Q4+Q5+Q6 )
            
        # Segunda Diagonal Principal
        for i in range(1,n):
            
            Q7 = -p( (x[i]+x[i+1])/2 ) / h
            Q8 = q( (x[i]+x[i+1])/2 ) / 2
            
            A[1].append( Q7+Q8 )
    
    else:
        A = [ [],[],[],[] ]
        
        for i in range(1,n):
            A[0].append(aij(i,i))
        for i in range(1,n-1):
            A[1].append(aij(i,i+1))
        for i in range(1,n-2):
            A[2].append(aij(i,i+2))
        for i in range(1,n-3):
            A[3].append(aij(i,i+3))
        
    return A


# Calcula aproximação de u em x
def v_barra(x,c):
    
    n = N
    
    soma = 0
    for i in range(n):
        soma += c[i]*g(i,x)
    
    return soma


# Calcula maior erro entre funções u e v
def erro(u,v):

    max = 0
    for ui,vi in zip(u,v):
        
        erro = abs(ui-vi)
        
        if erro>max:
            max = erro
        
    return max


# Splines cúbicos
def gauss_3_pontos(m,xi,i,j):
    
    h = H
     
    P1 = 5/9 * m( -0.6**(1/2), xi, i ,j)
    P2 = 8/9 * m( 0, xi, i , j)
    P3 = 5/9 * m( 0.6**(1/2), xi, i , j)
    
    return P1 + P2 + P3

def m_ap(u, xi, i, j):
    
    h = H
    x = X
    
    return p( u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g_derivada(i, u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g_derivada(j, u*h/2 + x[xi]/2 + x[xi+1]/2 )

def m_aq(u, xi, gi, gj):
    
    h = H
    x = X
    
    return q( u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g(gi, u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g(gj, u*h/2 + x[xi]/2 + x[xi+1]/2 )

def m_d(u , xi, bi, j):
    
    h = H
    x = X    
    
    
    return f( u*h/2 + x[xi]/2 + x[xi+1]/2 )*B( bi , u*h/2 + x[xi]/2 + x[xi+1]/2 )
    
    
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


# main      
def main():
    
    entrada = list(map(int,input(bcolors.BOLD+'\nDigite a, b e n no formato "a,b,n":     '+bcolors.ENDC).split(",")))
    linear = input(bcolors.BOLD+"Spline Linear (S/N):      "+bcolors.ENDC)
    
    if linear == "S" or linear == "s":
        linear = True
    else:
        linear = False
    
    # Definindo cte globais
    a = entrada[0]
    b = entrada[1]
    
    global N;   N = entrada[2]
    global H;   H = (b-a)/(N+1)
    global X;   X = [ a+i*H for i in range(N+2) ]    
    
    
    # Montando e Resolvendo Sistema Ac = D
    if linear is True:
        A = matrix_A()
        D = matrix_D()
        c = system_solve(A,D,1)
        
        # Adiciona c_0 e c_(n+1)
        c.insert(0,0)
        c.append(0)
        
        # Plot v_barra
        plt.scatter(X,c,color="r")
    

    else:
        A = matrix_A(linear = False)
        D = matrix_D(linear = False)
        c = system_solve(A,D,3)
        
        # Adiciona c_0 e c_(n+1)
        c.insert(0,0)
        c.append(0)
        
        # Plot v_barra
        linha = [i*0.0005 for i in range(2001)]
        v_b = [v_barra(x,c) for x in linha]
        plt.scatter(linha,v_b,color="r")


    # u(x)
    eixo_x = np.linspace(0,1,N+2)
    eixo_y = list(map(u,eixo_x))
    plt.plot(eixo_x,eixo_y)


    # Print Erro
    print(bcolors.BOLD+"\nErro máximo: "+bcolors.ENDC,erro(eixo_y,c),"\n")

    # Configurações do Plot
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend("uv")

    plt.show()

    
    
main()