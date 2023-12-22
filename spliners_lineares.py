from solve_SDP.main import system_solve
import matplotlib.pyplot as plt
import numpy as np

def matrix_D():
    
    x = X
    n = N
    h = H
    
    D = []
    
    
    for i in range(1,n):
        
        Q1 = (h/2)*(  f( (x[i]+x[i-1]) /2 )  )
        Q2 = (h/2)*(  f( (x[i]+x[i+1]) /2 )  )
    
        D.append( Q1+Q2 )
    
    return D
  

def matrix_A():

    n = N
    h = H
    x = X
    
    A = [ [],[] ]
    
    # Primeira Diagonal Principal
    for i in range(1,n):
        
        Q3 = p( (x[i]+x[i-1])/2 ) / h
        Q4 = h/4 * q( (x[i]+x[i-1])/2 )
        
        Q5 = p( (x[i]+x[i+1])/2 ) / h
        Q6 = h/4 * q( (x[i]+x[i+1])/2 )

        A[0].append( Q3+Q4+Q5+Q6 )
        
    # Segunda Diagonal Principal
    for i in range(n):
        
        Q7 = p( (x[i]+x[i+1])/2 ) / h
        Q8 = q( (x[i]+x[i+1])/2 ) / 2
        
        A[1].append( Q7+Q8 )
    
    return A


def f(x):
    return -1 #2*np.pi**2*np.sin(np.pi*x)

def p(x):
    return 1

def q(x):
    return 0 #np.pi**2


def v_barra(x,i,c1,c2):
    return (X[i+1]-x)/H*c1 + (x-X[i])/H*c2


def erro(u,v):
    x = X
    
    temp = []
    
    for xi in x:
        temp.append(abs(u(xi)-v(xi)))
        
    temp = 0
    for erro in temp:
        if erro>temp:
            temp = erro
        
    return temp

      
def main():
    
    global N
    global X
    global H
    
    entrada = [0,1,10]
    #list(map(int,input('Digite a, b e n no formato "a,b,n" \n\n').split(",")))
    
    
    a = entrada[0]
    b = entrada[1]
    N = entrada[2]
    
    
    H = (b-a)/(N)
    X = [ a+i*H for i in range(N+1) ]

    
    A = matrix_A()
    D = matrix_D()

    c = system_solve(A,D,1)
    
    c.insert(0,0)
    c.append(0)
    
    
    print("\nA = ",A)
    print("\nD = ",D)
    print("\nc = ",c,"\n")
    
    v = []

    
    for i in range(len(X)-2):
        
        v.append(v_barra(X[i],i,c[i],c[i+1]))
        
    v.append(0)
    
    C = list(map(lambda a: a,c))
    
    eixo_x = np.linspace(0,1,10)
    eixo_y = list(map(lambda a: 0.5*a*(1-a),eixo_x))#list(map(lambda a: np.sin(np.pi*a),eixo_x))
    plt.plot(eixo_x,eixo_y)
    
    plt.scatter(eixo_x,v,color="g")
    
    plt.scatter(X,C,color="r")


    
    plt.show()
    
    
main()