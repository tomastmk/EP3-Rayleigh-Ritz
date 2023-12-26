from solve_SDP.main import system_solve  #EP2
import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return x +((2-x)*np.exp(x))

def p(x):
    return np.exp(x)

def q(x):
    return np.exp(x)

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

def particiona_intervalo(a,b):
    '''Divide o intervalo em n+1 subintervalos'''

    x = []

    for i in range(N+1):
        xi = a + i*H
        x.append(xi)
    
    x.append(b)

    return x




# VETOR D:

def MV(i,j,u):
    '''Mudança de Variável (MV)'''
    
    x = ((X[j]-X[j-1])*u + X[j-1] + X[j])/2
    #x = (2*u - X[j] - X[j-1])/(X[j]-X[j-1])

    F1 = f(x)
    F2 = g(i,x)

    return F1*F2

def gauss_approx_3pts(i,j):

    integral = ( (5/9) *    MV(i,j,-(0.6)**(1/2)) + (8/9)*MV(i,j,0) + (5/9) * MV(i,j, (0.6)**(1/2)) )
    
    return ((X[j] - X[j-1])/2) * integral

def D():
    '''Calcula todos os coeficientes di
    
        - 1º for: itera os di
        - 2º for: calcula a integral de a até b, com subintegrais nos intervalos [xi-1, xi] para garantir a regularidade do integrando
    '''
    D_vector = []

    for i in range(N): #CONFERIR  
        sum = 0
        for j in range(1, N+2):
            sum += gauss_approx_3pts(i,j)
        D_vector.append(sum)

    return D_vector




# MATRIZ A:

def MV_A(i,j,k,u):

    x = ((X[k]-X[k-1])*u + X[k-1] + X[k])/2
    #x = (2*u - X[j] - X[j-1])/(X[j]-X[j-1])

    F1 = p(x)
    F2 = derivada_g(i,x)
    F3 = derivada_g(j,x)
    F4 = q(x)
    F5 = g(i,x)
    F6 = g(j,x)

    return F1*F2*F3 + F4*F5*F6

def gauss_approx_3pts_A(i,j,k):

    integral = ( (5/9) *    MV_A(i,j,k,-(0.6)**(1/2)) + (8/9)*MV_A(i,j,k, 0) + (5/9) * MV_A(i,j,k, (0.6)**(1/2)) )
    
    return ((X[k] - X[k-1])/2) * integral

def A():

    aii = [] #diagonal principal
    aij1 = [] #segunda diagonal principal
    aij2 = []
    aij3 = []

    for i in range(N):
        for j in range(N):
            sum = 0
            if abs(i-j) >= 4: 
                sum += 0
            elif i == j:
                for k in range(1, N+2):
                    sum += gauss_approx_3pts_A(i,j,k)
                aii.append(sum)
            elif j == i+1:
                for k in range(1, N+2):
                    sum += gauss_approx_3pts_A(i,j,k)
                aij1.append(sum)
            elif j == i+2:
                for k in range(1, N+2):
                    sum += gauss_approx_3pts_A(i,j,k)
                aij2.append(sum)
            elif j == i+3:
                for k in range(1, N+2):
                    sum += gauss_approx_3pts_A(i,j,k)
                aij3.append(sum)
            elif i == j+1:
                sum += 0
            elif i == j+2:
                sum += 0
            elif i == j+3:
                sum += 0

    return [aii, aij1, aij2, aij3]





# Aproximação v barra:

def v_barra(c, x):

    sum = 0
    for i in range(N):
        sum += (c[i]*g(i,x))

    return sum





def main():

    global H 
    global N 
    global a
    global b
    global X

    a, b, N = 0, 1, 10
    H = (b-a)/(N+1)

    X = particiona_intervalo(a,b)
    D_vector = D()
    print()
    A_matrix = A()
 
    print(D_vector)
    print("MATRIZ A")
    for linha in A_matrix:
        print(linha)

    c = system_solve(A_matrix,D_vector,1)

    print()
    print(c)

    eixo_x = np.linspace(0,1,100)
    plt.plot(eixo_x, (eixo_x-1)*(np.exp(-eixo_x) - 1))

    v = []
    for x in eixo_x:
        v.append(v_barra(c,x))

    plt.scatter(eixo_x, v, color='r')
    plt.show()

main()