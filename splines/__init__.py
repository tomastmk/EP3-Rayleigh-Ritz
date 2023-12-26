import numpy as np

entrada = [0,1,10]
#list(map(int,input('Digite a, b e n no formato "a,b,n" \n\n').split(",")))

a = entrada[0]
b = entrada[1]
N = entrada[2]

H = (b-a)/(N+1)
X = [ a+i*H for i in range(N+2) ]

def f(x):
    return x+(2-x)*np.exp(x) #2*np.pi**2*np.sin(np.pi*x)

def p(x):
    return np.exp(x)

def q(x):
    return np.exp(x) #np.pi**2
