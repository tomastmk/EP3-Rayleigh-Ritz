from solve_SDP.main import system_solve




def matrix_A():

    N = n
    H = h
    x = X
    
    A = []
    
    # Primeira Diagonal Principal
    for i in range(N):
        temp0 = p( (x[i-1]+x[i])/2 )/H**2
        temp1 = q( (x[i-1]+x[i])/2 )/4
        temp2 = p( (x[i+1]+x[i])/2 )/H**2
        temp3 = q( (x[i+1]+x[i])/2 )/4

        A.append(temp0+temp1+temp2+temp3)
        
    # Segunda Diagonal Principal
    for i in range(N-1):
        A.append(p( (x[i+1]+x[i])/2 )/H**2 + 3/4*q(  (x[i+1]+x[i])/2  ))
    
    return A


def matrix_D():
    
    x = X
    N = n    
    
    D = []
    
    for i in range(N):
        D.append( (f( (x[i-1]+X[i])/2 )  -  f( (x[i+1]+x[i])/2 ))/2 )
    
    return D
  


def f(x):
    return -1

def p(x):
    return 1

def q(x):
    return 0
      
        
def main():
    
    global n
    global X
    global h
    
    entrada = list(map(int,input('Digite a, b e n no formato "a,b,n" \n\n').split(",")))
    
    
    a = entrada[0]
    b = entrada[1]
    n = entrada[2]
    h = (b-a)/(n+1)
    X = [b-i*h for i in range(n+1)]
    

    
    A = matrix_A()
    D = matrix_D()
    
    c = system_solve(A,D,1)
    
main()