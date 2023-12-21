from solve_SDP.main import system_solve


def matrix_D():
    
    x = X
    n = N
    h = H
    
    D = []
    
    
    for i in range(n):
        
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
    for i in range(N):
        
        Q3 = p( (x[i]+x[i-1])/2 ) / h
        Q4 = h/4 * q( (x[i]+x[i-1])/2 )
        
        Q5 = p( (x[i]+x[i+1])/2 ) / h
        Q6 = h/4 * q( (x[i]+x[i+1])/2 )

        A[0].append( Q3+Q4+Q5+Q6 )
        
    # Segunda Diagonal Principal
    for i in range(n-1):
        
        Q7 = p( (x[i]+x[i+1])/2 ) / h
        Q8 = q( (x[i]+x[i+1])/2 ) / 2
        
        A[1].append( Q7+Q8)
    
    return A


def f(x):
    return -1

def p(x):
    return 1

def q(x):
    return 0
      
      
def main():
    
    global N
    global X
    global H
    
    entrada = [0,1,4]
    #list(map(int,input('Digite a, b e n no formato "a,b,n" \n\n').split(",")))
    
    
    a = entrada[0]
    b = entrada[1]
    N = entrada[2]-1
    
    
    H = (b-a)/(N+1)
    X = [ a+i*H for i in range(N+2) ]

    
    A = matrix_A()
    D = matrix_D()
    c = system_solve(A,D,1)
    
    i = int(input("Digite i:  "))
    v_barra_x_i = c[i]
    
    print("\nA = ",A)
    print("\nD = ",D)
    print("\nc = ",c,"\n"),

    print("reposta final:  ",v_barra_x_i,"\n")
    
main()