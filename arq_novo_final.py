def gauss_3_pontos(m,xi,i,j):
    
    h = H
     
    P1 = 5/9 * m( -0.6**(1/2), xi, i ,j)
    P2 = 8/9 * m( 0, xi, i )
    P3 = 5/9 * m( 0.6**(1/2), xi, i , j)
    
    return P1 + P2 + P3


def m_d(u , xi, bi, j):
    
    h = H
    x = X    
    
    return f( u*h/2 + x[xi]/2 + x[xi+1]/2 )*B( bi , u*h/2 + x[xi]/2 + x[xi+1]/2 )
    
    
def di(i):
    
    h = h
    n = N
    
    if i == 0:
        
        return h/2 * (gauss_3_pontos(m_d, 0, 0) + gauss_3_pontos(m_d, 1 , 0)) - 2*h*gauss_3_pontos(m_d,0,-1)
    
    elif i == 1:
        
        return h/2 * (gauss_3_pontos(m_d,0,1)+gauss_3_pontos(m_d,1,1)+gauss_3_pontos(m_d,2,1)-gauss_3_pontos(m_d,0,-1))
    
    elif 2 <= i <= n-1:
        
        return h/2 * (gauss_3_pontos(m_d,i-2,i)+gauss_3_pontos(m_d,i-1,i)+gauss_3_pontos(m_d,i,i)+gauss_3_pontos(m_d,i+1,i))
    
    elif i == n:
        
        return h/2 * (gauss_3_pontos(m_d,n-2,n)+gauss_3_pontos(m_d,n-1,n)+gauss_3_pontos(m_d,n,n)+gauss_3_pontos(m_d,n,n+2))

    elif i == n+1:
        
        return h/2 * (gauss_3_pontos(m_d,n-1,n)+gauss_3_pontos(m_d,n,n)) - 2*h*gauss_3_pontos(m_d,n,n+2)
    
    

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
    
    return q( u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g(i, u*h/2 + x[xi]/2 + x[xi+1]/2 ) * g(j, u*h/2 + x[xi]/2 + x[xi+1]/2 )


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
    
    if i == N:
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
