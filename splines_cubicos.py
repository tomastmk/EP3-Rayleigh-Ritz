def aij(i,j): # j = 0,1,2,3
    
    soma = 0
    for l in range(4-j): 
        
        var = x[i+1-l]+x[i+2-l]
        f1 = -(0.6)**(1/2)*h        
        
        
        P1 = 5/9 * p( (f1 + var)/2 ) * derivada_g( i, (f1 + var)/2 ) * derivada_g( i+j, (f1 + var)/2 )
        P2 = 8/9 * p( (var)/2 ) * derivada_g( i, (var)/2 ) * derivada_g( i+j, (var)/2 )
        P3 = 5/9 * p( (-f1+var)/2 ) * derivada_g( i, (-f1+var)/2 ) * derivada_g( i+j, (-f1 + var)/2 )
        
        Q1 = 5/9 * q( (f1+var)/2 ) * g( i, (f1+var)/2 ) * g( i+j, (f1+var)/2 )
        Q2 = 8/9 * q( (var)/2 ) * g( i, (var)/2 ) * g( i+j, (var)/2 )
        Q3 = 5/9 * q( (-f1+var)/2 ) * g( i, (-f1+var)/2 ) * g( i+j, (-f1+var)/2 )
        
        
        soma += P1+P2+P3 + Q1+Q2+Q3
        
    return h/2 * soma


def di(i):
    
    f1 = -(0.6)**(1/2)*h 
    soma = 0    
    
    if i == 0:
                        
        for l in range(0,2):
            var = x[i+l]+x[i+l+1]
        
            P1 = 5/9 * f( (f1 + var)/2 ) * B( 0, (f1 + var)/2 ) 
            P2 = 8/9 * f( (var)/2 ) * B( 0, (var)/2 ) 
            P3 = 5/9 * f( (-f1+var)/2 ) * B( 0, (-f1+var)/2 ) 
            
            soma += P1+P2+P3
        
        Q1 = 5/9 * f( (f1 + x[1] + x[0])/2 ) * B( -1, (f1 + x[1] + x[0])/2 ) 
        Q2 = 8/9 * f( (x[1] + x[0])/2 ) * B( -1, (x[1] + x[0])/2 ) 
        Q3 = 5/9 * f( (-f1 + x[1] + x[0])/2 ) * B( -1, (-f1 + x[1] + x[0])/2 ) 
        
        return h/2 * soma - 2*h*(Q1+Q2+Q3)
    
    if i == 1:
                        
        for l in range(0,3):
            var = x[i+l]+x[i+l+1]
        
            P1 = 5/9 * f( (f1 + var)/2 ) * B( 1, (f1 + var)/2 ) 
            P2 = 8/9 * f( (var)/2 ) * B( 1, (var)/2 ) 
            P3 = 5/9 * f( (-f1+var)/2 ) * B( 1, (-f1+var)/2 ) 
            
            soma += P1+P2+P3
        
        Q1 = 5/9 * f( (f1 + x[1] + x[0])/2 ) * B( -1, (f1 + x[1] + x[0])/2 ) 
        Q2 = 8/9 * f( (x[1] + x[0])/2 ) * B( -1, (x[1] + x[0])/2 ) 
        Q3 = 5/9 * f( (-f1 + x[1] + x[0])/2 ) * B( -1, (-f1 + x[1] + x[0])/2 ) 
        
        return h/2 * soma - 2*h*(Q1+Q2+Q3)
    
    
    elif 2<=i<=n-1:
                
        for l in range(-2,2):
            var = x[i-l]+x[i-l+1]
        
            P1 = 5/9 * f( (f1 + var)/2 ) * B( i, (f1 + var)/2 ) 
            P2 = 8/9 * f( (var)/2 ) * B( i, (var)/2 ) 
            P3 = 5/9 * f( (-f1+var)/2 ) * B( i, (-f1+var)/2 ) 
            
            soma += P1+P2+P3
        
        return h/2 * soma
    
    elif i == n:
                
        for l in range(-2,1):
            var = x[i-l]+x[i-l+1]
        
            P1 = 5/9 * f( (f1 + var)/2 ) * B( n, (f1 + var)/2 ) 
            P2 = 8/9 * f( (var)/2 ) * B( n, (var)/2 ) 
            P3 = 5/9 * f( (-f1+var)/2 ) * B( n, (-f1+var)/2 ) 
            
            soma += P1+P2+P3
        
        Q1 = 5/9 * f( (f1 + x[n] + x[n+1])/2 ) * B( n+2, (f1 + x[n] + x[n+1])/2 ) 
        Q2 = 8/9 * f( (x[n] + x[n+1])/2 ) * B( n+2, (x[n] + x[n+1])/2 ) 
        Q3 = 5/9 * f( (-f1 + x[n] + x[n+1])/2 ) * B( n+2, (-f1 + x[n] + x[n+1])/2 ) 
        
        return h/2 * (soma + Q1+Q2+Q3)
    
    elif i == n+1:
                
        for l in range(-1,1):
            var = x[n-l]+x[n-l+1]
        
            P1 = 5/9 * f( (f1 + var)/2 ) * B( n, (f1 + var)/2 ) 
            P2 = 8/9 * f( (var)/2 ) * B( n, (var)/2 ) 
            P3 = 5/9 * f( (-f1+var)/2 ) * B( n, (-f1+var)/2 ) 
            
            soma += P1+P2+P3
        
        Q1 = 5/9 * f( (f1 + x[n] + x[n+1])/2 ) * B( n+2, (f1 + x[n] + x[n+1])/2 ) 
        Q2 = 8/9 * f( (x[n] + x[n+1])/2 ) * B( n+2, (x[n] + x[n+1])/2 ) 
        Q3 = 5/9 * f( (-f1 + x[n] + x[n+1])/2 ) * B( n+2, (-f1 + x[n] + x[n+1])/2 ) 
        
        return h/2 * soma  -  2*h*(Q1+Q2+Q3)
    
    
      
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