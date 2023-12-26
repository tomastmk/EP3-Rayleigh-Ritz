import splines as var

def B(i,x):
    
    h = var.H
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
    h = var.H
    t = (x-a-i*h)/h
    
    if t<=-2:   return 0
    
    elif -2<t and t<=-1:   return 3/4 * (2+t)**2
        
    elif -1<t and t<=0:   return 3/4 * (2+t)**2 - 3*(1+t)**2
        
    elif 0<t and t<=1:    return -3/4 * (2-t)**2 + 3*(1-t)**2
        
    elif 1<t and t<=2:   return -3/4 * (2-t)**2
        
    else:   return 0


def g(i,x):
    
    n = var.N
    
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
    
    n = var.N
    
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
           