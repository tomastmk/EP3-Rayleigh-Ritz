def g(i,x):
    
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
            
def derivada_B(i,x):
    
    t = (x-a-i*h)/h
    
    if t<=-2:   return 0
    
    elif -2<t and t<=-1:   return 3/4 * (2+t)**2
        
    elif -1<t and t<=0:   return 3/4 * (2+t)**2 - 3*(1+t)**2
        
    elif 0<t and t<=1:    return -3/4 * (2-t)**2 + 3*(1-t)**2
        
    elif 1<t and t<=2:   return -3/4 * (2-t)**2
        
    else:   return 0

def m(i,u):
    
    
    
    x = (u * (b-a) + b + a)/2        
    
    Q1 = f(x)
    
    Q2 = g(i,x)
    
    return Q1*Q2

def di(i,x):
    
    integral = ( (5/9) * m(i,-(0,6)**(1/2)) + (8/9)*m(i,0) + (5/9) * m(i, (0,6)**(1/2)) )
    
    return integral
    
