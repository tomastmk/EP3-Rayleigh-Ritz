import sympy as sp


t2 = sp.symbols('t2')
t1 = sp.symbols('t1')


    
    

E1 = 1/4*((2-t1)**3-4*(1-t1)**3)-(2-t2)**3


sp.solve(E1,t1=4,t2=3)