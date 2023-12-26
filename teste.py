from splines.cubic_matrix import matrix_A,matrix_D
from splines.v_barra_cubico import v_barra
from splines.erro import erro

from solve_SDP.main import system_solve  #EP2
import numpy as np
import matplotlib.pyplot as plt
import splines as var


class bcolors:
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    RED = "\033[0;31m"
    BROWN = "\033[0;33m"
    PURPLE = "\033[0;35m"
    CYAN = "\033[0;36m"
    BOLD = "\033[1m"
    ENDC = "\033[0m"

N = var.N

D = matrix_D()
A = matrix_A()
c = system_solve(A,D,3)
c.insert(0,0)
c.append(0)










# u(x)
eixo_x = np.linspace(0,1,N+1)
eixo_y = list(map(lambda a: (a-1)*(np.exp(-a)-1),eixo_x))#list(map(lambda a: np.sin(np.pi*a),eixo_x))
plt.plot(eixo_x,eixo_y)

linha = []
for i in range(2001):
    linha.append(i*0.0005)
    
aprox = []
for x in linha:
    aprox.append(v_barra(x,c))

# Plot v_barra
plt.scatter(linha,aprox,color="r")
plt.legend("uv")

plt.xlabel('x')
plt.ylabel('y')

print(bcolors.BOLD+"\nErro máximo: "+bcolors.ENDC,erro(eixo_y,c),"\n")

plt.show()