
def erro(u,v):

    erros = []
    
    for ui,vi in zip(u,v):
        erro = abs(ui-vi)
        erros.append(erro)
        
    max = 0
    
    for erro in erros:
        if erro>max:
            max = erro
        
    return max
