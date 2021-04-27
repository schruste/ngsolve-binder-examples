import matplotlib.pyplot as plt
from numpy import nan

from IPython import display

def Draw1D(mesh, coefs, keep=False, n_p=2, figsize=(20,4)):
    """
        draw coefficient functions with matplotlib
    arguments:
        mesh: ngsolve.comp.Mesh
            Mesh (1D) on that functions should be drawn
        coefs: list of tuples, e.g. [(a,"a"),(b,"b")]
            first component of tuple is assumed to be a coefficient function
            second component is assumed to be a string that is used as a label
        n_p: int
            number of sampling points on every element (minimum is 2)
    """
    if n_p <= 2:
        n_p = 2
        
    eps = 1e-6 
    
    x_v = [p[0] for p in mesh.ngmesh.Points()]
    x_s = []
    f_s = {}

    miny = 1e99
    for f, name in coefs:
        f_s[name] = []
        
    x_s.append(nan)
    for f,name in coefs:
        f_s[name].append(nan)
        
    for el in mesh.ngmesh.Elements1D():
        left = mesh.ngmesh.Points()[el.points[0]][0]
        right = mesh.ngmesh.Points()[el.points[1]][0]
        for l in range(n_p):
            y = left + eps + (l / (n_p-1)) * (right - eps -left) 
            x_s.append(y)
            for f,name in coefs:
                ff = f(mesh(y))
                miny = min(miny,ff)
                f_s[name].append(ff)
                
        x_s.append(nan)
        for f,name in coefs:
            f_s[name].append(nan)

        
    # plt.clf()
    # display.display(plt.gcf())
    plt.figure(figsize=figsize)
    for f,name in coefs:
        plt.plot(x_s,f_s[name],label=name)
    plt.plot(x_v,[miny for v in x_v],'|',label='Knoten')
    plt.xlabel("x")
    plt.legend()
    plt.show()
    if keep:
        display.clear_output(wait=True)
    
