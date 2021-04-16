#from pyvista import set_plot_theme
#set_plot_theme('document')

from ngsolve import VTKOutput, CoefficientFunction, Mesh
# from IPython import get_ipython
# ipython = get_ipython()
# ipython.magic("matplotlib inline")
try:
    import pyvista as pv
    def PyVistaDraw(obj, mesh=None, name="dummy", sd=0, filename="tmp", elevation = 1, static=False, show_edges = True):
        if isinstance(obj, CoefficientFunction):
            gfu = obj
            if mesh == None:
                mesh = gfu.space.mesh
            if gfu.dim == 2:
                gfu = CoefficientFunction((gfu[0],gfu[1],0))

            vtkout = VTKOutput(mesh, [gfu], [name], 
                               filename = filename, subdivision=sd)
            vtkout.Do()
            plot_mesh = pv.read(filename+".vtk")

            
            p = pv.Plotter(notebook=True)

            if gfu.dim == 1 and elevation != 0:
                plot_mesh = plot_mesh.warp_by_scalar(name,normal=(0,0,1),factor=elevation)

            p.add_mesh(plot_mesh, scalars=name, label=name)

            if show_edges:
                edges = plot_mesh.extract_feature_edges(20)
                p.add_mesh(edges, color="black") 
            
            if gfu.dim > 1:
                glyphs = plot_mesh.glyph(orient=name, scale=name, factor=1, geom=pv.Arrow())
                p.add_mesh(glyphs, color="black")

            return p.show(title=name,use_panel=not static, auto_close=False)
        elif isinstance(obj,Mesh):
            mesh = obj
            vtkout = VTKOutput(mesh, [], [], 
                               filename = filename, subdivision=sd)
            vtkout.Do()
            plot_mesh = pv.read(filename+".vtk")
            p = pv.Plotter(notebook=True)
            edges = plot_mesh.extract_feature_edges(20)
            p.add_mesh(plot_mesh, show_edges=False) 
            p.add_mesh(edges, color="black") 
            return p.show(title=name,use_panel=not static, auto_close=False)
        else:
            raise Exception("obj-type is not viable")
        
    # Draw = PyVistaDraw
    print("* PyVista-Drawing with 'PyVistaDraw' available")
except ImportError:
    try:
        from ngsolve.webgui import Draw
    except ImportError:
        from netgen import gui
        from ngsolve import Draw
        print("Spawning NGSolve's GUI with standard 'Draw' function")
    print("Using NGSolve's 'Draw' from webgui!")

try:
    from ngsolve.webgui import Draw
    print("Using NGSolve's 'Draw' from webgui!")
except ImportError:
    from netgen import gui
    from ngsolve import Draw
    print("Spawning NGSolve's GUI with standard 'Draw' function")

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
    plt.plot(x_v,[miny for v in x_v],'|',label='vertices')
    plt.xlabel("x")
    plt.legend()
    plt.show()
    if keep:
        display.clear_output(wait=True)
    
print("* 1D Drawing available with patched 'Draw' or simply 'Draw1D'")

oldDraw = Draw
def Draw(cf,mesh,label,*args,**kwargs):
    if mesh.dim == 1:
        Draw1D(mesh,[(cf,label)],*args,**kwargs)
    else:
        oldDraw(cf,mesh,label,*args,**kwargs)
