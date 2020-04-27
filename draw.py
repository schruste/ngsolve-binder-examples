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
    Draw = PyVistaDraw
    print("Overwriting NGSolve's 'Draw' function with PyVista-Drawing")
    print("Reset to GUI Draw with:")
    print("  from  netgen import  gui # <- spawn GUI")
    print("  from ngsolve import Draw # <- take GUI version of Draw")
except ImportError:
    from netgen import gui
    from ngsolve import Draw
    print("Spawning NGSolve's GUI with standard 'Draw' function")


