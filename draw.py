from pyvista import set_plot_theme
set_plot_theme('document')

from ngsolve import VTKOutput

# from IPython import get_ipython
# ipython = get_ipython()
# ipython.magic("matplotlib inline")

import pyvista as pv
def PyVistaDraw(gfu, mesh, name, sd=0, filename="tmp"):
    if gfu.dim > 1:
        raise Exception("no vectorial things yet")
    vtkout = VTKOutput(mesh, [gfu], [name], 
                       filename = filename, subdivision=sd)
    vtkout.Do()
    plot_mesh = pv.read(filename+".vtk")
    p = pv.Plotter(notebook=True)
    p.add_mesh(plot_mesh.warp_by_scalar(name,normal=(0,0,1),factor=1), scalars=name, label=name) 
    return p.show(title=name,use_panel=True, auto_close=False)    
