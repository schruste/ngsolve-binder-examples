import matplotlib.pyplot as plt
from random import random
from math import *

def DrawSpecialFunction(N=1):
    if N < 1:
        print("WARNING: Can not draw. Please choose a larger N")
        return

    x_v = [i/N for i in range(0,N+1)]
    plt.figure(figsize=(20, 3))
    plt.plot(x_v,[0 for v in x_v],'|',label='0')

    n = 2*N
    
    x_s = [i/n for i in range(0,n+1)]
    f_s = [0 if i%2==0 else random()-0.5 for i in range(0,n+1)]
    
    plt.plot(x_s,f_s)
    plt.show()
    
def DrawErrorToSmoothFunction(N=10, smooth_function= "sin(2*pi*x)"):
    try:
        f = lambda x: eval(smooth_function)
        f(0.5)
    except:
        print("expression for smooth_function not valid!")
        return
    
    if N < 1:
        print("WARNING: Can not draw. Please choose a larger N")
        return

    plt.figure(figsize=(20, 5))
    plt.subplot(2, 1, 1)
    
    x_s = [i/N for i in range(0,N+1)]
    F_s = [f(x) for x in x_s]
    plt.plot(x_s,F_s,"x-")

    M = 10*N
    x_s = [i/M for i in range(0,M+1)]
    f_s = [f(x) for x in x_s]
    plt.plot(x_s,f_s)
    plt.legend(["Annäherung","Exakte Funktion"])

    
    plt.subplot(2, 1, 2)
    
    m = max([abs(q) for q in f_s]) * 0.2

    x_s = [i/M for i in range(0,M+1)]

    def f_s_approx(x):
        l = floor(x * N)
        r = ceil(x * N)
        if l==r:
            return F_s[l]
        else:
            c = (x*N - l)/(r - l)
            return c * F_s[r] + (1-c) * F_s[l]
    
    f_s = [f(x)-f_s_approx(x) for x in x_s]
    plt.plot(x_s,f_s)

    maxd = max([abs(q) for q in f_s])
    print("maximale Abweichung bei N = ", N, ": ", maxd)
    m = max(1.1*maxd,m)
    plt.ylim((-m,m))
    plt.legend(["Abweichung"])

    # f_s = [0 if i%2==0 else random()-0.5 for i in range(0,n+1)]
    
    plt.show()
    
from mesh1d import *
from draw1d import *

from ngsolve import H1, GridFunction
from ngsolve import x as X
from ngsolve import sin as Sin
from ngsolve import cos as Cos

def DrawBasisFunction(N=4, i=0, order=1):
   mesh1D = Mesh1D(N)
   fes = H1(mesh1D, order=order)
   gf = GridFunction(fes)
   
   gf.vec[:] = 0
   if i >= fes.ndof:
       print(" i is too large. Setting it to", fes.ndof-1)
       i = fes.ndof - 1
       
   gf.vec[i] = 1
   Draw1D(mesh1D, [(gf,"Basisfunktion "+str(i))], n_p=2*order**2)


def ApproximateWithFESpace(N=64, order=3, f="sin(3·pi·x)"):

    if f == "sin(3·pi·x)":
        f = Sin(3*pi*X)
    elif f == "sin(3·pi·x)·cos(5·pi·x)":
        f = Sin(3*pi*X)*Cos(5*pi*X)
    else:
        print("no function provided")
        return
        
    mesh1D = Mesh1D(N)
    try:
        f(mesh1D(0.5))
    except:
        print("expression for f not valid!")
        return
    
    fes = H1(mesh1D, order=order)
    gf = GridFunction(fes)
    gf.Set(f)
    Draw1D(mesh1D, [(gf,"FE Approximation"), (f,"exakte Funktion")], n_p=20)
    for i in range(fes.ndof):
        print("{:5.2f}".format(gf.vec[i]),end=" ")
        if (i>0 and i%18==17):
            print("")

def PartialApproximateWithFESpace(N=64, order=3, f="sin(3·pi·x)", j = 2):

    i = j
    
    if f == "sin(3·pi·x)":
        f = Sin(3*pi*X)
    elif f == "sin(3·pi·x)·cos(5·pi·x)":
        f = Sin(3*pi*X)*Cos(5*pi*X)
    else:
        print("no function provided")
        return
        
    mesh1D = Mesh1D(N)
    try:
        f(mesh1D(0.5))
    except:
        print("expression for f not valid!")
        return
    
    fes = H1(mesh1D, order=order)

    if (i > fes.ndof):
        print( "j is too large. Setting j = ", fes.ndof-1)
        i = fes.ndof
    
    gf = GridFunction(fes)
    gf.Set(f)
    for j in range(i,fes.ndof):
        gf.vec[j] = 0
    Draw1D(mesh1D, [(gf,"FE Approximation"), (f,"exakte Funktion")], n_p=20)
    for j in range(i):
        print("{:5.2f}".format(gf.vec[j]),end=" ")
        if (j>0 and j%18==17):
            print("")



from ngsolve import IfPos, BilinearForm, LinearForm, H1, GridFunction, grad, SymbolicBFI, SymbolicLFI, Integrate



def ComputeMatrixEntry2(N=8, order=1, k=1, i=0, j=0 ):
    mesh1D = Mesh1D(N)
    fes = H1(mesh1D, order=order) #, dirichlet=)

    gf1 = GridFunction(fes)
    gf2 = GridFunction(fes)

    gf1.vec[:] = 0; gf1.vec[i] = 1.
    gf2.vec[:] = 0; gf2.vec[j] = 1.

    Draw1D(mesh1D,[(gf1,"phi_i"),(gf2,"phi_j")],n_p=5*order**2, figsize=(12,3.5))
    Draw1D(mesh1D,[(grad(gf1)[0],"dphi_idx"),(grad(gf2)[0],"dphi_jdx")],n_p=5*order**2, figsize=(12,3.5))
    print("A[i,j] = ", Integrate(gf1.Deriv()*gf2.Deriv(),mesh1D,order=2*(order-1)))

from functools import partial
    
def ComputeMatrixEntry(N=8, order=1, k=1):
    options = {
        "i" : IntSlider(min=0, max=N*order, step=1, continuous_update=True, description='i', value=0),
        "j" : IntSlider(min=0, max=N*order, step=1, continuous_update=True, description='j', value=0),
    }
    ComputeMatrixEntry3 = partial(ComputeMatrixEntry2,N,order,k)
    out = interactive_output(ComputeMatrixEntry3, options)
    ui = HBox([options["i"], options["j"]])
    display.display(ui, out)
    

import scipy.sparse as sp

import matplotlib.pylab as plt

def Spy(N=8, order=1):
    mesh1D = Mesh1D(N)
    fes = H1(mesh1D, order=order) #, dirichlet=)
    u,v = fes.TnT()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v))
    a.Assemble()
    rows,cols,vals = a.mat.COO()
    A = sp.csr_matrix((vals,(rows,cols)))
    plt.figure(figsize=(7,7))
    plt.spy(A)
    plt.show()

def SpyDG(N=8, order=1):
    mesh1D = Mesh1D(N)
    fes = Discontinuous(H1(mesh1D, order=order)) #, dirichlet=)
    u,v = fes.TnT()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v))
    a.Assemble()
    rows,cols,vals = a.mat.COO()
    A = sp.csr_matrix((vals,(rows,cols)))
    plt.figure(figsize=(7,7))
    plt.spy(A)
    plt.show()
    
    


def Heat1DFEM( N=8,
               order=1, 
               k1 = 1, 
               k2 = 1, 
               Q1=0, 
               Q2=10, 
               boundary_condition_left = "Robin", 
               boundary_condition_right = "Dirichlet", 
               value_left = 0, 
               value_right = 1, 
               q_value_left = 0, 
               q_value_right = 1, 
               r_value_left = 0, 
               r_value_right = 1,
               intervalsize = 0.14):
    if (boundary_condition_left == "Neumann" or (boundary_condition_left == "Robin" and r_value_left==0)) and (boundary_condition_right == "Neumann" or (boundary_condition_right == "Robin" and r_value_right==0)):
        print("Temperatur ist nicht eindeutig bestimmt.")
        #return
    mesh1D = Mesh1D(N,interval=(0,intervalsize))
    dbnds = []
    if boundary_condition_left == "Dirichlet":
        dbnds.append(1)
        # print(True)
    if boundary_condition_right == "Dirichlet":
        dbnds.append(2)
    
    fes = H1(mesh1D, order=order, dirichlet=dbnds)
    gf = GridFunction(fes)
        
    if boundary_condition_left == "Dirichlet":
        gf.vec[0] = value_left
    
    if boundary_condition_right == "Dirichlet":
        gf.vec[N] = value_right

       
        
    Q = IfPos(X-0.5*intervalsize,Q2,Q1)    
    k = IfPos(X-0.5*intervalsize,k2,k1)    
        
    u,v = fes.TnT()    
    a = BilinearForm(fes)    
    a += SymbolicBFI(k * grad(u) * grad(v))
    if boundary_condition_left == "Robin":
        a += SymbolicBFI( r_value_left * u * v, definedon = mesh1D.Boundaries("left"))
    if boundary_condition_right == "Robin":
        a += SymbolicBFI( r_value_right * u * v, definedon = mesh1D.Boundaries("right") )
    a.Assemble()
    f = LinearForm(fes)
    f += SymbolicLFI(Q * v)
    f.Assemble()
    if boundary_condition_left == "Neumann":
        f.vec[0] +=  q_value_left
    elif boundary_condition_left == "Robin":
        f.vec[0] +=  r_value_left * value_left
    if boundary_condition_right == "Neumann":
        f.vec[N] +=  q_value_right
    elif boundary_condition_right == "Robin":
        f.vec[N] +=  r_value_right * value_right
        
    f.vec.data -= a.mat * gf.vec
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * f.vec
    Draw1D(mesh1D,[(gf,"u_h")],n_p=5*order**2)            
    
def align(q):
    interactive_plot = q
    interactive_plot.layout.height = '500px'
    return interactive_plot

from ipywidgets import interact, interactive_output, interact_manual, interactive, FloatSlider, IntSlider, Dropdown, HBox, VBox

def Heat1DExample():
    options = {
        "N" : Dropdown(description='N', index=1, options=(2, 3, 4, 8, 16, 32, 64), value=8),
        "order" : Dropdown(description='order', index=1, options=(1,2,3,4), value=1),
        "k1" : FloatSlider(min=1, max=100, step=0.1, continuous_update=False, description=r'\(k_1 / [W/m K]\)', value=50),
        "k2" : FloatSlider(min=1, max=100, step=0.1, continuous_update=False, description=r'\(k_2 / [W/m K]\)', value=50),
        "Q1" : FloatSlider(min=0, max=1e7, step=1e4, continuous_update=False, description=r'\(q_Q^1 / [W/m^3]\)', value=30000),
        "Q2" : FloatSlider(min=0, max=1e7, step=1e4, continuous_update=False, description=r'\(q_Q^2 / [W/m^3]\)', value=30000),
        "boundary_condition_left" : Dropdown(description='bnd.cnd.left', index=1, options=("Dirichlet", "Neumann", "Robin"), value= "Dirichlet"),
        "boundary_condition_right" : Dropdown(description='bnd.cnd.right', index=1, options=("Dirichlet", "Neumann", "Robin"), value= "Dirichlet"),
        "value_left" : FloatSlider(min=0, max=100, step=0.5, description=r'\(T_l/ [°C]\)', continuous_update=False,value=20),
        "value_right" : FloatSlider(min=0, max=100, step=0.5, description=r'\(T_r/ [°C]\)', continuous_update=False,value=20),
        "q_value_left" : FloatSlider(min=0, max=100, step=0.5, description=r'\(q_S^l/ [W/m^3]\)', continuous_update=False),
        "q_value_right" : FloatSlider(min=0, max=100, step=0.5, description=r'\(q_S^r/ [W/m^3]\)', continuous_update=False),
        "r_value_left" : FloatSlider(min=0, max=200, step=0.1, description=r'\(\alpha_l / [W\!/m^2 K]\)', continuous_update=False, value=1),
        "r_value_right" : FloatSlider(min=0, max=200, step=0.1, description=r'\(\alpha_{r\!} / [W\!/m^2 K]\)', continuous_update=False, value=1),
        "intervalsize" : FloatSlider(min=0.01, max=1, step=0.01, description=r'\(l / [m]\)', continuous_update=False, value=1),
    }
    
    out = interactive_output(Heat1DFEM, options)
    
    ui = VBox([options["intervalsize"],
               HBox([options["N"], options["order"]]),
               HBox([options["k1"], options["k2"]]),
               HBox([options["Q1"], options["Q2"]]),
               HBox([options["boundary_condition_left"], options["boundary_condition_right"]]),
               HBox([options["q_value_left"], options["q_value_right"]]),
               HBox([options["value_left"], options["value_right"]]),
               HBox([options["r_value_left"], options["r_value_right"]])])
    
    display.display(ui, out)

from ngsolve import *
from ngsolve.webgui import Draw
class DrawBasisFunction2DiDrawer:

    def __init__(self, gfu):
        self.first=True
        self.gfu = gfu
    def Draw(self, i):
        if i >= self.gfu.space.ndof:
            return None
        self.gfu.vec[:]=0
        self.gfu.vec[i]=1
        if self.first:
            self.scene = Draw(self.gfu,self.gfu.space.mesh,"basis_fct",deformation=True)
            self.first = False
        else:
            self.scene.Redraw()
            

def DrawBasisFunction2D(gfu):
    drawer = DrawBasisFunction2DiDrawer(gfu)
    options = {
        "i" : IntSlider(min=0, max=gfu.space.ndof-1, step=1, continuous_update=False, description='i', value=0),
    }
    out = interactive_output(drawer.Draw, options)
    ui = HBox([options["i"]])
    display.display(ui, out)
    drawer.Draw(0)

def DrawOneBasisFunction(gfu,i):
    gfu.vec[:]=0
    gfu.vec[i]=1
    scene = Draw(gfu,gfu.space.mesh,"basis_fct",deformation=True)
    
