from netgen.meshing import Mesh as NGMesh
from netgen.meshing import MeshPoint, Element1D, Element0D, Pnt
from ngsolve import Mesh as NGSMesh

def Mesh1D(elements, interval=(0,1), periodic=False):
    """
        generate an equidistant 1D mesh with N cells
    """
    N = elements
    left = interval[0]
    right = interval[1]
    mesh = NGMesh(dim=1)
    pids = []
    for i in range(N+1):
        pids.append (mesh.Add (MeshPoint(Pnt(left + i/N * (right-left), 0, 0))))
    for i in range(N):
        mesh.Add(Element1D([pids[i],pids[i+1]],index=1))
    mesh.Add (Element0D( pids[0], index=1))
    mesh.Add (Element0D( pids[N], index=2))
    mesh.SetBCName(0,"left")
    mesh.SetBCName(1,"right")
    if periodic:
        mesh.AddPointIdentification(pids[0],pids[N],1,2)
    ngsmesh = NGSMesh(mesh) 
    return ngsmesh
