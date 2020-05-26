import meshio

msh_file = "example02a.msh"

msh = meshio.read(msh_file)
meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]}))
meshio.write("mf.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                    cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))

from dolfin import *
mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
File("dolfinmesh.pvd").write(mesh)
File("dolfincellfunc.pvd").write(mf)
