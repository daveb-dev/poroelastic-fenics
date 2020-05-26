from fenics import *
import matplotlib.pyplot as plt
import numpy as np


mesh_file_name = "grid/example02a"
mesh = Mesh(mesh_file_name +".xml")
subdomains = MeshFunction("size_t", mesh, mesh_file_name+"_physical_region.xml")
bndry  = MeshFunction("size_t", mesh, mesh_file_name+"_facet_region.xml")


with XDMFFile('mesh.xdmf') as f:
    f.write(mesh)

with XDMFFile('facets.xdmf') as f:
    f.write(bndry)

with XDMFFile('subdomains.xdmf') as f:
    f.write(subdomains)


#plot(mesh)

#plt.show()



