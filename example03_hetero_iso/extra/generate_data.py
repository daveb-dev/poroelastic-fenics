from dolfin import *
import numpy as np
import os
from utilities import *
#dirname = os.path.dirname(os.path.abspath('__file__')) 

# Step 1: Store the data into xml files 


mesh_file_name = "grids/grid1/example01"
mesh = Mesh(mesh_file_name +".xml")
subdomains = MeshFunction("size_t", mesh, mesh_file_name+"_physical_region.xml")
bndry  = MeshFunction("size_t", mesh, mesh_file_name+"_facet_region.xml")



# Create mesh functions for phi, kxx, kyy
phi = MeshFunction("double", mesh, mesh.topology().dim())
kxx = MeshFunction("double", mesh, mesh.topology().dim())
kyy = MeshFunction("double", mesh, mesh.topology().dim())


p = gaussianField((0.2,0.4), mesh.num_cells())
# generate permeability
tau_ = 0.81
dp = 10e-6
perm =  p**3*(dp)**2/(tau_*72*(1.- p)**2)
perm = perm.astype('float32') 

phi.array()[:] = p[:]
kxx.array()[:] = perm[:]
kyy.array()[:] = perm[:]

for cell in cells(mesh):
    index = cell.index()
    phi[cell] = float(p[index])
    kxx[cell] = perm[index]
    kyy[cell] = perm[index]
    


 
# Store to file
mesh_file = File("data/"+mesh_file_name +".xml")
phi_file =  File("data/phi.xml.gz")
kxx_file =  File("data/kxx.xml.gz")
kyy_file =  File("data/kyy.xml.gz")

mesh_file << mesh
phi_file << phi
kxx_file << kxx
kyy_file << kyy

#plot(phi)
#plot(kxx)
#plot(kyy)

# Save vtk file
phi.rename("phi","phi")
kxx.rename("kxx", "kxx")
kyy.rename("kyy", "kyy")
File("data/phi.pvd") <<phi
File("data/kxx.pvd") << kxx
File("data/kyy.pvd") << kyy






    
    
    
 

# # 1.a Create mesh functions for c00, c01, c11
# phi = MeshFunction("double", mesh, 2)
# kxx = MeshFunction("double", mesh, 2)
# kyy = MeshFunction("double", mesh, 2)

# # # 1.b Iterate over mesh and set values
# # for cell in cells(mesh):
# #     if cell.midpoint().x() < 0.5:
# #                 c00[cell] = 1.0
# #                 c01[cell] = 0.3
# #                 c11[cell] = 2.0
# #     else:
# #                 c00[cell] = 3.0
# #                 c01[cell] = 0.5
# #                 c11[cell] = 4.0

# # Reading from text file
# #poro_txt = np.ones((8,8)).flatten() # read from the file
# poro_txt = np.ones(mesh.num_cells()) # read from the file
# # Create Mesh Function
# # Iterate through each cell and update the MeshFunction object
# # for cell in cells(mesh):
#     # Update MeshFunction object for each cell
#     #kxx[cell] = U[cell.index]                 

# for cell in cells(mesh):
#     phi[cell] = poro_txt[int(cell.index())]

# #1.c Store to file
# # mesh_file = File("mesh.xml.gz")
# # c00_file = File("c00.xml.gz")
# # c01_file = File("c01.xml.gz")
# # c11_file = File("c11.xml.gz")
    


# # mesh_file = XDMFFile(MPI.comm_world,"mesh.xdmf")
# # phi_out = File(MPI.comm_world,"phi.xml.gz")
# # phi_out << phi

# # phi_out = HDF5File(MPI.comm_world, "phi.h5", "w")

# # kx_out  = XDMFFile(MPI.comm_world, "kx.xdmf")
# # ky_out  = XDMFFile(MPI.comm_world, "ky.xdmf")



# # #mesh_file << mesh
# # phi_out.write(f, "/phi")
# # #phi_out.write_checkpoint(c00, "phi", XDMFFile.Encoding.HDF5, False)
# # #c01_file << c01
# # #c11_file << c11               

# # phi_out.close()
