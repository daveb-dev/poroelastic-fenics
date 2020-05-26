from dolfin import *
import numpy as np
import os
dirname = os.path.dirname(os.path.abspath( __file__ )) 

# Step 1: Store the data into xml files 

Lx = 365.76
Ly = 670.56
Nx = 60
Ny = 220
dx = Lx/Nx
dy = Ly/Ny
mesh = RectangleMesh.create([Point(0,0), Point(Lx,Ly)], [Nx, Ny], CellType.Type.quadrilateral)
#File(dirname+"/grids/mesh.xml") << mesh
V = FunctionSpace(mesh, "DQ", 0)


def create_SPE10_slice2D(Nx, Ny, x_shift = 0, y_shift = 0, z_shift = 0, perm_factor = 1.0):
    import os
    dirname = os.path.dirname(os.path.abspath( __file__ ))        
    kk = z_shift #layer
    lines = np.loadtxt(dirname + "/data/spe_phi.dat")
     
    single_line = np.reshape(lines, 60*220*85)

    phi = np.zeros((Nx,Ny))
    for j in range(0,Ny):
        for i in range(0,Nx):
            phi[i][j] = single_line[(i+x_shift)+(j+y_shift)*60+kk*220*60]


    # load data for permeability
    lines = np.loadtxt(dirname + "/data/spe_perm.dat")
    print("Using SPE10 permeability/porosity fields")
    print("Permeability factor is: ", perm_factor)
    lines =  lines*9.869233e-10*perm_factor #converting md to mm^2  

    single_line = np.reshape(lines, [3, 60*220*85]).transpose()
    kxx = np.zeros((Nx,Ny))
    for j in range(0,Ny):
        for i in range(0,Nx):
            kxx[i][j] = single_line[(i+x_shift)+(j+y_shift)*60+kk*220*60,0]
    
    kyy = np.zeros((Nx,Ny))
    for j in range(0,Ny):
        for i in range(0,Nx):
            kyy[i][j] = single_line[(i+x_shift)+(j+y_shift)*60+kk*220*60,1]
            
                    
    coords = V.tabulate_dof_coordinates()
    def coordtoij(x,y):
        i = np.floor(x / dx).astype(int)
        j = np.floor(y / dy).astype(int)
        return i,j
        
    i, j = coordtoij(coords[:, 0], coords[:, 1])  
            
    import pandas as pd
    df = pd.DataFrame({"phi":phi[i,j], "kxx":kxx[i, j], "kyy":kyy[i, j]})
    df.to_csv("data/data.csv", index=False)
    
    df2 = df/(1000*9.81)
    df2.to_csv("data/data2.csv", index=False)
    
    print("Data written succussfully...!")
    
    
if __name__ == "__main__":
    create_SPE10_slice2D(Nx,Ny)
    
    
    
 

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
