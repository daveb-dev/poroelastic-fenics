#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:21:15 2020
@author: ranjeet
"""

import dolfin as df

class MeshCreater:
    def __init__(self):
        self.mesh = None
        self.boundary_markers = None
        # Mark facets for Neumann BCs
        self.boundary_indices = {'left': 1,
                        'right': 2,
                        'top': 3,
                        'bottom': 4}
        self.cellType = None
    

    def UnitSqaureMesh(self, nx=16,ny=16, save=False):
        self.mesh = df.UnitSquareMesh(nx,ny)
        # Mark boundary subdomains
        left = df.CompiledSubDomain("near(x[0],0) && on_boundary")
        right = df.CompiledSubDomain("near(x[0],1.0) && on_boundary")
        bottom = df.CompiledSubDomain("near(x[1],0) && on_boundary")
        top = df.CompiledSubDomain("near(x[1],1.0) && on_boundary")
        
        self.__create_marker(left,right,bottom,top)
        
        
        
    def __create_marker(self, left,right,bottom,top):
        self.boundary_markers = df.MeshFunction("size_t", self.mesh, dim=1, value=0)
        left.mark(self.boundary_markers, self.boundary_indices["left"])
        right.mark(self.boundary_markers, self.boundary_indices["right"])
        top.mark(self.boundary_markers, self.boundary_indices["top"])
        bottom.mark(self.boundary_markers, self.boundary_indices["bottom"])
        
        

    def RectangleMesh(self, nx=16,ny=16, point2=(1,1), point1=(0,0), cellType="tri"): 
        "cellType = {'tri','quad'}"
        self.cellType = cellType
        x0,y0 = point1
        Lx, Ly = point2
        if cellType=="quad":
            self.mesh = df.RectangleMesh.create([df.Point(point1), df.Point(point2)], [nx, ny],df.CellType.Type.quadrilateral)
        else:
            self.mesh = df.RectangleMesh.create([df.Point(point1), df.Point(point2)], [nx, ny],df.CellType.Type.triangle)
    
        # Mark boundary subdomains
        left = df.CompiledSubDomain("near(x[0],x0) && on_boundary", x0=x0)
        right = df.CompiledSubDomain("near(x[0],Lx) && on_boundary",Lx=Lx)
        bottom = df.CompiledSubDomain("near(x[1],y0) && on_boundary",y0=y0)
        top = df.CompiledSubDomain("near(x[1],Ly) && on_boundary",Ly=Ly)
              
        self.__create_marker(left,right,bottom,top)
        
        
    def save_mesh(self, fileName="mesh"):
        
        if self.cellType=="quad":
            with df.XDMFFile( fileName+".xdmf") as f:
                f.write(mesh)
            with df.XDMFFile( fileName+"_facet.xdmf") as f:
                f.write(boundary_markers)
                
            f.close()
        else:
            df.File(fileName+".xml") << self.mesh
            df.File(fileName+"_facet.xml") << self.boundary_markers
        #     with df.XDMFFile("mesh_functions.xdmf") as f:
        #         f.write(boundary_markers)
        
    def print_bndry_info(self):
        print(self.boundary_indices)
        
    def get_mesh(self):
        return self.mesh
    
    def get_facet(self):
        return self.boundary_markers
        


# A simple test
#for i in range(1,5):
#    length = df.assemble(1.*ds(i))
#    print("side {:d} is of length {:4.2f}".format(i,length))
    
    
if __name__ == "__main__" :
    Nx = 30
    Ny = 60
    Nz = 85
    Lx = 182.88
    Ly = 182.88
    Lz = 51.816
    
    creater = MeshCreater()
    #creater.UnitSqaureMesh()
    creater.RectangleMesh(Nx,Ny, point2=(Lx,Ly))
    mesh = creater.mesh
    boundary_markers = creater.boundary_markers
    # Redefine element of area to include information about the markers
    ds = df.ds(domain=mesh, subdomain_data=boundary_markers)
    # dx = df.dx(domain=mesh)
    
    
    # # A simple test
    # mesh, markers, ds = mUnitSqaureMesh()
    # df.plot(mesh)
    
    for i in range(1,5):
        length = df.assemble(1.*ds(i))
        print("side {:d} is of length {:4.2f}".format(i,length))
        
    creater.save_mesh("spe10mesh")
    
    # # Reading Mesh
    # mesh2 = df.Mesh()
    # df.XDMFFile("spe10mesh.xdmf").read(mesh2)
    
    # # Read bndry marker
    # bndry = df.MeshFunction("size_t", mesh2, dim=1, value=0)
    # df.XDMFFile("spe10mesh_facet.xdmf").read(bndry)
    
    
        
        
    
        
    