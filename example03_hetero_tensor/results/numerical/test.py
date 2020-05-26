def nitsche1_poisson(N,p):
    'Classical (symmetric) Nitsche formulation.'
    function_name= 'Symmetric Nitsche Method'
    mesh = UnitSquareMesh(N, N)
    NT=mesh.num_cells()
    V = FunctionSpace(mesh, ’CG’, p)
    u = TrialFunction(V)
    v = TestFunction(V)

    beta_value = 2
    beta = Constant(beta_value)
    h_E = MinFacetEdgeLength(mesh)#mesh.ufl_cell().max_facet_edge_length
    n = FacetNormal(mesh)

    a = inner(grad(u), grad(v))*dx - inner(dot(grad(u), n), v)*ds -\
            inner(u, dot(grad(v), n))*ds + beta*h_E**-1*inner(u, v)*ds

    L = inner(f, v)*dx -\
            inner(u_exact , dot(grad(v), n))*ds + beta*h_E**-1*inner(u_exact , v)*ds

    uh = Function(V)
    solve(a == L, uh)

    # Save solution to file
    #file = File("poisson.pvd")
    #file << uh
    if N==64:
        viz1=plot(uh)
    viz1.write_png(’snitche’)
    # plot(uh, title=’numeric’)
    # plot(u_exact , mesh=mesh, title=’exact’)
    # interactive()

    # Compute norm of error
    E = FunctionSpace(mesh, ’DG’, 4)
    uh = interpolate(uh, E)
    u = interpolate(u_exact , E)
    e = uh - u

    norm_L2 = assemble(inner(e, e)*dx, mesh=mesh)
    norm_H10 = assemble(inner(grad(e), grad(e))*dx, mesh=mesh)
    norm_edge = assemble(beta*h_E**-1*inner(e, e)*ds)

    norm_H1 = norm_L2 + norm_H10 + norm_edge
    norm_L2 = sqrt(norm_L2)
    norm_H1 = sqrt(norm_H1)
    norm_H10 = sqrt(norm_H10)

    return Result(fun_name=function_name ,NT=NT,h=mesh.hmin(), L2=norm_L2 ,
            H1=norm_H1 , H10=norm_H10)


def nitsche2_poisson(N,p):
    'Unsymmetric Nitsche formulation.'
    function_name= 'Unsymmetric Nitche Method'
    mesh = UnitSquareMesh(N, N)
    NT=mesh.num_cells()

    V = FunctionSpace(mesh, ’CG’, p)
    u = TrialFunction(V)
    v = TestFunction(V)

    beta_value = 2
    beta = Constant(beta_value)
    h_E = MinFacetEdgeLength(mesh)#mesh.ufl_cell().max_facet_edge_length
    n = FacetNormal(mesh)

    a = inner(grad(u), grad(v))*dx - inner(dot(grad(u), n), v)*ds +\
            inner(u, dot(grad(v), n))*ds + beta*h_E**-1*inner(u, v)*ds

    L = inner(f, v)*dx +\
            inner(u_exact , dot(grad(v), n))*ds + beta*h_E**-1*inner(u_exact , v)*ds

    uh = Function(V)
    solve(a == L, uh)

    # plot(uh, title=’numeric’)
    # plot(u_exact , mesh=mesh, title=’exact’)
    # interactive()
    if N==64:
        viz1=plot(uh)
    viz1.write_png(’nnitche’)
    # Compute norm of error
    E = FunctionSpace(mesh, ’DG’, 4)
    uh = interpolate(uh, E)
    u = interpolate(u_exact , E)
    e = uh - u

    norm_L2 = assemble(inner(e, e)*dx, mesh=mesh)
    norm_H10 = assemble(inner(grad(e), grad(e))*dx, mesh=mesh)
    norm_edge = assemble(beta*h_E**-1*inner(e, e)*ds)
    norm_H1 = norm_L2 + norm_H10 + norm_edge

    norm_L2 = sqrt(norm_L2)
    norm_H1 = sqrt(norm_H1)
    norm_H10 = sqrt(norm_H10)

    return Result(fun_name=function_name ,NT=NT,h=mesh.hmin(), L2=norm_L2 ,H1=norm_H1 , H10=norm_H10)

#
#-----------------------------------------------------------------------------

methods = [classical_poisson , nitsche1_poisson , nitsche2_poisson]

#print ’The Method:{:d}’,format( method
for m in [1,2]: #[0,1,2]:
    for p in [1,2]:
        method = methods[m]

#print "The Method:{}".format(method)

#norm_type = ’H1’
R = method(N=4,p=p)
print "The method: {0} with degree of polynomial {1:}".format(R.fun_name ,p)
print "{0:>6s} {1:>6s} {2:>10s} {3:>7s} {4:>6s} {5:>7s}".format("Elem", "h", "L^2", "rate", "H^1","rate")
h_ = R.h
e_ = getattr(R, ’H1’)
eL2_= getattr(R,’L2’)
for N in [8, 16, 32, 64,128,256,512]:
    R = method(N,p=p)
h = R.h
NT = R.NT
e = getattr(R, ’H1’)
eL2 = getattr(R,’L2’)
rate = ln(e/e_)/ln(h/h_)
rateL2 = ln(eL2/eL2_)/ln(h/h_)
# print ’h error rate’
print({0:6d} {h:.3E} {eL2:.5E} {rateL2:.4f} {e:.5E} {rate:.4f}’.format(NT,h=h,eL2=eL2,rateL2=rateL2 , e=e, rate=rate))
