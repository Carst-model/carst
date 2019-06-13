from firedrake import *

n = 30
mesh = UnitSquareMesh(n, n)
nu_max = 0.007
t = 0

# We choose degree 2 continuous Lagrange polynomials. We also need a
# piecewise linear space for output purposes::

V = FunctionSpace(mesh, "CG", 1)

# We also need solution functions for the current and the next
# timestep. Note that, since this is a nonlinear problem, we don't
# define trial functions::

phi_ = Function(V, name="SedOld")
phi = Function(V, name="Sed")
surface = Function(V, name="surface")
limit = Function(V,name="limiter")
thickness = Function(V, name="thickness")
sea_level = Function(V, name="sea_level")
depth = Function(V, name="depth")
nu = Function(V, name="diff")

v = TestFunction(V)

# For this problem we need an initial condition::

x = SpatialCoordinate(mesh)
ic = project((1.0/(2.0*sqrt(2.*3.14*0.1*0.1)) * exp(-((x[0]-0.9)*(x[0]-0.9))/(2*0.1*0.1)))/5.0, V)
#ic = project(0.25*sin(x[0]*3.1459),V)
#ic = project(Expression(0.8),V)
#land = project(Expression('0.0'), V)
sl_time = Constant(0.5*cos((t/50)*180/3.14159))
sea_level = project(sl_time, V)

land = project(-sin(x[0]*3.14159), V)
outfile = File("land.pvd")
outfile.write(land)


#ic = project(x[0],V)
# We start with current value of u set to the initial condition, but we
# also use the initial condition as our starting guess for the next
# value of u::

phi_.assign(ic)
phi.assign(ic)

tiny = 1e-10
surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
thickness.interpolate(surface-land)
limit.interpolate((surface - land)/(surface-land+tiny))
depth.interpolate(sea_level - surface)
nu.project(nu_max*((2/(sqrt(2*3.14159))*exp(-0.5*depth**2))+0.2022))

# The timestep is set to produce an advective Courant number of
# around 1. Since we are employing backward Euler, this is stricter than
# is required for stability, but ensures good temporal resolution of the
# system's evolution::

timestep = 1.0/n
inflow = 0

# Here we finally get to define the residual of the equation. In the advection
# term we need to contract the test function :math:`v` with 
# :math:`(u\cdot\nabla)u`, which is the derivative of the velocity in the
# direction :math:`u`. This directional derivative can be written as
# ``dot(u,nabla_grad(u))`` since ``nabla_grad(u)[i,j]``:math:`=\partial_i u_j`.
# Note once again that for a nonlinear problem, there are no trial functions in
# the formulation. These will be created automatically when the residual
# is differentiated by the nonlinear solver::

F = (inner((phi - phi_)/timestep, v)
     + limit*nu*inner(grad(phi+land), grad(v)))*dx

# We now create an object for output visualisation::

outfile = File("diffusion.pvd")
tot = File("surface.pvd")
outfile.write(phi, limit, time=0)
sl = File("Sea_level.pvd")


# Finally, we loop over the timesteps solving the equation each time and
# outputting each result::

end = 10
while (t <= end):
    solve(F == 0, phi)#, bcs=[bc_land, bc_sea])
    phi_.assign(phi)
    t += timestep
    outfile.write(phi, limit, time=t)
    limit.interpolate((surface - land)/(surface-land+tiny))
    surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
    depth.interpolate(sea_level - surface)
    nu.interpolate(nu_max*((2/(sqrt(2*3.14159))*exp(-0.5*depth**2))+0.2022))
    thickness.interpolate(surface-land)
    sl_time.assign(0.5*cos((t/50)*180/3.14159))
    sea_level.interpolate(sl_time)
    tot.write(surface, thickness, nu, time=t)
    sl.write(sea_level, time=t)
