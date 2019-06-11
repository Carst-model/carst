from firedrake import *

n = 30
mesh = RectangleMesh(50, 25, 10000, 5000)
nu_max = 1.0
t = 0
output_time = 500
carb_pr = 5.0 # mm/yr

V = FunctionSpace(mesh, "CG", 1)

phi_ = Function(V, name="SedOld")
phi = Function(V, name="Sed")
surface = Function(V, name="surface")
limit = Function(V,name="limiter")
thickness = Function(V, name="thickness")
depth = Function(V, name="depth")
nu = Function(V, name="diff")
slf = Function(V, name="sea_level")
light_func = Function(V, name="light_attenuation")
layers = []

v = TestFunction(V)

x = SpatialCoordinate(mesh)
ic = project((20000*(1/(2*sqrt(2*3.14*250*250)) * exp(-((x[0]-6000)*(x[0]-6000))/(2*250*250)))) + (50000*(1/(2*sqrt(2*3.14*1000*1000)) * exp(-((x[0]-4000)*(x[0]-4000))/(2*1000*1000)))), V)
sl_time = Constant(25*sin(t/1000000/180*3.14159))
sea_level = Function(V, val=interpolate(sl_time,V), name="sea_level")
# if we don't do this, the sea level ouput has a random function name. 
# Not helpful...
slf.interpolate(sea_level)


land = project(100*(tanh(0.0005*(x[0]-6000))), V, name="Starting_topo")
outfile = File("land.pvd")
outfile.write(land)


#ic = project(x[0],V)
phi_.assign(ic)
phi.assign(ic)

tiny = 1e-10
surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
thickness.interpolate(surface-land)
limit.interpolate((surface - land)/(surface-land+tiny))
depth.interpolate(sea_level - surface)
nu.project(nu_max*((2/(sqrt(2*3.14159))*exp(-0.5*depth**2))+0.2022))

timestep = 50.0
inflow = 0

F = (inner((phi - phi_)/timestep, v)
     + limit*nu*inner(grad(phi+land), grad(v)))*dx

# We now create an object for output visualisation::

geometry = File("surfaces.pvd")
geometry.write(surface, phi, time=0)
layer_data = File("layer_data.pvd")
layer_data.write(nu, thickness, depth, time=0)
sl = File("Sea_level.pvd")
sl.write(slf, time=0)

end = 50000
while (t <= end):
    # grow carbs
    light_func.interpolate(1.0/(1+exp(-2*depth*25.0))* exp(-1.0*depth/10))
    solve(F == 0, phi)
    phi += carb_pr * light_func
    phi_.assign(phi)
    t += timestep
    limit.interpolate((surface - land)/(surface-land+tiny))
    surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
    depth.interpolate(sea_level - surface)
    nu.interpolate(nu_max*((2/(sqrt(2*3.14159))*exp(-0.5*((depth-5)/10)**2))+0.2022))
    thickness.interpolate(surface-land)
    sl_time.assign(25*sin(t/100000*180/3.14159))
    sea_level.interpolate(sl_time)
    slf.interpolate(sea_level)
    if t%output_time == 0:
        layer_data.write(nu, thickness, depth, time=t)
        geometry.write(surface, phi, time=t)
        sl.write(slf, time=t)
        layers.append([surface, phi, depth, nu, thickness])
