from firedrake import *

# Set numerical constants
# n = 30
START_TIME = 0
END_TIME = 20000
TIMESTEP = 50
MESH = RectangleMesh(50, 25, 10000, 5000)
NU_MAX = 1.0
OUTPUT_TIME = 500
OUTPUT_FOLDER = "output/"
INFLOW = 0
TINY = 1e-10
# FUNC_NAMES = (
    # "sed_old",
    # "phi",
    # "surface",
    # "limiter",
    # "thickness",
    # "depth",
    # "diff",
# )

# Initialise our blank functions
V = FunctionSpace(MESH, "CG", 1)
sed_old = Function(V, name="sed_old")
phi = Function(V, name="sed")
surface = Function(V, name="surface")
limiter = Function(V, name="limiter")
thickness = Function(V, name="thickness")
depth = Function(V, name="depth")
diff = Function(V, name="diff")
slf = Function(V, name="sea_level")
# funcs = {func_name: Function(V, name=func_name) for func_name in FUNC_NAMES}
layers = []

v = TestFunction(V)

x = SpatialCoordinate(MESH)
ic = project(
    (
        20000
        * (1 / (2*sqrt(2*3.14*250*250)))
        * exp(
            -((x[0]-6000) * (x[0]-6000))
            / (2*250*250)
        )
    )
    + (
        50000
        * (1 / (2*sqrt(2*3.14*1000*1000)))
        * exp(
            -((x[0]-4000) * (x[0]-4000))
            / (2*1000*1000)
        )
    ), V)

sl_time = Constant(25*sin(START_TIME/100000/180*3.14159))
sea_level = Function(V, val=interpolate(sl_time, V), name="sea_level")
# if we don't do this, the sea level ouput has a random function name.
# Not helpful...
slf.interpolate(sea_level)

land = project(100*(tanh(0.0005*(x[0]-6000))), V, name="Starting_topo")
outfile = File(OUTPUT_FOLDER + "land.pvd")
outfile.write(land)

sed_old.assign(ic)
phi.assign(ic)

surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
thickness.interpolate(surface-land)
limiter.interpolate((surface - land)/(surface-land+TINY))
depth.interpolate(sea_level - surface)
diff.project(NU_MAX*((2/(sqrt(2*3.14159))*exp(-0.5*depth**2))+0.2022))

F = (inner((phi - sed_old)/TIMESTEP, v)
     + limiter*diff*inner(grad(phi+land), grad(v)))*dx

# We now create an object for output visualisation::

geometry = File(OUTPUT_FOLDER + "surfaces.pvd")
geometry.write(surface, phi, time=0)
layer_data = File(OUTPUT_FOLDER + "layer_data.pvd")
layer_data.write(diff, thickness, depth, time=0)
sl = File(OUTPUT_FOLDER + "sea_level.pvd")
sl.write(slf, time=0)

current_time = START_TIME
while (current_time <= END_TIME):
    solve(F == 0, phi)
    sed_old.assign(phi)
    current_time += TIMESTEP
    limiter.interpolate((surface - land)/(surface-land+TINY))
    surface.interpolate(
        (
            ((land+phi) + land)
            + abs((land+phi)-land)
        )
        / 2
    )
    depth.interpolate(sea_level - surface)
    diff.interpolate(
        NU_MAX
        * (
            (2/(sqrt(2*3.14159))*exp(-0.5*((depth-5)/10)**2))
            + 0.2022
        )
    )
    thickness.interpolate(surface-land)
    sl_time.assign(25*sin(current_time/100000*180/3.14159))
    sea_level.interpolate(sl_time)
    slf.interpolate(sea_level)
    if current_time % OUTPUT_TIME == 0:
        layer_data.write(diff, thickness, depth, time=current_time)
        geometry.write(surface, phi, time=current_time)
        sl.write(slf, time=current_time)
        layers.append([surface, phi, depth, diff, thickness])
