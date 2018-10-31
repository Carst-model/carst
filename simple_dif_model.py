from firedrake import *
import sys

class SedimentModel:
    """A class to track a single sediment type
       on a non-uniform domain, using diffusion only
       to move sediment around
    """

    def __init__(self):
        """Set up sensible initial values
        """
        self.V = None
        self.mesh = None
        self.s0 = None
        self.end_time = 1
        self.dt = 1
        self.alpha = Constant(1)
        self.inflow_rate = Expression('0')
        self.output_time = 100000
        self.L = None
        self.s1 = None
        self.s = None

    def set_mesh(self,mesh):
        """Set which mesh to use and define function space"""
        self.mesh = mesh
        self.V = FunctionSpace(mesh, "CG", 1)

    def set_initial_conditions(self,topography,sediment):
        # Set the initial topography and sediment with
        # UFL expressions
        self.h0 = topography
        self.s0 = sediment

    def set_timestep(self,dt):
        self.dt = dt

    def set_end_time(self,time):
        self.end_time = time

    def init(self,plot_init=False):
        """Initialise the solvers, etc, and setup the problem
        """

        tiny = 1e-10

        # initial guess of solution
        s_1 = Function(self.V)
        self.s = Function(self.V)
        v = TestFunction(self.V)
        self.s_1.assign(self.s0)
        self.s.assign(self.s0)
        h = Function(V)
        h.assign(self.h0)
        f = Constant(0)
        
        # This limiter stops the diffuive process where the amount of sediment (s) is zero
        self.limit = (s_1 + abs(s_1))/(2*s_1 + tiny) # limiting term

        
        # RHS
        a = (s)*v*dx + limit*self.dt*inner(nabla_grad(s+h), nabla_grad(v))*self.alpha*dx
        # LHS
        self.L = (self.s_1 + self.dt*f)*v*dx - self.inflow_rate*v*ds
        
        F = a - L
        self.b = None
        self.A = assemble(a)   # assemble only once, before the time stepping


    def solve(self):
        """Solve the problem"""
        t = 0
        while t <= self.end_time:
            self.b = assemble(self.L)
            solve(self.A, self.s, self.b)

            self.b = assemble(self.L, self.b, bcs=bc)
            t += self.dt
            self.s_1.assign(self.s)

    def get_total_height_array(self):
        return self.s_1.vector().array()+self.h.vector().array()

    def get_total_height(self):
        return self.s_1+self.h

    def get_sed_height_array(self):
        return self.s_1.vector().array()

    def get_sed_height(self):
        return self.s_1

    def get_topographic_height_array(self):
        return self.h_1.vector().array()

    def get_topographic_height(self):
        return self.h_1

    def set_diffusion_coeff(self,coeff):
        self.alpha = Constant(coeff)

if __name__ == "__main__":
    
    #create a simple testcase
    model = SedimentModel()
    mesh = SquareMesh(10,10,1)
    model.set_mesh(mesh)
    init_cond = Expression('x/10.') # simple slope
    init_sed = Expression('1.0/(2.0*sqrt(2.*3.14*0.4*0.4)) * exp(-((x-6)*(x-6))/(2*0.4*0.4))') # this gives
    # a Gaussian bump centred on 6, with height of ~0.5m
    model.set_initial_conditions(init_cond,init_sed)
    model.set_end_time(10)
    model.set_diffusion_coeff(0.0001)
    model.solve()
    # answer should be 1 everywhere
    print(model.get_sed_height_array())


