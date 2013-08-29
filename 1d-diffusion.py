#! /usr/bin/python

from dolfin import *
import sys

N       = 32
kappa = 0.125
T       = 1
c       = 1.0
h = 1.0 / (N + 1)
dt = 5.0 * h  * h / kappa
iteration = 0

info("Initial parameters: N = {0}; kappa = {1}; h = {2}; dt = {3}".format(N, kappa, h, dt))
mesh = UnitIntervalMesh(N)
V = FunctionSpace(mesh, 'CG', 1)

u0 = Expression('((0.375 < x[0] && x[0] < 0.625)) ? 1 : 0')

def boundary(x, on_boundary):
  tol = 1E-15
  return abs(x[0]) < tol or \
         abs(x[0] - 1) < tol

bc = DirichletBC(V, Constant(0.0), boundary)

u_1 = interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

a = u*v*dx + kappa*dt*c*inner(nabla_grad(u), nabla_grad(v))*dx

A = assemble(a)

print A.array()

u = Function(V)
t = dt

viz_w = plot(u_1,
             range_min = 0.0,
             range_max = 1.0,
             rescale = False)
interactive()

while t <= T:
  L = u_1*v*dx - kappa*dt*(1-c)*inner(nabla_grad(u_1), nabla_grad(v))*dx
  b = assemble(L)
  bc.apply(A, b)
  solve(A, u.vector(), b, "bicgstab", "ilu")
  info("Iteration {0} (Time {1}): Min = {2}; Max = {3}".format(iteration, t, u.vector().min(), u.vector().max()))
  viz_w.plot(u)
  interactive()
  u_1.assign(u)
  t += dt
  iteration += 1
  c = 0.5
