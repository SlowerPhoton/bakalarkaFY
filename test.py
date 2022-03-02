import sympy as sp
from sympy.solvers.ode.systems import dsolve_system
from sympy import Function
e, ar, arp = sp.symbols("e,ar,arp", cls=Function)
t = sp.symbols("t")
k1, k2 = sp.symbols("k1,k2")

eq1 = sp.Eq(e(t).diff(t), k1*e(t)*ar(t) - k2*e(t)*ar(t)*arp(t))
eq2 = sp.Eq(ar(t).diff(t), -k1*e(t)*ar(t) + k2*e(t)*ar(t)*arp(t))
eq3 = sp.Eq(arp(t).diff(t), k1*e(t)*ar(t) - k2*e(t)*ar(t)*arp(t))

dsolve_system([eq1, eq2, eq3])
